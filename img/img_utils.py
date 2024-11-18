""" ImagesCombiner operates on a set of FIT images thru single line python statements

frames are first loaded into memory using the Images class methods.
then operations are applied in sequence on a single python line.

eg. : to generate a master bias image : 
master_bias = Images.from_fit(dir = "../CAPTURE/test01/", filter = "offset-*.fit") \
                    .trim('600, 600, 2700, 1400') \
                    .sigmaclip() 
                    .offset(1500 * u.adu)
                        
eg. : to reduce spectra raw frames : 
master_sciences = Images.from_fit(dir = "../CAPTURE/test01/", filter = "agdra-*.fit") \
                        .trim('600, 600, 2700, 1400' ) \
                        .reduce(master_bias, master_dark, master_flat, 'EXPTIME') \
                        .spec_align()
                  
"""
from typing import List
from pathlib import Path
import warnings
import numpy as np
import logging
from astropy.io import fits
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty
from astropy.visualization import astropy_mpl_style, quantity_support
from astropy.utils.exceptions import AstropyWarning
from ccdproc import combine, subtract_bias, subtract_dark, flat_correct
from ccdproc import trim_image, Combiner, ccd_process, cosmicray_median, cosmicray_lacosmic, create_deviation
from ccdproc import ImageFileCollection, gain_correct
from astropy.stats import mad_std
from astropy.nddata.blocks import block_reduce
from astropy import visualization as aviz

#import astroalign as aa
from scipy.signal import fftconvolve

from app.config import Config

warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', UserWarning)

class ImagesCombiner(object):

    """
    maintains images set array and file names
    max memory is used by ccdproc routines to avoid OOM exceptions when working with large set of big images
    """
    def __init__(self, images: List[CCDData], names: List[str], max_memory: float = 1e9):
        self._images: List[CCDData] = images
        self._image_names: List[str] = names
        self._memory_limit: float = max_memory
    """
    returns a specific image array thru its index
    """
    def __getitem__(self, i:int) -> np.ndarray:
        return self._images[i]

    """
    returns the number of frames loaded in this set
    """
    def __len__(self) -> int:
        return len(self._images)
    
    def get_image_names(self) -> List[str]:
        return self._image_names

    """
    returns the sum frame of frames loaded in this set
    """
    def sum(self) -> CCDData:
        logging.info(f'sum combine on {len(self._images)} images ...')
        return(combine(self._images,
                       method = 'sum',
                       dtype = np.float32,
                       mem_limit = self._memory_limit)
              )
   
    """
    returns the sigmaclip'ed frame of frames loaded in this set
    """
    def sigmaclip(self, low_thresh: int = 5, high_thresh: int = 5) -> CCDData:
        logging.info(f'sigmaclip combine on {len(self._images)} images ...')
        return(combine(self._images,
                       method = 'average',
                       sigma_clip = True, 
                       sigma_clip_low_thresh = low_thresh, 
                       sigma_clip_high_thresh = high_thresh,
                       sigma_clip_func = np.ma.median, 
                       signma_clip_dev_func = mad_std, 
                       dtype = np.float32,
                       mem_limit = self._memory_limit)
              )

    """
    returns the median frame of frames loaded in this set
    """
    def median(self) -> CCDData:
        logging.info(f'median combine on {len(self._images)} images ...')
        return (combine(self._images, 
                        method = 'median', 
                        dtype = np.float32, 
                        mem_limit = self._memory_limit)
               )

    """
    trim all frames loaded in this set
    trim_region is the rectangle : x1, y1, x2, y2 
    """
    def trim(self, trim_region: str | None):
        if trim_region is not None:
            for i in range(0, len(self._images)):
                self._images[i] = trim_image(self._images[i][eval(trim_region)[1]:eval(trim_region)[3], eval(trim_region)[0]:eval(trim_region)[2]])

            logging.info(f'{len(self._images)} images trimmed to ({trim_region})')
        else:
            logging.info('no trimming')
        return self

    """
    add a scalar to all frames loaded in this set
    """
    def offset(self, scalar):
        for i in range(0, len(self._images)):
            self._images[i] = CCDData(CCDData.add(self._images[i], scalar), header = self._images[i].header) #, unit = self._images[i].unit
            
        logging.info(f'{len(self._images)} images added by ({operand})')
        return self
    """
    substract a master bias frame to all frames loaded in this set
    """
    def bias_substract(self, frame):
        for i in range(0, len(self._images)):
            self._images[i] = subtract_bias(self._images[i], frame)
            
        logging.info(f'masterbias substracted to {len(self._images)} images')
        return self

    """
    substract a master dark frame to all frames loaded in this set
    """
    def dark_substract(self, frame, scale_exposure: bool = True, exposure = 'EXPTIME'):
        for i in range(0, len(self._images)):
            self._images[i] = subtract_dark(self._images[i], frame, scale = scale_exposure, exposure_time = exposure, exposure_unit = u.second)                
        
        logging.info(f'masterdark substracted to {len(self._images)} images')
        return self
    
    """
    divide a master flat frame to all frames loaded in this set
    """
    def flat_divide(self, frame):
        for i in range(0, len(self._images)):
            self._images[i] = flat_correct(ccd = self._images[i], flat = frame, min_value = None, norm_value = 10000 * u.adu)

        logging.info(f'masterflat divided to {len(self._images)} images')
        return self

    """
    process science frames
    masterdark frame is scaled according to science frame exposure duration
    """
    def reduce(self, master_bias, master_dark, master_flat, exposure_key = 'EXPTIME'):
        for i in range(0, len(self._images)):
            self._images[i] = ccd_process(ccd = self._images[i], 
                oscan = None, 
                gain_corrected = True, 
                trim = None, 
                error = False,
#                gain = camera_electronic_gain*u.electron/u.adu ,
#                readnoise = camera_readout_noise*u.electron,
                master_bias = master_bias,
                dark_frame = master_dark,
                master_flat = master_flat,
                exposure_key = exposure_key,
                exposure_unit = u.second,
                dark_scale = True)            

        logging.info(f'{len(self._images)} images reduced')
        return self

    """
    align a set of loaded frames - specific to spectra fields (fft based)
    """
    def spec_align(self, ref_image_index: int = 0):    
        ### Collect arrays and crosscorrelate all (except the first) with the first.
        logging.info('align: fftconvolve running...')
        nX, nY = self._images[ref_image_index].shape
        correlations = [fftconvolve(self._images[ref_image_index].data.astype('float32'),
                                    image[::-1, ::-1].data.astype('float32'),
                                    mode='same') 
                        for image in self._images[1:]]
    
        ### For each image determine the coordinate of maximum cross-correlation.
        logging.info('align: get max cross-correlation for every image...')
        shift_indices = [np.unravel_index(np.argmax(corr_array, axis=None), corr_array.shape) 
                         for corr_array in correlations]
        
        deltas = [(ind[0] - int(nX / 2), ind[1] - int(nY / 2)) for ind in shift_indices]
        logging.info('align: images deltas = ' + repr(deltas))
    
        ### Warn for ghost images if realignment requires shifting by more than
        ### 15% of the field size.
        x_frac = abs(max(deltas, key=lambda x: abs(x[0]))[0]) / nX
        y_frac = abs(max(deltas, key=lambda x: abs(x[1]))[1]) / nY
        t_frac = max(x_frac, y_frac)
        if t_frac > 0.15:
            logging.warning('align: shifting by {}% of the field size'.format(int(100 * t_frac)))
    
        ### Roll the images to realign them and return their median.
        logging.info('align: images realignement ...')
        realigned_images = [CCDData(np.roll(image, deltas[i], axis=(0, 1)).data.astype('float32'), unit = u.adu, header = image.header) 
                            for (i, image) in enumerate(self._images[1:])]

        ### do not forget the reference image
        realigned_images.append(CCDData(self._images[ref_image_index].data.astype('float32'), unit = u.adu, header = self._images[ref_image_index].header))
        logging.info('align: complete')
        return ImagesCombiner(realigned_images, [])

"""
Images class implements the file loader methods
"""
class Images(ImagesCombiner):
    def __init__(self, images: List[CCDData]):
        ImagesCombiner.__init__(self, images, [])

    """
    collect and sort (according to fit header 'date-obs') file names using a wildcard filter
    """
    @classmethod
    def find_files(cls, directory: str, files_filter: str, sort_key: str = 'date-obs') -> List[str]:
        ic = ImageFileCollection(directory, glob_include = files_filter)
        ic.sort([sort_key])
        return (ic.files_filtered(include_path=True))

    """
    load a set of FIT images
    create deviation metadata from the gain & readnoise provided
    """
    @classmethod
    def from_fit(cls, dir: str, filter: str, 
                 camera_electronic_gain: float = 1.2 * u.electron / u.adu, 
                 camera_readout_noise: float =  2.2 * u.electron):
        
        images = []
        names: List[str] = []

        for fp in Images.find_files(directory = dir, files_filter = filter):
            #images.append(create_deviation(CCDData.read(fp, unit = u.adu),
             #                              gain = camera_electronic_gain,
              #                             readnoise = camera_readout_noise,
               #                            disregard_nan = True
                #                          ))
            images.append(CCDData.read(fp, unit = u.adu))
            names.append(fp)
            logging.info(f'image : {fp} loaded')
        
        return ImagesCombiner(images=images, names=names)
        #return cls(images)

    @classmethod
    def from_fits(cls, imgs: list[str],
                 camera_electronic_gain: float = 1.2 * u.electron / u.adu, 
                 camera_readout_noise: float =  2.2 * u.electron):
        
        images = []
        names: List[str] = []
        for fp in imgs:
            #images.append(create_deviation(CCDData.read(fp, unit = u.adu),
             #                              gain = camera_electronic_gain,
              #                             readnoise = camera_readout_noise,
               #                            disregard_nan = True
                #                          ))
            images.append(CCDData.read(fp, unit = u.adu))
            names.append(fp)
            logging.info(f'image : {fp} loaded')
        
        return ImagesCombiner(images=images, names=names)
        #return cls(images)

    #@classmethod
    def reduce_images_ccdproc(self) -> CCDData | None:
        conf: Config = Config()
        TRIM_REGION = None
        EXPOSURE_KEY = 'EXPTIME'
        CAPTURE_DIR =  str(Path(self.get_image_names()[0]).absolute().parent) + '/'
        
        # read master frames
        master_bias = None
        master_dark = None
        master_flat = None

        try:
            if (bias_file := conf.get_str('pre_processing', 'master_offset')) is not None:
                master_bias = CCDData.read(CAPTURE_DIR + bias_file, unit = u.adu)
                logging.info(f"masterbias loaded")
        except Exception as e:
            logging.error(f"cannot read masterbias: {e}")

        try:
            if (dark_file := conf.get_str('pre_processing', 'master_dark')) is not None:
                master_dark = CCDData.read(CAPTURE_DIR + dark_file, unit = u.adu)
                logging.info(f"masterdark loaded")
        except Exception as e:
            logging.error(f"cannot read masterdark: {e}")

        try:
            if (flat_file := conf.get_str('pre_processing', 'master_flat')) is not None:
                master_flat = CCDData.read(CAPTURE_DIR + flat_file, unit = u.adu)
                logging.info(f"masterflat loaded")
        except Exception as e:
            logging.error(f"cannot read masterflat: {e}")

        ### reduce science frames
        try:
            master_sciences = self.reduce(
                                        master_bias = master_bias, 
                                        master_dark = master_dark, 
                                        master_flat = master_flat, 
                                        exposure_key = EXPOSURE_KEY
                                        )        
        except Exception as e:
            logging.error(f"unable to reduce data: {e}")
            return None

        ### combine frames (sum or median) & save master science frame
        return master_sciences.sum()
    
    #@classmethod
    def reduce_images_numpy(self, img_data: List[np.ndarray] , preprocess: bool = False) -> np.ndarray:
        logging.info('summing images data ...')

        _reduced_img = img_data[0]
        for img in img_data:
            _reduced_img = np.add(_reduced_img, img)

        return _reduced_img
