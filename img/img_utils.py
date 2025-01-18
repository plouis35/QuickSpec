""" 
utility routines to operate on a set of FIT images thru single line python statements

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
from astropy.utils.exceptions import AstropyWarning
from ccdproc import combine, subtract_bias, subtract_dark, flat_correct
from ccdproc import trim_image, Combiner, ccd_process, cosmicray_median, cosmicray_lacosmic, create_deviation
from ccdproc import ImageFileCollection, gain_correct
from astropy.stats import mad_std

from app.config import Config

warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', UserWarning)

class ImagesCombiner(object):
    def __init__(self, images: List[CCDData], names: List[str], max_memory: float = 2e9) -> None:
        """
        maintains images set array and file names
        max memory is used by ccdproc routines to avoid OOM exceptions when working with large set of big images
 
        Args:
            images (List[CCDData]): list of images data
            names (List[str]): list of images names
            max_memory (float, optional): max memory (in Bytes) to use when using CCDProc funcs. Defaults to 1e9.
        """        
        self._images: List[CCDData] = images
        self._image_names: List[str] = names
        self._memory_limit: float = max_memory

    def __getitem__(self, i:int) -> CCDData:
        """
        returns a specific image array thru its index

        Args:
            i (int): index of image to return

        Returns:
            CCDData: image data
        """        
        return self._images[i]

    def __len__(self) -> int:
        """
        returns the number of images loaded in this set

        Returns:
            int: nb of images
        """        
        return len(self._images)
    
    def get_image_names(self) -> List[str]:
        """
        returns all image names

        Returns:
            List[str]: list of names
        """        
        return self._image_names

    def sum(self) -> CCDData:
        """
        returns the sum of images loaded in this set

        Returns:
            CCDData: image sum'ed 
        """        
        logging.info(f'sum combine on {len(self._images)} images ...')
        return(combine(self._images,
                       method = 'sum',
                       dtype = np.float32,
                       mem_limit = self._memory_limit)
              )

    def sigmaclip(self, low_thresh: int = 5, high_thresh: int = 5) -> CCDData:
        """
        returns the median'ed - sigmaclip'ed of images loaded in this set

        Args:
            low_thresh (int, optional): any pixel greater than that will be rejected. Defaults to 5.
            high_thresh (int, optional): any pixel lower than that will be rejected. Defaults to 5.

        Returns:
            CCDData: image clipped
        """        
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

    def median(self) -> CCDData:
        """
        returns the median of images loaded in this set

        Returns:
            CCDData: image median'ed
        """        
        logging.info(f'median combine on {len(self._images)} images ...')
        return (combine(self._images, 
                        method = 'median', 
                        dtype = np.float32, 
                        mem_limit = self._memory_limit)
               )

    def trim(self, trim_region: str | None):
        """
        trim all frames loaded in this set

        Args:
            trim_region (str | None): x1, y1, x2, y2 rectangle to trim to

        Returns:
            self: images set updated
        """        
        if trim_region is not None:
            for i in range(0, len(self._images)):
                self._images[i] = trim_image(self._images[i][eval(trim_region)[1]:eval(trim_region)[3],
                                                              eval(trim_region)[0]:eval(trim_region)[2]])

            logging.info(f'{len(self._images)} images trimmed to ({trim_region})')
        else:
            logging.info('no trimming')
        return self

    def y_crop(self, y_ratio: float | None):
        """
        trim all frames loaded in this set on an y-axis percentage

        Args:
            y_ratio (float | None): percentage (0.1 .. 1.0) ratio from y-center

        Returns:
            self: images set updated
        """        
        if y_ratio is not None:
            x1 = 0
            y1 = round(self._images[0].shape[0] * y_ratio)
            x2 = self._images[0].shape[1]
            y2 = round(self._images[0].shape[0] - (self._images[0].shape[0] * y_ratio))

            for i in range(0, len(self._images)):
                self._images[i] = trim_image(self._images[i][y1:y2, x1:x2])

            logging.info(f"{len(self._images)} images y-cropped to {y_ratio}x: {x1=},{x2=},{y1=},{y2=}")
        else:
            logging.info('no y-cropping to do')
        return self


    def offset(self, scalar):
        """
        add a scalar to all frames loaded in this set

        Args:
            scalar (int): number to add to all images pixel values

        Returns:
            self: images set updated
        """        
        for i in range(0, len(self._images)):
            self._images[i] = CCDData(CCDData.add(self._images[i], scalar), 
                                      header = self._images[i].header) #, unit = self._images[i].unit
            
        logging.info(f'{len(self._images)} images added by ({scalar})')
        return self
    
    def bias_substract(self, frame):
        """
        substract a master bias frame to all frames loaded in this set

        Args:
            frame (CCDdata): bias frame to use

        Returns:
            self: images set updated
        """        
        for i in range(0, len(self._images)):
            self._images[i] = subtract_bias(self._images[i], frame)
            
        logging.info(f'masterbias substracted to {len(self._images)} images')
        return self

    def dark_substract(self, frame, scale_exposure: bool = True, exposure = 'EXPTIME'):
        """
        substract a master dark frame to all frames loaded in this set

        Args:
            frame (CCDdata): dark frame to use
            scale_exposure (bool):  If True, scale the dark frame by the exposure times. Default is True  

        Returns:
            self: images set updated
        """        
        for i in range(0, len(self._images)):
            self._images[i] = subtract_dark(self._images[i], frame, scale = scale_exposure, exposure_time = exposure, exposure_unit = u.second)                
        
        logging.info(f'masterdark substracted to {len(self._images)} images')
        return self
    
    def flat_divide(self, frame):
        """
        divide a master flat frame to all frames loaded in this set

        Args:
            frame (CCDdata): flat frame to use

        Returns:
            self: images set updated
        """        
        for i in range(0, len(self._images)):
            self._images[i] = flat_correct(ccd = self._images[i], flat = frame, min_value = None) #, norm_value = 10000 * u.adu)

        logging.info(f'masterflat divided to {len(self._images)} images')
        return self

    def reduce(self, master_bias, master_dark, master_flat, exposure_key = 'EXPTIME'):
        """
        process science frames
        masterdark frame is scaled according to science frame exposure duration

        Args:
            master_bias (CCDData): bias to use
            master_dark (CCDData): dark to use
            master_flat (CCDData): flat to use
            exposure_key (str, optional): FIT keyword to use for exposure. Defaults to 'EXPTIME'.
            
        Returns:
            self: images set updated
        """        
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

    def reduce_images(self) -> CCDData | None:
        """
        wrapper to CCDProc reduce routine - operates on all images set

        Returns:
            CCDData | None: sum of images reduced
        """        
        conf: Config = Config()
        TRIM_REGION = None
        EXPOSURE_KEY = 'EXPTIME'
        CAPTURE_DIR =  str(Path(self.get_image_names()[0]).absolute().parent) + '/'
        
        # read master frames
        master_bias = None
        master_dark = None
        master_flat = None

        y_ratio: float | None = conf.get_float('pre_processing','crop_auto')

        try:
            if (bias_file := conf.get_str('pre_processing', 'master_offset')) is not None:
                master_bias = CCDData.read(CAPTURE_DIR + bias_file, unit = u.adu)
                if y_ratio is not None:
                    x1 = 0
                    y1 = round(master_bias.shape[0] * y_ratio)
                    x2 = master_bias.shape[1]
                    y2 = round(master_bias.shape[0] - (master_bias.shape[0] * y_ratio))
                    master_bias = trim_image(master_bias[y1:y2, x1:x2])

                logging.info(f"masterbias loaded")
        except Exception as e:
            logging.error(f"cannot read masterbias: {e}")

        try:
            if (dark_file := conf.get_str('pre_processing', 'master_dark')) is not None:
                master_dark = CCDData.read(CAPTURE_DIR + dark_file, unit = u.adu)
                if y_ratio is not None:
                    x1 = 0
                    y1 = round(master_dark.shape[0] * y_ratio)
                    x2 = master_dark.shape[1]
                    y2 = round(master_dark.shape[0] - (master_dark.shape[0] * y_ratio))
                    master_dark = trim_image(master_dark[y1:y2, x1:x2])

                logging.info(f"masterdark loaded")
        except Exception as e:
            logging.error(f"cannot read masterdark: {e}")

        try:
            if (flat_file := conf.get_str('pre_processing', 'master_flat')) is not None:
                master_flat = CCDData.read(CAPTURE_DIR + flat_file, unit = u.adu)
                if y_ratio is not None:
                    x1 = 0
                    y1 = round(master_flat.shape[0] * y_ratio)
                    x2 = master_flat.shape[1]
                    y2 = round(master_flat.shape[0] - (master_flat.shape[0] * y_ratio))
                    master_flat = trim_image(master_flat[y1:y2, x1:x2])

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

        ### combine reduced frames
        #return master_sciences.median()
        return master_sciences.sum()

class Images(ImagesCombiner):
    """
    Images class implements the image file loader methods

    Args:
        ImagesCombiner (List[CCDData]): images set
    """    
    def __init__(self, images: List[CCDData]) -> None:
        ImagesCombiner.__init__(self, images, [])

    @classmethod
    def find_files(cls, directory: str, files_filter: str, sort_key: str = 'date-obs') -> List[str]:
        """
        collect and sort (according to fit header 'date-obs') file names using a wildcard filter

        Args:
            directory (str): path
            files_filter (str): wildcard filter
            sort_key (str, optional): FIT keyword to use for sorting. Defaults to 'date-obs'.

        Returns:
            List[str]: list if image names
        """        
        ic = ImageFileCollection(directory, glob_include = files_filter)
        ic.sort([sort_key])
        return (ic.files_filtered(include_path=True))

    @classmethod
    def from_fit(cls, dir: str, filter: str, 
                 camera_electronic_gain: float = 1.2 * u.electron / u.adu, 
                 camera_readout_noise: float =  2.2 * u.electron):
        """_summary_
        load a set of FIT images from a directory

        Args:
            dir (str): path
            filter (str): wildcard filter
            camera_electronic_gain (float, optional): . Defaults to 1.2*u.electron/u.adu.
            camera_readout_noise (float, optional): . Defaults to 2.2*u.electron.

        Returns:
            List (CCData): images set loaded
        """        
        
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

    @classmethod
    def from_fits(cls, imgs: list[str],
                 camera_electronic_gain: float = 1.2 * u.electron / u.adu, 
                 camera_readout_noise: float =  2.2 * u.electron):
        """
        load a set of FIT images from a list of names

        Args:
            imgs (list[str]): list of names
            camera_electronic_gain (float, optional): . Defaults to 1.2*u.electron/u.adu.
            camera_readout_noise (float, optional): . Defaults to 2.2*u.electron.

        Returns:
            _type_: _description_
        """        
        
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
