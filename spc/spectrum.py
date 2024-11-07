import logging
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from astropy.utils.exceptions import AstropyWarning
from astropy import units as u
from astropy.nddata import CCDData, NDData, StdDevUncertainty
from astropy.stats import mad_std
from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.table import QTable

from specutils.spectra.spectrum1d import Spectrum1D

from specreduce.tracing import FlatTrace, FitTrace
from specreduce.background import Background
from specreduce.extract import BoxcarExtract

from specreduce import WavelengthCalibration1D

from app.os_utils import OSUtils
from app.config import Config
from img.image import Image
from img.img_utils import show_image, reduce_images

class Spectrum(object):

    def __init__(self, axe_img: Axes, axe_spc: Axes) -> None:
        self.conf = Config()
        self.ax_img: Axes = axe_img
        self.ax_spc: Axes = axe_spc
        self.figure: Figure = axe_spc.get_figure()

    def do_calibration (self) -> bool:
        logging.info('extracting science spectra ...')

        self.ax_spc.clear()

        _reduced_img: CCDData | None = reduce_images(images=Image.img_names, preprocess=True)
        if Image.img_stacked is None: return False

        Image.img_stacked = _reduced_img.data
        
       # collect image stats
        vstd = Image.img_stacked.std()
        vmean = Image.img_stacked.mean()
        _min = Image.img_stacked.min()
        _max = Image.img_stacked.max()
        vmin = vmean - vstd
        vmax = vmean + vstd

        # display image
        logging.info (f"image stats : min = {_min}, max = {_max}, mean = {vmean}, std = {vstd}")
        show_image(image = Image.img_stacked,
                        fig = self.ax_img.get_figure(),
                        ax = self.ax_img,
                        show_colorbar = False, #show_colorbar, 
                        cmap = self.conf.get_str('window', 'colormap'))

        self.figure.canvas.draw_idle()

        master_science: np.ndarray = Image.img_stacked
        
        #trace_model : one of Chebyshev1D, Legendre1D, Polynomial1D, or Spline1D
        #peak_method : One of gaussian, centroid, or max. gaussian
        #trace_model : one of ~astropy.modeling.polynomial.Chebyshev1D,
        #  ~astropy.modeling.polynomial.Legendre1D, 
        # ~astropy.modeling.polynomial.Polynomial1D,
        #  or ~astropy.modeling.spline.Spline1D, optional The 1-D polynomial model used to fit the trace to the bins' peak pixels. S
        # pline1D models are fit with Astropy's 'SplineSmoothingFitter', while the other models are fit with the 'LevMarLSQFitter'.
        #  [default: models.Polynomial1D(degree=1)]

        _guess: float | None = self.conf.get_float('processing', 'trace_y_guess')
        try:
            sci_tr:FitTrace = FitTrace(master_science, 
                                    bins = self.conf.get_int('processing', 'trace_x_bins'), 
                                    trace_model=models.Chebyshev1D(degree=2), 
                                    peak_method = 'centroid',     #'gaussian', 
                                    window = self.conf.get_int('processing', 'trace_y_window'),
                                    guess =  _guess)
            
        except Exception as e:
            logging.error(f"unable to fit trace : {e}")
            return False

        logging.info(f'trace fitted : y = {sci_tr.trace}')

        _bg_separation = self.conf.get_int('processing', 'sky_y_offset')    #80
        _bg_width = self.conf.get_int('processing', 'sky_y_size')    #50
        _trace_width = self.conf.get_int('processing', 'trace_y_size')   #15
        try:
            bg: Background = Background.two_sided(master_science, sci_tr, separation=_bg_separation, width=_bg_width) 
        except Exception as e:
            logging.error(f"unable to fit background : {e}")
            return False
        
        logging.info('background extracted')
        try:
            extract = BoxcarExtract(master_science - bg, sci_tr, width = _trace_width)
        except Exception as e:
            logging.error(f"unable to extract background : {e}")
            return False

        logging.info('background substracted')

        try:
            sci_spectrum: Spectrum1D = extract.spectrum
        except Exception as e:
            logging.error(f"unable to extract science spectrum : {e}")
            return False

        logging.info('spectrum extracted')

        if self.conf.get_bool('processing', 'show_trace'):
            #self.ax_img.imshow(bg.bkg_wimage, origin='lower', aspect='auto', cmap=plt.cm.gray, alpha=0.2)
            #self.ax_img.imshow(sci_tr.image.data, origin='lower', aspect='auto', cmap=plt.cm.gray, alpha=0.2)
            self.ax_img.step(sci_spectrum.spectral_axis, sci_tr.trace , color='red', linestyle='dashed', linewidth = '0.3')
            self.ax_img.step(sci_spectrum.spectral_axis, sci_tr.trace + extract.width , color='green', linestyle='dashed', linewidth = '0.5')  #, alpha=0.2)
            self.ax_img.step(sci_spectrum.spectral_axis, sci_tr.trace - extract.width , color='green', linestyle='dashed', linewidth = '0.5')  #, alpha=0.2)

            self.ax_img.step(sci_spectrum.spectral_axis, sci_tr.trace + (_bg_separation) , color='blue', linewidth = '0.5', linestyle='dashed')  #, alpha=0.2)
            self.ax_img.step(sci_spectrum.spectral_axis, sci_tr.trace + (_bg_separation + _bg_width) , color='blue', linewidth = '0.5', linestyle='dashed')  #, alpha=0.2)

            self.ax_img.step(sci_spectrum.spectral_axis, sci_tr.trace - (_bg_separation + _bg_width) , color='blue', linewidth = '0.5', linestyle='dashed')  #, alpha=0.2)
            self.ax_img.step(sci_spectrum.spectral_axis, sci_tr.trace - (_bg_separation) , color='blue', linewidth = '0.5', linestyle='dashed')  #, alpha=0.2)
            
            #self.ax_spc.set_title('spectrum2D + background + trace fitted')

            #self.ax_spc.step(sci_spectrum.spectral_axis , sci_spectrum.flux, color='red', linewidth = '0.3')
            #self.ax_spc.set_xlabel('Pixels')
            #self.ax_spc.set_ylabel('ADU')

        logging.info('calibrating neon spectrum...')
        
        _wavelength = self.conf.get_str('processing', 'calib_x_wavelength')     #[6506.53, 6532.88, 6598.95, 6678.28, 6717.04]*u.AA
        _pixels = self.conf.get_str('processing', 'calib_x_pixel')     #[770, 1190, 2240, 3484, 4160]*u.pix

        #__wavelength = [6506.53, 6532.88, 6598.95, 6678.28, 6717.04]*u.AA
        #__pixels = [770, 1190, 2240, 3484, 4160]*u.pix

        wavelength = [float(x) for x in _wavelength.replace(',', '').split()]*u.AA
        pixels = [float(x) for x in _pixels.replace(',', '').split()]*u.pix
        #print(wavelength, pixels)
        logging.info(f"calibrating pixels set : {pixels}")
        logging.info(f"with wavelengths set : {wavelength}")

        #input_spectrum, matched_line_list=None, line_pixels=None, line_wavelengths=None, catalog=None, input_model=Linear1D(), fitter=None
        #fitter: ~astropy.modeling.fitting.Fitter, optional The fitter to use in optimizing the model fit. Defaults to
        #~astropy.modeling.fitting.LinearLSQFitter if the model to fit is linear
        #or ~astropy.modeling.fitting.LMLSQFitter if the model to fit is non-linear.
        cal = WavelengthCalibration1D(input_spectrum = sci_spectrum,
            line_wavelengths = wavelength,
            line_pixels = pixels,
            #matched_line_list = line_list,
            #input_model = models.Polynomial1D(degree=2), # .Linear1D(),
            #input_model = models.Polynomial1D(degree = 2),
            #fitter = fitting.LMLSQFitter()
            #fitter = fitting.LinearLSQFitter()
            )
        
        logging.info(f"residuals : {repr(cal.residuals)}")
        
        logging.info('neon spectrum calibrated')

        calibrated_spectrum: Spectrum1D = cal.apply_to_spectrum(sci_spectrum)
        logging.info('science spectrum calibrated')

        sci_mean_norm_region = calibrated_spectrum[6500 * u.AA: 6520 * u.AA].flux.mean()       # starEx2400 : high resolution
        final_spec = Spectrum1D(spectral_axis = calibrated_spectrum.wavelength, flux = calibrated_spectrum.flux / sci_mean_norm_region)  

        self.ax_spc.set_ylabel('Relative intensity')
        self.ax_spc.set_xlabel('Wavelength (Angstrom)')
        self.ax_spc.grid(color = 'grey', linestyle = '--', linewidth = 0.5)

        self.ax_spc.plot(final_spec.spectral_axis, final_spec.flux, color='red', linewidth = '0.4')
        self.figure.canvas.draw_idle()

        logging.info('calibration complete')


        plt.figure(figsize = (10,6))
        plt.step(final_spec.wavelength, final_spec.flux, color='black', linewidth = '0.6') #, where="mid")
        plt.xlabel('Wavelength (Ang)')
        plt.ylabel('ADU')
        #plt.ylim(-10000, 1e6)

        Spectrum.show_lines(ax = self.ax_spc, show_line = True)
        return True
    
    @staticmethod
    def show_lines(ax = None, show_line = True):
        """
        Show lines onto a plot. 
            
        Parameters
        ----------    
        ax : AxesSubplot
            The axis onto which the emission/absoption lines needs to be plotted.
            If ax = None, then the plotting function uses plt, rather than axis.
            
        show_line : bool
            Whether or not to draw vertical dashed lines. Default is True.
        
        Returns
        -------
        None
        
        """
        
        lines_to_display = [
            {"name" : 'Zero Order',"label" :'Zero', "lambda" :0.00}, 
            {"name" : 'Hydrogen', "label" : 'Hα',  "lambda" :656.2852}, 
            {"name" : 'Hydrogen', "label" : 'Hβ',  "lambda" :486.133 },
            {"name" : 'Hydrogen', "label" : 'Hγ',  "lambda" :434.047}, 
            {"name" : 'Hydrogen', "label" : 'Hδ',  "lambda" :410.174}, 
            {"name" : 'Hydrogen', "label" : 'Hε',  "lambda" :397.007}, 
            {"name" : 'Hydrogen', "label" : 'Hζ',  "lambda" :388.9049}, 
            {"name" : 'Hydrogen', "label" : 'Hη',  "lambda" :383.5384}, 
            {"name" : 'Hydrogen', "label" : 'Hθ',  "lambda" :379.75} , 
            {"name" : 'Iron', "label" : 'Fe',  "lambda" :527.04}, 
            {"name" : 'Iron', "label" : 'Fe',  "lambda" :516.89}, 
            {"name" : 'Iron', "label" : 'Fe',  "lambda" :495.76}, 
            {"name" : 'Iron', "label" : 'Fe',  "lambda" :466.81}, 
            {"name" : 'Iron', "label" : 'Fe',  "lambda" :438.36}, 
            {"name" : 'Iron', "label" : 'Fe',  "lambda" :430.79}, 
            {"name" : 'Magnesium', "label" : 'MgII',  "lambda" :448.11}, 
            {"name" : 'Magnesium', "label" : 'Mg',  "lambda" :518.36}, 
            {"name" : 'Magnesium', "label" : 'Mg',  "lambda" :517.27}, 
            {"name" : 'Magnesium', "label" : 'Mg',  "lambda" :516.73}, 
            {"name" : 'Neon', "label" : 'NeI',  "lambda" :585.249}, 
            {"name" : 'Neon', "label" : 'NeI',  "lambda" :588.189}, 
            #{"name" : 'Mercury', "label" : 'Hg',  "lambda" :404.656}, 
            #{"name" : 'Mercury', "label" : 'Hg',  "lambda" :435.833}, 
            #{"name" : 'Mercury', "label" : 'Hg',  "lambda" :546.074}, 
            #{"name" : 'Mercury', "label" : 'Hg',  "lambda" :576.960}, 
            #{"name" : 'Mercury', "label" : 'Hg',  "lambda" :578.966}, 
            {"name" : 'Sodium', "label" : 'NaI',  "lambda" :589.00}, 
            {"name" : 'Sodium', "label" : 'NaI',  "lambda" :589.59}, 
            {"name" : 'Oxygen', "label" : 'O1',  "lambda" :615.82}, 
            {"name" : 'Oxygen', "label" : 'O2',  "lambda" :627.77}, 
            {"name" : 'Oxygen', "label" : 'O2',  "lambda" :686.9}, 
            {"name" : 'Oxygen', "label" : 'O2',  "lambda" :718.6}, 
            {"name" : 'Oxygen', "label" : 'O2',  "lambda" :760.5}, 
            {"name" : 'Oxygen', "label" : 'O2',  "lambda" :898.77}, 
            {"name" : 'Oxygen', "label" : 'OIII',  "lambda" :495.9}, 
            {"name" : 'Oxygen', "label" : 'OIII',  "lambda" :500.69}, 
            {"name" : 'Water', "label" : 'H2O',  "lambda" :651.65}, 
            {"name" : 'Water', "label" : 'H2O',  "lambda" :694.07}, 
            {"name" : 'Water', "label" : 'H2O',  "lambda" :695.64}, 
            {"name" : 'Water', "label" : 'H2O',  "lambda" :698.90}, 
            {"name" : 'Calcium', "label" : 'Ca+ H',  "lambda" :396.85}, 
            {"name" : 'Calcium', "label" : 'Ca+ K',  "lambda" :393.37}, 
            {"name" : 'Helium', "label" : 'He I',  "lambda" :706.52}, 
            {"name" : 'Helium', "label" : 'He I',  "lambda" :667.82}, 
            {"name" : 'Helium', "label" : 'He I',  "lambda" :587.56}, 
            {"name" : 'Helium', "label" : 'He I',  "lambda" :501.57}, 
            {"name" : 'Helium', "label" : 'He I',  "lambda" :447.148}, 
            {"name" : 'Silicon', "label" : 'Si II',  "lambda" :634.71}, 
            {"name" : 'Silicon', "label" : 'Si II',  "lambda" :637.14}, 
            {"name" : 'Terbium', "label" : 'Tb',  "lambda" :487.7}, 
            {"name" : 'Terbium', "label" : 'Tb',  "lambda" :542.4}, 
            {"name" : 'Europium', "label" : 'Eu',  "lambda" :611.6}, 
            {"name" : 'Titanium', "label" : 'Ti+',  "lambda" :336.11}
        ]

        if (ax == None):
            ax = plt.gca()
            
        xbounds = ax.get_xbound()   # Getting the x-range of the plot     
        for ii in range(len(lines_to_display)):
            lam = lines_to_display[ii]['lambda'] * 10    # nm to AA
            if (lam > xbounds[0]) & (lam < xbounds[1]):
                ax.axvline(lam, 0.95, 1.0, color = 'yellow', lw = 1.0)
                ax.axvline(lam, color = 'yellow', lw = 0.8, linestyle = '--')
                trans = ax.get_xaxis_transform()
                ax.annotate(lines_to_display[ii]['label'], xy = (lam, 1.05), xycoords = trans, \
                        fontsize = 8, rotation = 90, color = 'yellow')