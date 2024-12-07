import logging
import numpy as np
import math
from pathlib import Path

import tkinter as tk
from tkinter import ttk

import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import matplotlib.colors as mpl_colors
import matplotlib.colorbar as cb
from matplotlib.text import Annotation
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends._backend_tk import NavigationToolbar2Tk
from matplotlib.colors import LinearSegmentedColormap

from astropy.utils.exceptions import AstropyWarning
from astropy import units as u
from astropy.nddata import CCDData, NDData, StdDevUncertainty
from astropy.stats import mad_std
from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.table import QTable
from astropy.nddata import NDDataRef

from specutils.spectra.spectrum1d import Spectrum1D
from specutils.manipulation import median_smooth, gaussian_smooth
from specutils.analysis import snr, snr_derived
from specutils.manipulation import FluxConservingResampler, LinearInterpolatedResampler, SplineInterpolatedResampler

from specreduce.tracing import FlatTrace, FitTrace
from specreduce.background import Background
from specreduce.extract import BoxcarExtract
from specreduce.fluxcal import FluxCalibration
from specreduce import WavelengthCalibration1D

from app.os_utils import OSUtils
from app.config import Config
from img.image import Image
from img.img_utils import Images
from spc.spc_utils import rgb

class Spectrum(object):

    def __init__(self, spc_frame: ttk.Frame, img_axe: Axes) -> None: #, axe_img: Axes, axe_spc: Axes) -> None:
        self.conf = Config()

        # assign spectrum instance variables
        self.science_spectrum: Spectrum1D = None
        self.science_trace: FitTrace = None
        self.showed_lines: bool = False
        self.colors = ('blue', 'red', 'green', 'orange', 'cyan')
        self.lines_color = 'yellow'    

        # create figure and axe
        self.spc_figure = Figure(figsize=(10, 3))
        self.spc_axe = self.spc_figure.add_subplot(111)
        self.img_axe = img_axe
        self.spc_axe.grid(color = 'grey', linestyle = '--', linewidth = 0.5)

        # create spectrum canvas to draw to
        self.spc_canvas = FigureCanvasTkAgg(self.spc_figure, spc_frame)
        self.spc_canvas.draw()

        # create customized toolbar
        spc_toolbar = CustomSpcToolbar(self.spc_canvas, spc_frame)

        # add new buttons
        self.clear_button = ttk.Button(spc_toolbar, text="Clear", command=self.clear_spectra)
        self.clear_button.pack(side=tk.LEFT, padx=5, pady=0)

        self.lines_button = ttk.Button(spc_toolbar, text="Show lines", command=self.show_lines)
        self.lines_button.pack(side=tk.LEFT, padx=5, pady=0)

        self.lines_button = ttk.Button(spc_toolbar, text="Colorize", command=self.colorize)
        self.lines_button.pack(side=tk.LEFT, padx=5, pady=0)

        spc_toolbar.update()
        spc_toolbar.pack(side=tk.TOP, fill=tk.X)
        self.spc_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def open_spectrum( self, spc_name: str) -> bool:
        # open spectrum data
        try:
            spec1d: Spectrum1D = Spectrum1D.read(spc_name)
            self.show_spectrum(spec1d, True)

        except Exception as e:
            logging.error(f"{e}")
            return False
        
        self.science_spectrum = spec1d
        return True

    def colorize(self) -> None:
        if self.science_spectrum is None:
            logging.info("please calibrate first")
            return
        
        min_wl = 3800 # self.science_spectrum.wavelength.value.min()
        max_wl = 7500 #self.science_spectrum.wavelength.value.max()
        bin_size = 50

        colors = plt.cm.jet((self.science_spectrum.wavelength.value - min_wl) / (max_wl - min_wl))
        for i in range(0, len(self.science_spectrum.wavelength) - 1, bin_size):
            self.spc_axe.fill_between(self.science_spectrum.wavelength.value[i:i+bin_size], 0, 
                            self.science_spectrum.flux.value[i:i+bin_size], color=colors[i]) #, alpha=0.3)
            
        self.spc_figure.canvas.draw_idle()

    def clear_spectra(self) -> None:
        self.spc_axe.clear()
        self.spc_axe.grid(color = 'grey', linestyle = '--', linewidth = 0.5)
        self.showed_lines: bool = False
        self.spc_figure.canvas.draw_idle()

    def show_spectrum(self, spectrum: Spectrum1D, calibrated: bool = False) -> None:
        # set legend titles
        if calibrated:
            self.spc_axe.set_ylabel('Relative intensity')
            self.spc_axe.set_xlabel('Wavelength (Angstrom)')
            def format_coord(x,y):
                return "Lambda: ({:.2f}, Intensity: {:.2f})".format(x,y)
            self.spc_axe.format_coord=format_coord
        else:
            self.spc_axe.set_xlabel('Pixels')
            self.spc_axe.set_ylabel('ADU')
            def format_coord(x,y):
                return "Pixel: ({:.0f}, ADU: {:.0f})".format(x,y)
            self.spc_axe.format_coord=format_coord
        
        # plot spectrum
        self.colors = np.roll(self.colors, 1) # rotate a new color every call 
        _color = self.colors[0]
        self.spc_axe.plot(spectrum.spectral_axis , spectrum.flux, color=_color, linewidth = '0.8')
        self.spc_axe.grid(color = 'grey', linestyle = '--', linewidth = 0.5)

        self.spc_figure.canvas.draw_idle()

    def do_trace(self, img_stacked:np.ndarray) -> bool:
        if img_stacked is None: 
            logging.error("please reduce image(s) before tracing")
            return False

        master_science: np.ndarray = img_stacked
        
        #trace_model : one of Chebyshev1D, Legendre1D, Polynomial1D, or Spline1D
        #peak_method : One of gaussian, centroid, or max. gaussian
        #trace_model : one of ~astropy.modeling.polynomial.Chebyshev1D,
        #  ~astropy.modeling.polynomial.Legendre1D, 
        # ~astropy.modeling.polynomial.Polynomial1D,
        #  or ~astropy.modeling.spline.Spline1D, optional The 1-D polynomial model used to fit the trace to the bins' peak pixels. S
        # pline1D models are fit with Astropy's 'SplineSmoothingFitter', while the other models are fit with the 'LevMarLSQFitter'.
        #  [default: models.Polynomial1D(degree=1)]

        try:
            self.science_trace:FitTrace = FitTrace(master_science, 
                                    bins=self.conf.get_int('processing', 'trace_x_bins'), 
                                    #trace_model=models.Chebyshev1D(degree=2), 
                                    trace_model=models.Polynomial1D(degree=2),
                                    peak_method='gaussian' , #centroid',     #'gaussian', 
                                    window=self.conf.get_int('processing', 'trace_y_window'),
                                    guess=self.conf.get_float('processing', 'trace_y_guess')
                                    )
            
        except Exception as e:
            logging.error(f"unable to fit trace : {e}")
            return False

        # trace spectrum fitted
        self.science_trace.shape
        self.img_axe.plot(self.science_trace.trace , color='red', linestyle='dashed', linewidth = '0.5')
        self.img_axe.get_figure().canvas.draw_idle()

        logging.info(f'trace fitted : y = {self.science_trace.trace}')

        return True


    def do_extract(self, img_stacked:np.ndarray) -> bool:

        if (img_stacked is None) or (self.science_trace is None): 
            logging.error("please fit trace before extracting spectrum")
            return False

        master_science: np.ndarray = img_stacked
        
        try:
            bg_trace: Background = Background.two_sided(master_science, 
                                                  self.science_trace, 
                                                  separation=self.conf.get_int('processing', 'sky_y_offset'), 
                                                  width=self.conf.get_int('processing', 'sky_y_size')  ) 
        except Exception as e:
            logging.error(f"unable to fit background : {e}")
            return False
        
        logging.info('background extracted')
                
        try:
            extracted_spectrum = BoxcarExtract(master_science - bg_trace, 
                                    self.science_trace, 
                                    width = self.conf.get_float('processing', 'trace_y_size') )
        except Exception as e:
            logging.error(f"unable to extract background : {e}")
            return False

        logging.info('background substracted')

        try:
            science_spectrum = extracted_spectrum.spectrum
        except Exception as e:
            logging.error(f"unable to extract science spectrum : {e}")
            return False

        self.science_spectrum = science_spectrum
        logging.info('spectrum extracted')

        self.clear_spectra()

        # trace sky zones
        self.img_axe.plot(self.science_spectrum.spectral_axis, self.science_trace.trace + extracted_spectrum.width , color='blue', 
                            linestyle='dashed', linewidth = '0.5')  #, alpha=0.2)
        self.img_axe.plot(self.science_spectrum.spectral_axis, self.science_trace.trace - extracted_spectrum.width , color='blue', 
                            linestyle='dashed', linewidth = '0.5')  #, alpha=0.2)

        self.img_axe.plot(self.science_spectrum.spectral_axis, self.science_trace.trace + (self.conf.get_int('processing', 'sky_y_offset')) , 
                            color='green', linewidth = '0.5', linestyle='dashed')  #, alpha=0.2)
        self.img_axe.plot(self.science_spectrum.spectral_axis, 
                            self.science_trace.trace + (self.conf.get_int('processing', 'sky_y_offset') + self.conf.get_int('processing', 'sky_y_size')) , 
                            color='green', linewidth = '0.5', linestyle='dashed')  #, alpha=0.2)
        self.img_axe.plot(self.science_spectrum.spectral_axis, 
                            self.science_trace.trace - (self.conf.get_int('processing', 'sky_y_offset') + self.conf.get_int('processing', 'sky_y_size')) , 
                            color='green', linewidth = '0.5', linestyle='dashed')  #, alpha=0.2)
        self.img_axe.plot(self.science_spectrum.spectral_axis, 
                            self.science_trace.trace - (self.conf.get_int('processing', 'sky_y_offset')) , 
                            color='green', linewidth = '0.5', linestyle='dashed')  #, alpha=0.2)
        
        self.img_axe.get_figure().canvas.draw_idle()

        self.show_spectrum(self.science_spectrum, False)

        return True
    
    def do_calibrate(self) -> bool:
        if self.science_spectrum is None: 
            logging.error(f'please extract spectrum trace before calibrating')
            return False
        
        logging.info('calibrating spectrum...')
        
        _wavelength = self.conf.get_str('processing', 'calib_x_wavelength')     #[6506.53, 6532.88, 6598.95, 6678.28, 6717.04]*u.AA
        _pixels = self.conf.get_str('processing', 'calib_x_pixel')     #[770, 1190, 2240, 3484, 4160]*u.pix

        wavelength = [float(x) for x in _wavelength.replace(',', '').split()]*u.AA
        pixels = [float(x) for x in _pixels.replace(',', '').split()]*u.pix
        logging.info(f"calibrating pixels set : {pixels}")
        logging.info(f"with wavelengths set : {wavelength}")

        try:
            #input_spectrum, matched_line_list=None, line_pixels=None, line_wavelengths=None, catalog=None, input_model=Linear1D(), fitter=None
            #fitter: ~astropy.modeling.fitting.Fitter, optional The fitter to use in optimizing the model fit. Defaults to
            #~astropy.modeling.fitting.LinearLSQFitter if the model to fit is linear
            #or ~astropy.modeling.fitting.LMLSQFitter if the model to fit is non-linear.
            cal = WavelengthCalibration1D(input_spectrum = self.science_spectrum,
                line_wavelengths = wavelength,
                line_pixels = pixels,
                #matched_line_list = line_list,
                #input_model = models.Linear1D(),
                #fitter = fitting.LMLSQFitter(),
                input_model = models.Polynomial1D(degree = 2),
                #fitter = fitting.LMLSQFitter()
                #fitter = fitting.LinearLSQFitter()
                )
        except Exception as e:
            logging.error(f"unable to calibraet spectrum : {e}")
            return False

        logging.info(f"residuals : {repr(cal.residuals)}")
        
        # calibrate science spectrum
        calibrated_spectrum: Spectrum1D = cal.apply_to_spectrum(self.science_spectrum)
        logging.info('spectrum calibrated')

        # normalize to 1 arround 6500 - 6520 Ang
        #sci_mean_norm_region = calibrated_spectrum[6500 * u.AA: 6520 * u.AA].flux.mean()       # starEx2400 : high resolution
        sci_mean_norm_region = calibrated_spectrum[6500 * u.AA: 6520 * u.AA].flux.mean()       # starEx2400 : high resolution
        normalized_spec = Spectrum1D(spectral_axis = calibrated_spectrum.wavelength, flux = calibrated_spectrum.flux / sci_mean_norm_region)  

        logging.info('spectrum normalized')

        # aply response file if any defined
        #sci_spectrum = FluxCalibration(normalized_spec, 1.0) #, airmass = airmass) 
        final_spec = normalized_spec

        try:
            if (respFile := self.conf.get_str('processing', 'response_file')) is not None:
                respFile = f"{self.conf.get_conf_directory()}/{respFile}"
                logging.info(f"opening {respFile}...")
                resp1d: Spectrum1D = Spectrum1D.read(respFile)
                _factor = int(resp1d.shape[0] / normalized_spec.shape[0])
                logging.info(f"{normalized_spec.shape[0]=}, {resp1d.shape[0]=}, {_factor=}")
                
                _resp1d_ndd = NDDataRef(resp1d)
                _resp1d_ndd.wcs = None
                _resp1d = _resp1d_ndd[::_factor].data
                final_spec = normalized_spec / _resp1d
                logging.info('response applied')
            else:
                logging.warning("no response file to apply")

        except Exception as e:
            logging.error(f"{e}")
            logging.error("no response file applied")

        # apply median smoothing (if any defined)
        smooth: int | None = self.conf.get_int('post_processing', 'median_smooth') 
        if smooth is not None:
            smooth_spec: Spectrum1D = median_smooth(final_spec, width = smooth) 
            self.science_spectrum = smooth_spec
            logging.info(f"median smooth applied={smooth}")
        else:
            self.science_spectrum = final_spec

        logging.info(f'snr = {snr_derived(self.science_spectrum):.1f}')

        #self.spc_axe.clear()
        self.clear_spectra()
        self.show_spectrum(self.science_spectrum, True)

        logging.info('calibration complete')
        return True
    
    def show_lines(self, ax = None, show_line = True) -> None:
        if ax is None: ax = self.spc_axe
                        
        xbounds = ax.get_xbound()   
        trans = ax.get_xaxis_transform()

        if self.showed_lines is False:
            for wave, elm in self.conf.config.items('lines'):
                _lambda = (float(wave) * 10)   # convert nm to ang
                if (_lambda > xbounds[0]) & (_lambda < xbounds[1]):
                        #ax.axvline(lam, 0.95, 1.0, color = self.lines_color, lw = 0.5)
                        ax.axvline(_lambda, color = self.lines_color, lw = 0.5, linestyle = '--', alpha=0.8)
                        #trans = ax.get_xaxis_transform()
                        ax.annotate(elm, xy=(_lambda, 1.05), xycoords=trans, fontsize=8, rotation=90, color=self.lines_color)
                        
            # show colorband
            self.showed_lines = True
        else:
            # clear lines
            for line in ax.lines: 
                if line.get_color() == self.lines_color:
                    line.remove() 

            # clear elements
            for elm in ax.get_children():
                if isinstance(elm, Annotation): 
                    elm.remove()
            
            # redraw spectrum
            #self.show_spectrum(self.sci_spectrum, True)
            self.showed_lines = False
                    
        self.spc_figure.canvas.draw_idle()

class CustomSpcToolbar(NavigationToolbar2Tk):
    def __init__(self, canvas, parent) -> None:
        # list of toolitems to add/modify to the toolbar, format is:
        # (
        #   text, # the text of the button (often not visible to users)
        #   tooltip_text, # the tooltip shown on hover (where possible)
        #   image_file, # name of the image for the button (without the extension)
        #   name_of_method, # name of the method in NavigationToolbar2 to call
        # )
        self.toolitems = (
            ('Home', 'Reset zoom to original view', 'home_large', 'home'),
            ('Back', 'Back to previous view', 'back_large', 'back'),
            ('Forward', 'Forward to next view', 'forward_large', 'forward'),
            (None, None, None, None),
            ('Pan', 'Pan axes with left mouse, zoom with right', 'move_large', 'pan'),
            ('Zoom', 'Zoom to rectangle', 'zoom_to_rect_large', 'zoom'),
            (None, None, None, None),
            ('Save', 'Save the figure', 'filesave_large', 'save_figure'), 
        )

        super().__init__(canvas = canvas, window = parent, pack_toolbar = True)

                
