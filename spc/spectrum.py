import logging
import numpy as np
#import random
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import matplotlib.colors as mpl_colors
import matplotlib.colorbar as cb
from matplotlib.text import Annotation
import matplotlib.gridspec as gridspec

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

    def __init__(self, axe_img: Axes, axe_spc: Axes) -> None:
        self.conf = Config()
        self.ax_img: Axes = axe_img
        self.ax_spc: Axes = axe_spc
        self.figure: Figure = axe_spc.get_figure()
        self.sci_spectrum: Spectrum1D = None
        self.showed_lines: bool = False
        self.colors = ('blue', 'red', 'green', 'orange', 'cyan')
        self.lines_color = 'yellow'    
        self.ax_spc.grid(color = 'grey', linestyle = '--', linewidth = 0.5)

    def open_spectrum( self, spc_name: str) -> None:
        # open spectrum data
        try:
            spec1d: Spectrum1D = Spectrum1D.read(spc_name)
            self.show_spectrum(spec1d, True)

        except Exception as e:
            logging.error(f"{e}")
            return
        
        self.sci_spectrum = spec1d

    def clear_spectrum(self) -> None:
        #self._figure.clear()
        #self._figure.clf()
        self.ax_spc.clear()
        self.figure.canvas.draw()

    def show_spectrum(self, spectrum: Spectrum1D, calibrated: bool = False) -> None:
        # show axes legend
        if calibrated:
            self.ax_spc.set_ylabel('Relative intensity')
            self.ax_spc.set_xlabel('Wavelength (Angstrom)')
        else:
            self.ax_spc.set_xlabel('Pixels')
            self.ax_spc.set_ylabel('ADU')
        
        # plot spectrum
        self.colors = np.roll(self.colors, 1) # pick a new color every call 
        _color = self.colors[0]
        self.ax_spc.plot(spectrum.spectral_axis , spectrum.flux, color=_color, linewidth = '0.8')
        self.ax_spc.grid(color = 'grey', linestyle = '--', linewidth = 0.5)

        self.figure.canvas.draw_idle()

    def do_extract(self, img_stacked:np.ndarray) -> bool:
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
            sci_trace:FitTrace = FitTrace(master_science, 
                                    bins=self.conf.get_int('processing', 'trace_x_bins'), 
                                    #trace_model=models.Chebyshev1D(degree=2), 
                                    trace_model=models.Polynomial1D(degree=2),
                                    peak_method='centroid',     #'gaussian', 
                                    window=self.conf.get_int('processing', 'trace_y_window'),
                                    guess=self.conf.get_float('processing', 'trace_y_guess')
                                    )
            
        except Exception as e:
            logging.error(f"unable to fit trace : {e}")
            return False

        logging.info(f'trace fitted : y = {sci_trace.trace}')
        
        try:
            bg: Background = Background.two_sided(master_science, 
                                                  sci_trace, 
                                                  separation=self.conf.get_int('processing', 'sky_y_offset'), 
                                                  width=self.conf.get_int('processing', 'sky_y_size')  ) 
        except Exception as e:
            logging.error(f"unable to fit background : {e}")
            return False
        
        logging.info('background extracted')
                
        try:
            extract = BoxcarExtract(master_science - bg, 
                                    sci_trace, 
                                    width = self.conf.get_float('processing', 'trace_y_size') )
        except Exception as e:
            logging.error(f"unable to extract background : {e}")
            return False

        logging.info('background substracted')

        try:
            sci_spectrum = extract.spectrum
        except Exception as e:
            logging.error(f"unable to extract science spectrum : {e}")
            return False

        self.sci_spectrum = sci_spectrum
        logging.info('spectrum extracted')

        self.ax_spc.clear()

        # trace spectrum zone
        self.ax_img.plot(self.sci_spectrum.spectral_axis, sci_trace.trace , color='red', 
                            linestyle='dashed', linewidth = '0.5')
        self.ax_img.plot(self.sci_spectrum.spectral_axis, sci_trace.trace + extract.width , color='green', 
                            linestyle='dashed', linewidth = '0.5')  #, alpha=0.2)
        self.ax_img.plot(self.sci_spectrum.spectral_axis, sci_trace.trace - extract.width , color='green', 
                            linestyle='dashed', linewidth = '0.5')  #, alpha=0.2)
        
        # trace sky zones
        self.ax_img.plot(self.sci_spectrum.spectral_axis, sci_trace.trace + (self.conf.get_int('processing', 'sky_y_offset')) , 
                            color='green', linewidth = '0.5', linestyle='dashed')  #, alpha=0.2)
        self.ax_img.plot(self.sci_spectrum.spectral_axis, 
                            sci_trace.trace + (self.conf.get_int('processing', 'sky_y_offset') + self.conf.get_int('processing', 'sky_y_size')) , 
                            color='green', linewidth = '0.5', linestyle='dashed')  #, alpha=0.2)
        self.ax_img.plot(self.sci_spectrum.spectral_axis, 
                            sci_trace.trace - (self.conf.get_int('processing', 'sky_y_offset') + self.conf.get_int('processing', 'sky_y_size')) , 
                            color='green', linewidth = '0.5', linestyle='dashed')  #, alpha=0.2)
        self.ax_img.plot(self.sci_spectrum.spectral_axis, 
                            sci_trace.trace - (self.conf.get_int('processing', 'sky_y_offset')) , 
                            color='green', linewidth = '0.5', linestyle='dashed')  #, alpha=0.2)
        
        self.show_spectrum(self.sci_spectrum, False)

        return True
    
    def do_calibrate(self) -> bool:
        if self.sci_spectrum is None: 
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
            cal = WavelengthCalibration1D(input_spectrum = self.sci_spectrum,
                line_wavelengths = wavelength,
                line_pixels = pixels,
                #matched_line_list = line_list,
                #input_model = models.Linear1D(),
                input_model = models.Polynomial1D(degree = 2),
                fitter = fitting.LMLSQFitter()
                #fitter = fitting.LinearLSQFitter()
                )
        except Exception as e:
            logging.error(f"unable to calibraet spectrum : {e}")
            return False

        logging.info(f"residuals : {repr(cal.residuals)}")
        
        # calibrate science spectrum
        calibrated_spectrum: Spectrum1D = cal.apply_to_spectrum(self.sci_spectrum)
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
            self.sci_spectrum = smooth_spec
            logging.info(f"median smooth applied={smooth}")
        else:
            self.sci_spectrum = final_spec

        logging.info(f'snr = {snr_derived(self.sci_spectrum):.1f}')

        self.ax_spc.clear()

        self.show_spectrum(self.sci_spectrum, True)

        logging.info('calibration complete')
        return True
    
    
    def do_clear(self) -> None:
        self.ax_spc.clear()
        #self.ax_img.clear()
        self.figure.canvas.draw_idle()

    def show_lines(self, ax = None, show_line = True):
        if ax is None: ax = self.ax_spc
                        
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
            #""""
            _band_size = int(self.sci_spectrum.wavelength[1].value - self.sci_spectrum.wavelength[0].value)
            if _band_size <= 1:
                _band_size = 1
            _min_flux = np.min(self.sci_spectrum.flux.value)
            _max_flux = np.max(self.sci_spectrum.flux.value)
            for i, _wave in zip(range(0, len(self.sci_spectrum.wavelength)), self.sci_spectrum.spectral_axis):
                if (i+1) % _band_size == 0:
                    _lambda = _wave.value
                    _flux = ((self.sci_spectrum.flux.value[i]) / (_max_flux - _min_flux))   # convert to 0..1 range
                    if (_flux < 0) or (math.isnan(_flux)):
                        _flux = 0
                    elif (_flux > 1):
                        _flux = 1
                    #self.ax_spc.axvline(_lambda, 0.95, 1.0, color = rgb(_lambda / 10.0), lw = _band_size, alpha=_flux)
                    if (_lambda > xbounds[0]) & (_lambda < xbounds[1]):
                        ax.annotate('||', xy=(_lambda - _band_size, 0.90), xycoords=trans, 
                                    fontsize=14, weight='bold', color=rgb(_lambda / 10.0), alpha=_flux)

            #"""
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
                    
        self.figure.canvas.draw_idle()

                
