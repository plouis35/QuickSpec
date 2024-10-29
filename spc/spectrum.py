import logging
import warnings
import psutil
import contextlib

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import matplotlib.pyplot as plt
from matplotlib.image import AxesImage
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from astropy.utils.exceptions import AstropyWarning
from astropy import units as u
from astropy.nddata import CCDData, NDData, StdDevUncertainty
from astropy.stats import mad_std
from astropy.io import fits
from astropy.modeling import models, fitting

from specutils.spectra.spectrum1d import Spectrum1D

from specreduce.tracing import FlatTrace, FitTrace
from specreduce.background import Background
from specreduce.extract import BoxcarExtract

from specreduce import WavelengthCalibration1D


from app.config import Config

warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', UserWarning)

class spec1d:

    def __init__(self, ax_spc: Axes, spec2d) -> None:
        self.conf = Config()
        spec1d.spec2d_cls = spec2d
        self._spec1d_axe: Axes = ax_spc
        self._spec2d_axe: Axes = spec2d.img_axe
        self._figure: Figure | None = ax_spc.get_figure()

        # create calib button
        self._run_calib_ax = self._figure.add_axes(rect=(0.935, 0.35, 0.06, 0.05)) # (left, bottom, width, height
        self._bt_run = Button(self._run_calib_ax,'Calib', color='k', hovercolor='grey')
        self._bt_run.on_clicked(self.run_calib)

    def run_calib(self, event) -> None:
        logging.info('extracting science spectra ...')

        master_science: NDData =  spec1d.spec2d_cls.img_stacked 

        sci_tr:FitTrace = FitTrace(master_science, bins = 64, 
                                  trace_model=models.Polynomial1D(degree=2), 
                                  peak_method = 'gaussian', 
                                  window = 50) #, guess=605) #, guess=407)

        logging.info('trace fitted')
        #trace_model : one of Chebyshev1D, Legendre1D, Polynomial1D, or Spline1D
        #peak_method : One of gaussian, centroid, or max. gaussian
        _bg_separation = 80
        _bg_width = 50
        _trace_width = 15
        bg: Background = Background.two_sided(master_science, sci_tr, separation=_bg_separation, width=_bg_width) 
        logging.info('background extracted')
        extract = BoxcarExtract(master_science - bg, sci_tr, width = _trace_width)
        logging.info('background substracted')
        sci_spectrum: Spectrum1D = extract.spectrum
        logging.info('spectrum extracted')

        self._spec2d_axe.imshow(bg.bkg_wimage, origin='lower', aspect='auto', cmap=plt.cm.gray, alpha=0.2)
        self._spec2d_axe.imshow(sci_tr.image.data, origin='lower', aspect='auto', cmap=plt.cm.gray, alpha=0.2)
        self._spec2d_axe.step(sci_spectrum.spectral_axis, sci_tr.trace , color='white', linewidth = '0.3')
        self._spec2d_axe.step(sci_spectrum.spectral_axis, sci_tr.trace + extract.width , color='white', linestyle='dashed', alpha=0.2)
        self._spec2d_axe.step(sci_spectrum.spectral_axis, sci_tr.trace - extract.width , color='white', linestyle='dashed', alpha=0.2)
        #self._spec1d_axe.set_title('spectrum2D + background + trace fitted')

        self._spec1d_axe.step(sci_spectrum.spectral_axis , sci_spectrum.flux, color='red', linewidth = '0.3')
        self._spec1d_axe.set_xlabel('Pixels')
        self._spec1d_axe.set_ylabel('ADU')

        logging.info('calibrating neon spectrum...')

        self._spec1d_axe.clear()

        wavelength = [6506.53, 6532.88, 6598.95, 6678.28, 6717.04]*u.AA
        pixels = [770, 1190, 2240, 3484, 4160]*u.pix

        #line_list = QTable([pixels, wavelength], names=["pixel_center", "wavelength"])
        #input_spectrum, matched_line_list=None, line_pixels=None, line_wavelengths=None, catalog=None, input_model=Linear1D(), fitter=None
        cal = WavelengthCalibration1D(input_spectrum = sci_spectrum,
            #matched_line_list = line_list,
            line_wavelengths = wavelength,
            line_pixels = pixels,
            #input_model = models.Polynomial1D(degree = 2),
            #fitter = fitting.LMLSQFitter()
            fitter = fitting.LinearLSQFitter()
            )
        
        print('residuals :', str(cal.residuals) )
        print('fitted ', str(cal.fitted_model) )

        logging.info('spectrum calibrated')

        calibrated_spectrum = cal.apply_to_spectrum(sci_spectrum)

        self._spec1d_axe.set_ylabel('Relative intensity')
        self._spec1d_axe.set_xlabel('Wavelength (ang)')

        self._spec1d_axe.step(calibrated_spectrum.spectral_axis, calibrated_spectrum.flux, color='green', linewidth = '0.3')

        #logging.info('residuals :', cal.residuals )
        #logging.info('fitted ', cal.fitted_model )






