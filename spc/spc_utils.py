"""
utility routines to trace, extract and calibrate 2D spectra
"""
import logging
import numpy as np

from specutils.spectra.spectrum1d import Spectrum1D
from specutils.manipulation import median_smooth, gaussian_smooth
from specutils.analysis import snr, snr_derived

from astropy.utils.exceptions import AstropyWarning
from astropy import units as u
from astropy.modeling import models, fitting
from astropy.nddata import NDDataRef
from astropy.nddata import CCDData

from specutils.spectra.spectrum1d import Spectrum1D

from specreduce.tracing import FlatTrace, FitTrace
from specreduce.background import Background
from specreduce.extract import BoxcarExtract
from specreduce.fluxcal import FluxCalibration
from specreduce import WavelengthCalibration1D

from app.config import Config

class SPCUtils(object):

    @staticmethod
    def trace_spectrum(img_stacked:CCDData) -> FitTrace | FlatTrace | None:
        """
        identify trace in a 2D spectrum image
        two methods implemented from specreduce: 
        1) 'flat' : just a straight line along the X-axis
        2) 'fit' : split the spectrum into bins and identify their y-position

        Args:
            img_stacked (np.ndarray): 2D image to scan

        Returns:
            FitTrace | FlatTrace |  None: specreduce Fit/flatTrace class - will be used to extract spectrum later
        """        
        conf = Config()

        if (_trace_method := conf.get_str('processing', 'trace_method')) is not None:
            if _trace_method == 'fit': 
                try:
                    if (_trace_model := conf.get_str('processing', 'trace_model')) is None:
                        _trace_model = "models.Polynomial1D(degree=2)"

                    if (_peak_method := conf.get_str('processing', 'peak_model')) is None:
                        _peak_method = "max"

                    science_trace = FitTrace(image=img_stacked, 
                                            bins=conf.get_int('processing', 'trace_x_bins'), 
                                            trace_model=eval(_trace_model),
                                            peak_method=_peak_method,
                                            window=conf.get_int('processing', 'trace_y_window'),
                                            guess=conf.get_float('processing', 'trace_y_guess')
                                            )
                    logging.info(f'trace fitted : trace model fitted = {science_trace.trace_model_fit}')

                except Exception as e:
                    logging.error(f"unable to fit trace : {e}")
                    return None

            elif _trace_method == 'flat': 
                if conf.get_float('processing', 'trace_y_guess') is None:
                    logging.error("please define a 'trace_y_guess' for flat trace mode")
                    return None
                try:
                    science_trace = FlatTrace(image=img_stacked, 
                                            trace_pos=conf.get_float('processing', 'trace_y_guess')
                                            )
                except Exception as e:
                    logging.error(f"unable to flat trace : {e}")
                    return None

            else: 
                logging.error(f"unknown trace method: {_trace_method}")
                return None
        else:
            logging.error("please define trace method in configuration file")
            return None
            
        logging.info(f'trace fitted : y = {science_trace.trace}')

        return science_trace



    @staticmethod
    def extract_spectrum(img_stacked:CCDData, science_trace: FitTrace) -> Spectrum1D | None:
        """
        extract 1D spectrum arround trace fitted
        Args:
            img_stacked (np.ndarray): 2D spectrum image
            science_trace (FitTrace): identified trace returned by trace_spectrum

        Returns:
            Spectrum1D | None: uncalibrated 1D spectrum 
        """      
        conf = Config()

        if conf.get_bool('processing', 'sky_substract') in (None, True):
            try:
                bg_trace: Background = Background.two_sided(img_stacked, 
                                                    science_trace, 
                                                    separation=conf.get_int('processing', 'sky_y_offset'), 
                                                    width=conf.get_int('processing', 'sky_y_size')  ) 
            except Exception as e:
                logging.error(f"unable to fit background : {e}")
                return None
            
            logging.info('background extracted')
                    
            try:
                _extract = BoxcarExtract(img_stacked - bg_trace, 
                                        science_trace,
                                        width = conf.get_float('processing', 'trace_y_size') )
            except Exception as e:
                logging.error(f"unable to extract spectrum : {e}")
                return None
            
            logging.info('background substracted')
        else:
            try:
                _extract = BoxcarExtract(img_stacked, 
                                        science_trace,
                                        width = conf.get_float('processing', 'trace_y_size') )
            except Exception as e:
                logging.error(f"unable to extract spectrum : {e}")
                return None
            
            logging.info("no background to substract")
        
        try:
            _extracted_spectrum: Spectrum1D = _extract.spectrum
        except Exception as e:
            logging.error(f"unable to extract science spectrum : {e}")
            return None
        
        return _extracted_spectrum


    @staticmethod
    def calibrate_spectrum(science_spectrum: Spectrum1D, pixels: str, wavelength: str) -> Spectrum1D | None:  
        """
        calibrate 1D spectrum according to a list of pixel / wavelength parameters 

        Args:
            science_spectrum (Spectrum1D): uncalibrated 1D spectrum
            pixels (str): list of pixel positions
            wavelength (str): list of corresponding wavelength values

        Returns:
            Spectrum1D | None: calibrated 1D spectrum
        """
        conf = Config()

        # convert args to pix/AA
        _wavelength = [float(x) for x in wavelength.replace(',', '').split()]*u.Unit('angstrom')
        _pixels = [float(x) for x in pixels.replace(',', '').split()]*u.Unit('pixel')
        logging.info(f"calibrating pixels set : {pixels}")
        logging.info(f"with wavelengths set : {wavelength}")

        if (_input_model := conf.get_str('processing', 'input_model')) is None:
            _input_model = "models.Polynomial1D(degree=2)"

        if (_fitter_model := conf.get_str('processing', 'fitter_model')) is None:
            _fitter_model = "fitting.LMLSQFitter()"

        try:
            calibration = WavelengthCalibration1D(input_spectrum = science_spectrum,
                line_wavelengths = _wavelength,
                line_pixels = _pixels,
                fitter = eval(_fitter_model),
                input_model = eval(_input_model)
                )
        except Exception as e:
            logging.error(f"unable to calibrate spectrum : {e}")
            return None
        
        # calibrate science spectrum
        calibrated_spectrum: Spectrum1D = calibration.apply_to_spectrum(science_spectrum)
        logging.info(f"spectrum calibrated - residuals : {repr(calibration.residuals)}")
        logging.info(f"spectrum calibrated - fitted model : {repr(calibration.fitted_model )}")

        # normalize to 1 
        # TODO: make sure the normalized range is within the calibration range !
        norm_region = calibrated_spectrum[6500 * u.Unit('angstrom'): 6520 * u.Unit('angstrom')].flux.mean() 
        normalized_spec = Spectrum1D(spectral_axis = calibrated_spectrum.wavelength, 
                                     flux = calibrated_spectrum.flux / norm_region)  
        logging.info('spectrum normalized to 1')

        return normalized_spec
    
    @staticmethod
    def apply_response(science_spectrum: Spectrum1D) -> Spectrum1D | None:  
        """
        divide calibrated 1D spectrum by a response file

        Args:
            science_spectrum (Spectrum1D): calibrated 1D spectrum

        Returns:
            Spectrum1D | None: response corrected 1D spectrum
        """        
        conf = Config()

        final_spec: Spectrum1D = None

        try:
            if (respFile := conf.get_str('processing', 'response_file')) is not None:
                respFile = f"{conf.get_conf_directory()}/{respFile}"
                logging.info(f"opening {respFile}...")
                resp1d: Spectrum1D = Spectrum1D.read(respFile)
                _factor = int(resp1d.shape[0] / science_spectrum.shape[0])
                logging.info(f"{science_spectrum.shape[0]=}, {resp1d.shape[0]=}, {_factor=}")
                
                _resp1d_ndd = NDDataRef(resp1d)
                _resp1d_ndd.wcs = None

                # rebin response to science spectrum axis size
                _resp1d = _resp1d_ndd[::_factor].data

                # apply response
                final_spec = science_spectrum / _resp1d
                logging.info('response applied')
            else:
                final_spec = science_spectrum
                logging.info("no response file to apply")

        except Exception as e:
            logging.error(f"{e}")
            logging.error("no response file applied")
        
        return final_spec


    @staticmethod
    def median_smooth(science_spectrum: Spectrum1D, smooth_width: int = 1) -> Spectrum1D | None:
        """
        apply a median filter to a spectrum

        Args:
            science_spectrum (Spectrum1D): spectrum to smooth
        Returns:
            Spectrum1D | None: smooth'ed spectrum
        """
        smooth_spec: Spectrum1D = median_smooth(science_spectrum, width=smooth_width) 
        science_spectrum = smooth_spec
        return science_spectrum

if __name__ == "__main__":
    pass
