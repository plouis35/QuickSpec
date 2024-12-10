import logging
import numpy as np

from specutils.spectra.spectrum1d import Spectrum1D
from specutils.manipulation import median_smooth, gaussian_smooth
from specutils.analysis import snr, snr_derived

from astropy.utils.exceptions import AstropyWarning
from astropy import units as u
from astropy.modeling import models, fitting
from astropy.nddata import NDDataRef

from specutils.spectra.spectrum1d import Spectrum1D

from specreduce.tracing import FlatTrace, FitTrace
from specreduce.background import Background
from specreduce.extract import BoxcarExtract
from specreduce.fluxcal import FluxCalibration
from specreduce import WavelengthCalibration1D

from app.config import Config

class SPCUtils(object):
    conf = Config()

    @staticmethod
    def trace_spectrum(img_stacked:np.ndarray) -> FitTrace | None:

      # peak_method : one of gaussian, centroid
      # trace_model : one of :
      #  ~astropy.modeling.polynomial.Chebyshev1D,
      #  ~astropy.modeling.polynomial.Legendre1D, 
      #  ~astropy.modeling.polynomial.Polynomial1D,
      #  ~astropy.modeling.spline.Spline1D, 
      #  optional The 1-D polynomial model used to fit the trace to the bins' peak pixels. 
      #  Spline1D models are fit with Astropy's 'SplineSmoothingFitter', while the other models are fit with the 'LevMarLSQFitter'.
      #  [default: models.Polynomial1D(degree=1)]

      try:
          science_trace:FitTrace = FitTrace(img_stacked, 
                                  bins=SPCUtils.conf.get_int('processing', 'trace_x_bins'), 
                                  trace_model=models.Chebyshev1D(degree=2), 
                                  #trace_model=models.Polynomial1D(degree=2),
                                  peak_method='centroid',     #'gaussian', 
                                  window=SPCUtils.conf.get_int('processing', 'trace_y_window'),
                                  guess=SPCUtils.conf.get_float('processing', 'trace_y_guess')
                                  )
          
      except Exception as e:
          logging.error(f"unable to fit trace : {e}")
          return None

      return science_trace



    @staticmethod
    def extract_spectrum(img_stacked:np.ndarray, science_trace: FitTrace) -> Spectrum1D | None:
        try:
            bg_trace: Background = Background.two_sided(img_stacked, 
                                                  science_trace, 
                                                  separation=SPCUtils.conf.get_int('processing', 'sky_y_offset'), 
                                                  width=SPCUtils.conf.get_int('processing', 'sky_y_size')  ) 
        except Exception as e:
            logging.error(f"unable to fit background : {e}")
            return None
        
        logging.info('background extracted')
                
        try:
            extracted_spectrum = BoxcarExtract(img_stacked - bg_trace, 
                                    science_trace, 
                                    width = SPCUtils.conf.get_float('processing', 'trace_y_size') )
        except Exception as e:
            logging.error(f"unable to extract background : {e}")
            return None

        logging.info('background substracted')

        try:
            extracted_spectrum = extracted_spectrum.spectrum
        except Exception as e:
            logging.error(f"unable to extract science spectrum : {e}")
            return None
        
        return extracted_spectrum


    @staticmethod
    def calibrate_spectrum(science_spectrum: Spectrum1D, pixels: str, wavelength: str) -> Spectrum1D | None:  
        """_summary_

        Args:
            science_spectrum (Spectrum1D): _description_
            pixels (str): _description_
            wavelength (str): _description_

        Returns:
            Spectrum1D | None: _description_
        """

        # convert args to pix/AA
        _wavelength = [float(x) for x in wavelength.replace(',', '').split()]*u.AA
        _pixels = [float(x) for x in pixels.replace(',', '').split()]*u.pix
        logging.info(f"calibrating pixels set : {pixels}")
        logging.info(f"with wavelengths set : {wavelength}")

        try:
            calibration = WavelengthCalibration1D(input_spectrum = science_spectrum,
                line_wavelengths = _wavelength,
                line_pixels = _pixels,
                input_model = models.Linear1D(),
                fitter = fitting.LinearLSQFitter(),
                #input_model = models.Polynomial1D(degree = 2),
                #fitter = fitting.LMLSQFitter(),
                )
        except Exception as e:
            logging.error(f"unable to calibrate spectrum : {e}")
            return None
        
        # calibrate science spectrum
        calibrated_spectrum: Spectrum1D = calibration.apply_to_spectrum(science_spectrum)
        logging.info(f"spectrum calibrated - residuals : {repr(calibration.residuals)}")

        # normalize to 1
        norm_region = calibrated_spectrum[6500 * u.AA: 6520 * u.AA].flux.mean() 
        normalized_spec = Spectrum1D(spectral_axis = calibrated_spectrum.wavelength, 
                                     flux = calibrated_spectrum.flux / norm_region)  
        logging.info('spectrum normalized to 1')

        return normalized_spec
    
    @staticmethod
    def apply_response(science_spectrum: Spectrum1D) -> Spectrum1D | None:  
        # apply response file if any defined

        final_spec: Spectrum1D = None

        try:
            if (respFile := SPCUtils.conf.get_str('processing', 'response_file')) is not None:
                respFile = f"{SPCUtils.conf.get_conf_directory()}/{respFile}"
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

                # compute airmass
                #target_coord = SkyCoord.from_name(TARGET)
                #target_time = Time(science_spectrum.header['DATE-OBS'])
                #obs_coord = EarthLocation(lon = OBS_LONGITUDE * u.deg, lat = OBS_LATITUDE * u.deg)
                #altaz = AltAz(obstime=target_time, location = obs_coord)
                #ZD = target_coord.transform_to(AltAz(obstime = target_time, location = obs_coord)).zen
                #airmass = 1.0 / np.cos(ZD)
                #logging.info(f"computed {ZD=}, {airmass=}")

                #cal_spec = FluxCalibration.airmass_cor(science_spectrum)
                #fc = FluxCalibration()
                #sci_spec = fc.mag2flux(science_spectrum)
                #cal_spec = fc.standard_sensfunc(standard=_resp1d)
                #final_spec = fc.apply_sensfunc(_resp1d)

                logging.info('response applied')
            else:
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
        logging.info(f"median smooth applied={smooth_width}")

        return science_spectrum

if __name__ == "__main__":
    pass
