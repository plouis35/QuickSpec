"""
SpectrumModel — spectra processing logic

Holds the state of the current spectrum pipeline and exposes all processing steps.
"""
import logging
import numpy as np
from pathlib import Path

from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.modeling import models, fitting
from specutils import Spectrum
from specutils.analysis import snr_derived
from specreduce.tracing import FlatTrace, FitTrace

from app.config import Config
import spc.spc_utils as spc_utils


class SpectrumModel:
    """
    Holds all spectrum processing state and business logic.

    Attributes:
        science_spectrum: the current 1D spectrum (None until extracted or loaded)
        science_trace:    the current trace (None until traced)
    """

    def __init__(self) -> None:
        self.conf = Config()
        self.science_spectrum: Spectrum | None = None
        self.science_trace: FitTrace | FlatTrace | None = None

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def reset(self) -> None:
        """Clear all processing state."""
        self.science_spectrum = None
        self.science_trace = None

    @staticmethod
    def get_object_name(header) -> str:
        """
        Extract object name from a FITS header.

        Args:
            header: FITS header (or None)

        Returns:
            str: object name, or 'no_name' if not found
        """
        if header is None:
            return 'no_name'
        if 'OBJNAME' in header:
            return header['OBJNAME']
        if 'OBJECT' in header:
            return header['OBJECT']
        return 'no_name'

    def is_calibrated(self) -> bool:
        """Return True if the current spectrum has a wavelength axis (not pixels)."""
        if self.science_spectrum is None:
            return False
        return self.science_spectrum.spectral_axis.unit != u.Unit('pixel')

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------

    def load_spectrum(self, path: str) -> Spectrum | None:
        """
        Read a 1D spectrum from a FITS or DAT file.

        Args:
            path (str): file path

        Returns:
            Spectrum | None: loaded spectrum, or None on error
        """
        suffix = Path(path).suffix.lower()

        if suffix == '.dat':
            try:
                arr = np.loadtxt(path)
                spec = Spectrum(
                    spectral_axis=arr[:, 0] * u.Unit('Angstrom'),
                    flux=arr[:, 1] * u.Unit('mJy'),
                )
                self.science_spectrum = spec
                return spec
            except Exception as e:
                logging.error(f"{path}: {e}")
                return None

        # FITS format
        try:
            hdus = fits.open(path)
            if len(hdus) == 1:
                spec = Spectrum.read(path)
            elif len(hdus) == 2:
                ccd = CCDData.read(path, hdu=1, unit=u.Unit('adu'))
                spec = Spectrum(
                    spectral_axis=ccd.data['wavelength'] * u.Unit('Angstrom'),
                    flux=ccd.data['flux'] * u.Unit('mJy'),
                )
            else:
                logging.error(f"{path}: no data found in HDUs")
                return None

            self.science_spectrum = spec
            return spec

        except Exception as e:
            logging.error(f"{path}: {e}")
            return None

    # ------------------------------------------------------------------
    # Processing steps
    # ------------------------------------------------------------------

    def trace(self, img_stacked: CCDData) -> FitTrace | FlatTrace | None:
        """
        Find the spectral trace in a 2D image.

        Args:
            img_stacked (CCDData): stacked 2D spectrum

        Returns:
            FitTrace | FlatTrace | None: computed trace, also stored in self.science_trace
        """
        if img_stacked is None:
            logging.error("please reduce image(s) before fitting trace")
            return None

        self.science_trace = None
        _method = self.conf.get_str('processing', 'trace_method')

        if _method is None:
            logging.error("please define trace method in configuration file")
            return None

        if _method == 'fit':
            _model_str = self.conf.get_str('processing', 'trace_model')
            _model = spc_utils._parse_model(
                _model_str, spc_utils.ALLOWED_TRACE_MODELS,
                models.Polynomial1D(degree=2), 'trace_model')
            _peak   = self.conf.get_str('processing', 'peak_model') or 'gaussian'
            _bins   = self.conf.get_int('processing', 'trace_x_bins') or 12
            _window = self.conf.get_int('processing', 'trace_y_window') or 50
            _guess  = self.conf.get_float('processing', 'trace_y_guess')

            result = spc_utils.trace_spectrum(
                img_stacked, mode='fit',
                bins=_bins, guess=_guess, window=_window,
                trace_model=_model, peak_method=_peak,
            )

        elif _method == 'flat':
            _guess = self.conf.get_float('processing', 'trace_y_guess')
            if _guess is None:
                logging.error("please define 'trace_y_guess' for flat trace mode")
                return None
            result = spc_utils.trace_spectrum(img_stacked, mode='flat', guess=_guess)

        else:
            logging.error(f"unknown trace method: {_method}")
            return None

        self.science_trace = result
        return result

    def extract(self, img_stacked: CCDData) -> Spectrum | None:
        """
        Extract a 1D spectrum from the 2D image using the current trace.

        Args:
            img_stacked (CCDData): stacked 2D spectrum

        Returns:
            Spectrum | None: uncalibrated 1D spectrum
        """
        if img_stacked is None or self.science_trace is None:
            logging.error("please fit trace before extracting spectrum")
            return None

        result = spc_utils.extract_spectrum(img_stacked, self.science_trace)
        if result is not None:
            self.science_spectrum = result
        return result

    def calibrate(self, img_stacked: CCDData) -> Spectrum | None:
        """
        Wavelength-calibrate the current spectrum.

        Args:
            img_stacked (CCDData): used only to retrieve the FITS header

        Returns:
            Spectrum | None: calibrated spectrum, with optional wavelength shift applied
        """
        if self.science_spectrum is None:
            logging.error("please extract spectrum before calibrating")
            return None

        _wavelength = self.conf.get_str('processing', 'calib_x_wavelength')
        _pixels     = self.conf.get_str('processing', 'calib_x_pixel')
        if _wavelength is None or _pixels is None:
            logging.error("please specify lines/waves index in config file")
            return None

        calibrated = spc_utils.calibrate_spectrum(self.science_spectrum, _pixels, _wavelength)
        if calibrated is None:
            return None

        shift = self.conf.get_float('post_processing', 'shift_wavelength') or 0.0
        logging.info(f"shifting wavelength by {shift} Å")
        result = Spectrum(
            spectral_axis=calibrated.spectral_axis + shift * u.Unit('Angstrom'),
            flux=calibrated.flux,
        )
        self.science_spectrum = result
        logging.info(f'snr = {snr_derived(self.science_spectrum):.1f}')
        return result

    def apply_response(self) -> Spectrum | None:
        """
        Apply an instrumental response correction to the current spectrum.

        Returns:
            Spectrum | None: response-corrected spectrum
        """
        if self.science_spectrum is None:
            logging.error("please calibrate spectrum before applying response file")
            return None

        result = spc_utils.apply_response(self.science_spectrum)
        if result is not None:
            self.science_spectrum = result
        return result

    def smooth(self) -> Spectrum | None:
        """
        Apply median smooth, crop and normalization to the current spectrum.

        Returns:
            Spectrum | None: post-processed spectrum
        """
        if self.science_spectrum is None:
            logging.error("please calibrate spectrum before smoothing")
            return None

        spectrum = self.science_spectrum

        # median smooth
        smooth_width = self.conf.get_int('post_processing', 'median_smooth')
        if smooth_width is not None:
            smoothed = spc_utils.apply_median_smooth(spectrum, smooth_width)
            if smoothed is not None:
                spectrum = smoothed
                logging.info(f"median smooth applied={smooth_width}")
        else:
            logging.info("no median smooth to apply")

        # crop
        crop_str = self.conf.get_str('post_processing', 'crop_region')
        if crop_str is not None:
            parsed = spc_utils._parse_float_list(crop_str, 'crop_region')
            if parsed is None or len(parsed) < 2:
                logging.error(f"invalid crop_region value: {crop_str}")
                return None
            mn, mx = parsed[0], parsed[1]
            logging.info(f"crop to ({mn}, {mx}) Å")
            spectrum = spectrum[mn * u.Unit('angstrom'): mx * u.Unit('angstrom')]
        else:
            logging.info("no crop regions defined")

        # normalize
        norm_str = self.conf.get_str('processing', 'normalized_region')
        if norm_str is None:
            mn_norm, mx_norm = 6500.0, 6520.0
        else:
            parsed = spc_utils._parse_float_list(norm_str, 'normalized_region')
            if parsed is None or len(parsed) < 2:
                logging.error(f"invalid normalized_region value: {norm_str}")
                return None
            mn_norm, mx_norm = parsed[0], parsed[1]

        logging.info(f"normalizing on ({mn_norm}, {mx_norm}) Å")
        try:
            norm_val = spectrum[
                mn_norm * u.Unit('angstrom'): mx_norm * u.Unit('angstrom')
            ].flux.mean()
            spectrum = Spectrum(
                spectral_axis=spectrum.wavelength,
                flux=spectrum.flux / norm_val,
            )
            logging.info('spectrum normalized to 1')
        except Exception as e:
            logging.error(f"unable to normalize spectrum: {e}")
            return None

        self.science_spectrum = spectrum
        return spectrum
