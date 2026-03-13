"""
Utility functions to trace, extract and calibrate 2D spectra.
"""
import logging

from specutils import Spectrum
from specutils.manipulation import median_smooth
from specutils.manipulation import FluxConservingResampler

from astropy import units as u
from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.nddata import CCDData

from specreduce.tracing import FlatTrace, FitTrace
from specreduce.background import Background
from specreduce.extract import BoxcarExtract
from specreduce import WavelengthCalibration1D

from app.config import Config

# ---------------------------------------------------------------------------
# Safe model / fitter resolution
# ---------------------------------------------------------------------------

ALLOWED_TRACE_MODELS: dict = {
    'models.Polynomial1D(degree=2)': models.Polynomial1D(degree=2),
    'models.Polynomial1D(degree=3)': models.Polynomial1D(degree=3),
    'models.Chebyshev1D(degree=2)':  models.Chebyshev1D(degree=2),
    'models.Chebyshev1D(degree=3)':  models.Chebyshev1D(degree=3),
    'models.Legendre1D(degree=2)':   models.Legendre1D(degree=2),
    'models.Legendre1D(degree=3)':   models.Legendre1D(degree=3),
}

ALLOWED_FITTERS: dict = {
    'fitting.LMLSQFitter()':     fitting.LMLSQFitter(),
    'fitting.LinearLSQFitter()': fitting.LinearLSQFitter(),
    'fitting.LevMarLSQFitter()': fitting.LevMarLSQFitter(),
}


def _parse_model(value: str | None, allowed: dict, default, label: str):
    """
    Safe replacement for eval() on model/fitter config values.
    """
    if value is None:
        return default
    if value not in allowed:
        logging.warning(f"unknown {label} '{value}', falling back to default")
        return default
    return allowed[value]


def _parse_float_list(value: str, label: str) -> list[float] | None:
    """
    Safely parse a comma-separated list of floats from a config string.
    """
    try:
        return [float(x.strip()) for x in value.split(',')]
    except (ValueError, AttributeError) as e:
        logging.error(f"invalid {label} value '{value}': {e}")
        return None


# ---------------------------------------------------------------------------
# Core processing functions
# ---------------------------------------------------------------------------

def trace_spectrum(
    img_stacked: CCDData,
    mode: str = 'fit',
    bins: int | None = None,
    guess: float | None = None,
    window: int | None = None,
    trace_model=None,
    peak_method: str | None = 'gaussian',
) -> FitTrace | FlatTrace | None:
    """
    Identify the trace in a 2D spectrum image.

    Two modes:
    - 'fit'  : automatic trace fitting
    - 'flat' : straight horizontal line at a fixed y position

    Args:
        img_stacked (CCDData): stacked 2D spectrum image
        mode (str): 'fit' or 'flat'
        bins (int | None): number of bins along dispersion axis (fit mode)
        guess (float | None): initial y guess for trace position
        window (int | None): half-window around guess (fit mode)
        trace_model: astropy model instance (fit mode)
        peak_method (str | None): 'gaussian', 'centroid', or 'max' (fit mode)

    Returns:
        FitTrace | FlatTrace | None
    """
    if trace_model is None:
        trace_model = models.Polynomial1D(degree=2)

    logging.info(f"trace args: {mode=}, {bins=}, {trace_model=}, {peak_method=}, {window=}, {guess=}")

    if mode == 'fit':
        try:
            science_trace = FitTrace(
                image=img_stacked,
                bins=bins,
                trace_model=trace_model,
                peak_method=peak_method,
                window=window,
                guess=guess,
            )
            logging.info(f'trace fitted: automatic model = {science_trace.trace_model_fit}')
            return science_trace
        except Exception as e:
            logging.error(f"unable to fit trace: {e}")
            return None

    elif mode == 'flat':
        try:
            science_trace = FlatTrace(image=img_stacked, trace_pos=guess)
            logging.info(f'trace fitted: flat model = {science_trace.trace}')
            return science_trace
        except Exception as e:
            logging.error(f"unable to flat trace: {e}")
            return None

    else:
        logging.error(f"unknown trace method: {mode}")
        return None


def extract_spectrum(
    img_stacked: CCDData,
    science_trace: FitTrace | FlatTrace,
) -> Spectrum | None:
    """
    Extract a 1D spectrum around the fitted trace.

    Args:
        img_stacked (CCDData): 2D spectrum image
        science_trace (FitTrace | FlatTrace): trace from trace_spectrum()

    Returns:
        Spectrum | None: uncalibrated 1D spectrum
    """
    conf = Config()

    if conf.get_bool('processing', 'sky_substract') in (None, True):
        try:
            bg_trace = Background.two_sided(
                img_stacked,
                science_trace,
                separation=conf.get_int('processing', 'sky_y_offset'),
                width=conf.get_int('processing', 'sky_y_size'),
            )
        except Exception as e:
            logging.error(f"unable to fit background: {e}")
            return None

        logging.info('background extracted')

        try:
            _extract = BoxcarExtract(
                img_stacked - bg_trace,
                science_trace,
                width=conf.get_float('processing', 'trace_y_size'),
            )
        except Exception as e:
            logging.error(f"unable to extract spectrum: {e}")
            return None

        logging.info('background subtracted')

    else:
        try:
            _extract = BoxcarExtract(
                img_stacked,
                science_trace,
                width=conf.get_float('processing', 'trace_y_size'),
            )
        except Exception as e:
            logging.error(f"unable to extract spectrum: {e}")
            return None

        logging.info("no background to subtract")

    try:
        return _extract.spectrum
    except Exception as e:
        logging.error(f"unable to extract science spectrum: {e}")
        return None


def calibrate_spectrum(
    science_spectrum: Spectrum,
    pixels: str,
    wavelength: str,
) -> Spectrum | None:
    """
    Calibrate a 1D spectrum using pixel/wavelength reference pairs.

    Args:
        science_spectrum (Spectrum): uncalibrated 1D spectrum
        pixels (str): comma-separated pixel positions
        wavelength (str): comma-separated wavelength values (Angstrom)

    Returns:
        Spectrum | None: calibrated and normalized 1D spectrum
    """
    conf = Config()

    _wavelength = [float(x) for x in wavelength.replace(',', '').split()] * u.Unit('angstrom')
    _pixels = [float(x) for x in pixels.replace(',', '').split()] * u.Unit('pixel')
    logging.info(f"calibrating pixels set: {pixels}")
    logging.info(f"with wavelengths set: {wavelength}")

    _input_model_str = conf.get_str('processing', 'input_model') or "models.Polynomial1D(degree=2)"
    _input_model = _parse_model(_input_model_str, ALLOWED_TRACE_MODELS, models.Polynomial1D(degree=2), 'input_model')

    _fitter_model_str = conf.get_str('processing', 'fitter_model') or "fitting.LMLSQFitter()"
    _fitter_model = _parse_model(_fitter_model_str, ALLOWED_FITTERS, fitting.LMLSQFitter(), 'fitter_model')

    try:
        calibration = WavelengthCalibration1D(
            input_spectrum=science_spectrum,
            line_wavelengths=_wavelength,
            line_pixels=_pixels,
            fitter=_fitter_model,
            input_model=_input_model,
        )
    except Exception as e:
        logging.error(f"unable to calibrate spectrum: {e}")
        return None

    try:
        calibrated = calibration.apply_to_spectrum(science_spectrum)
        logging.info(f"spectrum calibrated - residuals: {repr(calibration.residuals)}")
        logging.info(f"spectrum calibrated - fitted model: {repr(calibration.fitted_model)}")
    except Exception as e:
        logging.error(f"unable to apply calibration to spectrum: {e}")
        return None

    # normalize to 1
    if (_norm_region := conf.get_str('processing', 'normalized_region')) is None:
        _min_norm, _max_norm = 6500.0, 6520.0
    else:
        _parsed = _parse_float_list(_norm_region, 'normalized_region')
        if _parsed is None or len(_parsed) < 2:
            logging.error(f"invalid normalized_region value: {_norm_region}")
            return None
        _min_norm, _max_norm = _parsed[0], _parsed[1]

    logging.info(f"normalized regions to use = ({_min_norm}, {_max_norm})")

    try:
        norm_value = calibrated[
            _min_norm * u.Unit('angstrom'): _max_norm * u.Unit('angstrom')
        ].flux.mean()
        normalized = Spectrum(
            spectral_axis=calibrated.wavelength,
            flux=calibrated.flux / norm_value,
        )
        logging.info('spectrum normalized to 1')
        return normalized
    except Exception as e:
        logging.error(f"unable to normalize spectrum: {e}")
        return None


def apply_response(science_spectrum: Spectrum) -> Spectrum | None:
    """
    Divide a calibrated 1D spectrum by an instrumental response file.

    Args:
        science_spectrum (Spectrum): calibrated 1D spectrum

    Returns:
        Spectrum | None: response-corrected 1D spectrum
    """
    conf = Config()

    try:
        resp_file = conf.get_str('processing', 'response_file')
        if resp_file is None:
            logging.info("no response file to apply")
            return science_spectrum

        resp_path = f"{conf.get_conf_directory()}/{resp_file}"
        logging.info(f"opening {resp_path}...")

        hdus = fits.open(resp_path)
        if len(hdus) == 1:
            resp1d: Spectrum = Spectrum.read(resp_path)
        elif len(hdus) == 2:
            _spc = CCDData.read(resp_path, hdu=1, unit=u.Unit('adu'))
            resp1d = Spectrum(
                spectral_axis=_spc.data['wavelength'] * u.Unit('Angstrom'),
                flux=_spc.data['flux'] * u.Unit('mJy'),
            )
        else:
            logging.error(f"{resp_path}: no data found in HDUs")
            return None

        resampler = FluxConservingResampler(extrapolation_treatment='truncate')
        resp_resampled = resampler(resp1d, science_spectrum.spectral_axis)
        spec_resampled = resampler(science_spectrum, resp_resampled.spectral_axis)
        final_spec = spec_resampled / resp_resampled
        logging.info('response applied')
        return final_spec

    except Exception as e:
        logging.error(f"{e}")
        logging.error("no response file applied")
        return None


def apply_median_smooth(science_spectrum: Spectrum, smooth_width: int = 1) -> Spectrum | None:
    """
    Apply a median filter to a spectrum.

    Args:
        science_spectrum (Spectrum): spectrum to smooth
        smooth_width (int): kernel width in pixels

    Returns:
        Spectrum | None: smoothed spectrum
    """
    return median_smooth(science_spectrum, width=smooth_width)
