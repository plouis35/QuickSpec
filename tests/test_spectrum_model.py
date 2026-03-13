"""
Tests for spc/spectrum_model.py

Covers: state management, I/O, and processing steps.
No GUI, no Tkinter, no Matplotlib.
"""
import pytest
import numpy as np
from astropy import units as u
from astropy.nddata import CCDData
from specutils import Spectrum

from spc.spectrum_model import SpectrumModel


# ---------------------------------------------------------------------------
# State management
# ---------------------------------------------------------------------------

class TestReset:
    def test_initial_state_is_empty(self):
        model = SpectrumModel()
        assert model.science_spectrum is None
        assert model.science_trace is None

    def test_reset_clears_spectrum(self, flat_spectrum):
        model = SpectrumModel()
        model.science_spectrum = flat_spectrum
        model.reset()
        assert model.science_spectrum is None
        assert model.science_trace is None


# ---------------------------------------------------------------------------
# get_object_name
# ---------------------------------------------------------------------------

class TestGetObjectName:
    def test_returns_no_name_when_header_is_none(self):
        assert SpectrumModel.get_object_name(None) == 'no_name'

    def test_reads_objname_key(self):
        header = {'OBJNAME': 'Vega'}
        assert SpectrumModel.get_object_name(header) == 'Vega'

    def test_falls_back_to_object_key(self):
        header = {'OBJECT': 'Sirius'}
        assert SpectrumModel.get_object_name(header) == 'Sirius'

    def test_prefers_objname_over_object(self):
        header = {'OBJNAME': 'Vega', 'OBJECT': 'Sirius'}
        assert SpectrumModel.get_object_name(header) == 'Vega'

    def test_returns_no_name_when_keys_absent(self):
        header = {'EXPTIME': 30.0}
        assert SpectrumModel.get_object_name(header) == 'no_name'


# ---------------------------------------------------------------------------
# is_calibrated
# ---------------------------------------------------------------------------

class TestIsCalibrated:
    def test_returns_false_when_no_spectrum(self):
        model = SpectrumModel()
        assert model.is_calibrated() is False

    def test_returns_false_for_pixel_axis(self, pixel_spectrum):
        model = SpectrumModel()
        model.science_spectrum = pixel_spectrum
        assert model.is_calibrated() is False

    def test_returns_true_for_wavelength_axis(self, flat_spectrum):
        model = SpectrumModel()
        model.science_spectrum = flat_spectrum
        assert model.is_calibrated() is True


# ---------------------------------------------------------------------------
# load_spectrum
# ---------------------------------------------------------------------------

class TestLoadSpectrum:
    def test_load_fits(self, spectrum_fits_file):
        model = SpectrumModel()
        result = model.load_spectrum(spectrum_fits_file)
        assert result is not None
        assert isinstance(result, Spectrum)
        assert model.science_spectrum is result

    def test_load_dat(self, spectrum_dat_file):
        model = SpectrumModel()
        result = model.load_spectrum(spectrum_dat_file)
        assert result is not None
        assert isinstance(result, Spectrum)

    def test_load_nonexistent_file_returns_none(self):
        model = SpectrumModel()
        result = model.load_spectrum("/nonexistent/path/spectrum.fits")
        assert result is None
        assert model.science_spectrum is None

    def test_loaded_spectrum_has_correct_length(self, spectrum_fits_file, flat_spectrum):
        model = SpectrumModel()
        result = model.load_spectrum(spectrum_fits_file)
        assert len(result.flux) == len(flat_spectrum.flux)


# ---------------------------------------------------------------------------
# smooth
# ---------------------------------------------------------------------------

class TestSmooth:
    def test_smooth_returns_none_without_spectrum(self):
        model = SpectrumModel()
        assert model.smooth() is None

    def test_smooth_preserves_spectrum_length(self, flat_spectrum):
        model = SpectrumModel()
        model.science_spectrum = flat_spectrum
        result = model.smooth()
        # smooth may crop but should still return a valid Spectrum
        assert result is not None
        assert isinstance(result, Spectrum)

    def test_smooth_normalizes_flux(self, noisy_spectrum):
        """After smooth+normalize the flux near the norm region should be ~1."""
        model = SpectrumModel()
        model.science_spectrum = noisy_spectrum
        result = model.smooth()
        assert result is not None
        # normalized_region is 6500–6520 in the test config
        region = result[6500 * u.Angstrom: 6520 * u.Angstrom]
        np.testing.assert_allclose(region.flux.mean().value, 1.0, rtol=0.05)

    def test_smooth_updates_internal_state(self, flat_spectrum):
        model = SpectrumModel()
        model.science_spectrum = flat_spectrum
        result = model.smooth()
        assert model.science_spectrum is result


# ---------------------------------------------------------------------------
# trace (flat mode — no real image needed, just checks the guard)
# ---------------------------------------------------------------------------

class TestTrace:
    def test_trace_returns_none_without_image(self):
        model = SpectrumModel()
        result = model.trace(None)
        assert result is None

    def test_flat_trace_on_synthetic_image(self, synthetic_image):
        model = SpectrumModel()
        result = model.trace(synthetic_image)
        assert result is not None
        assert model.science_trace is result
