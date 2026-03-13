"""
Tests for spc/spc_utils.py

Covers: helpers internes (_parse_model, _parse_float_list),
et fonctions de traitement (apply_median_smooth, trace, extract).
"""
import pytest
import numpy as np
from astropy import units as u
from astropy.modeling import models
from specutils import Spectrum

from spc.spc_utils import (
    _parse_model,
    _parse_float_list,
    ALLOWED_TRACE_MODELS,
    ALLOWED_FITTERS,
    apply_median_smooth,
    trace_spectrum,
    extract_spectrum,
)


# ---------------------------------------------------------------------------
# _parse_float_list
# ---------------------------------------------------------------------------

class TestParseFloatList:
    def test_parses_two_values(self):
        result = _parse_float_list("6500, 6520", "test")
        assert result == pytest.approx([6500.0, 6520.0])

    def test_parses_without_spaces(self):
        result = _parse_float_list("100,200,300", "test")
        assert result == pytest.approx([100.0, 200.0, 300.0])

    def test_returns_none_on_invalid_input(self):
        result = _parse_float_list("not, a, number", "test")
        assert result is None

    def test_returns_none_on_empty_string(self):
        result = _parse_float_list("", "test")
        assert result is None


# ---------------------------------------------------------------------------
# _parse_model
# ---------------------------------------------------------------------------

class TestParseModel:
    def test_returns_known_model(self):
        result = _parse_model(
            'models.Polynomial1D(degree=2)',
            ALLOWED_TRACE_MODELS,
            models.Polynomial1D(degree=1),
            'trace_model',
        )
        assert isinstance(result, models.Polynomial1D)
        assert result.degree == 2

    def test_returns_default_for_unknown_key(self):
        default = models.Polynomial1D(degree=1)
        result = _parse_model(
            'models.UnknownModel()',
            ALLOWED_TRACE_MODELS,
            default,
            'trace_model',
        )
        assert result is default

    def test_returns_default_when_value_is_none(self):
        default = models.Polynomial1D(degree=1)
        result = _parse_model(None, ALLOWED_TRACE_MODELS, default, 'trace_model')
        assert result is default


# ---------------------------------------------------------------------------
# apply_median_smooth
# ---------------------------------------------------------------------------

class TestApplyMedianSmooth:
    def test_smooth_returns_spectrum(self, flat_spectrum):
        result = apply_median_smooth(flat_spectrum, smooth_width=5)
        assert isinstance(result, Spectrum)

    def test_smooth_preserves_length(self, noisy_spectrum):
        result = apply_median_smooth(noisy_spectrum, smooth_width=3)
        assert len(result.flux) == len(noisy_spectrum.flux)

    def test_smooth_reduces_noise(self, noisy_spectrum):
        """Smoothed spectrum should have lower std than the original."""
        smoothed = apply_median_smooth(noisy_spectrum, smooth_width=11)
        assert smoothed.flux.std() < noisy_spectrum.flux.std()

    def test_smooth_width_one_is_identity(self, flat_spectrum):
        result = apply_median_smooth(flat_spectrum, smooth_width=1)
        np.testing.assert_allclose(result.flux.value, flat_spectrum.flux.value, rtol=1e-5)


# ---------------------------------------------------------------------------
# trace_spectrum
# ---------------------------------------------------------------------------

class TestTraceSpectrum:
    def test_flat_trace_returns_result(self, synthetic_image):
        result = trace_spectrum(synthetic_image, mode='flat', guess=50)
        assert result is not None

    def test_fit_trace_on_synthetic_image(self, synthetic_image):
        result = trace_spectrum(
            synthetic_image, mode='fit',
            bins=10, guess=50, window=20,
        )
        assert result is not None

    def test_unknown_mode_returns_none(self, synthetic_image):
        result = trace_spectrum(synthetic_image, mode='unknown')
        assert result is None

    def test_flat_trace_centred_near_guess(self, synthetic_image):
        """Flat trace should sit exactly at the requested y position."""
        from specreduce.tracing import FlatTrace
        result = trace_spectrum(synthetic_image, mode='flat', guess=50)
        assert isinstance(result, FlatTrace)
        np.testing.assert_allclose(result.trace_pos, 50)


# ---------------------------------------------------------------------------
# extract_spectrum
# ---------------------------------------------------------------------------

class TestExtractSpectrum:
    def test_extract_returns_spectrum(self, synthetic_image):
        trace = trace_spectrum(synthetic_image, mode='flat', guess=50)
        result = extract_spectrum(synthetic_image, trace)
        assert result is not None
        assert isinstance(result, Spectrum)

    def test_extracted_length_matches_image_width(self, synthetic_image):
        trace = trace_spectrum(synthetic_image, mode='flat', guess=50)
        result = extract_spectrum(synthetic_image, trace)
        assert len(result.flux) == synthetic_image.data.shape[1]  # 200 columns
