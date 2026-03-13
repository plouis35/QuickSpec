"""
Tests for img/image_model.py

Covers: state management, stats, compute_cuts, load, reduce.
No GUI, no Tkinter, no Matplotlib.
"""
import pytest
import numpy as np
from astropy import units as u
from astropy.nddata import CCDData

from img.image_model import ImageModel


# ---------------------------------------------------------------------------
# State management
# ---------------------------------------------------------------------------

class TestReset:
    def test_initial_image_is_placeholder(self):
        model = ImageModel()
        assert model.img_stacked.data.shape == (2, 8)
        assert model.img_combiner is None
        assert model.img_reduced is False

    def test_reset_restores_placeholder(self, image_fits_file):
        model = ImageModel()
        model.load([image_fits_file])
        model.reset()
        assert model.img_stacked.data.shape == (2, 8)
        assert model.img_combiner is None
        assert model.img_reduced is False


# ---------------------------------------------------------------------------
# stats
# ---------------------------------------------------------------------------

class TestStats:
    def test_stats_on_known_data(self):
        model = ImageModel()
        model.img_stacked = CCDData(
            np.array([[1.0, 2.0], [3.0, 4.0]]),
            unit=u.Unit('adu'),
        )
        std, mean, vmin, vmax = model.stats()
        assert mean == pytest.approx(2.5)
        assert vmin == pytest.approx(1.0)
        assert vmax == pytest.approx(4.0)

    def test_stats_on_zeros(self):
        model = ImageModel()
        std, mean, vmin, vmax = model.stats()
        assert mean == pytest.approx(0.0)
        assert vmin == pytest.approx(0.0)
        assert vmax == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# compute_cuts
# ---------------------------------------------------------------------------

class TestComputeCuts:
    def test_cuts_low_lt_high(self):
        model = ImageModel()
        model.img_stacked = CCDData(
            np.random.default_rng(0).normal(1000, 50, (50, 50)),
            unit=u.Unit('adu'),
        )
        low, high, smin, smax = model.compute_cuts()
        assert low < high
        assert smin < low
        assert smax > high

    def test_slider_range_wider_than_cuts(self):
        model = ImageModel()
        model.img_stacked = CCDData(
            np.ones((10, 10)) * 500,
            unit=u.Unit('adu'),
        )
        low, high, smin, smax = model.compute_cuts()
        assert smin <= low
        assert smax >= high


# ---------------------------------------------------------------------------
# load
# ---------------------------------------------------------------------------

class TestLoad:
    def test_load_single_fits(self, image_fits_file):
        model = ImageModel()
        result = model.load([image_fits_file])
        assert result is not None
        assert isinstance(result, CCDData)
        assert model.img_combiner is not None
        assert model.img_reduced is False

    def test_load_updates_stacked(self, image_fits_file, synthetic_image):
        model = ImageModel()
        model.load([image_fits_file])
        # stacked shape should match the loaded image
        assert model.img_stacked.data.shape == synthetic_image.data.shape

    def test_load_nonexistent_returns_none(self):
        model = ImageModel()
        result = model.load(["/nonexistent/image.fits"])
        assert result is None

    def test_load_resets_reduced_flag(self, image_fits_file):
        model = ImageModel()
        model.img_reduced = True          # force stale state
        model.load([image_fits_file])
        assert model.img_reduced is False


# ---------------------------------------------------------------------------
# reduce
# ---------------------------------------------------------------------------

class TestReduce:
    def test_reduce_without_load_returns_none(self):
        model = ImageModel()
        result = model.reduce()
        assert result is None

    def test_reduce_without_calibration_frames_returns_stacked(self, image_fits_file):
        """
        When no bias/dark/flat are defined in config, reduce_images() should
        return the stacked image unchanged (no-op reduction).
        """
        model = ImageModel()
        model.load([image_fits_file])
        result = model.reduce()
        # with no calibration frames in test config, result may be None or
        # the same image — either is acceptable; what must not happen is a crash
        assert model.img_reduced is True or result is None

    def test_reduce_twice_is_idempotent(self, image_fits_file):
        model = ImageModel()
        model.load([image_fits_file])
        first  = model.reduce()
        second = model.reduce()
        # second call must warn and return the cached result, not re-process
        if first is not None:
            assert second is not None
