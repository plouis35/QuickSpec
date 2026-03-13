"""
conftest.py — shared pytest fixtures for QuickSpec tests.

All fixtures produce synthetic data (no real FITS files needed).
They are available in all test modules automatically.
"""
import os
import sys
import pytest
import numpy as np

from astropy import units as u
from astropy.nddata import CCDData
from specutils import Spectrum

# Make sure the project root is on sys.path when running pytest from any dir
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# ---------------------------------------------------------------------------
# Config isolation
# ---------------------------------------------------------------------------
@pytest.fixture(autouse=True)
def isolated_config(tmp_path, monkeypatch):
    """
    Write a minimal quickspec.ini in a temp dir and redirect Config there.

    Strategy:
      1. Patch the class-level _configDir / _configFile BEFORE any Config()
         call, so _initialize() reads our temp file.
      2. Reset _instance so the singleton is re-created fresh each test.
      3. Restore everything via monkeypatch + explicit cleanup in the teardown.
    """
    cfg_content = """\
[logger]
level = WARNING

[display]
theme = dark
contrast_level = 6
colorize_min_wl = 3800
colorize_max_wl = 7500
colorize_bin_size = 20

[pre_processing]
max_memory = 1e9

[processing]
trace_method = flat
trace_y_guess = 50
trace_y_size = 10
sky_substract = false
sky_y_offset = 20
sky_y_size = 10
calib_x_pixel = 100, 200, 300
calib_x_wavelength = 5000, 6000, 7000
normalized_region = 6500, 6520
input_model = models.Polynomial1D(degree=2)
fitter_model = fitting.LMLSQFitter()

[post_processing]
shift_wavelength = 0

[lines]
656.3 = Ha
486.1 = Hb
"""
    cfg_file = tmp_path / "quickspec.ini"
    cfg_file.write_text(cfg_content)

    from app.config import Config

    # Reset singleton then redirect class variables
    Config._instance = None
    monkeypatch.setattr(Config, '_configDir',  str(tmp_path))
    monkeypatch.setattr(Config, '_configFile', 'quickspec.ini')

    yield

    # Teardown: ensure next test gets a fresh singleton
    Config._instance = None


# ---------------------------------------------------------------------------
# Synthetic 1D spectrum fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def flat_spectrum():
    """Flat (flux=1.0) calibrated spectrum, 500 pts, 5000–7000 Å."""
    wave = np.linspace(5000, 7000, 500) * u.Angstrom
    flux = np.ones(500) * u.mJy
    return Spectrum(spectral_axis=wave, flux=flux)


@pytest.fixture
def noisy_spectrum():
    """Gaussian continuum + Hα emission line + noise, 800 pts, 4000–8000 Å."""
    rng = np.random.default_rng(42)
    wave = np.linspace(4000, 8000, 800) * u.Angstrom
    continuum = np.ones(800)
    line = 2.0 * np.exp(-0.5 * ((wave.value - 6563) / 5) ** 2)
    noise = rng.normal(0, 0.02, 800)
    flux = (continuum + line + noise) * u.mJy
    return Spectrum(spectral_axis=wave, flux=flux)


@pytest.fixture
def pixel_spectrum():
    """Uncalibrated spectrum with pixel axis (as produced by extract_spectrum)."""
    pixels = np.arange(200) * u.pixel
    flux = np.ones(200) * u.ct
    return Spectrum(spectral_axis=pixels, flux=flux)


# ---------------------------------------------------------------------------
# Synthetic 2D image fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def synthetic_image():
    """
    100×200 CCDData with a Gaussian spectral trace centred at row 50.
    """
    rng = np.random.default_rng(0)
    data = rng.normal(0, 5, (100, 200)).astype(np.float32)
    for row in range(100):
        data[row, :] += 1000 * np.exp(-0.5 * ((row - 50) / 4) ** 2)
    return CCDData(data, unit=u.Unit('adu'))


@pytest.fixture
def blank_image():
    """All-zeros 2×8 CCDData (same shape as the app's initial placeholder)."""
    return CCDData(np.zeros((2, 8)), unit=u.Unit('adu'))


# ---------------------------------------------------------------------------
# Temporary FITS / DAT file fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def spectrum_fits_file(tmp_path, flat_spectrum):
    """flat_spectrum written to a FITS file; returns the path."""
    path = str(tmp_path / "test_spectrum.fits")
    flat_spectrum.write(path, overwrite=True)
    return path


@pytest.fixture
def spectrum_dat_file(tmp_path, flat_spectrum):
    """flat_spectrum written to a 2-column DAT file; returns the path."""
    path = str(tmp_path / "test_spectrum.dat")
    data = np.column_stack([
        flat_spectrum.spectral_axis.value,
        flat_spectrum.flux.value,
    ])
    np.savetxt(path, data)
    return path


@pytest.fixture
def image_fits_file(tmp_path, synthetic_image):
    """synthetic_image written to a FITS file; returns the path."""
    path = str(tmp_path / "test_image.fits")
    synthetic_image.write(path, overwrite=True)
    return path
