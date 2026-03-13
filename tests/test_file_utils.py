"""
Tests for app/file_utils.py

Covers: classify_files() — routing DAT, 1D FITS, 2D FITS, bad files.
"""
import pytest
import numpy as np
from pathlib import Path
from astropy import units as u
from astropy.nddata import CCDData

from app.file_utils import classify_files


@pytest.fixture
def dat_file(tmp_path):
    path = tmp_path / "spec.dat"
    np.savetxt(str(path), np.column_stack([
        np.linspace(5000, 7000, 100),
        np.ones(100),
    ]))
    return str(path)


@pytest.fixture
def fits_2d_file(tmp_path, synthetic_image):
    path = str(tmp_path / "image.fits")
    synthetic_image.write(path, overwrite=True)
    return path


@pytest.fixture
def fits_1d_file(tmp_path, flat_spectrum):
    path = str(tmp_path / "spectrum.fits")
    flat_spectrum.write(path, overwrite=True)
    return path


class TestClassifyFiles:
    def test_dat_goes_to_spectra(self, dat_file):
        spectra, images = classify_files((dat_file,))
        assert dat_file in spectra
        assert dat_file not in images

    def test_2d_fits_goes_to_images(self, fits_2d_file):
        spectra, images = classify_files((fits_2d_file,))
        assert fits_2d_file in images
        assert fits_2d_file not in spectra

    def test_1d_fits_goes_to_spectra(self, fits_1d_file):
        spectra, images = classify_files((fits_1d_file,))
        assert fits_1d_file in spectra
        assert fits_1d_file not in images

    def test_mixed_files_split_correctly(self, dat_file, fits_2d_file, fits_1d_file):
        spectra, images = classify_files((dat_file, fits_2d_file, fits_1d_file))
        assert len(spectra) == 2   # dat + 1D fits
        assert len(images)  == 1   # 2D fits

    def test_nonexistent_file_is_skipped(self):
        spectra, images = classify_files(("/no/such/file.fits",))
        assert spectra == []
        assert images  == []

    def test_empty_input_returns_empty_lists(self):
        spectra, images = classify_files(())
        assert spectra == []
        assert images  == []
