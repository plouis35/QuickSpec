"""
file_utils — FITS/DAT file inspection helpers.
"""
import logging
from pathlib import Path


from astropy import units as u
from astropy.nddata import CCDData
from astropy.io import fits

def classify_files(paths: tuple[str, ...]) -> tuple[list[str], list[str]]:
    """
    Classify a list of selected file paths into 1D spectra and 2D images.

    Args:
        paths: tuple of file paths as returned by askopenfilenames()

    Returns:
        tuple:
            - spectrum_paths (list[str]): DAT files and 1D FITS spectra
            - image_paths (list[str]):    2D FITS images
    """
    spectrum_paths: list[str] = []
    image_paths: list[str] = []

    for path in paths:
        suffix = Path(path).suffix.lower()

        if suffix == '.dat':
            logging.info(f"{path} is a spectrum (2-column DAT)")
            spectrum_paths.append(path)
            continue

        # FITS file — inspect dimensionality
        try:
            hdu_count = len(fits.open(path))
            fit_data = CCDData.read(
                path,
                unit=u.dimensionless_unscaled,
                hdu=hdu_count - 1,
            )
            if fit_data.ndim == 1:
                logging.info(f"{path} is a spectrum (naxis=1, shape={fit_data.shape})")
                spectrum_paths.append(path)
            elif fit_data.ndim == 2:
                logging.info(f"{path} is a 2D image (naxis=2, shape={fit_data.shape})")
                image_paths.append(path)
            else:
                logging.error(f"{path} is not a supported format (naxis > 2)")

        except Exception as e:
            logging.error(f"{path}: {e}")

    return spectrum_paths, image_paths
