"""
ImageModel — processing modules for 2D spectral images.
"""
import logging
import warnings
import numpy as np

from astropy import units as u
from astropy.nddata import CCDData
from astropy.utils.exceptions import AstropyWarning

from app.config import Config
from img.img_utils import Images, ImagesCombiner

warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', UserWarning)

# Empty placeholder shape used before any image is loaded
_INIT_SHAPE: tuple[int, int] = (2, 8)


class ImageModel:
    """
    Holds the state of the current 2D image pipeline and exposes
    one method per processing step.

    Attributes:
        img_stacked:  the current combined/reduced CCDData frame
        img_combiner: the ImagesCombiner holding all loaded frames
        img_reduced:  True once reduce() has been called successfully
    """

    def __init__(self) -> None:
        self.conf = Config()
        self.img_stacked: CCDData = CCDData(np.zeros(_INIT_SHAPE), unit=u.Unit('adu'))
        self.img_combiner: ImagesCombiner | None = None
        self.img_reduced: bool = False

    # ------------------------------------------------------------------
    # State management
    # ------------------------------------------------------------------

    def reset(self) -> None:
        """Clear all image state back to initial placeholder."""
        self.img_stacked = CCDData(np.zeros(_INIT_SHAPE), unit=u.Unit('adu'))
        self.img_combiner = None
        self.img_reduced = False

    def stats(self) -> tuple[float, float, float, float]:
        """
        Compute basic statistics on the current stacked image.

        Returns:
            tuple: (std, mean, min, max)
        """
        data = self.img_stacked.data
        return data.std(), data.mean(), data.min(), data.max()

    def compute_cuts(self) -> tuple[float, float, float, float]:
        """
        Compute display cut levels and slider range from image statistics.

        Returns:
            tuple: (low_cut, high_cut, slider_min, slider_max)
        """
        v_std, v_mean, _, _ = self.stats()
        nb_sigma = self.conf.get_int('display', 'contrast_level') or 6

        low_cut  = v_mean - nb_sigma * v_std
        high_cut = v_mean + nb_sigma * v_std
        slider_min = low_cut  - 1 * nb_sigma * v_std
        slider_max = high_cut + 8 * nb_sigma * v_std
        return low_cut, high_cut, slider_min, slider_max

    # ------------------------------------------------------------------
    # Processing steps
    # ------------------------------------------------------------------

    def load(self, paths: list[str]) -> CCDData | None:
        """
        Load and sum a list of FITS images.

        Args:
            paths (list[str]): FITS file paths to load

        Returns:
            CCDData | None: summed image, also stored in self.img_stacked
        """
        max_memory = self.conf.get_float('pre_processing', 'max_memory') or 1e9
        y_crop     = self.conf.get_str('pre_processing', 'y_crop')

        try:
            combiner = Images.from_fits_by_names_list(
                imgs=paths, max_memory=max_memory
            ).y_crop(y_crop=y_crop)
            stacked = combiner.sum()
        except Exception as e:
            logging.error(f"load failed: {e}")
            return None

        self.img_stacked  = stacked.copy()
        self.img_combiner = combiner
        self.img_reduced  = False

        v_std, v_mean, v_min, v_max = self.stats()
        logging.info(f"image stats: min={v_min}, max={v_max}, mean={v_mean}, std={v_std}")
        return self.img_stacked

    def reduce(self) -> CCDData | None:
        """
        Reduce the loaded images using master calibration frames (bias/dark/flat).

        Returns:
            CCDData | None: reduced and combined image, or None on error
        """
        if self.img_combiner is None:
            logging.error("please load some images before reducing")
            return None

        if self.img_reduced:
            logging.warning("reduce already done — skipped")
            return self.img_stacked

        try:
            result = self.img_combiner.reduce_images()
        except Exception as e:
            logging.error(f"reduce failed: {e}")
            return None

        if result is None:
            logging.error("unable to reduce images")
            return None

        self.img_stacked = result.copy()
        self.img_reduced = True

        v_std, v_mean, v_min, v_max = self.stats()
        logging.info(f"image stats: min={v_min}, max={v_max}, mean={v_mean}, std={v_std}")
        return self.img_stacked
