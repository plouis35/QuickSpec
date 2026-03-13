"""
Image — thin controller wiring Image_Model <-> Image_View.

  - instantiate model and view
  - translate UI events (Clear button, sliders) into model/view calls
  - orchestrate load → display and reduce → display sequences
"""
import logging
import numpy as np

from tkinter import ttk
from astropy import units as u
from astropy.nddata import CCDData

from app.config import Config
from img.image_model import ImageModel
from img.image_view import ImageView


class Image:
    """Controller: coordinates ImageModel and ImageView."""

    def __init__(self, img_frame: ttk.Frame, bt_frame: ttk.Frame) -> None:
        """
        Args:
            img_frame (ttk.Frame): frame for the Matplotlib canvas
            bt_frame (ttk.Frame):  frame for buttons and sliders
        """
        self.conf = Config()
        self.model = ImageModel()
        self.view = ImageView(
            img_frame,
            bt_frame,
            on_clear=self.clear_image,
        )
        self.view.set_cut_callback(self._on_cut_changed)

        # show initial blank image
        self.view.show_image(self.model.img_stacked.data, show_colorbar=True)

    # ------------------------------------------------------------------
    # Public interface (called by Application / main.py)
    # ------------------------------------------------------------------

    @property
    def img_stacked(self) -> CCDData:
        """Current stacked/reduced image — accessed by Spectrum controller."""
        return self.model.img_stacked

    @property
    def img_axe(self):
        """Matplotlib Axes of the image panel — accessed by Spectrum controller."""
        return self.view.axe

    def clear_image(self) -> None:
        """Reset model state and clear the image panel."""
        self.model.reset()
        self.view.clear()
        self.view.show_image(self.model.img_stacked.data, show_colorbar=False)

    def load_images(self, paths: list[str]) -> bool:
        """
        Load, sum and display a list of FITS images.

        Args:
            paths (list[str]): FITS file paths

        Returns:
            bool: True on success
        """
        result = self.model.load(paths)
        if result is None:
            return False

        self.view.show_image(self.model.img_stacked.data, show_colorbar=True)
        self._refresh_cuts()
        return True

    def reduce_images(self) -> bool:
        """
        Reduce loaded images with calibration frames and redisplay.

        Returns:
            bool: True on success
        """
        result = self.model.reduce()
        if result is None:
            return False

        self.view.show_image(self.model.img_stacked.data, show_colorbar=True)
        self._refresh_cuts()
        return True

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _refresh_cuts(self) -> None:
        """Recompute and apply cut levels after a new image is displayed."""
        low_cut, high_cut, slider_min, slider_max = self.model.compute_cuts()
        self.view.set_slider_range(slider_min, slider_max, low_cut, high_cut)

    def _on_cut_changed(self, low_cut: float | None, high_cut: float | None) -> None:
        """
        Slider callback: apply new cut levels to the view.

        Args:
            low_cut (float | None):  new vmin, or None
            high_cut (float | None): new vmax, or None
        """
        self.view.apply_cut(low_cut, high_cut)
