"""
Image — controller code for Image_Model <-> Image_View.

  - instantiate model and view
  - translate UI events (Clear button, sliders) into model calls
"""
import logging
import numpy as np

from tkinter import ttk
from astropy import units as u
from astropy.nddata import CCDData

from app.config import Config
from img.image_model import ImageModel
from img.image_view import ImageView
from img import sa100_geometry as sa100

class Image:
    """Controller: coordinates ImageModel and ImageView."""

    def __init__(self, img_frame: ttk.Frame, bt_frame: ttk.Frame, on_sa100_ready=None) -> None:
        """
        Args:
            img_frame (ttk.Frame): frame for the Matplotlib canvas
            bt_frame (ttk.Frame):  frame for buttons and sliders
            on_sa100_ready (callable): optional callback fired after SA100 selection
        """
        self.conf = Config()
        self.model = ImageModel()
        self._on_sa100_ready = on_sa100_ready
        self.view = ImageView(
            img_frame,
            bt_frame,
            on_clear=self.clear_image,
            on_sa100_points=self._on_sa100_points,
        )
        self.view.set_cut_callback(self._on_cut_changed)

        # show initial blank image
        self.view.show_image(self.model.img_stacked.data, show_colorbar=True)

    # ------------------------------------------------------------------
    # Public
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


    # ------------------------------------------------------------------
    # SA100 mode — toggle + 2-click capture
    # ------------------------------------------------------------------
    def _on_sa100_points(self, x0: float, y0: float,
                          x1: float, y1: float, x2: float, y2: float) -> None:
        """
        Process the SA100 rectification pipeline after the user
        has clicked order 0 (x0, y0) and a spectrum point (x1, y1) and a sky zone (sky_offset)

        Steps:
          1. Compute dispersion angle
          2. Rotate image
          3. Find where order 0 lands after rotation
          4. Crop spectrum and sky strips
          5. Display rectified image
          6. Overlay zone markers
          7. Extract spectrum
        """

        img = self.model.img_stacked
        if img is None or img.data.shape == (2, 8):
            logging.error("SA100: please load an image first")
            return

        # find angle
        angle = sa100.compute_angle(x0, y0, x1, y1)

        # rotate accordingly
        rectified = sa100.rectify(img, angle)

        # identify order 0 location after rotation
        h, w = rectified.data.shape
        cx, cy = w / 2.0, h / 2.0
        x0r, y0r = sa100.rotate_point(x0, y0, cx, cy, angle)

        # crop
        half_h    = self.conf.get_int('sa100', 'half_height')    or 5
        sky_off   = int(y2 - y1)
        sky_half  = self.conf.get_int('sa100', 'sky_half_height') or 10

        spc_strip, sky_strip = sa100.crop_spectrum(
            rectified, x0r, y0r,
            half_height=half_h,
            sky_offset=sky_off,
            sky_half_height=sky_half,
        )

        # display rectified image 
        self.view.axe.clear()
        self.view.show_image(rectified.data, show_colorbar=True)
        self._refresh_cuts()

        # display spectrum zones
        self.view.draw_sa100_zones(
            x0=x0r, y0=y0r,
            x1=w, y1=y0r,
            half_height=half_h,
            sky_offset=sky_off,
            sky_half_height=sky_half,
        )

        # process to spectrum 1D 
        if self._on_sa100_ready is not None:
            self._on_sa100_ready(spc_strip, sky_strip, x0r)
        else:
            logging.info("SA100: rectified image ready — connect on_sa100_ready to extract")


 
    # ------------------------------------------------------------------
    # Private
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
