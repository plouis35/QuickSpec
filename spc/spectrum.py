"""
Spectrum — thin controller wiring Spectrum_Model <-> Spectrum_View.

  - instantiate model and view
  - translate UI events (button clicks) into model calls
  - pass model results to the view for display
  - read config values needed to bridge model ↔ view

"""
import logging
from pathlib import Path

from tkinter import ttk
from matplotlib.axes import Axes
from astropy.nddata import CCDData

from app.config import Config
from spc.spectrum_model import SpectrumModel
from spc.spectrum_view import SpectrumView


class Spectrum:
    """Controller: coordinates SpectrumModel and SpectrumView."""

    def __init__(self, spc_frame: ttk.Frame, img_axe: Axes) -> None:
        """
        Args:
            spc_frame (ttk.Frame): parent frame for the spectrum panel
            img_axe (Axes): 2D image axes (for trace/zone overlays)
        """
        self.conf = Config()
        self.img_axe = img_axe

        self.model = SpectrumModel()
        self.view = SpectrumView(
            spc_frame,
            on_clear=self.cb_clear,
            on_show_lines=self.cb_show_lines,
            on_colorize=self.cb_colorize,
        )

    # ------------------------------------------------------------------
    # Public interface (called by Application / main.py)
    # ------------------------------------------------------------------

    def open_spectrum(self, path: str) -> bool:
        """Load and display a 1D spectrum from a file."""
        spec = self.model.load_spectrum(path)
        if spec is None:
            return False
        self.view.show_spectrum(name=Path(path).stem, spectrum=spec, calibrated=True)
        return True

    def do_trace(self, img_stacked: CCDData) -> bool:
        """Find and display the spectral trace on the 2D image."""
        trace = self.model.trace(img_stacked)
        if trace is None:
            return False
        self.view.draw_trace_on_image(self.img_axe, trace.trace)
        return True

    def do_extract(self, img_stacked: CCDData) -> bool:
        """Extract 1D spectrum and display extraction zones."""
        spec = self.model.extract(img_stacked)
        if spec is None:
            return False

        trace_y_size = self.conf.get_int('processing', 'trace_y_size') or 15
        sky_sub = self.conf.get_bool('processing', 'sky_substract')
        if sky_sub in (None, True):
            sky_offset = self.conf.get_int('processing', 'sky_y_offset')
            sky_size   = self.conf.get_int('processing', 'sky_y_size')
        else:
            sky_offset = sky_size = None

        self.view.draw_extraction_zones(
            img_axe=self.img_axe,
            spectral_axis=spec.spectral_axis,
            trace=self.model.science_trace.trace,
            trace_y_size=trace_y_size,
            sky_y_offset=sky_offset,
            sky_y_size=sky_size,
        )
        self.view.clear()
        self.view.show_spectrum(
            name=SpectrumModel.get_object_name(img_stacked.header),
            spectrum=spec,
            calibrated=False,
        )
        return True

    def do_calibrate(self, img_stacked: CCDData) -> bool:
        """Wavelength-calibrate the spectrum and redisplay."""
        spec = self.model.calibrate(img_stacked)
        if spec is None:
            return False
        self.view.clear()
        self.view.show_spectrum(
            name=SpectrumModel.get_object_name(img_stacked.header),
            spectrum=spec,
            calibrated=True,
        )
        return True

    def do_response(self, img_stacked: CCDData) -> bool:
        """Apply instrumental response and redisplay."""
        spec = self.model.apply_response()
        if spec is None:
            return False
        self.view.clear()
        self.view.show_spectrum(
            name=SpectrumModel.get_object_name(img_stacked.header),
            spectrum=spec,
            calibrated=True,
        )
        return True

    def do_smooth(self, img_stacked: CCDData) -> bool:
        """Apply smooth/crop/normalize and redisplay."""
        spec = self.model.smooth()
        if spec is None:
            return False
        self.view.clear()
        self.view.show_spectrum(
            name=SpectrumModel.get_object_name(img_stacked.header),
            spectrum=spec,
            calibrated=True,
        )
        return True

    # ------------------------------------------------------------------
    # UI callbacks (wired to SpectrumView buttons)
    # ------------------------------------------------------------------

    def cb_clear(self) -> None:
        """Reset model state and clear the view."""
        self.model.reset()
        self.view.clear()

    def cb_show_lines(self) -> None:
        """Toggle spectral line markers."""
        if self.model.science_spectrum is None:
            logging.info("please calibrate first")
            return
        if not self.model.is_calibrated():
            logging.warning("spectrum needs to be calibrated first")
            return
        lines = [
            (float(wave) * 10, label)
            for wave, label in self.conf.config.items('lines')
        ]
        self.view.toggle_lines(lines)

    def cb_colorize(self) -> None:
        """Toggle rainbow colorization under the spectrum."""
        if self.model.science_spectrum is None:
            logging.info("please calibrate first")
            return
        if not self.model.is_calibrated():
            logging.warning("spectrum needs to be calibrated first")
            return

        spec     = self.model.science_spectrum
        min_wl   = self.conf.get_int('display', 'colorize_min_wl')   or 3800
        max_wl   = self.conf.get_int('display', 'colorize_max_wl')   or 7500
        bin_size = self.conf.get_int('display', 'colorize_bin_size') or 20

        self.view.toggle_colorize(
            wavelengths=spec.wavelength.value,
            flux=spec.flux.value,
            min_wl=min_wl,
            max_wl=max_wl,
            bin_size=bin_size,
        )
