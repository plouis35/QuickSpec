"""
ImageView — all UI rendering logic for the 2D image panel.

Depends on: Tkinter, Matplotlib.
"""
import logging

import numpy as np
import tkinter as tk
from tkinter import ttk

from matplotlib.image import AxesImage
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.colorbar import Colorbar
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends._backend_tk import NavigationToolbar2Tk

from app.config import Config


class ImageView:
    """
    Manages all Tkinter/Matplotlib widgets for the 2D image panel.

    Receives array data to display; knows nothing about how
    images were loaded or reduced.
    """

    def __init__(self, img_frame: ttk.Frame, bt_frame: ttk.Frame, on_clear) -> None:
        """
        Build the image panel (canvas, toolbar, sliders).

        Args:
            img_frame (ttk.Frame): frame for the Matplotlib canvas
            bt_frame (ttk.Frame):  frame for the control buttons / sliders
            on_clear (callable):   callback for the Clear button
        """
        self.conf = Config()
        self._colorbar: Colorbar | None = None
        self._axes_image: AxesImage | None = None
        self._cmap: str = 'grey'

        # Matplotlib figure + axes
        self.figure = Figure(figsize=(10, 3))
        self.axe: Axes = self.figure.add_subplot(111)

        # Tkinter canvas
        self.canvas = FigureCanvasTkAgg(self.figure, img_frame)
        self.canvas.draw()

        # Custom toolbar + Clear button + colormap selector
        self.toolbar = _CustomImgToolbar(self.canvas, img_frame)

        ttk.Button(self.toolbar, text="Clear", command=on_clear).pack(side=tk.LEFT, padx=5)

        _cmap_options = ["grey", "inferno", "nipy_spectral", "rainbow", "gnuplot"]
        _cmap_var = tk.StringVar(value=_cmap_options[0])
        ttk.OptionMenu(
            self.toolbar, _cmap_var, _cmap_options[0], *_cmap_options,
            command=self._on_cmap_change,
        ).pack(side=tk.LEFT, padx=5)

        self.toolbar.update()
        self.toolbar.pack(side=tk.TOP, fill=tk.X)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Cut-level sliders (packed into bt_frame, right-aligned)
        slider_frame = ttk.Frame(bt_frame)
        slider_frame.pack(side=tk.RIGHT, fill=tk.BOTH, padx=35, expand=True)

        # High cut slider
        high_frame = ttk.Frame(slider_frame)
        high_frame.pack(side=tk.TOP, fill=tk.X, padx=5)
        self._label_high = ttk.Label(high_frame, text="0")
        self._label_high.pack(side=tk.LEFT)
        self.slider_high = ttk.Scale(high_frame, from_=0, to=0, orient=tk.HORIZONTAL)
        self.slider_high.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.slider_high.bind("<ButtonRelease>", self._on_slider)

        # Low cut slider
        low_frame = ttk.Frame(slider_frame)
        low_frame.pack(side=tk.TOP, fill=tk.X, padx=5)
        self._label_low = ttk.Label(low_frame, text="0")
        self._label_low.pack(side=tk.LEFT)
        self.slider_low = ttk.Scale(low_frame, from_=0, to=0, orient=tk.HORIZONTAL)
        self.slider_low.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.slider_low.bind("<ButtonRelease>", self._on_slider)

        # Slider change callback (set by controller)
        self._on_cut_changed = None  # callable(low: float|None, high: float|None)

    def set_cut_callback(self, callback) -> None:
        """
        Register the callback that the controller wants called when sliders move.

        Args:
            callback (callable): fn(low: float | None, high: float | None)
        """
        self._on_cut_changed = callback

    # ------------------------------------------------------------------
    # Public draw methods
    # ------------------------------------------------------------------

    def show_image(self, data: np.ndarray, show_colorbar: bool = True) -> None:
        """
        Render a 2D array as an image on the axes.

        Args:
            data (np.ndarray): pixel data to display
            show_colorbar (bool): whether to show/refresh the colorbar
        """
        self.axe.axis('off')
        self.axe.set_yscale('linear')
        self.axe.format_coord = lambda x, y: f"(x, y): ({x:.0f}, {y:.0f})"

        self._axes_image = self.axe.imshow(
            data,
            origin='lower',
            interpolation='none',
            aspect='equal',
            cmap=self._cmap,
        )

        if show_colorbar:
            if self._colorbar is not None:
                try:
                    self._colorbar.remove()
                except Exception:
                    pass
            self._colorbar = self.figure.colorbar(
                self._axes_image, ax=self.axe, location='right', shrink=0.6
            )

        self.canvas.draw_idle()

    def clear(self) -> None:
        """Clear the image axes and reset toolbar state."""
        self.axe.clear()
        self.toolbar.update()
        self.canvas.draw_idle()

    def set_slider_range(self, slider_min: float, slider_max: float,
                         low_cut: float, high_cut: float) -> None:
        """
        Configure slider range and initial positions.

        Args:
            slider_min (float): minimum slider value
            slider_max (float): maximum slider value
            low_cut (float):  initial low cut position
            high_cut (float): initial high cut position
        """
        for slider in (self.slider_low, self.slider_high):
            slider.config(from_=slider_min, to=slider_max)
        self._apply_cut(low_cut=low_cut, high_cut=high_cut)

    def apply_cut(self, low_cut: float | None, high_cut: float | None) -> None:
        """
        Update the image display cut levels (public entry point).

        Args:
            low_cut (float | None):  new vmin, or None to leave unchanged
            high_cut (float | None): new vmax, or None to leave unchanged
        """
        self._apply_cut(low_cut, high_cut)
        self.canvas.draw_idle()

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _apply_cut(self, low_cut: float | None, high_cut: float | None) -> None:
        if self._axes_image is None:
            return
        try:
            if low_cut is not None:
                self._axes_image.norm.vmin = low_cut
                self._label_low.config(text=int(low_cut))
                self.slider_low.set(low_cut)
            if high_cut is not None:
                self._axes_image.norm.vmax = high_cut
                self._label_high.config(text=int(high_cut))
                self.slider_high.set(high_cut)
        except Exception as e:
            logging.error(f"cut level error: {e}")

    def _on_slider(self, event) -> None:
        """Translate a slider release event into a cut callback."""
        if self._on_cut_changed is None:
            return
        slider: ttk.Scale = event.widget
        if slider is self.slider_high:
            self._on_cut_changed(None, slider.get())
        elif slider is self.slider_low:
            self._on_cut_changed(slider.get(), None)

    def _on_cmap_change(self, cmap: str) -> None:
        """Apply a colormap change immediately."""
        self._cmap = cmap
        if self._axes_image is not None:
            self._axes_image.set_cmap(cmap)
            self.figure.canvas.draw()
        logging.info(f"colormap changed to {cmap}")


# ---------------------------------------------------------------------------
# Custom Matplotlib toolbar
# ---------------------------------------------------------------------------

class _CustomImgToolbar(NavigationToolbar2Tk):
    toolitems = (
        ('Home',    'Reset zoom to original view',               'home_large',       'home'),
        ('Back',    'Back to previous view',                     'back_large',       'back'),
        ('Forward', 'Forward to next view',                      'forward_large',    'forward'),
        (None, None, None, None),
        ('Pan',     'Pan axes with left mouse, zoom with right', 'move_large',       'pan'),
        ('Zoom',    'Zoom to rectangle',                         'zoom_to_rect_large', 'zoom'),
        (None, None, None, None),
        ('Save',    'Save the figure',                           'filesave_large',   'save_figure'),
    )

    def __init__(self, canvas, parent) -> None:
        super().__init__(canvas=canvas, window=parent, pack_toolbar=True)
