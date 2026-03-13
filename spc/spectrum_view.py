"""
SpectrumView — all UI rendering logic for the 1D spectrum panel.

Depends on: Tkinter, Matplotlib.
"""
import logging

import numpy as np
import tkinter as tk
from tkinter import ttk

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.text import Annotation
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends._backend_tk import NavigationToolbar2Tk

from astropy import units as u
from specutils import Spectrum

from app.config import Config


class SpectrumView:
    """
    Manages all Tkinter/Matplotlib widgets for the 1D spectrum panel.

    Receives Spectrum objects to display; knows nothing about how
    they were computed.
    """

    def __init__(self, spc_frame: ttk.Frame, on_clear, on_show_lines, on_colorize) -> None:
        """
        Build the spectrum panel.

        Args:
            spc_frame (ttk.Frame): parent Tkinter frame
            on_clear (callable):       callback → reset button
            on_show_lines (callable):  callback → show lines button
            on_colorize (callable):    callback → colorize button
        """
        self.conf = Config()
        self.showed_lines: bool = False
        self.showed_colorized: bool = False
        self._colorize_collection = None   # LineCollection handle for fast removal
        self.colors = ('blue', 'red', 'green', 'orange', 'cyan')
        self.spectrum_color = 'grey'
        self._current_spectrum: Spectrum | None = None   # kept for color cycling only

        # theme-dependent line color
        theme = self.conf.get_str('display', 'theme')
        if theme == 'dark':
            self.lines_color = 'white'
        elif theme == 'light':
            self.lines_color = 'black'
        else:
            logging.warning(f"unsupported color theme: {theme!r}")
            self.lines_color = 'white'

        # Matplotlib figure + axes
        self.figure = Figure(figsize=(10, 3))
        self.axe: Axes = self.figure.add_subplot(111)
        self.axe.grid(color='grey', linestyle='--', linewidth=0.5)

        # Tkinter canvas
        self.canvas = FigureCanvasTkAgg(self.figure, spc_frame)
        self.canvas.draw()

        # Custom toolbar + extra buttons
        self.toolbar = _CustomSpcToolbar(self.canvas, spc_frame)

        ttk.Button(self.toolbar, text="Clear",      command=on_clear).pack(side=tk.LEFT, padx=5)
        ttk.Button(self.toolbar, text="Show lines", command=on_show_lines).pack(side=tk.LEFT, padx=5)
        ttk.Button(self.toolbar, text="Colorize",   command=on_colorize).pack(side=tk.LEFT, padx=5)

        self.toolbar.update()
        self.toolbar.pack(side=tk.TOP, fill=tk.X)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    # ------------------------------------------------------------------
    # Public draw methods
    # ------------------------------------------------------------------

    def clear(self) -> None:
        """Clear the spectrum axes and reset display state."""
        self.axe.clear()
        self.axe.grid(color='grey', linestyle='--', linewidth=0.5)
        self.showed_lines = False
        self.showed_colorized = False
        self._colorize_collection = None
        self._current_spectrum = None
        self.toolbar.update()
        self.figure.canvas.draw_idle()

    def show_spectrum(self, name: str, spectrum: Spectrum, calibrated: bool = False) -> None:
        """
        Plot a 1D spectrum on the axes.

        Args:
            name (str): legend label (underscores replaced to avoid matplotlib hiding it)
            spectrum (Spectrum): spectrum to plot
            calibrated (bool): True → wavelength axis, False → pixel axis
        """
        name = name.replace('_', '-')

        # axis labels & coord formatter
        if calibrated:
            self.axe.set_ylabel('Relative intensity')
            self.axe.set_xlabel('Wavelength (Angstrom)')
            self.axe.format_coord = lambda x, y: f"Lambda: ({x:.2f}, Intensity: {y:.2f})"
        else:
            self.axe.set_xlabel('Pixels')
            self.axe.set_ylabel('ADU')
            self.axe.format_coord = lambda x, y: f"Pixel: ({x:.0f}, ADU: {y:.0f})"

        # color cycling: grey for the "main" science spectrum, rotating colors for overlays
        if spectrum is self._current_spectrum:
            color = self.spectrum_color
        else:
            self.colors = tuple(self.colors[-1:]) + tuple(self.colors[:-1])
            color = self.colors[0]
            self._current_spectrum = spectrum

        self.axe.plot(spectrum.spectral_axis, spectrum.flux,
                      label=name, color=color, linewidth=0.8)
        self.axe.grid(color='grey', linestyle='--', linewidth=0.5)
        self.axe.legend()
        self.figure.canvas.draw_idle()

    def toggle_lines(self, lines: list[tuple[float, str]]) -> None:
        """
        Toggle spectral line markers on/off.

        Args:
            lines: list of (wavelength_angstrom, label) pairs to draw
        """
        if self.showed_lines:
            # remove existing markers
            for line in self.axe.lines:
                if line.get_color() == self.lines_color:
                    line.remove()
            for child in self.axe.get_children():
                if isinstance(child, Annotation):
                    child.remove()
            self.showed_lines = False
        else:
            xmin, xmax = self.axe.get_xbound()
            trans = self.axe.get_xaxis_transform()
            for wavelength, label in lines:
                if xmin < wavelength < xmax:
                    self.axe.axvline(wavelength, color=self.lines_color,
                                     lw=0.5, linestyle='--', alpha=0.8)
                    self.axe.annotate(label, xy=(wavelength, 1.05),
                                      xycoords=trans, fontsize=8,
                                      rotation=90, color=self.lines_color)
            self.showed_lines = True

        self.figure.canvas.draw_idle()

    def toggle_colorize(self, wavelengths, flux, min_wl: int, max_wl: int, bin_size: int) -> None:
        """
        Toggle the rainbow colorization under the spectrum curve.

        Uses a single LineCollection instead of N fill_between calls —
        typically 20× faster for spectra with hundreds of points.

        Args:
            wavelengths: array of wavelength values (Angstrom, no units)
            flux:        array of flux values (no units)
            min_wl (int): min wavelength for color mapping
            max_wl (int): max wavelength for color mapping
            bin_size (int): unused — kept for API compatibility
        """
        if self.showed_colorized:
            # Remove the existing LineCollection
            if self._colorize_collection is not None:
                self._colorize_collection.remove()
                self._colorize_collection = None
            self.showed_colorized = False
        else:
            y0 = flux.min()
            # Build one vertical segment per wavelength point: (x, y0) → (x, flux)
            segments = np.array([
                [[w, y0], [w, f]]
                for w, f in zip(wavelengths, flux)
            ])
            # Map each wavelength to a color via turbo colormap
            norm_wl = np.clip((wavelengths - min_wl) / (max_wl - min_wl), 0, 1)
            colors = plt.cm.turbo(norm_wl)

            lc = LineCollection(segments, colors=colors, linewidths=1.5, alpha=1.0)
            self._colorize_collection = self.axe.add_collection(lc)
            self.showed_colorized = True

        self.figure.canvas.draw_idle()

    def draw_trace_on_image(self, img_axe: Axes, trace_data) -> None:
        """
        Overlay the fitted trace on the 2D image axes.

        Args:
            img_axe (Axes): the 2D image matplotlib axes
            trace_data: trace array (e.g. FitTrace.trace)
        """
        # remove any existing trace lines
        for child in img_axe.get_children():
            if isinstance(child, Line2D):
                child.remove()
        img_axe.plot(trace_data, color='red', linestyle='dashed', linewidth=0.8)
        img_axe.get_figure().canvas.draw_idle()

    def draw_extraction_zones(
        self,
        img_axe: Axes,
        spectral_axis,
        trace,
        trace_y_size: int,
        sky_y_offset: int | None,
        sky_y_size: int | None,
    ) -> None:
        """
        Draw extraction and sky zone markers on the 2D image axes.

        Args:
            img_axe (Axes): the 2D image matplotlib axes
            spectral_axis: x values (pixel positions along dispersion)
            trace: y trace array
            trace_y_size (int): half-width of the extraction box
            sky_y_offset (int | None): sky band offset from trace
            sky_y_size (int | None): sky band half-width
        """
        kw = dict(linestyle='dashed', linewidth=0.5)
        img_axe.plot(spectral_axis, trace + trace_y_size, color='blue', **kw)
        img_axe.plot(spectral_axis, trace - trace_y_size, color='blue', **kw)

        if sky_y_offset is not None and sky_y_size is not None:
            for sign in (+1, -1):
                img_axe.plot(spectral_axis, trace + sign * sky_y_offset, color='green', **kw)
                img_axe.plot(spectral_axis, trace + sign * (sky_y_offset + sky_y_size),
                             color='green', **kw)

        img_axe.get_figure().canvas.draw_idle()


# ---------------------------------------------------------------------------
# Custom Matplotlib toolbar
# ---------------------------------------------------------------------------

class _CustomSpcToolbar(NavigationToolbar2Tk):
    toolitems = (
        ('Home',    'Reset zoom to original view',            'home_large',       'home'),
        ('Back',    'Back to previous view',                  'back_large',       'back'),
        ('Forward', 'Forward to next view',                   'forward_large',    'forward'),
        (None, None, None, None),
        ('Pan',     'Pan axes with left mouse, zoom with right', 'move_large',   'pan'),
        ('Zoom',    'Zoom to rectangle',                      'zoom_to_rect_large', 'zoom'),
        (None, None, None, None),
        ('Save',    'Save the figure',                        'filesave_large',   'save_figure'),
    )

    def __init__(self, canvas, parent) -> None:
        super().__init__(canvas=canvas, window=parent, pack_toolbar=True)
