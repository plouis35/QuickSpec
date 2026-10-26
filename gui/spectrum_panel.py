"""PySide6 and Matplotlib presentation for the 1D spectrum panel."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.collections import LineCollection
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.colors import to_rgb
from PySide6.QtWidgets import QPushButton, QVBoxLayout, QWidget
from specutils import Spectrum

from app.config import Config


class SpectrumPanel(QWidget):
    """Manage the Qt and Matplotlib widgets for the 1D spectrum panel."""

    def __init__(
        self, parent: QWidget, on_clear, on_show_lines, on_colorize
    ) -> None:
        super().__init__(parent)
        self.conf = Config()
        self.showed_lines = False
        self.showed_colorized = False
        self._colorize_collection: LineCollection | None = None
        self._line_markers = []
        self._line_annotations = []
        self.colors = ('blue', 'red', 'green', 'orange', 'cyan')
        self.spectrum_color = 'grey'
        self._current_spectrum: Spectrum | None = None

        self.figure = Figure(figsize=(10, 3))
        self.axe: Axes = self.figure.add_subplot(111)
        self.axe.grid(color='grey', linestyle='--', linewidth=0.5)
        self.canvas = FigureCanvasQTAgg(self.figure)
        #self.toolbar = NavigationToolbar2QT(self.canvas, self)
        self.toolbar = _CustomImgToolbar(self.canvas, self)
        for label, callback in (
            ("Clear", on_clear),
            ("Show lines", on_show_lines),
            ("Colorize", on_colorize),
        ):
            button = QPushButton(label, self)
            button.clicked.connect(callback)
            self.toolbar.addWidget(button)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas, 1)

    def clear(self) -> None:
        self.axe.clear()
        self.axe.grid(color='grey', linestyle='--', linewidth=0.5)
        self.showed_lines = False
        self.showed_colorized = False
        self._colorize_collection = None
        self._line_markers.clear()
        self._line_annotations.clear()
        self._current_spectrum = None
        self.toolbar.update()
        self.canvas.draw_idle()

    def show_spectrum(
        self, name: str, spectrum: Spectrum, calibrated: bool = False
    ) -> None:
        name = name.replace('_', '-')
        if calibrated:
            self.axe.set_ylabel('Relative intensity')
            self.axe.set_xlabel('Wavelength (Angstrom)')
            self.axe.format_coord = (
                lambda x, y: f"Lambda: ({x:.2f}, Intensity: {y:.2f})"
            )
        else:
            self.axe.set_xlabel('Pixels')
            self.axe.set_ylabel('ADU')
            self.axe.format_coord = lambda x, y: f"Pixel: ({x:.0f}, ADU: {y:.0f})"

        if spectrum is self._current_spectrum:
            color = self.spectrum_color
        else:
            self.colors = tuple(self.colors[-1:]) + tuple(self.colors[:-1])
            color = self.colors[0]
            self._current_spectrum = spectrum

        self.axe.plot(
            spectrum.spectral_axis,
            spectrum.flux,
            label=name,
            color=color,
            linewidth=0.8,
        )
        self.axe.grid(color='grey', linestyle='--', linewidth=0.5)
        self.axe.legend()
        self.canvas.draw_idle()

    def toggle_lines(self, lines: list[tuple[float, str]]) -> None:
        if self.showed_lines:
            for artist in self._line_markers + self._line_annotations:
                artist.remove()
            self._line_markers.clear()
            self._line_annotations.clear()
            self.showed_lines = False
        else:
            xmin, xmax = self.axe.get_xbound()
            transform = self.axe.get_xaxis_transform()
            line_color = self._contrast_color()
            for wavelength, label in lines:
                if xmin < wavelength < xmax:
                    marker = self.axe.axvline(
                        wavelength,
                        color=line_color,
                        lw=0.5,
                        linestyle='--',
                        alpha=0.8,
                    )
                    annotation = self.axe.annotate(
                        label,
                        xy=(wavelength, 1.05),
                        xycoords=transform,
                        fontsize=8,
                        rotation=90,
                        color=line_color,
                    )
                    self._line_markers.append(marker)
                    self._line_annotations.append(annotation)
            self.showed_lines = True
        self.canvas.draw_idle()

    def _contrast_color(self) -> str:
        red, green, blue = to_rgb(self.axe.get_facecolor())
        luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
        return 'black' if luminance > 0.5 else 'white'

    def toggle_colorize(
        self, wavelengths, flux, min_wl: int, max_wl: int, bin_size: int
    ) -> None:
        if self.showed_colorized:
            if self._colorize_collection is not None:
                self._colorize_collection.remove()
                self._colorize_collection = None
            self.showed_colorized = False
        else:
            y0 = flux.min()
            segments = np.array(
                [[[w, y0], [w, f]] for w, f in zip(wavelengths, flux)]
            )
            norm_wavelength = np.clip(
                (wavelengths - min_wl) / (max_wl - min_wl), 0, 1
            )
            colors = plt.cm.turbo(norm_wavelength)
            collection = LineCollection(
                segments, colors=colors, linewidths=1.5, alpha=1.0
            )
            self._colorize_collection = self.axe.add_collection(collection)
            self.showed_colorized = True
        self.canvas.draw_idle()

    def draw_trace_on_image(self, img_axe: Axes, trace_data) -> None:
        for child in img_axe.get_children():
            if isinstance(child, Line2D):
                child.remove()
        img_axe.plot(
            trace_data, color='red', linestyle='dashed', linewidth=0.8
        )
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
        line_options = dict(linestyle='dashed', linewidth=0.5)
        img_axe.plot(
            spectral_axis, trace + trace_y_size, color='blue', **line_options
        )
        img_axe.plot(
            spectral_axis, trace - trace_y_size, color='blue', **line_options
        )
        if sky_y_offset is not None and sky_y_size is not None:
            for sign in (+1, -1):
                img_axe.plot(
                    spectral_axis,
                    trace + sign * sky_y_offset,
                    color='green',
                    **line_options,
                )
                img_axe.plot(
                    spectral_axis,
                    trace + sign * (sky_y_offset + sky_y_size),
                    color='green',
                    **line_options,
                )
        img_axe.get_figure().canvas.draw_idle()


# ---------------------------------------------------------------------------
# Custom Matplotlib toolbar
# ---------------------------------------------------------------------------

class _CustomImgToolbar(NavigationToolbar2QT):
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
        super().__init__(canvas=canvas) #, window=parent, pack_toolbar=True)

