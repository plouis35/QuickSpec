"""PySide6 and Matplotlib presentation for the 2D image panel."""
import logging

import numpy as np
from matplotlib.axes import Axes
from matplotlib.colorbar import Colorbar
from matplotlib.figure import Figure
from matplotlib.image import AxesImage
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from PySide6.QtCore import QSignalBlocker, Qt
from PySide6.QtWidgets import (
    QComboBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSlider,
    QVBoxLayout,
    QWidget,
)
from skimage.measure import centroid

from app.config import Config


class ImagePanel(QWidget):
    """Manage the Qt and Matplotlib widgets for the 2D image panel."""

    def __init__(self, parent: QWidget, on_clear, on_sa100_points) -> None:
        super().__init__(parent)
        self.conf = Config()
        self._colorbar: Colorbar | None = None
        self._axes_image: AxesImage | None = None
        self._cmap = 'grey'
        self._on_sa100_points = on_sa100_points
        self._sa100_active = False
        self._sa100_clicks: list[tuple[float, float]] = []
        self._sa100_markers = []
        self._click_cid = None
        self._on_cut_changed = None
        self._slider_min = 0.0
        self._slider_max = 1.0

        self.figure = Figure(figsize=(10, 3))
        self.axe: Axes = self.figure.add_subplot(111)
        self.canvas = FigureCanvasQTAgg(self.figure)
        #self.toolbar = NavigationToolbar2QT(self.canvas, self)
        self.toolbar = _CustomImgToolbar(self.canvas, self)
        self._add_toolbar_controls(on_clear)

        self.slider_low, self._label_low = self._create_cut_slider()
        self.slider_high, self._label_high = self._create_cut_slider()
        self.cut_controls = QWidget(self)
        low_controls = QHBoxLayout()
        low_controls.setContentsMargins(0, 0, 0, 0)
        low_controls.addWidget(QLabel("Low"))
        low_controls.addWidget(self.slider_low, 1)
        low_controls.addWidget(self._label_low)
        high_controls = QHBoxLayout()
        high_controls.setContentsMargins(0, 0, 0, 0)
        high_controls.addWidget(QLabel("High"))
        high_controls.addWidget(self.slider_high, 1)
        high_controls.addWidget(self._label_high)
        cut_layout = QVBoxLayout(self.cut_controls)
        cut_layout.setContentsMargins(0, 0, 0, 0)
        cut_layout.addWidget(self.toolbar)
        cut_layout.addLayout(high_controls, 1)
        cut_layout.addLayout(low_controls, 1)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(4)
        #layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas, 1)

        self.slider_low.valueChanged.connect(
            lambda value: self._slider_moved(self.slider_low, value)
        )
        self.slider_high.valueChanged.connect(
            lambda value: self._slider_moved(self.slider_high, value)
        )

    def _add_toolbar_controls(self, on_clear) -> None:
        controls = QWidget(self)
        layout = QHBoxLayout(controls)
        layout.setContentsMargins(4, 0, 4, 0)
        clear_button = QPushButton("Clear", controls)
        clear_button.clicked.connect(on_clear)
        layout.addWidget(clear_button)
        self._cmap_selector = QComboBox(controls)
        self._cmap_selector.addItems(
            ["Grey", "Inferno", "Nipy_spectral", "Rainbow", "Gnuplot"]
        )
        self._cmap_selector.setMaximumWidth(150)
        self._cmap_selector.currentTextChanged.connect(self._on_cmap_change)
        layout.addWidget(self._cmap_selector)
        self.toolbar.addWidget(controls)

    def _create_cut_slider(self) -> tuple[QSlider, QLabel]:
        slider = QSlider(Qt.Orientation.Horizontal, self)
        slider.setRange(0, 10000)
        label = QLabel("0", self)
        return slider, label

    def set_cut_callback(self, callback) -> None:
        self._on_cut_changed = callback

    def show_image(self, data: np.ndarray, show_colorbar: bool = True) -> None:
        self.axe.axis('off')
        self.axe.set_yscale('linear')
        self.axe.format_coord = lambda x, y: f"(x, y): ({x:.0f}, {y:.0f})"
        self._axes_image = self.axe.imshow(
            data,
            origin='lower',
            interpolation='none',
            aspect='equal',
            cmap=self._cmap.lower(),
        )
        if show_colorbar:
            if self._colorbar is not None:
                self._colorbar.remove()
            self._colorbar = self.figure.colorbar(
                self._axes_image, ax=self.axe, location='right', shrink=0.6
            )
        self.canvas.draw_idle()

    def clear(self) -> None:
        self._cancel_sa100()
        if self._colorbar is not None:
            self._colorbar.remove()
            self._colorbar = None
        self.axe.clear()
        self._axes_image = None
        self.toolbar.update()
        self.canvas.draw_idle()

    def set_slider_range(
        self, slider_min: float, slider_max: float, low_cut: float, high_cut: float
    ) -> None:
        if slider_max <= slider_min:
            padding = max(abs(slider_min) * 0.01, 1.0)
            slider_min -= padding
            slider_max += padding
        self._slider_min = slider_min
        self._slider_max = slider_max
        self._apply_cut(low_cut, high_cut)

    def apply_cut(self, low_cut: float | None, high_cut: float | None) -> None:
        self._apply_cut(low_cut, high_cut)
        self.canvas.draw_idle()

    def draw_sa100_zones(
        self,
        x0: float,
        y0: float,
        x1: float,
        y1: float,
        half_height: int,
        sky_offset: int,
        sky_half_height: int,
    ) -> None:
        for artist in getattr(self, '_zone_artists', []):
            artist.remove()
        self._zone_artists = []
        line_options = dict(linewidth=0.8, linestyle='--')
        for offset in (+half_height, -half_height):
            self._zone_artists.append(
                self.axe.axhline(y0 + offset, color='red', **line_options)
            )
        for offset in (
            sky_offset + sky_half_height,
            sky_offset - sky_half_height,
        ):
            self._zone_artists.append(
                self.axe.axhline(y0 + offset, color='blue', **line_options)
            )
        self._zone_artists.append(
            self.axe.axvline(x0, color='yellow', linewidth=0.8, linestyle=':')
        )
        self.canvas.draw_idle()

    def toggle_sa100(self) -> None:
        if self._sa100_active:
            self._cancel_sa100()
        else:
            self._start_sa100()

    def _start_sa100(self) -> None:
        self._sa100_active = True
        self._sa100_clicks = []
        self.canvas.setCursor(Qt.CursorShape.CrossCursor)
        self._click_cid = self.canvas.mpl_connect(
            'button_press_event', self._on_canvas_click
        )
        logging.info(
            "SA100 mode ON — click 0 order, then a point on the spectrum, "
            "and finally a point for the sky"
        )

    def _cancel_sa100(self) -> None:
        self._sa100_active = False
        self._sa100_clicks = []
        self.canvas.unsetCursor()
        if self._click_cid is not None:
            self.canvas.mpl_disconnect(self._click_cid)
            self._click_cid = None
        for artist in self._sa100_markers:
            artist.remove()
        self._sa100_markers.clear()
        self.canvas.draw_idle()

    def _on_canvas_click(self, event) -> None:
        search_half_size = 20
        if event.inaxes != self.axe or not event.dblclick:
            return
        if event.xdata is None or event.ydata is None or self._axes_image is None:
            return

        x, y = int(event.xdata), int(event.ydata)
        raw_data = np.asarray(self._axes_image.get_array())
        y_min = max(0, y - search_half_size)
        y_max = min(raw_data.shape[0], y + search_half_size)
        x_min = max(0, x - search_half_size)
        x_max = min(raw_data.shape[1], x + search_half_size)
        crop = raw_data[y_min:y_max, x_min:x_max]
        if crop.size == 0:
            return
        image_subtracted = np.maximum(0, crop - np.median(crop))
        click_index = len(self._sa100_clicks)

        if click_index < 2:
            centroid_y, centroid_x = centroid(image_subtracted)
            x, y = x_min + centroid_x, y_min + centroid_y
            color = 'r+' if click_index == 0 else 'b+'
        else:
            color = 'g+'

        self._sa100_clicks.append((x, y))
        marker, = self.axe.plot(
            x, y, color, markersize=18, markeredgewidth=2
        )
        self._sa100_markers.append(marker)
        self.canvas.draw_idle()
        logging.info(f"SA100: point {click_index + 1} at ({x:.1f}, {y:.1f})")

        if len(self._sa100_clicks) == 3:
            points = self._sa100_clicks.copy()
            self._cancel_sa100()
            if self._on_sa100_points is not None:
                self._on_sa100_points(*points[0], *points[1], *points[2])

    def _apply_cut(
        self, low_cut: float | None, high_cut: float | None
    ) -> None:
        if self._axes_image is None:
            return
        if low_cut is not None:
            self._axes_image.norm.vmin = low_cut
            self._label_low.setText(str(int(low_cut)))
            with QSignalBlocker(self.slider_low):
                self.slider_low.setValue(self._slider_position(low_cut))
        if high_cut is not None:
            self._axes_image.norm.vmax = high_cut
            self._label_high.setText(str(int(high_cut)))
            with QSignalBlocker(self.slider_high):
                self.slider_high.setValue(self._slider_position(high_cut))

    def _slider_position(self, value: float) -> int:
        span = self._slider_max - self._slider_min
        if span <= 0:
            return 0
        return round((value - self._slider_min) / span * 10000)

    def _slider_moved(self, slider: QSlider, position: int) -> None:
        if self._on_cut_changed is None:
            return
        
        value = self._slider_min + position / 10000 * (
            self._slider_max - self._slider_min
        )

        if slider is self.slider_low:
            self._on_cut_changed(value, None)
        else:
            self._on_cut_changed(None, value)

    def _on_cmap_change(self, cmap: str) -> None:
        self._cmap = cmap
        if self._axes_image is not None:
            self._axes_image.set_cmap(cmap.lower())
            self.canvas.draw_idle()
        logging.info(f"colormap changed to {cmap}")

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

