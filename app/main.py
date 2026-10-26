"""
main application class:
- initialize logging and configuration
- create the GUI
- instantiate image and spectrum processors
- monitor configuration changes and new FITS files
"""
import logging
from pathlib import Path
import functools

import matplotlib.pyplot as plt
from PySide6.QtCore import QObject, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QKeySequence, QPalette
from PySide6.QtWidgets import (
    QApplication,
    QFileDialog,
    QHBoxLayout,
    QMainWindow,
    QMenu,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from app.logger import LogHandler
import app.os_utils as os_utils
from app.file_utils import classify_files
from app.config import Config
from app.file_monitor import FileMonitor
from gui.image_panel import ImagePanel
from gui.spectrum_panel import SpectrumPanel
from processing.image_processor import ImageProcessor
from processing.spectrum_processor import SpectrumProcessor
from img import sa100_geometry as sa100


class _MonitorSignals(QObject):
    new_fits = Signal(str)
    config_changed = Signal()


class Application(QMainWindow):

    def __init__(self, app_name: str, app_version: str) -> None:
        """
        initialize GUI components
        start watchdog timer

        Args:
            app_name (str): app name 
            app_version (str): app version
        """        
        super().__init__()
        self.conf: Config = Config()
        LogHandler().initialize()

        self._apply_theme()
        plt.rcParams['figure.constrained_layout.use'] = True
        
        self.setWindowTitle(f"{app_name} v{app_version}")
        self.resize(1100, 800)
        self.setMinimumSize(800, 600)
        self.app_name = app_name
        self.app_version = app_version

        self._monitor_signals = _MonitorSignals(self)
        self._monitor_signals.new_fits.connect(self._on_new_fits_file)
        self._monitor_signals.config_changed.connect(self._on_config_file_changed)
        self.file_monitor = FileMonitor(
            on_new_fits=self._monitor_signals.new_fits.emit,
            on_config_changed=self._monitor_signals.config_changed.emit,
        )
        self._config_reload_timer = QTimer(self)
        self._config_reload_timer.setSingleShot(True)
        self._config_reload_timer.setInterval(1000)
        self._config_reload_timer.timeout.connect(self._reload_config)

        # create image and spectrum panels
        self.create_panels()
        self.create_menus()

        self._start_file_monitor()

        # and display major packages versions installed
        logging.info(f"{app_name} v{app_version} started")
        os_utils.show_versions()

    def _apply_theme(self) -> None:
        app = QApplication.instance()
        if app is None:
            raise RuntimeError("QApplication must be created before Application")

        theme = self.conf.get_str('display', 'theme')
        if theme not in ('dark', 'light'):
            if theme is not None:
                logging.error(f"unsupported color theme: {theme!r}")
            theme = 'dark'

        app.setStyle("Fusion")
        if theme == 'dark':
            plt.style.use('dark_background')
            palette = QPalette()
            palette.setColor(QPalette.ColorRole.Window, QColor(53, 53, 53))
            palette.setColor(QPalette.ColorRole.WindowText, QColor(240, 240, 240))
            palette.setColor(QPalette.ColorRole.Base, QColor(35, 35, 35))
            palette.setColor(QPalette.ColorRole.AlternateBase, QColor(53, 53, 53))
            palette.setColor(QPalette.ColorRole.Text, QColor(240, 240, 240))
            palette.setColor(QPalette.ColorRole.Button, QColor(53, 53, 53))
            palette.setColor(QPalette.ColorRole.ButtonText, QColor(240, 240, 240))
            palette.setColor(QPalette.ColorRole.Highlight, QColor(42, 130, 218))
            palette.setColor(
                QPalette.ColorRole.HighlightedText, QColor(255, 255, 255)
            )
        elif theme == 'light':
            plt.style.use('grayscale')
            palette = QPalette()
            palette.setColor(QPalette.ColorRole.Window, QColor(240, 240, 240))
            palette.setColor(QPalette.ColorRole.WindowText, QColor(20, 20, 20))
            palette.setColor(QPalette.ColorRole.Base, QColor(255, 255, 255))
            palette.setColor(QPalette.ColorRole.AlternateBase, QColor(245, 245, 245))
            palette.setColor(QPalette.ColorRole.Text, QColor(20, 20, 20))
            palette.setColor(QPalette.ColorRole.Button, QColor(240, 240, 240))
            palette.setColor(QPalette.ColorRole.ButtonText, QColor(20, 20, 20))
            palette.setColor(QPalette.ColorRole.Highlight, QColor(42, 130, 218))
            palette.setColor(
                QPalette.ColorRole.HighlightedText, QColor(255, 255, 255)
            )

        app.setPalette(palette)
        self._theme = theme

    def create_panels(self) -> None:
        """
        Create the shared image and spectrum display panels.
        """
        central_widget = QWidget(self)
        layout = QVBoxLayout(central_widget)
        layout.setContentsMargins(6, 6, 6, 6)
        self.controls_layout = QHBoxLayout()
        layout.addLayout(self.controls_layout)

        self.image_processor = ImageProcessor()
        self.spectrum_processor = SpectrumProcessor()
        self.image_panel = ImagePanel(
            central_widget,
            on_clear=self.cb_clear_image,
            on_sa100_points=self._on_sa100_points,
        )
        self.image_panel.set_cut_callback(self._on_cut_changed)
        self.image_panel.show_image(
            self.image_processor.img_stacked.data, show_colorbar=True
        )
        self.controls_layout.addWidget(self.image_panel.cut_controls, 1)
        self.spectrum_panel = SpectrumPanel(
            central_widget,
            on_clear=self.cb_clear_spectrum,
            on_show_lines=self.cb_show_lines,
            on_colorize=self.cb_colorize,
        )
        splitter = QSplitter(Qt.Orientation.Vertical, central_widget)
        splitter.addWidget(self.image_panel)
        splitter.addWidget(self.spectrum_panel)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 1)
        layout.addWidget(splitter, 1)
        self.setCentralWidget(central_widget)

    def create_menus(self) -> None:
        """
        Create the application menu bar.
        """
        menu_bar = self.menuBar()
        file_menu = menu_bar.addMenu("File")
        load_action = file_menu.addAction("Load", self.cb_open_files)
        load_action.setShortcut(QKeySequence("Ctrl+O"))
        file_menu.addSeparator()
        exit_action = file_menu.addAction(
            "Exit", lambda _checked=False: self.close()
        )
        exit_action.setShortcut(QKeySequence("Ctrl+X"))

        def run_step(name: str, action) -> None:
            logging.info(f"step '{name}' started...")
            action()

        def add_entries(menu: QMenu, entries: list) -> list:
            actions = []
            for entry in entries:
                if entry is None:
                    menu.addSeparator()
                else:
                    label, action = entry
                    actions.append(menu.addAction(
                        label,
                        lambda checked=False, n=label, a=action: run_step(n, a),
                    ))
            return actions

        longslit_menu = menu_bar.addMenu("Longslit")
        longslit_actions = add_entries(longslit_menu, [
            ("Process auto",              self.cb_run_all),
            None,
            ("Find spectrum",             self.cb_trace_spectrum),
            ("Extract spectrum",          self.cb_extract_spectrum),
            ("Calibrate spectrum",        self.cb_calibrate_spectrum),
            ("Apply response",            self.cb_apply_response),
            ("Smooth, crop & normalize",  self.cb_smooth_spectrum),
        ])
        longslit_actions[0].setShortcut(QKeySequence("Ctrl+P"))

        slitless_menu = menu_bar.addMenu("Slitless")
        slitless_actions = add_entries(slitless_menu, [
            ("3-click process",           self.cb_slitless),
        ])
        slitless_actions[0].setShortcut(QKeySequence("Ctrl+Q"))
        
    def set_cursor(self, cursor: str = '') -> None:
        """
        Set the busy cursor while a synchronous processing action is running.

        Args:
            cursor: non-empty to show the busy cursor; empty to restore default
        """        
        if cursor:
            self.setCursor(Qt.CursorShape.BusyCursor)# .WaitCursor)
        else:
            self.unsetCursor()
        QApplication.processEvents()

    def set_title(self, title: str = '') -> None:
        """
        Set the title while retaining the application name and version.

        Args:
            title: additional title text.
        """        
        self.setWindowTitle(f"{self.app_name} v{self.app_version} - {title}")
        
    # local callbacks for buttons
    @staticmethod
    def run_long_operation(func):
        """
        decorator for callback buttons
        display 'waiting' cursor while processing
        (does not work well on Windows platforms ...)
        """        
        @functools.wraps(func)
        def wrap(self, *args, **kwargs):
            self.set_cursor("watch")
            try:
                return func(self, *args, **kwargs)
            finally:
                self.set_cursor()
        return wrap

    @run_long_operation
    def cb_run_all(self) -> bool:
        logging.info('run all started...')
        for action in (
                        self.cb_trace_spectrum, 
                        self.cb_extract_spectrum, 
                        self.cb_calibrate_spectrum,
                        self.cb_apply_response,
                        self.cb_smooth_spectrum): 
             if action() is not True:
                 logging.error('run all aborted')
                 return False
             
        return True
    

    @run_long_operation
    def cb_trace_spectrum(self) -> bool:
        trace = self.spectrum_processor.trace(self.image_processor.img_stacked)
        if trace is None:
            return False
        self.spectrum_panel.draw_trace_on_image(
            self.image_panel.axe, trace.trace
        )
        return True

    @run_long_operation
    def cb_extract_spectrum(self) -> bool:
        image = self.image_processor.img_stacked
        spectrum = self.spectrum_processor.extract(image)
        if spectrum is None:
            return False

        trace = self.spectrum_processor.science_trace
        if trace is None:
            logging.error("spectrum trace is unavailable after extraction")
            return False
        trace_y_size = self.conf.get_int('processing', 'trace_y_size') or 15
        sky_sub = self.conf.get_bool('processing', 'sky_substract')
        if sky_sub in (None, True):
            sky_offset = self.conf.get_int('processing', 'sky_y_offset')
            sky_size = self.conf.get_int('processing', 'sky_y_size')
        else:
            sky_offset = sky_size = None

        self.spectrum_panel.draw_extraction_zones(
            img_axe=self.image_panel.axe,
            spectral_axis=spectrum.spectral_axis,
            trace=trace.trace,
            trace_y_size=trace_y_size,
            sky_y_offset=sky_offset,
            sky_y_size=sky_size,
        )
        self._show_spectrum(
            name=SpectrumProcessor.get_object_name(image.header),
            spectrum=spectrum,
            calibrated=False,
        )
        return True

    @run_long_operation
    def cb_calibrate_spectrum(self) -> bool:
        spectrum = self.spectrum_processor.calibrate(self.image_processor.img_stacked)
        if spectrum is None:
            return False
        self._show_spectrum(
            name=SpectrumProcessor.get_object_name(
                self.image_processor.img_stacked.header
            ),
            spectrum=spectrum,
            calibrated=True,
        )
        return True

    @run_long_operation
    def cb_apply_response(self) -> bool:
        spectrum = self.spectrum_processor.apply_response()
        if spectrum is None:
            return False
        self._show_spectrum(
            name=SpectrumProcessor.get_object_name(
                self.image_processor.img_stacked.header
            ),
            spectrum=spectrum,
            calibrated=True,
        )
        return True

    @run_long_operation
    def cb_smooth_spectrum(self) -> bool:
        spectrum = self.spectrum_processor.smooth()
        if spectrum is None:
            return False
        self._show_spectrum(
            name=SpectrumProcessor.get_object_name(
                self.image_processor.img_stacked.header
            ),
            spectrum=spectrum,
            calibrated=True,
        )
        return True

    @run_long_operation
    def cb_slitless(self) -> bool:
        # restore stacked image
        self.image_panel.clear()
        self.image_panel.show_image(
            self.image_processor.img_stacked.data, show_colorbar=True
        )
        self._refresh_cuts()

        # process sa100 images
        self.image_panel.toggle_sa100()
        return True

    @run_long_operation
    def cb_open_files(self) -> bool:
        """
        Prompt for input files and dispatch them to the processing components.

        Returns:
            bool: True when at least one file was processed
        """
        paths, _selected_filter = QFileDialog.getOpenFileNames(
            self,
            'Select image(s) or spectrum(s)',
            self.conf.get_conf_directory(),
            'FITS and DAT files (*.fit *.fts *.fits *.dat)',
        )

        if not paths:
            return False

        self.conf.set_conf_directory(os_utils.get_path_directory(path=paths[0]))
        self._start_file_monitor()

        spectrum_paths, image_paths = classify_files(paths)

        for spc_path in spectrum_paths:
            spectrum = self.spectrum_processor.load_spectrum(spc_path)
            if spectrum is not None:
                self._show_spectrum(
                    name=Path(spc_path).stem,
                    spectrum=spectrum,
                    calibrated=True,
                    clear=False,
                )

        if image_paths:
            self.cb_clear_image()
            self._load_images(image_paths)

        return True

    def _on_sa100_points(
        self, x0: float, y0: float, x1: float, y1: float, x2: float, y2: float
    ) -> None:
        result = sa100.rectify_sa100(self.image_processor.img_stacked, 
                                     x0, y0, x1, y1, y2)
        if result is None:
            return
        (
            rectified,
            spectrum_strip,
            sky_strip,
            x0_rotated,
            y0_rotated,
            half_height,
            sky_offset,
            sky_half_height,
        ) = result

        self.image_panel.clear()
        self.image_panel.show_image(rectified.data, show_colorbar=True)
        self._refresh_cuts()
        self.image_panel.draw_sa100_zones(
            x0=x0_rotated,
            y0=y0_rotated,
            x1=rectified.data.shape[1],
            y1=y0_rotated,
            half_height=half_height,
            sky_offset=sky_offset,
            sky_half_height=sky_half_height,
        )
        self._extract_sa100(spectrum_strip, sky_strip, x0_rotated)

    def _extract_sa100(self, spectrum_strip, sky_strip, x0_rotated: float) -> None:
        """
        Extract and display the spectrum after the user selects SA100 reference points.
        """
        self.set_cursor("watch")
        try:
            spectrum = self.spectrum_processor.sa100_extract(
                spectrum_strip, sky_strip, x0_rotated
            )
            if spectrum is not None:
                self._show_spectrum(
                    name=SpectrumProcessor.get_object_name(spectrum_strip.header),
                    spectrum=spectrum,
                    calibrated=True,
                )
        finally:
            self.set_cursor()

    def cb_clear_image(self) -> None:
        self.image_processor.reset()
        self.image_panel.clear()
        self.image_panel.show_image(
            self.image_processor.img_stacked.data, show_colorbar=False
        )

    def _load_images(self, paths: list[str]) -> bool:
        image = self.image_processor.load(paths)
        if image is None:
            return False
        self.image_panel.show_image(image.data, show_colorbar=True)
        self._refresh_cuts()
        return True

    def cb_clear_spectrum(self) -> None:
        self.spectrum_processor.reset()
        self.spectrum_panel.clear()

    def _show_spectrum(
        self, name: str, spectrum, calibrated: bool, clear: bool = True
    ) -> None:
        if clear:
            self.spectrum_panel.clear()
        self.spectrum_panel.show_spectrum(
            name=name, spectrum=spectrum, calibrated=calibrated
        )

    def cb_show_lines(self) -> None:
        self.conf.reload_if_changed()
        if self.spectrum_processor.science_spectrum is None:
            logging.info("please calibrate first")
            return
        if not self.spectrum_processor.is_calibrated():
            logging.warning("spectrum needs to be calibrated first")
            return
        lines = [
            (float(wavelength) * 10, label)
            for wavelength, label in self.conf.config.items('lines')
        ]
        self.spectrum_panel.toggle_lines(lines)

    def cb_colorize(self) -> None:
        if self.spectrum_processor.science_spectrum is None:
            logging.info("please calibrate first")
            return
        if not self.spectrum_processor.is_calibrated():
            logging.warning("spectrum needs to be calibrated first")
            return

        spectrum = self.spectrum_processor.science_spectrum
        self.spectrum_panel.toggle_colorize(
            wavelengths=spectrum.wavelength.value,
            flux=spectrum.flux.value,
            min_wl=self.conf.get_int('display', 'colorize_min_wl') or 3800,
            max_wl=self.conf.get_int('display', 'colorize_max_wl') or 7500,
            bin_size=self.conf.get_int('display', 'colorize_bin_size') or 20,
        )

    def _refresh_cuts(self) -> None:
        low, high, slider_min, slider_max = self.image_processor.compute_cuts()
        self.image_panel.set_slider_range(slider_min, slider_max, low, high)

    def _on_cut_changed(
        self, low_cut: float | None, high_cut: float | None
    ) -> None:
        self.image_panel.apply_cut(low_cut, high_cut)

    def _start_file_monitor(self) -> None:
        config_path = Path(self.conf.config_path).resolve()
        self.file_monitor.start(str(config_path.parent), str(config_path))

    def _on_config_file_changed(self) -> None:
        self._config_reload_timer.start()

    def _reload_config(self) -> None:
        if self.conf.reload():
            previous_theme = self._theme
            self._apply_theme()
            if self._theme != previous_theme and hasattr(self, 'spectrum_panel'):
                self.spectrum_panel.clear()

    def closeEvent(self, event) -> None:
        self.file_monitor.stop()
        event.accept()

    def _on_new_fits_file(self, filepath: str) -> None:
        """
        Callback fired by the watchdog thread when a new FITS file appears.
        Queues processing on the Qt main thread via a Qt signal.

        Args:
            filepath (str): absolute path of the new FITS file
        """
        self._process_new_fits(filepath)

    def _process_new_fits(self, filepath: str) -> None:
        """
        Load and process a newly detected FITS file (runs on the Qt main thread).

        Args:
            filepath (str): absolute path of the new FITS file
        """
        if self.conf.get_bool('processing', 'auto_process') not in (None, True):
            return
        logging.info(f"new FIT file detected: {filepath}")
        self.set_cursor("watch")
        try:
            self.cb_clear_image()
            logging.info(f"loading {filepath}...")
            if not self._load_images([filepath]):
                return
            logging.info(f"processing {filepath}...")
            self.cb_run_all()
        finally:
            self.set_cursor()
