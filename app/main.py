"""
main application class:
- initialize logging and configuration
- create tkinter main GUI
- instantiate 2D image and 1D spectrum classes
- start a watchdog thread (via the 'watchdog' library) to monitor new FITS files
"""
import logging
import os
from pathlib import Path
import functools

import tkinter as tk
from tkinter import ttk
from tkinter.filedialog import askopenfilenames

import matplotlib.pyplot as plt

from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler, FileCreatedEvent

from app.logger import LogHandler
import app.os_utils as os_utils
from app.file_utils import classify_files
from app.config import Config
from img.image import Image
from spc.spectrum import Spectrum

class Application(tk.Tk):

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

        # set window style from Microsoft azure
        self.tk.call("source", "azure.tcl")

        # set theme and colors
        if (_theme := self.conf.get_str('display', 'theme')) in (None, 'dark'):
            _tk_theme = 'dark'
            _mpl_theme = 'dark_background'

        elif _theme == 'light':
            _tk_theme = 'light'
            _mpl_theme = 'grayscale'
        else:
            logging.error(f"unsupported color theme: {_theme=}")
            _tk_theme = 'dark'
            _mpl_theme = 'dark_background'
            
        plt.style.use(_mpl_theme)        
        self.tk.call("set_theme", _tk_theme)
        self.tk.call('set', 'tk_strictMotif', '1')

        plt.rcParams['figure.constrained_layout.use'] = True
        
        self.title(f"{app_name} v{app_version}")
        self.app_name = app_name
        self.app_version = app_version

        # create image and spectrum panels
        self.create_panels()
        self.create_buttons()

        # start watchdog for automatic FITS file processing
        self._start_watchdog()

        # and display major packages versions installed
        logging.info(f"{app_name} v{app_version} started")
        os_utils.show_versions()

    def create_panels(self) -> None:
        """
        create tkinter panels for 2D and 1D spectrum display
        instantiate image and spectrum classes to manage them
        """        
        # create top frame to hold buttons and sliders
        self.bt_frame = ttk.Frame(self)
        self.bt_frame.pack(side=tk.TOP, fill=tk.X, pady=5)

        # create two frames to hold image and spectrum
        paned_window = ttk.PanedWindow(self, orient=tk.VERTICAL)
        paned_window.pack(fill=tk.BOTH, expand=True)
        img_frame = ttk.Frame(paned_window)
        spc_frame = ttk.Frame(paned_window)
        paned_window.add(child=img_frame)
        paned_window.add(child=spc_frame)

        # intialize image axe
        self._image = Image(img_frame, self.bt_frame)

        # initialize spectrum axe
        self._spectrum = Spectrum(spc_frame, self._image.img_axe)

    def create_buttons(self) -> None:
        """
        creates buttons to load, process_all, process_step_by_step and SA100 mode
        """        
        bt_load = ttk.Button(self.bt_frame, text="Load & reduce", command=self.cb_open_files) 
        bt_load.pack(side=tk.LEFT, padx=5, pady=0)

        bt_run = ttk.Button(self.bt_frame, text="Run all", command=self.cb_run_all)
        bt_run.pack(side=tk.LEFT, padx=5, pady=0)

        _step_options = [
                                    "Find spectrum",
                                    "Extract spectrum", 
                                    "Calibrate spectrum", 
                                    "Apply response",
                                    "Smooth, crop & normalize"
                                    ]
        bt_step_default = "Run step"
        _var = tk.StringVar(value=bt_step_default)

        _step_map = {
            _step_options[0]: self.cb_trace_spectrum,
            _step_options[1]: self.cb_extract_spectrum,
            _step_options[2]: self.cb_calibrate_spectrum,
            _step_options[3]: self.cb_apply_response,
            _step_options[4]: self.cb_smooth_spectrum,
        }

        def cb_run_step(selected_step: tk.StringVar) -> None:
            logging.info(f"step {selected_step} started...")
            if (action := _step_map.get(selected_step)) is not None:
                action()
            else:
                logging.warning(f"unknown step: {selected_step}")
            _var.set(bt_step_default)

        bt_steps = ttk.OptionMenu(self.bt_frame, _var, bt_step_default, *(_step_options), command = cb_run_step)
        bt_steps.pack(side=tk.LEFT, padx=5, pady=0)

        bt_sa100 = ttk.Button(self.bt_frame, text="Slitless", command=self.cb_slitless) 
        bt_sa100.pack(side=tk.LEFT, padx=5, pady=0)


    def set_cursor(self, cursor: str = '') -> None:
        """
        set cursor icon to 'hourglass' mode
        NOTE: works well on Linux and macOS - not on Windows...

        Args:
            mode (str) : either 'watch' (hourglass) or '' (back to default)
        """        
        self.config(cursor=cursor)
        self.update()

    def set_title(self, title: str = '') -> None:
        """
        set window title

        Args:
            title (str, optional): new title. Defaults to ''.
        """        
        self.title(f"{self.app_name} v{self.app_version} - {title}")
        
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
            retcode = func(self, *args, **kwargs)
            self.set_cursor()    
            return retcode
        return wrap

    @run_long_operation
    def cb_run_all(self) -> bool:
        logging.info('run all started...')
        #for action in ( self.cb_reduce_images, 
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
        return self._spectrum.do_trace(self._image.img_stacked)

    @run_long_operation
    def cb_extract_spectrum(self) -> bool:
        return self._spectrum.do_extract(self._image.img_stacked)

    @run_long_operation
    def cb_calibrate_spectrum(self) -> bool:
        return self._spectrum.do_calibrate(self._image.img_stacked)

    @run_long_operation
    def cb_apply_response(self) -> bool:
        return self._spectrum.do_response(self._image.img_stacked)

    @run_long_operation
    def cb_smooth_spectrum(self) -> bool:
        return self._spectrum.do_smooth(self._image.img_stacked)

    @run_long_operation
    def cb_slitless(self) -> bool:
        return True

    @run_long_operation
    def cb_open_files(self) -> bool:
        """
        Prompt user to select files, classify them, and dispatch to
        Image or Spectrum controllers.

        Returns:
            bool: True when at least one file was processed
        """
        #self.tk.call('set', 'tk_strictMotif', '1')

        paths = askopenfilenames(
            title='Select image(s) or spectrum(s)',
            filetypes=[
                ("fits files", '*.fit'),
                ("fits files", "*.fts"),
                ("fits files", "*.fits"),
                ("dat files", "*.dat"),
            ],
        )
        #self.tk.call('set', 'tk_strictMotif', '0')

        if not paths:
            return False

        self.conf.set_conf_directory(os_utils.get_path_directory(path=paths[0]))

        spectrum_paths, image_paths = classify_files(paths)

        for spc_path in spectrum_paths:
            self._spectrum.open_spectrum(spc_path)

        if image_paths:
            self._image.clear_image()
            self._image.load_images(image_paths)

        return True

    def _start_watchdog(self) -> None:
        """
        Start a watchdog Observer thread that reacts to new FITS file creation
        in the currently monitored directory.
        """
        self._watchdog_observer: Observer | None = None

        auto_process = self.conf.get_bool('processing', 'auto_process')
        if auto_process not in (None, True):
            logging.debug("watchdog disabled by configuration")
            return

        watch_path = os_utils.get_current_path()
        if watch_path == '.':
            logging.debug("watchdog: no directory selected yet, skipping start")
            return

        handler = _FitsEventHandler(on_new_file=self._on_new_fits_file)
        self._watchdog_observer = Observer()
        self._watchdog_observer.schedule(handler, path=watch_path, recursive=False)
        self._watchdog_observer.start()
        logging.info(f"watchdog started on: {watch_path}")

    def _stop_watchdog(self) -> None:
        """Stop the watchdog Observer thread if running."""
        if getattr(self, '_watchdog_observer', None) is not None:
            self._watchdog_observer.stop()
            self._watchdog_observer.join()
            self._watchdog_observer = None
            logging.info("watchdog stopped")

    def _on_new_fits_file(self, filepath: str) -> None:
        """
        Callback fired by the watchdog thread when a new FITS file appears.
        Schedules the actual processing on the Tk main thread via after().

        Args:
            filepath (str): absolute path of the new FITS file
        """
        logging.info(f"new FIT file detected: {filepath}")
        # Marshal back to Tk thread — never touch GUI from a background thread
        self.after(0, lambda: self._process_new_fits(filepath))

    def _process_new_fits(self, filepath: str) -> None:
        """
        Load and process a newly detected FITS file (runs on Tk main thread).

        Args:
            filepath (str): absolute path of the new FITS file
        """
        self.set_cursor("watch")
        self._image.clear_image()
        logging.info(f"loading {filepath}...")
        self._image.load_images([filepath])
        logging.info(f"processing {filepath}...")
        self.cb_run_all()
        self.set_cursor()



class _FitsEventHandler(FileSystemEventHandler):
    """
    Watchdog event handler that fires a callback when a new FITS file appears.
    Runs in a background thread — must NOT touch the Tk GUI directly.
    """
    FITS_SUFFIXES = {'.fit', '.fts', '.fits'}

    def __init__(self, on_new_file) -> None:
        """
        Args:
            on_new_file (callable): called with the new file path (str)
        """
        super().__init__()
        self._on_new_file = on_new_file

    def on_created(self, event: FileCreatedEvent) -> None:
        if event.is_directory:
            return
        if Path(event.src_path).suffix.lower() in self.FITS_SUFFIXES:
            self._on_new_file(event.src_path)
