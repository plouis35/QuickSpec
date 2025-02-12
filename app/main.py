"""
main application class:
- initialize logging and configuration
- creates tkinter main GUI
- instantiates 2D image and 1D spectrum classes
- creates a watchdog routine to monitor new file creation for automatic spectrum processing
"""
import logging
import time
from pathlib import Path
import functools
import os

import tkinter as tk
from tkinter import ttk
from tkinter.filedialog import askopenfilenames

import matplotlib.pyplot as plt

from astropy import units as u
from astropy.nddata import CCDData

from app.logger import LogHandler
from app.os_utils import os_utils

from app.config import Config
from img.image import Image
from app.os_utils import os_utils
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
        plt.rcParams['figure.constrained_layout.use'] = True
        
        self.title(f"{app_name} v{app_version}")
        self.app_name = app_name
        self.app_version = app_version

        # create image and spectrum panels
        self.create_panels()
        self.create_buttons()

        # create a watchdog timer to capture new files creation
        self._last_timer: float = time.time()
        self.after_idle(self.watch_files)

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
        creates global buttons to load, process_all and process_step_by_step spectra
        """        
        bt_load = ttk.Button(self.bt_frame, text="Load", command=self.cb_open_files) 
        bt_load.pack(side=tk.LEFT, padx=5, pady=0)

        bt_run = ttk.Button(self.bt_frame, text="Run all", command=self.cb_run_all)
        bt_run.pack(side=tk.LEFT, padx=5, pady=0)

        _step_options = ["Reduce image(s)", 
                                    "Trace spectrum",
                                    "Extract spectrum", 
                                    "Calibrate spectrum", 
                                    "Apply response",
                                    "Smooth spectrum"]
        bt_step_default = "Run step"
        _var = tk.StringVar(value=bt_step_default)

        def cb_run_step(selected_step: tk.StringVar) -> None:
            logging.info(f"step {selected_step} started...")
            if selected_step == _step_options[0]: self.cb_reduce_images()
            elif selected_step == _step_options[1]: self.cb_trace_spectrum()
            elif selected_step == _step_options[2]: self.cb_extract_spectrum()
            elif selected_step == _step_options[3]: self.cb_calibrate_spectrum()
            elif selected_step == _step_options[4]: self.cb_apply_response()
            elif selected_step == _step_options[5]: self.cb_smooth_spectrum()
            else: pass 
            _var.set(bt_step_default)

        bt_steps = ttk.OptionMenu(self.bt_frame, _var, bt_step_default, *(_step_options), command = cb_run_step)
        bt_steps.pack(side=tk.LEFT, padx=5, pady=0)

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
        for action in ( self.cb_reduce_images, 
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
    def cb_reduce_images(self) -> bool:
        return self._image.reduce_images()

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

    def cb_open_files(self) -> bool:
        """
        request user to select file(s) to open
        open and display selected files according to type (2D or 1D spectra)

        Returns:
            bool: True when successfull
        """        
        # get list of files to open
        path = askopenfilenames(title='Select image(s) or spectrum(s)',
                            filetypes=[("fits files", '*.fit'), 
                                       ("fits files", "*.fts"), 
                                       ("fits files", "*.fits")
                                       ],
                            )
        if path == '': return False

        # create new config file if not existing
        self.conf.set_conf_directory(os_utils.get_path_directory(path=path[0]))

        # check header to either load an 2D image or a 1D spectrum
        img_names: list[str] = []
        for img_name in path:
            fit_data: CCDData = CCDData.read(img_name, unit=u.dimensionless_unscaled)
            if fit_data.ndim == 1:
                # this is a spectrum fit file: display directly
                logging.info(f"{img_name} is a spectrum (naxis = 1, shape = {fit_data.shape})")
                self._spectrum.open_spectrum(img_name)

            elif fit_data.ndim == 2:
                # this is a fit image - keep it to load them all together
                logging.info(f"{img_name} is a fit image (naxis = 2, shape = {fit_data.shape})")
                img_names.append(img_name)
            else:
                # not supported fit format
                logging.error(f"{img_name} is not a supported fit format (naxis > 2)")

        # load and stack selected images 
        if len(img_names) > 0:
            self._image.clear_image()
            self.set_cursor("watch")
            self._image.load_images(img_names)
            self.set_cursor()    
        
        # set window title according to images names
        #if len(img_names) > 1:
         #   self.set_title(f"[{Path(path[0]).stem} .. {Path(path[-1]).stem}]")
        #elif len(img_names) == 1:
          #  self.set_title(f"[{Path(path[0]).stem}]")
        #else:
         #   pass

        return True

    def watch_files(self):
        """
        watchdog monitor routine
        used to check whether a new file is created under selected directory
        if so, process the new file using run_all callback
        TODO: raise errors if the new file is a 1D spectrum (ignored for now)
        """        

        logging.debug(f"watchdog current time is : {time.strftime('%H:%M:%S', time.localtime())}")
        auto_process: bool | None = self.conf.get_bool('processing', 'auto_process')
        if auto_process in (None, True):
            logging.debug(f"watchdog enabled")
            if ((path := os_utils.get_current_path()) != '.'):
                new_file = os_utils.list_files(path, '*.fit*')[0]
                if os.path.getmtime(f"{path}/{new_file}") >= self._last_timer:
                    logging.info(f"new FIT file detected: {new_file}")
                    self._last_timer = time.time()

                    # load and process the new file
                    self.set_cursor("watch")

                    self._image.clear_image()
                    logging.info(f"loading {new_file}...")
                    self._image.load_images([f"{path}/{new_file}"])
                    logging.info(f"processing {new_file}...")
                    self.cb_run_all()
                    
                    self.set_cursor()
        else:
            logging.debug(f"watchdog disabled")

        # restart watchdog every second
        self.after(1000, self.watch_files) 

