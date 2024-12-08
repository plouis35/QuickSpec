import logging
import time
from pathlib import Path
import functools
import os

import tkinter as tk
from tkinter import ttk
from tkinter.filedialog import askopenfilenames

import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends._backend_tk import NavigationToolbar2Tk

from astropy import units as u
from astropy.nddata import CCDData

from app.logger import LogHandler
from app.os_utils import OSUtils

from app.config import Config
from img.image import Image
from app.os_utils import OSUtils
from spc.spectrum import Spectrum

class Application(tk.Tk):

    def __init__(self, app_name: str, app_version: str) -> None:
        super().__init__()

        self.iconbitmap(r'./quickspec.ico')
        self.title(f"{app_name} - {app_version}")
        self.app_name = app_name
        self.app_version = app_version
        self.conf: Config = Config()
        LogHandler().initialize()
        self.create_panels()
        self.create_buttons()

        # create a timer to capture new files created
        self.last_timer: float = time.time()
        self.after_idle(self.watch_files)

    def create_panels(self) -> None:
        plt.style.use('dark_background')        
        plt.rcParams['figure.constrained_layout.use'] = True
        
        # create top frame to hold buttons and sliders
        self.bt_frame = ttk.Frame(self)
        self.bt_frame.pack(side=tk.TOP, fill=tk.X)

        # create two frames to hold image and spectrum
        paned_window = ttk.PanedWindow(self, orient=tk.VERTICAL)
        paned_window.pack(fill=tk.BOTH, expand=True)
        img_frame = ttk.Frame(paned_window)
        spc_frame = ttk.Frame(paned_window)
        paned_window.add(child=img_frame) #, weight=1)
        paned_window.add(child=spc_frame) #, weight=1)

        # intialize image axe
        self._image = Image(img_frame, self.bt_frame)

        # initialize spectrum axe
        self._spectrum = Spectrum(spc_frame, self._image.img_axe)

    def create_buttons(self) -> None:
        bt_load = ttk.Button(self.bt_frame, text="Load", command=self.cb_open_files) 
        bt_load.pack(side=tk.LEFT, padx=5, pady=0)

        bt_run = ttk.Button(self.bt_frame, text="Run all", command=self.cb_run_all)
        bt_run.pack(side=tk.LEFT, padx=5, pady=0)

        _step_options = ["Reduce image(s)", "Find spectrum", "Extract spectrum", "Calibrate spectrum"]
        bt_step_default = "Run step"
        _var = tk.StringVar(value=bt_step_default)

        # TODO : use callbacks instead
        def cb_run_step(selected_step: tk.StringVar) -> None:
            logging.info(f"step {selected_step} started...")
            if selected_step == _step_options[0]: self.cb_reduce_images()
            elif selected_step == _step_options[1]: self.cb_trace_spectrum()
            elif selected_step == _step_options[2]: self.cb_extract_spectrum()
            elif selected_step == _step_options[3]: self.cb_calibrate_spectrum()
            else: pass 
            _var.set(bt_step_default)

        bt_steps = ttk.OptionMenu(self.bt_frame, _var, bt_step_default, *(_step_options), command = cb_run_step)
        bt_steps.pack(side=tk.LEFT, padx=5, pady=0)

    # local callbacks for buttons
    @staticmethod
    def run_long_operation(func):
        @functools.wraps(func)
        def wrap(self, *args, **kwargs):
            self.config(cursor="watch")
            self.update()
            retcode = func(self, *args, **kwargs)
            self.config(cursor="")    
            self.update()
            return retcode
        return wrap

    @run_long_operation
    def cb_run_all(self) -> bool:
        for action in ( self.cb_reduce_images, 
                        self.cb_trace_spectrum, 
                        self.cb_extract_spectrum, 
                        self.cb_calibrate_spectrum): 
             if action() is not True:
                 logging.error('run all aborted')
                 return False
        return True
    
    @run_long_operation
    def cb_reduce_images(self) -> bool:
        #self._image.clear_image()
        return self._image.reduce_images()

    @run_long_operation
    def cb_trace_spectrum(self) -> bool:
        self._spectrum.clear_spectra()
        return self._spectrum.do_trace(self._image.img_stacked.data)

    @run_long_operation
    def cb_extract_spectrum(self) -> bool:
        self._spectrum.clear_spectra()
        return self._spectrum.do_extract(self._image.img_stacked.data)

    @run_long_operation
    def cb_calibrate_spectrum(self) -> bool:
        self._spectrum.clear_spectra()
        return self._spectrum.do_calibrate()

    def cb_open_files(self) -> bool:
        # get list of files to open
        path = askopenfilenames(title='Select image(s) or spectrum(s)',
                            filetypes=[("fits files", '*.fit'), 
                                       ("fits files", "*.fts"), 
                                       ("fits files", "*.fits")
                                       ],
                            )
        if path == '': return False

        # create new config file if not existing
        self.conf.set_conf_directory(OSUtils.get_path_directory(path=path[0]))

        # check header to either load an image 2D or a spectrum 1D
        img_names: list[str] = []
        for img_name in path:
            fit_data: CCDData = CCDData.read(img_name, unit=u.dimensionless_unscaled)
            if fit_data.ndim == 1:
                # this is a spectrum fit file: display directly
                logging.info(f"{img_name} is a spectrum (naxis = 1)")
                self._spectrum.open_spectrum(img_name)

            elif fit_data.ndim == 2:
                # this is a fit image - keep it to load them all together
                logging.info(f"{img_name} is a fit image (naxis = 2)")
                img_names.append(img_name)
            else:
                # not supported fit format
                logging.error(f"{img_name} is not a supported fit format (naxis > 2)")

        # clear existing images
        self._image.clear_image()

        # show wait cursor while loading images
        if len(img_names) > 0:
            self.config(cursor="watch")
            self.update()
            self._image.load_images(img_names)
            self.config(cursor="")    
        
        # update names in title
        if len(img_names) > 0:
            self.title(f"{self.app_name} - {self.app_version} [{Path(path[0]).stem}..{Path(path[-1]).stem}]")
        else:
            self.title(f"{self.app_name} - {self.app_version} [{Path(path[0]).stem}]")

        return True

    def watch_files(self):
        #logging.debug(f"watcher time is : {time.strftime("%H:%M:%S", time.localtime())}")

        if ((path := OSUtils.get_current_path()) != '.'):
            new_file = OSUtils.list_files(path, '*.fit*')[0]
            if os.path.getmtime(f"{path}/{new_file}") >= self.last_timer:
                logging.info(f"new FIT detected: {new_file}")
                self.last_timer = time.time()

        self.after(1000, self.watch_files) 
