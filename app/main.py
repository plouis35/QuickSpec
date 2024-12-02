import logging
import time
from pathlib import Path

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
        self.title(f"{app_name} - {app_version}")
        self.app_name = app_name
        self.app_version = app_version
        self.conf: Config = Config()
        LogHandler().initialize()
        self.create_panels()
        self.create_buttons()

        # create a timer to capture new files created
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
        bt_load = ttk.Button(self.bt_frame, text="Load", command=self.open_files)
        bt_load.pack(side=tk.LEFT, padx=5, pady=0)

        bt_run = ttk.Button(self.bt_frame, text="Run all", command=self.run_all)
        bt_run.pack(side=tk.LEFT, padx=5, pady=0)

        bt_reduce = ttk.Button(master=self.bt_frame, text="Reduce", command=self.reduce_images)
        bt_reduce.pack(side=tk.LEFT, padx=5, pady=0)

        bt_extract = ttk.Button(master=self.bt_frame, text="Extract", command=self.extract_spectrum)
        bt_extract.pack(side=tk.LEFT, padx=5, pady=0)

        bt_calibrate = ttk.Button(master=self.bt_frame, text="Calibrate", command=self.calibrate_spectrum)
        bt_calibrate.pack(side=tk.LEFT, padx=5, pady=0)

        bt_lines = ttk.Button(master=self.bt_frame, text="Show lines", command=self._spectrum.show_lines)
        bt_lines.pack(side=tk.LEFT, padx=5, pady=0)

        #bt_clear = ttk.Button(master=frame, text="Clear", command=self._spectrum.do_clear)
        #bt_clear.pack(side=tk.LEFT, padx=5, pady=0)

    # local callbacks for buttons
    def run_all(self) -> None:
        for action in ( self.reduce_images, 
                        self.extract_spectrum, 
                        self.calibrate_spectrum): action()
    
    def extract_spectrum(self) -> None:
        self._spectrum.clear_spectrum()
        self._spectrum.do_extract(self._image.img_stacked.data)

    def reduce_images(self) -> None:
        #self._image.clear_image()
        self._image.reduce_images()

    def calibrate_spectrum(self) -> None:
        self._spectrum.clear_spectrum()
        self._spectrum.do_calibrate()

    def open_files(self) -> None:
        # get list of files to open
        path = askopenfilenames(title='Select image(s) or spectra',
                            filetypes=[("fits files", '*.fit'), 
                                       ("fits files", "*.fts"), 
                                       ("fits files", "*.fits")
                                       ],
                            )
        if path == '': return 

        # create new config file if not existing
        self.conf.set_conf_directory(OSUtils.get_path_directory(path=path[0]))

        # check file type to call proper load routine 
        img_names: list[str] = []
        for img_name in path:
            fit_data: CCDData = CCDData.read(img_name, unit=u.dimensionless_unscaled)
            if fit_data.ndim == 1:
                # this is a spectrum fit file: display directly
                logging.info(f"{img_name} is a spectrum (naxis = 1)")
                self._spectrum.open_spectrum(img_name)

            elif fit_data.ndim == 2:
                # this is a fit image - keep it to load it later
                logging.info(f"{img_name} is a fit image (naxis = 2)")
                img_names.append(img_name)
            else:
                # not supported fit format
                logging.error(f"{img_name} is not a supported fit format (naxis > 2)")

        if len(img_names) > 0:
            self._image.load_images(img_names)
        
        self.title(f"{self.app_name} - {self.app_version} [{Path(path[0]).stem}...]")

        return 

    def watch_files(self):
        #logging.info(f"watcher time is : {time.strftime("%H:%M:%S", time.localtime())}")
        self.after(1000, self.watch_files) 
