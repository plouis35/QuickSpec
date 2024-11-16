import logging
import threading
import time

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
from app.config import Config
from app.os_utils import OSUtils

from img.image import Image
from spc.spectrum import Spectrum

class Application(tk.Tk):

    def __init__(self, app_name: str, app_version: str) -> None:
        """_summary_

        Args:
            app_name (str): _description_
            app_version (str): _description_
        """
        super().__init__()
        self.title(app_name)
        self.app_name = app_name

        self.conf = Config()
        LogHandler().set()
        OSUtils.show_versions()
        #self.create_watcher()
        self.create_gui()


    def create_gui(self) -> None:
        plt.style.use('dark_background')        
        plt.rcParams['figure.constrained_layout.use'] = True
        self.geometry(self.conf.get_str('display', 'geometry'))
        self.wm_state(self.conf.get_str('display', 'state'))

        # create a unique figure with 2 axes to hold image & spectrum horizontaly stacked
        self.figure = Figure(figsize=(5, 4)) #, dpi=100)
        self.axe_img = self.figure.add_subplot(211)
        self.axe_spc = self.figure.add_subplot(212)
        logging.debug('figure created')

        # create top frame to hold buttons and sliders
        frame = ttk.Frame(self)
        frame.pack(side=tk.TOP, fill=tk.X)

        # intialize img axe
        self._image = Image(self.axe_img, self.axe_spc, frame)

        # initialize spec axe
        self._spectrum = Spectrum(self.axe_img, self.axe_spc)

        # create buttons
        bt_load = ttk.Button(frame, text="Load", command=self.open_files)
        bt_load.pack(side=tk.LEFT, padx=5, pady=0)

        bt_run = ttk.Button(frame, text="Run all", command=self.run_all)
        bt_run.pack(side=tk.LEFT, padx=5, pady=0)

        bt_reduce = ttk.Button(master=frame, text="Reduce", command=self._image.reduce_images)
        bt_reduce.pack(side=tk.LEFT, padx=5, pady=0)

        bt_extract = ttk.Button(master=frame, text="Extract", command=self._spectrum.do_extract)
        bt_extract.pack(side=tk.LEFT, padx=5, pady=0)

        bt_calibrate = ttk.Button(master=frame, text="Calibrate", command=self._spectrum.do_calibrate)
        bt_calibrate.pack(side=tk.LEFT, padx=5, pady=0)

        bt_response = ttk.Button(master=frame, text="Response", command=self._spectrum.do_response)
        bt_response.pack(side=tk.LEFT, padx=5, pady=0)

        # create canvas
        canvas = FigureCanvasTkAgg(self.figure, self)
        canvas.draw()

        # create toolbar
        toolbar = NavigationToolbar2Tk(canvas, self, pack_toolbar=False)
        toolbar.children['!button4'].pack_forget()      # ugly... should use another method to remove the conf button.
        #toolbar.config(background='grey')        
        toolbar.config(height=100)
        toolbar.update()

        # pack all widgets
        toolbar.pack(side=tk.BOTTOM, fill=tk.X)
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    # define callbacks for buttons
    def run_all(self) -> None:
        for action in ( self._image.reduce_images, 
                        self._spectrum.do_extract, 
                        self._spectrum.do_calibrate, 
                        self._spectrum.do_response ):
            if action() is False:
                logging.warning(f"runall aborted")
                return
    

    def open_files(self) -> None:
        # create openfile dialog
        path = askopenfilenames(title='Select image(s) or a directory for watch mode',
                            initialdir = self.conf.get_str('files', 'initial_directory'),
                            filetypes=[("fits files", '*.fit'), 
                                       ("fits files", "*.fts"), 
                                       ("fits files", "*.fits")
                                       ],
                            )
        #path = (['albireo-10.fit'])
        if path == '': 
            return 

        # check file type to call proper load routine 
        img_names: list[str] = []
        for img_name in path:
            fit_data: CCDData = CCDData.read(img_name, unit=u.dimensionless_unscaled)
            if fit_data.ndim == 1:
                # this is a spectrum fit file: display directly
                logging.info(f"{img_name} is a spectrum (naxis = 1)")
                self._spectrum.open_spectrum(img_name)

            elif fit_data.ndim == 2:
                # this is a fit image - keep it to stack it later
                logging.info(f"{img_name} is a fit image (naxis = 2)")
                img_names.append(img_name)
            else:
                # not supported fit format
                logging.error(f"{img_name} is not a supported fit format (naxis > 2)")

        self._image.load_image(img_names)
        return 
    
    def create_watcher(self) -> None:
        def watch() -> None:
            while True:
                try:
                    self.title(f"{self.app_name} [{OSUtils.get_memory_used()}MB]")
                    time.sleep(1)
                except Exception as e:
                    #logging.warning(f"Watch thread : {e}")
                    logging.info("Watch thread terminated")
                    # exit thread whatever exception is (during destroy window)
                    break

        self.watch_thread = threading.Thread(target = watch, name='quickspec_watch_thread')
        self.watch_thread.start()
        logging.info('watch thread started - tid = {}'.format(self.watch_thread))          

