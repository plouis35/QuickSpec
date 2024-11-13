import logging
import threading
import time

import tkinter as tk
from tkinter import ttk

import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends._backend_tk import NavigationToolbar2Tk

from app.logger import LogHandler
from app.config import Config
from app.os_utils import OSUtils

from img.image import Image
from spc.spectrum import Spectrum

class Application(tk.Tk):

    def __init__(self, app_name: str, app_version: str) -> None:
        super().__init__()

        self.conf = Config()
        LogHandler().set()
        OSUtils.log_versions()
        #self.tk.call('tk', 'scaling', '-displayof', '.', 2)
        self.geometry(self.conf.get_str('window', 'geometry'))
        self.title(app_name)
        self.app_name = app_name
        self.wm_state(self.conf.get_str('window', 'state'))
        #self.create_watcher()
        self.create_gui()

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


    def create_gui(self) -> None:
        plt.rcParams['figure.constrained_layout.use'] = True
        plt.style.use(self.conf.get_str('window', 'theme'))        

        # create a single figure for both image & spectrum horizontaly packed
        self.figure = Figure(figsize=(5, 4)) #, dpi=100)

        self.axe_img = self.figure.add_subplot(211)
        self.axe_spc = self.figure.add_subplot(212)
        logging.debug('figure created')

        # create top frame to hold buttons and sliders
        frame = ttk.Frame(self)
        frame.pack(side=tk.TOP, fill=tk.X)

        # intialize img tools
        self._image = Image(self.axe_img, self.axe_spc, frame)

        # initialize spec tools
        self._spectrum = Spectrum(self.axe_img, self.axe_spc)

        # create buttons
        bt_load = ttk.Button(frame, text="Load", command=self._image.open_image)
        bt_load.pack(side=tk.LEFT, padx=5, pady=0)
        bt_run = ttk.Button(frame, text="Run all", command=self._spectrum.do_run_all)
        bt_run.pack(side=tk.LEFT, padx=5, pady=0)
        bt_reduce = ttk.Button(master=frame, text="Reduce", command=self._spectrum.do_reduce)
        bt_reduce.pack(side=tk.LEFT, padx=5, pady=0)
        bt_extract = ttk.Button(master=frame, text="Extract", command=self._spectrum.do_extract)
        bt_extract.pack(side=tk.LEFT, padx=5, pady=0)
        bt_calibrate = ttk.Button(master=frame, text="Calibrate", command=self._spectrum.do_calibrate)
        bt_calibrate.pack(side=tk.LEFT, padx=5, pady=0)
        bt_response = ttk.Button(master=frame, text="Response", command=self._spectrum.do_response)
        bt_response.pack(side=tk.LEFT, padx=5, pady=0)

        #bt_clear = ttk.Button(master=frame, text="Clear", command=self._image.clear_image)
        #bt_clear.pack(side=tk.LEFT, padx=5, pady=0)

        # create canvas 
        canvas = FigureCanvasTkAgg(self.figure, self)
        canvas.draw()

        # create toolbar
        toolbar = NavigationToolbar2Tk(canvas, self, pack_toolbar=False)
        toolbar.children['!button4'].pack_forget()
        #toolbar.config(background='grey')        
        toolbar.update()

        # pack all
        toolbar.pack(side=tk.BOTTOM, fill=tk.X)
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)


