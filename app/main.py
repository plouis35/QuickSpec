import logging
import threading
import time
import os
import numpy as np

import tkinter as tk
from tkinter import ttk

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.image import AxesImage
from matplotlib.widgets import Button
from matplotlib.widgets import RangeSlider
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
                self.title(f"{self.app_name} [{OSUtils.get_memory_used()}MB]")
                time.sleep(1)

        self.watch_thread = threading.Thread(target = watch, name='EASYSPEC_watch_thread')
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

        # intialize img tools
        self._image = Image(self.axe_img, self.axe_spc)

        # initialize spec tools
        self._spectrum = Spectrum(self.axe_img, self.axe_spc)

        # callback for cuts sliders    
        def update_sliders(slider, label):
            label.config(text=f"[{slider.get():.0f}]")
            if self._image.image is not None:           
                if slider == slider_low:
                    self._image.update_image(slider.get(), None)
                elif slider == slider_high:
                    self._image.update_image(None, slider.get())
                else:
                    logging.error("internal : unknown slider event on repr({slider})")
                                
        # create button bar and sliders
        frame = ttk.Frame(self)
        frame.pack(side=tk.TOP, fill=tk.X)

        bt_load = ttk.Button(frame, text="Load", command=self._image.open_image)
        bt_load.pack(side=tk.LEFT, padx=5, pady=0)
        bt_run = ttk.Button(frame, text="Process", command=self._spectrum.do_calibration)
        bt_run.pack(side=tk.LEFT, padx=5, pady=0)

        slider_frame = ttk.Frame(frame)
        slider_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        # Slider high
        slider_high_frame = ttk.Frame(slider_frame)
        slider_high_frame.pack(side=tk.TOP, fill=tk.X, padx=50, pady=0)
        slider_high_value = ttk.Label(slider_high_frame, text="[65535]")
        slider_high_value.pack(side=tk.LEFT)
        #slider_high_label = ttk.Label(slider_high_frame, text="0")
        #slider_high_label.pack(side=tk.LEFT)
        slider_high = ttk.Scale(slider_high_frame, from_=0, to=65535, orient=tk.HORIZONTAL) #, showvalue=False, resolution = 100)
        slider_high.pack(side=tk.LEFT, fill=tk.X, expand=True)
        slider_high.set(65535)
        #slider_high_max_label = ttk.Label(slider_high_frame, text="100")
        #slider_high_max_label.pack(side=tk.LEFT)
        slider_high.bind("<Motion>", lambda event: update_sliders(slider_high, slider_high_value))

        # Slider low
        slider_low_frame = ttk.Frame(slider_frame)
        slider_low_frame.pack(side=tk.TOP, fill=tk.X, padx=50, pady=0)
        slider_low_value = ttk.Label(slider_low_frame, text="[0]")
        slider_low_value.pack(side=tk.LEFT)
        #slider_low_label = ttk.Label(slider_low_frame, text="0")
        #slider_low_label.pack(side=tk.LEFT)
        slider_low = ttk.Scale(slider_low_frame, from_=0, to=65535, orient=tk.HORIZONTAL) #, showvalue=False, resolution=100)   #, resolution = 100)
        slider_low.pack(side=tk.LEFT, fill=tk.X, expand=True)
        slider_low.set(0)
        #slider_low_max_label = ttk.Label(slider_low_frame, text="100")
        #slider_low_max_label.pack(side=tk.LEFT)
        slider_low.bind("<Motion>", lambda event: update_sliders(slider_low, slider_low_value))
        
        # create canvas 
        canvas = FigureCanvasTkAgg(self.figure, self)

        # create customized toolbar
        class CustomToolbar(NavigationToolbar2Tk):

            def __init__(self, canvas, parent) -> None:
                # list of toolitems to add/modify to the toolbar, format is:
                # (
                #   text, # the text of the button (often not visible to users)
                #   tooltip_text, # the tooltip shown on hover (where possible)
                #   image_file, # name of the image for the button (without the extension)
                #   name_of_method, # name of the method in NavigationToolbar2 to call
                # )
                # this is enforced by MPL lib
                self.toolitems = (
                    ('Home', 'Reset zoom to original view', 'home', 'home'),
                    ('Back', 'Back to previous view', 'back', 'back'),
                    ('Forward', 'Forward to next view', 'forward', 'forward'),
                    (None, None, None, None),
                    ('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan'),
                    ('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'),
                    (None, None, None, None),     #("Configure", "Change parameters", "subplots", "edit_config"),                    
                    ('Save', 'Save the figure', 'filesave', 'save_figure'), 
                )

                super().__init__(canvas = canvas, window = parent, pack_toolbar = True)

        toolbar = CustomToolbar(canvas, self)
        toolbar.update()
        logging.debug('toolbar updated')

        # pack all
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
