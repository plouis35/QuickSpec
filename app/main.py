import logging
import threading
import time
import numpy as np

import tkinter as tk
from tkinter import ttk, Menu
import tkinter.font as tkFont
#import customtkinter as ctk

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.image import AxesImage
from matplotlib.widgets import Button
from matplotlib.widgets import RangeSlider


from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends._backend_tk import NavigationToolbar2Tk

from app.logger import LogHandler
from app.config import Config
from img.img_tools import ImgTools
from spc.spc_tools import Calibration
from app.os_utils import OSUtils


class Application(tk.Tk):

    def __init__(self, app_name: str, app_version: str) -> None:
        super().__init__()

        self.conf = Config()
        LogHandler().set()
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
        default_font = tkFont.nametofont("TkDefaultFont")
        default_font.configure(size=18)
        self.option_add("*Font", default_font)
        plt.style.use(self.conf.get_str('window', 'theme'))        

        # create a single figure for both image & spectrum horizontaly packed
        self.figure, (self.axe_img, self.axe_spc) = plt.subplots(2, 1, figsize=(12, 8)) #, num="QuickSpec")
        #plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.50)   
        logging.debug('figure created')

        # intialize img tools
        ImgTools(self.axe_img, self.axe_spc)

        # initialize spec tools
        _spcTools = Calibration(self.axe_img, self.axe_spc)
        
        # create buttons
        self._button_frame = ttk.Frame(self)
        self._button_frame.pack(side='top', fill='x')

        self.__btn_load = ttk.Button(self._button_frame, text="Load", command=ImgTools.open_image)
        self.__btn_load.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)

        self.__btn_run = ttk.Button(self._button_frame, text="Run", command=_spcTools.do_calibration)
        self.__btn_run.grid(row=0, column=1, sticky="nsew", padx=10, pady=10)

        self._button_frame.grid_columnconfigure(0, weight=1)
        self._button_frame.grid_columnconfigure(1, weight=1)

        logging.debug('buttons created')

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
        canvas.get_tk_widget().pack(side = tk.TOP, fill = tk.BOTH, expand = True)
        canvas.draw() #_idle()

