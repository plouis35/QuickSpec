import logging
import warnings
import numpy as np
import tkinter as tk
from tkinter import ttk
import tkinter.messagebox as msg
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.figure import Figure

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg #, NavigationToolbar2Tk
from matplotlib.backends._backend_tk import NavigationToolbar2Tk

from app.config import Config
from app.logger import LogHandler
#from app.toolbar import CustomToolbar  
from img.image import spec2d
from spc.spectrum import spec1d

class Application(tk.Tk):
    def __init__(self, version: str) -> None:
        super().__init__()
        self.title(f"QuickSpec")
        
        # read config 
        self.conf = Config()

        # catch useless warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)

        # set theme & geometry
        plt.rcParams['figure.constrained_layout.use'] = True
        plt.style.use(self.conf.get_config('window', 'theme'))
        self.geometry(self.conf.get_config('window', 'geometry'))
        self.wm_state(self.conf.get_config('window', 'state'))

        # create & configure logger frame
        log_frame = tk.Frame(self, bg='gray20')
        log_frame.pack(side=tk.BOTTOM, fill=tk.X, expand=False)
        log_text = tk.Text(log_frame, state='disabled', height=self.conf.get_config('logger', 'nb_lines'), bg='black', fg='white')
        log_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        log_text.configure(font = (self.conf.get_config('logger', 'font')))

        scrollbar = ttk.Scrollbar(log_frame, command=log_text.yview)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        log_text['yscrollcommand'] = scrollbar.set

        self.logger = logging.getLogger()
        self.logger.setLevel(self.conf.get_config('logger', 'level'))
        handler = LogHandler(log_text)
        self.logger.addHandler(handler)
        logging.info('logger started')
 
        # create a single frame for both image & spectrum horizontaly packed
        figure = Figure(figsize=(5, 4))
        ax_img = figure.add_subplot(211)        # nrows/ncols/index
        ax_spc = figure.add_subplot(212)
        logging.debug('figure created')

        # create image widget
        self._image = spec2d(ax_img)
        logging.debug('image frame created')

        # create spectrum widget
        self._spectrum = spec1d(ax_spc, self._image)
        logging.debug('spectrum frame created')
        
        # create canvas 
        canvas = FigureCanvasTkAgg(figure, self)

        # create toolbar
        toolbar = NavigationToolbar2Tk(canvas, self, pack_toolbar = False)
        toolbar.pack(side = tk.TOP, fill = tk.X)

        # pack all
        canvas.get_tk_widget().pack(side = tk.TOP, fill = tk.BOTH, expand = True)
        canvas.draw_idle()
