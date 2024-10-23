import logging
import warnings
import numpy as np
import tkinter as tk
from tkinter import ttk
import tkinter.messagebox as msg
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

from app.config import Config
from app.logger import LogHandler
from app.toolbar import CustomToolbar  
from img.image import spec2d
from spc.spectrum import spec1d


class Application(tk.Tk):
    def __init__(self, version):
        super().__init__()
        self.title(f"QuickSpec")
        
        # read config 
        self.conf = Config().config

        # catch useless warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)

        # set theme & geometry
        plt.rcParams['figure.constrained_layout.use'] = True
        plt.style.use(self.conf['window']['theme'])
        self.geometry(self.conf['window']['geometry'])
        self.wm_state(self.conf['window']['state'])

        # create & configure logger frame
        log_frame = tk.Frame(self, bg='gray20')
        log_frame.pack(side=tk.BOTTOM, fill=tk.X, expand=False)
        log_text = tk.Text(log_frame, state='disabled', height=self.conf['logger']['nb_lines'], bg='black', fg='white')
        log_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        log_text.configure(font = (self.conf['logger']['font']))

        scrollbar = ttk.Scrollbar(log_frame, command=log_text.yview)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        log_text['yscrollcommand'] = scrollbar.set

        self.logger = logging.getLogger()
        self.logger.setLevel(self.conf['logger']['level'])
        handler = LogHandler(log_text)
        self.logger.addHandler(handler)
        logging.info('logger started')

        # pre-create toolbar frame on top
        toolbar_frame = tk.Frame(self) #, height = 80, bg='white')
        toolbar_frame.pack(side=tk.TOP, fill=tk.X)
        logging.debug('toolbar created')
 
        # create a single frame for both image & spectrum horizontaly packed
        figure = Figure(figsize=(5, 4))
        ax_img = figure.add_subplot(211)
        ax_spc = figure.add_subplot(212)
        logging.debug('figure created')

        # create image widget
        self._image = spec2d(ax_img)
        logging.debug('image frame created')

        # create spectrum widget
        self._spectrum = spec1d(ax_spc)
        logging.debug('spectrum frame created')

        # create canvas for customized toolbar
        canvas = FigureCanvasTkAgg(figure, self)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        logging.debug('canvas created')

        # and finally update toolbar
        toolbar = CustomToolbar(canvas, toolbar_frame)
        toolbar.update()
        logging.debug('toolbar updated')


