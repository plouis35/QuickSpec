import logging
import threading
import time
import numpy as np

import tkinter as tk
from tkinter import ttk, Menu

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.image import AxesImage
from matplotlib.widgets import Button

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends._backend_tk import NavigationToolbar2Tk

from app.logger import LogHandler
from app.config import Config
from img.img_utils import ImgTools
from spc.calibration import Calibration
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
        #plt.rcParams['font.size'] = 5
        plt.style.use(self.conf.get_str('window', 'theme'))        

        # create a single figure for both image & spectrum horizontaly packed
        self.figure, (self.axe_img, self.axe_spc) = plt.subplots(2, 1, figsize=(8, 6)) #, num="QuickSpec")
        plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.99)      
        logging.debug('figure created')

        # intialize imgTools
        ImgTools(self.axe_img, self.axe_spc)

        # create load button
        #self._load_image_ax = self.figure.add_axes(rect=(0.08, 0.54, 0.06, 0.05)) # (left, bottom, width, height
        #self._bt_load_image = Button(self._load_image_ax, 'Load', color='k') #, hovercolor='grey')
        #self._bt_load_image.on_clicked(ImgTools.open_image)

        # initialize calibration
        _spcTools = Calibration(self.axe_img, self.axe_spc)
        
        # create menus
        class MenuBar(tk.Menu):
            def __init__(self, root) -> None:
                super().__init__(root)
                show_all = tk.BooleanVar()
                show_all.set(True)
                show_done = tk.BooleanVar()
                show_not_done = tk.BooleanVar()

                self.fileMenu = tk.Menu(self, tearoff=False)
                self.fileMenu.add_command(label ='Open...', command = ImgTools.open_image)
                self.fileMenu.entryconfigure('Open...', accelerator='Command+O')
                self.fileMenu.add_command(label ='Open recent...', command = None) 
                self.fileMenu.add_separator()
                self.fileMenu.add_command(label ='Exit', command = root.destroy)

                self.processMenu = tk.Menu(self, tearoff=False)
                self.processMenu.add_command(label ='Run all', command = _spcTools.do_calibration)
                self.processMenu.add_separator()
                self.processMenu.add_checkbutton(label="Watch mode", onvalue=1, offvalue=0, variable=show_all)
                self.processMenu.add_separator()
                self.processMenu.add_command(label ='Trace', command = None)
                self.processMenu.add_command(label ='Calibrate', command = None)
                self.processMenu.add_command(label ='Apply response', command = None)
                
                self.viewMenu = tk.Menu(self, tearoff=False)
                self.viewMenu.add_checkbutton(label ='Show common lines',  onvalue=1, offvalue=0, variable=show_all)
                self.viewMenu.add_checkbutton(label ='Show atm lines',  onvalue=1, offvalue=0, variable=show_all)

                self.aboutMenu = tk.Menu(self, tearoff=False)
                self.aboutMenu.add_command(label =f"About...", command = None)

                self.add_cascade(label='File',menu=self.fileMenu)
                self.add_cascade(label='Process',menu=self.processMenu)
                self.add_cascade(label='View',menu=self.viewMenu)
                self.add_cascade(label='About',menu=self.aboutMenu)

        #container = tk.Frame(self)
        #container.pack(fill=tk.BOTH,expand=True)
        root_menu=MenuBar(self)
        self.config(menu=root_menu)

        # create run button
        #self._run_calib_ax: Axes = self.figure.add_axes(rect=(0.15, 0.54, 0.06, 0.05)) # (left, bottom, width, height
        #self._bt_run = Button(self._run_calib_ax, 'Run', color='k') #, hovercolor='grey')
        #self._bt_run.on_clicked(_spcTools.do_calibration)

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
        canvas.draw_idle()

