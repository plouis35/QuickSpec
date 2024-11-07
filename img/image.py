import logging
import numpy as np
import warnings
import os

import matplotlib.pyplot as plt
from matplotlib.image import AxesImage
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.colorbar import Colorbar

from astropy import visualization as aviz
from astropy.nddata.blocks import block_reduce
from astropy.nddata.utils import Cutout2D

from astropy import units as u
from astropy.nddata import CCDData, NDData, StdDevUncertainty
from astropy.stats import mad_std
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning

from ccdproc import Combiner, combine, subtract_bias, subtract_dark, flat_correct
from ccdproc import trim_image, Combiner, ccd_process, cosmicray_median

import tkinter as tk
from tkinter import ttk
#import customtkinter as ctk

from tkinter.filedialog import askopenfilenames

from app.config import Config
from img.img_utils import show_image, reduce_images


warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', UserWarning)

class Image(object):

    img_stacked: np.ndarray = np.zeros((2, 8))
    img_names:tuple[str, ...] = None

    def __init__(self, axe_img: Axes, axe_spc: Axes) -> None:
        self.conf: Config = Config()
        self._colorbar: Colorbar = None
        self._ax_img: Axes = axe_img
        self._ax_spc: Axes = axe_spc
        self._figure: Figure = axe_img.get_figure()
        self.image: AxesImage = None

        # create cuts range slider
        def format_coord(x,y) -> str:
            return f'x={x:.0f}, y={y:.0f}'

        self._ax_img.format_coord = format_coord
        
        #plt.subplots_adjust(bottom=0.3)

        # show a dummy (zeros) image to start            
        show_image(image = Image.img_stacked,
                        fig = self._figure,
                        ax = self._ax_img,
                        show_colorbar=False, 
                        cmap = 'Grays')

    def update_image(self, val) -> None:
        # Update the image's colormap
        self.image.norm.vmin = val[0]
        self.image.norm.vmax = val[1]
        
        # update image cuts levels
        self.image.set_clim(val[0], val[1])

        # uodate colorbar
        #self._colorbar.update_normal(self.image)

        # Redraw the figure to ensure it updates
        self._figure.canvas.draw_idle()


    def open_image(self) -> None:
        # create openfile dialog
        self.path: tuple[str, ...] | Literal[''] = askopenfilenames(title='Select image(s) or a directory for watch mode',
                            initialdir = self.conf.get_str('files', 'initial_directory'),
                            defaultextension = self.conf.get_str('files', 'file_types')
                            )
        #path = (['albireo-10.fit'])
        if self.path == '': 
            return 
        else: 
            self.load_image(self.path)
            return
        
    def load_image(self, path: tuple[str, ...]) -> None:
        # cleanup previous images
        self._ax_img.clear()

        # do not recreate colorbar
        if (self._colorbar is None):
            show_colorbar = True
        else:
            show_colorbar = False

        """"    
        _img_count: int = 0
        _img_data: list[np.ndarray] = []
        _img_names: list[str] = []
        
        # read new images into memory
        for img_name in path:
            logging.info(f"loading {img_name}...")
            _img_data.append(CCDData.read(img_name, unit = 'adu').data)
            _img_names.append(img_name)
            _img_count += 1

        Image.img_names = _img_names.copy()

        # stack images
        Image.img_stacked = _img_data[0]
        for i in range(1, _img_count):
            Image.img_stacked = np.add(Image.img_stacked, _img_data[i])

        del _img_data
        """

        Image.img_names = path
        Image.img_stacked = reduce_images(images=Image.img_names, preprocess=False).data

       # collect image stats
        vstd = Image.img_stacked.std()
        vmean = Image.img_stacked.mean()
        _min = Image.img_stacked.min()
        _max = Image.img_stacked.max()
        vmin = vmean - vstd
        vmax = vmean + vstd

        # display image
        logging.info (f"image stats : min = {_min}, max = {_max}, mean = {vmean}, std = {vstd}")
        show_image(image = Image.img_stacked,
                        fig = self._ax_img.get_figure(),
                        ax = self._ax_img,
                        show_colorbar = False, #show_colorbar, 
                        cmap = self.conf.get_str('window', 'colormap'))

        self._figure.canvas.draw_idle()

