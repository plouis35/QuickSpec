import logging
import numpy as np
import warnings

import matplotlib.pyplot as plt
from matplotlib.image import AxesImage
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.colorbar import Colorbar

from astropy import visualization as aviz
from astropy.nddata.blocks import block_reduce

from astropy import units as u
from astropy.nddata import CCDData, NDData, StdDevUncertainty
from astropy.stats import mad_std
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning

from ccdproc import Combiner, combine, subtract_bias, subtract_dark, flat_correct
from ccdproc import trim_image, Combiner, ccd_process, cosmicray_median

import tkinter as tk
from tkinter import ttk

from tkinter.filedialog import askopenfilenames

from app.config import Config
from img.img_utils import reduce_images

warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', UserWarning)

class Image(object):

    img_stacked: np.ndarray = None
    img_names:tuple[str, ...] = None
    img_colorbar: Colorbar = None

    def __init__(self, axe_img: Axes, axe_spc: Axes, frame: ttk.Frame) -> None:
        self.conf: Config = Config()
        self._ax_img: Axes = axe_img
        self._ax_spc: Axes = axe_spc
        self._figure: Figure = axe_img.get_figure()
        self.image: AxesImage = None
        Image.img_stacked = np.zeros((2, 8))
        
    # create sliders
        def update_slider(event):
            slider: ttk.Scale = event.widget
            if slider == self.slider_high:
                self.update_image(None, slider.get())
            if slider == self.slider_low:
                self.update_image(slider.get() , None)

        slider_frame = ttk.Frame(frame)
        slider_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        slider_high_frame = ttk.Frame(slider_frame)
        slider_high_frame.pack(side=tk.TOP, fill=tk.X, padx=50, pady=0)
        self.slider_high_value = ttk.Label(slider_high_frame, text="0")
        self.slider_high_value.pack(side=tk.LEFT)
        #slider_high_label = ttk.Label(slider_high_frame, text="0")
        #slider_high_label.pack(side=tk.LEFT)
        self.slider_high = ttk.Scale(slider_high_frame, from_=0, to=0, orient=tk.HORIZONTAL) #, showvalue=False, resolution = 100)
        self.slider_high.pack(side=tk.LEFT, fill=tk.X, expand=True)
        #self.slider_high.set(65535)
        #slider_high_max_label = ttk.Label(slider_high_frame, text="100")
        #slider_high_max_label.pack(side=tk.LEFT)
        self.slider_high.bind("<ButtonRelease>", update_slider) 

        slider_low_frame = ttk.Frame(slider_frame)
        slider_low_frame.pack(side=tk.TOP, fill=tk.X, padx=50, pady=0)
        self.slider_low_value = ttk.Label(slider_low_frame, text="0")
        self.slider_low_value.pack(side=tk.LEFT)
        #slider_low_label = ttk.Label(slider_low_frame, text="0")
        #slider_low_label.pack(side=tk.LEFT)
        self.slider_low = ttk.Scale(slider_low_frame, from_=0, to=0, orient=tk.HORIZONTAL) #, showvalue=False, resolution=100)   #, resolution = 100)
        self.slider_low.pack(side=tk.LEFT, fill=tk.X, expand=True)
        #self.slider_low.set(0)
        #slider_low_max_label = ttk.Label(slider_low_frame, text="100")
        #slider_low_max_label.pack(side=tk.LEFT)
        self.slider_low.bind("<ButtonRelease>", update_slider) 
        
        # show a dummy (zeros) image to start            
        self.image = self.show_image(image = Image.img_stacked,
                        fig_img = self._figure,
                        ax_img = self._ax_img,
                        show_colorbar=True, 
                        cmap = self.conf.get_str('window', 'colormap')) # type: ignore

    def clear_image(self) -> None:
        self._figure.clear()
        #self._figure.clf()
        #self._ax_img.clear()
        self._figure.canvas.draw()

    def update_image(self, low_cut: float | None, high_cut: float | None) -> None:
        # Update the image's colormap and cuts
        try:
            if (low_cut is not None) and (self.image.norm.vmax > low_cut): 
                self.image.norm.vmin = low_cut 
                self.slider_low_value.config(text=int(low_cut))
            if (high_cut is not None) and  (self.image.norm.vmin < high_cut): 
                self.image.norm.vmax = high_cut 
                self.slider_high_value.config(text=int(high_cut))

        except Exception as e:
            logging.error({e})
        
        self._figure.canvas.draw()

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
        #self._ax_img.clear()
        #self.image.remove()

        Image.img_names = path
        try:
            _imgs: CCDData | None = reduce_images(images=Image.img_names, preprocess=False)
        except Exception as e:
            logging.error(f"{e}")
            return 

        Image.img_stacked = _imgs.data

        # collect image stats
        vstd = Image.img_stacked.std()
        vmean = Image.img_stacked.mean()
        _min = Image.img_stacked.min()
        _max = Image.img_stacked.max()

        # update sliders positions
        self.slider_low.config(from_=_min)
        self.slider_low.config(to=_max)

        self.slider_high.config(from_=_min)
        self.slider_high.config(to=_max)

        nb_sigma = 5
        low_cut = vmean - (nb_sigma * vstd)
        high_cut = vmean + (nb_sigma * vstd)
        self.slider_low.set(low_cut)
        self.slider_high.set(high_cut)

        # display image
        logging.info (f"image stats : min = {_min}, max = {_max}, mean = {vmean}, std = {vstd}")
        self.image = self.show_image(image = Image.img_stacked,
                        fig_img = self._figure,
                        ax_img = self._ax_img,
                        show_colorbar = True, 
                        cmap = self.conf.get_str('window', 'colormap'))
        
        self.update_image(low_cut, high_cut)
        #self._figure.canvas.draw_idle()

    def show_image( self, image,
                    cmap: str,
                    show_colorbar: bool,
                    fig_img: Figure, 
                    ax_img: Axes) -> AxesImage:
        
        ax_img.axis('off')
        ax_img.set_yscale('linear')

        #norm = aviz.ImageNormalize(image,
         #                       interval=aviz.AsymmetricPercentileInterval(5, 95),
          #                      stretch=aviz.LinearStretch(), clip=True)
        #scale_args = dict(norm=norm)

        img: AxesImage = ax_img.imshow(
            image, 
            origin = 'lower', 
            #vmin = self.slider.val[0],   #image.min(), 
            #vmax = self.slider.val[1], #image.max(), 
            interpolation='none',
            aspect = 'equal',
            cmap = cmap,
            #**scale_args
        )

        def format_coord(x,y):
            return "x={:.0f}, y={:.0f} ->".format(x,y)
                    
        ax_img.format_coord=format_coord

        if show_colorbar:
            if Image.img_colorbar is not None:
                Image.img_colorbar.remove()
            
            Image.img_colorbar = fig_img.colorbar(img, ax = ax_img, location='right', shrink=0.6)

        return img
