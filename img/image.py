import logging
import numpy as np
import warnings
from pathlib import Path

import numpy as np

import tkinter as tk
from tkinter import ttk

import matplotlib.pyplot as plt
from matplotlib.image import AxesImage
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.colorbar import Colorbar


from astropy import units as u
from astropy.nddata import CCDData, NDData, StdDevUncertainty
from astropy.stats import mad_std
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
from astropy import visualization as aviz
from astropy.nddata.blocks import block_reduce

from ccdproc import Combiner, combine, subtract_bias, subtract_dark, flat_correct
from ccdproc import trim_image, Combiner, ccd_process, cosmicray_median

from app.config import Config
from img.img_utils import Images

warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', UserWarning)

class Image(object):

    img_stacked: np.ndarray = None
    img_names:list[str] = []
    img_count = 0
    img_colorbar: Colorbar = None

    def __init__(self, axe_img: Axes, axe_spc: Axes, frame: ttk.Frame) -> None:
        self.conf: Config = Config()
        self._ax_img: Axes = axe_img
        self._ax_spc: Axes = axe_spc
        self._figure: Figure = axe_img.get_figure()
        self.image: AxesImage = None
        Image.img_stacked = np.zeros((2, 8))
        
        # create sliders
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
        self.slider_high.bind("<ButtonRelease>", self.update_slider) 

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
        self.slider_low.bind("<ButtonRelease>", self.update_slider) 
        
        # show a dummy (zeros) image to start            
        self.image = self.show_image(image = Image.img_stacked,
                        fig_img = self._figure,
                        ax_img = self._ax_img,
                        show_colorbar=True, 
                        cmap = 'grey')  # #self.conf.get_str('display', 'colormap')) # type: ignore

    def update_slider(self, event) -> None:
        slider: ttk.Scale = event.widget
        if slider == self.slider_high:
            self.update_image(None, slider.get())
        if slider == self.slider_low:
            self.update_image(slider.get() , None)

    def clear_image(self) -> None:
        self._figure.clear()
        #self._figure.clf()
        #self._ax_img.clear()
        self._figure.canvas.draw()

    def stats_image(self) -> tuple[float, float, float, float]:
        v_std = Image.img_stacked.std()
        v_mean = Image.img_stacked.mean()
        v_min = Image.img_stacked.min()
        v_max = Image.img_stacked.max()
        return v_std, v_mean, v_min, v_max

    def update_image(self, low_cut: float | None = None, high_cut: float | None = None) -> None:
        # collect image stats
        if (low_cut is None) and (high_cut is None):
            v_std, v_mean, v_min, v_max = self.stats_image()

            # update sliders positions
            nb_sigma = self.conf.get_int('display', 'contrast_level')
            low_cut = v_mean - (nb_sigma * v_std)
            high_cut = v_mean + (nb_sigma * v_std)   

            self.slider_low.config(from_=v_min / nb_sigma)
            self.slider_low.config(to=v_max / nb_sigma)

            self.slider_high.config(from_=v_min / nb_sigma)
            self.slider_high.config(to=v_max / nb_sigma)

        # Update the image's colormap and cuts
        try:
            if (low_cut is not None) : #and (self.image.norm.vmax > low_cut): 
                self.image.norm.vmin = low_cut 
                self.slider_low_value.config(text=int(low_cut))
                self.slider_low.set(low_cut)
            if (high_cut is not None) : #and  (self.image.norm.vmin < high_cut): 
                self.image.norm.vmax = high_cut 
                self.slider_high_value.config(text=int(high_cut))
                self.slider_high.set(high_cut)

        except Exception as e:
            logging.error({e})
        

        self._figure.canvas.draw_idle()


    def load_image(self, path: list[str]) -> None:
        Image.img_names = path

        img_count = 0
        img_data = []

        # read images into memory
        for img_name in path:
            fit_data: CCDData = CCDData.read(img_name, unit=u.adu)
            if fit_data.ndim == 2:
                # this is a fit image - load it
                img_data.append(fit_data.data)
                logging.info(f"{img_name} loaded")
                img_count += 1
            else:
                # not supported fit format
                logging.error(f"{img_name} is not a supported fit format (naxis := 2)")

        # stack images
        Image.img_stacked = Images.reduce_images_numpy(img_data, False)

        #logging.info(f"{low_cut=}, {high_cut=}")

        # display image
        self.image = self.show_image(image = Image.img_stacked,
                        fig_img = self._figure,
                        ax_img = self._ax_img,
                        show_colorbar = True, 
                        cmap = self.conf.get_str('display', 'colormap'))
        
        self.update_image()
        self._ax_img.set_title(f"{Path(path[0]).stem}...", fontsize=10, loc='left') 

        v_std, v_mean, v_min, v_max = self.stats_image()
        logging.info (f"image stats: min={v_min}, max={v_max}, mean={v_mean}, std={v_std}")


    def show_image( self, image,
                    cmap: str,
                    show_colorbar: bool,
                    fig_img: Figure, 
                    ax_img: Axes) -> AxesImage:
        
        ax_img.axis('off')
        ax_img.set_yscale('linear')

        #norm = aviz.ImageNormalize(image,
         #                       interval=aviz.AsymmetricPercentileInterval(percl, percu),
          #                      stretch=aviz.LinearStretch(), clip=True)
        #scale_args = dict(norm=norm)

        img: AxesImage = ax_img.imshow(
            image, 
            origin = 'lower', 
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
    
    def reduce_images(self) -> None:
        _img_reduced = Images.reduce_images_ccdproc(Image.img_names, preprocess=True)

        if _img_reduced is not None:
            Image.img_stacked = _img_reduced.data
            self.image.set_data(Image.img_stacked.data)
            self.update_image()
        
            v_std, v_mean, v_min, v_max = self.stats_image()
            logging.info (f"image stats: min={v_min}, max={v_max}, mean={v_mean}, std={v_std}")
        else:
            logging.error("unable to reduce images")



