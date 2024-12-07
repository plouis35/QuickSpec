import logging
import numpy as np
import warnings
from pathlib import Path
import contextlib

import numpy as np

import tkinter as tk
from tkinter import ttk

import matplotlib.pyplot as plt
from matplotlib.image import AxesImage
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.colorbar import Colorbar
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends._backend_tk import NavigationToolbar2Tk

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
from img.img_utils import Images, ImagesCombiner

warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', UserWarning)

class Image(object):
    def __init__(self, img_frame: ttk.Frame, bt_frame: ttk.Frame) -> None: #, axe_img: Axes, axe_spc: Axes, frame: ttk.Frame) -> None:
        self.conf: Config = Config()

        # assign image instance variables
        self.image: AxesImage = None
        self.img_stacked: CCDData = CCDData(np.zeros((2,8)), unit=u.adu)
        self.img_names:list[str] = []
        self.img_count = 0
        self.img_colorbar: Colorbar = None
        self.img_combiner: ImagesCombiner = None
        self.img_cmap: str = 'grey'

        # create figure and axe for image
        self.img_figure = Figure(figsize=(10, 3))
        self.img_axe = self.img_figure.add_subplot(111)

        # create canvas to draw to
        self.img_canvas = FigureCanvasTkAgg(self.img_figure, img_frame) #img_frame)
        self.img_canvas.draw()

        # create customized toolbar
        img_toolbar = CustomImgToolbar(self.img_canvas, img_frame)

        # add new buttons
        self.clear_button = ttk.Button(img_toolbar, text="Clear", command=self.clear_image)
        self.clear_button.pack(side=tk.LEFT, padx=5, pady=0)

        _cmap_options = ["grey", "inferno", "magma", "plasma"]
        bt_cmap_default = "grey"
        _var = tk.StringVar(value=bt_cmap_default)

        def cb_cmap_option(cmap: tk.StringVar) -> None:
            self.image.set_cmap(str(cmap))
            self.img_figure.canvas.draw()
            self.img_cmap = str(cmap)
            logging.info(f"colormap changed to {cmap}")
            
        bt_cmap = ttk.OptionMenu(img_toolbar, _var, bt_cmap_default, *(_cmap_options), command = cb_cmap_option)
        bt_cmap.pack(side=tk.LEFT, padx=5, pady=0)

        # pack buttons
        img_toolbar.update()
        img_toolbar.pack(side=tk.TOP, fill=tk.X)
        self.img_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True) #side=tk.TOP, 

        # now create sliders
        slider_frame = ttk.Frame(bt_frame)
        slider_frame.pack(side=tk.RIGHT, fill=tk.BOTH, padx=35, pady=0, expand=True)

        slider_high_frame = ttk.Frame(slider_frame)
        slider_high_frame.pack(side=tk.TOP, fill=tk.X, padx=5, pady=0)
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
        slider_low_frame.pack(side=tk.TOP, fill=tk.X, padx=5, pady=0)
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
        
        # and finally show a dummy (zeros) image         
        self.image = self.show_image(image = self.img_stacked,
                        fig_img = self.img_figure,
                        ax_img = self.img_axe,
                        show_colorbar=True, 
                        #cmap = 'grey'
                        )

    def update_slider(self, event) -> None:
        slider: ttk.Scale = event.widget
        if slider == self.slider_high:
            self.update_image(None, slider.get())
        if slider == self.slider_low:
            self.update_image(slider.get() , None)

    def clear_image(self) -> None:
        self.img_stacked: CCDData = CCDData(np.zeros((2,8)), unit=u.adu)
        self.img_names:list[str] = []
        self.img_count = 0
        self.img_combiner: ImagesCombiner = None
        #self.image.set_data(self.img_stacked)
        #self.img_axe.clear()
        self.image = self.show_image(image = self.img_stacked,
                fig_img = self.img_figure,
                ax_img = self.img_axe,
                show_colorbar=False, 
                #cmap = 'grey'
                )

        self.img_figure.canvas.draw_idle()

    def stats_image(self) -> tuple[float, float, float, float]:
        v_std = self.img_stacked.data.std()
        v_mean = self.img_stacked.data.mean()
        v_min = self.img_stacked.data.min()
        v_max = self.img_stacked.data.max()
        return v_std, v_mean, v_min, v_max

    def update_image(self, low_cut: float | None = None, high_cut: float | None = None) -> None:
        # collect image stats
        if (low_cut is None) and (high_cut is None):
            v_std, v_mean, v_min, v_max = self.stats_image()

            # update sliders positions
            nb_sigma = 6 #self.conf.get_int('display', 'contrast_level')
            low_cut = v_mean - (nb_sigma * v_std)
            high_cut = v_mean + (nb_sigma * v_std)   

            min_cut = low_cut - (2 * nb_sigma * v_std)
            max_cut = high_cut + (8 * nb_sigma * v_std)
            
            self.slider_low.config(from_=min_cut)
            self.slider_low.config(to=max_cut)

            self.slider_high.config(from_=min_cut)
            self.slider_high.config(to=max_cut)

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
        
        self.img_figure.canvas.draw_idle()


    def load_images(self, path: list[str]) -> bool:
        _img_reduced: CCDData
        _img_combiner: ImagesCombiner

        try:
            _img_combiner = Images.from_fits(imgs=path).y_crop(y_ratio=self.conf.get_float('pre_processing','crop_auto'))
            _img_reduced = _img_combiner.sum()

        except Exception as e:
            logging.error(f"{e}")
            return False
        
        self.img_stacked = _img_reduced.copy()    
        self.img_combiner = _img_combiner
    
        # display image
        self.image = self.show_image(image = self.img_stacked.data,
                        fig_img = self.img_figure,
                        ax_img = self.img_axe,
                        show_colorbar = True, 
                        #cmap = self.conf.get_str('display', 'colormap')
                        )
        
        v_std, v_mean, v_min, v_max = self.stats_image()
        logging.info (f"image stats: min={v_min}, max={v_max}, mean={v_mean}, std={v_std}")

        self.update_image()

        return True

    
    def reduce_images(self) -> bool:
        _img_reduced: CCDData | None
        _img_combiner: ImagesCombiner = self.img_combiner

        if self.img_combiner is None:
            logging.error('please load some images before reducing')
            return False
        
        try:
            _img_reduced = _img_combiner.reduce_images()
            
        except Exception as e:
            logging.error(f"{e}")
            return False
        
        if _img_reduced is not None:
            self.img_stacked = _img_reduced.copy()    
            #self.image.set_data(self.img_stacked.data)
            #self.update_image()
        
            v_std, v_mean, v_min, v_max = self.stats_image()
            logging.info (f"image stats: min={v_min}, max={v_max}, mean={v_mean}, std={v_std}")
        else:
            logging.error("unable to reduce images")
            return False
    
        # display image
        self.image = self.show_image(image = self.img_stacked,
                        fig_img = self.img_figure,
                        ax_img = self.img_axe,
                        show_colorbar = True, 
                        #cmap = 'grey' #cmap = self.conf.get_str('display', 'colormap')
                        )
        
        self.update_image()
        return True

    def show_image( self, image,
                    #cmap: str,
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
            cmap = self.img_cmap,
            #**scale_args
        )

        if show_colorbar:
            if self.img_colorbar is not None:
                try:
                    self.img_colorbar.remove()
                except:
                    pass
            
            self.img_colorbar = fig_img.colorbar(img, ax = ax_img, location='right', shrink=0.6)

        return img


class CustomImgToolbar(NavigationToolbar2Tk):
    def __init__(self, canvas, parent) -> None:
        # list of toolitems to add/modify to the toolbar, format is:
        # (
        #   text, # the text of the button (often not visible to users)
        #   tooltip_text, # the tooltip shown on hover (where possible)
        #   image_file, # name of the image for the button (without the extension)
        #   name_of_method, # name of the method in NavigationToolbar2 to call
        # )
        self.toolitems = (
            ('Home', 'Reset zoom to original view', 'home_large', 'home'),
            ('Back', 'Back to previous view', 'back_large', 'back'),
            ('Forward', 'Forward to next view', 'forward_large', 'forward'),
            (None, None, None, None),
            ('Pan', 'Pan axes with left mouse, zoom with right', 'move_large', 'pan'),
            ('Zoom', 'Zoom to rectangle', 'zoom_to_rect_large', 'zoom'),
            (None, None, None, None),
            ('Save', 'Save the figure', 'filesave_large', 'save_figure'), 
        )

        super().__init__(canvas = canvas, window = parent, pack_toolbar = True)
