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

warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', UserWarning)

class Image(object):

    img_stacked: np.ndarray = np.zeros((2, 8))
    img_names:list[str] = None

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
        
        plt.subplots_adjust(bottom=0.3)

        # show a dummy (zeros) image to start            
        self.show_image(image = Image.img_stacked,
                        #percl = 0,
                        #percu = 99.5,
                        fig = self._figure,
                        ax = self._ax_img,
                        show_colorbar=True, 
                        cmap = 'Grays')

    def update_image(self, val) -> None:
        # Update the image's colormap
        self.image.norm.vmin = val[0]
        self.image.norm.vmax = val[1]
        
        # update image cuts levels
        self.image.set_clim(val[0], val[1])

        # uodate colorbar
        self._colorbar.update_normal(self.image)

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
            
        _img_count: int = 0
        _img_data: list[np.ndarray] = []
        _img_name: list[str] = []
        
        # read new images into memory
        for img_name in path:
            logging.info(f"loading {img_name}...")
            _img_data.append(CCDData.read(img_name, unit = 'adu').data)
            _img_name.append(img_name)
            _img_count += 1

        Image.img_names = _img_name.copy()

        # stack images
        Image.img_stacked = _img_data[0]
        for i in range(1, _img_count):
            Image.img_stacked = np.add(Image.img_stacked, _img_data[i])

        del _img_data

       # collect image stats
        vstd = Image.img_stacked.std()
        vmean = Image.img_stacked.mean()
        _min = Image.img_stacked.min()
        _max = Image.img_stacked.max()
        vmin = vmean - vstd
        vmax = vmean + vstd

        # display image
        logging.info (f"image stats : min = {_min}, max = {_max}, mean = {vmean}, std = {vstd}")
        self.image = self.show_image(image = Image.img_stacked,
                        #percl = 0,
                        #percu = 99.5,
                        fig = self._ax_img.get_figure(),
                        ax = self._ax_img,
                        show_colorbar = show_colorbar, 
                        cmap = self.conf.get_str('window', 'colormap'))

        self._figure.canvas.draw_idle()

    def show_image( self, image,
                    percl = 99.5,
                    percu = None,
                    is_mask = False,
                    figsize= (10, 10),
                    cmap = 'viridis',
                    log = False,
                    clip = True,
                    show_colorbar = True,
                    show_ticks = False,
                    fig = None, ax = None, input_ratio= None) -> AxesImage:
        """
        Show an image in matplotlib with some basic astronomically-appropriat stretching.
        from : https://github.com/astropy/ccd-reduction-and-photometry-guide/blob/main/notebooks/convenience_functions.py

        Parameters
        ----------
        image
            The image to show
        percl : number
            The percentile for the lower edge of the stretch (or both edges if ``percu`` is None)
        percu : number or None
            The percentile for the upper edge of the stretch (or None to use ``percl`` for both)
        figsize : 2-tuple
            The size of the matplotlib figure in inches

        Returns:
            matplotlib.image.AxesImage
        """
        if percu is None:
            percu = percl
            percl = 100 - percl

        if (fig is None and ax is not None) or (fig is not None and ax is None):
            raise ValueError('Must provide both "fig" and "ax" '
                             'if you provide one of them')
        elif fig is None and ax is None:
            if figsize is not None:
                # Rescale the fig size to match the image dimensions, roughly
                image_aspect_ratio = image.shape[0] / image.shape[1]
                figsize = (max(figsize) * image_aspect_ratio, max(figsize))

            fig, ax = plt.subplots(1, 1, figsize=figsize, layout = 'constrained')

        # To preserve details we should *really* downsample correctly and
        # not rely on matplotlib to do it correctly for us (it won't).

        # So, calculate the size of the figure in pixels, block_reduce to
        # roughly that,and display the block reduced image.

        # Thanks, https://stackoverflow.com/questions/29702424/how-to-get-matplotlib-figure-size
        fig_size_pix = fig.get_size_inches() * fig.dpi

        ratio:float = (image.shape // fig_size_pix).max()

        if ratio < 1:
            ratio = 1

        ratio = input_ratio or ratio

        reduced_data = block_reduce(image, ratio)

        if not is_mask:
            # Divide by the square of the ratio to keep the flux the same in the
            # reduced image. We do *not* want to do this for images which are
            # masks, since their values should be zero or one.
             reduced_data = reduced_data / ratio**2

        # Of course, now that we have downsampled, the axis limits are changed to
        # match the smaller image size. Setting the extent will do the trick to
        # change the axis display back to showing the actual extent of the image.
        extent = (0, image.shape[1], 0, image.shape[0])

        if log:
            stretch = aviz.LogStretch()
        else:
            stretch = aviz.LinearStretch()

        norm = aviz.ImageNormalize(reduced_data,
                                   interval=aviz.AsymmetricPercentileInterval(percl, percu),
                                   stretch=stretch, clip=clip)

        if is_mask:
            # The image is a mask in which pixels should be zero or one.
            # block_reduce may have changed some of the values, so reset here.
            reduced_data = reduced_data > 0
            # Set the image scale limits appropriately.
            scale_args = dict(vmin=0, vmax=1)
        else:
            scale_args = dict(norm=norm)

        im = ax.imshow(reduced_data, origin='lower',
                       cmap=cmap,
                       extent=extent,
                       aspect='equal',
                       **scale_args)
        
        ax.set_title(' ')
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        fig.set_label(' ')
        
        if show_colorbar:
            # I haven't a clue why the fraction and pad arguments below work to make
            # the colorbar the same height as the image, but they do....unless the image
            # is wider than it is tall. Sticking with this for now anyway...
            # Thanks: https://stackoverflow.com/a/26720422/3486425
    #        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            self._colorbar: Colorbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
            
            # In case someone in the future wants to improve this:
            # https://joseph-long.com/writing/colorbars/
            # https://stackoverflow.com/a/33505522/3486425
            # https://matplotlib.org/mpl_toolkits/axes_grid/users/overview.html#colorbar-whose-height-or-width-in-sync-with-the-master-axes
            
        if not show_ticks:
            ax.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False)

        return im

    @staticmethod
    def reduce_images() -> None:
        path: list[str] = Image.img_names 
        _img_count: int = 0
        _img_data: list[CCDData] = []
        _img_name: list[str] = []
        
        # read new images into memory
        for img_name in path:
            logging.info(f"loading {img_name}...")
            _img_data.append(CCDData.read(img_name, unit = 'adu'))
            _img_name.append(img_name)
            _img_count += 1

        # read master frames
        #capture_dir: str = os.path.basename(_img_name[0]) 
        capture_dir = '/Users/papa/Documents/ASTRO/CAPTURES/20240826_Void'
        logging.info(f"loading masterbias...")
        master_bias = CCDData.read(capture_dir + '/_offset.fit', unit = u.adu)
        logging.info(f"loading masterdark...")
        master_dark = CCDData.read(capture_dir + '/_dark.fit', unit = u.adu)
        logging.info(f"loading masterflat...")
        master_flat = CCDData.read(capture_dir + '/_flat.fit', unit = u.adu)

        # stack images
        for i in range(0, _img_count):
            logging.info(f"processing {_img_name[i]}...")
            _img_data[i] = ccd_process(ccd = _img_data[i], 
                oscan = None, 
                gain_corrected = True, 
                trim = None, 
                error = False,
#                gain = camera_electronic_gain*u.electron/u.adu ,
#                readnoise = camera_readout_noise*u.electron,
                master_bias = master_bias,
                dark_frame = master_dark,
                master_flat = master_flat,
                exposure_key = 'EXPTIME',
                exposure_unit = u.second,
                dark_scale = True)       

        logging.info(f"summing images...")
        Image.img_stacked = _img_data[0].data
        for i in range(1, _img_count):
            np.add(Image.img_stacked, _img_data[i])

        del _img_data


