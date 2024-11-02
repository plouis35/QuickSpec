import logging
import numpy as np
import psutil
import warnings

import matplotlib.pyplot as plt
from matplotlib.image import AxesImage
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.colorbar import Colorbar
from matplotlib.widgets import RangeSlider


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
from tkinter.filedialog import askopenfilenames

#from tkinter import Tk, filedialog, Frame, Button, Canvas


from app.config import Config

warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', UserWarning)

class ImgTools(object):

    def __init__(self, axe_img: Axes, axe_spc: Axes) -> None:
        ImgTools.conf: Config = Config()
        ImgTools._colorbar: Colorbar = None
        ImgTools._ax_img: Axes = axe_img
        ImgTools._ax_spc: Axes = axe_spc
        ImgTools._figure: Figure = axe_spc.get_figure()

        # public variables
        ImgTools.img_stacked: np.ndarray = np.zeros((2, 8))

        # create cuts range slider
        def format_coord(x,y) -> str:
            return f'x={x:.0f}, y={y:.0f}'

        ImgTools._ax_img.format_coord = format_coord
        
        ImgTools._slider_ax: Axes = ImgTools._figure.add_axes((0.082, 0.97, 0.55, 0.01))  # (left, bottom, width, height)
        ImgTools._slider: RangeSlider = RangeSlider(ImgTools._slider_ax, "Cuts: ",
                                   orientation = 'horizontal',
                                   valstep = 10,
                                   dragging = True,
                                   valinit = (0, 65535),
                                   valmin = 0,
                                   valmax = 65535)
        logging.debug('slider frame created')

        # show a dummy (zeros) image to start            
        ImgTools.show_image(image = ImgTools.img_stacked,
                        #percl = 0,
                        #percu = 99.5,
                        fig = ImgTools._figure,
                        ax = ImgTools._ax_img,
                        show_colorbar=True, 
                        cmap = 'Grays')

    @staticmethod
    def update_image(val) -> None:
        # Update the image's colormap
        ImgTools._img.norm.vmin = val[0]
        ImgTools._img.norm.vmax = val[1]
        
        # update image cuts levels
        ImgTools._img.set_clim(val[0], val[1])

        # uodate colorbar
        ImgTools._colorbar.update_normal(ImgTools._img)

        # Redraw the figure to ensure it updates
        ImgTools._figure.canvas.draw_idle()


    @staticmethod
    def open_image() -> None:
        # create openfile dialog

        path = askopenfilenames(title='Select image(s) or a directory for watch mode',
                            initialdir = ImgTools.conf.get_str('files', 'initial_directory'),
                            defaultextension = ImgTools.conf.get_str('files', 'file_types')
                            )
        #path = (['albireo-10.fit'])
        if path == '': 
            return 
        else: 
            ImgTools.load_image(path)
            return
        

 
    @staticmethod
    def load_image(path: tuple[str, ...]) -> None:

        # cleanup previous images
        ImgTools._ax_img.clear()

        # do not recreate colorbar
        if (ImgTools._colorbar is None):
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
#            _img_data.append(CCDData.read(img_name).data)
            _img_name.append(img_name)
            _img_count += 1

        # stack images
        ImgTools.img_stacked: np.ndarray = _img_data[0]
        for i in range(1, _img_count):
            ImgTools.img_stacked = np.add(ImgTools.img_stacked, _img_data[i])

        del _img_data

       # collect image stats
        vstd = ImgTools.img_stacked.std()
        vmean = ImgTools.img_stacked.mean()
        _min = ImgTools.img_stacked.min()
        _max = ImgTools.img_stacked.max()
        vmin = vmean - vstd
        vmax = vmean + vstd

        # display image
        logging.info (f"image stats : min = {_min}, max = {_max}, mean = {vmean}, std = {vstd}")
        ImgTools._img = ImgTools.show_image(image = ImgTools.img_stacked,
                        #percl = 0,
                        #percu = 99.5,
                        fig = ImgTools._ax_img.get_figure(),
                        ax = ImgTools._ax_img,
                        show_colorbar = show_colorbar, 
                        cmap = ImgTools.conf.get_str('window', 'colormap'))

        # update colorbar
        ImgTools._colorbar.update_normal(ImgTools._img)

        # update slider
        ImgTools._img.norm.vmax = vmax
        ImgTools._img.norm.vmin = vmin

        ImgTools._slider.valmin = vmin / 2 #_min 
        ImgTools._slider.valmax = vmax * 2
        ImgTools._slider.valinit = (vmin, vmax)
        ImgTools._slider_ax.set_xlim(vmin, vmax * 2)
        ImgTools._slider.reset()

        ImgTools._slider.on_changed(ImgTools.update_image)
        ImgTools._figure.canvas.draw_idle()

        logging.info(f"total memory used = {int(psutil.Process().memory_info().rss / (1024 * 1024))}MB")

    @staticmethod
    def show_image( image,
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
            ImgTools._colorbar: Colorbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
            
            # In case someone in the future wants to improve this:
            # https://joseph-long.com/writing/colorbars/
            # https://stackoverflow.com/a/33505522/3486425
            # https://matplotlib.org/mpl_toolkits/axes_grid/users/overview.html#colorbar-whose-height-or-width-in-sync-with-the-master-axes
            
        if not show_ticks:
            ax.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False)

        return im
