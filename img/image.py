import logging
import warnings
import psutil
import numpy as np

import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from tkinter.filedialog import askopenfilename

import matplotlib.pyplot as plt
from matplotlib.figure import Axes

from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.widgets import RangeSlider

from app.config import Config

from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty
from astropy.stats import mad_std
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning

from ccdproc import Combiner, combine, subtract_bias, subtract_dark, flat_correct
from ccdproc import trim_image, Combiner, ccd_process, cosmicray_median

#from align_combine import align_and_combine

from astropy import visualization as aviz
from astropy.nddata.blocks import block_reduce
from astropy.nddata.utils import Cutout2D

from astropy.modeling import models
from specreduce import tracing, background, extract

class spec2d:

    # set class variables
    _fig_axe = None         # main axe passed to class
    _img_name = []
    _img_axe = None         # showimage returned axe
    _img_data = []
    _img_count = 0
    _slider = None
    _slider_ax = None
    _img_stacked = None
    _colorbar = None
    _clear_image_ax = None
    _bt_clear_image = None

    
    
    def __init__(self, ax: Axes):
        warnings.simplefilter('ignore', category=AstropyWarning)
        warnings.simplefilter('ignore', UserWarning)

        spec2d.conf = Config().config
        
        # change default label format (x, y)
        def format_coord(x,y):
            return f'x={x:.0f}, y={y:.0f}'

        spec2d._fig_axe = ax
        spec2d._fig_axe.format_coord=format_coord
        
        # Create the RangeSlider
        spec2d._slider_ax = ax.get_figure().add_axes([0.052, 0.98, 0.75, 0.015])
        # (left, bottom, width, height)
        spec2d._slider = RangeSlider(self._slider_ax, "Cuts: ",
                                   orientation = 'horizontal',
                                   valstep = 10,
                                   dragging = True,
                                   valinit = [0, 65535], #[vmin, vmax],
                                   valmin = 0,
                                   valmax = 65535) #vmax * 1.5)

        # create a dummy (zeros) image to start            
        spec2d._img_axe = spec2d.show_image(image = np.zeros((2, 8)),
                        #percl = 0,
                        #percu = 99.5,
                        fig = spec2d._fig_axe.get_figure(),
                        ax = spec2d._fig_axe,
                        show_colorbar=True, 
                        cmap = 'Grays') #spec2d.conf['window']['colormap'])

        """
        # 
        spec2d._clear_image_ax = ax.get_figure().add_axes([0.045, 0.90, 0.05, 0.025])
        spec2d._bt_clear_image = Button(spec2d._clear_image_ax,
                                        'Clear',
                                        color='k',
                                        hovercolor='r')
        spec2d._bt_clear_image.on_clicked(spec2d.clear_image)
        """
        
    @staticmethod
    def clear_image(event):
        logging.info("clear frame")

    @staticmethod
    def update(val):
        # Update the image's colormap
        spec2d._img_axe.norm.vmin = val[0]
        spec2d._img_axe.norm.vmax = val[1]
        
        # update image cuts levels
        spec2d._img_axe.set_clim([val[0], val[1]])

        # uodate colorbar
        spec2d._colorbar.update_normal(spec2d._img_axe)

        # Redraw the figure to ensure it updates
        spec2d._fig_axe.get_figure().canvas.draw_idle()

    @staticmethod
    def load_image():
        # create openfile dialog
        path = askopenfilename(title='Select image(s)',
                               initialdir = spec2d.conf['files']['initial_directory'],
                               defaultextension = spec2d.conf['files']['file_types'],
                               multiple = True)
        if path == '': return

        # cleanup previous images
        spec2d._fig_axe.clear()

        # do not recreate colorbar
        if (spec2d._img_axe is None):
            show_colorbar = True
        else:
            show_colorbar = False
            
        spec2d._img_count = 0
        spec2d._img_data = []
        spec2d._img_name = []
        
        # read new images into memory
        for img_name in path:
            logging.info(f"loading {img_name}...")
            spec2d._img_data.append(CCDData.read(img_name, unit = u.adu).data)
            spec2d._img_name.append(img_name)
            spec2d._img_count += 1

        # stack images
        spec2d._img_stacked = spec2d._img_data[0]
        for i in range(1, spec2d._img_count):
            spec2d._img_stacked = np.add(spec2d._img_stacked, spec2d._img_data[i])

       # collect image stats
        vstd = spec2d._img_stacked.std()
        vmean = spec2d._img_stacked.mean()
        _min = spec2d._img_stacked.min()
        _max = spec2d._img_stacked.max()
        vmin = vmean - vstd
        vmax = vmean + vstd

        # display image
        logging.info (f"image stats : min = {_min}, max = {_max}, mean = {vmean}, std = {vstd}")
        spec2d._img_axe = spec2d.show_image(image = spec2d._img_stacked,
                        #percl = 0,
                        #percu = 99.5,
                        fig = spec2d._fig_axe.get_figure(),
                        ax = spec2d._fig_axe,
                        show_colorbar = show_colorbar, 
                        cmap = spec2d.conf['window']['colormap'])

        # update colorbar
        spec2d._colorbar.update_normal(spec2d._img_axe)

        # update slider
        spec2d._img_axe.norm.vmin = vmin
        spec2d._img_axe.norm.vmax = vmax

        spec2d._slider.valmin = _min
        spec2d._slider.valmax = vmax * 2
        spec2d._slider.valinit = [vmin, vmax]
        spec2d._slider_ax.set_xlim(vmin, vmax * 2)
        spec2d._slider.reset()

        spec2d._slider.on_changed(spec2d.update)
        spec2d._fig_axe.get_figure().canvas.draw_idle()

        logging.info(f"memory used = {int(psutil.Process().memory_info().rss / (1024 * 1024))}MB")

    @staticmethod
    def show_image(image,
                   percl=99.5, percu=None, is_mask=False,
                   figsize=(10, 10),
                   cmap='viridis', log=False, clip=True,
                   show_colorbar=True, show_ticks=False,
                   fig=None, ax=None, input_ratio=None):
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

        ratio = (image.shape // fig_size_pix).max()

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
        extent = [0, image.shape[1], 0, image.shape[0]]

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
            #if (spec2d._colorbar is not None): spec2d._colorbar.remove()            
            spec2d._colorbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
            
            # In case someone in the future wants to improve this:
            # https://joseph-long.com/writing/colorbars/
            # https://stackoverflow.com/a/33505522/3486425
            # https://matplotlib.org/mpl_toolkits/axes_grid/users/overview.html#colorbar-whose-height-or-width-in-sync-with-the-master-axes
            
        if not show_ticks:
            ax.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False)

        return im
