import logging
import numpy as np
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
    def __init__(self, ax: Axes):
        self.conf = Config().config
        self._ax = ax

        def format_coord(x,y):
            return f'x={x:.0f}, y={y:.0f}'
            
        self._ax.format_coord=format_coord
        
        #self._img = np.random.rand(4000, 5000)

        self._img = CCDData.read('albireo-10.fit', unit = u.adu).data
        self._img_ax = self.show_image(image = self._img,
                        #percl = 0,
                        #percu = 99.5,
                        fig = self._ax.get_figure(),
                        ax = self._ax,
                        #show_colorbar=False,
                        cmap = self.conf['window']['colormap'])

        self._img_ax.norm.vmin = self._img.min()
        self._img_ax.norm.vmax = self._img.max()

        # Create the RangeSlider
        #plt.subplots_adjust(top=0.25)
        self._slider_ax = self._img_ax.get_figure().add_axes([0.065, 0.95, 0.77, 0.015])
        # (left, bottom, width, height)
        self._slider = RangeSlider(self._slider_ax, "Cuts: ",
                                   orientation = 'horizontal',
                                   valinit = [self._img.min(), self._img.max()],
                                   valmin = 0,
                                   valmax = self._img.max() * 1.0)

        self._slider.on_changed(self.update)

    def update(self, val):
        # Update the image's colormap
        self._img_ax.norm.vmin = val[0]
        self._img_ax.norm.vmax = val[1]

        self._img_ax.set_clim([val[0], val[1]])

        # Redraw the figure to ensure it updates
        self._img_ax.get_figure().canvas.draw_idle()
       
    def show_image(self, image,
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
            fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
            
            # In case someone in the future wants to improve this:
            # https://joseph-long.com/writing/colorbars/
            # https://stackoverflow.com/a/33505522/3486425
            # https://matplotlib.org/mpl_toolkits/axes_grid/users/overview.html#colorbar-whose-height-or-width-in-sync-with-the-master-axes

        if not show_ticks:
            ax.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False)

        return im
