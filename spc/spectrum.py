import logging
import numpy as np

import tkinter as tk
from tkinter import ttk

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.text import Annotation
from matplotlib.lines import Line2D

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends._backend_tk import NavigationToolbar2Tk
from matplotlib.colors import LinearSegmentedColormap

from astropy.utils.exceptions import AstropyWarning
from astropy import units as u

from specutils.spectra.spectrum1d import Spectrum1D
from specutils.analysis import snr, snr_derived

from specreduce.tracing import FlatTrace, FitTrace

from app.config import Config
from spc.spc_utils import SPCUtils

class Spectrum(object):

    def __init__(self, spc_frame: ttk.Frame, img_axe: Axes) -> None: #, axe_img: Axes, axe_spc: Axes) -> None:
        self.conf = Config()

        # declare/assign instance variables
        self.science_spectrum: Spectrum1D = None
        self.science_trace: FitTrace = None
        self.showed_lines: bool = False
        self.showed_colorized: bool = False
        self.colors = ('blue', 'red', 'green', 'orange', 'cyan')
        self.lines_color = 'yellow'    
        self.spectrum_color = 'grey'

        # create figure and axe
        self.spc_figure = Figure(figsize=(10, 3))
        self.spc_axe: Axes = self.spc_figure.add_subplot(111)
        self.img_axe: Axes = img_axe
        self.spc_axe.grid(color = 'grey', linestyle = '--', linewidth = 0.5)

        # create spectrum canvas to draw to
        self.spc_canvas = FigureCanvasTkAgg(self.spc_figure, spc_frame)
        self.spc_canvas.draw()

        # create customized toolbar
        self.spc_toolbar = CustomSpcToolbar(self.spc_canvas, spc_frame)

        # add new buttons
        self.clear_button = ttk.Button(self.spc_toolbar, text="Clear", command=self.reset_spectra)
        self.clear_button.pack(side=tk.LEFT, padx=5, pady=0)

        self.lines_button = ttk.Button(self.spc_toolbar, text="Show lines", command=self.show_lines)
        self.lines_button.pack(side=tk.LEFT, padx=5, pady=0)

        self.lines_button = ttk.Button(self.spc_toolbar, text="Colorize", command=self.colorize)
        self.lines_button.pack(side=tk.LEFT, padx=5, pady=0)

        self.spc_toolbar.update()
        self.spc_toolbar.pack(side=tk.TOP, fill=tk.X)
        self.spc_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def open_spectrum( self, spc_name: str) -> bool:
        # open spectrum data
        try:
            spec1d: Spectrum1D = Spectrum1D.read(spc_name)
            self.show_spectrum(spec1d, True)

        except Exception as e:
            logging.error(f"{e}")
            return False
        
        self.science_spectrum = spec1d
        return True

    def colorize(self) -> None:
        if self.science_spectrum is None:
            logging.info("please calibrate first")
            return
        
        min_wl = 3800 
        max_wl = 7500
        bin_size = 21 # size (AA) of each colored slice under spectrum
        alpha = 1.0

        if self.showed_colorized:
            # show 'black' color under spectrum
            for i in range(0, len(self.science_spectrum.wavelength) - 1, bin_size):
                self.spc_axe.fill_between(self.science_spectrum.wavelength.value[i:i+bin_size], 0, 
                                self.science_spectrum.flux.value[i:i+bin_size], color='black', alpha=alpha)
            self.showed_colorized = False
        else:
            # show 'rainbow' color under spectrum
            colors = plt.cm.turbo((self.science_spectrum.wavelength.value - min_wl) / (max_wl - min_wl))            
            for i in range(0, len(self.science_spectrum.wavelength) - 0, bin_size):
                self.spc_axe.fill_between(x=self.science_spectrum.wavelength.value[i:i+bin_size], 
                                          y1=0, y2=self.science_spectrum.flux.value[i:i+bin_size], 
                                          color=colors[i], alpha=alpha)
            self.showed_colorized = True
            
        self.spc_figure.canvas.draw_idle()

    def reset_spectra(self) -> None:
        self.science_spectrum = None
        self.science_trace = None
        self.clear_spectra()

    def clear_spectra(self) -> None:
        self.spc_axe.clear()
        self.spc_axe.grid(color = 'grey', linestyle = '--', linewidth = 0.5)
        self.showed_lines = False
        self.showed_colorized = False
        self.spc_toolbar.update()
        self.spc_figure.canvas.draw_idle()

    def show_spectrum(self, spectrum: Spectrum1D, calibrated: bool = False) -> None:
        # set legend titles
        if calibrated:
            self.spc_axe.set_ylabel('Relative intensity')
            self.spc_axe.set_xlabel('Wavelength (Angstrom)')
            def format_coord(x,y):
                return "Lambda: ({:.2f}, Intensity: {:.2f})".format(x,y)
            self.spc_axe.format_coord=format_coord
        else:
            self.spc_axe.set_xlabel('Pixels')
            self.spc_axe.set_ylabel('ADU')
            def format_coord(x,y):
                return "Pixel: ({:.0f}, ADU: {:.0f})".format(x,y)
            self.spc_axe.format_coord=format_coord
        
        # plot spectrum
        if spectrum == self.science_spectrum:
            _color = self.spectrum_color
        else:
            self.colors = np.roll(self.colors, 1) # rotate a new color every call 
            _color = self.colors[0]

        self.spc_axe.plot(spectrum.spectral_axis , spectrum.flux, color=_color, linewidth = '0.8')
        self.spc_axe.grid(color = 'grey', linestyle = '--', linewidth = 0.5)

        self.spc_figure.canvas.draw_idle()

    def do_trace(self, img_stacked:np.ndarray) -> bool:
        if img_stacked is None: 
            logging.error("please reduce image(s) before fitting trace")
            return False

        science_trace = SPCUtils.trace_spectrum(img_stacked)
        if science_trace is None:
            return False
        
        self.science_trace = science_trace

        # trace fitted spectrum
        for elm in self.img_axe.get_children():
                if isinstance(elm, Line2D): 
                    elm.remove()

        self.img_axe.plot(self.science_trace.trace , color='red', linestyle='dashed', linewidth = '0.5')
        self.img_axe.get_figure().canvas.draw_idle()

        return True


    def do_extract(self, img_stacked:np.ndarray) -> bool:

        if (img_stacked is None) or (self.science_trace is None): 
            logging.error("please fit trace before extracting spectrum")
            return False

        extracted_spectrum = SPCUtils.extract_spectrum(img_stacked, self.science_trace)
        if extracted_spectrum is None:
            return False
        
        self.science_spectrum = extracted_spectrum
        logging.info('spectrum extracted')

        #self.clear_spectra()

        # trace sky zones
        self.img_axe.plot(self.science_spectrum.spectral_axis, self.science_trace.trace + self.conf.get_int('processing', 'trace_y_size'), 
                color='blue', linestyle='dashed', linewidth = '0.5')  #, alpha=0.2)
        self.img_axe.plot(self.science_spectrum.spectral_axis, self.science_trace.trace - self.conf.get_int('processing', 'trace_y_size'), 
                color='blue', linestyle='dashed', linewidth = '0.5')  #, alpha=0.2)

        self.img_axe.plot(self.science_spectrum.spectral_axis, self.science_trace.trace + (self.conf.get_int('processing', 'sky_y_offset')) , 
                color='green', linewidth = '0.5', linestyle='dashed')  #, alpha=0.2)
        self.img_axe.plot(self.science_spectrum.spectral_axis, self.science_trace.trace - (self.conf.get_int('processing', 'sky_y_offset')) , 
                color='green', linewidth = '0.5', linestyle='dashed')  #, alpha=0.2)

        self.img_axe.plot(self.science_spectrum.spectral_axis, 
                self.science_trace.trace + (self.conf.get_int('processing', 'sky_y_offset') + self.conf.get_int('processing', 'sky_y_size')) , 
                color='green', linewidth = '0.5', linestyle='dashed')  #, alpha=0.2)
        self.img_axe.plot(self.science_spectrum.spectral_axis, 
                self.science_trace.trace - (self.conf.get_int('processing', 'sky_y_offset') + self.conf.get_int('processing', 'sky_y_size')) , 
                color='green', linewidth = '0.5', linestyle='dashed')  #, alpha=0.2)
        
        self.img_axe.get_figure().canvas.draw_idle()

        # show extracted spectrum
        self.clear_spectra()
        self.show_spectrum(self.science_spectrum, False)

        return True
    
    def do_calibrate(self) -> bool:
        if self.science_spectrum is None: 
            logging.error(f'please extract spectrum trace before calibrating')
            return False
        
        logging.info('calibrating spectrum...')
        
        _wavelength = self.conf.get_str('processing', 'calib_x_wavelength')  
        _pixels = self.conf.get_str('processing', 'calib_x_pixel') 

        if (_wavelength is None) or (_pixels is None):
            logging.error("please specify lines/waves index in config file")
            return False
        
        calibrated_spectrum = SPCUtils.calibrate_spectrum(self.science_spectrum, _pixels, _wavelength)
        if calibrated_spectrum is None:
            return False        

        self.science_spectrum = calibrated_spectrum

        # show derived-SNR
        logging.info(f'snr = {snr_derived(self.science_spectrum):.1f}')

        # show spectrum
        self.clear_spectra()
        self.show_spectrum(self.science_spectrum, True)

        logging.info('calibration complete')
        return True


    def do_response(self) -> bool:
        if self.science_spectrum is None: 
            logging.error(f'please calibrate spectrum before applying response file')
            return False
        
        final_spec: Spectrum1D = SPCUtils.apply_response(self.science_spectrum) 
        if final_spec is None:
            return False        

        self.science_spectrum = final_spec

        # show spectrum
        self.clear_spectra()
        self.show_spectrum(self.science_spectrum, True)

        logging.info('response applied')
        return True


    def do_smooth(self) -> bool:
        if self.science_spectrum is None: 
            logging.error(f'please calibrate spectrum trace before smoothing')
            return False
        
        # apply median smoothing (if any defined)
        smooth: int | None = self.conf.get_int('post_processing', 'median_smooth') 
        if smooth is not None:
            smooth_spec: Spectrum1D = SPCUtils.median_smooth(self.science_spectrum, smooth) 
            if smooth_spec is not None:
                self.science_spectrum = smooth_spec
                logging.info(f"median smooth applied={smooth}")
        else:
            self.science_spectrum = self.science_spectrum
            logging.info(f"no median smooth to apply")

        # show spectrum
        self.clear_spectra()
        self.show_spectrum(self.science_spectrum, True)

        return True

    def show_lines(self, ax = None, show_line = True) -> None:
        if ax is None: ax = self.spc_axe
                        
        xbounds = ax.get_xbound()   
        trans = ax.get_xaxis_transform()

        if self.showed_lines is False:
            for wave, elm in self.conf.config.items('lines'):
                _lambda = (float(wave) * 10)   # convert nm to ang
                if (_lambda > xbounds[0]) & (_lambda < xbounds[1]):
                        #ax.axvline(lam, 0.95, 1.0, color = self.lines_color, lw = 0.5)
                        ax.axvline(_lambda, color = self.lines_color, lw = 0.5, linestyle = '--', alpha=0.8)
                        #trans = ax.get_xaxis_transform()
                        ax.annotate(elm, xy=(_lambda, 1.05), xycoords=trans, fontsize=8, rotation=90, color=self.lines_color)
                        
            # show colorband
            self.showed_lines = True
        else:
            # clear lines
            for line in ax.lines: 
                if line.get_color() == self.lines_color:
                    line.remove() 

            # clear elements
            for elm in ax.get_children():
                if isinstance(elm, Annotation): 
                    elm.remove()
            
            # redraw spectrum
            #self.show_spectrum(self.sci_spectrum, True)
            self.showed_lines = False
                    
        self.spc_figure.canvas.draw_idle()

class CustomSpcToolbar(NavigationToolbar2Tk):
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

                
