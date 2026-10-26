"""
Interactive 'removal' of deviant lines from a 1D spectrum
New values are interpolated from a mouse selected region boundary * random value to simulate noise
"""
import sys
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('qtagg')
import matplotlib.pyplot as plt
from PySide6.QtWidgets import QApplication, QFileDialog
from matplotlib.widgets import Button, TextBox, SpanSelector
from specutils import Spectrum1D
from astropy.io import fits
import astropy.units as u

app = QApplication.instance() or QApplication(sys.argv)

# globals
file_path: str = ''
spectrum: Spectrum1D = Spectrum1D(spectral_axis=[0] * u.Unit('Angstrom'), flux=[0] * u.Unit('adu') )
fits_header: dict()
noise_sigma: float = 0.005 # initial 'noise' level (sigma to random.normal)
history = [spectrum.flux.value.copy()]

# Display spectrum
plt.style.use('dark_background')  
fig, ax = plt.subplots()
fig.suptitle(Path(file_path).stem)
ax.set_title('')
ax.grid(color = 'grey', linestyle = '--', linewidth = 0.5)
ax.set_xlabel("Wavelength (Å)")
ax.set_ylabel("Relative flux")
line, = ax.plot(spectrum.spectral_axis.value, spectrum.flux.value) #, label="Spectrum")
#ax.set_xlim(spectrum.spectral_axis.value[0], spectrum.spectral_axis.value[-1])

# keep space for
plt.subplots_adjust(bottom=0.3)

# callback functions to buttons
def cb_undo(event) -> None:
    global spectrum
    if len(history) <= 1:
        return
    history.pop()
    spectrum = Spectrum1D(spectral_axis=spectrum.spectral_axis, flux=history[-1] * u.Unit("adu"))
    line.set_ydata(spectrum.flux.value)
    fig.canvas.draw()

def cb_save_fits(event) -> None:
    if file_path and fits_header:
        hdu = fits.PrimaryHDU(data=spectrum.flux.value, header=fits_header)
        hdu.header['CRVAL1'] = spectrum.spectral_axis.value[0]
        hdu.header['CDELT1'] = spectrum.spectral_axis.value[1] - spectrum.spectral_axis.value[0]
        hdu.header['NAXIS1'] = len(spectrum.flux.value)
        hdu.writeto(f"{Path(file_path).stem}_corrected{Path(file_path).suffix}", overwrite=True)
        sys.exit(0)
    else:
        print(f"FATAL: no file selected")
        return

def cb_load_fits(event) -> None:
    global spectrum
    global history
    global fits_header
    global file_path
    global line

    file_path, _selected_filter = QFileDialog.getOpenFileName(
        None,
        'Select spectrum',
        '',
        'FITS files (*.fit *.fts *.fits)',
    )

    if not file_path:
        print(f"FATAL: no file selected")
        return

    try:
        spectrum = Spectrum1D.read(file_path)
        line, = ax.plot(spectrum.spectral_axis.value, spectrum.flux.value) #, label="Spectrum")
        ax.set_xlim(spectrum.spectral_axis.value[0], spectrum.spectral_axis.value[-1])
        fig.canvas.draw()
    except Exception as e:
        print(f"FATAL: {file_path} : {e}")
        return

    # save original fits header
    fits_header=spectrum.meta['header']

    # save initial undo history 
    history = [spectrum.flux.value.copy()]

def cb_update_sigma(text: str) -> None:
    global noise_sigma
    try:
        noise_sigma = float(text)
        fig.canvas.draw()
    except ValueError:
        pass

# create buttons & sigma entry
ax_load = plt.axes([0.05, 0.05, 0.15, 0.075]) #left, bottom, width, height
btn_load = Button(ax_load, 'Load', color='blue')
btn_load.on_clicked(cb_load_fits)

ax_undo = plt.axes([0.25, 0.05, 0.15, 0.075]) #left, bottom, width, height
btn_undo = Button(ax_undo, 'Undo', color='blue')
btn_undo.on_clicked(cb_undo)

ax_save = plt.axes([0.45, 0.05, 0.15, 0.075]) #left, bottom, width, height
btn_save = Button(ax_save, 'Save & exit', color='blue')
btn_save.on_clicked(cb_save_fits)

ax_sigma = plt.axes([0.86, 0.05, 0.1, 0.075]) #left, bottom, width, height
text_sigma = TextBox(ax_sigma, 'Noise: ', initial=str(noise_sigma), color='blue')
text_sigma.on_submit(cb_update_sigma)


# correction func
def apply_interpolation(xmin: int, xmax: int, sigma: float) -> None:
    """
    generate new values interpolated from xmin/xmax selected zone
    add a random multiplier to every value to simulate noise

    Args:
        xmin (int): lower x axis selection
        xmax (int): upper y axis selection
        sigma (float): 'noise' simulated level
    """    
    # extract selected zone into a sub-spectrum
    wavelength_mask = (spectrum.spectral_axis.value >= xmin) & (spectrum.spectral_axis.value <= xmax)
    x_mask = spectrum.spectral_axis.value[wavelength_mask]

    # extract values to be interpolated from
    wavelength_left = np.where(spectrum.spectral_axis.value < xmin)[0][-1]
    wavelength_right = np.where(spectrum.spectral_axis.value > xmax)[0][0]
    flux_left = spectrum.flux.value[wavelength_left]
    flux_right = spectrum.flux.value[wavelength_right]

    # interpolate
    interp_values = np.interp(x_mask, [xmin, xmax], [flux_left, flux_right])

    # add noise to interpolated values
    noise = np.random.normal(0, sigma, size=len(interp_values))
    spectrum.flux.value[wavelength_mask] = interp_values + noise

    # update history for undo
    history.append(spectrum.flux.value.copy())

    # update graph
    line.set_ydata(spectrum.flux.value)
    fig.canvas.draw()

# define range to select 
def on_select(xmin: int, xmax: int):
    apply_interpolation(xmin, xmax, noise_sigma)

span = SpanSelector(ax, on_select, 'horizontal', useblit=True, props=dict(alpha=0.5, facecolor='red'))

# show everything and wait for buttons
plt.show()
