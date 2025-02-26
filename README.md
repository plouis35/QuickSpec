# QuickSpec
![Alt text](./screenshot_01.PNG)

Quick-look spectra reduction &amp; calibration tool - aimed at demonstrating spectroscopy science to public or students

:warning: does not produce science-quality spectra - please use instead:

- [specinti](https://solex.astrosurf.com/specinti1_fr.html)
- [demetra](https://www.shelyak.com/logiciel/logiciel-demetra/)
- [spcaudace](http://spcaudace.free.fr)
- [EasyAstro](https://github.com/plouis35/EasyAstro.git)

# Functions
- images viewer with cuts level, zoom, pan, registration and auto sum
- spectra colorized viewer with selected lines display
- one-click images reducer, registration, spectrum extract, calibration, response and normalisation processing
- new image detection for automatic processing during capture session
- single parameter file with changes dynamically applied

# Current limits
- only FITS format images can be read
- only FITS or DAT (2-columns) format spectra can be read
- generated spectra can only be saved as images

# Pre-requisites
- DOF (bias, dark, flat) and response reference files are to be generated by another app

# Installation
from binaries:

download and extract the following ZIP file
- [Windows](https://<to_come...>) (not generated yet)
- [MacOS](https://<to_come...>) (not generated yet)
- [Linux](https://<to_come...>) (not generated yet)

then go to your extracted location and execute:```quickspec[.exe]``` (can take several minutes for the first execution...)

from sources:

- download and install a python interpreter: [miniforge](https://github.com/conda-forge/miniforge)
- update (or create a new) a python environment: ```$ conda env update --file environment.yml```
- clone QuickSpec sources: ```$ git clone https://github.com/plouis35/QuickSpec.git``` (or download zip)
- go to cloned directory: ```$ cd quickspec```
- install requirements: ```$ pip install -r requirements.txt```
- and run Quickspec : ```$ python quickspec.py```

# Usage
- move the QuickSpec console window so that logging information is always visible
- press ```load``` button to navigate to and select spectra and/or calibration images; they will automatically be summed
- use a text editor to open and adjust the parameter file created under that images directory (quickspec.ini); it is recommended to open the calibration images first in order to define the proper x-pixel / wavelength tuples
- press ```run all``` button to generate the calibrated spectrum
- press ```colorize``` and ```show lines``` to explain spectrum properties
- for demonstration purposes (and debugging), do not hesite to use the ```run step``` button
- press ```clear``` or any panel to reset loaded data

# Configuration file reference (quickspec.ini)
```
[logger]
level = INFO                                    # DEBUG, INFO, WARNING, ERROR

[display]
theme = dark                                    # dark, light

[pre_processing]
#y_crop = 0.5, 0.3                              # y-relative center, y-relative size arround center
master_offset = _offset.fit                     # masterbias - generated by another app
master_dark = _dark.fit                         # masterdark - generated by another app
master_flat = _flat.fit                         # masterflat - generated by another app
max_memory = 1e9                                # max bytes allocatable during processing
auto_process = Yes                              # detect and process new file under current selected directory
registration = No                              # align all spectra images to the first image

[processing]
trace_method = fit                              # fit (fit a curve), flat (horizontal line)
peak_model = gaussian                           # max, gaussian, centroid
trace_model = models.Polynomial1D(degree=2)     # models.polynomial.Chebyshev1D(), models.polynomial.Legendre1D(), models.Spline1D()
#trace_y_guess = 1695                           # fixed y
trace_y_size = 15                               # spectrum bin size
trace_y_window = 50                             # fit mode window size to search 
trace_x_bins = 12                               # nb of split sections to fit trace
sky_substract = Yes                             # remove sky background ?
sky_y_size = 60                                 # y-size of sky bands
sky_y_offset = 60                               # offset to spectrum trace

calib_x_pixel = 770, 1190, 2240, 3520, 4160                             # pixels x-position of calibration lines
calib_x_wavelength = 6506.53, 6532.88, 6598.95, 6678.28, 6717.04        # corresponding wavelength in Angstrom
input_model = models.Polynomial1D(degree=2)                             # degree changes allowed only

normalized_region = 6500, 6520                  # region to use for normalization
response_file = _rep.fits                       # response spectrum - generated by another app

[post_processing]
#shift_wavelength = 5.5                         # wavelength shift in Angstrom
#median_smooth = 7                              # size of median kernel to apply 

[lines]                                         # WARNING: special character not allowed (such as greek chars...)
0.00 = Zero
656.28 = H
486.13 = H
434.04 = H
410.17 = H
397.00 = H 
388.90 = H 
383.53 = H
379.75 = H
527.04 = Fe 
516.89 = Fe 
495.76 = Fe 
466.81 = Fe 
438.36 = Fe 
430.79 = Fe 
448.11 = MgII
518.36 = Mg 
517.27 = Mg 
516.73 = Mg 
585.24 = NeI
588.18 = NeI
589.00 = NaI
589.59 = NaI
615.82 = O1 
627.77 = O2 
686.90 = O2 
718.60 = O2 
760.50 = O2 
898.77 = O2 
495.90 = OIII
500.69 = OIII
651.65 = H2O
694.07 = H2O
695.64 = H2O
698.90 = H2O
396.85 = CaII
393.37 = CaII
706.52 = He
667.82 = He
587.56 = He
501.57 = He
447.14 = He
634.71 = SiII
637.14 = SiII
487.70 = Tb 
542.40 = Tb 
611.60 = Eu 
336.11 = TiII
```

# Some spectrograph calibration configurations: 

```
# HR (starex2400)
calib_x_pixel = 770, 1190, 2240, 3520, 4160
calib_x_wavelength = 6506.53, 6532.88, 6598.95, 6678.28, 6717.04
Ò
# LR (dados200)
calib_x_pixel: 61, 683, 940, 1400, 1540
calib_x_wavelength: 4333.56, 5400.56, 5852.49, 6678.28, 6929.47

# MR (dados900)
calib_x_pixel = 349.965, 376.590, 441.225, 487.131, 539.054, 608.398,  728.015
calib_x_wavelength = 4500.9, 4524.7, 4582.7, 4624.3, 4671.2, 4734.1, 4843.3

```