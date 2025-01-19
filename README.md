# QuickSpec
Quick-look spectra reduction &amp; calibration tool - aimed at demonstrating spectroscopy to public
> :warning: does not produce science-quality spectra - please use instead:
> - [specinti](https://solex.astrosurf.com/specinti1_fr.html)
> - [demetra](https://www.shelyak.com/logiciel/logiciel-demetra/)
> - [spcaudace](http://spcaudace.free.fr)

# Usage
> [Youtube (French)](https://)

# Installation
from binaries:
> download and extract the following ZIP file:
> - [Windows](https://)
> - [MacOS](https://)
> - [Linux](https://)

> then go to your extracted location and execute:```quickspec[.exe]``` (can take several minutes for the first execution...)

from sources:
> - download and install a python interpreter: [miniforge](https://github.com/conda-forge/miniforge)
> - create a python environment: ```$ conda env create ...```
> - clone QuickSpec sources: ```$ git clone https://github.com/plouis35/QuickSpec.git``` (or download zip)
> - go to cloned directory: ```$ cd quickspec```
> - install requirements: ```$ pip install -r requirements.txt```
> - and run Quickspec : ```$ python quickspec.py```

# Configuration file reference: 

```
[logger]
level = INFO

[display]
theme = dark            # light
contrast_level = 6      # from 1 to 9 : 

[pre_processing]
crop_auto = 0.5, 0.3    # y central relative (0.0 = bottom, 1.0 = top), crop ratio (0.0 = 0%, 1.0 = 100%)
master_offset = _offset.fit
master_dark = _dark.fit
master_flat = _flat.fit

[processing]
#auto_process = False
trace_model = models.Chebyshev1D(degree=2)
#trace_model = models.Spline1D(degree=2)
#trace_model = models.Polynomial1D(degree=2)            # default
#trace_model = models.Legendre1D(degree=2)
#trace_y_guess = 1695
trace_y_size = 15
trace_y_window = 50
trace_x_bins = 12
sky_y_size = 140
sky_y_offset = 120

# HR (starex2400)
calib_x_pixel = 770, 1190, 2240, 3520, 4160
calib_x_wavelength = 6506.53, 6532.88, 6598.95, 6678.28, 6717.04
response_file = _rep.fits

# LR (alpy600)
#calib_x_pixel = 959, 1645, 2130
#calib_x_wavelength = 4200.7, 5852.4, 6965.4
#response_file = _rep.fits

# LR (dados200)
#calib_x_pixel: 61, 683, 940, 1400, 1540
#calib_x_wavelength: 4333.56, 5400.56, 5852.49, 6678.28, 6929.47
#response_file = _rep.fits

# MR (dados900)
#calib_x_pixel: 61, 683, 940, 1400, 1540
#calib_x_wavelength: 4333.56, 5400.56, 5852.49, 6678.28, 6929.47
#response_file = _rep.fits

[post_processing]
#median_smooth = 7

[lines]
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
