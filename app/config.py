import os
import time
import configparser
import logging
import functools
from typing import LiteralString

class Config(object):
    """
    Singleton-type class to manage configuration ini file
    Handles exceptions (i.e. returns '' instead of raising except)
    Is designed to be initiated by several modules
    Returns:
        config: ConfigParser
    """

    # private class variables
    _instance = None
    _configDir: str = './'
    _configFile: str = 'quickspec.ini'
    _last_time_check: float = time.time()

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(Config, cls).__new__(cls)
            cls._instance._initialize()
        return cls._instance

    def _initialize(self) -> None:        
        self._config_path: str = os.path.join(Config._configDir, Config._configFile)
        self.config = configparser.ConfigParser(inline_comment_prefixes=('#', ';')) 

        if not os.path.isfile(self._config_path):
            self._new_config()

        self.config.read(self._config_path)

    def set_conf_directory(self, new_dir: str) -> None:
        Config._configDir = new_dir
        self._initialize()

    @staticmethod
    def _check_changes(func):

        @functools.wraps(func)
        def wrap(self, *args, **kwargs):
            # check if config file has been modified - if so re-read it
            if Config._last_time_check <= os.path.getmtime(self._config_path): 
                logging.warning("config has changed - reloading it...")
                self.config.read(self._config_path)
                Config._last_time_check = time.time()

            return func(self, *args, **kwargs)
        return wrap

    @_check_changes
    def get_str(self, section, key) -> str | None:
        try:
            return self.config.get(section, key)
        except (configparser.NoSectionError, configparser.NoOptionError, ValueError) as e:
            logging.warning(f"configuration key not found: {section}.{key} {e} in {self._config_path} file")
            return None

    @_check_changes
    def get_int(self, section, key) -> int | None:
        try:
            return self.config.getint(section, key)
        except (configparser.NoSectionError, configparser.NoOptionError, ValueError) as e:
            logging.warning(f"configuration key not found: {section}.{key} {e} in {self._config_path} file")
            return None

    @_check_changes
    def get_float(self, section, key) -> float | None:
        try:
            return self.config.getfloat(section, key)
        except (configparser.NoSectionError, configparser.NoOptionError, ValueError) as e:
            logging.warning(f"configuration key not found: {section}.{key} {e} in {self._config_path} file")
            return None

    @_check_changes
    def get_bool(self, section, key) -> bool | None:
        try:
            return self.config.getboolean(section, key)
        except (configparser.NoSectionError, configparser.NoOptionError, ValueError) as e:
            logging.warning(f"configuration key not found: {section}.{key} {e} in {self._config_path} file")
            return None

    def save(self) -> None:
        with open(self._config_path, 'w') as configfile:
            self.config.write(configfile)

    def check_section(self, section: str) -> None:
        if not self.config.has_section(section):
            self.config.add_section(section)

    def _new_config(self) -> None:
        # static config template to write to first opening of a directory
        contents: LiteralString = f"""[logger]
level = INFO

[display]
colormap = inferno                  # magma
contrast_level = 6                  # contrast 'magic' level (high=1 ... low=10)

[pre_processing]
#crop_region = 0, 1000, 5496, 2500  # x1, y1, x2, y2 
master_offset = masterbias.fit      #_offset.fit
master_dark = masterdark.fit        #_dark.fit
master_flat = masterflat.fit        #_flat.fit
master_science = masterscience.fit
cosmics_cleanup = no

[processing]
#trace_y_guess = 1695
trace_y_size = 15
trace_y_window = 50
trace_x_bins = 64
sky_y_size = 140
sky_y_offset = 120
# HR (starex2400)
calib_x_pixel = 770, 1190, 2240, 3520, 4160
calib_x_wavelength = 6506.53, 6532.88, 6598.95, 6678.28, 6717.04
calib_poly_order = 2
# LR (dados200 or alpy600)
#calib_x_pixel = 770, 1190, 2240, 3520, 4160
#calib_x_wavelength = 6506.53, 6532.88, 6598.95, 6678.28, 6717.04
#calib_poly_order = 2
reponse_file = masterresponse.fits

[post_processing]
#median_smooth = 7

[lines]
0.00 = Zero
656.2852 = Hα 
486.133 = Hβ 
434.047 = Hγ 
410.174 = Hδ 
397.007 = Hε 
388.9049 = Hζ 
383.5384 = Hη 
379.75 = Hθ 
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
585.249 = NeI
588.189 = NeI
589.00 = NaI
589.59 = NaI
615.82 = O1 
627.77 = O2 
686.9 = O2 
718.6 = O2 
760.5 = O2 
898.77 = O2 
495.9 = OIII,
500.69 = OIII
651.65 = H2O
694.07 = H2O
695.64 = H2O
698.90 = H2O
396.85 = Ca+ H
393.37 = Ca+ K
706.52 = He I
667.82 = He I
587.56 = He I
501.57 = He I
447.148 = He I
634.71 = Si II
637.14 = Si II
487.7 = Tb 
542.4 = Tb 
611.6 = Eu 
336.11 = Ti+

"""
        logging.warning(f"creating a new config file: {self._config_path}")
        with open(self._config_path, 'w') as configfile:
            configfile.write(contents + '\n')

# test
if __name__ == "__main__":
    conf = Config()

    for key, value in conf.config.items('lines'):
        print(f"\tfor key {key} -> {value} (value)")

    """"
    from astropy import units as u
    from astropy.table import QTable

    __wavelength = [6506.53, 6532.88, 6598.95, 6678.28, 6717.04]*u.AA
    __pixels = [770, 1190, 2240, 3484, 4160]*u.pix
    print(' ok -> ', __wavelength, __pixels)

    _wavelength = conf.get_str('processing', 'calib_x_wavelength') #[6506.53, 6532.88, 6598.95, 6678.28, 6717.04]*u.AA
    _pixels = conf.get_str('processing', 'calib_x_pixel')    #[770, 1190, 2240, 3484, 4160]*u.pix
    print('raw -> ', _wavelength, _pixels)

    #___wavelength = _wavelength
    #___wavelength = _wavelength.replace(',', '').split()
    ____wavelength = [float(x) for x in _wavelength.replace(',', '').split()]*u.AA
    print('cnv -> ', ____wavelength)

    #line_list = QTable([_pixels, _wavelength], names=["pixel_center", "wavelength"]) #, dtype=[u.pix, u.AA])

    wavelength = _wavelength
    pixels = _pixels
    print('nok -> ', wavelength, pixels)
    """
