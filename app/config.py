"""
Singleton-type class to manage configuration INI files
Hide exceptions related to keyword/section missing (i.e. returns None instead of raising exception)
"""
import os
import configparser
import logging
import functools
from typing import LiteralString

class Config(object):
    # private 
    _instance = None
    _configDir: str = './'
    _configFile: str = 'quickspec.ini'

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(Config, cls).__new__(cls)
            cls._instance._initialize()
        return cls._instance

    def _initialize(self) -> None:    
        """
        configure INI file syntax
        create a new INI file if not exist
        """         
        self._config_path: str = os.path.join(Config._configDir, Config._configFile)
        if not os.path.isfile(self._config_path):
            self._new_config()
        self.config = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))
        self._config_mtime_ns = None
        self.reload()

    def set_conf_directory(self, new_dir: str) -> None:
        Config._configDir = new_dir
        self._initialize()

    def get_conf_directory(self) -> str:
        return Config._configDir

    @property
    def config_path(self) -> str:
        return self._config_path

    def reload(self) -> bool:
        """Load the current configuration file and remember its modification time."""
        try:
            config_mtime_ns = os.stat(self._config_path).st_mtime_ns
            config = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))
            with open(self._config_path, encoding='utf-8') as configfile:
                config.read_file(configfile)
        except (OSError, configparser.Error) as e:
            logging.error(f"unable to reload configuration {self._config_path}: {e}")
            return False

        self.config = config
        self._config_mtime_ns = config_mtime_ns
        logging.debug(f"configuration reloaded: {self._config_path}")
        return True

    def reload_if_changed(self) -> bool:
        """Reload the configuration when its file timestamp has changed."""
        try:
            config_mtime_ns = os.stat(self._config_path).st_mtime_ns
        except OSError as e:
            logging.error(f"unable to check configuration {self._config_path}: {e}")
            return False
        if config_mtime_ns == self._config_mtime_ns:
            return False
        return self.reload()

    @staticmethod
    def check_changes(func):
        """
        check if config file has been modified
        decorate all config read routines below

        Args:
            func (_type_): routine to decorate

        Returns:
            _type_: routine return code
        """        
        @functools.wraps(func)
        def wrap(self, *args, **kwargs):
            self.reload_if_changed()
            return func(self, *args, **kwargs)
        return wrap

    @check_changes
    def get_str(self, section, key) -> str | None:
        try:
            return self.config.get(section, key)
        except (configparser.NoSectionError, configparser.NoOptionError, ValueError) as e:
            logging.debug(f"configuration key not found: {section}.{key} {e} in {self._config_path} file")
            return None

    @check_changes
    def get_int(self, section, key) -> int | None:
        try:
            return self.config.getint(section, key)
        except (configparser.NoSectionError, configparser.NoOptionError, ValueError) as e:
            logging.debug(f"configuration key not found: {section}.{key} {e} in {self._config_path} file")
            return None

    @check_changes
    def get_float(self, section, key) -> float | None:
        try:
            return self.config.getfloat(section, key)
        except (configparser.NoSectionError, configparser.NoOptionError, ValueError) as e:
            logging.debug(f"configuration key not found: {section}.{key} {e} in {self._config_path} file")
            return None

    @check_changes
    def get_bool(self, section, key) -> bool | None:
        try:
            return self.config.getboolean(section, key)
        except (configparser.NoSectionError, configparser.NoOptionError, ValueError) as e:
            logging.debug(f"configuration key not found: {section}.{key} {e} in {self._config_path} file")
            return None

    def save(self) -> None:
        with open(self._config_path, 'w', encoding="utf-8") as configfile:
            self.config.write(configfile)

    def check_section(self, section: str) -> None:
        if not self.config.has_section(section):
            self.config.add_section(section)

    def _new_config(self) -> None:
        """
        generates a new config file template under selected directory
        NOTE: be sure README.md is updated accordingly
        """        
        # static config template to write to first opening of a directory
        contents: LiteralString = f"""
[logger]
level = INFO                                    # DEBUG, INFO, WARNING, ERROR

[display]
theme = dark                                    # dark, light

[sa100]
distance_mm = 35.0                               # SA100-to-sensor distance (mm)
pixel_size_um = 4.65                             # sensor pixel size (µm)
order_zero_offset_px = 0                         # fine pixel offset if needed (px)
half_height = 5                                  # spectrum half binning height (px)
sky_half_height = 10                             # sky half height (px)

[pre_processing]
#y_crop = 0.5, 0.3                              # y-relative center, y-relative size arround center
master_offset = _offset.fits                     # masterbias - generated by another app
master_dark = _dark.fits                         # masterdark - generated by another app
master_flat = _flat.fits                         # masterflat - generated by another app
max_memory = 1e9                                # max bytes allocated during processing
auto_process = No                              # detect and process new file under current selected directory
registration = Yes

[processing]
trace_method = fit                              # fit (automatic detection and trace), flat (horizontal trace)
peak_model = gaussian                           # max, gaussian, centroid
trace_model = models.Polynomial1D(degree=2)     # models.Chebyshev1D(), models.Legendre1D(), models.Spline1D()
#trace_y_guess = 1695                           # fixed y
trace_y_size = 15                               # spectrum bin size
trace_y_window = 50                             # fit mode window size to search 
trace_x_bins = 22                               # nb of split sections to fit trace
sky_substract = Yes                             # remove sky background ?
sky_y_size = 60                                 # y-size of sky bands
sky_y_offset = 60                               # offset to spectrum trace

calib_x_pixel = 770, 1190, 2240, 3520, 4160                             # pixels x-position of calibration lines
calib_x_wavelength = 6506.53, 6532.88, 6598.95, 6678.28, 6717.04        # corresponding wavelength in Angstrom

normalized_region = 6550, 6600                  # region to use for normalization
response_file = _rep.fits                       # response spectrum - generated by another app

[post_processing]
#median_smooth = 7                              # size of median kernel to apply 
#shift_wavelength = 3.5
#crop_region = 4200, 7150

[lines]                                         # ATTENTION: special character not allowed (such as greek chars...)
0.00 = Zero
#652.38 = atm
#653.67 = atm
#654.39 = atm
#654.86 = atm
#655.26 = atm
#657.21 = atm
#657.48 = atm
#658.08 = atm
#658.65 = atm
#659.44 = atm
656.28 = H
486.13 = H
434.04 = H
410.17 = H
397.00 = H 
388.90 = H 
383.53 = H
379.75 = H
#527.04 = Fe 
#516.89 = Fe 
#495.76 = Fe 
#466.81 = Fe 
#438.36 = Fe 
#430.79 = Fe 
#448.11 = MgII
518.36 = Mg 
517.27 = Mg 
516.73 = Mg 
#585.24 = NeI
#588.18 = NeI
589.00 = NaI
589.59 = NaI
#615.82 = O1 
#627.77 = O2 
#686.90 = O2 
#718.60 = O2 
#760.50 = O2 
#898.77 = O2 
#495.90 = OIII
#500.69 = OIII
#651.65 = H2O
#694.07 = H2O
#695.64 = H2O
#698.90 = H2O
396.85 = CaII
393.37 = CaII
#706.52 = He
#667.82 = He
#587.56 = He
#501.57 = He
#447.14 = He
#634.71 = SiII
#637.14 = SiII
#487.70 = Tb 
#542.40 = Tb 
#611.60 = Eu 
#336.11 = TiII
"""
        logging.debug(f"creating a new config file: {self._config_path}")
        with open(self._config_path, 'w', encoding="utf-8") as configfile:
            configfile.write(contents + '\n')
