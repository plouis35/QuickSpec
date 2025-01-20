"""
Singleton-type class to manage configuration INI files
Hide exceptions related to keyword/section missing (i.e. returns None instead of raising exception)
"""
import os
import time
import configparser
import logging
import functools
from typing import LiteralString

class Config(object):
    # private 
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
        """
        configure INI file syntax
        create a new INI file if not exist
        """         
        self._config_path: str = os.path.join(Config._configDir, Config._configFile)
        self.config = configparser.ConfigParser(inline_comment_prefixes=('#', ';')) 

        if not os.path.isfile(self._config_path):
            self._new_config()

        self.config.read(self._config_path)

    def set_conf_directory(self, new_dir: str) -> None:
        Config._configDir = new_dir
        self._initialize()

    def get_conf_directory(self) -> str:
        return Config._configDir

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
            # check if config file has been modified - if so re-read it
            if Config._last_time_check <= os.path.getmtime(self._config_path): 
                logging.warning("config has been modified - reloading it...")
                self.config.read(self._config_path)
                Config._last_time_check = time.time()

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
        """        
        # static config template to write to first opening of a directory
        contents: LiteralString = f"""[logger]
level = INFO

[display]
theme = dark

[pre_processing]
master_offset = _offset.fit
master_dark = _dark.fit
master_flat = _flat.fit

[processing]
#trace_y_guess = 1695
trace_y_size = 15
trace_y_window = 50
trace_x_bins = 12
sky_y_size = 140
sky_y_offset = 120

calib_x_pixel = 770, 1190, 2240, 3520, 4160
calib_x_wavelength = 6506.53, 6532.88, 6598.95, 6678.28, 6717.04
response_file = _rep.fits

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
"""
        logging.info(f"creating a new config file: {self._config_path}")
        with open(self._config_path, 'w', encoding="utf-8") as configfile:
            configfile.write(contents + '\n')
