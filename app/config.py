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
    _configFile: str = 'quickspec.ini'
    _configDir: str = './'
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

        contents: LiteralString = f"""[logger]
level = INFO
nb_lines = 5                        # lines displayed
font = 'Console', 14                # font, size

[window]
theme = dark_background             # 'classic', 'dark_background', 'fast', 'fivethirtyeight', 'ggplot', 'grayscale'
geometry = 1200x800+330+133         # widthxheight+x+y
state = normal                      # normal, iconic, withdrawn, or zoomed
colormap = magma                    #
"""
        with open(self._config_path, 'w') as configfile:
            configfile.write(contents + '\n')

# test
if __name__ == "__main__":
    conf = Config()

    from astropy import units as u
    from astropy.table import QTable

    """"
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
