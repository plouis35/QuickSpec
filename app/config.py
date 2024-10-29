import os
import configparser
import logging
from typing import LiteralString

class Config:
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

    def get_config(self, section, key) -> str:
        try:
            return self.config.get(section, key)
        except (configparser.NoSectionError, configparser.NoOptionError, ValueError):
            logging.error(f"configuration error : {section}.{key} key not found in {self._config_path} file")
            return ""

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

if __name__ == "__main__":
    config_manager = Config()
    log_level = config_manager.get_config('logger', 'leel')
    print(f'log level: {log_level}')
