"""
logging routines
output to console and a new log file created under app directory
"""
import logging
from colorlog import ColoredFormatter
from datetime import datetime
from app.config import Config

class LogHandler(logging.Handler):
    def __init__(self) -> None:
        self.conf = Config()
        
        LOG_LEVEL: str | None = self.conf.get_str('logger', 'level')

        logging.root.setLevel(LOG_LEVEL)
        console_formatter = ColoredFormatter("%(log_color)s%(asctime)s %(levelname)s - %(message)s%(reset)s ",
            "%H:%M:%S")
        file_formatter = logging.Formatter(
            "%(asctime)s : %(levelname)s : [%(filename)s:%(lineno)s - %(funcName)s()] : %(message)s",
            "%Y-%m-%d %H:%M:%S")
        self.stream = logging.StreamHandler()
        self.stream.setLevel(LOG_LEVEL)
        self.stream.setFormatter(console_formatter)
        self.log = logging.getLogger()
        self.log.setLevel(LOG_LEVEL)

        self.file_handler = logging.FileHandler(f"quickspec_{datetime.now().strftime('%d-%m-%Y_%H%M%S')}.log")
        self.file_handler.setLevel(LOG_LEVEL)
        self.file_handler.setFormatter(file_formatter)

    def initialize(self) -> None:
        """
        add console and file handlers
        to be called once when app starts
        other modules juste import logging and use logging.warn/info/error to log messages
        """        
        self.log.addHandler(self.stream)
        self.log.addHandler(self.file_handler)
