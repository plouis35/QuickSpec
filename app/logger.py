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
        CONSOLE_LOGFORMAT = "%(log_color)s%(asctime)s %(levelname)s - %(message)s%(reset)s "
        FILE_LOGFORMAT = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

        logging.root.setLevel(LOG_LEVEL)
        formatter = ColoredFormatter(CONSOLE_LOGFORMAT)
        self.stream = logging.StreamHandler()
        self.stream.setLevel(LOG_LEVEL)
        self.stream.setFormatter(formatter)
        self.log = logging.getLogger()
        self.log.setLevel(LOG_LEVEL)

        self.file_handler = logging.FileHandler(f"quickspec_{datetime.now().strftime('%d-%m-%Y_%H%M%S')}.log")
        self.file_handler.setLevel(LOG_LEVEL)
        self.file_handler.setFormatter(FILE_LOGFORMAT)

    def initialize(self) -> None:
        """
        add console and file handlers
        to be called once when app starts
        other modules juste import logging and use logging.warn/info/error to log messages
        """        
        self.log.addHandler(self.stream)
        self.log.addHandler(self.file_handler)
