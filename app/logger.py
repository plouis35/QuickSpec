import logging
from colorlog import ColoredFormatter
from app.config import Config

class LogHandler(logging.Handler):
    def __init__(self) -> None:
        self.conf = Config()
        
        LOG_LEVEL = self.conf.get_str('logger', 'level')
        LOGFORMAT = format = "%(log_color)s%(asctime)s %(levelname)s - %(message)s%(reset)s "

        logging.root.setLevel(LOG_LEVEL)
        formatter = ColoredFormatter(LOGFORMAT)
        self.stream = logging.StreamHandler()
        self.stream.setLevel(LOG_LEVEL)
        self.stream.setFormatter(formatter)
        self.log = logging.getLogger()
        self.log.setLevel(LOG_LEVEL)

    def set(self) -> None:
        self.log.addHandler(self.stream)
