import logging
import psutil

class OSUtils(object):

    @staticmethod
    def get_memory_used() -> float:
        return int(psutil.Process().memory_info().rss / (1024 * 1024))
    
    @staticmethod
    def log_memory_used() -> None:
        logging.info(f"current memory used = {OSUtils.get_memory_used()}MB")
