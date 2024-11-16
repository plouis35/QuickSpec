import logging
import sys, os
import platform
from pathlib import Path
import psutil
from importlib.metadata import version  

class OSUtils(object):

    @staticmethod
    def get_memory_used() -> float:
        return int(psutil.Process().memory_info().rss / (1024 * 1024))
    
    @staticmethod
    def log_memory_used() -> None:
        logging.info(f"current memory used = {OSUtils.get_memory_used()}MB")

    @staticmethod
    def show_versions() -> None:
        logging.info("Versions installed: ")
        for module in ['numpy','astropy', 'specutils', 'specreduce', 'ccdproc']:
            try:
                logging.info(f"{module} = v{version(module)}")
            except Exception as e:
                logging.info(f"{e}")

    @staticmethod
    def get_path_directory(path: str) -> str:
        return str(Path(path).absolute().parent)

    # test
if __name__ == "__main__":
    OSUtils.show_versions()
    print(OSUtils.get_path_directory('/Users/papa/Documents/ASTRO/CAPTURES/202408id/_offset.fit'))

