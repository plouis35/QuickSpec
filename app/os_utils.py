import logging
import sys, os
import platform
import time
from pathlib import Path
#import psutil
from importlib.metadata import version  
import fnmatch, pathlib

class OSUtils(object):
    # private to class
    _current_path: str = '.'

    #@staticmethod
    #def get_memory_used() -> float:
        #return int(psutil.Process().memory_info().rss / (1024 * 1024))
    
    @staticmethod
    def log_memory_used() -> None:
        logging.info(f"current memory used = {OSUtils.get_memory_used()}MB")

    @staticmethod
    def show_versions() -> None:
        logging.info("Versions installed: ")
        for module in ['matplotlib', 'numpy','astropy', 'specutils', 'specreduce', 'ccdproc']:
            try:
                print(f"{module} = v{version(module)}")
                #logging.info(f"{module} = v{version(module)}")
            except Exception as e:
                logging.info(f"{e}")

    @staticmethod
    def get_path_directory(path: str) -> str:
        OSUtils._current_path = str(Path(path).absolute().parent)
        return OSUtils._current_path

    @staticmethod
    def get_current_path() -> str:
        return OSUtils._current_path

    @staticmethod
    def list_files(path: str, name: str = '*') -> list:
        """
        returns a list of files under a path, reverse sorted by last modified time

        Args:
            path (str): path to scan 
            name (str, optional): name filter. Defaults to '*'.

        Returns:
            list(str): list of files (sorted)
        """      
        image_types = ['fit', 'fts', 'fits']
        return fnmatch.filter((str(i).split(os.sep)[-1] for i in sorted(
                    pathlib.Path(path).iterdir(), 
                    key = os.path.getmtime, 
                    reverse = True)
                ) , name)


# test
if __name__ == "__main__":
    OSUtils.show_versions()
    #print(OSUtils.get_path_directory('/Users/papa/Documents/ASTRO/CAPTURES/202408id/_offset.fit'))
    #print(OSUtils.list_files('/Users/papa/Documents/ASTRO/CAPTURES/20231007_Void', '*.fit*'))

