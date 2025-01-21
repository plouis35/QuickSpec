import logging
import os, sys, platform
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
    
    #@staticmethod
    #def log_memory_used() -> None:
        #logging.info(f"current memory used = {OSUtils.get_memory_used()}MB")

    @staticmethod
    def show_versions() -> None:
        """
        show some packages version
        """        
        from tkinter import TkVersion

        logging.info("versions running: ")
        logging.info(f"platform = {platform.system()}, release = {platform.release()}")
        logging.info(f"python = {sys.version}")
        logging.info(f"tkinter = {TkVersion}")
        for module in ['matplotlib', 'numpy', 'astropy', 'specutils', 'specreduce', 'ccdproc']:
            try:
                logging.info(f"{module} = {version(module)}")
            except Exception as e:
                logging.info(f"{e}")

    @staticmethod
    def get_path_directory(path: str) -> str:
        """
        get absolute current selected directory

        Args:
            path (str): relative path

        Returns:
            str: _description_
        """        
        OSUtils._current_path = str(Path(path).absolute().parent)
        return OSUtils._current_path

    @staticmethod
    def get_current_path() -> str:
        """
        get current path

        Returns:
            str: current path
        """        
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
    

