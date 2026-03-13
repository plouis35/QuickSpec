"""
os_utils — OS and file system utility functions.
"""
import logging
import os
import sys
import platform
from pathlib import Path
from importlib.metadata import version

# Current working directory — updated when the user opens a file
_current_path: str = '.'


def get_path_directory(path: str) -> str:
    """
    Derive and store the parent directory of a given file path.

    Args:
        path (str): any file path

    Returns:
        str: absolute parent directory
    """
    global _current_path
    _current_path = str(Path(path).absolute().parent)
    return _current_path


def get_current_path() -> str:
    """
    Return the currently monitored directory.

    Returns:
        str: absolute path, or '.' if none selected yet
    """
    return _current_path


def show_versions() -> None:
    """Log platform and key package versions at startup."""
    from tkinter import TkVersion

    logging.info("versions installed:")
    logging.info(f"platform = {platform.system()}, release = {platform.release()}")
    logging.info(f"python = {sys.version}")
    logging.info(f"tkinter = {TkVersion}")
    for module in ['matplotlib', 'numpy', 'scipy', 'astropy', 'specutils', 'specreduce', 'ccdproc']:
        try:
            logging.info(f"{module} = {version(module)}")
        except Exception as e:
            logging.info(f"{e}")
