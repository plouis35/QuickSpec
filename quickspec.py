"""_summary_
"""
__app__ = 'QuickSpec'
__version__ = '0.1'
import logging
from app.logger import LogHandler
from app.main import Application

if __name__ == "__main__":
    LogHandler().set()
    logging.info(f"{__app__} v{__version__} starting...")
    Application().run()
    logging.info(f"{__app__} ending...")
