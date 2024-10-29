__version__ = '0.1'
import logging
from app.main import Application

if __name__ == "__main__":
    app = Application(__version__)
    logging.info(f"QuickSpec version {__version__} ready")
    app.mainloop()
