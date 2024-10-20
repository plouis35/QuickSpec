__version__ = '0.1'

import os, sys
import logging
import tkinter as tk
from tkinter import ttk
from app.main import Application

if __name__ == "__main__":
    app = Application(__version__)
    logging.info("QuickSpec started")
    app.mainloop()
