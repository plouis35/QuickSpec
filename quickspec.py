import os, sys
import logging
import tkinter as tk
from tkinter import ttk
from app.main import Application

__version__ = '0.1'

if __name__ == "__main__":
    app = Application(__version__)
    logging.info("QuickSpec started")
    app.mainloop()
