"""
app starter class
"""
__app__ = 'QuickSpec'
__version__ = '0.5'

from app.main import Application

if __name__ == "__main__":
    Application(__app__, __version__).mainloop()  
