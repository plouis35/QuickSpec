"""Application entry point."""
__app__ = 'QuickSpec'
__version__ = '0.9'

from PySide6.QtWidgets import QApplication

from app.main import Application

if __name__ == "__main__":
    app = QApplication([])
    app.setApplicationName(__app__)
    window = Application(__app__, __version__)
    window.show()
    raise SystemExit(app.exec())
