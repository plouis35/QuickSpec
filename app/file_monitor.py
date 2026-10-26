"""Filesystem monitoring for configuration edits and newly captured FITS files."""
import logging
from pathlib import Path

from watchdog.events import FileCreatedEvent, FileSystemEvent, FileSystemEventHandler
from watchdog.observers import Observer


class _QuickSpecEventHandler(FileSystemEventHandler):
    FITS_SUFFIXES = {'.fit', '.fts', '.fits'}

    def __init__(self, config_path: Path, on_new_fits, on_config_changed) -> None:
        super().__init__()
        self.config_path = config_path.resolve()
        self.on_new_fits = on_new_fits
        self.on_config_changed = on_config_changed

    def on_created(self, event: FileCreatedEvent) -> None:
        if not event.is_directory:
            self._handle_new_file(Path(event.src_path))

    def on_modified(self, event: FileSystemEvent) -> None:
        if not event.is_directory:
            self._handle_config(Path(event.src_path))

    def on_moved(self, event: FileSystemEvent) -> None:
        if event.is_directory:
            return
        destination = getattr(event, 'dest_path', None)
        if destination is not None:
            self._handle_new_file(Path(destination))

    def _handle_config(self, path: Path) -> None:
        if path.resolve() == self.config_path:
            self.on_config_changed()

    def _handle_new_file(self, path: Path) -> None:
        if path.resolve() == self.config_path:
            self.on_config_changed()
        elif path.suffix.lower() in self.FITS_SUFFIXES:
            self.on_new_fits(str(path))


class FileMonitor:
    """Watch a directory for configuration edits and new FITS images."""

    def __init__(self, on_new_fits, on_config_changed) -> None:
        self.on_new_fits = on_new_fits
        self.on_config_changed = on_config_changed
        self._observer: Observer | None = None

    def start(self, directory: str, config_path: str) -> None:
        self.stop()
        path = Path(directory).resolve()
        if not path.is_dir():
            logging.error(f"unable to monitor missing directory: {path}")
            return

        handler = _QuickSpecEventHandler(
            Path(config_path), self.on_new_fits, self.on_config_changed
        )
        self._observer = Observer()
        self._observer.schedule(handler, path=str(path), recursive=False)
        self._observer.start()
        logging.info(f"monitoring files in: {path}")

    def stop(self) -> None:
        if self._observer is not None:
            self._observer.stop()
            self._observer.join()
            self._observer = None
            logging.info("file monitoring stopped")
