from pathlib import Path

from watchdog.events import FileCreatedEvent, FileModifiedEvent

from app.file_monitor import _QuickSpecEventHandler


def test_monitor_dispatches_config_changes_and_new_fits(tmp_path):
    config_path = tmp_path / 'quickspec.ini'
    config_path.touch()
    new_fits = tmp_path / 'capture.FITS'
    config_changes = []
    fits_files = []
    handler = _QuickSpecEventHandler(
        config_path,
        on_new_fits=fits_files.append,
        on_config_changed=lambda: config_changes.append(True),
    )

    handler.on_modified(FileModifiedEvent(str(config_path)))
    handler.on_created(FileCreatedEvent(str(new_fits)))

    assert config_changes == [True]
    assert fits_files == [str(new_fits)]
