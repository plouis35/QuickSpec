import io
import logging

from colorlog import ColoredFormatter

from app import logger as logger_module


def test_initialize_replaces_default_console_handler(tmp_path, monkeypatch):
    root_logger = logging.getLogger()
    monkeypatch.setattr(root_logger, 'handlers', [])
    plain_output = io.StringIO()
    formatted_output = io.StringIO()
    logging.basicConfig(stream=plain_output, level=logging.INFO)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        logger_module,
        'Config',
        lambda: type('TestConfig', (), {'get_str': lambda *_: 'INFO'})(),
    )

    handler = logger_module.LogHandler()
    handler.stream = logging.StreamHandler(formatted_output)
    handler.stream.setLevel(logging.INFO)
    handler.stream.setFormatter(
        ColoredFormatter(
            "%(log_color)s%(asctime)s %(levelname)s - %(message)s%(reset)s ",
            "%H:%M:%S",
        )
    )
    handler.file_handler.close()
    handler.file_handler = logging.FileHandler(tmp_path / 'test.log')
    handler.file_handler.setLevel(logging.INFO)
    handler.file_handler.setFormatter(logging.Formatter('%(message)s'))
    try:
        handler.initialize()
        logging.info("single formatted message")

        assert plain_output.getvalue() == ''
        assert "INFO - single formatted message" in formatted_output.getvalue()
        assert sum(
            type(existing) is logging.StreamHandler
            for existing in root_logger.handlers
        ) == 1
    finally:
        for existing in tuple(root_logger.handlers):
            root_logger.removeHandler(existing)
            existing.close()
