from pathlib import Path

from app.config import Config


def test_configuration_defines_supported_application_themes():
    config = Config()

    assert config.get_str('display', 'theme') == 'dark'


def test_reload_if_changed_loads_updated_values():
    config = Config()
    config_file = Path(config.config_path)
    config_file.write_text(
        config_file.read_text().replace('contrast_level = 6', 'contrast_level = 9')
    )

    assert config.reload_if_changed() is True
    assert config.get_int('display', 'contrast_level') == 9
    assert config.reload_if_changed() is False


def test_invalid_reload_keeps_last_valid_configuration():
    config = Config()
    config_file = Path(config.config_path)
    config_file.write_text('[display\ninvalid')

    assert config.reload_if_changed() is False
    assert config.get_int('display', 'contrast_level') == 6
