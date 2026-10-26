import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest
from PySide6.QtWidgets import QApplication, QWidget
from pytest import approx

from gui.image_panel import ImagePanel
from gui.spectrum_panel import SpectrumPanel
from app import main as main_module


def test_image_panel_renders_with_qt_canvas(qapp, synthetic_image):
    parent = QWidget()
    panel = ImagePanel(parent, on_clear=lambda: None, on_sa100_points=None)

    panel.show_image(synthetic_image.data)
    panel.set_slider_range(-10, 10, -5, 5)

    assert panel.canvas.parent() is panel
    assert panel._axes_image is not None
    assert panel._axes_image.get_array().shape == synthetic_image.data.shape
    assert len(panel.figure.axes) == 2
    assert panel.cut_controls.parent() is panel

    cuts = []
    def on_cut_changed(low, high):
        cuts.append((low, high))
        panel.apply_cut(low, high)

    panel.set_cut_callback(on_cut_changed)
    panel.slider_low.setValue(3000)
    assert cuts[0][0] == approx(-4)
    assert cuts[0][1] is None
    assert panel._axes_image.norm.vmin == approx(-4)
    panel.slider_high.setValue(8000)
    assert cuts[1][0] is None
    assert cuts[1][1] == approx(6)
    assert panel._axes_image.norm.vmax == approx(6)

    panel.clear()
    panel.show_image(synthetic_image.data)
    assert len(panel.figure.axes) == 2
    parent.close()


def test_spectrum_panel_renders_with_qt_canvas(qapp, flat_spectrum):
    parent = QWidget()
    panel = SpectrumPanel(
        parent,
        on_clear=lambda: None,
        on_show_lines=lambda: None,
        on_colorize=lambda: None,
    )

    panel.show_spectrum("test", flat_spectrum, calibrated=True)

    assert panel.canvas.parent() is panel
    assert len(panel.axe.lines) == 1
    assert panel.axe.get_xlabel() == "Wavelength (Angstrom)"
    assert panel._contrast_color() == "black"
    panel.axe.set_facecolor("black")
    assert panel._contrast_color() == "white"
    parent.close()


def test_application_has_classic_menus(qapp, monkeypatch):
    monkeypatch.setattr(main_module.LogHandler, 'initialize', lambda _self: None)
    monkeypatch.setattr(main_module.os_utils, 'show_versions', lambda: None)
    window = main_module.Application("QuickSpec", "test")
    try:
        assert window._theme == 'dark'
        assert qapp.palette().color(
            qapp.palette().ColorRole.Window
        ).red() == 53
        assert main_module.plt.rcParams['figure.facecolor'] == (
            main_module.plt.style.library['dark_background']['figure.facecolor']
        )
        window.conf.config.set('display', 'theme', 'light')
        window._apply_theme()
        assert window._theme == 'light'
        assert qapp.palette().color(
            qapp.palette().ColorRole.Window
        ).red() == 240
        assert main_module.plt.rcParams['figure.facecolor'] == (
            main_module.plt.style.library['grayscale']['figure.facecolor']
        )
        window.conf.config.set('display', 'theme', 'dark')
        window._apply_theme()
        menus = {
            action.text(): action.menu()
            for action in window.menuBar().actions()
        }
        assert list(menus) == ["File", "Longslit", "Slitless"]
        assert [action.text() for action in menus["File"].actions()] == [
            "Load",
            "",
            "Exit",
        ]
        assert menus["File"].actions()[0].shortcut().toString() == "Ctrl+O"
        assert menus["File"].actions()[2].shortcut().toString() == "Ctrl+X"
        longslit_actions = menus["Longslit"].actions()
        assert longslit_actions[0].text() == "Process auto"
        assert longslit_actions[0].shortcut().toString() == "Ctrl+P"
        assert longslit_actions[1].isSeparator()
        assert [action.text() for action in longslit_actions[2:]] == [
            "Find spectrum",
            "Extract spectrum",
            "Calibrate spectrum",
            "Apply response",
            "Smooth, crop & normalize",
        ]
        assert [action.text() for action in menus["Slitless"].actions()] == [
            "3-click process"
        ]
        assert menus["Slitless"].actions()[0].shortcut().toString() == "Ctrl+Q"

        monitor_stopped = []
        monkeypatch.setattr(
            window.file_monitor,
            "stop",
            lambda: monitor_stopped.append(True),
        )
        window.show()
        menus["File"].actions()[2].trigger()
        qapp.processEvents()
        assert not window.isVisible()
        assert monitor_stopped == [True]
    finally:
        window.close()

@pytest.fixture(scope="session")
def qapp():
    return QApplication.instance() or QApplication([])
