"""Transducer ▸ Device panel: the device list shows the study's library device, and none with
custom geometry (not the list's first entry)."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")


@pytest.fixture
def panel(schema, prestus_root):
    from PyQt6.QtWidgets import QApplication

    from prestus_nexus.equipment import EquipmentLibrary
    from prestus_nexus.model import ConfigModel
    from prestus_nexus.ui.transducer_panel import DevicePanel
    app = QApplication.instance() or QApplication([])
    m = ConfigModel(schema, EquipmentLibrary(prestus_root))
    p = DevicePanel(m)
    p.app = app
    return p


def test_custom_geometry_shows_no_device(panel):
    panel.model.reset("transducer.serial", 0)
    panel.refresh(0)
    assert panel.rb_custom.isChecked()
    assert panel.device.currentIndex() == -1 and not panel.device.isEnabled()


def test_switching_to_the_library_picks_and_sets_a_device(panel):
    panel.model.reset("transducer.serial", 0)
    panel.refresh(0)
    panel.rb_lib.setChecked(True)
    assert panel.device.currentIndex() == 0
    assert panel.model.get("transducer.serial", 0) == panel.device.itemData(0)


def test_library_device_is_shown(panel):
    serial = "CTX_500_026"
    panel.model.set("transducer.serial", serial, 0)
    panel.refresh(0)
    assert panel.rb_lib.isChecked() and panel.device.currentData() == serial
    panel.model.reset("transducer.serial", 0)               # back to custom: the device clears
    panel.refresh(0)
    assert panel.device.currentIndex() == -1
