"""Transducer tab ▸ Equipment (advanced): the sub-tab, the form and saving into an own library."""
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")

from PyQt6.QtCore import QObject, pyqtSignal  # noqa: E402
from PyQt6.QtWidgets import QApplication, QMessageBox  # noqa: E402

from prestus_nexus import yamlio  # noqa: E402
from prestus_nexus.equipment import EquipmentLibrary  # noqa: E402
from prestus_nexus.model import ConfigModel  # noqa: E402
from prestus_nexus.session import Session  # noqa: E402
from prestus_nexus.ui.main_window import EQUIPMENT, MainWindow  # noqa: E402


class FakeBridge(QObject):
    state_changed = pyqtSignal(str)
    log = pyqtSignal(str)

    def __init__(self):
        super().__init__()
        self.state, self.error, self.calls, self.issues = "ready", "", [], []

    def stop(self):
        pass

    def call(self, fn, args=None, on_done=None, on_error=None, on_progress=None):
        self.calls.append((fn, args))
        if (args or {}).get("action") == "layout":
            g = args["matrix"]["matrix_shape"]["define_here"]["grid_shape"]
            n = g["rect"]["elem_n_row"] * g["rect"]["elem_n_col"] if g["type"] == "rect" else 7
            on_done({"x_mm": [float(i) for i in range(n)], "y_mm": [0.0] * n, "z_mm": [0.0] * n, "n": n,
                     "error": ""})
            return
        on_done({"issues": list(self.issues)})


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def window(app, schema, prestus_root, tmp_path, monkeypatch):
    fake = FakeBridge()
    monkeypatch.setattr(MainWindow, "bridge", lambda self: fake)
    monkeypatch.setattr(QMessageBox, "question", lambda *a, **k: QMessageBox.StandardButton.Yes)
    m = ConfigModel(schema, EquipmentLibrary(prestus_root))
    m.set("transducer.serial", "CTX_250_001")
    m.save(tmp_path / "config_study.yaml")
    w = MainWindow(m, Session(tmp_path / "session.yaml"), prestus_root)
    w.fake = fake
    yield w
    w.dirty = False
    w.close()


def _tab_visible(window) -> bool:
    tab = window.tab_pages["transducer"]
    return tab.subtabs.isTabVisible(tab.subtabs.indexOf(window.pages[EQUIPMENT]))


def test_equipment_subtab_follows_advanced(window):
    tab = window.tab_pages["transducer"]
    tab.advanced.setChecked(False)
    assert not _tab_visible(window)
    tab.advanced.setChecked(True)
    assert _tab_visible(window)


def test_copy_geometry_and_save_new_transducer(window, tmp_path):
    page = window.pages[EQUIPMENT]
    page.refresh()
    page.start.setCurrentText("CTX_250_001")
    assert page.serial.text() == ""                      # a copy needs a new serial
    assert page.idd.text() and page.od.text()
    assert not page.save_btn.isEnabled()                 # no serial yet
    page.serial.setText("CTX_250_777")
    assert page.save_btn.isEnabled()
    lib = tmp_path / "own_lib"
    lib.mkdir()
    page.folder.setCurrentText(str(lib))
    page.save()
    assert window.fake.calls and window.fake.calls[-1][0] == "nexus_equipment"
    entry = yamlio.load_file(lib / "tran-CTX_250_777_transducer.yaml")
    assert entry["serial"] == "CTX_250_777" and entry["transducer"]["annular"]["elem_n"] >= 1


def test_add_driving_system_to_shipped_transducer(window, tmp_path):
    page = window.pages[EQUIPMENT]
    page.refresh()
    page.keep_serial.setChecked(True)
    page.start.setCurrentText("CTX_250_001")
    assert page.serial.text() == "CTX_250_001"
    n_before = len(page.rows)
    page.add_combo()
    row = page.rows[-1]
    row.ds.setCurrentText("MY_NEW_DS")
    assert row.is_new() and row.new_box.isVisibleTo(page)
    row.ds_name.setText("Lab driver")
    lib = tmp_path / "own_lib2"
    lib.mkdir()
    page.folder.setCurrentText(str(lib))
    page.save()
    entry = yamlio.load_file(lib / "tran-CTX_250_001_transducer.yaml")
    assert "MY_NEW_DS" in entry["combos"] and len(entry["combos"]) == n_before + 1
    ds = yamlio.load_file(lib / "ds-MY_NEW_DS_driver.yaml")
    assert ds["type"] == "driving_system" and ds["name"] == "Lab driver"


def test_never_saves_into_shipped_equipment(window, prestus_root):
    page = window.pages[EQUIPMENT]
    page.refresh()
    page.start.setCurrentText("CTX_250_001")
    page.serial.setText("CTX_250_888")
    page.folder.setCurrentText(str(prestus_root / "config" / "equipment"))
    page.save()
    assert not (prestus_root / "config" / "equipment" / "tran-CTX_250_888_transducer.yaml").exists()
    assert "shipped config/equipment is not changed" in page.msg.toPlainText()


def test_matrix_transducer_layout_from_matlab_and_save(window, tmp_path):
    page = window.pages[EQUIPMENT]
    page.refresh()
    page.start.setCurrentIndex(0)                         # blank
    page.kind.setCurrentText("matrix")
    page.serial.setText("MX_TEST")
    page.m_grid.setCurrentText("rect")
    page.r_nrow.setValue(4)
    page.r_ncol.setValue(5)
    assert not page.save_btn.isEnabled()                  # layout not known yet
    page._request_layout()
    assert window.fake.calls[-1][1]["action"] == "layout"
    assert page.save_btn.isEnabled()
    assert "20 elements" in page.msg.toPlainText()
    lib = tmp_path / "own_mx"
    lib.mkdir()
    page.folder.setCurrentText(str(lib))
    page.save()
    entry = yamlio.load_file(lib / "tran-MX_TEST_transducer.yaml")
    assert entry["n_elem"] == 20 and entry["transducer"]["matrix"]["matrix_shape"]["define_here"]["grid_shape"][
        "rect"]["elem_n_col"] == 5


def test_shipped_save_for_maintainers(window, tmp_path, monkeypatch):
    page = window.pages[EQUIPMENT]
    shipped = tmp_path / "equipment"
    (shipped / "ds").mkdir(parents=True)
    monkeypatch.setattr(type(page), "shipped_dir", lambda self: shipped)
    page.refresh()
    page.start.setCurrentText("CTX_250_001")
    page.serial.setText("CTX_250_555")
    page.shipped.setChecked(True)
    page.save()
    assert (shipped / "CTX_250_555.yaml").is_file()


def test_view_geometry_from_calibration(window):
    from prestus_nexus.ui.main_window import CALIBRATION
    cal = window.show_special(CALIBRATION)
    window.tab_pages["transducer"].advanced.setChecked(False)
    cal._view_geometry()
    page = window.pages[EQUIPMENT]
    assert _tab_visible(window) and page.start.currentData() == "CTX_250_001"
    assert page.idd.text() and page.serial.text() == ""      # shown as a copy: nothing is changed


@pytest.mark.parametrize("kind", ["annular", "matrix"])
def test_preview_frame_holds_while_editing_radius_and_exit_plane(window, kind):
    """Changing the curvature radius or the exit plane distance keeps the side view's scale and
    apex on screen; loading another transducer, a resize, or an overflow refits it."""
    page = window.pages[EQUIPMENT]
    page.refresh()
    page.start.setCurrentText("CTX_250_001")
    page.kind.setCurrentText(kind)
    sch = page.schematic
    sch.resize(600, 300)
    frame = lambda: (sch.grab(), sch._view[2:])[1]     # noqa: E731  (paint, then the held frame)
    held = frame()
    for spin, delta in ((page.rc, 5.0), (page.rc, -5.0), (page.dep, 4.0), (page.dep, -4.0)):
        spin.setValue(spin.value() + delta)
        assert frame() == held, (spin is page.rc, delta)
    page.rc.setValue(page.rc.value() * 6)              # far past the frame: refits (negative control)
    assert frame() != held
    page.start.setCurrentText("CTX_250_001")            # a loaded transducer refits
    refit = frame()
    sch.resize(800, 300)                                # a resize refits
    assert frame() != refit
    img = sch.grab()                                    # export painting is fitted, frame untouched
    assert not img.isNull()


def test_preview_refits_when_the_drawing_shrinks_to_under_half(window):
    page = window.pages[EQUIPMENT]
    page.refresh()
    page.start.setCurrentText("CTX_250_001")
    sch = page.schematic
    sch.resize(600, 300)
    sch.grab()
    page.rc.setValue(page.rc.value() * 4)               # zoom out (refit)
    page.dep.setValue(page.rc.value() - 3)
    sch.grab()
    held = sch._view[2]
    page.rc.setValue(page.rc.value() / 4)               # back to a small drawing: refits, zooms in
    page.dep.setValue(page.rc.value() - 3)
    sch.grab()
    assert sch._view[2] > held
