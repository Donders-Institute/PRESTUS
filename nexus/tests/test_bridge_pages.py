"""Calibration and Localite pages against recorded MATLAB answers (fixtures/bridge, made by
make_bridge_fixtures.m), so they are tested on exactly what nexus_bridge sends."""

import json
import os
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")

from PyQt6.QtCore import QObject, Qt, pyqtSignal  # noqa: E402
from PyQt6.QtWidgets import QApplication, QFileDialog, QMessageBox  # noqa: E402

from prestus_nexus.equipment import EquipmentLibrary  # noqa: E402
from prestus_nexus.model import ConfigModel  # noqa: E402
from prestus_nexus.session import Session  # noqa: E402
from prestus_nexus.ui.main_window import CALIBRATION, LOCALITE, MainWindow  # noqa: E402

FIX = Path(__file__).parent / "fixtures" / "bridge"


def fixture(name: str):
    return json.loads((FIX / f"{name}.json").read_text(encoding="utf-8"))


class FakeBridge(QObject):
    state_changed = pyqtSignal(str)
    log = pyqtSignal(str)

    def __init__(self):
        super().__init__()
        self.state, self.error, self.calls = "ready", "", []
        self.fail: dict[str, tuple[str, str]] = {}

    def stop(self):
        pass

    def call(self, fn, args=None, on_done=None, on_error=None, on_progress=None):
        from prestus_nexus.bridge import BridgeError
        args = args or {}
        self.calls.append((fn, args))
        action = args.get("action")
        if action in self.fail:
            on_error(BridgeError(*self.fail[action]))
            return
        if fn == "nexus_calibration":
            name = {"open": "calibration_open", "load": "calibration_open",
                    "view": "calibration_view_predicted" if args.get("predicted") else "calibration_view_measured",
                    "fit": "calibration_fit_static" if args.get("mode") == "static" else "calibration_fit_per_depth"
                    }.get(action)
            out = fixture(name) if name else {"file": str(Path(args.get("out_dir", "")) / "entry.yaml"),
                                              "path": args.get("path", ""), "what": "fit"}
        else:
            out = fixture(f"neuronav_{action}") if action != "sessions" else {}
        if on_progress is not None:
            on_progress("working")
        on_done(out)


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def window(app, schema, prestus_root, tmp_path, monkeypatch):
    fake = FakeBridge()
    monkeypatch.setattr(MainWindow, "bridge", lambda self: fake)
    m = ConfigModel(schema, EquipmentLibrary(prestus_root))
    m.set("transducer.serial", "CTX_250_001")
    m.set("transducer.combo.ds_serial", "SC_105_010")
    m.set("transducer.focal_distance_ep", 47.5)
    m.set("path.localite_raw", str(tmp_path / "raw"))
    m.set("path.localite_post", str(tmp_path / "post"))
    (tmp_path / "raw" / "sub-001" / "ses-01").mkdir(parents=True)
    m.save(tmp_path / "config_study.yaml")
    w = MainWindow(m, Session(tmp_path / "session.yaml"), prestus_root)
    w.fake = fake
    yield w
    w.dirty = False
    w.close()


def test_calibration_page_opens_the_config_transducer_and_fits(window):
    page = window.show_special(CALIBRATION)
    calls = window.fake.calls
    assert calls[0][1]["action"] == "open" and calls[0][1]["serial"] == "CTX_250_001"
    assert calls[0][1]["driving_system"] == "SC_105_010" and calls[0][1]["library_path"] == ""
    assert calls[1][1]["action"] == "view" and calls[1][1]["predicted"] is False   # profiles loaded
    assert page.depths.text() == "47.5"                         # the study's own depth
    assert page.measured.count() == 12 and "CTX_250_001_SC_105_010" in page.coverage.text()
    assert page.table.rowCount() == 0                           # no fit yet, measured view
    assert page.ax_r2.get_ylabel().startswith("R²")

    page.predict()
    assert calls[-1][1]["depths"] == [47.5] and page.ax_r2.get_ylabel().startswith("Focus")
    assert page.table.rowCount() == 3 and page.table.horizontalHeaderItem(2).text() == "Entry"
    assert page.show_depth.itemText(0) == "All depths"
    page.show_depth.setCurrentIndex(1)                          # one depth: R² in the legend
    assert any("R²" in t.get_text() for t in page.ax_prof.get_legend().get_texts())

    for i in range(page.measured.count()):
        page.measured.item(i).setCheckState(Qt.CheckState.Unchecked)
    page.fit()
    assert "Select at least one" in page.log.toPlainText().splitlines()[-1]
    page.measured.item(4).setCheckState(Qt.CheckState.Checked)
    page.fit()
    fit_args = calls[-1][1]
    assert fit_args["action"] == "fit" and fit_args["measured_depths"] == [31.11]
    assert fit_args["settings"]["isppa_wcm2"] == 30 and fit_args["free_params"] == ["phase", "amp"]
    assert page.has_fit and page.table.horizontalHeaderItem(0).text() == "Element"

    page.mode.setCurrentIndex(1)                                # static: one row per element
    page.fit()
    assert page.table.horizontalHeaderItem(2).text() == "Gain"
    assert "c_steer" in page.log.toPlainText()


def test_calibration_depths_follow_the_setting_and_the_placement(window, monkeypatch):
    from prestus_nexus.ui.main_window import VIEWER
    viewer = window.pages[VIEWER]
    monkeypatch.setattr(viewer, "target_depth_ep", lambda ti: 52.04)
    page = window.show_special(CALIBRATION)
    assert page.depths.text() == "47.5, 52"                     # setting and placement target
    assert "transducer 1 target 52 mm" in page.log.toPlainText()

    window.model.set("transducer.focal_distance_ep", 45.0)      # config changed: follows
    page.refresh()
    assert page.depths.text() == "45, 52" and window.fake.calls[-1][1]["depths"] == [45.0, 52.0]

    page.depths.setText("30, 40")                               # typed by hand: kept
    window.model.set("transducer.focal_distance_ep", 50.0)
    page.refresh()
    assert page.depths.text() == "30, 40"
    page.use_study_depths()                                     # the button refills on request
    assert page.depths.text() == "50, 52"


def test_calibration_save_sets_library_path_and_k_wave_script(window, tmp_path, monkeypatch):
    page = window.show_special(CALIBRATION)
    page.fit()
    lib = tmp_path / "mylib"
    monkeypatch.setattr(QFileDialog, "getExistingDirectory", lambda *a, **k: str(lib))
    page.save()
    assert window.model.get("transducer.library_path") == str(lib)
    assert window.fake.calls[-1][1] == {"action": "save", "session": page.session_id,
                                        "out_dir": str(lib), "allow_shipped": False}
    assert "Saved" in page.log.toPlainText()

    monkeypatch.setattr(QMessageBox, "exec", lambda self: 0)
    monkeypatch.setattr(QMessageBox, "clickedButton", lambda self: None)
    page.kwave_check()
    script = tmp_path / "calibration" / "kwave_drive_CTX_250_001_SC_105_010.m"
    assert window.fake.calls[-1][1]["action"] == "write_fit" and script.is_file()
    text = script.read_text()
    assert "kwave_calibrate_drive(fit.entry, fit.prof, ko)" in text and str(lib) in text
    assert "calibration_save_entry(lib, fit" in text


def test_calibration_page_reopens_after_matlab_restart(window):
    page = window.show_special(CALIBRATION)
    window.fake.fail["view"] = ("No calibration session", "nexus_calibration:noSession")
    page.predict()
    assert page.opened is None
    del window.fake.fail["view"]
    n = len(window.fake.calls)
    page.predict()
    assert [c[1]["action"] for c in window.fake.calls[n:]] == ["open", "view"]


def test_localite_page_detect_write_and_use(window, monkeypatch):
    page = window.show_special(LOCALITE)
    assert page.sub.currentText() == "sub-001" and page.ses.currentText() == "ses-01"
    assert "Raw exports" in page.summary.text()
    page.detect()
    args = window.fake.calls[-1][1]
    assert args["action"] == "detect" and args["config"].endswith("config_study.yaml")
    assert page.table.rowCount() == 2 and page.table.item(0, 3).text() == "LEFT PUL"
    assert "1.0 mm off axis" in page.table.item(0, 4).text()
    assert page.write_btn.isEnabled() and page.raw.rowCount() == 2

    page.table.item(0, 3).setText("L_PUL")
    page.table.item(1, 0).setCheckState(Qt.CheckState.Unchecked)
    page.write()
    tm = window.fake.calls[-1][1]["target_map"]
    assert tm == [{"series_index": 1, "transducer_id": 0, "target_name": "L_PUL", "color": ""}]
    assert "Wrote 2 target(s)" in page.status.text()            # the recorded answer wrote both

    page.use_in_study()
    assert window.model.get("placement.localite.session") == "ses-01"
    assert window.model.get("placement.localite.target") == "L_PUL"

    window.model.set("placement.localite.voxel_size", 1.0)      # stale detection: no write
    page.refresh()
    assert not page.write_btn.isEnabled()


def test_localite_batch_saves_label_map_and_runs_script(window, monkeypatch):
    from prestus_nexus.ui.localite_page import BatchDialog
    from prestus_nexus.ui.run_page import RunPage
    ran = []
    monkeypatch.setattr(RunPage, "run_script", lambda self, task, script: ran.append((task, script)) or True)
    dlg = BatchDialog(window)
    assert dlg.sessions.rowCount() == 1 and dlg.labels.rowCount() == 2
    dlg.labels.item(0, 4).setText("L_PUL")
    dlg.write_mode.setCurrentText("named")
    dlg.run()
    assert window.model.get("neuronav.target_labels") == [{"target": "L_PUL", "labels": ["LEFT PUL"]}]
    assert window.model.get("neuronav.batch_write") == "named"
    task, script = ran[0]
    assert "1 sessions" in task and script.name == "ingest_study.m"
    text = script.read_text()
    assert "selection = {'sub-001', 'ses-01'};" in text and "neuronav_ingest_batch(parameters, selection" in text
