"""Calibration and Localite pages against recorded MATLAB answers (fixtures/bridge, made by
make_bridge_fixtures.m), so they are tested on exactly what nexus_bridge sends."""

import json
import os
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")

from PyQt6.QtCore import QObject, Qt, pyqtSignal  # noqa: E402
from PyQt6.QtWidgets import QApplication, QFileDialog, QMessageBox, QSplitter  # noqa: E402

from prestus_nexus.equipment import EquipmentLibrary  # noqa: E402
from prestus_nexus.model import ConfigModel  # noqa: E402
from prestus_nexus.session import Session  # noqa: E402
from prestus_nexus.ui.main_window import CALIBRATION, LOCALITE, PCT, MainWindow  # noqa: E402

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
        self.handlers: dict = {}            # fn -> args -> result, for facades without fixtures

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
        if fn in self.handlers:
            out = self.handlers[fn](args)
        elif fn == "nexus_calibration":
            name = {"open": "calibration_open", "load": "calibration_open",
                    "view": "calibration_view_predicted" if args.get("predicted") else "calibration_view_measured",
                    "fit": "calibration_fit_per_depth"
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
    assert calls == [] and page.depths.text() == "47.5"         # MATLAB starts on the first action
    assert page.coverage.text().startswith("Library entry: empirical") and page.entry_plot.arrays["cal_range"]
    assert "CTX_250_001" not in page.coverage.text()            # the selector names the device
    page.evaluate_measured()
    assert calls[0][1]["action"] == "open" and calls[0][1]["serial"] == "CTX_250_001"
    assert calls[0][1]["driving_system"] == "SC_105_010" and calls[0][1]["library_path"] == ""
    assert calls[1][1]["action"] == "view" and calls[1][1]["predicted"] is False   # profiles loaded
    assert page.depths.text() == "47.5"                         # the study's own depth
    assert page.measured.count() == 12 and page.coverage.text().startswith("Library entry:")
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


def test_calibration_depths_follow_the_setting_and_the_placement(window, monkeypatch):
    from prestus_nexus.ui.main_window import VIEWER
    viewer = window.pages[VIEWER]
    monkeypatch.setattr(viewer, "target_depth_ep", lambda ti: 52.04)
    page = window.show_special(CALIBRATION)
    assert page.depths.text() == "47.5, 52"                     # setting and placement target
    assert "transducer 1 target 52 mm" in page.log.toPlainText()
    page.predict()                                              # opens MATLAB (lazy)

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


def test_calibration_report_opens_the_entry_report(window, tmp_path, monkeypatch):
    # Report: enabled with a library entry; asks nexus_calibration for the entry's HTML
    # report and opens it; disabled, saying why, without an entry
    from PyQt6.QtGui import QDesktopServices
    page = window.show_special(CALIBRATION)
    page.refresh()
    assert page.entry is not None and page.report_btn.isEnabled()
    html = tmp_path / "CTX_250_001_SC_105_010_report.html"
    html.write_text("<html></html>")
    opened = []
    monkeypatch.setattr(QDesktopServices, "openUrl", lambda url: opened.append(url.toLocalFile()))
    calls = window.fake.calls
    fake_call = window.fake.call

    def call(fn, args=None, on_done=None, on_error=None, on_progress=None):
        if (args or {}).get("action") == "report":
            calls.append((fn, args))
            return on_done({"file": str(html), "library_file": "entry.yaml"})
        return fake_call(fn, args, on_done, on_error, on_progress)
    monkeypatch.setattr(window.fake, "call", call)
    page.open_report()
    assert calls[-1] == ("nexus_calibration", {"action": "report", "session": page.session_id})
    assert opened == [str(html)] and html.name in page.log.toPlainText()

    page.entry, page.info = None, {}                    # negative control: no entry
    page._update_buttons()
    assert not page.report_btn.isEnabled() and "save a fit first" in page.report_btn.toolTip()


def test_calibration_steps_sit_below_the_plots_and_fit_the_minimum_width(window):
    page = window.show_special(CALIBRATION)
    window.show()
    split = page.findChild(QSplitter)
    assert split.orientation() == Qt.Orientation.Vertical
    assert page.canvas.minimumHeight() >= 400               # two axes never squash
    window.resize(window.minimumWidth(), 700)
    QApplication.processEvents()
    for b in (page.load_btn, page.predict_btn, page.fit_btn, page.save_btn, page.kw_btn):
        right = b.mapTo(page, b.rect().topRight()).x()
        assert right <= page.width(), b.text()              # nothing clipped at the narrowest window


def test_cards_named_by_their_tab_drop_their_title(window):
    cards = {p.title.text(): p for tab in window.tab_pages.values() for p in tab.sections}
    assert cards["Study"].headless                              # same name as the Study tab
    assert cards["Where to run"].headless                       # same name as its sub-tab
    assert not cards["Pipeline stages"].headless                # shares Run ▸ Stages with Report
    assert not cards["Report"].headless
    assert not cards["Targeting analysis"].headless             # shares the Study tab
    assert not cards["Calibration conventions"].headless        # sits inside the Calibration page
    study = cards["Study"]
    assert study.title.isHidden() and study.toggle.isHidden() and not study.doc_button.isHidden()
    study.set_auto_collapse(lambda: True)                       # nothing to collapse to
    assert not study.body.isHidden()


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
    assert "Wrote 2 targets" in page.status.text()            # the recorded answer wrote both

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


def test_calibration_advanced_and_ring_options_reach_matlab(window, tmp_path, monkeypatch):
    from prestus_nexus.ui import calibration_page as cp
    page = window.show_special(CALIBRATION)
    planes = tmp_path / "planes.zip"
    planes.write_bytes(b"")

    def advanced_exec(dlg):                                     # preset, then two edits
        dlg.pre.setCurrentText("oneil_phase_only")
        dlg.nfm.setValue(0.2)
        dlg.grid.setText("auto")
        return 1

    def ring_exec(dlg):
        dlg.nsub.setValue(6)
        dlg.file.setText(str(planes))
        dlg.kind.setCurrentText("transverse")
        return 1
    monkeypatch.setattr(cp.AdvancedDialog, "exec", advanced_exec)
    monkeypatch.setattr(cp.RingOptionsDialog, "exec", ring_exec)
    page._advanced()
    assert (page.fwd.currentText(), page.huber.value(), page.fp_amp.isChecked()) == ("oneil", 10, False)
    assert page.advanced == {"Optimizer": "fexminimize", "Lambda": 0, "NearFieldPenalty": 0,
                             "NearFieldMargin": 0.2}
    assert page.extra == {"GridDims": "auto"} and page.adv_btn.text() == "Advanced (5)"
    page.ring.setCurrentIndex(page.ring.findData("fit"))
    page._ring_options()
    assert page.rpo_btn.text() == "Options (3)"
    page.fit()
    s = window.fake.calls[-1][1]["settings"]
    assert s["advanced"]["NearFieldMargin"] == 0.2 and s["huber_delta"] == 10
    ro = s["ring_profile_options"]
    assert ro["NSub"] == 6 and ro["OffAxis"] == str(planes) and ro["OffAxisOptions"]["Kind"] == "transverse"
    assert window.fake.calls[-1][1]["ring_profile"] == "fit"


def test_pct_page_predicts_and_prefers_the_default_run(window, tmp_path, monkeypatch):
    import nibabel as nib
    import numpy as np
    m = window.model
    m.set("subject_id", 3)
    m.set("path.pct", str(tmp_path / "pct"))
    m.set("path.sim", str(tmp_path / "sim"))
    m.set("simulation.medium", "layered")
    (tmp_path / "pct" / "sub-003").mkdir(parents=True)
    nib.save(nib.Nifti1Image(np.zeros((4, 4, 4), np.float32), np.eye(4)),
             str(tmp_path / "pct" / "sub-003" / "pseudoCT.nii.gz"))
    m.save(tmp_path / "config_study.yaml")
    window.dirty = False
    monkeypatch.setattr(MainWindow, "matlab_path", lambda self: "/fake/matlab")
    cache = tmp_path / "sim" / "sub-003" / "cache"
    prop = {"min": 1800, "max": 2200, "histogram": {"counts": [3, 5], "edges": [1800, 2000, 2200]}}
    const = {"min": 4.7, "max": 4.7, "histogram": {"counts": 8, "edges": [4.2, 5.2]}}

    def predict(args):
        cache.mkdir(parents=True, exist_ok=True)
        out = cache / "sub-003_layered_pct_prediction.json"
        out.write_text(json.dumps({"predicted": True, "source": "native", "n_voxels": 8,
                                   "properties": {"density": prop, "alpha_coeff": const},
                                   "reference": {"skull": {"density": 1850}}}))
        return {"file": str(out), "source": "native", "n_voxels": 8}
    window.fake.handlers["nexus_pct"] = predict

    page = window.show_special(PCT)
    assert page.predict_btn.isEnabled() and page.stats_file is None
    page.predict()
    fn, args = window.fake.calls[-1]
    assert fn == "nexus_pct" and args == {"action": "predict", "config": str(m.path)}
    assert page.stats_kind == "predicted" and "predicted from the CT (native voxels" in page.stats_title.text()
    att = page.axes[3]
    assert len(att.patches) >= 1 and any("constant 4.7" in t.get_text() for t in att.texts)
    assert page.fig.legends and "Bound variants (homogeneous skull)" in \
        [t.get_text() for t in page.fig.legends[0].get_texts()]

    (cache / "sub-003_layered_medium_stats.json").write_text(json.dumps(
        {"n_voxels": 8, "properties": {"density": prop}, "reference": {"skull": {"density": 1850}}}))
    page._load_stats()
    assert page.stats_kind == "real" and "mapped skull of the default run" in page.stats_title.text()
    assert "a prediction also exists" in page.stats_note.text()

    window.fake.fail["predict"] = ("No pseudoCT.nii.gz found", "pct_predict_skull_properties:notFound")
    page.predict()
    assert "Prediction failed: No pseudoCT.nii.gz found" in page.stats_note.text()
    assert page.predict_btn.isEnabled()


def test_pct_flow_steps_and_curves(window, tmp_path, monkeypatch):
    m = window.model
    m.set("simulation.medium", "layered")
    m.set("pct.enabled", 1)
    m.set("pct.skull_mapping", "wiesinger")
    monkeypatch.setattr(MainWindow, "matlab_path", lambda self: "/fake/matlab")
    hu = [-200, 1000, 3000]
    window.fake.handlers["nexus_pct"] = lambda args: {
        "hu": hu, "density": [1000, 1600, 2700], "sound_speed": [1500, 2300, 3600],
        "alpha_coeff": [0.5, 8, 20], "alpha_power": [1, 1, 1],
        "mapping": {"density": "k-plan", "sound_speed": "k-plan", "attenuation": "k-plan",
                    "density_calibration": "kplan_default"},
        "reference": {"skull": {"density": 1850, "sound_speed": 2800, "alpha_coeff": 13}}}
    flow = window.pct_flow
    flow.refresh()
    texts = [s.text.text() for s in flow.steps]
    assert [s.glyph for s in flow.steps] == ["mri", "histogram", "line_map", "ct", "table"]
    assert "UTE / PETRA image" in texts[0] and "HU = -2000 × UTE + 2042" in texts[2]
    assert flow.ute_canvas.isVisibleTo(flow) and flow.ute_fig.axes[0].lines
    flow.request_curves(force=True)
    fn, args = window.fake.calls[-1]
    assert fn == "nexus_pct" and args["action"] == "curves" and args["pct"]["skull_mapping"] == "wiesinger"
    assert len(flow.prop_fig.axes) == 3 and "Density (k-plan)" in flow.prop_fig.axes[0].get_title(loc="left")
    n = len(window.fake.calls)
    flow.refresh()                                     # same settings: drawn from the cache
    flow.request_curves(force=True)
    assert len(window.fake.calls) == n

    m.set("path.ct_pattern", "sub-%1$03d_CT.nii.gz")
    flow.refresh()
    assert "Measured CT" in flow.steps[0].text.text() and not flow.steps[2].isVisibleTo(flow)


def test_calibration_any_library_transducer(window, tmp_path, monkeypatch):
    """calibration_gui's lists: any library transducer and driving system, saved to a chosen
    folder without changing the config; the k-Wave check with its own simulation config."""
    page = window.show_special(CALIBRATION)
    page.any_lib.setChecked(True)
    assert page.lib_serial.isVisibleTo(page) and not page.tx.isVisibleTo(page)
    page.lib_serial.setCurrentIndex(page.lib_serial.findData("CTX_250_001"))
    page.lib_ds.setCurrentIndex(page.lib_ds.findData("SC_105_010"))
    assert page.current() == (-1, "CTX_250_001", "SC_105_010")
    assert page.depths.text() == ""                              # no study depths for a library device
    page.fit()
    calls = window.fake.calls
    assert calls[0][1]["action"] == "open" and calls[0][1]["serial"] == "CTX_250_001"
    lib = tmp_path / "own_cal"
    monkeypatch.setattr(QFileDialog, "getExistingDirectory", lambda *a, **k: str(lib))
    page.save()
    assert calls[-1][1]["action"] == "save" and calls[-1][1]["out_dir"] == str(lib)
    assert window.model.get("transducer.library_path") in (None, "", [])    # config unchanged
    assert page.library_path() == str(lib)

    monkeypatch.setattr(type(page), "_library_sim_config", lambda self: "")
    monkeypatch.setattr(QMessageBox, "exec", lambda self: 0)
    monkeypatch.setattr(QMessageBox, "clickedButton", lambda self: None)
    page.kwave_check()
    script = lib / "calibration" / "kwave_drive_CTX_250_001_SC_105_010.m"
    assert script.is_file() and "SimConfig" not in script.read_text()  # PRESTUS defaults
    page.any_lib.setChecked(False)
    assert page.tx.isVisibleTo(page) and page.current()[0] == 0
