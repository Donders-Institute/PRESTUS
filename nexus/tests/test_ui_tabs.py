"""Tabbed main window: tabs, advanced options, navigation, timing panel, Localite import."""

import json
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")

from PyQt6.QtWidgets import QApplication  # noqa: E402

from prestus_nexus.model import ConfigModel  # noqa: E402
from prestus_nexus.session import Session  # noqa: E402
from prestus_nexus.ui.main_window import STUDY, MainWindow  # noqa: E402


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def window(app, schema, prestus_root, tmp_path):
    m = ConfigModel(schema)
    m.set("path.sim", str(tmp_path / "sim"))
    m.save(tmp_path / "config_t.yaml")
    w = MainWindow(m, Session(tmp_path / "s.yaml"), prestus_root)
    yield w
    w.dirty = False
    w.close()


def _row(page, path):
    return next(r for r in page.rows if r.node.path == path)


def test_tabs_and_sections(window):
    keys = [window.tabs.widget(i).key for i in range(window.tabs.count())]
    assert keys == ["setup", "transducer", "placement", "simulation", "study", "run", "results"]
    assert window._section_tab["timing"] == "simulation" and window._section_tab["study"] == "study"
    # platform is shown on the Run tab (HPC card), the backends on the Simulation tab
    assert any(r.node.path == "platform" for r in window.pages["hpc"].rows)
    assert any(r.node.path == "acoustic_backend" for r in window.pages["simulation"].rows)
    assert not any(r.node.path == "platform" for r in window.pages["general"].rows)


def test_advanced_options_toggle_and_override(window):
    tab = window.tab_pages["simulation"]
    window.tabs.setCurrentWidget(tab)
    grid = window.pages["grid"]
    assert not _row(grid, "grid.pml_size").block.isVisibleTo(grid)       # advanced: hidden
    assert _row(grid, "grid.resolution_mm").block.isVisibleTo(grid)      # basic: shown
    assert "advanced hidden" in tab.hidden.text()
    tab.advanced.setChecked(True)
    assert _row(grid, "grid.pml_size").block.isVisibleTo(grid)
    tab.advanced.setChecked(False)
    window.model.set("grid.pml_size", 20)                                 # set advanced -> stays visible
    tab.refresh()
    assert _row(grid, "grid.pml_size").block.isVisibleTo(grid)


def test_card_advanced_button_opens_one_card(window):
    tab = window.tab_pages["simulation"]
    window.tabs.setCurrentWidget(tab)
    grid, timing = window.pages["grid"], window.pages["timing"]
    assert grid.adv_button.isVisibleTo(grid) and "advanced" in grid.adv_button.text()
    grid.adv_button.click()
    assert grid.show_advanced and not timing.show_advanced and not tab.advanced.isChecked()
    assert _row(grid, "grid.pml_size").block.isVisibleTo(grid)
    assert not _row(timing, "timing.pt_timestep").block.isVisibleTo(timing)
    assert grid.adv_button.text() == "Hide advanced"
    grid.adv_button.click()
    assert not _row(grid, "grid.pml_size").block.isVisibleTo(grid)


def test_tune_for_goal_reveals_and_marks_keys(window):
    tab = window.tab_pages["simulation"]
    window.tabs.setCurrentWidget(tab)
    grid = window.pages["grid"]
    assert tab.tune.isVisibleTo(tab) and not (res := window.tab_pages["results"]).tune.isVisibleTo(res)
    tab.tune.setCurrentIndex(tab.tune.findData("speed"))
    assert window.tune_goal == "speed"
    assert window.tab_pages["setup"].tune.currentData() == "speed"       # shared by all tabs
    assert _row(grid, "grid.pml_size").block.isVisibleTo(grid)          # advanced, tagged speed
    assert not _row(grid, "grid.use_kWaveArray").block.isVisibleTo(grid)  # advanced, accuracy only
    assert "speed" in _row(grid, "grid.resolution_mm").label.text()     # basic fields are marked too
    assert not grid.show_advanced
    window.set_tune_goal(None)
    assert not _row(grid, "grid.pml_size").block.isVisibleTo(grid)
    assert "speed" not in _row(grid, "grid.resolution_mm").label.text()


def test_goto_reveals_advanced_field(window):
    window.goto("headmodel.smooth_fwhm_mm")
    tab = window.tab_pages["simulation"]
    assert window.tabs.currentWidget() is tab
    page = window.pages["headmodel"]
    assert page.show_advanced and _row(page, "headmodel.smooth_fwhm_mm").block.isVisibleTo(page)


def test_labels_and_descriptions(window):
    row = _row(window.pages["timing"], "timing.pri")
    assert "Pulse repetition interval" in row.label.text()
    assert "Key: timing.pri" in row.label.toolTip()


def test_timing_panel_export_and_import(window, tmp_path):
    for k, v in {"timing.pd": 0.02, "timing.pri": 0.2, "timing.ptd": 1, "transducer.freq_hz": 5e5}.items():
        window.model.set(k, v)
    panel = window.timing_panel
    panel.refresh()
    out = panel.export_settings(tmp_path / "calc.json")
    data = json.loads(out.read_text())
    assert data["Pulse Duration [ms]"] == 20 and data["Fundamental Frequency [MHz]"] == 0.5
    data["Pulse Duration [ms]"] = 30
    data["Isppa [W/cm2]"] = 25
    (tmp_path / "calc2.json").write_text(json.dumps(data))
    assert panel.import_settings(tmp_path / "calc2.json", apply_optional=True)
    assert window.model.get("timing.pd") == pytest.approx(0.03)
    assert window.model.get("transducer.target_isppa_wcm2") == 25
    assert panel.diagram.p.pd == pytest.approx(0.03)


def test_localite_dialog_builds_bilateral_runs(window, tmp_path):
    from prestus_nexus.ui.localite_import import LocaliteImportDialog
    post = tmp_path / "post"
    post.mkdir()
    (post / "neuronav_ingest_summary.csv").write_text(
        "sub,ses,target\nsub-001,ses-02,R_MD\nsub-001,ses-02,L_MD\nsub-002,ses-02,L_MD\n")
    window.model.set("path.localite_post", str(post))
    window.model.set("placement.mode", "localite")
    window.model.add_transducer()
    dlg = LocaliteImportDialog(window)
    assert dlg.g_session.isChecked()                        # 2 transducers -> per-session runs
    rows = dlg.result_manifest.rows
    assert [r["placement.localite.target"] for r in rows] == ["[L_MD, R_MD]"]
    assert "sub-002" in dlg.notes.text()
    dlg.g_target.setChecked(True)
    assert len(dlg.result_manifest.rows) == 3
    study = window.show_special(STUDY)
    assert window.tabs.currentWidget().key == "study" and study.isVisibleTo(window)
