"""Tabbed main window: tabs, advanced options, navigation, timing panel, Localite import."""

import json
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")

from PyQt6.QtWidgets import QApplication  # noqa: E402

from prestus_nexus.model import ConfigModel  # noqa: E402
from prestus_nexus.session import Session  # noqa: E402
from prestus_nexus.uncertainty import PROPERTIES  # noqa: E402
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
    assert keys == ["start", "setup", "transducer", "placement", "simulation", "study", "run", "results"]
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


def test_tissue_properties_locked_by_default(window, monkeypatch):
    from PyQt6.QtWidgets import QMessageBox
    page = window.pages["medium_properties"]
    row = _row(page, "medium_properties.skull.sound_speed")
    assert page.locked and not row.editor.isEnabled() and not row.reset.isEnabled()
    assert not window.pages["grid"].lock_button.isVisibleTo(window.pages["grid"])
    monkeypatch.setattr(QMessageBox, "warning", lambda *a, **k: QMessageBox.StandardButton.Cancel)
    page.lock_button.click()                                              # declined: stays locked
    assert page.locked and not row.editor.isEnabled()
    monkeypatch.setattr(QMessageBox, "warning", lambda *a, **k: QMessageBox.StandardButton.Yes)
    page.lock_button.click()
    assert not page.locked and row.editor.isEnabled()
    page.lock_button.click()                                              # relocking needs no confirmation
    assert page.locked and not row.editor.isEnabled()


def test_uncertainty_card_lists_variants(window):
    card = window.uncertainty_card
    window.tabs.setCurrentWidget(window.tab_pages["simulation"])
    assert not card.isVisibleTo(window)                                   # uncertainty off
    window.model.set("simulation.uncertainty", True)
    window.tab_pages["simulation"].refresh()
    assert card.isVisibleTo(window)
    assert not card.sets.errors
    n_tissues = len(card.sets.sets["default"])
    assert card.table.rowCount() == n_tissues * len(PROPERTIES)
    if card.canvases:
        assert len(card.canvases["brain"].figure.axes) == 7               # one column per plotted property
        assert len(card.canvases["skull"].figure.axes) == 7


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


def test_placement_modes_greyed_without_segmentation(window, tmp_path):
    anat = tmp_path / "anat"
    anat.mkdir()
    (anat / "sub-003_T1w.nii.gz").write_bytes(b"")
    window.model.set("subject_id", 3)
    window.model.set("path.anat", str(anat))
    page = window.pages["placement"]
    page.refresh()
    ed = _row(page, "placement.mode").editor
    items = {v: ed.combo.model().item(i) for i, v in enumerate(ed._items)}
    assert items["manual"].isEnabled()
    assert not items["heuristic"].isEnabled()                 # no m2m, no SimNIBS
    assert "unavailable" in items["heuristic"].text()
    assert "simnibs_bin_path" in items["plantus"].toolTip()
    window.model.set("placement.mode", "heuristic")             # current value stays selectable
    page.refresh()
    ed = _row(page, "placement.mode").editor
    assert ed.combo.model().item(ed._items.index("heuristic")).isEnabled()


def test_plantus_actions_run_placement_only(window, monkeypatch):
    from prestus_nexus.ui.main_window import RUN, VIEWER, PlantusActions
    from prestus_nexus.ui.run_page import RunPage
    ran = []
    monkeypatch.setattr(RunPage, "run_script", lambda self, task, script: ran.append((task, script)) or True)
    acts = window.pages["placement"].findChild(PlantusActions)
    window.model.set("placement.mode", "manual")
    acts.refresh()
    assert acts.isHidden()
    window.model.set("placement.mode", "plantus")
    acts.refresh()
    assert not acts.isHidden()
    acts._run(True)
    task, script = ran[0]
    assert task == PlantusActions.TASK and script.name == "plantus_t.m"
    assert "'placement.plantus.interactive', 1;" in script.read_text()
    shown = []
    monkeypatch.setattr(type(window), "show_special", lambda self, key: shown.append(key))
    window.pages[RUN].finished.emit(PlantusActions.TASK, True)
    assert shown == [VIEWER] and not acts._waiting
def test_choice_notes_never_become_windows(window, app):
    """A choice's plain-text note (overlay choice_math '> …') stays inside its field: replacing it
    on refresh used to leave each old note as a top-level window (the pop-up seen in advanced mode)."""
    from PyQt6.QtWidgets import QApplication
    window.model.set("thermal.cem43_iso", 1)
    window.show()
    tab = window.tab_pages["sonication"]
    window.tabs.setCurrentWidget(tab)
    for on in (True, False, True):
        tab.advanced.setChecked(on)
        window._changed()
        QApplication.processEvents()
    window.model.set("thermal.cem43_iso", 0)
    window._changed()
    QApplication.processEvents()
    strays = [w for w in QApplication.topLevelWidgets() if w.isVisible() and w is not window]
    assert not strays, [type(w).__name__ for w in strays]


