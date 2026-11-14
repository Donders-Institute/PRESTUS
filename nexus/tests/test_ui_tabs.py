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
from prestus_nexus.ui.uncertainty_view import FILE_KEYS  # noqa: E402


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
    assert keys == ["start", "setup", "transducer", "placement", "medium", "sonication", "solver",
                    "study", "run", "results"]
    assert window._section_tab["timing"] == "sonication" and window._section_tab["study"] == "study"
    assert window._section_tab["medium_properties"] == "medium" and window._section_tab["grid"] == "solver"
    # the medium choice is the Medium tab's header row; the uncertainty switch has its own card under the properties
    assert window.pages["__medium__"] is window.medium_selector
    assert not any(r.node.path == "simulation.medium" for p in ("simulation", "medium_properties")
                   for r in window.pages[p].rows)
    unc = next(r for r in window.pages["uncertainty"].rows if r.node.path == "simulation.uncertainty")
    assert not unc.node.advanced and window._section_tab["uncertainty"] == "medium"
    tab = window.tab_pages["medium"]
    assert [tab.subtabs.tabText(i) for i in range(tab.subtabs.count())] == \
        ["Head model", "Acoustic properties", "Uncertainty", "Skull CT"]
    assert tab.subtab_of(window.pages["segmentation"]) is window.medium_subtabs["head"]
    assert tab.subtab_of(window.pages["medium_properties"]) is window.medium_subtabs["acoustic"]
    assert tab.subtab_of(window.pages["uncertainty"]) is window.medium_subtabs["uncertainty"]
    assert tab.subtab_of(window.pages["pct"]) is window.medium_subtabs["pct"]
    assert any(r.node.path == "pct.enabled" for r in window.pages["pct"].rows)
    # platform is shown on the Run tab (HPC card), the backends on the Solver card
    assert any(r.node.path == "platform" for r in window.pages["hpc"].rows)
    assert any(r.node.path == "acoustic_backend" for r in window.pages["simulation"].rows)
    assert not any(r.node.path == "platform" for r in window.pages["general"].rows)


def test_advanced_options_toggle_and_override(window):
    tab = window.tab_pages["solver"]
    window.tabs.setCurrentWidget(tab)
    grid = window.pages["grid"]
    assert not _row(grid, "grid.pml_size").block.isVisibleTo(grid)       # advanced: hidden
    assert _row(grid, "grid.resolution_mm").block.isVisibleTo(grid)      # basic: shown
    assert tab.advanced.text().startswith("Advanced (")                     # hidden count in the toggle
    tab.advanced.setChecked(True)
    assert _row(grid, "grid.pml_size").block.isVisibleTo(grid)
    tab.advanced.setChecked(False)
    window.model.set("grid.pml_size", 20)                                 # set advanced -> stays visible
    tab.refresh()
    assert _row(grid, "grid.pml_size").block.isVisibleTo(grid)


def test_tissue_properties_locked_by_default(window, monkeypatch):
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import QMessageBox
    from prestus_nexus.ui import style
    tab = window.tab_pages["medium"]
    window.tabs.setCurrentWidget(tab)
    tab.subtabs.setCurrentWidget(window.medium_subtabs["acoustic"])
    tab.advanced.setChecked(False)
    page, panel = window.pages["medium_properties"], window.medium_panel
    assert page.isVisibleTo(tab)                                          # shown without advanced
    assert not _row(page, "medium_properties.skull.sound_speed").block.isVisibleTo(page)
    assert panel.table.isVisibleTo(page)
    row, col = panel.props.index("sound_speed"), [t.name for t in panel.tissues].index("skull")
    cell = lambda: panel.table.item(row, col)
    editable = lambda: bool(cell().flags() & Qt.ItemFlag.ItemIsEditable)
    assert page.locked and not editable()
    assert cell().foreground().color().name() == style.active()["ink4"].lower()   # greyed
    assert not page.lock_button.isVisibleTo(page)                          # unlocking is advanced
    assert not window.pages["grid"].lock_button.isVisibleTo(window.pages["grid"])
    page.set_show_advanced(True)
    assert page.lock_button.isVisibleTo(page)
    monkeypatch.setattr(QMessageBox, "warning", lambda *a, **k: QMessageBox.StandardButton.Cancel)
    page.lock_button.click()                                              # declined: stays locked
    assert page.locked and not editable()
    monkeypatch.setattr(QMessageBox, "warning", lambda *a, **k: QMessageBox.StandardButton.Yes)
    page.lock_button.click()
    assert not page.locked and editable()
    cell().setText("2900")                                                # edit a cell
    assert window.model.get("medium_properties.skull.sound_speed") == 2900
    assert window.model.is_modified("medium_properties.skull.sound_speed")
    cell().setText("")                                                    # empty: back to default
    assert not window.model.is_overridden("medium_properties.skull.sound_speed")
    page.lock_button.click()                                              # relocking needs no confirmation
    assert page.locked and not editable()


def test_property_table_rows(window):
    panel = window.medium_panel
    tissues = [t.name for t in panel.tissues]
    # thermal-only rows are shaded, acoustic rows are not
    shade = lambda prop: panel.table.item(panel.props.index(prop), 0).background().style()
    assert shade("thermal_conductivity") != shade("sound_speed")
    assert shade("absorption_coeff") == shade("perfusion")
    # absorption: the field a tissue does not set is shown as implied, on the attenuation power law
    col = tissues.index("skull")
    coeff = panel.table.item(panel.props.index("absorption_coeff"), col)
    a = window.model.get("medium_properties.skull.alpha_coeff")
    f = window.model.get("medium_properties.skull.absorption_fraction")
    assert coeff.text() == f"({f * a:.4g})" and coeff.font().italic()
    window.model.set("medium_properties.skull.absorption_coeff", 3.78)
    panel.refresh()
    frac = panel.table.item(panel.props.index("absorption_fraction"), col)
    assert frac.text() == f"({3.78 / a:.4g})"


def test_medium_selector_controls_subtabs(window):
    tab, panel = window.tab_pages["medium"], window.medium_panel
    window.tabs.setCurrentWidget(tab)
    vis = lambda name: tab.subtabs.isTabVisible(tab.subtabs.indexOf(window.medium_subtabs[name]))
    assert vis("head") and vis("acoustic") and vis("uncertainty") and vis("pct")   # layered
    window.medium_selector.combo.setCurrentIndex(window.medium_selector.combo.findData("water"))
    window.medium_selector._picked(window.medium_selector.combo.currentIndex())
    assert window.model.get("simulation.medium") == "water"
    assert not vis("head") and vis("acoustic") and not vis("uncertainty") and not vis("pct")
    assert tab.subtabs.currentWidget() is window.medium_subtabs["acoustic"]
    shown = [t.name for i, t in enumerate(panel.tissues) if not panel.table.isColumnHidden(i)]
    assert shown == ["water"]
    window.model.set("simulation.medium", "phantom")
    tab.refresh()
    assert vis("head") and not vis("pct")
    assert all(not panel.table.isColumnHidden(i) for i in range(len(panel.tissues)))
    assert "Phantom" in panel.note.text()
    window.goto("simulation.medium")
    assert window.tabs.currentWidget() is tab


def test_goto_opens_pct_subtab(window):
    window.model.set("pct.enabled", 1)
    window.goto("pct.mapping_density")
    tab = window.tab_pages["medium"]
    assert tab.subtabs.currentWidget() is window.medium_subtabs["pct"]
    assert "Density" in window.pct_formulas.view.toPlainText()


def test_dependent_fields(window):
    m = window.model
    assert not m.is_visible("path.localite_post")                         # placement.mode manual
    m.set("placement.mode", "localite")
    assert m.is_visible("path.localite_post") and not m.is_visible("path.localite_export")
    m.set("placement.mode", "heuristic")
    assert m.is_visible("path.localite_export") and not m.is_visible("path.localite_raw")
    assert not m.is_visible("pct.mapping_density")                        # pct.enabled 0
    m.set("pct.enabled", 1)
    assert m.is_visible("pct.mapping_density") and not m.is_visible("pct.custom_slope")
    m.set("pct.skull_mapping", "custom")
    assert m.is_visible("pct.custom_slope")
    assert not m.is_visible("thermal.log_thermal_updates") and not m.is_visible("thermal_backend")
    m.set("modules.run_heating_sims", 1)
    assert m.is_visible("thermal.log_thermal_updates") and m.is_visible("thermal_backend")
    assert not m.is_visible("time_reversal.placement_search.max_angle_deg")
    m.set("time_reversal.enabled", True)
    assert m.is_visible("time_reversal.placement_search.enabled")
    assert not m.is_visible("time_reversal.placement_search.max_angle_deg")
    m.set("time_reversal.placement_search.enabled", True)
    assert m.is_visible("time_reversal.placement_search.max_angle_deg")
    assert m.is_visible("time_reversal.pose_workflow.enabled")


def test_uncertainty_card_lists_variants(window):
    card = window.uncertainty_card
    window.tabs.setCurrentWidget(window.tab_pages["medium"])
    window.tab_pages["medium"].subtabs.setCurrentWidget(window.medium_subtabs["uncertainty"])
    assert not card.isVisibleTo(window)                                   # uncertainty off
    window.model.set("simulation.uncertainty", True)
    window.tab_pages["medium"].refresh()
    assert card.isVisibleTo(window)
    assert not card.sets.errors
    n_tissues = len(card.sets.sets["default"])
    assert card.table.rowCount() == n_tissues * len(PROPERTIES)
    plots = card.canvases                                                 # the bounds plot: Uncertainty sub-tab
    if plots:
        assert len(plots["brain"].figure.axes) == 7                       # one column per plotted property
        assert len(plots["skull"].figure.axes) == 7


def test_card_advanced_button_opens_one_card(window):
    tab = window.tab_pages["solver"]
    window.tabs.setCurrentWidget(tab)
    grid, timing = window.pages["grid"], window.pages["timing"]
    assert grid.adv_button.isVisibleTo(grid) and grid.adv_button.text().isdigit()   # icon + hidden count
    grid.adv_button.click()
    assert grid.show_advanced and not timing.show_advanced and not tab.advanced.isChecked()
    assert _row(grid, "grid.pml_size").block.isVisibleTo(grid)
    assert not _row(timing, "timing.pt_timestep").block.isVisibleTo(timing)
    assert grid.adv_button.text() == ""                                  # shown: icon only
    grid.adv_button.click()
    assert not _row(grid, "grid.pml_size").block.isVisibleTo(grid)


def test_tune_for_is_hidden(window):
    tab = window.tab_pages["solver"]
    window.tabs.setCurrentWidget(tab)
    tab.advanced.setChecked(True)
    grid = window.pages["grid"]
    assert not tab.tune.isVisibleTo(tab) and not tab.tune_label.isVisibleTo(tab)
    assert "speed" not in _row(grid, "grid.resolution_mm").label.text()   # no goal tags

def test_goto_reveals_advanced_field(window):
    window.goto("headmodel.smooth_fwhm_mm")
    tab = window.tab_pages["medium"]
    assert window.tabs.currentWidget() is tab
    page = window.pages["head_smoothing"]                       # the Smoothing card
    assert page.show_advanced and _row(page, "headmodel.smooth_fwhm_mm").block.isVisibleTo(page)


def test_medium_cards_regrouped(window):
    rows = lambda sec: [r.node.path for r in window.pages[sec].rows]
    assert rows("head_smoothing")[0] == "headmodel.smooth_method"
    assert "headmodel.smooth_properties" in rows("head_smoothing")
    assert rows("grid_refinement")[0] == "headmodel.csf_expansion"
    assert not any(p.startswith("headmodel.smooth") or p == "headmodel.csf_expansion" for p in rows("headmodel"))
    layers = window.pages["layers"]
    assert layers.lockable and layers.locked
    assert window.model.schema.section_is_advanced("layers")
    assert "charm" in layers.lock_note and layers.icon_name == "layers"
    tab = window.tab_pages["medium"]
    window.tabs.setCurrentWidget(tab)
    window.model.set("simulation.medium", "phantom")
    tab.refresh()
    assert layers.isVisibleTo(tab)                     # a phantom is defined by its layers


def test_labels_and_descriptions(window):
    row = _row(window.pages["timing"], "timing.pri")
    assert "Pulse repetition interval" in row.label.text()
    assert "timing.pri" in row.label.toolTip()


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
def test_coupling_on_sonication_tab(window):
    tab = window.tab_pages["sonication"]
    window.tabs.setCurrentWidget(tab)
    page = window.pages["coupling"]
    assert window._section_tab["coupling"] == "sonication"
    assert any(r.node.path == "simulation.transducer_coupling" for r in page.rows)
    assert not any(r.node.path == "simulation.transducer_coupling" for r in window.pages["simulation"].rows)
    assert any(r.node.path == "simulation.debug" for r in window.pages["modules"].rows)   # Run tab
    tab.refresh()
    assert not page.isVisibleTo(tab)                           # one transducer: nothing to couple
    window.pages["transducer"]._add_transducer()
    window._changed()
    tab.refresh()
    assert page.isVisibleTo(tab) and window.coupling_view.isVisibleTo(page)
    assert window.coupling_link.isVisibleTo(window.pages["transducer"]) or \
        not window.pages["transducer"].isVisible()
    window.coupling_view.toggle.setChecked(True)
    assert "async" in window.coupling_view.text.text()
    window.goto("simulation.transducer_coupling")
    assert window.tabs.currentWidget() is tab


def test_run_tab_subtabs(window):
    from prestus_nexus.ui.main_window import REVIEW, RUN
    tab = window.tab_pages["run"]
    assert [tab.subtabs.tabText(i) for i in range(tab.subtabs.count())] == \
        ["Stages", "Where to run", "Plan", "Review YAML", "Launch", "Jobs"]
    assert tab.subtab_of(window.pages["modules"]) is tab.subtabs.widget(0)
    assert tab.subtab_of(window.pages["hpc"]) is tab.subtabs.widget(1)
    assert window.show_special(RUN) is tab.subtabs.currentWidget()
    assert window.show_special(REVIEW) is tab.subtabs.currentWidget()
    window.goto("hpc.timelimit")
    assert tab.subtabs.currentWidget() is tab.subtabs.widget(1)


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


def test_uncertainty_switch_turns_on(window):
    """Ticking the switch on the Uncertainty analysis card: the plot appears, the variant table and
    files sit in the card's closed advanced fold,
    and the window refresh (I/O preview of a non-single mode) completes."""
    window.model.set("simulation.medium", "layered")
    window.goto("simulation.uncertainty")
    page = window.pages["uncertainty"]
    row = next(r for r in page.rows if r.node.path == "simulation.uncertainty")
    row.editor.box.click()
    assert window.model.get("simulation.uncertainty") is True
    assert sum(r.block.isVisibleTo(page) for r in page.rows) == 1      # the switch only
    card = window.uncertainty_card
    files = [r for r in page.rows if r.node.path in FILE_KEYS]
    assert len(files) == 4 and not card.advanced.isVisibleTo(card)        # files and table: folded
    assert all(r.block.parentWidget() is card.advanced for r in files)
    tab = window.tab_pages["medium"]
    assert tab.subtabs.currentWidget() is window.medium_subtabs["uncertainty"]   # its own sub-tab
    assert not card.isHidden()
    window.goto("simulation.uncertainty_temp_high")                      # goto opens the fold
    assert card.advanced.isVisibleTo(card)
    assert all(r.block.isVisibleTo(card) for r in files)
    row.editor.box.click()
    assert window.model.get("simulation.uncertainty") is False
    assert window.uncertainty_card.isHidden()
