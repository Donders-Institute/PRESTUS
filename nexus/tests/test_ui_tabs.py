"""Tabbed main window: tabs, advanced options, navigation, timing panel, Localite import."""

import json
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")

from PyQt6.QtWidgets import QApplication  # noqa: E402

from prestus_nexus.model import ConfigModel  # noqa: E402
from prestus_nexus.session import Session  # noqa: E402
from prestus_nexus.uncertainty import OVERRIDES_KEY, PROPERTIES  # noqa: E402
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
    sel = window.placement_selector                              # the Placement tab's header row
    assert not any(r.node.path == "placement.mode" for r in window.pages["placement"].rows)
    sel.refresh()
    combo = sel.combo
    items = {combo.itemData(i): combo.model().item(i) for i in range(combo.count())}
    assert items["manual"].isEnabled()
    assert not items["heuristic"].isEnabled()                 # no m2m, no SimNIBS
    assert "unavailable" in items["heuristic"].text()
    assert "SimNIBS bin folder" in items["plantus"].toolTip()
    window.model.set("placement.mode", "heuristic")             # current value stays selectable
    sel.refresh()
    assert combo.model().item(combo.findData("heuristic")).isEnabled()
    assert "Cannot run yet" in sel.note.text()


def test_plantus_actions_run_placement_only(window, monkeypatch):
    from prestus_nexus.ui.main_window import RUN, VIEWER, PlantusActions
    from prestus_nexus.ui.plantus_page import TASK
    from prestus_nexus.ui.run_page import RunPage
    ran = []
    monkeypatch.setattr(RunPage, "run_script", lambda self, task, script: ran.append((task, script)) or True)
    acts = window.pages["placement"].findChild(PlantusActions)
    window.model.set("placement.mode", "manual")
    acts.refresh()
    assert acts.isHidden()
    window.model.set("placement.mode", "plantus")
    window.model.set("placement.plantus.vertex", 12)
    acts.refresh()
    assert not acts.isHidden()
    acts._run(True)
    task, script = ran[0]
    assert task == TASK and script.name == "plantus_t.m"
    assert "'placement.plantus.interactive', 1;" in script.read_text()
    assert window.model.get("placement.plantus.vertex") == []      # maps recomputed: choice cleared
    shown = []
    monkeypatch.setattr(type(window), "show_special", lambda self, key: shown.append(key))
    window.pages[RUN].finished.emit(TASK, True)
    assert shown == [VIEWER] and not acts._waiting


def test_plantus_page_places_chosen_vertex(window, monkeypatch, tmp_path):
    from prestus_nexus.ui.main_window import PLANTUS, RUN
    from prestus_nexus.ui.plantus_page import TASK
    from prestus_nexus.ui.run_page import RunPage
    import test_plantus
    ran = []
    monkeypatch.setattr(RunPage, "run_script", lambda self, task, script: ran.append((task, script)) or True)
    page = window.pages[PLANTUS]
    window.model.set("path.sim", str(tmp_path / "sim"))
    window.model.set("placement.plantus.target_name", "M1")
    window.model.set("placement.mode", "plantus")
    page.refresh()
    assert page.stack.is_empty()                                    # no maps yet: Compute maps offered
    path = test_plantus.write_metrics(window.model)
    page.refresh()
    assert page.maps is not None and page.maps.path == path
    assert page.clicked == page.maps.selected
    page.select(5)
    page.place()
    task, script = ran[0]
    assert task == TASK and "'placement.plantus.interactive', 0;" in script.read_text()
    assert window.model.get("placement.plantus.vertex") == 5
    window.pages[RUN].finished.emit(TASK, True)
    assert not page._waiting


def test_plantus_page_submits_job_when_run_tab_submits(window, monkeypatch, tmp_path):
    from prestus_nexus.ui.main_window import PLANTUS, RUN
    from prestus_nexus.ui.run_page import RunPage
    sent = []
    monkeypatch.setattr(RunPage, "submits_to_cluster", lambda self: True)
    monkeypatch.setattr(RunPage, "submit_partial", lambda self, name, what, overrides, extra="": sent.append(name) or True)
    monkeypatch.setattr(RunPage, "run_script", lambda self, *a: pytest.fail("ran locally"))
    window.model.set("placement.mode", "plantus")
    page = window.pages[PLANTUS]
    page.compute()
    assert sent == ["plantus"] and page._watching and not page._waiting
    window.pages[RUN].jobs_changed.emit()                          # job update: reloads quietly


def test_transducer_tab_cards(window):
    from prestus_nexus.model import TX
    tab = window.tab_pages["transducer"]
    window.tabs.setCurrentWidget(tab)
    tx, geo, drive = (window.pages[s] for s in ("transducer", "tx_geometry", "tx_drive"))
    assert [p.section for p in tab.sections[:3]] == ["transducer", "tx_geometry", "tx_drive"]
    # one transducer selector: the Geometry and Drive cards follow the Transducer card
    assert tx.tx_bar.isVisibleTo(tx) and not geo.tx_bar.isVisibleTo(geo) and not drive.tx_bar.isVisibleTo(drive)
    assert any(r.node.path == f"{TX}.target_isppa_wcm2" for r in drive.rows)
    assert any(r.node.path == f"{TX}.annular.elem_phase_deg" for r in drive.rows)
    assert not any(r.node.path == f"{TX}.annular.elem_phase_deg" for r in geo.rows)
    assert any(r.node.path == f"{TX}.annular.elem_od_mm" for r in geo.rows)
    # the calibration coverage plot is on the Calibration sub-tab only
    assert not hasattr(window, "cal_cov")
    tx._add_transducer()                                                  # selects the new transducer
    tab.refresh()
    assert geo.ti == 1 and drive.ti == 1
    tx.tx_combo.setCurrentIndex(0)
    tab.refresh()
    # geometry from the library: the card stays listed, collapsed, until opened by hand
    window.model.set(f"{TX}.type", "annular")
    window._changed()
    assert not geo.body.isHidden() and geo.isVisibleTo(tab)
    if window.model.library is not None and window.model.library.lookup("CTX_500_026") is not None:
        window.model.set(f"{TX}.serial", "CTX_500_026")
        window._changed()
        assert geo.body.isHidden() and geo.isVisibleTo(tab) and "library" in geo.count.text()
        geo.toggle.click()
        window._changed()
        assert not geo.body.isHidden()


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
        ["Stages", "Where to run", "Plan", "Review YAML", "Launch", "Jobs", "Live view"]
    assert not tab.subtabs.isTabVisible(6)          # Live view: only once a run has frames
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
    """Ticking the switch on the Uncertainty analysis card: the plot appears, the variant table is
    folded and locked; unlocked, an edit is bold red and warned,
    and the window refresh (I/O preview of a non-single mode) completes."""
    window.model.set("simulation.medium", "layered")
    window.goto("simulation.uncertainty")
    page = window.pages["uncertainty"]
    row = next(r for r in page.rows if r.node.path == "simulation.uncertainty")
    row.editor.box.click()
    assert window.model.get("simulation.uncertainty") is True
    assert [r.node.path for r in page.rows if r.block.isVisibleTo(page)] == \
        ["simulation.uncertainty", "simulation.uncertainty_fov"]       # the switch and the variants' grid
    from PyQt6.QtCore import Qt
    from prestus_nexus.ui import style
    card = window.uncertainty_card
    assert not card.values.isVisibleTo(card)                              # table: folded, plot only
    hidden = next(r for r in page.rows if r.node.path == OVERRIDES_KEY)
    assert not hidden.block.isVisibleTo(window)                           # edited in the table only
    tab = window.tab_pages["medium"]
    assert tab.subtabs.currentWidget() is window.medium_subtabs["uncertainty"]   # its own sub-tab
    assert not card.isHidden()
    window.goto(OVERRIDES_KEY)                                            # goto opens the fold
    assert card.values.isVisibleTo(card)
    # locked like the property table: greyed, nothing editable
    cell = _uncertainty_cell(card, "Brain", "Sound speed", "Temp-H")
    assert card.locked and not cell.flags() & Qt.ItemFlag.ItemIsEditable
    assert not card.warning.isVisibleTo(card)
    card.set_locked(False)
    cell = _uncertainty_cell(card, "Brain", "Sound speed", "Temp-H")
    assert cell.flags() & Qt.ItemFlag.ItemIsEditable
    bound = card.sets.bound("temp_high", "brain", "sound_speed")
    cell.setText("1600")                                                  # an edit: bold red, warned
    assert window.model.get(OVERRIDES_KEY) == {"temp_high": {"brain": {"sound_speed": 1600}}}
    cell = _uncertainty_cell(card, "Brain", "Sound speed", "Temp-H")
    assert cell.font().bold() and cell.foreground().color().name() == style.active()["red"].lower()
    assert card.warning.isVisibleTo(card) and "1 bound value edited" in card.warning.text()
    assert not _uncertainty_cell(card, "Brain", "Perfusion", "Trans-H").flags() & Qt.ItemFlag.ItemIsEditable
    cell.setText(f"{bound}")                                              # the bound again: edit removed
    assert window.model.get(OVERRIDES_KEY) in (None, [], {})
    assert not card.warning.isVisibleTo(card)
    row.editor.box.click()
    assert window.model.get("simulation.uncertainty") is False
    assert window.uncertainty_card.isHidden()


def _uncertainty_cell(card, tissue, prop, variant):
    t = card.table
    col = [t.horizontalHeaderItem(c).text() for c in range(t.columnCount())].index(variant)
    row = next(r for r in range(t.rowCount())
               if t.item(r, 0).text() == tissue and t.item(r, 1).text() == prop)
    return t.item(row, col)


def _segmentation(tmp_path, labels=(1, 2, 3, 5, 7, 8), report=True):
    """A small m2m folder: final_tissues.nii.gz with these labels, and charm_report.html."""
    import nibabel as nib
    import numpy as np
    m2m = tmp_path / "seg" / "m2m_sub-001"
    m2m.mkdir(parents=True)
    seg = np.zeros((12, 12, 12), dtype=np.uint8)
    for i, label in enumerate(labels):
        seg[i] = label
    nib.save(nib.Nifti1Image(seg, np.eye(4)), str(m2m / "final_tissues.nii.gz"))
    if report:
        (m2m / "charm_report.html").write_text("<html><body>CHARM report</body></html>")
    return m2m


def test_charm_report_card(window, tmp_path, monkeypatch):
    pytest.importorskip("nibabel")
    from PyQt6.QtGui import QDesktopServices
    tab, card = window.tab_pages["medium"], window.charm_report_card
    seg = window.pages["segmentation"]
    window.model.set("subject_id", 1)
    window.model.set("path.seg", str(tmp_path / "seg"))
    tab.refresh()
    assert card.isHidden()                            # no segmentation report yet
    m2m = _segmentation(tmp_path)
    tab.refresh()
    assert not card.isHidden() and card.report == m2m / "charm_report.html"
    assert card.parentWidget() is seg.parentWidget()  # its own card, beside the segmentation card
    column = seg.parentWidget().layout()
    assert column.indexOf(card) > column.indexOf(seg)
    assert not card.button.isChecked() and card.body.isHidden()     # closed until opened
    opened = []
    monkeypatch.setattr(QDesktopServices, "openUrl", lambda url: opened.append(url.toLocalFile()))
    window.simnibs_panel.refresh()
    window.simnibs_panel.report.click()               # 'Open charm report': the browser
    assert opened == [(m2m / "charm_report.html").as_posix()]
    window.model.set("simulation.medium", "water")
    tab.refresh()
    assert card.isHidden()


def test_charm_report_css_follows_palette():
    from prestus_nexus.ui import style
    from prestus_nexus.ui.charm_report import STYLE_ID, report_css, style_js
    for palette in (style.C, style.C_DARK):
        css = report_css(palette)
        assert palette["surface"] in css and palette["ink"] in css and palette["accent"] in css
        assert "font-size: 13px" in css and "font-size: 14px" in css      # card text, card titles
    assert STYLE_ID in style_js(style.C) and "textContent" in style_js(style.C)


def test_property_table_greys_tissues_without_voxels(window, tmp_path):
    pytest.importorskip("nibabel")
    from PyQt6.QtCore import Qt
    from prestus_nexus.ui import style
    _segmentation(tmp_path)                           # charm labels: no combined skull (4)
    window.model.set("subject_id", 1)
    window.model.set("path.seg", str(tmp_path / "seg"))
    panel = window.medium_panel
    panel.refresh()
    panel.wait_census()
    names = [t.name for t in panel.tissues]
    head = lambda t: panel.table.horizontalHeaderItem(names.index(t)).text()
    ink4 = style.active()["ink4"].lower()
    cell = lambda t: panel.table.item(panel.props.index("sound_speed"), names.index(t))
    assert panel.census_state == "done"
    assert not panel.census["skull"].modelled and panel.census["skull_cortical"].modelled
    tip = lambda t: panel.table.horizontalHeaderItem(names.index(t)).toolTip()
    # the volume is in the header's tooltip only; the header names why a tissue is not modelled
    assert "no voxels" in head("skull") and head("skin") == "Skin" and "mL" in tip("skin")
    assert "Skull (single layer)" in panel.note.text()
    panel.set_locked(False)
    assert cell("skull").foreground().color().name() == ink4          # greyed though unlocked
    assert cell("skin").foreground().color().name() != ink4
    assert cell("skull").flags() & Qt.ItemFlag.ItemIsEditable          # still editable
    window.model.set("pct.enabled", 1)                # pCT: one skull layer, split ones merged
    panel.refresh()
    assert panel.census["skull"].modelled
    assert "pCT" in head("skull_cortical") and cell("skull_cortical").foreground().color().name() == ink4


def test_label_column_and_tabs_follow_the_window(window):
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import QGridLayout
    from prestus_nexus.ui.section import LABEL_MIN_WIDTH, LABEL_WIDTH, label_width
    assert label_width(400) == LABEL_MIN_WIDTH and label_width(2000) == LABEL_WIDTH
    assert LABEL_MIN_WIDTH < label_width(800) < LABEL_WIDTH
    page = window.pages["grid"]
    grids = [g for g in page.findChildren(QGridLayout) if g.property("label_inset") is not None]
    assert grids
    page.resize(600, page.height())
    page._fit_label_column()
    assert all(g.columnMinimumWidth(0) == label_width(600) - g.property("label_inset") for g in grids)
    assert all(r.block.maximumWidth() == label_width(600) + 40 for r in page.rows)
    page.resize(1200, page.height())
    page._fit_label_column()
    assert all(g.columnMinimumWidth(0) == LABEL_WIDTH - g.property("label_inset") for g in grids)
    for bar in (window.tabs.tabBar(), window.tab_pages["medium"].subtabs.tabBar()):
        assert bar.expanding() and bar.elideMode() == Qt.TextElideMode.ElideNone
        assert not bar.usesScrollButtons()                                # full labels, none hidden


def test_narrow_tabs_show_icons_except_quickstart(app, window):
    from prestus_nexus.ui import icons
    from prestus_nexus.ui.tabbar import tab_widget
    from PyQt6.QtWidgets import QWidget
    w = tab_widget()
    bar = w.tabBar()
    bar.set_compactable()
    w.addTab(QWidget(), "QuickStart")
    bar.set_tab(0, "QuickStart")                                    # no short icon: keeps its label
    for i, (label, icon) in enumerate((("Transducer", "tx_annular"), ("Placement", "target"),
                                       ("Medium  2", "layers"), ("Sonication", "pulse")), start=1):
        w.addTab(QWidget(), label)
        bar.set_tab(i, label, tooltip="No problems found", short_icon=icons.icon(icon),
                    short_text="2" if label.startswith("Medium") else "")
    w.resize(2000, 300)
    w.show()
    app.processEvents()
    assert not bar.is_compact() and bar.tabText(1) == "Transducer"
    assert bar.minimumSizeHint().width() < bar.full_labels_width()   # the window may get narrower
    w.resize(bar.minimumSizeHint().width() + 10, 300)
    app.processEvents()
    assert bar.is_compact()
    assert bar.tabText(0) == "QuickStart" and bar.tabText(1) == "" and bar.tabText(3) == "2"
    assert not bar.tabIcon(1).isNull() and "Transducer" in bar.tabToolTip(1)
    w.resize(2000, 300)
    app.processEvents()
    assert not bar.is_compact() and bar.tabText(1) == "Transducer"
    assert bar.tabToolTip(1) == "No problems found"
    w.close()
    main = window.tabs.tabBar()                                      # the main tabs: compactable
    assert main.minimumSizeHint().width() < main.full_labels_width()


# ── Find setting (ui/field_search.py) ─────────────────────────────────

def test_find_setting_matches_label_key_and_help(window):
    from prestus_nexus.ui.field_search import entries, search
    items = entries(window.model.schema, window.display_section, window.can_goto)
    assert items and all(window.can_goto(e.node.path) for e in items)
    by_label = search(items, "k-wave backend")
    assert by_label[0].node.path == "simulation.code_type"
    assert search(items, "code_type")[0].node.path == "simulation.code_type"      # the key
    assert any(e.node.path == "simulation.code_type" for e in search(items, "apple silicon"))  # help/choice label
    assert search(items, "") == [] and search(items, "zzqqxx") == []


def test_find_setting_opens_the_field(window):
    from PyQt6.QtCore import Qt
    window.find_setting()
    dlg = window._field_search
    dlg.edit.setText("water_baseline_grid")
    item = dlg.list.currentItem()
    assert item is not None and item.data(Qt.ItemDataRole.UserRole) == "grid.water_baseline_grid"
    dlg._open(item)
    assert not dlg.isVisible()
    assert window.tabs.currentWidget().key == window.model.schema.tab_of(window.display_section("grid.water_baseline_grid"))


def test_settings_subtab_only_while_it_shows_something(window):
    """Every tab: the 'Settings' sub-tab is shown exactly when it has a visible section or a card;
    Results has none, so it opens on its results."""
    checked = 0
    for key, tab in window.tab_pages.items():
        if tab.subtabs is None or tab.settings_page is None:
            continue
        tab.refresh()
        i = tab.subtabs.indexOf(tab.settings_page)
        has = bool(tab.cards) or any(not p.isHidden() for p in tab.sections)
        assert tab.subtabs.isTabVisible(i) == has, key
        if not has:
            assert tab.subtabs.currentIndex() != i, key
        checked += 1
    assert checked
    results = window.tab_pages["results"]
    assert not results.subtabs.isTabVisible(results.subtabs.indexOf(results.settings_page))
