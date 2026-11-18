"""Main-window wiring, offscreen: menus, undo, the transducers overview, guarded saves."""

import os
import stat

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")

from PyQt6.QtWidgets import QApplication  # noqa: E402

from prestus_nexus.model import ConfigModel  # noqa: E402
from prestus_nexus.session import Session  # noqa: E402
from prestus_nexus.ui import errors  # noqa: E402
from prestus_nexus.ui.main_window import MainWindow  # noqa: E402


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def window(app, schema, prestus_root, tmp_path):
    w = MainWindow(ConfigModel(schema), Session(tmp_path / "session.yaml"), prestus_root)
    yield w
    w.dirty = False
    w.close()


def _menu(window, title):
    for action in window.menuBar().actions():
        if action.text().replace("&", "") == title:
            return action.menu()
    raise AssertionError(f"no {title} menu")


def test_menus_and_names(window):
    titles = [a.text().replace("&", "") for a in window.menuBar().actions()]
    assert titles[:2] == ["File", "Edit"] and "Help" in titles
    tools = [a.text() for a in _menu(window, "Tools").actions()]
    for name in ("Localite", "PlanTUS", "TUS entry", "Calibration", "Usage statistics…"):
        assert name in tools
    shortcuts = {a.shortcut().toString() for a in _menu(window, "Help").actions()}
    assert "F1" in shortcuts


def test_undo_reverts_a_change(window):
    assert not window.undo_action.isEnabled()
    window.model.set("subject_id", 7)
    assert window.undo_action.isEnabled()
    window.undo_action.trigger()
    assert window.model.get("subject_id") != 7
    assert window.redo_action.isEnabled()
    window.redo_action.trigger()
    assert window.model.get("subject_id") == 7


def test_undo_covers_the_study_table(window):
    """Study-tab table edits and config edits share one linear history."""
    from prestus_nexus.ui.main_window import STUDY
    page = window.show_special(STUDY)
    page.add_subjects("1-3")
    window.model.set("subject_id", 7)
    page.add_column("condition_b")
    assert window.undo_action.text() == "Undo add column"
    window.undo_action.trigger()
    assert "condition_b" not in page.manifest.columns
    assert page.table.columnCount() == len(page.manifest.columns)     # the table was refilled
    window.undo_action.trigger()
    assert window.model.get("subject_id") != 7 and len(page.manifest.rows) == 3
    window.undo_action.trigger()
    assert not page.manifest.rows and page.table.rowCount() == 0
    window.redo_action.trigger()
    assert [r["subject_id"] for r in page.manifest.rows] == ["1", "2", "3"]


def test_transducers_overview_on_placement_settings(window):
    """The overview sits on the Placement card (no sub-tab of its own), shown with 2+
    transducers; its table and the card's transducer selector pick the same transducer."""
    tab = window.tab_pages["placement"]
    place, overview = window.pages["placement"], window.transducers_page
    assert tab.subtabs.indexOf(overview) < 0                  # no Transducers sub-tab
    assert place.extras.indexOf(overview) >= 0
    tab.refresh()
    assert overview.isHidden() and overview.details.isHidden() and not place.tx_bar.isHidden()
    overview.add_btn.click()                                   # add from the overview
    tab.refresh()
    assert window.model.n_transducers == 2
    assert not overview.isHidden() and not overview.details.isHidden() and place.tx_bar.isHidden()
    overview.table.selectRow(0)
    overview.table.selectRow(1)
    assert place.ti == 1                                       # the table selects what is edited
    place.tx_combo.setCurrentIndex(0)
    assert overview.table.currentRow() == 0                    # and follows the selector
    overview.table.selectRow(1)
    overview.remove_btn.click()
    tab.refresh()
    assert window.model.n_transducers == 1 and overview.isHidden()


@pytest.mark.skipif(os.name == "nt", reason="POSIX permissions")
def test_failed_save_keeps_file_and_window(window, tmp_path, monkeypatch):
    shown = []
    monkeypatch.setattr(errors, "show_error", lambda *a, **k: shown.append(a))
    folder = tmp_path / "ro"
    folder.mkdir()
    study = folder / "config_x.yaml"
    study.write_text("subject_id: 1\n")
    folder.chmod(stat.S_IRUSR | stat.S_IXUSR)       # no new files: the atomic write must fail
    try:
        window.model.path = study
        window.model.set("subject_id", 3)
        window.dirty = True
        assert window.save() is False
        assert shown and window.dirty
        assert study.read_text() == "subject_id: 1\n"
    finally:
        folder.chmod(stat.S_IRWXU)



def test_one_current_transducer_across_selectors(window):
    """Picking a transducer on any selector (overview table, Placement card, Transducer card,
    Viewer) makes it the current one everywhere; adding or removing moves it too."""
    from prestus_nexus.ui.main_window import VIEWER
    m = window.model
    for _ in range(2):
        m.add_transducer()
    window.tab_pages["placement"].refresh()
    place, tx, overview = window.pages["placement"], window.pages["transducer"], window.transducers_page
    viewer = window.pages.get(VIEWER)
    selectors = {
        "table": lambda ti: overview.table.selectRow(ti),
        "placement card": lambda ti: place.tx_combo.setCurrentIndex(ti),
        "transducer card": lambda ti: tx.tx_combo.setCurrentIndex(ti),
    }
    if viewer is not None:
        selectors["viewer"] = lambda ti: viewer.tx_combo.setCurrentIndex(ti)
    for name, pick in selectors.items():
        for ti in (2, 0, 1):
            pick(ti)
            assert place.ti == tx.ti == ti, name
            assert overview.table.currentRow() == ti, name
            assert f"Transducer {ti + 1}" in overview.editing.text(), name
            if viewer is not None:
                assert viewer.ti == ti, name
    tx.tx_combo.setCurrentIndex(1)
    window.pages["tx_drive"].ti == 1 or window.pages["tx_drive"].refresh()
    assert window.pages["tx_drive"].ti == 1             # follows the Transducer card
    overview.add_btn.click()
    assert place.ti == tx.ti == 3 and m.n_transducers == 4


def test_card_subtitles_are_on_hover(window):
    """Every card's subtitle (overlay 'sub') is the hover text of its title, not a line of text."""
    from PyQt6.QtWidgets import QLabel
    from prestus_nexus.ui.section import SectionPage
    pages = [p for p in window.pages.values() if isinstance(p, SectionPage)]
    with_sub = [p for p in pages if window.model.schema.section_subtitle(p.section)]
    assert with_sub
    for p in with_sub:
        sub = window.model.schema.section_subtitle(p.section)
        assert p.title.toolTip() == sub, p.section
        shown = [w.text() for w in p.findChildren(QLabel) if w.objectName() == "CardSub"]
        assert sub not in shown, p.section


def test_units_follow_the_labels(window):
    """Every field with a unit shows it in its label, 'Name (unit)', and no unit column."""
    from PyQt6.QtWidgets import QLabel
    from prestus_nexus.ui.section import SectionPage
    seen = 0
    for page in window.pages.values():
        if not isinstance(page, SectionPage):
            continue
        for row in page.rows:
            if row.node.unit:
                seen += 1
                assert f"({row.node.unit})" in row.label.text().replace("&lt;", "<"), row.node.path
                beside = [w.text() for w in row.container.findChildren(QLabel)]
                assert row.node.unit not in beside, row.node.path
    assert seen > 10


def test_sonication_without_thermal_stage(window):
    """Thermal stage off: the thermal model and sequence cards are hidden and one line says
    why, with a link that turns heating on; the pulse protocol stays."""
    window.goto("timing.pd")
    tab = window.tab_pages["sonication"]
    window.model.set("modules.run_heating_sims", 0)
    tab.refresh()
    assert not window.heating_off_note.isHidden()
    assert window.pages["thermal"].isHidden() and window.sequence_card.isHidden()
    assert not window.timing_panel.isHidden()
    window.heating_off_note.linkActivated.emit("#on")
    tab.refresh()
    assert window.model.get("modules.run_heating_sims") == 1
    assert window.heating_off_note.isHidden()
    assert not window.pages["thermal"].isHidden() and not window.sequence_card.isHidden()


def test_pipeline_stages_have_icons(window):
    """Every stage row on Run > Stages shows its icon between the checkbox and the label."""
    from PyQt6.QtWidgets import QLabel
    from prestus_nexus.ui.section import STAGE_ICONS
    from prestus_nexus.ui import icons
    page = window.pages["modules"]
    assert set(STAGE_ICONS.values()) <= set(icons.names())
    rows = [r for r in page.rows if r.node.path in STAGE_ICONS]
    assert rows
    for r in rows:
        pics = [w for w in r.block.findChildren(QLabel) if w.pixmap() is not None and not w.pixmap().isNull()]
        assert pics, r.node.path


def test_placement_shows_only_the_active_method_page(window):
    """Of Localite, Heuristic, PlanTUS, and TUS entry only the method in use has a sub-tab;
    one opened from the Tools menu shows while it is open; Viewer always shows."""
    from prestus_nexus.ui.main_window import HEURISTIC, LOCALITE, PLANTUS, TUS_ENTRY, VIEWER
    tab = window.tab_pages["placement"]
    sub = tab.subtabs
    shown = lambda key: sub.isTabVisible(sub.indexOf(window.pages[key]))   # noqa: E731
    methods = {"localite": LOCALITE, "heuristic": HEURISTIC, "plantus": PLANTUS, "tus_entry": TUS_ENTRY}
    for mode, key in methods.items():
        window.model.set("placement.mode", mode)
        sub.setCurrentWidget(window.pages[VIEWER])
        tab.refresh()
        assert shown(key) and shown(VIEWER), mode
        assert not any(shown(k) for k in methods.values() if k != key), mode
    window.model.set("placement.mode", "manual")
    sub.setCurrentWidget(window.pages[VIEWER])
    tab.refresh()
    assert not any(shown(k) for k in methods.values())
    window.show_special(PLANTUS)                  # Tools menu: shown while open
    tab.refresh()
    assert shown(PLANTUS)
    sub.setCurrentWidget(window.pages[VIEWER])
    tab.refresh()
    assert not shown(PLANTUS)


def test_review_yaml_marks_values_that_differ(window):
    """Review YAML: each 'key: value' line maps to its key (transducer list items by index);
    lines that differ from the default are tinted; one equal to the default is not."""
    from prestus_nexus.ui.main_window import ReviewPage, yaml_line_keys
    text = "subject_id: 3\ntiming:\n  pd: 0.02\ntransducer:\n  - serial: A\n    combo:\n      ds_serial: B\n  - serial: C\n"
    keys = [(p, ti) for _n, p, ti in yaml_line_keys(text)]
    assert keys == [("subject_id", 0), ("timing.pd", 0), ("transducer.serial", 0),
                    ("transducer.combo.ds_serial", 0), ("transducer.serial", 1)]
    m = window.model
    m.set("subject_id", 3)
    m.overrides.setdefault("grid", {})["resolution_mm"] = 0.5      # equal to the default
    page = ReviewPage(window)
    page.refresh()
    by_path = {info[0]: info[2] for info in page._lines.values()}
    assert by_path["subject_id"] is True and by_path["grid.resolution_mm"] is False
    assert len(page.text.extraSelections()) == sum(by_path.values())
