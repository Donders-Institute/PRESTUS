"""Main window layout: workflow tabs spanning the window, Placement opening on the Viewer, positions
edited per transducer on Placement, the serial edited only in the Device panel."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")

from PyQt6.QtWidgets import QApplication  # noqa: E402

from prestus_nexus.equipment import EquipmentLibrary  # noqa: E402
from prestus_nexus.model import ConfigModel  # noqa: E402
from prestus_nexus.session import Session  # noqa: E402
from prestus_nexus.ui.main_window import MainWindow  # noqa: E402


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def window(app, schema, prestus_root, tmp_path):
    w = MainWindow(ConfigModel(schema, EquipmentLibrary(prestus_root)), Session(tmp_path / "s.yaml"),
                   prestus_root)
    yield w
    w.dirty = False
    w.close()


def test_tabs_span_the_window(window):
    bar = window.tabs.tabBar()
    assert bar.expanding()
    assert bar.label(0) == "QuickStart"
    assert not bar.label(1)[0].isdigit()                           # plain labels, no step numbers


def test_tab_problems_colour_the_label_and_count_by_advanced(window):
    """A tab with problems keeps its plain label in the status colour; the count sits next
    to the tab's Advanced button."""
    bar = window.tabs.tabBar()
    flagged = 0
    for i in range(window.tabs.count()):
        tab = window.tabs.widget(i)
        if not hasattr(tab, "issue_badge"):
            continue
        n = sum(1 for iss in window.issues if window._issue_tab(iss) == tab.key
                and iss.level in ("error", "warning"))
        assert not any(ch.isdigit() for ch in bar.label(i))
        assert bar.tabIcon(i).isNull() or bar.is_compact()          # no status dot
        assert (bar.label_color(i) is not None) == bool(n)
        assert tab.issue_badge.isVisibleTo(tab) == bool(n)
        if n:
            flagged += 1
            assert tab.issue_badge.text().split()[-1].isdigit()
    assert flagged                          # an untitled study has problems (no subject, paths)


def test_placement_opens_on_settings_first(window):
    tab = window.tab_pages["placement"]
    subs = tab.subtabs
    assert subs.tabText(0) == "Settings"
    assert subs.widget(0) is tab.settings_page
    assert subs.currentWidget() is tab.settings_page
    assert subs.tabText(1) == "Viewer"


def test_positions_are_edited_per_transducer_on_placement(window):
    page = window.pages["placement"]
    assert page.has_tx and page.device_panel is None
    window.model.add_transducer()
    page.refresh()
    assert page.tx_combo.count() == 2
    page.tx_combo.setCurrentIndex(1)
    row = next(r for r in page.rows if r.node.path == "transducer.focus_pos")
    row.editor.edit.setText("1, 2, 3")
    row.editor._commit()
    assert window.model.get("transducer.focus_pos", 1) == [1, 2, 3]
    assert window.model.get("transducer.focus_pos", 0) in ([], None)


def test_serial_is_edited_only_in_the_device_panel(window):
    tx = window.pages["transducer"]
    tx.refresh()
    shown = {r.node.path for r in tx.rows if tx._row_shown(r)}
    assert "transducer.serial" not in shown and "transducer.combo.ds_serial" not in shown
    assert tx.device_panel is not None
