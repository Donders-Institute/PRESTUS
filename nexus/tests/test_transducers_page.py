"""Transducers overview (Placement ▸ Settings, or a page of its own), offscreen."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")

from PyQt6.QtWidgets import QApplication  # noqa: E402

from prestus_nexus.model import ConfigModel  # noqa: E402
from prestus_nexus.ui.transducers_page import TransducersPage  # noqa: E402


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def _model(schema):
    m = ConfigModel(schema)
    m.set("transducer.trans_pos", [90, 100, 150])
    m.set("transducer.focus_pos", [100, 100, 100])
    m.set("transducer.freq_hz", 500000)
    m.add_transducer()
    m.set("transducer.trans_pos", [110, 100, 150], 1)
    m.set("transducer.focus_pos", [100, 100, 100], 1)
    return m


def test_one_transducer_shows_the_empty_state(app, schema):
    page = TransducersPage(ConfigModel(schema))
    page.refresh()
    assert page.stack.is_empty()


def test_table_tilts_and_mode_change(app, schema):
    m = _model(schema)
    page = TransducersPage(m)
    page.refresh()
    assert not page.stack.is_empty()
    assert page.table.rowCount() == 2
    t1 = page.table.item(0, 6).text()
    assert t1 == page.table.item(1, 6).text() and t1.startswith("11.3")    # atan(10/50)
    assert "Grid z" in page.orient_text.text()
    page.mode_editor.edited.emit("transducer_axis")
    assert m.get("grid.mode") == "transducer_axis"
    assert page.table.item(0, 6).text().startswith("0.0")
    assert page.table.item(1, 6).text().startswith("22.6")
    page.mode_editor.edited.emit("ras_plus")
    assert page.table.item(0, 6).text().startswith("11.3")
    rows = []
    page.edit_in_viewer.connect(rows.append)
    page.table.selectRow(1)
    page.edit_btn.click()
    assert rows == [1]


def test_async_layout_row_follows_coupling(app, schema):
    m = _model(schema)
    page = TransducersPage(m)
    page.refresh()
    assert page.async_row.isHidden()
    m.set("simulation.transducer_coupling", "async")
    page.refresh()
    assert not page.async_row.isHidden()
    page.async_editor.edited.emit("per_transducer")
    assert m.get("grid.async_layout") == "per_transducer"


def test_coarse_tilted_beams_are_flagged(app, schema):
    m = _model(schema)
    m.set("transducer.trans_pos", [150, 100, 150], 1)     # 56 deg apart: each ~28 deg off z
    m.set("grid.resolution_mm", 1.0)
    page = TransducersPage(m)
    page.refresh()
    assert not page.orient_note.isHidden()
    assert "points per wavelength" in page.orient_note.text()


def test_head_sketch_beside_the_table(app, schema):
    """The sketch shows every positioned transducer, left of the other one here, and hides
    when no transducer has a position."""
    m = _model(schema)
    page = TransducersPage(m)
    page.refresh()
    heads = page.sketch.heads
    assert len(heads) == 2 and not page.sketch.isHidden()
    assert heads[0][0][0] < heads[1][0][0]               # x 90 is left of x 110 (RAS+)
    assert not page.sketch.grab().isNull()
    for ti in range(2):
        m.reset("transducer.trans_pos", ti)
    page.refresh()
    assert page.sketch.isHidden()


def test_compact_table(app, schema, prestus_root):
    """Numbered rows in the sketch colours, the device's short name, units in the headers."""
    from prestus_nexus.equipment import EquipmentLibrary
    from prestus_nexus.ui import style
    m = _model(schema)
    m.library = EquipmentLibrary(prestus_root)
    m.set("transducer.serial", "CTX_500_026")
    page = TransducersPage(m)
    page.refresh()
    t = page.table
    heads = [t.horizontalHeaderItem(c).text() for c in range(t.columnCount())]
    assert heads[0] == "" and heads[1] == "Transducer" and heads[2] == "Position (vox)"
    assert [t.item(r, 0).text() for r in range(2)] == ["1", "2"]
    assert t.item(1, 0).background().color().name().lower() == style.series(1).lower()
    assert t.item(0, 1).text() == "CTX-500 #026" and "CTX_500_026" in t.item(0, 1).toolTip()
    assert t.item(0, 2).text() == "90, 100, 150"              # no brackets, unit in the header
    assert m.transducer_title(1) == "Transducer 2: CTX-500 #026"
