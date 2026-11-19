"""Stages card (section.CHECK_FIRST): each stage's checkbox sits before its label; rows without a
checkbox keep an empty slot so the labels line up. Every other card keeps the label | input layout."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")

from PyQt6.QtWidgets import QApplication  # noqa: E402

from prestus_nexus.model import ConfigModel  # noqa: E402
from prestus_nexus.session import Session  # noqa: E402
from prestus_nexus.ui import style  # noqa: E402
from prestus_nexus.ui.fields import CheckEditor  # noqa: E402
from prestus_nexus.ui.section import CHECK_FIRST, CHECK_SLOT, SectionPage  # noqa: E402


@pytest.fixture(scope="module")
def pages(schema, prestus_root, tmp_path_factory):
    app = QApplication.instance() or QApplication([])
    from prestus_nexus.ui.main_window import MainWindow
    style.apply(app, prestus_root)
    w = MainWindow(ConfigModel(schema), Session(tmp_path_factory.mktemp("s") / "session.yaml"), prestus_root)
    w.resize(1200, 900)
    w.show()
    app.processEvents()
    yield w.findChildren(SectionPage)
    w.dirty = False
    w.close()


def _line(row) -> list:
    """Widgets of the row's label line (checkbox slot, label), left to right."""
    line = row.block.layout().itemAt(0).layout()
    return [line.itemAt(i).widget() for i in range(line.count())]


def test_stage_checkbox_precedes_label(pages):
    stage_pages = [p for p in pages if p.section in CHECK_FIRST]
    assert stage_pages, "no Stages card built"
    for p in stage_pages:
        checks = [r for r in p.rows if isinstance(r.editor, CheckEditor)]
        assert checks, f"{p.section}: no checkbox rows"
        for r in checks:
            assert r.block.isAncestorOf(r.editor), f"{r.node.path}: checkbox not in the label column"
            line = _line(r)          # checkbox, the stage's icon, label
            assert line[0] is r.editor and line[-1] is r.label and len(line) == 3, \
                f"{r.node.path}: checkbox not before its icon and label"
        # rows without a checkbox keep an empty slot of the same width, so the labels line up
        for r in p.rows:
            slot, icon, label = _line(r)
            assert label is r.label and slot.width() == CHECK_SLOT == slot.minimumWidth(), r.node.path
            assert icon.width() == 18, r.node.path


def test_other_cards_keep_input_after_label(pages):
    others = [r for p in pages if p.section not in CHECK_FIRST for r in p.rows
              if isinstance(r.editor, CheckEditor)]
    assert others, "negative control has no checkbox rows"
    for r in others:
        assert r.container.isAncestorOf(r.editor), f"{r.node.path}: checkbox moved off its input column"
