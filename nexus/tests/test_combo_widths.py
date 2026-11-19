"""Every drop-down in the Nexus shows its longest item whole: style.apply fits each non-editable
QComboBox (style.fit_combo) with slack for platform fonts that render wider than the test's."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")

from PyQt6.QtWidgets import (QApplication, QComboBox, QStyle, QStyleOptionComboBox,  # noqa: E402
                             QTabWidget)

from prestus_nexus.model import ConfigModel  # noqa: E402
from prestus_nexus.session import Session  # noqa: E402
from prestus_nexus.ui import style  # noqa: E402

MIN_SLACK = 4      # px of the item area left over by the longest item (style.COMBO_SLACK is 8)


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture(scope="module")
def window(app, schema, prestus_root, tmp_path_factory):
    from prestus_nexus.ui.main_window import MainWindow
    style.apply(app, prestus_root)
    w = MainWindow(ConfigModel(schema), Session(tmp_path_factory.mktemp("s") / "session.yaml"), prestus_root)
    w.resize(1000, 800)        # narrow: neighbours that stretch squeeze the combos most
    w.show()
    yield w
    w.dirty = False
    w.close()


def _text_need(c: QComboBox) -> int:
    fm = c.fontMetrics()
    need = max(fm.horizontalAdvance(c.itemText(i)) for i in range(c.count()))
    if any(not c.itemIcon(i).isNull() for i in range(c.count())):
        need += c.iconSize().width() + 4
    return need


def _item_area(c: QComboBox) -> int:
    opt = QStyleOptionComboBox()
    c.initStyleOption(opt)
    return c.style().subControlRect(QStyle.ComplexControl.CC_ComboBox, opt,
                                    QStyle.SubControl.SC_ComboBoxEditField, c).width()


def _visible_combos(app, window) -> dict:
    """Each tab and sub-tab in turn: the combos it shows, laid out at the window's size."""
    seen = {}

    def collect():
        app.processEvents()
        for c in window.findChildren(QComboBox):
            if c.isVisible() and not c.isEditable() and c.count():
                seen[id(c)] = c
    for i in range(window.tabs.count()):
        window.tabs.setCurrentIndex(i)
        collect()
        for sub in window.tabs.widget(i).findChildren(QTabWidget):
            for j in range(sub.count()):
                sub.setCurrentIndex(j)
                collect()
    return seen


def test_every_combo_shows_its_longest_item(app, window):
    combos = _visible_combos(app, window)
    assert len(combos) >= 10, "the walk over tabs found too few combos to mean anything"
    clipped = []
    for c in combos.values():
        need = _text_need(c)
        if style.combo_fit_width(c) > style.COMBO_MAX_FIT:
            continue                                   # long paths may clip by design
        if _item_area(c) < need + MIN_SLACK:
            longest = max((c.itemText(i) for i in range(c.count())), key=len)
            clipped.append(f"{type(c.parentWidget()).__name__}: {longest!r} "
                           f"({_item_area(c)} px for {need} px)")
    assert not clipped, "combos that clip their longest item:\n" + "\n".join(clipped)


def test_combo_refits_when_items_change(app):
    style.apply(app)
    c = QComboBox()
    c.addItems(["Free water", "Head model"])
    c.show()
    app.processEvents()
    c.setItemText(0, "Localite neuronavigation  (needs setup)")   # relabelled after show
    app.processEvents()
    assert _item_area(c) >= _text_need(c) + MIN_SLACK
    c.close()
