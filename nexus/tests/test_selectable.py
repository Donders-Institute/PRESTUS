"""Label text can be selected and copied everywhere in the Nexus."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")

from PyQt6.QtCore import Qt  # noqa: E402
from PyQt6.QtGui import QPixmap  # noqa: E402
from PyQt6.QtWidgets import QApplication, QLabel  # noqa: E402

from prestus_nexus.ui import selectable  # noqa: E402

SEL = Qt.TextInteractionFlag.TextSelectableByMouse
LINKS = Qt.TextInteractionFlag.LinksAccessibleByMouse


@pytest.fixture(scope="module")
def app():
    a = QApplication.instance() or QApplication([])
    selectable.install(a)
    selectable.install(a)                   # idempotent
    return a


def test_plain_label_becomes_selectable(app):
    lbl = QLabel("ISPPA 30 W/cm²")
    lbl.show()
    app.processEvents()
    assert lbl.textInteractionFlags() & SEL
    lbl.close()


def test_link_label_keeps_links(app):
    lbl = QLabel('<a href="x">go</a>')
    lbl.setTextInteractionFlags(LINKS)
    lbl.show()
    app.processEvents()
    assert lbl.textInteractionFlags() & SEL
    assert lbl.textInteractionFlags() & LINKS
    lbl.close()


def test_pixmap_label_left_alone(app):
    lbl = QLabel()
    pm = QPixmap(8, 8)
    pm.fill(Qt.GlobalColor.red)
    lbl.setPixmap(pm)
    lbl.show()
    app.processEvents()
    assert not lbl.textInteractionFlags() & SEL
    lbl.close()
