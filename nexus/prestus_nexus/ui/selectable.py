"""Every label's text can be selected with the mouse and copied (Ctrl/Cmd+C or right-click ▸ Copy).

One application-wide event filter adds ``TextSelectableByMouse`` to each ``QLabel`` when Qt first
polishes it, so labels built anywhere (pages, dialogs, message boxes) need no per-widget call.
Existing flags such as ``LinksAccessibleByMouse`` are kept. Labels that show only a pixmap or
movie are left alone.
"""

from __future__ import annotations

from PyQt6.QtCore import QEvent, QObject, Qt
from PyQt6.QtWidgets import QApplication, QLabel

_SELECT = Qt.TextInteractionFlag.TextSelectableByMouse


def make_selectable(label: QLabel) -> None:
    """Add mouse selection to *label*, keeping its other interaction flags."""
    if label.pixmap() is not None and not label.pixmap().isNull():
        return
    if label.movie() is not None:
        return
    flags = label.textInteractionFlags()
    if flags & _SELECT:
        return
    if flags == Qt.TextInteractionFlag.NoTextInteraction:
        flags = _SELECT
    else:
        flags |= _SELECT
    label.setTextInteractionFlags(flags)


class _SelectableLabels(QObject):
    def eventFilter(self, obj, event) -> bool:          # noqa: N802 (Qt API)
        if event.type() == QEvent.Type.Polish and isinstance(obj, QLabel):
            make_selectable(obj)
        return False


def install(app: QApplication) -> None:
    """Make all labels of *app* selectable; installing twice is a no-op."""
    if getattr(app, "_selectable_labels", None) is not None:
        return
    f = _SelectableLabels(app)
    app.installEventFilter(f)
    app._selectable_labels = f
