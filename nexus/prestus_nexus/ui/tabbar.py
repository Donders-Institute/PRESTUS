"""A tab bar that can draw chosen tabs' labels in bold (e.g. the sub-tab of the placement
method in use). QTabBar has no per-tab font, so the labels are painted here; tab shapes are
the style's own."""

from __future__ import annotations

from PyQt6.QtCore import QSize
from PyQt6.QtGui import QFontMetrics
from PyQt6.QtWidgets import QStyle, QStyleOptionTab, QStylePainter, QTabBar, QTabWidget


class EmphasisTabBar(QTabBar):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._bold: set[int] = set()

    def bold_tabs(self) -> set[int]:
        return set(self._bold)

    def set_bold(self, indices) -> None:
        indices = {i for i in indices if 0 <= i < self.count()}
        if indices == self._bold:
            return
        self._bold = indices
        for i in range(self.count()):           # tab widths follow the text: re-lay out
            self.setTabText(i, self.tabText(i))
        self.update()

    def _bold_font(self):
        f = self.font()
        f.setBold(True)
        return f

    def tabSizeHint(self, index: int) -> QSize:       # noqa: N802
        size = super().tabSizeHint(index)
        if index in self._bold:
            text = self.tabText(index)
            extra = QFontMetrics(self._bold_font()).horizontalAdvance(text) - \
                QFontMetrics(self.font()).horizontalAdvance(text)
            size.setWidth(size.width() + max(0, extra))
        return size

    def paintEvent(self, event) -> None:               # noqa: N802
        if not self._bold:
            super().paintEvent(event)
            return
        painter = QStylePainter(self)
        order = [i for i in range(self.count()) if i != self.currentIndex()] + \
            ([self.currentIndex()] if self.currentIndex() >= 0 else [])     # current on top
        for i in order:
            if not self.isTabVisible(i):
                continue
            opt = QStyleOptionTab()
            self.initStyleOption(opt, i)
            painter.drawControl(QStyle.ControlElement.CE_TabBarTabShape, opt)
            painter.setFont(self._bold_font() if i in self._bold else self.font())
            painter.drawControl(QStyle.ControlElement.CE_TabBarTabLabel, opt)


def tab_widget() -> QTabWidget:
    """A document-mode QTabWidget with an EmphasisTabBar (set before any tab is added)."""
    w = QTabWidget()
    w.setTabBar(EmphasisTabBar(w))
    w.setDocumentMode(True)
    return w
