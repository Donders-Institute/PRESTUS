"""One clear panel in place of a view that has nothing to show yet: what is missing, why,
and one button that fixes it. Pages swap it with their content in a QStackedWidget
(``EmptyStack``), so blank frames and placeholder axes never appear.

Text styles come from the app stylesheet's object names (CardTitle, CardSub, Card), so
light and dark mode need nothing here."""

from __future__ import annotations

from typing import Callable

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (QFrame, QHBoxLayout, QLabel, QPushButton, QSizePolicy,
                             QStackedWidget, QVBoxLayout, QWidget)


class EmptyState(QWidget):
    """A centred card: title (what is missing), body (why / what to do), one action."""

    def __init__(self, title: str = "", body: str = "", parent: QWidget | None = None):
        super().__init__(parent)
        self._action: Callable[[], None] | None = None
        outer = QVBoxLayout(self)
        outer.addStretch(1)
        row = QHBoxLayout()
        row.addStretch(1)
        card = QFrame()
        card.setObjectName("Card")
        card.setMaximumWidth(560)
        card.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        v = QVBoxLayout(card)
        v.setContentsMargins(28, 22, 28, 22)
        v.setSpacing(8)
        self.title = QLabel(title)
        self.title.setObjectName("CardTitle")
        self.title.setWordWrap(True)
        self.title.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        v.addWidget(self.title)
        self.body = QLabel(body)
        self.body.setObjectName("CardSub")
        self.body.setWordWrap(True)
        self.body.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        self.body.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        v.addWidget(self.body)
        brow = QHBoxLayout()
        brow.addStretch(1)
        self.button = QPushButton()
        self.button.setProperty("primary", True)
        self.button.clicked.connect(self._clicked)
        brow.addWidget(self.button)
        brow.addStretch(1)
        v.addSpacing(4)
        v.addLayout(brow)
        row.addWidget(card, 1)
        row.addStretch(1)
        outer.addLayout(row)
        outer.addStretch(2)
        self.button.hide()

    def set(self, title: str, body: str = "", action: str | None = None,
            callback: Callable[[], None] | None = None, tooltip: str = "") -> None:
        """Show this message; ``action``/``callback`` add the one button (else none)."""
        self.title.setText(title)
        self.body.setText(body)
        self.body.setVisible(bool(body))
        self._action = callback
        self.button.setText(action or "")
        self.button.setToolTip(tooltip)
        self.button.setVisible(bool(action and callback))

    def _clicked(self) -> None:
        if self._action is not None:
            self._action()


class EmptyStack(QStackedWidget):
    """Page 0: the content; page 1: an EmptyState. ``show_empty``/``show_content`` switch."""

    def __init__(self, content: QWidget, parent: QWidget | None = None):
        super().__init__(parent)
        self.content = content
        self.empty = EmptyState()
        self.addWidget(content)
        self.addWidget(self.empty)

    def show_empty(self, title: str, body: str = "", action: str | None = None,
                   callback: Callable[[], None] | None = None, tooltip: str = "") -> None:
        self.empty.set(title, body, action, callback, tooltip)
        self.setCurrentWidget(self.empty)

    def show_content(self) -> None:
        self.setCurrentWidget(self.content)

    def is_empty(self) -> bool:
        return self.currentWidget() is self.empty
