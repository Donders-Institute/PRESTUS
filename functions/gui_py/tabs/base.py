"""Base class for all PRESTUS GUI tabs."""

from __future__ import annotations
from PyQt6.QtWidgets import (
    QWidget, QScrollArea, QVBoxLayout, QFormLayout, QGridLayout,
    QLabel, QFrame, QSizePolicy,
)
from PyQt6.QtCore import Qt
import style
from widgets import make_section, make_note, make_label, make_hline


class BaseTab(QScrollArea):
    """Scrollable container that hosts a form widget."""

    def __init__(self, registry: dict, parent=None):
        super().__init__(parent)
        self.registry = registry
        self.setWidgetResizable(True)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.setStyleSheet("border: none; background: white;")

        self._inner = QWidget()
        self._inner.setStyleSheet(f"background:{style.BG_PANEL};")
        self._layout = QVBoxLayout(self._inner)
        self._layout.setContentsMargins(24, 16, 24, 24)
        self._layout.setSpacing(2)
        self.setWidget(self._inner)

    def _add_section(self, text: str):
        self._layout.addSpacing(8)
        self._layout.addWidget(make_section(text))
        self._layout.addSpacing(2)

    def _add_note(self, text: str):
        self._layout.addWidget(make_note(text))

    def _add_row(self, label: str, widget, note: str = ""):
        row = QWidget()
        row.setStyleSheet("background: transparent;")
        hl = QGridLayout(row)
        hl.setContentsMargins(0, 1, 0, 1)
        hl.setColumnStretch(1, 1)
        hl.setColumnMinimumWidth(0, 210)

        lbl = make_label(label)
        hl.addWidget(lbl, 0, 0)
        hl.addWidget(widget, 0, 1)
        if note:
            n = make_note(note)
            hl.addWidget(n, 1, 1)
        self._layout.addWidget(row)

    def _add_row_wide(self, widget, note: str = ""):
        """Full-width widget (no label column)."""
        self._layout.addWidget(widget)
        if note:
            self._layout.addWidget(make_note(note))

    def _add_spacer(self):
        self._layout.addSpacing(6)

    def _finish(self):
        self._layout.addStretch()
