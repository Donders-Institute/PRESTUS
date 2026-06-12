"""Reusable form-building helpers for the PRESTUS Python GUI.

All widgets that bind to a config parameter are registered in a shared
``registry`` dict (tag -> QWidget) owned by the main app.  The same tag
conventions as the MATLAB GUI are used (dotted paths like
``transducer.annular.elem_n``).
"""

from __future__ import annotations
from typing import Any

from PyQt6.QtWidgets import (
    QWidget, QLabel, QLineEdit, QDoubleSpinBox, QComboBox, QCheckBox,
    QPushButton, QHBoxLayout, QVBoxLayout, QFileDialog, QSizePolicy,
    QFrame,
)
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QDoubleValidator

import style


# ── Sentinel for "not set" (NaN equivalent) ──────────────────────────────────
_NAN_PLACEHOLDER = "NaN"


def _nan_or_float(text: str):
    text = text.strip()
    if text == "" or text.lower() == "nan":
        return None
    try:
        return float(text)
    except ValueError:
        return None


# ── Section header ────────────────────────────────────────────────────────────

def make_section(text: str) -> QLabel:
    lbl = QLabel("  " + text.upper())
    lbl.setProperty("class", "section")
    lbl.setStyleSheet(
        f"background:{style.ACCENT_LT}; color:{style.ACCENT_DK};"
        f"font-weight:bold; font-size:{style.FS_SM}px; padding:3px 8px;"
    )
    lbl.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
    lbl.setMinimumHeight(22)
    return lbl


def make_note(text: str) -> QLabel:
    lbl = QLabel(text)
    lbl.setWordWrap(True)
    lbl.setStyleSheet(f"color:{style.TEXT_SUB}; font-size:{style.FS_SM}px;")
    lbl.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
    return lbl


def make_label(text: str) -> QLabel:
    mandatory = text.startswith("*")
    clean = text.lstrip("* ").strip()
    if mandatory:
        lbl = QLabel(f'<span style="color:{style.MANDATORY}">*</span> {clean}')
        lbl.setTextFormat(Qt.TextFormat.RichText)
    else:
        lbl = QLabel(clean)
    lbl.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
    lbl.setStyleSheet(f"font-size:{style.FS}px;")
    return lbl


# ── BrowseEdit: text field + "···" browse button ──────────────────────────────

class BrowseEdit(QWidget):
    """A QLineEdit paired with a browse button for file/directory selection."""

    def __init__(self, tag: str, browse_type: str = "dir",
                 registry: dict | None = None, mandatory: bool = False,
                 parent=None):
        super().__init__(parent)
        self._tag = tag
        self._browse_type = browse_type

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)

        self.edit = QLineEdit()
        if mandatory:
            self.edit.setStyleSheet(
                f"background:#fff0f0; border:1px solid {style.MANDATORY};"
                f"border-radius:4px; padding:3px 6px;"
            )
        layout.addWidget(self.edit, stretch=1)

        btn = QPushButton("···")
        btn.setProperty("class", "browse")
        btn.setStyleSheet(
            f"background:#dde3f0; color:{style.TEXT};"
            "border-radius:4px; padding:4px 8px; min-width:36px;"
        )
        btn.setFixedWidth(42)
        btn.clicked.connect(self._browse)
        layout.addWidget(btn)

        if registry is not None:
            registry[tag] = self

    @property
    def tag(self):
        return self._tag

    def value(self) -> str:
        return self.edit.text()

    def set_value(self, v):
        self.edit.setText(str(v) if v is not None else "")

    def _browse(self):
        if self._browse_type == "dir":
            p = QFileDialog.getExistingDirectory(self, "Select folder")
            if p:
                self.edit.setText(p)
        elif self._browse_type == "file":
            p, _ = QFileDialog.getOpenFileName(self, "Select file")
            if p:
                self.edit.setText(p)
        else:
            p, _ = QFileDialog.getOpenFileName(
                self, "Select YAML", filter="YAML (*.yaml *.yml);;All files (*)")
            if p:
                self.edit.setText(p)


# ── FloatEdit: numeric field that accepts NaN ─────────────────────────────────

class FloatEdit(QLineEdit):
    """Single-line numeric field; empty/NaN → returns None."""

    def __init__(self, tag: str, default=None, registry: dict | None = None,
                 parent=None):
        super().__init__(parent)
        self._tag = tag
        self.setAlignment(Qt.AlignmentFlag.AlignRight)
        if default is None:
            self.setPlaceholderText(_NAN_PLACEHOLDER)
        else:
            self.setText(str(default))
        if registry is not None:
            registry[tag] = self

    @property
    def tag(self):
        return self._tag

    def value(self):
        return _nan_or_float(self.text())

    def set_value(self, v):
        if v is None or (isinstance(v, float) and v != v):  # NaN check
            self.clear()
        else:
            self.setText(str(v))


# ── TextEdit: plain string field ──────────────────────────────────────────────

class TextEdit(QLineEdit):
    def __init__(self, tag: str, default: str = "", mandatory: bool = False,
                 registry: dict | None = None, parent=None):
        super().__init__(parent)
        self._tag = tag
        self.setText(default)
        if mandatory:
            self.setStyleSheet(
                f"background:#fff0f0; border:1px solid {style.MANDATORY};"
                "border-radius:4px; padding:3px 6px;"
            )
        if registry is not None:
            registry[tag] = self

    @property
    def tag(self):
        return self._tag

    def value(self) -> str:
        return self.text()

    def set_value(self, v):
        self.setText(str(v) if v is not None else "")


# ── Dropdown ──────────────────────────────────────────────────────────────────

class Dropdown(QComboBox):
    def __init__(self, tag: str, items: list[str], default: str | None = None,
                 registry: dict | None = None, parent=None):
        super().__init__(parent)
        self._tag = tag
        self.addItems(items)
        if default and default in items:
            self.setCurrentText(default)
        if registry is not None:
            registry[tag] = self

    @property
    def tag(self):
        return self._tag

    def value(self) -> str:
        return self.currentText()

    def set_value(self, v):
        if v is not None:
            idx = self.findText(str(v))
            if idx >= 0:
                self.setCurrentIndex(idx)
            else:
                # Add it if not present (e.g. serial from saved config)
                self.addItem(str(v))
                self.setCurrentText(str(v))


# ── Checkbox ──────────────────────────────────────────────────────────────────

class Checkbox(QCheckBox):
    def __init__(self, tag: str, label: str = "", default: bool = False,
                 registry: dict | None = None, parent=None):
        super().__init__(label, parent)
        self._tag = tag
        self.setChecked(default)
        if registry is not None:
            registry[tag] = self

    @property
    def tag(self):
        return self._tag

    def value(self) -> bool:
        return self.isChecked()

    def set_value(self, v):
        self.setChecked(bool(v))


# ── XYZWidget: three numeric fields side-by-side for 3D coordinates ───────────

class XYZWidget(QWidget):
    """Three FloatEdit fields for X/Y/Z, tagged as base_tag_0/1/2."""

    def __init__(self, base_tag: str, defaults=(None, None, None),
                 registry: dict | None = None, parent=None):
        super().__init__(parent)
        self._base_tag = base_tag
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)
        self._fields = []
        for i, axis in enumerate(("X", "Y", "Z")):
            lbl = QLabel(axis)
            lbl.setStyleSheet(f"color:{style.TEXT_SUB}; font-size:{style.FS_SM}px;")
            lbl.setFixedWidth(14)
            layout.addWidget(lbl)
            tag = f"{base_tag}_{i}"
            fe = FloatEdit(tag, defaults[i] if defaults[i] is not None else None)
            layout.addWidget(fe, stretch=1)
            self._fields.append(fe)
            if registry is not None:
                registry[tag] = fe

    def value(self):
        return [f.value() for f in self._fields]

    def set_value(self, v):
        if v is None:
            return
        for i, val in enumerate(v[:3]):
            self._fields[i].set_value(val)


# ── CSVArrayEdit: comma-separated list of numbers ─────────────────────────────

class CSVArrayEdit(QLineEdit):
    """Displays/edits a numeric array as a comma-separated string."""

    def __init__(self, tag: str, default: str = "", registry: dict | None = None,
                 parent=None):
        super().__init__(parent)
        self._tag = tag
        self.setText(default)
        if registry is not None:
            registry[tag] = self

    @property
    def tag(self):
        return self._tag

    def value(self):
        raw = self.text().strip()
        if not raw:
            return None
        try:
            return [float(x) for x in raw.replace(";", ",").split(",") if x.strip()]
        except ValueError:
            return None

    def set_value(self, v):
        if v is None:
            self.clear()
        elif isinstance(v, (list, tuple)):
            self.setText(", ".join(str(x) for x in v))
        elif isinstance(v, str):
            self.setText(v)
        else:
            self.setText(str(v))


# ── Horizontal divider ────────────────────────────────────────────────────────

def make_hline() -> QFrame:
    line = QFrame()
    line.setFrameShape(QFrame.Shape.HLine)
    line.setStyleSheet(f"color:{style.BORDER};")
    return line
