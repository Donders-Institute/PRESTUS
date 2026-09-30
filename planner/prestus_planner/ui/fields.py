"""Editors for single config values, chosen by schema kind."""

from __future__ import annotations

from typing import Any

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (QCheckBox, QComboBox, QFileDialog, QHBoxLayout, QLineEdit,
                             QToolButton, QWidget)

from .. import yamlio
from ..checks import check_value
from ..model import is_empty
from ..schema import Node

TRISTATE = [("auto (~)", None), ("true", True), ("false", False)]


def _num_text(value: Any) -> str:
    """Scalar for display: integral floats without '.0' (300000.0 -> 300000)."""
    if isinstance(value, float) and value.is_integer() and abs(value) < 1e15:
        return str(int(value))
    return yamlio.format_scalar(value)


def display_text(value: Any) -> str:
    """Text shown in a line edit for a value (plain strings unquoted)."""
    if isinstance(value, str):
        return value
    if isinstance(value, list):
        if any(isinstance(v, (list, dict)) for v in value):
            return yamlio._flow(value)          # nested lists: explicit YAML
        return ", ".join(_num_text(v) for v in value)
    if isinstance(value, dict):
        return yamlio.dump(value).strip()
    return _num_text(value)


class FieldEditor(QWidget):
    """Base class. Emits ``edited(value)`` with a parsed, valid value."""

    edited = pyqtSignal(object)

    def __init__(self, node: Node, parent: QWidget | None = None):
        super().__init__(parent)
        self.node = node
        self.baseline: Any = node.default
        self._layout = QHBoxLayout(self)
        self._layout.setContentsMargins(0, 0, 0, 0)

    def set_value(self, value: Any, baseline: Any) -> None:
        raise NotImplementedError

    def set_error(self, message: str | None) -> None:
        pass


class TextEditor(FieldEditor):
    """Line edit that parses its text according to the node kind."""

    def __init__(self, node: Node, parent: QWidget | None = None):
        super().__init__(node, parent)
        self.edit = QLineEdit()
        self.edit.setClearButtonEnabled(node.kind in ("str", "dir", "file", "path", "strlist"))
        self._layout.addWidget(self.edit, 1)
        self.edit.editingFinished.connect(self._commit)
        self.edit.textEdited.connect(lambda _: self.set_error(None))
        if node.kind in ("dir", "file", "path"):
            browse = QToolButton(text="…")
            browse.setToolTip("Browse")
            browse.clicked.connect(self._browse)
            self._layout.addWidget(browse)
        placeholder = {"vec3": "x, y, z", "numbers": "value or v1, v2, …",
                       "strlist": "item1; item2", "yaml": "YAML value"}.get(node.kind)
        if placeholder:
            self.edit.setPlaceholderText(placeholder)
        self._shown = ""

    def set_value(self, value: Any, baseline: Any) -> None:
        self.baseline = baseline
        text = "; ".join(value) if self.node.kind == "strlist" and isinstance(value, list) \
            else display_text(value)
        self._shown = text
        if self.edit.text() != text:
            self.edit.setText(text)
        self.set_error(None)

    def set_error(self, message: str | None) -> None:
        self.edit.setStyleSheet("QLineEdit { border: 1px solid #d9534f; }" if message else "")
        self.edit.setToolTip(message or "")

    def parse(self, text: str) -> Any:
        kind = self.node.kind
        text = text.strip()
        if not text:
            # Empty input restores the baseline's kind of emptiness ([] vs '')
            return self.baseline if is_empty(self.baseline) else \
                ("" if kind in ("str", "dir", "file", "path") else [])
        if kind in ("str", "dir", "file", "path"):
            return text
        if kind == "strlist":
            return [p.strip() for p in text.replace(",", ";").split(";") if p.strip()]
        if kind == "yaml":
            return yamlio.load(text)
        # numeric kinds: accept '1, 2, 3', '[1, 2, 3]', '1 2 3' or a scalar
        body = text if text.startswith("[") else "[" + ", ".join(text.replace(",", " ").split()) + "]"
        value = yamlio.load(body)
        if kind in ("number", "float", "int") or (kind == "numbers" and len(value) == 1
                                                   and not text.startswith("[")):
            if len(value) != 1:
                raise ValueError("expected a single number")
            value = value[0]
        if kind == "int" and isinstance(value, float) and value.is_integer():
            value = int(value)
        return value

    def _commit(self) -> None:
        text = self.edit.text()
        if text == self._shown:
            return
        try:
            value = self.parse(text)
        except Exception as exc:          # YAML or number parse error
            self.set_error(f"Cannot parse: {exc}")
            return
        problem = check_value(self.node, value)
        if problem:
            self.set_error(problem)
            return
        self._shown = text
        self.edited.emit(value)

    def _browse(self) -> None:
        start = self.edit.text() or ""
        if self.node.kind == "dir":
            chosen = QFileDialog.getExistingDirectory(self, self.node.path, start)
        elif self.node.kind == "file":
            chosen, _ = QFileDialog.getOpenFileName(self, self.node.path, start)
        else:
            chosen = QFileDialog.getExistingDirectory(self, self.node.path, start)
        if chosen:
            self.edit.setText(chosen)
            self._commit()


class ChoiceEditor(FieldEditor):
    def __init__(self, node: Node, parent: QWidget | None = None):
        super().__init__(node, parent)
        self.combo = QComboBox()
        self._layout.addWidget(self.combo, 1)
        self._items: list[Any] = []
        self._setting = False
        self.combo.currentIndexChanged.connect(self._changed)

    def _options(self, value: Any) -> list[tuple[str, Any]]:
        if self.node.kind == "tristate":
            return list(TRISTATE)
        opts = [(str(c), c) for c in self.node.choices or []]
        if not is_empty(value) and value not in (self.node.choices or []):
            opts.append((f"{value}  (not recognised)", value))
        if is_empty(value):
            opts.insert(0, ("(not set)", value))
        return opts

    def set_value(self, value: Any, baseline: Any) -> None:
        self.baseline = baseline
        self._setting = True
        opts = self._options(value)
        self._items = [v for _, v in opts]
        self.combo.clear()
        self.combo.addItems([label for label, _ in opts])
        for i, v in enumerate(self._items):
            if type(v) is type(value) and v == value:
                self.combo.setCurrentIndex(i)
                break
        self._setting = False

    def _changed(self, index: int) -> None:
        if self._setting or index < 0:
            return
        self.edited.emit(self._items[index])


class CheckEditor(FieldEditor):
    """Checkbox; 'flag' writes 0/1, 'bool' writes true/false."""

    def __init__(self, node: Node, parent: QWidget | None = None):
        super().__init__(node, parent)
        self.box = QCheckBox()
        self._layout.addWidget(self.box)
        self._layout.addStretch(1)
        self._setting = False
        self.box.toggled.connect(self._changed)

    def set_value(self, value: Any, baseline: Any) -> None:
        self.baseline = baseline
        self._setting = True
        self.box.setChecked(bool(value) if not is_empty(value) else False)
        self._setting = False

    def _changed(self, checked: bool) -> None:
        if self._setting:
            return
        self.edited.emit((1 if checked else 0) if self.node.kind == "flag" else checked)


def make_editor(node: Node, parent: QWidget | None = None) -> FieldEditor:
    if node.kind == "tristate" or node.choices:
        return ChoiceEditor(node, parent)
    if node.kind in ("flag", "bool"):
        return CheckEditor(node, parent)
    return TextEditor(node, parent)
