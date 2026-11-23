"""Editors for single config values, chosen by schema kind."""

from __future__ import annotations

import html
import re
from typing import Any

from PyQt6.QtCore import QEvent, Qt, pyqtSignal
from PyQt6.QtWidgets import (QSizePolicy, QCheckBox, QComboBox, QFileDialog, QHBoxLayout, QLabel, QLineEdit,
                             QToolTip, QVBoxLayout, QWidget)

from .. import labels, yamlio
from ..checks import check_value
from ..model import is_empty
from ..schema import Node
from . import icons, style
from .tooltip import choice_tip

TRISTATE = [("auto (~)", None), ("true", True), ("false", False)]
PATH_KINDS = ("dir", "file", "path")


def full_value_tip(text: str) -> str:
    """Hover text with a path field's whole value; a ':'-separated search path
    (LD_LIBRARY_PATH) gets one folder per line."""
    parts = text.split(":")
    if len(parts) > 1 and all(not p or re.match(r"[/$~]", p) for p in parts):
        return "<br>".join(html.escape(p) for p in parts)
    return text


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
        if node.kind in PATH_KINDS:
            self.edit.installEventFilter(self)      # hover: the whole value when it is clipped
            browse = icons.PathButton("file" if node.kind == "file" else "dir")
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
            if self.node.kind in PATH_KINDS and not self.edit.hasFocus():
                self.edit.setCursorPosition(0)      # show the start of a long path, not its tail
        self.set_error(None)

    def eventFilter(self, obj, event) -> bool:     # noqa: N802
        """A path field whose value does not fit shows all of it on hover (an error message,
        when there is one, keeps the tooltip)."""
        if obj is self.edit and event.type() == QEvent.Type.ToolTip and not self.edit.toolTip():
            text = self.edit.text()
            room = self.edit.contentsRect().width() - 8     # the frame's text margins
            if text and self.edit.fontMetrics().horizontalAdvance(text) > room:
                QToolTip.showText(event.globalPos(), full_value_tip(text), self.edit)
                return True
        return super().eventFilter(obj, event)

    def set_error(self, message: str | None) -> None:
        self.edit.setStyleSheet(f"QLineEdit {{ border: 1px solid {style.active()['red']}; }}" if message else "")
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
        except Exception as exc:  # reported: next to the field (YAML or number parse error)
            self.set_error(f"Cannot parse: {exc}")
            return
        problem = check_value(self.node, value)
        if problem:
            self.set_error(problem)
            return
        self.set_error(None)            # a value that arrives without typing (Browse…) clears it too
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


def _clear_layout(layout) -> None:
    """Remove a layout's widgets for good, at once: hidden first, so that a widget leaving
    its parent can never show as a window of its own (it did: a stray note window per
    refresh), then deleted."""
    while layout.count():
        w = layout.takeAt(0).widget()
        if w is not None:
            w.hide()
            w.setParent(None)
            w.deleteLater()


FIGURE_PX = 120        # width of a choice figure; its caption sits beside it


class ChoiceEditor(FieldEditor):
    def __init__(self, node: Node, parent: QWidget | None = None):
        super().__init__(node, parent)
        self.combo = QComboBox()
        self.combo.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)  # as wide as text fields
        box = QVBoxLayout()
        box.setContentsMargins(0, 0, 0, 0)
        box.setSpacing(2)
        box.addWidget(self.combo)
        self.note = QLabel()            # why the selected choice cannot run yet (set_item_status)
        self.note.setWordWrap(True)
        self.note.setVisible(False)
        box.addWidget(self.note)
        self.math = QVBoxLayout()       # overlay 'choice_math': the selected choice's formula
        self.math.setContentsMargins(0, 4, 0, 0)
        self.math.setSpacing(2)
        box.addLayout(self.math)
        self.figures = QVBoxLayout()    # overlay 'choice_figures': one row per schematic,
        self.figures.setContentsMargins(0, 4, 0, 0)   # its caption beside it
        self.figures.setSpacing(8)
        box.addLayout(self.figures)
        self._layout.addLayout(box, 1)
        self._items: list[Any] = []
        self._setting = False
        self.combo.currentIndexChanged.connect(self._changed)

    def _options(self, value: Any) -> list[tuple[str, Any]]:
        if self.node.kind == "tristate":
            return list(TRISTATE)
        opts = [(labels.choice(self.node, c), c) for c in self.node.choices or []]
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
        if self.node.choice_icons:
            from . import icons
            for i, v in enumerate(self._items):
                ic = icons.for_choice(self.node, v)
                if ic is not None:
                    self.combo.setItemIcon(i, ic)
        for i, v in enumerate(self._items):
            if self.node.kind != "tristate" and not is_empty(v):
                self.combo.setItemData(i, choice_tip(self.node, v), Qt.ItemDataRole.ToolTipRole)
        for i, v in enumerate(self._items):
            if type(v) is type(value) and v == value:
                self.combo.setCurrentIndex(i)
                break
        self._setting = False
        self._show_math()

    def _show_math(self) -> None:
        """The formula lines of the selected choice (overlay 'choice_math'), rendered."""
        _clear_layout(self.math)
        self._show_figures()
        if not self.node.choice_math or self.combo.currentIndex() < 0:
            return
        lines = self.node.choice_math.get(self._items[self.combo.currentIndex()]) or []
        from PyQt6.QtGui import QPixmap
        from .pct_formulas import render_math
        for line in [lines] if isinstance(lines, str) else lines:
            label = QLabel(self)            # parented at once: never a window of its own
            text = line[2:] if line.startswith("> ") else None      # '> ': a plain-text note
            img = render_math(line, style.active()["ink"], 11) if text is None else None
            if img is not None:
                label.setPixmap(QPixmap.fromImage(img))
            else:
                label.setText(text if text is not None else line)   # LaTeX source as fallback
                label.setWordWrap(True)
                label.setObjectName("CardSub")
            self.math.addWidget(label)

    def _show_figures(self) -> None:
        """The figures of the selected choice (overlay 'choice_figures'), each with its caption."""
        _clear_layout(self.figures)
        if not self.node.choice_figures or self.combo.currentIndex() < 0:
            return
        value = self._items[self.combo.currentIndex()]
        names = self.node.choice_figures.get(value, []) if isinstance(value, str) else []
        if not names:
            return
        from . import choice_icons, icons
        for name in names:
            cell = QWidget(self)
            row = QHBoxLayout(cell)
            row.setContentsMargins(0, 0, 0, 0)
            row.setSpacing(12)
            pic = QLabel()
            pic.setPixmap(icons.figure(name, FIGURE_PX))
            pic.setFixedSize(pic.sizeHint())
            text = QLabel(choice_icons.caption(name))
            text.setWordWrap(True)
            text.setObjectName("CardSub")
            text.setMinimumWidth(FIGURE_PX)
            text.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
            cell.setMinimumHeight(pic.sizeHint().height())   # wrapped text must not squeeze it
            row.addWidget(pic, 0, Qt.AlignmentFlag.AlignVCenter)
            row.addWidget(text, 1, Qt.AlignmentFlag.AlignVCenter)
            self.figures.addWidget(cell)

    def _changed(self, index: int) -> None:
        if self._setting or index < 0:
            return
        self._show_math()
        self.edited.emit(self._items[index])

    def set_item_status(self, status: dict) -> None:
        """Mark choices by availability (checks.ModeStatus per value): 'unavailable'
        items are disabled, 'pending' ones get a note; the reason is the tooltip.
        The current value always stays selectable."""
        model = self.combo.model()
        current = self._items[self.combo.currentIndex()] if self.combo.currentIndex() >= 0 else None
        for i, v in enumerate(self._items):
            st = status.get(v) if isinstance(v, str) else None
            item = model.item(i) if hasattr(model, "item") else None
            if st is None or item is None:
                continue
            label = labels.choice(self.node, v)
            if st.state == "unavailable":
                label += "  (unavailable)"
                item.setEnabled(v == current)
            elif st.state == "pending":
                label += "  (needs setup)"
            item.setText(label)
            item.setToolTip(choice_tip(self.node, v, st.reason or "All inputs found"))
        lines = []
        st = status.get(current) if isinstance(current, str) else None
        if st is not None and st.state != "ready" and st.reason:
            colour = style.active()["red" if st.state == "unavailable" else "amber"]
            what = "Cannot run yet" if st.state == "unavailable" else "Still to set up"
            lines.append(f"<span style='color:{colour}'>{what}:</span> {html.escape(st.reason)}")
        # Disabled choices show no tooltip on most platforms: say here why they are unavailable
        for v in self._items:
            other = status.get(v) if isinstance(v, str) and v != current else None
            if other is not None and other.state == "unavailable" and other.reason:
                lines.append(f"<span style='color:{style.active()['ink3']}'>{html.escape(labels.choice(self.node, v))} "
                             f"unavailable: {html.escape(other.reason)}</span>")
        self.note.setText("<br>".join(lines))
        self.note.setVisible(bool(lines))


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
    if node.path == "transducer.target_isppa_wcm2":    # line edit + "Sweep…" (multi-ISPPA)
        from .isppa_sweep import IsppaSweepEditor
        return IsppaSweepEditor(node, parent)
    if node.kind == "tristate" or node.choices:
        return ChoiceEditor(node, parent)
    if node.kind in ("flag", "bool"):
        return CheckEditor(node, parent)
    return TextEditor(node, parent)
