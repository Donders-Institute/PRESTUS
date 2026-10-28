"""The form for one config section, shown as a card on a tab of the main window."""

from __future__ import annotations

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QCursor
from PyQt6.QtWidgets import (QComboBox, QFrame, QGridLayout, QHBoxLayout, QLabel, QMessageBox, QPushButton,
                             QScrollArea, QToolButton, QToolTip, QVBoxLayout, QWidget)

from .. import labels
from ..checks import check_value, placement_mode_status
from ..model import TX, ConfigModel
from ..schema import Node
from . import style
from .fields import FieldEditor, display_text, make_editor
from .transducer_panel import DevicePanel

LABEL_WIDTH = 290
# Inputs stay near their label instead of stretching across a wide window: values and
# choices get a short field, paths and lists a longer one (both still shrink when narrow).
FIELD_MAX_WIDTH = {"dir": 760, "file": 760, "path": 760, "str": 560, "strlist": 560, "yaml": 760, "numbers": 420}
FIELD_MAX_DEFAULT = 380


def _tooltip(node: Node) -> str:
    """Full help for hover: description, key, default."""
    title = labels.label(node)
    parts = [f"<b>{title}</b>"]
    if node.help:
        parts.append(node.help.replace("\n", "<br>"))
    parts.append(f"<span style='color:{style.active()['ink4']}'>Key: {node.path}</span>")
    if not node.is_group:
        parts.append(f"<span style='color:{style.active()['ink4']}'>Default: {display_text(node.default) or '(empty)'}</span>")
    if node.tune:
        parts.append(f"<span style='color:{style.active()['ink4']}'>Tuning: affects {', '.join(node.tune)}</span>")
    if node.extra:
        parts.append(f"<span style='color:{style.active()['ink4']}'>Not in config_default.yaml; read by PRESTUS code.</span>")
    return "<br>".join(parts)


class _Row:
    def __init__(self, node: Node, label: QLabel, container: QWidget, editor: FieldEditor,
                 reset: QToolButton, block: QWidget):
        self.node, self.label, self.container, self.editor, self.reset = \
            node, label, container, editor, reset
        self.block = block


class _Group:
    def __init__(self, node: Node, frame: QFrame, rows: list, groups: list):
        self.node, self.frame, self.rows, self.groups = node, frame, rows, groups


class SectionPage(QFrame):
    """Card with the form for all keys of one top-level section."""

    changed = pyqtSignal()
    advanced_toggled = pyqtSignal()     # the card's own 'advanced' button was used

    def __init__(self, model: ConfigModel, section: str, nodes: list[Node],
                 title: str, parent: QWidget | None = None, subtitle: str = ""):
        super().__init__(parent)
        self.setObjectName("Card")
        self.model = model
        self.section = section
        self.ti = 0
        self.show_advanced = False
        self.goal: str | None = None        # 'Tune for' goal: its keys are shown even when advanced
        self.rows: list[_Row] = []
        self.groups: list[tuple[Node, QFrame]] = []
        self._top_rows: list[_Row] = []
        self._top_groups: list[_Group] = []

        self.vbox = QVBoxLayout(self)
        self.vbox.setContentsMargins(16, 12, 16, 14)
        head = QHBoxLayout()
        self.title = QLabel(title)
        self.title.setObjectName("CardTitle")
        head.addWidget(self.title)
        self.count = QLabel()
        self.count.setObjectName("CardSub")
        head.addWidget(self.count)
        head.addStretch(1)
        self.adv_button = QToolButton()
        self.adv_button.setToolTip("Show or hide this card's advanced settings: the ones most studies "
                                   "leave at their defaults.")
        self.adv_button.clicked.connect(self._toggle_card_advanced)
        head.addWidget(self.adv_button)
        self.toggle = QToolButton(text="Hide")
        self.toggle.setToolTip("Collapse / expand")
        self.toggle.clicked.connect(self._toggle_body)
        head.addWidget(self.toggle)
        self.vbox.addLayout(head)
        if subtitle:
            sub = QLabel(subtitle)
            sub.setObjectName("CardSub")
            sub.setWordWrap(True)
            self.vbox.addWidget(sub)

        self.body = QWidget()
        self.body_layout = QVBoxLayout(self.body)
        self.body_layout.setContentsMargins(0, 6, 0, 0)
        self.vbox.addWidget(self.body)

        self.device_panel: DevicePanel | None = None
        if section == TX:
            self.body_layout.addLayout(self._transducer_bar())
            self.device_panel = DevicePanel(model)
            self.device_panel.changed.connect(self._panel_changed)
            self.body_layout.addWidget(self.device_panel)
        self.extras = QVBoxLayout()
        self.body_layout.addLayout(self.extras)

        grid = self._new_grid()
        self.body_layout.addLayout(grid)
        self._top_rows, self._top_groups = self._fill(grid, self.body_layout, nodes)
        self.refresh()

    # ── construction ────────────────────────────────────────────────────
    @staticmethod
    def _new_grid() -> QGridLayout:
        g = QGridLayout()
        g.setHorizontalSpacing(14)
        g.setVerticalSpacing(8)
        g.setColumnMinimumWidth(0, LABEL_WIDTH)
        g.setColumnStretch(1, 1)
        return g

    def _fill(self, grid: QGridLayout, outer: QVBoxLayout, nodes: list[Node]):
        rows, groups = [], []
        leaves = [n for n in nodes if not n.is_group]
        for node in leaves:
            rows.append(self._add_row(grid, node))
        for node in [n for n in nodes if n.is_group]:
            frame = QFrame()
            frame.setObjectName("SubGroup")
            frame.setToolTip(_tooltip(node))
            v = QVBoxLayout(frame)
            v.setContentsMargins(12, 8, 12, 10)
            t = QLabel(labels.label(node))
            t.setObjectName("SubGroupTitle")
            t.setToolTip(_tooltip(node))
            v.addWidget(t)
            inner = self._new_grid()
            inner.setColumnMinimumWidth(0, LABEL_WIDTH - 14)
            v.addLayout(inner)
            sub_rows, sub_groups = self._fill(inner, v, node.children)
            outer.addWidget(frame)
            self.groups.append((node, frame))
            groups.append(_Group(node, frame, sub_rows, sub_groups))
        return rows, groups

    def _add_row(self, grid: QGridLayout, node: Node) -> _Row:
        tip = _tooltip(node)
        block = QWidget()
        bl = QVBoxLayout(block)
        bl.setContentsMargins(0, 0, 0, 0)
        bl.setSpacing(1)
        label = QLabel(labels.label(node))
        label.setObjectName("FieldLabel")
        label.setToolTip(tip)
        label.setWordWrap(True)
        bl.addWidget(label)

        editor = make_editor(node)
        editor.setToolTip(tip)
        editor.setEnabled(not node.derived)
        container = QWidget()
        h = QHBoxLayout(container)
        h.setContentsMargins(0, 0, 0, 0)
        h.addWidget(editor, 1)
        if node.unit:
            u = QLabel(node.unit)
            u.setObjectName("CardSub")
            h.addWidget(u)
        helpb = QToolButton(text="?")
        helpb.setToolTip("What is this?")
        helpb.clicked.connect(lambda _=False, t=tip: QToolTip.showText(QCursor.pos(), t, helpb))
        h.addWidget(helpb)
        reset = QToolButton(text="↺")
        reset.setToolTip("Reset to default")
        h.addWidget(reset)
        container.setMaximumWidth(FIELD_MAX_WIDTH.get(node.kind, FIELD_MAX_DEFAULT))
        r = grid.rowCount()
        grid.addWidget(block, r, 0, Qt.AlignmentFlag.AlignTop)
        grid.addWidget(container, r, 1, Qt.AlignmentFlag.AlignTop)
        row = _Row(node, label, container, editor, reset, block)
        editor.edited.connect(lambda value, rr=row: self._edited(rr, value))
        reset.clicked.connect(lambda _=False, rr=row: self._reset(rr))
        self.rows.append(row)
        return row

    def add_extra(self, widget: QWidget) -> None:
        """A widget shown above the fields (diagram, panel, …)."""
        self.extras.addWidget(widget)

    def _toggle_body(self) -> None:
        self.body.setVisible(not self.body.isVisible())
        self.toggle.setText("Hide" if self.body.isVisible() else "Show")

    def _transducer_bar(self) -> QHBoxLayout:
        bar = QHBoxLayout()
        self.tx_combo = QComboBox()
        self.tx_combo.currentIndexChanged.connect(self._select_transducer)
        add = QPushButton("+")
        add.setToolTip("Add a transducer")
        add.setProperty("compact", True)
        add.clicked.connect(self._add_transducer)
        self.tx_remove = QPushButton("−")
        self.tx_remove.setToolTip("Remove this transducer")
        self.tx_remove.setProperty("compact", True)
        self.tx_remove.clicked.connect(self._remove_transducer)
        self.tx_note = QLabel()
        self.tx_note.setObjectName("CardSub")
        bar.addWidget(QLabel("Editing:"))
        bar.addWidget(self.tx_combo)
        bar.addWidget(add)
        bar.addWidget(self.tx_remove)
        bar.addWidget(self.tx_note, 1)
        return bar

    # ── transducer handling ─────────────────────────────────────────────
    def _update_tx_combo(self) -> None:
        self.tx_combo.blockSignals(True)
        self.tx_combo.clear()
        self.tx_combo.addItems([f"Transducer {i + 1}" for i in range(self.model.n_transducers)])
        self.ti = min(self.ti, self.model.n_transducers - 1)
        self.tx_combo.setCurrentIndex(self.ti)
        self.tx_combo.blockSignals(False)
        self.tx_remove.setEnabled(self.model.n_transducers > 1)
        self.tx_note.setText("Unset values fall back to transducer 1 (PRESTUS merge rule)"
                             if self.ti > 0 else "")

    def _select_transducer(self, index: int) -> None:
        if index >= 0:
            self.ti = index
            self.refresh()
            self.changed.emit()

    def _add_transducer(self) -> None:
        self.ti = self.model.add_transducer()
        self.refresh()
        self.changed.emit()

    def _remove_transducer(self) -> None:
        self.model.remove_transducer(self.ti)
        self.ti = max(0, self.ti - 1)
        self.refresh()
        self.changed.emit()

    def _panel_changed(self) -> None:
        self.refresh()
        self.changed.emit()

    # ── editing ─────────────────────────────────────────────────────────
    def _edited(self, row: _Row, value) -> None:
        self.model.set(row.node.path, value, self.ti)
        self.refresh()
        self.changed.emit()

    def _reset(self, row: _Row) -> None:
        self.model.reset(row.node.path, self.ti)
        self.refresh()
        self.changed.emit()

    # ── visibility ──────────────────────────────────────────────────────
    def set_show_advanced(self, show: bool) -> None:
        self.show_advanced = show
        self.refresh()

    def _toggle_card_advanced(self) -> None:
        self.set_show_advanced(not self.show_advanced)
        self.advanced_toggled.emit()

    def set_goal(self, goal: str | None) -> None:
        self.goal = goal or None
        self.refresh()

    def _matches(self, node: Node) -> bool:
        return self.goal is not None and self.goal in node.tune

    def _ti_for(self, node: Node) -> int:
        return self.ti if node.path.startswith(TX + ".") or node.path == TX else 0

    def _row_shown(self, row: _Row) -> bool:
        path, ti = row.node.path, self._ti_for(row.node)
        if not self.model.is_visible(path, ti):
            return False
        return self.show_advanced or not row.node.advanced or self.model.is_overridden(path, ti) \
            or self._matches(row.node)

    def _group_shown(self, g: _Group) -> bool:
        ti = self._ti_for(g.node)
        if not self.model.is_visible(g.node.path, ti):
            return False
        any_child = any(self._row_shown(r) for r in g.rows) or any(self._group_shown(s) for s in g.groups)
        if not any_child:
            return False
        return self.show_advanced or not g.node.advanced or self._has_override(g) or self._has_match(g)

    def _has_match(self, g: _Group) -> bool:
        return any(self._matches(r.node) for r in g.rows) or any(self._has_match(s) for s in g.groups)

    def matches_goal(self) -> bool:
        """The card holds a key tagged with the current 'Tune for' goal (in the current modes)."""
        return any(self._matches(r.node) and self.model.is_visible(r.node.path, self._ti_for(r.node))
                   for r in self.rows)

    def _has_override(self, g: _Group) -> bool:
        ti = self._ti_for(g.node)
        return any(self.model.is_overridden(r.node.path, ti) for r in g.rows) or \
            any(self._has_override(s) for s in g.groups)

    def has_overrides(self) -> bool:
        return any(self.model.is_overridden(r.node.path, self._ti_for(r.node)) for r in self.rows)

    def _advanced_rows(self) -> list[_Row]:
        """Advanced fields that apply in the current modes and hold no value (the ones the
        'advanced' switches show and hide)."""
        return [r for r in self.rows if r.node.advanced
                and self.model.is_visible(r.node.path, self._ti_for(r.node))
                and not self.model.is_overridden(r.node.path, self._ti_for(r.node))]

    def hidden_advanced(self) -> int:
        """Fields hidden only because they are advanced."""
        if self.show_advanced:
            return 0
        return sum(1 for r in self._advanced_rows() if not self._matches(r.node))

    def _tune_tags(self, node: Node) -> str:
        """' · speed, accuracy' after the label; the 'Tune for' goal highlighted."""
        if not node.tune or not (self.show_advanced or self.goal):
            return ""
        schema = self.model.schema
        tags = [f"<b style='color:{style.active()['accent']}'>{schema.tune_label(g).lower()}</b>" if g == self.goal
                else schema.tune_label(g).lower() for g in node.tune]
        return f" <span style='color:{style.active()['ink4']}'>· {', '.join(tags)}</span>"

    def refresh(self) -> None:
        for i in range(self.extras.count()):
            w = self.extras.itemAt(i).widget()
            if w is not None and hasattr(w, "refresh"):
                w.refresh()
        if self.section == TX:
            self._update_tx_combo()
            if self.device_panel is not None:
                self.device_panel.refresh(self.ti)

        def walk(rows: list[_Row], groups: list[_Group]) -> None:
            for g in groups:
                g.frame.setVisible(self._group_shown(g))
                walk(g.rows, g.groups)
            for row in rows:
                visible = self._row_shown(row)
                row.block.setVisible(visible)
                row.container.setVisible(visible)
                if not visible:
                    continue
                path, ti = row.node.path, self._ti_for(row.node)
                value = self.model.get(path, ti)
                row.editor.set_value(value, self.model.baseline(path, ti))
                row.editor.set_error(check_value(row.node, value))
                if path == "placement.mode" and hasattr(row.editor, "set_item_status"):
                    row.editor.set_item_status(placement_mode_status(self.model))
                name = labels.label(row.node)
                tags = self._tune_tags(row.node)
                if self.model.is_modified(path, ti):
                    row.label.setText(f"<b>{name}</b> <span style='color:#1f6fb2'>●</span>{tags}")
                elif self.model.from_library(path, ti):
                    row.label.setText(f"{name} <span style='color:#9aa7b4'>· from library</span>{tags}")
                else:
                    row.label.setText(f"{name}{tags}")
                row.reset.setVisible(self.model.is_overridden(path, ti))
        walk(self._top_rows, self._top_groups)
        n_set = sum(1 for r in self.rows if self.model.is_overridden(r.node.path, self._ti_for(r.node)))
        self.count.setText(f"· {n_set} set" if n_set else "")
        n_hidden = self.hidden_advanced()
        self.adv_button.setVisible(n_hidden > 0 or (self.show_advanced and bool(self._advanced_rows())))
        self.adv_button.setText("Hide advanced" if self.show_advanced else f"+ {n_hidden} advanced")

    def focus_key(self, path: str, ti: int | None = None) -> None:
        """Scroll to (and focus) the row for a key; reveals advanced fields if needed."""
        if ti is not None and self.section == TX and ti != self.ti:
            self.ti = ti
            self.refresh()
        for row in self.rows:
            if row.node.path != path:
                continue
            if not row.container.isVisible() and row.node.advanced and not self.show_advanced:
                self.set_show_advanced(True)
            if not self.body.isVisible():
                self._toggle_body()
            w = self.parentWidget()
            while w is not None and not isinstance(w, QScrollArea):
                w = w.parentWidget()
            if isinstance(w, QScrollArea):
                w.ensureWidgetVisible(row.container)
            row.editor.setFocus(Qt.FocusReason.OtherFocusReason)
            return
