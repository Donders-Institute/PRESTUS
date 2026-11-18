"""The form for one config section, shown as a card on a tab of the main window."""

from __future__ import annotations

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (QComboBox, QFrame, QGridLayout, QHBoxLayout, QLabel, QMessageBox, QPushButton,
                             QScrollArea, QSizePolicy, QToolButton, QVBoxLayout, QWidget)

from .. import docs, labels
from ..checks import check_value, placement_mode_status
from ..model import TX, ConfigModel, is_empty
from ..schema import Node
from . import icons, style
from .fields import CheckEditor, FieldEditor, display_text, make_editor
from .transducer_panel import DevicePanel

LABEL_WIDTH = 290       # label column of a wide card; narrower cards scale it down (label_width)
LABEL_MIN_WIDTH = 170
LABEL_SHARE = 0.3       # of the card width

def _label_html(node: Node) -> str:
    """The field label with its unit after it: 'Time step within a train (s)'."""
    import html
    name = labels.rich_label(node)
    return f"{name} ({html.escape(node.unit, quote=False)})" if node.unit else name
CHECK_FIRST = {"modules"}   # sections whose checkboxes sit before the label (a list of stages to tick)
# an icon per pipeline stage, between its checkbox and its label (icons.names())
STAGE_ICONS = {
    "modules.run_transducer_placement": "target", "modules.run_grid_setup": "grid",
    "modules.run_medium_setup": "layers", "modules.run_source_setup": "tx_annular",
    "modules.run_water_baseline": "baseline_beam", "modules.run_acoustic_sims": "wave",
    "modules.run_acoustic_analysis": "chart", "modules.run_heating_sims": "thermometer",
    "modules.run_thermal_analysis": "histogram", "modules.run_nifti_creation": "file",
    "modules.run_posthoc_water_sims": "medium_water", "modules.generate_report": "report",
    "modules.run_targeting_analysis": "aim_corrected", "modules.generate_group_report": "people",
    "modules.segmentation_only": "segmentation", "simulation.debug": "wrench",
}
CHECK_SLOT = 24             # width of that checkbox slot; other rows of the section keep it empty
PANEL_KEYS = {f"{TX}.serial", f"{TX}.combo.ds_serial"}   # the Device panel edits these


def label_width(card_width: int) -> int:
    """Label column for a card `card_width` wide: LABEL_SHARE of it, within
    [LABEL_MIN_WIDTH, LABEL_WIDTH]; labels wrap, so a narrow window keeps room for the inputs."""
    return max(LABEL_MIN_WIDTH, min(LABEL_WIDTH, int(card_width * LABEL_SHARE)))


def mark_label_grid(grid: QGridLayout, inset: int = 0) -> None:
    """Let the card resize column 0 of `grid` with its label column; `inset` is how much
    narrower it is (a grid inside an indented frame), so its inputs start where the card's do."""
    grid.setProperty("label_inset", inset)
    grid.setColumnMinimumWidth(0, LABEL_WIDTH - inset)


def _tooltip(node: Node) -> str:
    """Full help for hover: description, choices, key, default (ui/tooltip.py)."""
    from .tooltip import html_for
    return html_for(node, labels.label(node), None if node.is_group else display_text(node.default))


def _help_popup(node: Node, schema, anchor: QWidget) -> None:
    """The ? button: the same help in a popup whose doc links can be clicked, with the key's
    documentation page (docs.key_doc)."""
    from .tooltip import html_for, show_help
    html = html_for(node, labels.label(node), None if node.is_group else display_text(node.default),
                    links=True, doc=docs.key_doc(node.path, schema))
    show_help(html, anchor)


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
    transducer_selected = pyqtSignal(int)   # the transducer this card now edits
    advanced_toggled = pyqtSignal()     # the card's own 'advanced' button was used

    def __init__(self, model: ConfigModel, section: str, nodes: list[Node],
                 title: str, parent: QWidget | None = None, subtitle: str = "", locked: bool = False,
                 icon: str | None = None):
        super().__init__(parent)
        self._label_w = LABEL_WIDTH
        self.setObjectName("Card")
        self.model = model
        self.section = section
        self.ti = 0
        self.show_advanced = False
        self.folded: set[str] = set()       # set_fold: rows another control on the card edits
        self.moved: dict = {}               # move_rows: path -> reveal() of the widget now holding the row
        self.fold_open = False
        self.goal: str | None = None        # 'Tune for' goal: its keys are shown even when advanced
        self.locked = locked                # read-only until unlocked (overlay 'locked')
        self.lockable = locked
        self.lock_note = ""                  # overlay 'lock_note': why the card is locked
        self.table = None                   # set_table: the section's own keys edited in one widget
        self.leader: "SectionPage | None" = None    # follow_transducer: card whose selector this follows
        self.collapse_if = None             # set_auto_collapse: callable -> start collapsed when True
        self.collapse_note = ""
        self._user_toggled = False
        self.rows: list[_Row] = []
        self.groups: list[tuple[Node, QFrame]] = []
        self._top_rows: list[_Row] = []
        self._top_groups: list[_Group] = []

        self.vbox = QVBoxLayout(self)
        self.vbox.setContentsMargins(16, 12, 16, 14)
        head = QHBoxLayout()
        self.icon_name = icon               # overlay 'sections: <name>: icon': a glyph left of the title
        self.icon = QLabel()
        self.icon.setVisible(bool(icon))
        if icon:
            self.icon.setFixedSize(22, 22)
        self.retint()
        head.addWidget(self.icon)
        self.title = QLabel(title)
        self.title.setObjectName("CardTitle")
        head.addWidget(self.title)
        self.count = QLabel()
        self.count.setObjectName("CardSub")
        head.addWidget(self.count)
        head.addStretch(1)
        self.doc_ref = docs.section_doc(section)      # overlay 'sections: <name>: doc'
        self.doc_button = QToolButton()
        self.doc_button.setText("↗")
        self.doc_button.setToolTip("Open documentation")
        self.doc_button.setAccessibleName(f"Open documentation: {title}")
        self.doc_button.clicked.connect(lambda _=False: docs.open_doc(self.doc_ref))
        self.doc_button.setVisible(bool(self.doc_ref))
        head.addWidget(self.doc_button)
        self.adv_button = QToolButton()
        self.adv_button.setIcon(icons.icon("advanced"))
        self.adv_button.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.adv_button.setToolTip("Show or hide this card's advanced settings: the ones most studies "
                                   "leave at their defaults.")
        self.adv_button.clicked.connect(self._toggle_card_advanced)
        head.addWidget(self.adv_button)
        self.lock_button = QToolButton()
        self.lock_button.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.lock_button.clicked.connect(self._toggle_lock)
        self.lock_button.setVisible(locked)
        head.addWidget(self.lock_button)
        self.toggle = QToolButton()
        self.toggle.setIcon(icons.icon("collapse"))
        self.toggle.setToolTip("Collapse this card")
        self.toggle.clicked.connect(self._toggle_body)
        head.addWidget(self.toggle)
        self.vbox.addLayout(head)
        self.headless = False               # set_headless: the tab already names this card
        self.subtitle = None
        if subtitle:                        # what the card is for: on hover, not as a line of text
            self.subtitle_text = subtitle
            for w in (self.title, self.icon):
                if w is not None:
                    w.setToolTip(subtitle)

        self.body = QWidget()
        self.body_layout = QVBoxLayout(self.body)
        self.body_layout.setContentsMargins(0, 6, 0, 0)
        self.vbox.addWidget(self.body)

        self.device_panel: DevicePanel | None = None
        # Per-transducer keys shown on another card (positions on Placement) get the same
        # transducer selector, so each transducer of a study stays editable there.
        self.has_tx = section == TX or any(n.path.startswith(TX + ".") for n in nodes)
        if self.has_tx:
            self.tx_bar = QWidget()
            self.tx_bar.setLayout(self._transducer_bar())
            self.body_layout.addWidget(self.tx_bar)
        if section == TX:
            self.device_panel = DevicePanel(model)
            self.device_panel.changed.connect(self._panel_changed)
            self.body_layout.addWidget(self.device_panel)
        self.extras = QVBoxLayout()
        self.body_layout.addLayout(self.extras)

        grid = self._new_grid()
        self.body_layout.addLayout(grid)
        self._top_rows, self._top_groups = self._fill(grid, self.body_layout, nodes)
        self.extras_below = QVBoxLayout()
        self.body_layout.addLayout(self.extras_below)
        self.refresh()

    def resizeEvent(self, event) -> None:                # noqa: N802
        super().resizeEvent(event)
        self._fit_label_column()

    def _fit_label_column(self) -> None:
        """Scale the label column with the card (label_width): its grids, the label blocks and
        panels that mark their grids (mark_label_grid)."""
        w = label_width(self.width())
        if w == self._label_w:
            return
        self._label_w = w
        for g in self.findChildren(QGridLayout):
            inset = g.property("label_inset")
            if inset is not None:
                g.setColumnMinimumWidth(0, w - int(inset))
        for row in self.rows:
            if not getattr(row, "lead", False):
                row.block.setMaximumWidth(w + 40)

    # ── construction ────────────────────────────────────────────────────
    @staticmethod
    def _new_grid() -> QGridLayout:
        g = QGridLayout()
        g.setHorizontalSpacing(14)
        g.setVerticalSpacing(8)
        mark_label_grid(g)
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
            t = QLabel(labels.rich_label(node))
            t.setTextFormat(Qt.TextFormat.RichText)
            t.setObjectName("SubGroupTitle")
            t.setToolTip(_tooltip(node))
            v.addWidget(t)
            inner = self._new_grid()
            mark_label_grid(inner, 14)
            v.addLayout(inner)
            # children another card shows (overlay 'group') stay there
            children = [c for c in node.children if not c.group or c.group == self.section]
            sub_rows, sub_groups = self._fill(inner, v, children)
            outer.addWidget(frame)
            self.groups.append((node, frame))
            groups.append(_Group(node, frame, sub_rows, sub_groups))
        return rows, groups

    def _add_row(self, grid: QGridLayout, node: Node) -> _Row:
        tip = _tooltip(node)
        block = QWidget()
        bl = QVBoxLayout(block)
        bl.setContentsMargins(0, 5, 0, 0)      # the label's first line on the input's text line
        bl.setSpacing(1)
        label = QLabel(_label_html(node))
        label.setTextFormat(Qt.TextFormat.RichText)       # escaped; I<sub>SPPA</sub> (labels.rich)
        label.setObjectName("FieldLabel")
        label.setToolTip(tip)
        label.setWordWrap(True)

        editor = make_editor(node)
        editor.setToolTip(tip)
        editor.setEnabled(not node.derived)
        # Stages: the checkbox leads its label; the other rows keep the slot empty so labels align
        check_first = self.section in CHECK_FIRST
        lead = check_first and isinstance(editor, CheckEditor)
        if check_first:
            line = QHBoxLayout()
            line.setContentsMargins(0, 0, 0, 0)
            line.setSpacing(4)
            slot = editor if lead else QWidget()
            slot.setFixedWidth(CHECK_SLOT)
            line.addWidget(slot, 0, Qt.AlignmentFlag.AlignTop)
            pic = QLabel()                      # the stage's icon (empty slot keeps labels aligned)
            pic.setFixedSize(18, 18)
            if node.path in STAGE_ICONS:
                pic.setPixmap(icons.icon(STAGE_ICONS[node.path], style.active()["ink2"]).pixmap(16, 16))
            line.addWidget(pic, 0, Qt.AlignmentFlag.AlignTop)
            line.addWidget(label, 1)
            bl.addLayout(line)
        else:
            bl.addWidget(label)
        container = QWidget()
        container.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        h = QHBoxLayout(container)
        h.setContentsMargins(0, 0, 0, 0)
        if lead:
            h.addStretch(1)                 # unit, ? and ↺ stay right-aligned like the other rows
        else:
            h.addWidget(editor, 1)
        # every row has the same columns: input | ? | ↺, each top-aligned; the unit is in the label
        helpb = QToolButton()
        helpb.setIcon(icons.icon("help", style.active()["ink3"]))
        helpb.setToolTip("What is this?")
        helpb.clicked.connect(lambda _=False, n=node: _help_popup(n, self.model.schema, helpb))
        h.addWidget(helpb, 0, Qt.AlignmentFlag.AlignTop)
        reset = QToolButton()
        reset.setIcon(icons.icon("reset", style.active()["ink3"]))
        reset.setToolTip("Reset to default")
        sp = reset.sizePolicy()
        sp.setRetainSizeWhenHidden(True)        # help icons stay in one column, set or not
        reset.setSizePolicy(sp)
        h.addWidget(reset, 0, Qt.AlignmentFlag.AlignTop)
        # The label column keeps its width; the input fills the rest of the card, so a wide
        # window gives longer fields rather than empty space beside them.
        if lead:                                # a stage: no input beside it, so one line
            label.setWordWrap(False)
        else:
            block.setMaximumWidth(self._label_w + 40)
        r = grid.rowCount()
        grid.addWidget(block, r, 0, Qt.AlignmentFlag.AlignTop)
        grid.addWidget(container, r, 1, Qt.AlignmentFlag.AlignTop)   # full width: the editors expand
        row = _Row(node, label, container, editor, reset, block)
        row.lead = lead
        editor.edited.connect(lambda value, rr=row: self._edited(rr, value))
        reset.clicked.connect(lambda _=False, rr=row: self._reset(rr))
        self.rows.append(row)
        return row

    def add_extra(self, widget: QWidget, below: bool = False) -> None:
        """A widget shown above the fields (diagram, panel, …), or below them."""
        (self.extras_below if below else self.extras).addWidget(widget)

    def move_rows(self, paths, grid: QGridLayout, reveal=None) -> None:
        """Show the rows of `paths` in `grid` on another widget (e.g. an advanced fold on a
        companion card); this page still fills, checks and resets them. `reveal()` opens that
        widget before focus_key scrolls to one of them."""
        for row in self.rows:
            if row.node.path in paths:
                r = grid.rowCount()
                grid.addWidget(row.block, r, 0, Qt.AlignmentFlag.AlignTop)
                grid.addWidget(row.container, r, 1, Qt.AlignmentFlag.AlignTop)
                self.moved[row.node.path] = reveal
        self.refresh()

    def set_table(self, widget: QWidget) -> None:
        """Edit the section's own keys in `widget` (below the fields) instead of one row each;
        keys other sections moved here keep their rows. The widget gets set_locked(bool),
        refresh(), focus(path) and a `changed` signal."""
        self.table = widget
        self.extras_below.addWidget(widget)
        widget.changed.connect(self._panel_changed)
        self.refresh()

    def _tabulated(self, node: Node) -> bool:
        return self.table is not None and node.section == self.section

    def follow_transducer(self, leader: "SectionPage") -> None:
        """Edit the transducer `leader` selects, without a selector of its own (cards on one
        tab share the Transducer card's selector)."""
        self.leader = leader
        self.tx_bar.setVisible(False)
        self.refresh()

    def set_fold(self, paths, open_: bool | None = None) -> None:
        """Hide the rows of `paths` while the fold is closed: values another control on the card
        edits (e.g. the timing diagram). `open_` opens or closes the fold."""
        self.folded = set(paths)
        if open_ is not None:
            self.fold_open = bool(open_)
        self.refresh()

    def set_auto_collapse(self, predicate, note: str = "") -> None:
        """Start collapsed while predicate() holds (e.g. geometry from the equipment library),
        until the card is opened or closed by hand; `note` follows the title meanwhile."""
        if self.headless:                   # no title to collapse to
            return
        self.collapse_if, self.collapse_note = predicate, note
        self.refresh()

    def set_headless(self) -> None:
        """Drop the icon, title, count, subtitle and collapse button when the tab or sub-tab
        label and intro already name the card (it is the page's only card, or has the tab's
        name). The documentation, advanced and lock buttons stay."""
        self.headless = True
        for w in (self.icon, self.title, self.count, self.subtitle, self.toggle):
            if w is not None:
                w.hide()
        self.collapse_if = None
        self._set_body(True)

    def _toggle_body(self) -> None:
        self._user_toggled = True
        self._set_body(self.body.isHidden())         # isHidden: also right while the card is off-screen

    def _set_body(self, shown: bool) -> None:
        self.body.setVisible(shown)
        self.toggle.setIcon(icons.icon("collapse" if shown else "expand"))
        self.toggle.setToolTip("Collapse this card" if shown else "Expand this card")

    def _transducer_bar(self) -> QHBoxLayout:
        bar = QHBoxLayout()
        bar.setContentsMargins(0, 0, 0, 0)         # flush with the field labels
        self.tx_combo = QComboBox()
        self.tx_combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)   # shows the serial
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
        for i in range(self.model.n_transducers):
            self.tx_combo.addItem(self.model.transducer_title(i))
            self.tx_combo.setItemData(i, self.model.transducer_tooltip(i), Qt.ItemDataRole.ToolTipRole)
        self.ti = min(self.ti, self.model.n_transducers - 1)
        self.tx_combo.setCurrentIndex(self.ti)
        self.tx_combo.setToolTip(self.model.transducer_tooltip(self.ti))
        self.tx_combo.blockSignals(False)
        self.tx_remove.setEnabled(self.model.n_transducers > 1)
        self.tx_note.setText("Unset values fall back to transducer 1 (PRESTUS merge rule)"
                             if self.ti > 0 else "")

    def _select_transducer(self, index: int) -> None:
        if index >= 0:
            self.ti = index
            self.refresh()
            self.transducer_selected.emit(self.ti)
            self.changed.emit()

    def _add_transducer(self) -> None:
        self.ti = self.model.add_transducer()
        self.refresh()
        self.transducer_selected.emit(self.ti)
        self.changed.emit()

    def _remove_transducer(self) -> None:
        self.model.remove_transducer(self.ti)
        self.ti = max(0, self.ti - 1)
        self.refresh()
        self.transducer_selected.emit(self.ti)
        self.changed.emit()

    def show_transducer(self, ti: int) -> None:
        """Edit transducer ti (the window's current transducer), without announcing it."""
        if self.has_tx and self.leader is None and 0 <= ti < self.model.n_transducers and ti != self.ti:
            self.ti = ti
            self.refresh()

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

    def _toggle_lock(self) -> None:
        if self.locked:
            note = self.lock_note or (
                "These are the acoustic medium property values that simulated pressures and "
                "temperatures depend on. Change them only for a documented reason (e.g. a sensitivity "
                "analysis) and always report the values you used.")
            answer = QMessageBox.warning(
                self, f"Unlock {self.title.text()}", f"{note}\n\nUnlock for editing?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
                QMessageBox.StandardButton.Cancel)
            if answer != QMessageBox.StandardButton.Yes:
                return
        self.locked = not self.locked
        self.refresh()

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

    def retint(self) -> None:
        """Redraw the card icon in the active palette (after a theme change)."""
        if self.icon_name:
            self.icon.setPixmap(icons.icon(self.icon_name, style.active()["accent"], 20).pixmap(20, 20))

    def _row_shown(self, row: _Row) -> bool:
        if self._tabulated(row.node):
            return False
        path, ti = row.node.path, self._ti_for(row.node)
        if row.node.derived and is_empty(self.model.get(path, ti)):
            return False                    # set by the pipelines (overlay 'derived'); shown once it holds a value
        if self.device_panel is not None and path in PANEL_KEYS:
            return False                    # edited in the Device panel above the fields
        if not self.model.is_visible(path, ti):
            return False
        if path in self.folded and not self.fold_open and not self._matches(row.node) \
                and not check_value(row.node, self.model.get(path, ti)):
            return False                    # folded away; an invalid value stays visible
        return self.show_advanced or not row.node.advanced or self.model.is_overridden(path, ti) \
            or self._matches(row.node)

    def _group_shown(self, g: _Group) -> bool:
        if self._tabulated(g.node):
            return False
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
        if not node.tune or not self.goal:              # only while a 'Tune for' goal is set
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
        for i in range(self.extras_below.count()):
            w = self.extras_below.itemAt(i).widget()
            if w is not None and w is not self.table and hasattr(w, "refresh"):
                w.refresh()
        if self.table is not None:
            self.table.set_locked(self.locked)
            self.table.refresh()
        if self.leader is not None:
            self.ti = min(self.leader.ti, self.model.n_transducers - 1)
        if self.has_tx:
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
                name = _label_html(row.node)
                tags = self._tune_tags(row.node)
                if self.model.is_modified(path, ti):
                    row.label.setText(f"<b>{name}</b> <span style='color:{style.active()['accent']}'>●</span>{tags}")
                elif self.model.from_library(path, ti):
                    row.label.setText(f"{name} <span style='color:{style.active()['ink4']}'>· from library</span>{tags}")
                else:
                    row.label.setText(f"{name}{tags}")
                row.reset.setVisible(self.model.is_overridden(path, ti))
                # the lock covers the section's own keys, not those moved here from other sections
                locked = self.locked and row.node.section == self.section
                row.editor.setEnabled(not row.node.derived and not locked)
                row.reset.setEnabled(not locked)
        walk(self._top_rows, self._top_groups)
        n_set = sum(1 for r in self.rows if self.model.is_overridden(r.node.path, self._ti_for(r.node)))
        self.count.setText(f"· {n_set} set" if n_set else "")
        if self.collapse_if is not None:
            collapse = bool(self.collapse_if())
            if not self._user_toggled:
                self._set_body(not collapse)
            if collapse and self.collapse_note and self.body.isHidden():
                self.count.setText(f"{self.count.text()}  · {self.collapse_note}".strip())
        n_hidden = self.hidden_advanced()
        self.adv_button.setVisible(n_hidden > 0 or (self.show_advanced and bool(self._advanced_rows()))
                                   or (self.lockable and self.locked))
        # the advanced icon with the number of hidden settings (greyed); highlighted while shown
        c = style.active()
        self.adv_button.setIcon(icons.icon("advanced", c["accent"] if self.show_advanced else None))
        self.adv_button.setText("" if self.show_advanced or not n_hidden else str(n_hidden))
        self.adv_button.setStyleSheet(f"QToolButton {{ color: {c['ink4']}; }}")
        self.adv_button.setToolTip(
            "Hide this card's advanced settings" if self.show_advanced else
            f"Show this card's {n_hidden} advanced setting{'s' if n_hidden != 1 else ''}: the ones most "
            "studies leave at their defaults" if n_hidden else "Show this card's advanced settings")
        # unlocking is an advanced action; a card unlocked once keeps its 'Lock' button
        self.lock_button.setVisible(self.lockable and (self.show_advanced or not self.locked))
        self.lock_button.setText("Unlock" if self.locked else "Lock")
        self.lock_button.setIcon(icons.icon("lock" if self.locked else "unlock"))
        self.lock_button.setToolTip("Read-only: defaults most studies should not change. Click to edit."
                                    if self.locked else "Make this card read-only again.")

    def focus_key(self, path: str, ti: int | None = None) -> None:
        """Scroll to (and focus) the row for a key; reveals advanced fields if needed."""
        if ti is not None and self.has_tx and ti != self.ti:
            self.ti = ti
            self.refresh()
        if self.table is not None and path.startswith(self.section + ".") and self.table.focus(path):
            if self.body.isHidden():
                self._set_body(True)
            w = self.parentWidget()
            while w is not None and not isinstance(w, QScrollArea):
                w = w.parentWidget()
            if isinstance(w, QScrollArea):
                w.ensureWidgetVisible(self.table)
            return
        if self.device_panel is not None and path in PANEL_KEYS:
            if self.body.isHidden():
                self._set_body(True)
            w = self.parentWidget()
            while w is not None and not isinstance(w, QScrollArea):
                w = w.parentWidget()
            if isinstance(w, QScrollArea):
                w.ensureWidgetVisible(self.device_panel)
            return
        for row in self.rows:
            if row.node.path != path:
                continue
            if not row.container.isVisible() and row.node.advanced and not self.show_advanced:
                self.set_show_advanced(True)
            if self.moved.get(path) is not None:
                self.moved[path]()
            if self.body.isHidden():
                self._set_body(True)
            w = self.parentWidget()
            while w is not None and not isinstance(w, QScrollArea):
                w = w.parentWidget()
            if isinstance(w, QScrollArea):
                w.ensureWidgetVisible(row.container)
            row.editor.setFocus(Qt.FocusReason.OtherFocusReason)
            return
