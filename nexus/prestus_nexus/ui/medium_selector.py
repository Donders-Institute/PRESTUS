"""A tab's header row for the choice that decides what the tab shows: simulation.medium on
the Medium tab (water: tissue properties of water only; phantom: tissue layers and all
properties; layered: everything, including the CT-based skull mapping), placement.mode on the Placement
tab (PLACEMENT_NOTE)."""

from __future__ import annotations

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import QComboBox, QFrame, QHBoxLayout, QLabel, QToolButton

from .. import labels
from . import icons, style
from .tooltip import choice_tip

KEY = "simulation.medium"
NOTE = {
    "water": "Free field in water: only the water properties apply.",
    "layered": "Multi-layer head model based on segmentation.",
    "phantom": "Designed test object made of tissue layers.",
}
PLACEMENT_KEY = "placement.mode"
PLACEMENT_NOTE = {
    "manual": "Transducer and focus coordinates set by hand, on the Viewer or as numbers.",
    "localite": "A pose recorded with Localite neuronavigation: ingest it on the Localite sub-tab.",
    "heuristic": "A search from an MNI target over scalp positions: compare and choose candidates on the "
                 "Heuristic sub-tab.",
    "mni": "MNI coordinates mapped to this subject and snapped to the scalp.",
    "plantus": "PlanTUS scores every skin position: compute and choose on the PlanTUS sub-tab.",
    "tus_entry": "TUS_entry picks the flattest scalp entry within the distance window: run it on the "
                 "TUS entry sub-tab.",
}


# placement methods whose sub-tab lists what is still to set up (setup_note): the selector only
# points there, so the header row stays short; other methods keep their (short) reason inline
PLACEMENT_SUBTAB = {"plantus": "PlanTUS", "tus_entry": "TUS entry"}


def _what(st) -> tuple[str, str]:
    colour = style.active()["red" if st.state == "unavailable" else "amber"]
    return colour, "Cannot run yet" if st.state == "unavailable" else "Still to set up"


def setup_note(st) -> str:
    """What a placement method still needs (checks.ModeStatus), one item per line; '' when ready."""
    if st is None or st.state == "ready" or not st.reason:
        return ""
    colour, what = _what(st)
    items = "".join(f"<li>{r[:1].upper()}{r[1:]}</li>" for r in st.reason.split("; ") if r)
    return f"<span style='color:{colour}'><b>{what}</b></span><ul style='margin:2px 0'>{items}</ul>"


class MediumSelector(QFrame):
    changed = pyqtSignal()

    def __init__(self, model, parent=None, key: str = KEY, title: str = "Medium",
                 notes: dict[str, str] | None = None, status=None):
        super().__init__(parent)
        self.setObjectName("Card")
        self.model = model
        self.key, self.title = key, title
        self.notes = NOTE if notes is None else notes
        self.status = status            # model -> {choice: checks.ModeStatus}: availability marks
        node = model.schema.get(key)
        h = QHBoxLayout(self)
        h.setContentsMargins(14, 8, 14, 8)
        self.label = QLabel(title)
        self.label.setObjectName("CardTitle")
        h.addWidget(self.label)
        self.combo = QComboBox()
        for choice in (node.choices if node and node.choices else ["water", "layered", "phantom"]):
            ic = icons.for_choice(node, choice)
            if ic is not None:
                self.combo.addItem(ic, labels.choice(node, choice), choice)
            else:
                self.combo.addItem(labels.choice(node, choice), choice)
            self.combo.setItemData(self.combo.count() - 1, choice_tip(node, choice), Qt.ItemDataRole.ToolTipRole)
        if node is not None:
            from .fields import display_text
            from .tooltip import html_for
            self.combo.setToolTip(html_for(node, title, display_text(node.default)))
        self._help_tip = self.combo.toolTip()       # the key's help; the choice's note goes above it
        self.combo.activated.connect(self._picked)
        h.addWidget(self.combo)
        self.reset = QToolButton(text="↺")
        self.reset.setToolTip("Reset to default")
        self.reset.clicked.connect(self._reset)
        h.addWidget(self.reset)
        self.note = QLabel()
        self.note.setObjectName("CardSub")
        self.note.setWordWrap(True)
        h.addWidget(self.note, 1)
        self.refresh()

    def _picked(self, index: int) -> None:
        self.model.set(self.key, self.combo.itemData(index))
        self.refresh()
        self.changed.emit()

    def _reset(self) -> None:
        self.model.reset(self.key)
        self.refresh()
        self.changed.emit()

    def refresh(self) -> None:
        value = self.model.get(self.key)
        i = self.combo.findData(value)
        self.combo.blockSignals(True)
        if i < 0 and value not in (None, ""):
            self.combo.addItem(str(value), value)
            i = self.combo.count() - 1
        self.combo.setCurrentIndex(max(i, 0))
        self.combo.blockSignals(False)
        modified = self.model.is_modified(self.key)
        self.label.setText(f"{self.title} <span style='color:{style.active()['accent']}'>●</span>"
                           if modified else self.title)
        self.reset.setVisible(self.model.is_overridden(self.key))
        # what the choice means: on hover over the drop-down; only its status is shown
        desc = self.notes.get(str(value), "")
        self.combo.setToolTip((f"<p>{desc}</p>" if desc else "") + self._help_tip)
        note = ""
        if self.status is not None:
            st = self._mark(self.status(self.model), value)
            tip = setup_note(st)
            if tip:
                colour, what = _what(st)
                sub = PLACEMENT_SUBTAB.get(str(value)) if self.key == PLACEMENT_KEY else None
                if sub and st.state == "pending":     # details on the sub-tab: one sentence here
                    note = f"See the {sub} tab for necessary settings."
                else:
                    where = f" (see the {sub} tab)" if sub else f": {st.reason}"
                    note = f"<span style='color:{colour}'>{what}</span>{where}"
            self.note.setToolTip(tip)
        self.note.setText(note)
        self.note.setVisible(bool(note))

    def _mark(self, status: dict, current):
        """Choices marked by availability, as the form's choice fields do (fields.set_item_status):
        'unavailable' ones are disabled (the current value stays selectable), 'pending' ones
        noted; the reason is the tooltip. Returns the current choice's status."""
        model = self.combo.model()
        for i in range(self.combo.count()):
            v = self.combo.itemData(i)
            st = status.get(v) if isinstance(v, str) else None
            item = model.item(i) if hasattr(model, "item") else None
            if st is None or item is None:
                continue
            node = self.model.schema.get(self.key)
            label = labels.choice(node, v) + {"unavailable": "  (unavailable)", "pending": "  (needs setup)"}.get(st.state, "")
            item.setEnabled(st.state != "unavailable" or v == current)
            item.setText(label)
            item.setToolTip(choice_tip(node, v, st.reason or "All inputs found"))
        return status.get(current) if isinstance(current, str) else None

    def focus_key(self, path: str, ti: int | None = None) -> None:
        self.combo.setFocus()
