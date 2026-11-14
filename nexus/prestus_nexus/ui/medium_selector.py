"""A tab's header row for the choice that decides what the tab shows: simulation.medium on
the Medium tab (water: tissue properties of water only; phantom: tissue layers and all
properties; layered: everything, including the CT-based skull mapping), placement.mode on the Placement
tab (PLACEMENT_NOTE)."""

from __future__ import annotations

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import QComboBox, QFrame, QHBoxLayout, QLabel, QToolButton

from . import icons, style

KEY = "simulation.medium"
NOTE = {
    "water": "Free field in water: only the water properties apply.",
    "layered": "Head model from the SimNIBS segmentation: head processing, tissue properties and "
               "optional skull mapping from a measured CT or a UTE pseudo-CT.",
    "phantom": "A designed test object: the segmentation is mapped onto tissues without "
               "preprocessing; only the tissues its layers assign are modelled.",
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
                self.combo.addItem(ic, str(choice), choice)
            else:
                self.combo.addItem(str(choice), choice)
        if node is not None:
            from .fields import display_text
            from .tooltip import html_for
            self.combo.setToolTip(html_for(node, title, display_text(node.default)))
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
        note = self.notes.get(str(value), "")
        if self.status is not None:
            st = self._mark(self.status(self.model), value)
            if st is not None and st.state != "ready" and st.reason:
                colour = style.active()["red" if st.state == "unavailable" else "amber"]
                what = "Cannot run yet" if st.state == "unavailable" else "Still to set up"
                note = f"{note}<br><span style='color:{colour}'>{what}:</span> {st.reason}"
        self.note.setText(note)

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
            label = str(v) + {"unavailable": "  (unavailable)", "pending": "  (needs setup)"}.get(st.state, "")
            item.setEnabled(st.state != "unavailable" or v == current)
            item.setText(label)
            item.setToolTip(st.reason or "All inputs found")
        return status.get(current) if isinstance(current, str) else None

    def focus_key(self, path: str, ti: int | None = None) -> None:
        self.combo.setFocus()
