"""Run profile picker: named hpc.* bundles (run_presets.yaml), shown above the HPC
section's fields. It shows the profile the study's settings match (or that none does);
applying a preset sets those fields on the open study -- it does not submit anything;
see doc_hpc.md for the actual srun/sbatch workflow."""

from __future__ import annotations

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import QComboBox, QHBoxLayout, QLabel, QPushButton, QVBoxLayout, QWidget

from .. import presets

CUSTOM = "Current settings (no profile)"


class RunProfilePanel(QWidget):
    changed = pyqtSignal()

    def __init__(self, window):
        super().__init__()
        self.window = window
        self.presets = presets.load_presets()

        v = QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 6)
        row = QHBoxLayout()
        row.addWidget(QLabel("Run profile"))
        self.combo = QComboBox()
        self.combo.addItem(CUSTOM)
        for p in self.presets:
            self.combo.addItem(p.name)
        self.combo.setToolTip("The profile the settings below match. Choose another and Apply "
                              "to set its values.")
        self.combo.currentIndexChanged.connect(self._show_description)
        row.addWidget(self.combo, 1)
        self.apply_btn = QPushButton("Apply")
        self.apply_btn.setToolTip("Sets the fields below from this profile. Never submits a job.")
        self.apply_btn.clicked.connect(self._apply)
        row.addWidget(self.apply_btn)
        v.addLayout(row)
        self.description = QLabel()
        self.description.setWordWrap(True)
        self.description.setObjectName("CardSub")
        v.addWidget(self.description)
        self.refresh()

    def _preset(self, index: int) -> presets.Preset | None:
        return self.presets[index - 1] if 1 <= index <= len(self.presets) else None

    def refresh(self) -> None:
        """Select the profile the study's settings match, or 'no profile'."""
        match = presets.current(self.window.model, self.presets)
        self.combo.blockSignals(True)
        self.combo.setCurrentIndex(self.presets.index(match) + 1 if match else 0)
        self.combo.blockSignals(False)
        self._show_description(self.combo.currentIndex())

    def showEvent(self, event) -> None:        # noqa: N802
        super().showEvent(event)
        self.refresh()

    def _show_description(self, index: int) -> None:
        p = self._preset(index)
        if p is None:
            self.description.setText("The settings below match no profile. Choose one and Apply "
                                     "to use its values.")
            self.apply_btn.setEnabled(False)
            return
        applied = presets.matches(self.window.model, p)
        self.description.setText(p.description + ("" if applied else "  (Not applied: press Apply.)"))
        self.apply_btn.setEnabled(not applied)

    def _apply(self) -> None:
        p = self._preset(self.combo.currentIndex())
        if p is not None:
            presets.apply(self.window.model, p)
            self.changed.emit()
            self.refresh()
