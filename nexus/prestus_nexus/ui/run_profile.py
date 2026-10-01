"""Run profile picker: named hpc.* bundles (run_presets.yaml), shown above the HPC
section's fields. Applying a preset sets those fields on the open study -- it does not
submit anything; see doc_hpc.md for the actual srun/sbatch workflow."""

from __future__ import annotations

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import QComboBox, QHBoxLayout, QLabel, QPushButton, QVBoxLayout, QWidget

from .. import presets


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
        for p in self.presets:
            self.combo.addItem(p.name)
        self.combo.currentIndexChanged.connect(self._show_description)
        row.addWidget(self.combo, 1)
        apply_btn = QPushButton("Apply")
        apply_btn.setToolTip("Sets the hpc.* fields below from this preset. Never submits a job.")
        apply_btn.clicked.connect(self._apply)
        row.addWidget(apply_btn)
        v.addLayout(row)
        self.description = QLabel()
        self.description.setWordWrap(True)
        self.description.setObjectName("CardSub")
        v.addWidget(self.description)
        self._show_description(self.combo.currentIndex())

    def _show_description(self, index: int) -> None:
        if 0 <= index < len(self.presets):
            self.description.setText(self.presets[index].description)

    def _apply(self) -> None:
        index = self.combo.currentIndex()
        if 0 <= index < len(self.presets):
            presets.apply(self.window.model, self.presets[index])
            self.changed.emit()
