"""Tab 5 – Medium Properties (editable table)."""

from __future__ import annotations
from PyQt6.QtWidgets import QTableWidget, QTableWidgetItem, QLabel
from PyQt6.QtCore import Qt
from tabs.base import BaseTab
import param_io
import style


class TabMedium(BaseTab):
    def __init__(self, registry: dict, parent=None):
        super().__init__(registry, parent)

        self._add_note(
            "Acoustic and thermal tissue properties. Leave cells at defaults unless you "
            "need to override them."
        )

        col_names = [
            "Sound Speed\n[m/s]", "Density\n[kg/m³]", "α coeff\n[dB/cm/MHz]",
            "α power", "Therm. cond.", "Heat cap.", "Perfusion", "Absorb. frac.",
        ]

        self.table = QTableWidget(len(param_io.TISSUES), len(param_io.MEDIUM_PROPS))
        self.table.setHorizontalHeaderLabels(col_names)
        self.table.setVerticalHeaderLabels(param_io.TISSUES)
        self.table.horizontalHeader().setDefaultSectionSize(110)
        self.table.verticalHeader().setDefaultSectionSize(28)
        self.table.setStyleSheet(
            f"QTableWidget {{ background:{style.BG_EDIT}; gridline-color:{style.BORDER}; }}"
            f"QHeaderView::section {{ background:{style.BG_FIG}; border:1px solid {style.BORDER}; "
            f"padding:4px; font-size:{style.FS_SM}px; }}"
        )
        self.table.setMinimumHeight(220)

        for ti, tissue in enumerate(param_io.TISSUES):
            defaults = param_io.MEDIUM_DEFAULTS[tissue]
            for pi, val in enumerate(defaults):
                item = QTableWidgetItem(str(val))
                item.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                self.table.setItem(ti, pi, item)

        self._layout.addWidget(self.table)
        self._add_note(
            "Tip: tissue layers can be excluded in the Advanced tab by modifying the "
            "active layers list."
        )
        self._finish()
