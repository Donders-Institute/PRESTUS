"""Time-reversal pose table viewer: the poses a pose search scored
(<prefix>_timereversal_poses_T<i>.csv, time_reversal_stage.m), best first, with a
"Use this pose…" button that starts a follow-up run at the selected pose."""

from __future__ import annotations

from pathlib import Path

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (QAbstractItemView, QComboBox, QHBoxLayout, QLabel, QPushButton,
                             QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget)

from .. import followup as fu

# (column, header) shown, in this order; reciprocity_score_poses.m names
COLUMNS = (
    ("pose", "Pose"),
    ("kind", "Kind"),
    ("p_corrected_rel", "Target pressure (rel.)"),
    ("p_geometric_rel", "Geometric aim (rel.)"),
    ("polar_deg", "Tilt from planned [°]"),
    ("azimuth_deg", "Azimuth [°]"),
    ("aim_offset_mm", "Aim offset [mm]"),
    ("residual_shift_mm", "Residual shift [mm]"),
    ("skin_mm", "Skin gap [mm]"),
    ("complete", "Complete"),
)


def _fmt(v) -> str:
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        return "" if v != v else f"{v:.3g}"
    return "" if v is None else str(v)


class PoseTable(QWidget):
    """Pose tables of one run. ``use_pose(dict)``: "Use this pose…" on the selected row."""

    use_pose = pyqtSignal(object)
    selection_changed = pyqtSignal(object)      # the selected row (dict) or None

    def __init__(self, show_use: bool = True, parent: QWidget | None = None):
        super().__init__(parent)
        self._tables: list[Path] = []
        self.rows: list[dict] = []
        v = QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 0)
        bar = QHBoxLayout()
        self.which = QComboBox()
        self.which.setToolTip("One pose table per transducer")
        self.which.currentIndexChanged.connect(self._load_table)
        bar.addWidget(self.which)
        self.note = QLabel()
        self.note.setObjectName("CardSub")
        self.note.setWordWrap(True)
        bar.addWidget(self.note, 1)
        self.use = QPushButton("Use this pose…")
        self.use.setToolTip("Create a follow-up config that simulates the selected pose "
                            "(manual placement in RAS mm)")
        self.use.setEnabled(False)
        self.use.clicked.connect(lambda: self.selected() is not None and self.use_pose.emit(self.selected()))
        self.use.setVisible(show_use)
        bar.addWidget(self.use)
        v.addLayout(bar)
        self.table = QTableWidget()
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.itemSelectionChanged.connect(self._selected)
        self.table.cellDoubleClicked.connect(lambda *_: self.use.isVisible() and self.use.click())
        v.addWidget(self.table, 1)

    def load(self, tables: list[Path]) -> None:
        self._tables = list(tables)
        self.which.blockSignals(True)
        self.which.clear()
        self.which.addItems([f"Transducer {p.stem.rsplit('_T', 1)[-1]}" for p in self._tables])
        self.which.blockSignals(False)
        self.which.setVisible(len(self._tables) > 1)
        self._load_table(0)

    def _load_table(self, index: int) -> None:
        path = self._tables[index] if 0 <= index < len(self._tables) else None
        try:
            self.rows = fu.read_pose_table(path) if path is not None else []
        except (OSError, ValueError) as exc:
            self.rows = []
            self.note.setText(f"Cannot read {path.name}: {exc}")
        else:
            ras = bool(self.rows) and all(k in self.rows[0] for k in fu.POSE_RAS_COLUMNS)
            self.note.setText(
                "" if path is None else
                f"{len(self.rows)} poses, best first by predicted target pressure (relative to the "
                "planned pose). Predictions rank poses; verify the chosen one with a simulation."
                + ("" if ras or not self.rows else
                   " This table has no RAS positions (search run without a T1 mapping), so its "
                   "poses cannot be used as a placement."))
        cols = [(k, h) for k, h in COLUMNS if not self.rows or k in self.rows[0]]
        self.table.clear()
        self.table.setColumnCount(len(cols))
        self.table.setHorizontalHeaderLabels([h for _, h in cols])
        self.table.setRowCount(len(self.rows))
        for r, row in enumerate(self.rows):
            for c, (k, _) in enumerate(cols):
                self.table.setItem(r, c, QTableWidgetItem(_fmt(row.get(k))))
        self.table.resizeColumnsToContents()
        if self.rows:
            self.table.selectRow(0)
        self._selected()

    def selected(self) -> dict | None:
        r = self.table.currentRow()
        return self.rows[r] if 0 <= r < len(self.rows) and self.table.selectedItems() else None

    def _selected(self) -> None:
        row = self.selected()
        ok = row is not None and all(k in row for k in fu.POSE_RAS_COLUMNS)
        self.use.setEnabled(ok)
        self.selection_changed.emit(row if ok else None)
