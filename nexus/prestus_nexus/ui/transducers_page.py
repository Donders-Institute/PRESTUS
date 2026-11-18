"""Transducers overview (shown with more than one transducer), embedded in the Placement card on
Placement ▸ Settings (TransducersPage(embedded=True)), or as a page of its own: one row
per transducer, the grid orientation under grid.mode with each beam's tilt, the async grid
mode, and what the single timing block means for several transducers.

The geometry is multitransducer.py (ports of mean_beam_axis.m / preproc_head.m). The table
is read-only: positions are edited in the placement viewer, where clicks, snapping and
undo already work per transducer ("Edit in placement viewer")."""

from __future__ import annotations

import html

from typing import Callable

import numpy as np
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import (QAbstractItemView, QFrame, QGridLayout, QHBoxLayout, QHeaderView,
                             QLabel, QPushButton, QScrollArea, QTableWidget, QTableWidgetItem,
                             QVBoxLayout, QWidget)

from .. import multitransducer as mt
from .. import placement as pl
from ..model import ConfigModel, is_empty
from . import style
from .empty_state import EmptyStack
from .readout import muted, status
from .head_sketch import HeadSketch

# Compact columns: the number (in the transducer's colour, as on the head sketch), the device,
# and values without units; units and full names are in the headers and their tooltips.
COLUMNS = ("", "Transducer", "Position", "Focus", "Depth EP (mm)", "ISPPA (W/cm²)", "Tilt (°)")
COLUMN_TIPS = ("Transducer number and colour (as on the head sketch)",
               "Device: equipment-library serial (driving system in the tooltip)",
               "Transducer position: bowl apex (transducer.trans_pos)",
               "Focus position: the target (transducer.focus_pos)",
               "Focal depth setting from the exit plane (transducer.focal_distance_ep); "
               "'auto' = from the positions",
               "Target free-water ISPPA (transducer.target_isppa_wcm2); 'drive' = the drive level",
               "Tilt of the beam from the grid z axis (grid.mode)")
UNIT = {"voxel": "vox", "ras": "mm RAS", "mni": "mm MNI"}

TIMING_NOTE = (
    "Coherent coupling fires all transducers together with the <b>timing</b> block. Async "
    "coupling switches each transducer's heat source on its own schedule "
    "(<b>timing.async_schedule</b>, thermal_schedule.m): by default each uses the timing "
    "block and they take turns within each PRI ('pulses'); 'trains' alternates pulse trains; "
    "'custom' takes a start delay and timing overrides per transducer. The diagram under "
    "Sonication › Timing shows who fires when; its Sonication sequence card sets up stances "
    "fired one after another.")


def _fmt_vec(v, nd: int = 1) -> str:
    if isinstance(v, (list, tuple, np.ndarray)) and len(v) == 3:
        try:
            return ", ".join(f"{float(x):.{nd}f}".rstrip("0").rstrip(".") for x in v)
        except (TypeError, ValueError):
            pass
    return "—"


def _num(v) -> float | None:
    return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else None


class _Card(QFrame):
    """A card; ``plain``: a titled block without its own frame (inside another card)."""

    def __init__(self, title: str, sub: str = "", parent: QWidget | None = None, plain: bool = False,
                 first: bool = False):
        super().__init__(parent)
        if not plain:
            self.setObjectName("Card")
        self.body = QVBoxLayout(self)
        self.body.setContentsMargins(*(((0, 0 if first else 10, 0, 4)) if plain else (14, 12, 14, 12)))
        self.body.setSpacing(8)
        t = QLabel(title)
        t.setObjectName("CardTitle")
        if sub:                            # what the block is for: on hover
            t.setToolTip(sub)
        if plain:
            f = t.font()
            f.setBold(True)
            t.setFont(f)
        self.body.addWidget(t)


class _Note(QLabel):
    """A plain tinted box (no accent bar) for a warning or a remark."""

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.setWordWrap(True)
        self.setTextFormat(Qt.TextFormat.RichText)
        self.hide()

    def show_lines(self, lines: list[str], level: str = "warn") -> None:
        if not lines:
            self.hide()
            return
        c = style.active()
        bg = c["thermal_soft"] if level == "warn" else c["accent_soft"]
        self.setStyleSheet(f"QLabel {{ background: {bg}; border: 1px solid {c['line']}; "
                           f"border-radius: 6px; padding: 6px 8px; color: {c['ink']}; }}")
        self.setText("<br>".join(lines))
        self.show()


class TransducersPage(QWidget):
    """The transducers overview. ``viewer_getter`` returns the PlacementPage (or None); when it
    has the T1 loaded, positions come from it (also for Localite / heuristic placement) and
    are taken into the T1's RAS+ voxel grid as PRESTUS does."""

    changed = pyqtSignal()
    edit_in_viewer = pyqtSignal(int)      # transducer index: MainWindow shows the Viewer sub-tab
    navigate = pyqtSignal(str)            # config key to show (MainWindow.goto)
    selected = pyqtSignal(int)            # row picked: the transducer the fields below edit
    add_requested = pyqtSignal()
    remove_requested = pyqtSignal(int)

    def __init__(self, model: ConfigModel, viewer_getter: Callable[[], object] | None = None,
                 parent: QWidget | None = None, embedded: bool = False):
        """``embedded``: inside the Placement card on Settings (no scroll area, blocks without
        frames, hidden with one transducer); the table then selects the transducer to edit."""
        super().__init__(parent)
        self.model = model
        self.viewer_getter = viewer_getter
        self.embedded = embedded
        self._filling = False
        self.current = 0
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        if embedded:
            self.stack = None
            v = outer
            v.setSpacing(4)
        else:
            content = QWidget()
            self.stack = EmptyStack(content)
            outer.addWidget(self.stack)
            scroll = QScrollArea()
            scroll.setWidgetResizable(True)
            scroll.setFrameShape(QFrame.Shape.NoFrame)
            inner = QWidget()
            scroll.setWidget(inner)
            cv = QVBoxLayout(content)
            cv.setContentsMargins(0, 0, 0, 0)
            cv.addWidget(scroll)
            v = QVBoxLayout(inner)
            v.setSpacing(12)

        # (a) one row per transducer
        table_card = _Card("Transducers", "One row per transducer: select one to edit its values "
                                          "below. Positions are also edited in the placement "
                                          "viewer; the device and ISPPA target on the Transducer "
                                          "tab." if embedded else
                                          "One row per transducer. Positions are edited in the "
                                          "placement viewer; the device and ISPPA target on the "
                                          "Transducer tab.", plain=embedded, first=True)
        self.table = QTableWidget(0, len(COLUMNS))
        self.table.setHorizontalHeaderLabels(COLUMNS)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.table.cellDoubleClicked.connect(lambda r, _c: self.edit_in_viewer.emit(r))
        self.table.currentCellChanged.connect(self._row_changed)
        trow = QHBoxLayout()
        trow.addWidget(self.table, 1)
        self.sketch = HeadSketch(150)            # where the transducers sit, seen from above
        trow.addWidget(self.sketch, 0, Qt.AlignmentFlag.AlignTop)
        table_card.body.addLayout(trow)
        row = QHBoxLayout()
        self.edit_btn = QPushButton("Edit in placement viewer")
        self.edit_btn.setToolTip("Show the selected transducer in the Viewer, where clicks and the "
                                 "Place actions edit its positions")
        self.edit_btn.clicked.connect(lambda: self.edit_in_viewer.emit(max(0, self.table.currentRow())))
        row.addWidget(self.edit_btn)
        self.add_btn = QPushButton("Add transducer")
        self.add_btn.setToolTip("Add a transducer; it starts from transducer 1's values")
        self.add_btn.clicked.connect(self.add_requested.emit)
        self.remove_btn = QPushButton("Remove transducer")
        self.remove_btn.setToolTip("Remove the selected transducer")
        self.remove_btn.clicked.connect(lambda: self.remove_requested.emit(max(0, self.table.currentRow())))
        for b in (self.add_btn, self.remove_btn):
            b.setVisible(embedded)
            row.addWidget(b)
        row.addStretch(1)
        table_card.body.addLayout(row)
        self.space_note = QLabel()
        self.space_note.setObjectName("CardSub")
        self.space_note.setWordWrap(True)
        table_card.body.addWidget(self.space_note)
        # embedded: which transducer the fields right below edit
        self.editing = QLabel()
        self.editing.setTextFormat(Qt.TextFormat.RichText)
        self.editing.setVisible(embedded)
        table_card.body.addWidget(self.editing)
        v.addWidget(table_card)

        # (b) grid orientation
        orient = _Card("Grid orientation", "PRESTUS rotates the head so one direction becomes the "
                                           "grid z axis (grid.mode). A beam along z is voxelised "
                                           "cleanly; an oblique beam crosses the grid at an angle.",
                       plain=embedded)
        from .fields import ChoiceEditor
        node = model.schema.get("grid.mode")
        self.mode_editor = ChoiceEditor(node)
        self.mode_editor.edited.connect(lambda val: self._set("grid.mode", val))
        grid = QGridLayout()
        grid.setColumnStretch(1, 1)
        grid.addWidget(QLabel("Grid orientation"), 0, 0, Qt.AlignmentFlag.AlignTop)
        grid.addWidget(self.mode_editor, 0, 1)
        orient.body.addLayout(grid)
        self.orient_text = QLabel()
        self.orient_text.setWordWrap(True)
        self.orient_text.setTextFormat(Qt.TextFormat.RichText)
        orient.body.addWidget(self.orient_text)
        self.orient_note = _Note()
        orient.body.addWidget(self.orient_note)
        # embedded: grid orientation and coupling go below the Placement card's fields
        # (self.details), so the table sits right above the positions it selects for
        self.details = QWidget()
        dv = QVBoxLayout(self.details)
        dv.setContentsMargins(0, 0, 0, 0)
        dv.setSpacing(4)
        (dv if embedded else v).addWidget(orient)

        # (c) coupling: the async grid mode (the coupling card itself lives on the Simulation tab)
        coupling = _Card("Coupling", "How the transducers are combined (simulation.transducer_coupling).",
                         plain=embedded)
        self.coupling_text = QLabel()
        self.coupling_text.setWordWrap(True)
        self.coupling_text.setTextFormat(Qt.TextFormat.RichText)
        coupling.body.addWidget(self.coupling_text)
        self.async_row = QWidget()
        ag = QGridLayout(self.async_row)
        ag.setContentsMargins(0, 0, 0, 0)
        ag.setColumnStretch(1, 1)
        ag.addWidget(QLabel("Async grid mode"), 0, 0, Qt.AlignmentFlag.AlignTop)
        self.async_editor = ChoiceEditor(model.schema.get("grid.async_layout"))
        self.async_editor.edited.connect(lambda val: self._set("grid.async_layout", val))
        ag.addWidget(self.async_editor, 0, 1)
        coupling.body.addWidget(self.async_row)
        self.coupling_box = coupling.body           # MainWindow may add the coupling card here
        (dv if embedded else v).addWidget(coupling)

        # (d) timing
        timing = _Card("Timing")
        t = QLabel(TIMING_NOTE)
        t.setWordWrap(True)
        t.setTextFormat(Qt.TextFormat.RichText)
        timing.body.addWidget(t)
        if not embedded:                 # on Settings, the Sonication tab says it
            v.addWidget(timing)
            v.addStretch(1)
        self._orientation: dict = {}

    # ── model ───────────────────────────────────────────────────────────
    def _set(self, key: str, value) -> None:
        self.model.set(key, value)
        self.changed.emit()
        self.refresh()

    def _space(self) -> str:
        return str(self.model.get("placement.manual.coordinate_space") or "voxel")

    def _viewer(self):
        v = self.viewer_getter() if self.viewer_getter is not None else None
        return v if v is not None and getattr(v, "t1", None) is not None else None

    def transducers(self) -> tuple[list[dict], str]:
        """Positions for the orientation and the space they are in (for the note)."""
        v = self._viewer()
        if v is not None:
            t1 = v.t1
            return (mt.model_transducers(self.model, v._positions_for,
                                         lambda r: t1.ras_to_canon(np.atleast_2d(r))[0]),
                    "the T1's RAS+ voxel grid, as PRESTUS computes it")
        space = self._space()
        how = {"voxel": "configured voxel positions (assumes a RAS+-stored T1; open the Viewer to use the image)",
               "ras": "configured RAS mm positions (exact for isotropic T1 voxels)",
               "mni": "configured MNI mm positions (approximate; the subject's head is not aligned to MNI)"}
        return mt.model_transducers(self.model), how.get(space, space)

    # ── refresh ─────────────────────────────────────────────────────────
    def _row_changed(self, row: int, *_args) -> None:
        if self._filling or row < 0:
            return
        self.current = row
        self.selected.emit(row)

    def select(self, ti: int) -> None:
        """Show transducer ti as the selected row (without emitting ``selected``)."""
        self.current = ti
        if 0 <= ti < self.table.rowCount() and self.table.currentRow() != ti:
            self._filling = True
            self.table.selectRow(ti)
            self._filling = False
        self._editing_text()

    def _editing_text(self) -> None:
        if self.embedded and 0 <= self.current < self.model.n_transducers:
            self.editing.setText(f"<b>Editing {html.escape(self.model.transducer_title(self.current))}</b>")

    def refresh(self) -> None:
        m = self.model
        n = m.n_transducers
        if self.embedded:
            self.setVisible(n > 1)
            self.details.setVisible(n > 1)
            if n <= 1:
                return
        elif n <= 1:
            self.stack.show_empty("One transducer",
                                  "This tab compares several transducers: their beams, the grid "
                                  "orientation, and the coupling. Add a transducer in the Viewer.",
                                  "Open the Viewer", lambda: self.edit_in_viewer.emit(0))
            return
        if self.stack is not None:
            self.stack.show_content()
        txs, how = self.transducers()
        o = mt.orientation(m, txs)
        self._orientation = o
        self._fill_table(txs, o)
        self.sketch.set_heads(mt.schematic_top_view(txs), how)
        self.space_note.setText(f"Tilts are computed from {how}.")
        self.space_note.hide()                       # on hover over the Tilt header instead
        tilt = self.table.horizontalHeaderItem(len(COLUMNS) - 1)
        if tilt is not None:
            tilt.setToolTip(f"{COLUMN_TIPS[-1]}. Computed from {how}.")
        self.mode_editor.set_value(m.get("grid.mode"), m.baseline("grid.mode"))
        self._fill_orientation(o)
        self._fill_coupling()

    def _fill_table(self, txs: list[dict], o: dict) -> None:
        m = self.model
        eff = m.effective_transducers()
        unit = UNIT.get(self._space(), "")
        self._filling = True
        self.table.setRowCount(len(eff))
        for col, (name, tip) in enumerate(zip(COLUMNS, COLUMN_TIPS)):
            head = QTableWidgetItem(f"{name} ({unit})" if col in (2, 3) and unit else name)
            head.setToolTip(tip)
            self.table.setHorizontalHeaderItem(col, head)
        for ti, t in enumerate(eff):
            serial = t.get("serial")
            lib = m.library
            dev = (lib.short_name(serial) if lib is not None else str(serial)) if not is_empty(serial) \
                else str(t.get("type") or "custom")
            ds = (t.get("combo") or {}).get("ds_serial")
            if not is_empty(ds) and lib is not None:
                full = (lib.driving_systems.get(str(ds)) or {}).get("name") or ""
                ds = f"{ds} ({full})" if full and full != ds else str(ds)
            bowl = pl.Bowl.from_transducer(t)
            ep = pl.resolve_focal_distance_ep(t, bowl)
            isppa = t.get("target_isppa_wcm2")
            if isinstance(isppa, list):
                isppa_s = ", ".join(f"{x:g}" for x in isppa if _num(x) is not None) if isppa else "drive"
            else:
                isppa_s = f"{isppa:g}" if _num(isppa) is not None else "drive"
            tilt = o["tilts"][ti] if ti < len(o["tilts"]) else None
            cells = [str(ti + 1), dev, _fmt_vec(t.get("trans_pos")), _fmt_vec(t.get("focus_pos")),
                     f"{ep:.3g}" if ep is not None else "auto", isppa_s,
                     f"{tilt:.1f}" if tilt is not None else "—"]
            for col, text in enumerate(cells):
                item = QTableWidgetItem(text)
                if col == 0:                       # the transducer's colour, as on the sketch
                    item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    item.setBackground(QColor(style.series(ti)))
                    item.setForeground(QColor("#ffffff"))
                    f = item.font()
                    f.setBold(True)
                    item.setFont(f)
                    item.setToolTip(m.transducer_title(ti))
                elif col == 1:
                    item.setToolTip(m.transducer_tooltip(ti) + (f"\nDriving system: {ds}" if not is_empty(ds) else ""))
                elif col == len(cells) - 1 and ti in o.get("coarse_tilted", []):
                    item.setToolTip("Tilted beam on a coarse grid; see Grid orientation")
                self.table.setItem(ti, col, item)
        self.current = min(self.current, len(eff) - 1)
        self.table.selectRow(self.current)
        self._filling = False
        self.remove_btn.setEnabled(len(eff) > 1)
        self._editing_text()
        self._fit_table()

    def _fit_table(self) -> None:
        """Table as tall as its rows, plus the horizontal scroll bar when the columns overflow
        (the head sketch beside it takes width), so no row hides behind the bar."""
        t = self.table
        h = t.horizontalHeader().height() + sum(t.rowHeight(r) for r in range(t.rowCount()))
        cols = sum(t.columnWidth(c) for c in range(t.columnCount()))
        if cols > t.viewport().width():
            h += t.horizontalScrollBar().sizeHint().height()
        t.setFixedHeight(h + 2 * t.frameWidth())

    def resizeEvent(self, event) -> None:  # noqa: N802 (Qt API)
        super().resizeEvent(event)
        self._fit_table()

    def _fill_orientation(self, o: dict) -> None:
        lines = []
        warn = []
        if o["error"]:
            self.orient_text.setText(status(o["error"], "bad" if "removed" in o["error"] else "warn"))
            self.orient_note.show_lines([])
            return
        z = o["z"]
        lines.append(f"Grid z: {_fmt_vec(z, 3)}"
                     + (muted("  (T1 superior axis: no rotation)") if o["mode"] == "ras_plus" else ""))
        tilts = ", ".join(f"{i + 1}: {a:.1f}°" if a is not None else f"{i + 1}: —"
                          for i, a in enumerate(o["tilts"]))
        if self.embedded:                   # the table's Tilt column shows the tilts
            self.orient_text.setText("")
            self.orient_text.hide()
            self.mode_editor.setToolTip("<br>".join(lines))
        else:
            lines.append(f"Beam tilt from grid z: {tilts}")
            self.orient_text.setText("<br>".join(lines))
        if o["fell_back"]:
            warn.append(status(f"The beams are too spread to average (mean-axis norm below "
                               f"grid.focal_axis_min_norm = {o['min_norm']:g}); PRESTUS aligns z with "
                               "transducer 1 instead and warns. Consider async coupling with "
                               "per_transducer grids.", "warn"))
        if o["coarse_tilted"]:
            ppw = o["ppw"]
            which = ", ".join(str(i + 1) for i in o["coarse_tilted"])
            warn.append(status(f"Transducer {which}: tilted ≥ {mt.TILT_WARN_DEG:g}° on a grid with "
                               f"{ppw:.1f} points per wavelength in water.", "warn")
                        + "<br>" + muted(
                            f"Guidance, not a rule: one test lost 16% of free-water ISPPA for a beam "
                            f"41° off axis at 2.5 points per wavelength (2 mm, 300 kHz; "
                            f"doc_multitransducer_async.md). The Nexus flags tilts from "
                            f"{mt.TILT_WARN_DEG:g}° below {mt.PPW_WARN:g} points per wavelength. A finer "
                            "grid, or async coupling with per_transducer grids, keeps each beam "
                            "on its own axis."))
        self.orient_note.show_lines(warn)

    def _fill_coupling(self) -> None:
        m = self.model
        coupling = m.get("simulation.transducer_coupling") or "coherent"
        if coupling == "async":
            self.coupling_text.setText(
                "<b>async</b>: one acoustic run per transducer, combined afterwards (no interference).")
            self.async_row.show()
            self.async_editor.set_value(m.get("grid.async_layout"),
                                        m.baseline("grid.async_layout"))
        else:
            self.coupling_text.setText(
                f"<b>{coupling}</b>: all transducers in one run on one grid, oriented by grid.mode "
                "above. " + muted("The async grid mode applies to async coupling only."))
            self.async_row.hide()
