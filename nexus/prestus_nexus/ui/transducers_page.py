"""Transducers sub-tab of the Placement tab (shown with more than one transducer): one row
per transducer, the grid orientation under grid.mode with each beam's tilt, the async grid
mode, and what the single timing block means for several transducers.

The geometry is multitransducer.py (ports of mean_beam_axis.m / preproc_head.m). The table
is read-only: positions are edited in the placement viewer, where clicks, snapping and
undo already work per transducer ("Edit in placement viewer")."""

from __future__ import annotations

from typing import Callable

import numpy as np
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (QAbstractItemView, QFrame, QGridLayout, QHBoxLayout, QHeaderView,
                             QLabel, QPushButton, QScrollArea, QTableWidget, QTableWidgetItem,
                             QVBoxLayout, QWidget)

from .. import multitransducer as mt
from .. import placement as pl
from ..model import ConfigModel, is_empty
from . import style
from .empty_state import EmptyStack
from .readout import muted, status

COLUMNS = ("Transducer", "Device", "Transducer position", "Focus position", "Focal depth (EP)",
           "Target ISPPA", "Tilt from grid z")
UNIT = {"voxel": "vox", "ras": "mm RAS", "mni": "mm MNI"}

TIMING_NOTE = (
    "PRESTUS has one <b>timing</b> block for all transducers. Coherent coupling fires them "
    "together with that timing. Async coupling heats with the sum of the transducers' "
    "intensities under that one timing (combine_async_intensity.m, thermal_parameters.m), "
    "so the time-averaged heat source is Σ Iᵢ × pd/pri: right when each transducer fires with "
    "the timing block's duty cycle in its own, non-overlapping window; too high when the "
    "transducers share the timing block's on-time between them. The firing schedule itself "
    "is not modelled.")


def _fmt_vec(v, nd: int = 1) -> str:
    if isinstance(v, (list, tuple, np.ndarray)) and len(v) == 3:
        try:
            return "[" + ", ".join(f"{float(x):.{nd}f}".rstrip("0").rstrip(".") for x in v) + "]"
        except (TypeError, ValueError):
            pass
    return "—"


def _num(v) -> float | None:
    return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else None


class _Card(QFrame):
    def __init__(self, title: str, sub: str = "", parent: QWidget | None = None):
        super().__init__(parent)
        self.setObjectName("Card")
        self.body = QVBoxLayout(self)
        self.body.setContentsMargins(14, 12, 14, 12)
        self.body.setSpacing(8)
        t = QLabel(title)
        t.setObjectName("CardTitle")
        self.body.addWidget(t)
        if sub:
            s = QLabel(sub)
            s.setObjectName("CardSub")
            s.setWordWrap(True)
            self.body.addWidget(s)


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
    """Placement ▸ Transducers. ``viewer_getter`` returns the PlacementPage (or None); when it
    has the T1 loaded, positions come from it (also for Localite / heuristic placement) and
    are taken into the T1's RAS+ voxel grid as PRESTUS does."""

    changed = pyqtSignal()
    edit_in_viewer = pyqtSignal(int)      # transducer index: MainWindow shows the Viewer sub-tab
    navigate = pyqtSignal(str)            # config key to show (MainWindow.goto)

    def __init__(self, model: ConfigModel, viewer_getter: Callable[[], object] | None = None,
                 parent: QWidget | None = None):
        super().__init__(parent)
        self.model = model
        self.viewer_getter = viewer_getter
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
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
        table_card = _Card("Transducers", "One row per transducer. Positions are edited in the "
                                          "placement viewer; the device and ISPPA target on the "
                                          "Transducer tab.")
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
        table_card.body.addWidget(self.table)
        row = QHBoxLayout()
        self.edit_btn = QPushButton("Edit in placement viewer")
        self.edit_btn.setToolTip("Show the selected transducer in the Viewer, where clicks and the "
                                 "Place actions edit its positions")
        self.edit_btn.clicked.connect(lambda: self.edit_in_viewer.emit(max(0, self.table.currentRow())))
        row.addWidget(self.edit_btn)
        row.addStretch(1)
        table_card.body.addLayout(row)
        self.space_note = QLabel()
        self.space_note.setObjectName("CardSub")
        self.space_note.setWordWrap(True)
        table_card.body.addWidget(self.space_note)
        v.addWidget(table_card)

        # (b) grid orientation
        orient = _Card("Grid orientation", "PRESTUS rotates the head so one direction becomes the "
                                           "grid z axis (grid.mode). A beam along z is voxelised "
                                           "cleanly; an oblique beam crosses the grid at an angle.")
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
        v.addWidget(orient)

        # (c) coupling: the async grid mode (the coupling card itself lives on the Simulation tab)
        coupling = _Card("Coupling", "How the transducers are combined (simulation.transducer_coupling).")
        self.coupling_text = QLabel()
        self.coupling_text.setWordWrap(True)
        self.coupling_text.setTextFormat(Qt.TextFormat.RichText)
        coupling.body.addWidget(self.coupling_text)
        self.async_row = QWidget()
        ag = QGridLayout(self.async_row)
        ag.setContentsMargins(0, 0, 0, 0)
        ag.setColumnStretch(1, 1)
        ag.addWidget(QLabel("Async grid mode"), 0, 0, Qt.AlignmentFlag.AlignTop)
        self.async_editor = ChoiceEditor(model.schema.get("simulation.async_grid_mode"))
        self.async_editor.edited.connect(lambda val: self._set("simulation.async_grid_mode", val))
        ag.addWidget(self.async_editor, 0, 1)
        coupling.body.addWidget(self.async_row)
        self.coupling_box = coupling.body           # MainWindow may add the coupling card here
        v.addWidget(coupling)

        # (d) timing
        timing = _Card("Timing")
        t = QLabel(TIMING_NOTE)
        t.setWordWrap(True)
        t.setTextFormat(Qt.TextFormat.RichText)
        timing.body.addWidget(t)
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
    def refresh(self) -> None:
        m = self.model
        n = m.n_transducers
        if n <= 1:
            self.stack.show_empty("One transducer",
                                  "This tab compares several transducers: their beams, the grid "
                                  "orientation and the coupling. Add a transducer in the Viewer.",
                                  "Open the Viewer", lambda: self.edit_in_viewer.emit(0))
            return
        self.stack.show_content()
        txs, how = self.transducers()
        o = mt.orientation(m, txs)
        self._orientation = o
        self._fill_table(txs, o)
        self.space_note.setText(f"Tilts are computed from {how}.")
        self.mode_editor.set_value(m.get("grid.mode"), m.baseline("grid.mode"))
        self._fill_orientation(o)
        self._fill_coupling()

    def _fill_table(self, txs: list[dict], o: dict) -> None:
        m = self.model
        eff = m.effective_transducers()
        unit = UNIT.get(self._space(), "")
        self.table.setRowCount(len(eff))
        for ti, t in enumerate(eff):
            dev = " / ".join(str(x) for x in (t.get("serial"), (t.get("combo") or {}).get("ds_serial"))
                             if not is_empty(x)) or str(t.get("type") or "custom")
            bowl = pl.Bowl.from_transducer(t)
            ep = pl.resolve_focal_distance_ep(t, bowl)
            isppa = t.get("target_isppa_wcm2")
            if isinstance(isppa, list):
                isppa_s = ", ".join(f"{x:g}" for x in isppa if _num(x) is not None) + " W/cm²" if isppa else "drive level"
            else:
                isppa_s = f"{isppa:g} W/cm²" if _num(isppa) is not None else "drive level"
            tilt = o["tilts"][ti] if ti < len(o["tilts"]) else None
            cells = [m.transducer_title(ti), dev,
                     f"{_fmt_vec(t.get('trans_pos'))} {unit}" if _fmt_vec(t.get('trans_pos')) != "—" else "—",
                     f"{_fmt_vec(t.get('focus_pos'))} {unit}" if _fmt_vec(t.get('focus_pos')) != "—" else "—",
                     f"{ep:.1f} mm" if ep is not None else "from positions",
                     isppa_s,
                     f"{tilt:.1f}°" if tilt is not None else "—"]
            for col, text in enumerate(cells):
                item = QTableWidgetItem(text)
                if col == len(cells) - 1 and ti in o.get("coarse_tilted", []):
                    item.setToolTip("Tilted beam on a coarse grid; see Grid orientation")
                self.table.setItem(ti, col, item)
        h = self.table.horizontalHeader().height() + sum(self.table.rowHeight(r) for r in range(len(eff)))
        self.table.setFixedHeight(h + 2 * self.table.frameWidth())

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
            self.async_editor.set_value(m.get("simulation.async_grid_mode"),
                                        m.baseline("simulation.async_grid_mode"))
        else:
            self.coupling_text.setText(
                f"<b>{coupling}</b>: all transducers in one run on one grid, oriented by grid.mode "
                "above. " + muted("The async grid mode applies to async coupling only."))
            self.async_row.hide()
