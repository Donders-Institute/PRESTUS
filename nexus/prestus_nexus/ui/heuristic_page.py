"""Placement ▸ Heuristic: choose among the heuristic's candidate positions.

transducer_positioning.m writes every candidate it evaluated to the placement run folder
(tpars_sub-XXX_<target>.csv: position, distance to the target and the beam criteria) and the
one it chose (sub-XXX_<target>.txt, after the criteria of tp_select_heuristic_position.m).
This page lists them, plots distance against the share of the beam path through excluded
tissue, and places a chosen candidate: placement.heuristic.candidate, used as it is instead
of the criteria, in a placement-only run that reuses the candidate list.
"""

from __future__ import annotations

import numpy as np
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (QAbstractItemView, QHBoxLayout, QHeaderView, QLabel, QPushButton,
                             QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget)

from .. import launcher
from ..heuristic_import import Candidate, candidates, selected
from ..model import is_empty
from ..outputs import dir_placement_run
from . import style
from .empty_state import EmptyStack
from .partial_run import RunFollower, run_partial

TASK = "Heuristic placement"
COLUMNS = (("idx", "Candidate", "{:.0f}"), ("dist_to_target", "Distance to target [vox]", "{:.1f}"),
           ("prop_intersect", "Beam through excluded tissue", "{:.3f}"),
           ("mean_dist_skin", "Skin distance, mean [vox]", "{:.2f}"),
           ("var_dist_skin", "Skin distance, var", "{:.2f}"),
           ("mean_dist_skull", "Skull distance, mean [vox]", "{:.2f}"),
           ("var_dist_skull", "Skull distance, var", "{:.2f}"))

try:
    from matplotlib.figure import Figure

    from .mpl_theme import ThemedCanvas
    HAVE_MPL = True
except ImportError:          # pragma: no cover
    HAVE_MPL = False


def run_info(model) -> tuple | None:
    """(run folder, subject id, target name) of the heuristic placement, or None."""
    run = dir_placement_run(model, "heuristic")
    sid = model.get("subject_id")
    target = model.get("placement.heuristic.target_name")
    if run is None or is_empty(target) or not isinstance(sid, (int, float)) or isinstance(sid, bool):
        return None
    return run, int(sid), str(target)


def heuristic_choice(cands: list[Candidate], sel) -> int | None:
    """idx of the candidate the heuristic placed (sub-XXX_<target>.txt; its transducer voxel is
    the candidate's, rounded), or None."""
    if sel is None or sel.trans_vox is None or not cands:
        return None
    d = [float(np.linalg.norm(np.round(c.trans_vox) - np.asarray(sel.trans_vox))) for c in cands]
    k = int(np.argmin(d))
    return cands[k].idx if d[k] <= 1.0 else None


def passes_intersection(c: Candidate, criterion) -> bool:
    """The heuristic's first filter (criterion_intersection, before it widens it)."""
    return not isinstance(criterion, (int, float)) or c.prop_intersect < criterion


class HeuristicPage(QWidget):
    def __init__(self, window, parent: QWidget | None = None):
        super().__init__(parent)
        self.window = window
        self.cands: list[Candidate] = []
        self.choice: int | None = None          # the heuristic's own pick
        self.current: int | None = None         # the candidate shown in the readout
        self._loaded = None
        self._stale = False
        self.follower = RunFollower(window, TASK, self._run_done)

        v = QVBoxLayout(self)
        intro = QLabel("The heuristic evaluates candidate positions on the scalp around the MNI target and "
                       "picks the closest one whose beam crosses little excluded tissue. Click a candidate in "
                       "the plot or the table to see its criteria, then <b>Place here</b> to use it instead.")
        intro.setObjectName("CardSub")
        intro.setWordWrap(True)
        v.addWidget(intro)
        bar = QHBoxLayout()
        self.compute_btn = QPushButton("Compute candidates")
        self.compute_btn.setToolTip("Run the heuristic placement again: evaluate the candidates and place its "
                                    "choice (clears placement.heuristic.candidate)")
        self.compute_btn.clicked.connect(self.compute)
        bar.addWidget(self.compute_btn)
        viewer = QPushButton("Open in Viewer")
        viewer.clicked.connect(self._viewer)
        bar.addWidget(viewer)
        bar.addStretch(1)
        v.addLayout(bar)

        self._content = QWidget()
        body = QVBoxLayout(self._content)
        body.setContentsMargins(0, 0, 0, 0)
        top = QHBoxLayout()
        self.canvas = None
        if HAVE_MPL:
            self.fig = Figure(figsize=(6, 3.6), layout="constrained")
            self.canvas = ThemedCanvas(self.fig, export_name="heuristic_candidates")
            self.canvas.setMinimumHeight(280)
            self.canvas.mpl_connect("button_press_event", self._clicked)
            top.addWidget(self.canvas, 3)
        side = QVBoxLayout()
        self.readout = QLabel()
        self.readout.setTextFormat(Qt.TextFormat.RichText)
        self.readout.setWordWrap(True)
        self.readout.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.readout.setMinimumWidth(260)
        side.addWidget(self.readout)
        self.place_btn = QPushButton("Place here")
        self.place_btn.setToolTip("Set placement.heuristic.candidate to this candidate and place it (Run tab; "
                                  "the candidate list is reused, not recomputed)")
        self.place_btn.clicked.connect(self.place)
        side.addWidget(self.place_btn)
        self.auto_btn = QPushButton("Use the heuristic's choice")
        self.auto_btn.setToolTip("Clear placement.heuristic.candidate and place the heuristic's choice")
        self.auto_btn.clicked.connect(self.use_heuristic)
        side.addWidget(self.auto_btn)
        side.addStretch(1)
        self.status = QLabel()
        self.status.setObjectName("CardSub")
        self.status.setWordWrap(True)
        side.addWidget(self.status)
        top.addLayout(side, 2)
        body.addLayout(top, 2)
        self.table = QTableWidget(0, len(COLUMNS))
        self.table.setHorizontalHeaderLabels([label for _k, label, _f in COLUMNS])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSortingEnabled(True)
        self.table.itemSelectionChanged.connect(self._table_selected)
        self.table.verticalHeader().setVisible(False)
        body.addWidget(self.table, 3)
        self.stack = EmptyStack(self._content)
        v.addWidget(self.stack, 1)

    @property
    def model(self):
        return self.window.model

    # ── data ────────────────────────────────────────────────────────────
    def chosen(self) -> int | None:
        v = self.model.get("placement.heuristic.candidate")
        return int(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else None

    def refresh(self) -> None:
        if self.window.isVisible() and not self.isVisible():
            self._stale = True
            return
        is_mode = self.model.get("placement.mode") == "heuristic"
        self.compute_btn.setEnabled(is_mode and not self.follower.waiting)
        info = run_info(self.model)
        if info is None:
            self.stack.show_empty("No heuristic run folder yet",
                                  "Set the subject, path.sim (or path.placement) and "
                                  "placement.heuristic.target_name (Settings).")
            return
        run, sid, target = info
        csv = run / f"tpars_sub-{sid:03d}_{target}.csv"
        try:
            key = (csv, csv.stat().st_mtime)
        except OSError:
            key = None
        if key is None:
            body = ("Compute candidates runs the heuristic for this subject and target (as a cluster job "
                    "where the Run tab submits there); its candidates then appear here. It needs the "
                    "SimNIBS segmentation." if is_mode else
                    "Choose 'heuristic' as the placement method (top of this tab) to run it.")
            self.stack.show_empty("No candidates yet", f"{body}<br>Looked for {csv}",
                                  "Compute candidates" if is_mode else None, self.compute if is_mode else None)
            return
        if key != self._loaded:
            self.cands = candidates(run, sid, target)
            self.choice = heuristic_choice(self.cands, selected(run, sid, target))
            self._loaded = key
            ids = {c.idx for c in self.cands}
            self.current = self.chosen() if self.chosen() in ids else self.choice
            self._fill_table()
        self.stack.show_content()
        self._plot()
        self._update_side()

    def showEvent(self, event) -> None:        # noqa: N802
        super().showEvent(event)
        if self._stale:
            self._stale = False
            self.refresh()

    def _fill_table(self) -> None:
        crit = self.model.get("placement.heuristic.criterion_intersection")
        self.table.blockSignals(True)
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(self.cands))
        for r, cand in enumerate(self.cands):
            for col, (key, _label, fmt) in enumerate(COLUMNS):
                value = getattr(cand, key)
                item = QTableWidgetItem()
                # a number (sorts numerically), rounded as the column shows it
                digits = int(fmt[3]) if fmt[3].isdigit() else 0
                item.setData(Qt.ItemDataRole.DisplayRole, int(value) if key == "idx" else round(float(value), digits))
                item.setData(Qt.ItemDataRole.UserRole, cand.idx)
                if not passes_intersection(cand, crit):
                    item.setForeground(Qt.GlobalColor.gray)
                if cand.idx == self.choice:
                    f = item.font()
                    f.setBold(True)
                    item.setFont(f)
                self.table.setItem(r, col, item)
        self.table.setSortingEnabled(True)
        self.table.sortItems(1)
        self.table.blockSignals(False)

    # ── plot ────────────────────────────────────────────────────────────
    def _plot(self) -> None:
        if self.canvas is None:
            return
        self.fig.clear()
        ax = self.fig.add_subplot()
        c = style.active()
        if not self.cands:
            self.canvas.draw_idle()
            return
        crit = self.model.get("placement.heuristic.criterion_intersection")
        x = np.array([k.dist_to_target for k in self.cands])
        y = np.array([k.prop_intersect for k in self.cands])
        ok = np.array([passes_intersection(k, crit) for k in self.cands])
        ax.scatter(x[~ok], y[~ok], s=14, color=c["ink4"], label="beam crosses too much excluded tissue")
        ax.scatter(x[ok], y[ok], s=14, color=c["accent"], label="passes the intersection criterion")
        if isinstance(crit, (int, float)):
            ax.axhline(crit, color=c["ink3"], linestyle="--", linewidth=0.8)
        marks = [(self.choice, "*", 160, c["green"], "heuristic's choice"),
                 (self.chosen(), "D", 60, c["ink"], "placement.heuristic.candidate"),
                 (self.current, "o", 140, "none", "selected here")]
        for idx, marker, size, colour, label in marks:
            k = next((q for q in self.cands if q.idx == idx), None)
            if k is None:
                continue
            ax.scatter([k.dist_to_target], [k.prop_intersect], s=size, marker=marker,
                       facecolors=colour if colour != "none" else "none",
                       edgecolors=c["red"] if colour == "none" else colour, linewidths=1.6, label=label,
                       zorder=5)
        ax.set_xlabel("Distance to target [voxels]")
        ax.set_ylabel("Beam path through excluded tissue")
        ax.legend(fontsize="small", loc="best")
        self.canvas.draw_idle()

    def _clicked(self, event) -> None:
        if event.inaxes is None or not self.cands:
            return
        ax = event.inaxes
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        d = [((k.dist_to_target - event.xdata) / (x1 - x0)) ** 2 + ((k.prop_intersect - event.ydata) / (y1 - y0)) ** 2
             for k in self.cands]
        self.select(self.cands[int(np.argmin(d))].idx)

    def _table_selected(self) -> None:
        items = self.table.selectedItems()
        if items:
            self.select(int(items[0].data(Qt.ItemDataRole.UserRole)))

    def select(self, idx: int | None) -> None:
        if idx is None or idx == self.current:
            return
        self.current = idx
        for r in range(self.table.rowCount()):
            it = self.table.item(r, 0)
            if it is not None and int(it.data(Qt.ItemDataRole.UserRole)) == idx:
                self.table.blockSignals(True)
                self.table.selectRow(r)
                self.table.blockSignals(False)
                break
        self._plot()
        self._update_side()

    def _update_side(self) -> None:
        k = next((q for q in self.cands if q.idx == self.current), None)
        if k is None:
            self.readout.setText("Click a candidate in the plot or the table.")
        else:
            crit = self.model.get("placement.heuristic.criterion_intersection")
            rows = [("Candidate", f"{k.idx}" + ("  (the heuristic's choice)" if k.idx == self.choice else "")),
                    ("Transducer [voxel]", "[" + " ".join(f"{v:.0f}" for v in k.trans_vox) + "]"),
                    ("Target [voxel]", "[" + " ".join(f"{v:.0f}" for v in k.targ_vox) + "]")]
            rows += [(label, fmt.format(getattr(k, key))) for key, label, fmt in COLUMNS[1:]]
            if not passes_intersection(k, crit):
                rows.append(("Note", f"above criterion_intersection ({crit}); the heuristic would skip it "
                                     "unless it widens the criterion"))
            self.readout.setText("<table>" + "".join(
                f"<tr><td>{a}</td><td style='padding-left:8px'><b>{b}</b></td></tr>" for a, b in rows)
                + "</table>")
        is_mode = self.model.get("placement.mode") == "heuristic"
        busy = self.follower.waiting
        self.place_btn.setEnabled(is_mode and k is not None and not busy)
        self.auto_btn.setEnabled(is_mode and self.chosen() is not None and not busy)
        cfg = self.chosen()
        self.status.setText(
            (f"placement.heuristic.candidate = {cfg}" if cfg is not None else
             "placement.heuristic.candidate = [] (the heuristic's choice)")
            + ("" if is_mode else "<br>placement.mode is not 'heuristic': placing is off."))

    # ── runs ────────────────────────────────────────────────────────────
    def _run(self, candidate: int | None, recompute: bool) -> None:
        if self.model.get("placement.mode") != "heuristic":
            return
        self.model.set("placement.heuristic.candidate", [] if candidate is None else int(candidate))
        overrides = launcher.placement_overrides("heuristic") + [
            ("io.overwrite_files", "always" if recompute else "never")]   # recompute or reuse tpars
        how = run_partial(self.window, TASK, "heuristic", "the heuristic placement", overrides,
                          launcher.PLACEMENT_ONLY)
        self.follower.start(how)
        self._update_side()

    def place(self) -> None:
        self._run(self.current, recompute=False)

    def use_heuristic(self) -> None:
        self._run(None, recompute=False)

    def compute(self) -> None:
        self._run(None, recompute=True)

    def _run_done(self, ok) -> None:
        if ok is not None:
            self.window.statusBar().showMessage(
                "Heuristic placement finished: the Viewer shows it." if ok else
                "Heuristic placement failed: see the Run tab's log.", 10000)
        if self.isVisible():
            self.refresh()

    def _viewer(self) -> None:
        from .main_window import VIEWER
        self.window.show_special(VIEWER)
