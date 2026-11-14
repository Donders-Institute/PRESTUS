"""Results tab: acoustic/thermal NIfTI viewer, the self-contained HTML report, the
output CSV, and a side-by-side comparison of every run in the subject's output folder.

One output folder holds several runs, told apart by their output affix (study
conditions, uncertainty variants, async per-transducer sub-runs and joint result,
sequential steps, ...). results.discover_runs lists them; the Run selector picks
one, and every sub-tab shows that run's files. The open config's own run is
selected by default (for async coupling: the joint result).

Showing the tab rescans the folder but reloads nothing unless the runs' files changed
(*Refresh* forces it). Only the visible sub-tab loads its files; the others load when
they are opened, so the T1, the maps and the HTML report are read on request."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PyQt6.QtCore import QUrl
from PyQt6.QtGui import QDesktopServices
from PyQt6.QtWidgets import (QComboBox, QHBoxLayout, QLabel, QPushButton, QTableWidget,
                             QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget)

from .. import placement as pl
from .. import results as rs
from ..model import is_empty
from ..outputs import dir_output
from .empty_state import EmptyStack
from .heat_view import HeatPanel
from .ortho_view import OrthoViewer, cached_volume

# Marker colours per transducer (1, 2, ...), drawn on the result maps.
TX_COLOURS = ("#ffffff", "#5ad1ff", "#9dff8a", "#ffd35a")


class _OverlayPanel(QWidget):
    """A viewer with a dropdown of the maps the selected run has."""

    def __init__(self, labels: tuple[str, ...], kind: str = "", set_t1=None,
                 parent: QWidget | None = None):
        super().__init__(parent)
        self._labels = labels
        self._kind = kind                   # 'acoustic' / 'thermal', for the empty state
        self._set_t1 = set_t1               # opens the T1 setting (MainWindow.goto)
        self._run: rs.Run | None = None

        v = QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 0)
        tb = QHBoxLayout()
        self.map_label = QLabel("Map")
        tb.addWidget(self.map_label)
        self.combo = QComboBox()
        self.combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
        self.combo.currentTextChanged.connect(self._load_overlay)
        tb.addWidget(self.combo)
        self.peak = QPushButton("Go to peak")
        self.peak.setToolTip("Centre the slices on the maximum of the shown map")
        self.peak.clicked.connect(lambda: self.viewer.centre_on(self.viewer.overlay_peak()))
        tb.addWidget(self.peak)
        tb.addStretch(1)
        self.note = QLabel()
        self.note.setObjectName("CardSub")
        tb.addWidget(self.note)
        v.addLayout(tb)
        self.viewer = OrthoViewer()
        button = getattr(getattr(self.viewer, "canvas", None), "export_button", None)
        if button is not None:                  # sub-XXX_<medium><affix>_<map>.svg …
            button.view = f"results_{kind or 'map'}"
            button.filename_hint = lambda: (f"{self._run.stem}_{self.combo.currentText()}"
                                            if self._run is not None else button.view)
        v.addWidget(self.viewer, 1)

    def load(self, run: rs.Run | None, anat: Path | None, markers) -> None:
        self._run = run
        self.viewer.set_anatomy(anat)
        self.viewer.set_markers(markers)
        current = self.combo.currentText()
        self.combo.blockSignals(True)
        self.combo.clear()
        available = [lbl for lbl in self._labels if run is not None and lbl in run.maps]
        self.combo.addItems(available)
        if current in available:
            self.combo.setCurrentText(current)
        self.combo.blockSignals(False)
        for w in (self.map_label, self.combo, self.peak):      # nothing to pick without maps
            w.setVisible(bool(available))
        self.note.setText("" if available or run is None else
                          f"This run wrote no {self._kind or 'such'} maps.")
        if run is None:
            self.viewer.set_empty_message("No run selected")
        elif anat is None:
            self.viewer.set_empty_message(
                "No T1 image to show the maps on",
                "The maps are drawn over the subject's T1 image; set path.anat and "
                "path.t1_pattern on Setup ▸ Paths.",
                "Set T1 image…" if self._set_t1 is not None else None, self._set_t1)
        self._load_overlay(self.combo.currentText(), centre=True)

    def _load_overlay(self, label: str, centre: bool = False) -> None:
        path = self._run.maps.get(label) if self._run is not None else None
        unit, zero_based, threshold = rs.MAP_DISPLAY.get(label, ("", False, 0))
        self.viewer.set_overlay(path, unit=unit, zero_based=zero_based, threshold_pct=threshold)
        if centre:
            self.viewer.when_loaded(lambda: self.viewer.centre_on(self.viewer.overlay_peak()))


class ResultsPage(QWidget):
    """Run selector over Acoustic / Thermal / Heating / Report / Data Table / Compare sub-tabs."""

    def __init__(self, window, parent: QWidget | None = None):
        super().__init__(parent)
        self.window = window
        self.runs: list[rs.Run] = []
        self._all_runs: list[rs.Run] = []   # including async per-transducer runs
        self._parts: list[rs.Run] = []      # selected run + its async per-transducer runs
        self._fingerprint = None            # folder state the page shows
        self._part: rs.Run | None = None    # run the sub-tabs show
        self._stale: set[QWidget] = set()   # sub-tabs not yet loaded for self._part
        self._marker_cache: tuple | None = None

        v = QVBoxLayout(self)
        title = QLabel("Results")
        title.setObjectName("CardTitle")
        v.addWidget(title)

        tb = QHBoxLayout()
        refresh = QPushButton("Refresh")
        refresh.setToolTip("Rescan the output folder and reload the shown files")
        refresh.clicked.connect(lambda: self.refresh(force=True))
        tb.addWidget(refresh)
        self.run_label = QLabel("Run")
        tb.addWidget(self.run_label)
        self.run_combo = QComboBox()
        self.run_combo.setMinimumWidth(320)
        self.run_combo.setToolTip("Runs found in this subject's output folder, by output affix")
        self.run_combo.currentIndexChanged.connect(self._show_run)
        tb.addWidget(self.run_combo)
        self.parts = QComboBox()
        self.parts.setToolTip("Async coupling: the joint result or one transducer's own field")
        self.parts.currentIndexChanged.connect(self._show_part)
        tb.addWidget(self.parts)
        open_html = QPushButton("Open HTML report")
        open_html.setToolTip("Open the selected report in the system's web browser")
        open_html.clicked.connect(self._open_html)
        tb.addWidget(open_html)
        open_folder = QPushButton("Open output folder")
        open_folder.setToolTip("Show this subject's output folder in the file manager")
        open_folder.clicked.connect(self._open_folder)
        tb.addWidget(open_folder)
        # run selector and open buttons do nothing until there are runs: hidden till then
        self._run_controls = (self.run_label, self.run_combo, open_html, open_folder)
        tb.addStretch(1)
        self.status = QLabel("No output loaded.")
        self.status.setObjectName("CardSub")
        tb.addWidget(self.status)
        v.addLayout(tb)

        self.subtabs = QTabWidget()
        self.stack = EmptyStack(self.subtabs)
        v.addWidget(self.stack, 1)
        acoustic = tuple(k for k in rs.MAPS if k not in rs.THERMAL_MAPS)
        set_t1 = (lambda: window.goto("path.anat")) if hasattr(window, "goto") else None
        self.acoustic = _OverlayPanel(acoustic, "acoustic", set_t1)
        self.subtabs.addTab(self.acoustic, "Acoustic")
        self.thermal = _OverlayPanel(rs.THERMAL_MAPS, "thermal", set_t1)
        self.subtabs.addTab(self.thermal, "Thermal")
        self.heat = HeatPanel()
        self.subtabs.addTab(self.heat, "Heating")
        self.report_tab, self.web, self.report_combo = self._build_report_tab()
        self.subtabs.addTab(self.report_tab, "Report")
        self.table = QTableWidget()
        self.subtabs.addTab(self.table, "Data Table")
        self.compare = QTableWidget()
        self.compare.setToolTip("One row per run; double-click a row to show that run")
        self.compare.cellDoubleClicked.connect(lambda row, _c: self.run_combo.setCurrentIndex(row))
        self.subtabs.addTab(self.compare, "Compare runs")
        self.subtabs.currentChanged.connect(lambda _i: self._load_visible())

    def _build_report_tab(self):
        w = QWidget()
        vl = QVBoxLayout(w)
        vl.setContentsMargins(0, 0, 0, 0)
        row = QHBoxLayout()
        row.addWidget(QLabel("Report"))
        combo = QComboBox()
        combo.setToolTip("The selected run's report, or any other report in the folder "
                         "(e.g. the sequential or multi-ISPPA summary)")
        combo.currentIndexChanged.connect(self._load_report)
        row.addWidget(combo, 1)
        vl.addLayout(row)
        try:
            from PyQt6.QtWebEngineWidgets import QWebEngineView
            web = QWebEngineView()
            vl.addWidget(web, 1)
            return w, web, combo
        except ImportError:
            note = QLabel("PyQt6-WebEngine not installed — use 'Open HTML report' to view it "
                          "in a browser instead.")
            note.setWordWrap(True)
            vl.addWidget(note)
            return w, None, combo

    # ── model-derived paths ─────────────────────────────────────────────
    def _base(self) -> Path | None:
        base = dir_output(self.window.model)
        return base if base is not None and base.is_dir() else None

    def _anat(self) -> Path | None:
        m = self.window.model
        anat, pattern, sid = m.get("path.anat"), m.get("path.t1_pattern"), m.get("subject_id")
        if is_empty(anat) or not isinstance(sid, (int, float)) or isinstance(sid, bool):
            return None
        return pl.find_t1(anat, pattern or "", int(sid))

    def refresh(self, force: bool = False) -> None:
        """Rescan the output folder; rebuild the page only when the runs or their files
        changed since the last scan (or when forced)."""
        base = self._base()
        if base is None:
            self._fingerprint = None
            self.runs, self._all_runs = [], []
            self._show_empty(dir_output(self.window.model))
            return
        all_runs = rs.runs_for(self.window.model)
        target = rs.default_run(all_runs, self.window.model)
        target_affix = all_runs[target].affix if target >= 0 else ""
        m = self.window.model
        config_markers = repr([(m.get("transducer.trans_pos", ti), m.get("transducer.focus_pos", ti))
                               for ti in range(m.n_transducers)]
                              + [m.get("placement.manual.coordinate_space")])
        fingerprint = (str(base), target_affix, tuple(_run_state(r) for r in all_runs),
                       tuple(p.name for p in rs.other_reports(base)), str(self._anat()),
                       config_markers)
        if not force and fingerprint == self._fingerprint:
            return
        same_target = self._fingerprint is not None and self._fingerprint[:2] == fingerprint[:2]
        shown = self._current()
        self._fingerprint = fingerprint
        self._marker_cache = None
        self._all_runs = all_runs
        self.runs = [r for r in all_runs if r.kind != "transducer"]
        self.run_combo.blockSignals(True)
        self.run_combo.clear()
        self.run_combo.addItems([r.label for r in self.runs])
        self.run_combo.blockSignals(False)
        self._fill_compare()
        if not self.runs:
            self._show_run(-1)
            self._show_empty(base)
            return
        self.stack.show_content()
        for w in self._run_controls:
            w.show()
        # Keep the run on show when only its files changed; else the config's own run.
        keep = shown.affix if same_target and shown is not None else target_affix
        idx = next((i for i, r in enumerate(self.runs) if r.affix == keep),
                   next((i for i, r in enumerate(self.runs) if r.affix == target_affix), 0))
        self.run_combo.blockSignals(True)
        self.run_combo.setCurrentIndex(idx)
        self.run_combo.blockSignals(False)
        self._show_run(idx)
        self.status.setText(f"{len(self.runs)} run(s) in {base}")

    def _show_empty(self, base: Path | None) -> None:
        """No runs to show: one panel that says so and leads to the Run tab."""
        sid = self.window.model.get("subject_id")
        sub = f"sub-{int(sid):03d}" if isinstance(sid, (int, float)) and not isinstance(sid, bool) else "this subject"
        where = f" PRESTUS writes them to {base}." if base is not None else \
            " Set path.sim (where PRESTUS writes its outputs) on Setup ▸ Paths."
        self.stack.show_empty(f"No results for {sub} yet",
                              "Run the pipeline on the Run tab; its maps, report and tables "
                              "appear here." + where, "Go to Run", self._goto_run)
        for w in self._run_controls:
            w.hide()
        self.parts.hide()
        self.status.setText("")

    def _goto_run(self) -> None:
        from .main_window import RUN
        self.window.show_special(RUN)

    # ── one run ─────────────────────────────────────────────────────────
    def _current(self) -> rs.Run | None:
        i = self.run_combo.currentIndex()
        return self.runs[i] if 0 <= i < len(self.runs) else None

    def _show_run(self, _index: int = 0) -> None:
        run = self._current()
        self.parts.blockSignals(True)
        self.parts.clear()
        self._parts = []
        if run is not None:
            self._parts = [run] + rs.related_runs(self._all_runs, run)
            if len(self._parts) > 1:
                self.parts.addItems(["Joint (both transducers)"] +
                                    [f"Transducer {r.transducer} alone" for r in self._parts[1:]])
        self.parts.setVisible(len(self._parts) > 1)
        self.parts.blockSignals(False)
        self._show_part(0)

    def _show_part(self, index: int) -> None:
        self._part = self._parts[index] if 0 <= index < len(self._parts) else self._current()
        self._stale = {self.acoustic, self.thermal, self.heat, self.report_tab}
        self._fill_reports(self._part)
        self._load_table(self._part)
        self._load_visible()

    def _thermal_run(self) -> rs.Run | None:
        """Async per-transducer fields have no thermal stage: show the joint run's."""
        run = self._part
        return self._current() if run is not None and run.kind == "transducer" else run

    def _load_visible(self) -> None:
        """Load the shown sub-tab's files if it has not shown the current run yet."""
        w = self.subtabs.currentWidget()
        if w not in self._stale:
            return
        self._stale.discard(w)
        if w is self.acoustic or w is self.thermal:
            anat = self._anat()
            markers = self._markers(self._part, anat)
            if w is self.acoustic:
                self.acoustic.load(self._part, anat, markers)
            else:
                self.thermal.load(self._thermal_run(), anat, markers)
        elif w is self.heat:
            self.heat.load(self._thermal_run())
        elif w is self.report_tab:
            self._load_report()

    def _markers(self, run: rs.Run | None, anat: Path | None):
        """Transducer (square) and focus (x) of every transducer of the run, as its
        placement stage recorded them, else as the open config sets them (manual voxels)."""
        if run is None or anat is None:
            return []
        key = (run.affix, str(anat))
        if self._marker_cache is not None and self._marker_cache[0] == key:
            return self._marker_cache[1]
        try:
            t1 = cached_volume(anat) or pl.Volume(anat)     # the header is enough here
        except Exception:
            return []
        positions = rs.placement_positions(run)
        if not positions and self.window.model.get("placement.manual.coordinate_space") in (None, "", "voxel"):
            positions = [(self.window.model.get("transducer.trans_pos", ti),
                          self.window.model.get("transducer.focus_pos", ti))
                         for ti in range(self.window.model.n_transducers)]
            if run.kind == "transducer" and run.transducer and run.transducer <= len(positions):
                positions = [positions[run.transducer - 1]]
        offset = (run.transducer - 1) if run.kind == "transducer" and run.transducer else 0
        markers = []
        for i, (trans, focus) in enumerate(positions):
            colour = TX_COLOURS[(i + offset) % len(TX_COLOURS)]
            for pos, marker in ((trans, "s"), (focus, "x")):
                if isinstance(pos, list) and len(pos) == 3 and all(isinstance(c, (int, float)) for c in pos):
                    markers.append((pl.prestus_voxel_to_ras(t1, pos), marker, colour))
        self._marker_cache = (key, markers)
        return markers

    def _fill_reports(self, run: rs.Run | None) -> None:
        base = self._base()
        reports = rs.other_reports(base) if base is not None else []
        self.report_combo.blockSignals(True)
        self.report_combo.clear()
        for p in reports:
            self.report_combo.addItem(p.name, str(p))
        if run is not None and run.report is not None:
            self.report_combo.setCurrentIndex(max(0, self.report_combo.findData(str(run.report))))
        self.report_combo.blockSignals(False)
        self._load_report()

    def _load_report(self, _index: int = 0) -> None:
        if self.subtabs.currentWidget() is not self.report_tab:
            self._stale.add(self.report_tab)       # load when the Report tab is opened
            return
        path = self.report_combo.currentData()
        if path and self.web is not None:
            self.web.setUrl(QUrl.fromLocalFile(path))

    def _load_table(self, run: rs.Run | None) -> None:
        metrics = rs.read_metrics(run.csv if run is not None else None)
        self.table.clear()
        self.table.setColumnCount(2)
        self.table.setRowCount(len(metrics))
        self.table.setHorizontalHeaderLabels(["Metric", "Value"])
        for ri, (k, value) in enumerate(metrics.items()):
            self.table.setItem(ri, 0, QTableWidgetItem(k))
            self.table.setItem(ri, 1, QTableWidgetItem(value))
        self.table.resizeColumnsToContents()

    def _fill_compare(self) -> None:
        cols = rs.COMPARE_COLUMNS
        self.compare.clear()
        self.compare.setColumnCount(len(cols) + 1)
        self.compare.setRowCount(len(self.runs))
        self.compare.setHorizontalHeaderLabels(["Run"] + [label for _, label in cols])
        for ri, run in enumerate(self.runs):
            metrics = rs.read_metrics(run.csv)
            self.compare.setItem(ri, 0, QTableWidgetItem(run.label))
            for ci, (key, _) in enumerate(cols, start=1):
                self.compare.setItem(ri, ci, QTableWidgetItem(_fmt(metrics.get(key, ""))))
        self.compare.resizeColumnsToContents()

    def _open_html(self) -> None:
        path = self.report_combo.currentData()
        if path:
            _open_in_os(Path(path))

    def _open_folder(self) -> None:
        base = self._base()
        if base is not None:
            _open_in_os(base)


def _run_state(run: rs.Run) -> tuple:
    """A run's files with their modification times, to tell whether anything changed."""
    files = [run.csv, run.report, run.heat, *run.maps.values()]
    state = []
    for p in files:
        try:
            state.append((p.name, p.stat().st_mtime_ns))
        except (AttributeError, OSError):
            state.append(None)
    return (run.medium, run.affix, tuple(state))


def _fmt(value: str) -> str:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return value or ""
    return "" if np.isnan(x) else f"{x:.4g}"


def _open_in_os(path: Path) -> None:
    """Open a file or folder with the system's default app (macOS, Windows and Linux)."""
    QDesktopServices.openUrl(QUrl.fromLocalFile(str(path)))
