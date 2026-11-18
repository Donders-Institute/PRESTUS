"""Results tab: a reporting overview of the run (ITRUSST items and non-significant-risk
check, itrusst_report.py), the acoustic/thermal NIfTI maps (Spatial), the thermal time
courses, the self-contained HTML report, the output CSV, and a side-by-side comparison of
every run in the subject's output folder. The Overview opens first.

One output folder holds several runs, told apart by their output affix (study
conditions, uncertainty variants, async per-transducer sub-runs and joint result,
sequential steps, ...). results.discover_runs lists them; the Run selector picks
one, and every sub-tab shows that run's files. The open config's own run is
selected by default (for async coupling: the joint result).

Showing the tab rescans the folder but reloads nothing unless the runs' files changed
(*Refresh* forces it). Only the visible sub-tab loads its files; the others load when
they are opened, so the T1, the maps and the HTML report are read on request."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from html import escape

from PyQt6.QtCore import QUrl, pyqtSignal
from PyQt6.QtGui import QDesktopServices, QGuiApplication
from PyQt6.QtWidgets import (QComboBox, QFileDialog, QHBoxLayout, QLabel, QMessageBox, QPushButton,
                             QTableWidget,
                             QTableWidgetItem, QTabWidget, QTextBrowser, QVBoxLayout, QWidget)

from ..labels import plural
from .. import focal_widths as fwm
from .. import followup as fu
from .. import itrusst_report as ir
from .. import placement as pl
from .. import results as rs
from ..model import is_empty
from ..outputs import dir_output
from .tabbar import EmphasisTabBar, fit_to_width
from . import style
from .empty_state import EmptyStack
from .heat_view import HeatPanel
from .ortho_view import OrthoViewer, cached_volume
from .pose_table import PoseTable

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


class OverviewPanel(QWidget):
    """The run's parameters as ITRUSST asks TUS studies to report them, with the
    non-significant-risk check; each row names its source (itrusst_report.py).

    Focal widths of runs whose CSV has no focal_fwhm_* columns are measured on their T1w
    intensity map (focal_widths.py) on a worker thread; the rows read 'measuring…' until
    the result arrives, and the overview is then rebuilt."""

    STATUS = {"within": ("ok", "within"), "near": ("warn", "near limit"), "above": ("bad", "above")}
    _widths_done = pyqtSignal(object, object)          # (token, {stem: FocalWidths | None})
    _executor = ThreadPoolExecutor(max_workers=1, thread_name_prefix="nexus-focal-widths")

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.report: ir.Report | None = None
        self._args = None                                # (model, run, runs) last shown
        self._widths: dict = {}
        self._token = 0
        self._widths_done.connect(self._on_widths)
        v = QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 0)
        bar = QHBoxLayout()
        self.summary = QLabel()
        self.summary.setWordWrap(True)
        bar.addWidget(self.summary, 1)
        self.copy = QPushButton("Copy as Markdown")
        self.copy.setToolTip("Copy the overview as Markdown tables, e.g. for a methods section or a lab notebook")
        self.copy.clicked.connect(self._copy)
        bar.addWidget(self.copy)
        v.addLayout(bar)
        self.view = QTextBrowser()
        self.view.setOpenLinks(False)
        v.addWidget(self.view, 1)

    def load(self, model, run: rs.Run | None, runs: list[rs.Run]) -> None:
        self._token += 1
        self._args = (model, run, runs)
        self._widths = {}
        todo = []
        if run is not None:
            for r in ir.width_runs(run, runs):
                inputs = fwm.inputs_for(r)
                axis = tuple(float(a) for a in (*inputs[2], *inputs[3]))
                hit = fwm.cached(inputs[0], inputs[1], axis)
                self._widths[r.stem] = hit if hit is not None else ir.PENDING
                if hit is None:
                    todo.append((r.stem, inputs))
        if todo:
            token = self._token
            fut = self._executor.submit(_measure_widths, todo)
            fut.add_done_callback(lambda f: self._widths_done.emit(token, f))
        self._show()

    def _on_widths(self, token, fut) -> None:
        if token != self._token or self._args is None:
            return                                       # another run is shown by now
        try:
            found = fut.result()
        except Exception:                                # noqa: BLE001 — shown as not available
            found = {}
        for stem, value in list(self._widths.items()):
            if value == ir.PENDING:
                self._widths[stem] = found.get(stem)
        self._show()

    def _show(self) -> None:
        model, run, runs = self._args
        self.report = ir.build(model, run, runs, self._widths) if run is not None else None
        self.copy.setEnabled(self.report is not None)
        if self.report is None:
            self.summary.setText("No run selected.")
            self.view.clear()
            return
        from .readout import status
        r = self.report
        if r.verdict == "within":
            head = status("Within ITRUSST non-significant-risk limits", "ok")
        elif r.verdict == "review":
            head = status(f"{plural(len(r.above), 'metric')} above their ITRUSST NSR reference level: "
                          f"{escape(', '.join(r.above))}. Inspect NSR adherence.", "bad")
        else:
            head = status("No ITRUSST-limited metric in this run's outputs", "muted")
        self.summary.setText(f"<b>{escape(r.run_label)}</b> ({escape(r.medium)})<br>{head}")   # verdict on its own line
        self.view.setHtml(self.html(r))

    def html(self, r: ir.Report) -> str:
        from .readout import status
        c = style.active()
        out = [f"<style>td {{ padding: 2px 10px 2px 0; }} th {{ text-align: left; color: {c['ink3']}; "
               f"font-weight: normal; padding: 2px 10px 2px 0; }} h3 {{ margin: 12px 0 4px 0; }}</style>"]
        for sec in r.sections:
            out.append(f"<h3>{escape(sec.title)}</h3>")
            if sec.note:
                out.append(f"<p style='color:{c['ink3']}'>{escape(sec.note)}</p>")
            if not sec.rows:
                continue
            limits = any(row.limit for row in sec.rows)
            head = "<tr><th>Item</th><th>Value</th>" + ("<th>Limit</th><th>Status</th>" if limits else "") \
                + "<th>Source</th></tr>"
            out.append(f"<table cellspacing='0' width='100%'>{head}")
            for row in sec.rows:
                value = escape(row.value)
                if row.value == ir.NA:
                    value = f"<span style='color:{c['ink4']}'>{value}</span>"
                cells = f"<td>{escape(row.label)}</td><td><b>{value}</b></td>"
                if limits:
                    st = self.STATUS.get(row.status)
                    cells += f"<td>{escape(row.limit)}</td><td>{status(st[1], st[0]) if st else ''}</td>"
                cells += f"<td style='color:{c['ink3']}'>{escape(row.source)}</td>"
                out.append(f"<tr>{cells}</tr>")
            out.append("</table>")
        return "".join(out)

    def _copy(self) -> None:
        if self.report is not None:
            QGuiApplication.clipboard().setText(ir.to_markdown(self.report))


def _measure_widths(todo) -> dict:
    """{stem: FocalWidths or None} for (stem, focal_widths.inputs_for) pairs (worker thread)."""
    out = {}
    for stem, (path, mask, trans, focus) in todo:
        try:
            out[stem] = fwm.measure(path, mask, trans, focus)
        except Exception:                                # noqa: BLE001 — unreadable map
            out[stem] = None
    return out


class ResultsPage(QWidget):
    """Run selector over Overview / Spatial (Acoustic, Thermal) / Heating / Report /
    Data table / Poses (time-reversal pose search only) / Compare sub-tabs.

    "Follow-up…" starts a new config from the selected run (followup_dialog); the run header
    shows its lineage (follow-up configs in the study folder, results.lineage_links).
    ``followup_created(path)`` is emitted with the saved follow-up config."""

    followup_created = pyqtSignal(object)

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
        self.followup_btn = QPushButton("Follow-up…")
        self.followup_btn.setToolTip("Create a new config that starts from this run's results: "
                                     "another intensity, timing or placement, a time-reversal pose, "
                                     "or a sequential run with heat carry-over")
        self.followup_btn.clicked.connect(lambda: self.open_followup())
        tb.addWidget(self.followup_btn)
        # run selector and open buttons do nothing until there are runs: hidden till then
        self._run_controls = (self.run_label, self.run_combo, open_html, open_folder, self.followup_btn)
        tb.addStretch(1)
        self.status = QLabel("No output loaded.")
        self.status.setObjectName("CardSub")
        tb.addWidget(self.status)
        v.addLayout(tb)
        self.lineage = QLabel()
        self.lineage.setObjectName("CardSub")
        self.lineage.setWordWrap(True)
        self.lineage.hide()
        v.addWidget(self.lineage)

        self.subtabs = QTabWidget()
        self.subtabs.setTabBar(EmphasisTabBar(self.subtabs))
        fit_to_width(self.subtabs.tabBar())
        self.stack = EmptyStack(self.subtabs)
        v.addWidget(self.stack, 1)
        self.overview = OverviewPanel()
        self.subtabs.addTab(self.overview, "Overview")
        acoustic = tuple(k for k in rs.MAPS if k not in rs.THERMAL_MAPS)
        set_t1 = (lambda: window.goto("path.anat")) if hasattr(window, "goto") else None
        self.spatial = QTabWidget()                 # the maps on the T1
        self.spatial.setTabBar(EmphasisTabBar(self.spatial))
        self.spatial.setDocumentMode(True)
        fit_to_width(self.spatial.tabBar())
        self.acoustic = _OverlayPanel(acoustic, "acoustic", set_t1)
        self.spatial.addTab(self.acoustic, "Acoustic")
        self.thermal = _OverlayPanel(rs.THERMAL_MAPS, "thermal", set_t1)
        self.spatial.addTab(self.thermal, "Thermal")
        self.spatial.currentChanged.connect(lambda _i: self._load_visible())
        self.subtabs.addTab(self.spatial, "Spatial")
        self.heat = HeatPanel(filename_hint=lambda: (f"{self._run.stem}_heating" if getattr(self, "_run", None)
                                                    is not None else "heating"))
        self.subtabs.addTab(self.heat, "Heating")
        self.report_tab, self.web, self.report_combo = self._build_report_tab()
        self.subtabs.addTab(self.report_tab, "Report")
        self.table = QTableWidget()
        self.subtabs.addTab(self.table, "Data table")
        self.poses = PoseTable()
        self.poses.use_pose.connect(lambda row: self.open_followup("pose_from_time_reversal", row))
        self.subtabs.addTab(self.poses, "Poses")
        self.subtabs.setTabVisible(self.subtabs.indexOf(self.poses), False)
        self.subtabs.setTabToolTip(self.subtabs.indexOf(self.poses),
                                   "Poses scored by the time-reversal pose search")
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
        combo.setToolTip("The selected run's report, or any other report of the subject (e.g. the "
                         "placement, sequential, or multi-ISPPA report), or a group report in path.sim")
        combo.currentIndexChanged.connect(self._load_report)
        row.addWidget(combo, 1)
        self.pdf_btn = QPushButton("Export report as PDF…")
        self.pdf_btn.setToolTip("Save the shown report as a PDF")
        self.pdf_btn.clicked.connect(self._export_pdf)
        self.pdf_btn.setEnabled(False)
        row.addWidget(self.pdf_btn)
        vl.addLayout(row)
        self._report_layout = vl
        return w, None, combo

    def _web(self):
        """The report view, created when the Report tab first shows a report (QtWebEngine is
        imported before the QApplication in app.py, as Qt requires; the view itself is only
        built on demand). None without PyQt6-WebEngine."""
        if self.web is None and not getattr(self, "_no_web", False):
            try:
                from PyQt6.QtWebEngineWidgets import QWebEngineView
            except ImportError:
                self._no_web = True
                note = QLabel("PyQt6-WebEngine not installed — use 'Open HTML report' to view it "
                              "in a browser instead.")
                note.setWordWrap(True)
                self._report_layout.addWidget(note)
                return None
            self.web = QWebEngineView()
            self._report_layout.addWidget(self.web, 1)
        return self.web

    def _export_pdf(self) -> None:
        web = self.web
        src = self.report_combo.currentData()
        if web is None or not src:
            return
        start = str(Path(src).with_suffix(".pdf"))
        chosen, _ = QFileDialog.getSaveFileName(self, "Export report as PDF", start, "PDF (*.pdf)")
        if not chosen:
            return
        target = Path(chosen)

        def done(data) -> None:
            from ..fsio import write_bytes_atomic
            try:
                if not data:
                    raise OSError("the report could not be printed")
                write_bytes_atomic(target, bytes(data))
            except OSError as exc:
                QMessageBox.warning(self, "PDF not saved", f"Could not save {target}:\n{exc}")
                return
            self.status.setText(f"Saved {target.name}")
        web.page().printToPdf(done)

    def show_view(self, widget: QWidget) -> None:
        """Bring a sub-tab forward, also one inside Spatial (Acoustic, Thermal)."""
        if self.spatial.indexOf(widget) >= 0:
            self.spatial.setCurrentWidget(widget)
            self.subtabs.setCurrentWidget(self.spatial)
        else:
            self.subtabs.setCurrentWidget(widget)
        self._load_visible()

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
        self.status.setText(f"{plural(len(self.runs), 'run')} in {base}")

    def _show_empty(self, base: Path | None) -> None:
        """No runs to show: one panel that says so and leads to the Run tab."""
        sid = self.window.model.get("subject_id")
        sub = f"sub-{int(sid):03d}" if isinstance(sid, (int, float)) and not isinstance(sid, bool) else "this subject"
        where = f" PRESTUS writes them to {base}." if base is not None else \
            " Set path.sim (where PRESTUS writes its outputs) on Setup ▸ Paths."
        self.stack.show_empty(f"No results for {sub} yet",
                              "Run the pipeline on the Run tab; its maps, report, and tables "
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
        self._stale = {self.overview, self.acoustic, self.thermal, self.heat, self.report_tab,
                       self.poses}
        self._fill_reports(self._part)
        self._load_table(self._part)
        self._show_lineage(self._current())
        tables = fu.pose_tables(self._part) if self._part is not None else []
        self.subtabs.setTabVisible(self.subtabs.indexOf(self.poses), bool(tables))
        self._pose_tables = tables
        self.followup_btn.setEnabled(self._part is not None and self.window.model.path is not None)
        self._load_visible()

    # ── follow-ups and lineage ──────────────────────────────────────────
    def _show_lineage(self, run: rs.Run | None) -> None:
        """'Follow-up of …' / 'Follow-ups: …' for the selected run, from the lineage headers of
        the configs in the study folder."""
        path = self.window.model.path
        if run is None or path is None:
            self.lineage.hide()
            return
        links = rs.lineage_links(path.parent)
        up, down = rs.lineage_for(run, links, path)
        parts = []
        if up is not None:
            kind = fu.TYPES[up.type].label if up.type in fu.TYPES else up.type
            parent = Path(up.parent_config).name if up.parent_config else "?"
            label = rs.describe_affix(up.parent_affix)[0]
            parts.append(f"Follow-up ({escape(kind)}) of run <b>{escape(label)}</b> in {escape(parent)}, "
                         f"config {escape(up.child_config.name)}")
        if down:
            items = ", ".join(f"{escape(rs.describe_affix(d.child_affix)[0])} ({escape(d.child_config.name)})"
                              for d in down)
            parts.append(f"Follow-ups of this run: {items}")
        self.lineage.setText("<br>".join(parts))      # parent and children on separate lines
        self.lineage.setVisible(bool(parts))

    def open_followup(self, ftype: str | None = None, pose: dict | None = None):
        """The follow-up dialog for the selected run (or async part)."""
        run = self._part or self._current()
        if run is None:
            return None
        if self.window.model.path is None:
            QMessageBox.information(self, "Save first", "Save this config first: a follow-up refers to it.")
            return None
        from .followup_dialog import FollowupDialog
        dlg = FollowupDialog(self.window, run, ftype=ftype, pose=pose, parent=self)
        dlg.created.connect(self.followup_created.emit)
        dlg.created.connect(lambda _p: self._show_lineage(self._current()))
        dlg.open()
        return dlg

    def _thermal_run(self) -> rs.Run | None:
        """Async per-transducer fields have no thermal stage: show the joint run's."""
        run = self._part
        return self._current() if run is not None and run.kind == "transducer" else run

    def _load_visible(self) -> None:
        """Load the shown sub-tab's files if it has not shown the current run yet."""
        w = self.subtabs.currentWidget()
        if w is self.spatial:
            w = self.spatial.currentWidget()
        if w not in self._stale:
            return
        self._stale.discard(w)
        if w is self.overview:
            self.overview.load(self.window.model, self._part, self._all_runs)
        elif w is self.acoustic or w is self.thermal:
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
        elif w is self.poses:
            self.poses.load(getattr(self, "_pose_tables", []))

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
            ti = i + offset                     # hover text: the device the open config names
            tip = (self.window.model.transducer_tooltip(ti) if ti < self.window.model.n_transducers
                   else f"Transducer {ti + 1}")
            for pos, marker in ((trans, "s"), (focus, "x")):
                if isinstance(pos, list) and len(pos) == 3 and all(isinstance(c, (int, float)) for c in pos):
                    markers.append((pl.prestus_voxel_to_ras(t1, pos), marker, colour,
                                    tip if marker == "s" else f"{tip} · focus"))
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
        web = self._web() if path else self.web
        if path and web is not None:
            web.setUrl(QUrl.fromLocalFile(path))
        self.pdf_btn.setEnabled(bool(path) and web is not None)

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
