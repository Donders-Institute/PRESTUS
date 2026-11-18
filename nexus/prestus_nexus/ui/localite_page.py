"""Placement tab ▸ Localite: detect the recorded series of a session, name them, write the
ingested markers, and use them in the study (the steps of neuronav_gui, run in MATLAB
through the bridge with the study config's own settings).

Batch ingestion of many sessions runs as a script on the Run tab; the MNI view maps the
batch summary with each subject's SimNIBS warp.
"""

from __future__ import annotations

import html

import datetime as _dt
from pathlib import Path

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import (QCheckBox, QComboBox, QDialog, QFormLayout, QGridLayout,
                             QHBoxLayout, QHeaderView, QLabel, QLineEdit, QPlainTextEdit,
                             QPushButton, QSplitter, QTableWidget, QTableWidgetItem, QTabWidget,
                             QToolButton, QVBoxLayout, QWidget)

from ..labels import plural
from .. import localite as lc
from ..bridge import BridgeError, as_floats, as_list
from . import style
from .empty_state import EmptyStack
from .matlab_status import MatlabStatus
from .tabbar import EmphasisTabBar

try:
    from .mpl_theme import ThemedCanvas as FigureCanvasQTAgg   # drawn in the Nexus theme
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D
except ImportError:                                    # pragma: no cover
    FigureCanvasQTAgg = None

PALETTE = style.SERIES          # colour-blind-safe categorical set (series, MNI groups)
SERIES_COLS = ["Write", "#", "Coil", "Target name", "Name from", "Colour", "Pulses",
               "Pos (mm)", "Pos SD", "Rot (rad)", "Rot SD", "Flag"]


def _item(text, editable: bool = False) -> QTableWidgetItem:
    it = QTableWidgetItem("" if text is None else str(text))
    if not editable:
        it.setFlags(it.flags() & ~Qt.ItemFlag.ItemIsEditable)
    return it


def _num(v, fmt: str = "{:.2f}") -> str:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return ""
    return "" if f != f else fmt.format(f)


class LocalitePage(QWidget):
    def __init__(self, window, parent: QWidget | None = None):
        super().__init__(parent)
        self.window = window
        self.result: dict | None = None        # last 'detect' answer
        self.snapshot = None                   # (sub, ses, settings) at detect time
        self.written: dict | None = None

        outer = QVBoxLayout(self)
        outer.setContentsMargins(12, 8, 12, 8)
        head = QHBoxLayout()
        intro = QLabel("Detect the recorded series of a session, check and name them, and write the "
                       "ingested markers (path.localite_post). The settings are the study config's; "
                       "edit them on the Settings tab.")
        intro.setObjectName("CardSub")
        intro.setWordWrap(True)
        head.addWidget(intro, 1)
        # MATLAB starts only on an action that needs it (detect, write, batch, MNI view):
        # showing this page reads folders only.
        bridge = window.bridge()
        self.matlab_status = MatlabStatus(bridge)
        head.addWidget(self.matlab_status)
        bridge.state_changed.connect(self._idle_status)     # after MatlabStatus._update
        self._idle_status(getattr(bridge, "state", ""))
        outer.addLayout(head)

        # one setting per line: the two folders always, the rest behind a 'more' link
        self.summary = QLabel()
        self.summary.setWordWrap(True)
        self.summary.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse |
                                             Qt.TextInteractionFlag.LinksAccessibleByMouse)
        self._summary_open = False
        self.summary.linkActivated.connect(self._toggle_summary)
        srow = QHBoxLayout()
        srow.addWidget(self.summary, 1)
        edit = QPushButton("Edit settings")
        edit.clicked.connect(lambda: window.goto("path.localite_raw"))
        srow.addWidget(edit)
        outer.addLayout(srow)

        row = QHBoxLayout()
        row.addWidget(QLabel("Subject"))
        self.sub = QComboBox(editable=True)
        self.sub.setMinimumWidth(110)
        self.sub.currentTextChanged.connect(self._sub_changed)
        row.addWidget(self.sub)
        row.addWidget(QLabel("Session"))
        self.ses = QComboBox(editable=True)
        self.ses.setMinimumWidth(90)
        self.ses.currentTextChanged.connect(lambda _t: self._inputs_changed())
        row.addWidget(self.ses)
        rescan = QToolButton(text="↻")
        rescan.setToolTip("List the subjects and sessions under the raw root again")
        rescan.clicked.connect(self._fill_subjects)
        row.addWidget(rescan)
        self.detect_btn = QPushButton("1. Detect series")
        self.detect_btn.setProperty("primary", True)
        self.detect_btn.clicked.connect(self.detect)
        row.addWidget(self.detect_btn)
        self.write_btn = QPushButton("2. Write outputs")
        self.write_btn.setEnabled(False)
        self.write_btn.clicked.connect(self.write)
        row.addWidget(self.write_btn)
        row.addStretch(1)
        self.use_btn = QPushButton("Use in this study")
        self.use_btn.setToolTip("Set placement.localite.session and .target to this session and the "
                                "written target(s)")
        self.use_btn.setEnabled(False)
        self.use_btn.clicked.connect(self.use_in_study)
        row.addWidget(self.use_btn)
        imp = QPushButton("Import into study…")
        imp.setToolTip("Study manifest rows from all ingested sessions and targets")
        imp.clicked.connect(self._import)
        row.addWidget(imp)
        batch = QPushButton("Ingest all sessions…")
        batch.setToolTip("Scan all subjects/sessions, map planned-target labels, and ingest many sessions")
        batch.clicked.connect(self.open_batch)
        row.addWidget(batch)
        mni = QPushButton("MNI view…")
        mni.setToolTip("Written targets of all subjects and sessions in MNI space (batch summary + SimNIBS warps)")
        mni.clicked.connect(self.open_mni)
        row.addWidget(mni)
        outer.addLayout(row)
        self.status = QLabel()
        self.status.setWordWrap(True)
        outer.addWidget(self.status)

        split = QSplitter(Qt.Orientation.Vertical)
        self.table = QTableWidget(0, len(SERIES_COLS))
        self.table.setHorizontalHeaderLabels(SERIES_COLS)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.table.itemChanged.connect(self._series_edited)
        split.addWidget(self.table)
        tabs = QTabWidget()
        tabs.setTabBar(EmphasisTabBar(tabs))               # the Nexus's tab dividers
        if FigureCanvasQTAgg is not None:
            self.fig_pos = Figure(figsize=(6, 4), layout="constrained")
            self.ax_pos = self.fig_pos.add_subplot(projection="3d")
            self.canvas_pos = FigureCanvasQTAgg(self.fig_pos, export_name="localite_positions")
            self.pos_stack = EmptyStack(self.canvas_pos)
            tabs.addTab(self.pos_stack, "Positions (3-D)")
            self.fig_dev = Figure(figsize=(6, 4), layout="constrained")
            self.ax_dev = self.fig_dev.add_subplot()
            self.canvas_dev = FigureCanvasQTAgg(self.fig_dev, export_name="localite_deviation")
            self.dev_stack = EmptyStack(self.canvas_dev)
            tabs.addTab(self.dev_stack, "Deviation per pulse")
        else:
            self.canvas_pos = self.canvas_dev = None
            self.pos_stack = self.dev_stack = None
        self.raw = QTableWidget(0, 5)
        self.raw.setHorizontalHeaderLabels(["Type", "File", "Folder (relative)", "Size (kB)", "Status"])
        self.raw.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        tabs.addTab(self.raw, "Raw marker files")
        self.output = QPlainTextEdit(readOnly=True)
        tabs.addTab(self.output, "Output")
        self.log = QPlainTextEdit(readOnly=True)
        tabs.addTab(self.log, "Log")
        self.tabs = tabs
        split.addWidget(tabs)
        split.setSizes([220, 420])
        outer.addWidget(split, 1)
        self._show_plot_empty()

    @property
    def model(self):
        return self.window.model

    def _idle_status(self, state: str) -> None:
        """While MATLAB is off, say which action starts it (nothing starts it on show)."""
        if state == "off":
            c = style.active()
            self.matlab_status.label.setText(f"<span style='color:{c['ink3']}'>●</span> MATLAB starts "
                                             "when you detect series")

    def _show_plot_empty(self) -> None:
        """Before a detection the plots have nothing to show: say what to do instead."""
        if self.pos_stack is None:
            return
        if not lc.raw_root(self.model):
            args = ("No Localite exports folder set",
                    "Set path.localite_raw to the folder with the sub-*/ses-* raw exports.",
                    "Set Localite folder…", lambda: self.window.goto("path.localite_raw"))
        elif not getattr(self, "sessions", None):
            args = ("No sessions found",
                    "No sub-*/ses-* folders under path.localite_raw. Check the folder.",
                    "Check Localite folder…", lambda: self.window.goto("path.localite_raw"))
        else:
            args = ("No series detected yet",
                    "Pick a subject and session, then detect their recorded series to plot the "
                    "pulse positions and their spread. This starts MATLAB (15–60 s the first time).",
                    "Detect series", self.detect)
        self.pos_stack.show_empty(*args)
        self.dev_stack.show_empty(*args)

    # ── inputs ──────────────────────────────────────────────────────────
    SUMMARY_SHOWN = 2       # settings_summary rows always shown (the folders)

    def _show_summary(self) -> None:
        rows = lc.settings_summary(self.model)
        shown = rows if self._summary_open else rows[:self.SUMMARY_SHOWN]
        lines = "".join(f"<tr><td>{k}</td><td style='padding-left:10px'>{html.escape(str(v))}</td></tr>"
                        for k, v in shown)
        n_more = len(rows) - self.SUMMARY_SHOWN
        toggle = ("<a href='fold'>Fewer settings</a>" if self._summary_open
                  else f"<a href='fold'>{plural(n_more, 'more setting')}</a>") if n_more > 0 else ""
        self.summary.setText(f"<table>{lines}</table>{toggle}")

    def _toggle_summary(self, _href: str = "") -> None:
        self._summary_open = not self._summary_open
        self._show_summary()

    def refresh(self) -> None:
        self._show_summary()
        if self.sub.count() == 0:
            self._fill_subjects()
        self._check_stale()

    def _fill_subjects(self) -> None:
        self.sessions = lc.list_sessions(lc.raw_root(self.model))
        cur = self.sub.currentText()
        self.sub.blockSignals(True)
        self.sub.clear()
        self.sub.addItems(list(self.sessions))
        if cur:
            self.sub.setCurrentText(cur)
        self.sub.blockSignals(False)
        self._sub_changed(self.sub.currentText())
        if not self.sessions:
            self._status("No sub-*/ses-* folders under the raw root (path.localite_raw).", "amber")
        if self.result is None:
            self._show_plot_empty()

    def _sub_changed(self, sub: str) -> None:
        cur = self.ses.currentText()
        self.ses.blockSignals(True)
        self.ses.clear()
        self.ses.addItems(getattr(self, "sessions", {}).get(sub, []))
        want = self.model.get("placement.localite.session")
        label = f"ses-{int(want):02d}" if isinstance(want, (int, float)) and not isinstance(want, bool) else str(want or "")
        for candidate in (cur, label):
            if candidate and self.ses.findText(candidate) >= 0:
                self.ses.setCurrentText(candidate)
                break
        self.ses.blockSignals(False)
        self._inputs_changed()

    def _session(self) -> tuple[str, str]:
        ses = self.ses.currentText().strip()
        if ses.isdigit():
            ses = f"ses-{int(ses):02d}"
        return self.sub.currentText().strip(), ses

    def _inputs_changed(self) -> None:
        self._check_stale()

    def _check_stale(self) -> None:
        if self.snapshot is None:
            return
        if self.snapshot != (self._session(), lc.settings_snapshot(self.model)):
            self.snapshot = None
            self.write_btn.setEnabled(False)
            self._status("Inputs or settings changed: detect the series again before writing.", "ink3")

    def _status(self, text: str, colour: str = "ink2") -> None:
        from . import style
        self.status.setText(f"<span style='color:{style.active()[colour]}'>{text}</span>")

    def _log(self, header: str, text: str = "") -> None:
        self.log.appendPlainText(f"── {header}  ({_dt.datetime.now():%H:%M:%S})")
        if text.strip():
            self.log.appendPlainText(text.rstrip())

    def _saved_config(self) -> str | None:
        """Detect/write read the saved study config, as a run would."""
        if (self.window.dirty or self.model.path is None) and not self.window.save():
            return None
        return str(self.model.path)

    def _failed(self, what: str, err: BridgeError) -> None:
        self.detect_btn.setEnabled(True)
        self._status(f"✗ {err.message}", "red")
        self._log(f"{what} failed", err.report or err.message)

    # ── detect / write ──────────────────────────────────────────────────
    def detect(self) -> None:
        sub, ses = self._session()
        if not sub or not ses:
            self._status("Choose a subject and session.", "red")
            return
        cfg = self._saved_config()
        if cfg is None:
            return
        self.detect_btn.setEnabled(False)
        self.write_btn.setEnabled(False)
        self._status(f"⟳ Detecting series of {sub} {ses}…", "accent")
        snap = ((sub, ses), lc.settings_snapshot(self.model))
        self.window.bridge().call(
            "nexus_neuronav", {"action": "detect", "config": cfg, "sub": sub, "ses": ses},
            on_done=lambda out: self._detected(out, snap), on_error=lambda e: self._failed("Detect", e))

    def _detected(self, out: dict, snap) -> None:
        self.detect_btn.setEnabled(True)
        self.result = out
        self.written = None
        self.use_btn.setEnabled(False)
        self._log(f"Detect {out['sub']} {out['ses']}", out.get("log", ""))
        series = [s for s in as_list(out.get("series"))]
        out["series"] = series
        self._fill_series(series)
        self._fill_raw(as_list(out.get("raw_files")))
        self._plot()
        n_ok = sum(1 for s in series if not s.get("empty"))
        if n_ok == 0:
            self._status("✗ No usable series found: see Log and Raw marker files.", "red")
            return
        self.snapshot = snap
        self.write_btn.setEnabled(True)
        self._status(f"✓ {n_ok} series detected ({out.get('markertype', '')}, "
                     f"{', '.join(as_list(out.get('sources')))}). Names: {out.get('n_from_coil_map', 0)} from "
                     f"the coil table, {out.get('n_from_targets', 0)} from planned targets. Check them, then write.",
                     "green")

    def _fill_series(self, series: list[dict]) -> None:
        t = self.table
        t.blockSignals(True)
        t.setRowCount(len(series))
        for r, s in enumerate(series):
            use = QTableWidgetItem()
            use.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled)
            use.setCheckState(Qt.CheckState.Checked if s.get("tracked") and not s.get("empty")
                              else Qt.CheckState.Unchecked)
            t.setItem(r, 0, use)
            cells = [s["index"], s.get("coil", ""), s.get("name", "") if not s.get("empty") else "",
                     s.get("name_from", "") if not s.get("empty") else "no data", "",
                     s.get("n_pulses", 0), _num(s.get("pos_dev_mean")), _num(s.get("pos_dev_std")),
                     _num(s.get("rot_dev_mean"), "{:.3f}"), _num(s.get("rot_dev_std"), "{:.3f}"), s.get("flag", "")]
            for c, val in enumerate(cells, start=1):
                t.setItem(r, c, _item(val, editable=c in (3, 5) and not s.get("empty")))
        t.blockSignals(False)

    def _fill_raw(self, rows: list[dict]) -> None:
        self.raw.setRowCount(len(rows))
        for r, f in enumerate(rows):
            for c, key in enumerate(("type", "file", "folder", "kb", "status")):
                self.raw.setItem(r, c, _item(f.get(key, "")))

    def _series_edited(self, item: QTableWidgetItem) -> None:
        if item.column() in (0, 5):
            self._plot()

    def _table_values(self) -> tuple[list[bool], list[str], list[str]]:
        t = self.table
        ticked = [t.item(r, 0).checkState() == Qt.CheckState.Checked for r in range(t.rowCount())]
        names = [t.item(r, 3).text() for r in range(t.rowCount())]
        colours = [t.item(r, 5).text() for r in range(t.rowCount())]
        return ticked, names, colours

    def write(self) -> None:
        self._check_stale()
        if self.snapshot is None or self.result is None:
            return
        cfg = self._saved_config()
        if cfg is None:
            return
        sub, ses = self.snapshot[0]
        ticked, names, colours = self._table_values()
        tm = lc.target_map(self.result["series"], ticked, names, colours)
        self.write_btn.setEnabled(False)
        self._status("⟳ Writing outputs…", "accent")
        self.window.bridge().call(
            "nexus_neuronav", {"action": "write", "config": cfg, "sub": sub, "ses": ses, "target_map": tm},
            on_done=lambda out: self._wrote(out, tm), on_error=self._write_failed)

    def _write_failed(self, err: BridgeError) -> None:
        self.write_btn.setEnabled(self.snapshot is not None)
        self._failed("Write", err)

    def _wrote(self, out: dict, tm: list[dict]) -> None:
        self.write_btn.setEnabled(True)
        self.written = dict(out, targets=[e["target_name"] for e in tm])
        self.output.setPlainText(f"{out['xml']}\n{out['json']}\n\n{out.get('json_text', '')}")
        self.tabs.setCurrentWidget(self.output)
        self._log("Write", f"{out['xml']}")
        self.use_btn.setEnabled(True)
        if out.get("kept_existing"):
            self._status(f"Outputs already exist and io.overwrite_files is 'never': kept {Path(out['xml']).parent}",
                         "amber")
        else:
            self._status(f"✓ Wrote {plural(out['n_targets'], 'target')} to {Path(out['xml']).parent}", "green")

    def use_in_study(self) -> None:
        if not self.written or self.snapshot is None:
            return
        changes = lc.use_in_study(self.model, self.snapshot[0][1], self.written["targets"])
        if self.model.get("placement.mode") != "localite":
            self._status("Set; note that placement.mode is not 'localite', so the run does not use them yet.",
                         "amber")
        self.window._changed()
        self.window.statusBar().showMessage(
            "Set " + ", ".join(f"{k} = {v}" for k, v in changes), 8000)

    def _import(self) -> None:
        from .main_window import STUDY
        self.window.show_special(STUDY).import_neuronav()

    # ── plots ───────────────────────────────────────────────────────────
    def _plot(self) -> None:
        if self.canvas_pos is None:
            return
        ax, axd = self.ax_pos, self.ax_dev
        ax.clear()
        axd.clear()
        series = (self.result or {}).get("series") or []
        if not any(not s.get("empty") for s in series):
            self._show_plot_empty()
            if self.result is not None:
                args = ("No usable series in this session",
                        "The detected series have no pulses; see the Log and Raw marker files tabs.")
                self.pos_stack.show_empty(*args)
                self.dev_stack.show_empty(*args)
            return
        self.pos_stack.show_content()
        self.dev_stack.show_content()
        ticked, names, colours = self._table_values() if self.table.rowCount() == len(series) else ([], [], [])
        gum = (self.result or {}).get("markertype") == "GUMMarkers"
        for r, s in enumerate(series):
            if s.get("empty"):
                continue
            c = colours[r].strip() if r < len(colours) else ""
            if not (c.startswith("#") and QColor(c).isValid()):
                c = PALETTE[r % len(PALETTE)]
            pts = [as_floats(p) for p in as_list(s.get("positions"))]
            label = f"{s['index']} · {names[r] if r < len(names) else s.get('name', '')}"
            if pts:
                xs, ys, zs = zip(*pts)
                ax.scatter(xs, ys, zs, s=4, color=c, alpha=0.7, label=label)
            m = as_floats(s.get("mean_pos"))
            d = as_floats(s.get("beam_dir"))
            if len(m) == 3:
                ax.scatter([m[0]], [m[1]], [m[2]], s=80, facecolors="none", edgecolors=c, linewidths=1.5)
                if len(d) == 3 and not gum:
                    ax.quiver(m[0], m[1], m[2], *(30 * x for x in d), color=c, lw=1.5, arrow_length_ratio=0.15)
                    depth = s.get("depth_mm")
                    if isinstance(depth, (int, float)) and depth == depth and s.get("name"):
                        e = [m[i] + depth * d[i] for i in range(3)]
                        ax.plot([m[0], e[0]], [m[1], e[1]], [m[2], e[2]], ":", color=c, lw=1)
            dev = as_floats(s.get("pos_dev"))
            axd.plot(range(1, len(dev) + 1), dev, ".-", color=c, lw=1.1, ms=5, label=label)
        pal = style.active()
        has_target = False
        for cand in as_list((self.result or {}).get("candidates")):
            xyz = as_floats(cand.get("xyz"))
            if len(xyz) == 3:
                has_target = True
                ax.scatter([xyz[0]], [xyz[1]], [xyz[2]], marker="x", color=pal["ink"], s=30)
                ax.text(xyz[0], xyz[1], xyz[2], "  " + str(cand.get("description", "")),
                        fontsize=style.PLOT_PT["small"], color=pal["ink2"])
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
        ax.set_zlabel("z [mm]")
        axd.set_xlabel("Pulse within series")
        axd.set_ylabel("Distance from series mean [mm]")
        axd.grid(True, alpha=0.3)
        if series:
            # series by colour, then what each mark means (in the theme's ink)
            handles, _labels = ax.get_legend_handles_labels()
            ink = pal["ink3"]
            handles += [Line2D([], [], ls="none", marker=".", ms=6, color=ink, label="pulse"),
                        Line2D([], [], ls="none", marker="o", ms=9, mfc="none", mec=ink, mew=1.5,
                               label="series mean")]
            if not gum:
                handles.append(Line2D([], [], color=ink, lw=1.5, label="beam (30 mm)"))
            if has_target:
                handles.append(Line2D([], [], ls="none", marker="x", ms=6, color=pal["ink"],
                                      label="planned target"))
            ax.legend(handles=handles, fontsize=style.PLOT_PT["small"], loc="upper left")
            if axd.get_legend_handles_labels()[0]:
                axd.legend(fontsize=style.PLOT_PT["small"])
        self.canvas_pos.draw_idle()
        self.canvas_dev.draw_idle()

    # ── batch / MNI ─────────────────────────────────────────────────────
    def open_batch(self) -> None:
        cfg = self._saved_config()
        if cfg is None:
            return
        BatchDialog(self.window, self).exec()

    def open_mni(self) -> None:
        cfg = self._saved_config()
        if cfg is None:
            return
        MniDialog(self.window, self).exec()


class BatchDialog(QDialog):
    """Scan all sessions (neuronav_scan_sessions), edit the label map (neuronav.target_labels),
    pick sessions and run neuronav_ingest_batch as a script on the Run tab."""

    def __init__(self, window, parent=None):
        super().__init__(parent)
        self.window = window
        self.setWindowTitle("Ingest Localite sessions")
        self.resize(1050, 680)
        v = QVBoxLayout(self)
        intro = QLabel("Uses the study config's settings. Tick sessions, check the label map, run. The "
                       "batch runs on the Run tab (<tt>neuronav_ingest_batch</tt>); afterwards, "
                       "<b>Import into study…</b> turns the written targets into study runs.")
        intro.setWordWrap(True)
        v.addWidget(intro)
        top = QHBoxLayout()
        scan = QPushButton("↻ Scan sessions")
        scan.setProperty("primary", True)
        scan.clicked.connect(self.scan)
        top.addWidget(scan)
        for text, mode in (("All", "all"), ("None", "none"), ("Not yet ingested", "new")):
            b = QPushButton(text)
            b.clicked.connect(lambda _=False, m=mode: self.select(m))
            top.addWidget(b)
        self.count = QLabel()
        top.addWidget(self.count, 1)
        stem = QPushButton("Group labels by stem")
        stem.setToolTip("Unmapped labels: target = label without a trailing _<number>")
        stem.clicked.connect(lambda: self.label_targets("stem"))
        top.addWidget(stem)
        asis = QPushButton("Use labels as found")
        asis.clicked.connect(lambda: self.label_targets("label"))
        top.addWidget(asis)
        v.addLayout(top)
        split = QSplitter()
        self.sessions = QTableWidget(0, 8)
        self.sessions.setHorizontalHeaderLabels(["Include", "Subject", "Session", "Trigger", "GUM",
                                                 "EntryTarget", "Planned labels", "Ingested"])
        self.sessions.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        split.addWidget(self.sessions)
        self.labels = QTableWidget(0, 5)
        self.labels.setHorizontalHeaderLabels(["Label (as recorded)", "Sessions", "Subjects", "Stem", "Target name"])
        self.labels.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        split.addWidget(self.labels)
        v.addWidget(split, 1)
        bottom = QHBoxLayout()
        bottom.addWidget(QLabel("Write"))
        self.write_mode = QComboBox()
        self.write_mode.addItems(["all", "named"])
        self.write_mode.setToolTip("all: every tracked series · named: only series named by the coil "
                                   "table or a planned target (neuronav.batch_write)")
        self.write_mode.setCurrentText(str(window.model.get("neuronav.batch_write") or "all"))
        bottom.addWidget(self.write_mode)
        bottom.addStretch(1)
        self.run_btn = QPushButton("▶ Save label map and run batch")
        self.run_btn.setProperty("primary", True)
        self.run_btn.clicked.connect(self.run)
        bottom.addWidget(self.run_btn)
        close = QPushButton("Close")
        close.clicked.connect(self.reject)
        bottom.addWidget(close)
        v.addLayout(bottom)
        self.scan()

    def scan(self) -> None:
        self.count.setText("⟳ Scanning…")
        self.window.bridge().call(
            "nexus_neuronav", {"action": "scan", "config": str(self.window.model.path)},
            on_done=self._scanned, on_error=lambda e: self.count.setText(f"✗ {e.message}"))

    def _scanned(self, out: dict) -> None:
        sess, labs = as_list(out.get("sessions")), as_list(out.get("labels"))
        self.sessions.setRowCount(len(sess))
        for r, s in enumerate(sess):
            usable = (s.get("n_trigger") or 0) + (s.get("n_gum") or 0) > 0
            inc = QTableWidgetItem()
            inc.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled)
            inc.setCheckState(Qt.CheckState.Checked if usable else Qt.CheckState.Unchecked)
            self.sessions.setItem(r, 0, inc)
            for c, key in enumerate(("sub", "ses", "n_trigger", "n_gum", "n_entrytarget", "labels"), start=1):
                self.sessions.setItem(r, c, _item(s.get(key, "")))
            self.sessions.setItem(r, 7, _item("✓" if s.get("ingested") else ""))
        self.labels.setRowCount(len(labs))
        for r, lab in enumerate(labs):
            for c, key in enumerate(("label", "n_sessions", "n_subjects", "stem")):
                self.labels.setItem(r, c, _item(lab.get(key, "")))
            self.labels.setItem(r, 4, _item(lab.get("target", lab.get("label", "")), editable=True))
        n_sub = len({s.get("sub") for s in sess})
        self.count.setText(f"{plural(len(sess), 'session')} in {plural(n_sub, 'subject')}; "
            f"{plural(len(labs), 'distinct planned-target label')}.")

    def select(self, mode: str) -> None:
        for r in range(self.sessions.rowCount()):
            usable = int(self.sessions.item(r, 3).text() or 0) + int(self.sessions.item(r, 4).text() or 0) > 0
            ingested = bool(self.sessions.item(r, 7).text())
            on = {"all": usable, "none": False, "new": usable and not ingested}[mode]
            self.sessions.item(r, 0).setCheckState(Qt.CheckState.Checked if on else Qt.CheckState.Unchecked)

    def label_targets(self, mode: str) -> None:
        for r in range(self.labels.rowCount()):
            label, target = self.labels.item(r, 0).text(), self.labels.item(r, 4).text()
            if mode == "label":
                self.labels.item(r, 4).setText(label)
            elif target == label:                      # keep curated entries
                self.labels.item(r, 4).setText(self.labels.item(r, 3).text())

    def label_map(self) -> list[dict]:
        """neuronav.target_labels: only labels whose target differs from the label (as neuronav_gui)."""
        groups: dict[str, list[str]] = {}
        for r in range(self.labels.rowCount()):
            label = self.labels.item(r, 0).text()
            target = self.labels.item(r, 4).text().strip() or label
            if target != label:
                groups.setdefault(target, []).append(label)
        return [{"target": t, "labels": labs} for t, labs in groups.items()]

    def selection(self) -> list[tuple[str, str]]:
        return [(self.sessions.item(r, 1).text(), self.sessions.item(r, 2).text())
                for r in range(self.sessions.rowCount())
                if self.sessions.item(r, 0).checkState() == Qt.CheckState.Checked]

    def run(self) -> None:
        from .. import launcher
        from .main_window import RUN
        sel = self.selection()
        if not sel:
            self.count.setText("Tick at least one session.")
            return
        m = self.window.model
        if lc.post_root(m) == "":
            self.count.setText("✗ Set path.localite_post first (Settings tab).")
            return
        labels = self.label_map()
        if labels or m.get("neuronav.target_labels"):
            m.set("neuronav.target_labels", labels)
        m.set("neuronav.batch_write", self.write_mode.currentText())
        self.window._changed()
        if not self.window.save():
            return
        script = launcher.write_ingest_script(m.path, self.window.prestus_root, sel)
        page = self.window.show_special(RUN)
        if page.run_script(f"Localite ingestion ({len(sel)} sessions)", script):
            self.accept()


class MniDialog(QDialog):
    """Targets of the batch summary in MNI space (neuronav_summary_to_mni), filterable."""

    def __init__(self, window, parent=None):
        super().__init__(parent)
        self.window = window
        self.rows: list[dict] = []
        self.head: tuple[list, list] | None = None
        self.setWindowTitle("Localite targets in MNI space")
        self.resize(1150, 760)
        m = window.model
        g = QHBoxLayout(self)
        left = QWidget()
        left.setMaximumWidth(380)
        form = QFormLayout(left)
        post = lc.post_root(m)
        self.csv = QLineEdit(str(Path(post) / "neuronav_ingest_summary.csv") if post else "")
        form.addRow("Summary CSV", self.csv)
        self.seg = QLineEdit(str(m.get("path.seg") or ""))
        self.seg.setToolTip("path.seg: m2m_<sub>/final_tissues.nii.gz and toMNI/Conform2MNI_nonl.nii.gz")
        form.addRow("SimNIBS root", self.seg)
        self.tpl = QLineEdit()
        self.tpl.setPlaceholderText("SimNIBS MNI152 (found automatically)")
        form.addRow("Head template", self.tpl)
        go = QPushButton("▶ Load && map to MNI")
        go.setProperty("primary", True)
        go.clicked.connect(self.load)
        form.addRow(go)
        self.status = QLabel("Map the summary to see its targets.")
        self.status.setWordWrap(True)
        form.addRow(self.status)
        self.show_trans = QCheckBox("Transducer", checked=True)
        self.show_focus = QCheckBox("Focus", checked=True)
        self.show_beam = QCheckBox("Beam axis", checked=True)
        self.show_head = QCheckBox("Head outline", checked=True)
        box = QGridLayout()
        for i, w in enumerate((self.show_trans, self.show_focus, self.show_beam, self.show_head)):
            w.toggled.connect(self.plot)
            box.addWidget(w, i // 2, i % 2)
        form.addRow(box)
        self.color_by = QComboBox()
        self.color_by.addItems(["target", "sub", "ses"])
        self.color_by.currentIndexChanged.connect(self.plot)
        form.addRow("Colour by", self.color_by)
        self.filters = {}
        for key, title in (("sub", "Subjects"), ("ses", "Sessions"), ("target", "Targets")):
            t = QTableWidget(0, 2)
            t.setHorizontalHeaderLabels(["Show", title])
            t.horizontalHeader().setStretchLastSection(True)
            t.itemChanged.connect(lambda _i: self.plot())
            t.setMaximumHeight(130)
            self.filters[key] = t
            form.addRow(t)
        g.addWidget(left)
        if FigureCanvasQTAgg is not None:
            self.fig = Figure(figsize=(7, 6), layout="constrained")
            self.ax = self.fig.add_subplot(projection="3d")
            self.canvas = FigureCanvasQTAgg(self.fig, export_name="localite_mni_targets")
            self.plot_stack = EmptyStack(self.canvas)
            self.plot_stack.show_empty("No targets mapped yet",
                                       "Map the batch summary to MNI space to see every written "
                                       "target. This starts MATLAB and loads each subject's warp.",
                                       "Load and map to MNI", self.load)
            g.addWidget(self.plot_stack, 1)
        else:
            self.canvas = None
            g.addWidget(QLabel("Install matplotlib to see the MNI view."), 1)

    def load(self) -> None:
        self.status.setText("⟳ Mapping to MNI… (the first subject loads its warp)")
        args = {"action": "mni", "config": str(self.window.model.path), "csv": self.csv.text().strip(),
                "seg": self.seg.text().strip(), "head": True}
        if self.tpl.text().strip():
            args["template"] = self.tpl.text().strip()
        self.window.bridge().call("nexus_neuronav", args, on_done=self._loaded,
                                  on_error=lambda e: self.status.setText(f"✗ {e.message}"),
                                  on_progress=lambda msg: self.status.setText(f"⟳ {msg}"))

    def _loaded(self, out: dict) -> None:
        self.rows = as_list(out.get("rows"))
        faces, verts = as_list(out.get("faces")), as_list(out.get("vertices"))
        self.head = (faces, verts) if faces and verts else None
        if out.get("template"):
            self.tpl.setText(out["template"])
        for key, t in self.filters.items():
            names = sorted({str(r.get(key, "")) for r in self.rows})
            t.blockSignals(True)
            t.setRowCount(len(names))
            for i, n in enumerate(names):
                it = QTableWidgetItem()
                it.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled)
                it.setCheckState(Qt.CheckState.Checked)
                t.setItem(i, 0, it)
                t.setItem(i, 1, _item(n))
            t.blockSignals(False)
        ok = [r for r in self.rows if not r.get("mni_status")]
        bad = sorted({f"{r.get('sub')} ({r.get('mni_status')})" for r in self.rows if r.get("mni_status")})
        self.status.setText(f"✓ {plural(len(ok), 'target')} of {plural(len({r.get('sub') for r in ok}), 'subject')} mapped."
                            + (f" Not mapped: {'; '.join(bad)}" if bad else ""))
        if self.canvas is not None:
            if ok:
                self.plot_stack.show_content()
            else:
                self.plot_stack.show_empty("No target could be mapped",
                                           "See the status on the left for each subject's reason.")
        self.plot()

    def _shown(self, key: str) -> set[str]:
        t = self.filters[key]
        return {t.item(i, 1).text() for i in range(t.rowCount())
                if t.item(i, 0).checkState() == Qt.CheckState.Checked}

    def plot(self) -> None:
        if self.canvas is None:
            return
        ax = self.ax
        ax.clear()
        if self.show_head.isChecked() and self.head is not None:
            faces = [[int(i) - 1 for i in as_floats(f)] for f in self.head[0]]
            verts = [as_floats(p) for p in self.head[1]]
            xs, ys, zs = zip(*verts)
            ax.plot_trisurf(xs, ys, zs, triangles=faces, color="#dcc7b8", alpha=0.12, linewidth=0)
        shown = {k: self._shown(k) for k in self.filters}
        use = [r for r in self.rows if not r.get("mni_status")
               and all(str(r.get(k, "")) in shown[k] for k in shown)]
        key = self.color_by.currentText()
        groups = sorted({str(r.get(key, "")) for r in use})
        for i, grp in enumerate(groups):
            c = PALETTE[i % len(PALETTE)]
            rs = [r for r in use if str(r.get(key, "")) == grp]
            tr = [(r["trans_mni_x"], r["trans_mni_y"], r["trans_mni_z"]) for r in rs]
            fo = [(r["focus_mni_x"], r["focus_mni_y"], r["focus_mni_z"]) for r in rs]
            if self.show_beam.isChecked():
                for a, b in zip(tr, fo):
                    if None not in a and None not in b:
                        ax.plot(*zip(a, b), "-", color=c, lw=0.8, alpha=0.6)
            label = f"{grp} ({len(rs)})"
            if self.show_trans.isChecked():
                pts = [p for p in tr if None not in p]
                if pts:
                    ax.scatter(*zip(*pts), marker="o", s=30, color=c, label=label)
                    label = None
            if self.show_focus.isChecked():
                pts = [p for p in fo if None not in p]
                if pts:
                    ax.scatter(*zip(*pts), marker="D", s=26, color=c, label=label)
        ax.set_xlabel("x MNI, +R [mm]")
        ax.set_ylabel("y MNI, +A [mm]")
        ax.set_zlabel("z MNI, +S [mm]")
        ax.set_title(f"Targets in MNI space ({len(use)} shown)", fontsize=style.PLOT_PT["title"])
        if groups:
            # groups by colour, then what each mark means (in the theme's ink)
            handles, _labels = ax.get_legend_handles_labels()
            ink = style.active()["ink3"]
            if self.show_trans.isChecked():
                handles.append(Line2D([], [], ls="none", marker="o", ms=6, color=ink, label="transducer"))
            if self.show_focus.isChecked():
                handles.append(Line2D([], [], ls="none", marker="D", ms=5, color=ink, label="focus"))
            if self.show_beam.isChecked():
                handles.append(Line2D([], [], color=ink, lw=0.8, label="beam axis"))
            ax.legend(handles=handles, fontsize=style.PLOT_PT["small"], loc="upper left")
        self.canvas.draw_idle()
