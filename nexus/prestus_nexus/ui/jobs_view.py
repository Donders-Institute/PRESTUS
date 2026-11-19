"""Cluster jobs: the job table and the selected job's logs (output, errors, report).

The jobs themselves (list, scheduler, SSH host, 30 s polling, cancel) live on the
RunPage; a JobsView only shows them. One view is the Run section's "Jobs" side tab;
"Pop out" opens another in its own window (JobsWindow), which can be maximised or
moved to a second screen. Each view follows its selected job's log every 10 s; logs are
read on a worker thread (over SSH that takes seconds), one read at a time per view.
"""

from __future__ import annotations

import copy
import html
import time

from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QBrush, QColor
from PyQt6.QtWidgets import (QAbstractItemView, QCheckBox, QHBoxLayout, QLabel, QPlainTextEdit,
                             QPushButton, QSplitter, QTabWidget, QToolButton, QTreeWidget,
                             QTreeWidgetItem, QVBoxLayout, QWidget)

from .. import jobs as hpc_jobs
from . import style
from .progress_panel import ProgressPanel
from .tabbar import EmphasisTabBar

JOB_COLUMNS = ["Job ID", "Name", "State", "Progress", "Elapsed / limit", "Left", "Where", "Peak memory", "Exit"]
PROGRESS_COL = 3
LOG_TABS = (("output", "Output"), ("error", "Errors"), ("report", "Report"))
FONT_PT = 11                      # default log font size; A-/A+ change it (session: job_log_font_pt)


def progress_cell(j: hpc_jobs.Job, now: float | None = None) -> str:
    """'acoustic 45 %, 3 min left' style summary of the job's progress ('' if none)."""
    p = hpc_jobs.job_progress(j)
    if p is None or not j.running:
        return ""
    now = time.time() if now is None else now
    est = p.estimate(now)
    whole = p.overall(j.expected, now)
    warn = " · exceeds time limit" if hpc_jobs.timeout_risk(j, whole[1] if whole else est.eta_s) else ""
    stage = p.solver.stage if p.solver is not None else p.stage
    if not stage:
        return est.outer.split(" · ")[0] if est.outer else ""     # e.g. a drive check: 'Depth 2 of 5'
    if est.outer:
        stage = f"{est.outer.split(' · ')[0]}: {stage}"
    text = stage.replace("_", " ")
    if est.fraction is not None:
        text += f" {'~' if est.estimated else ''}{100 * est.fraction:.0f} %"
    if est.eta_s is not None:
        from ..progress import format_duration
        text += f", {format_duration(est.eta_s)} left"
    return text + warn


def timeout_text(j: hpc_jobs.Job, eta_s: float | None) -> str:
    """The card's warning for a job that will probably not finish in its time limit."""
    if not hpc_jobs.timeout_risk(j, eta_s):
        return ""
    from ..progress import format_duration
    return (f"About {format_duration(eta_s)} of work left, but only {j.time_left} before the time "
            "limit: the job will probably end in TIMEOUT. Raise hpc.timelimit for the next run.")


def job_row(j: hpc_jobs.Job) -> list[str]:
    timing = f"{j.elapsed} / {j.time_limit}" if j.time_limit else j.elapsed
    where = f"{j.partition}: {j.where}" if j.partition and j.where else (j.where or j.partition)
    name = f"↳ {j.name}" if j.parent else j.name        # submitted by another job (its parent)
    return [j.job_id, name, j.state, progress_cell(j), timing, j.time_left, where, j.max_rss, j.exit_code]


class JobsView(QWidget):
    def __init__(self, page, popout: bool = False, parent: QWidget | None = None):
        super().__init__(parent)
        self.page = page
        self.popout = popout
        v = QVBoxLayout(self)
        m = 8 if popout else 0
        v.setContentsMargins(m, m, m, m)

        head = QHBoxLayout()
        self.label = QLabel("Cluster jobs")
        head.addWidget(self.label)
        head.addStretch(1)
        self.b_refresh = QPushButton("Refresh")
        self.b_refresh.clicked.connect(page.refresh_jobs)
        self.b_cancel = QPushButton("Cancel selected")
        self.b_cancel.clicked.connect(self.cancel)
        self.b_clear = QPushButton("Clear finished")
        self.b_clear.setToolTip("Remove the jobs that have ended from this list (the jobs and their "
                                "outputs are not touched)")
        self.b_clear.clicked.connect(page.clear_finished_jobs)
        head.addWidget(self.b_refresh)
        head.addWidget(self.b_cancel)
        head.addWidget(self.b_clear)
        if not popout:
            self.b_popout = QPushButton("Pop out ↗")
            self.b_popout.setToolTip("Open the jobs and their logs in a separate window "
                                     "(also: double-click a job)")
            self.b_popout.clicked.connect(lambda: page.pop_out_jobs(self.selected_id()))
            head.addWidget(self.b_popout)
        v.addLayout(head)

        self.split = QSplitter(Qt.Orientation.Vertical)
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(JOB_COLUMNS)
        self.tree.setRootIsDecorated(False)
        self.tree.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.tree.setMinimumHeight(110)       # about four jobs
        self.tree.itemSelectionChanged.connect(self._selected)
        if not popout:
            self.tree.itemDoubleClicked.connect(lambda item, _c: page.pop_out_jobs(item.text(0)))
        self.split.addWidget(self.tree)

        logs = QWidget()
        lv = QVBoxLayout(logs)
        lv.setContentsMargins(0, 4, 0, 0)
        self.path_label = QLabel("Select a job to see its output, errors, and report.")
        self.path_label.setWordWrap(True)
        self.path_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        lv.addWidget(self.path_label)
        # The selected job's progress (and live frame) from its output log (progress.py)
        self.progress_panel = ProgressPanel()
        self.progress_panel.live_requested.connect(lambda: self._show_selected_progress(select_live=True))
        lv.addWidget(self.progress_panel)
        self.tabs = QTabWidget()
        self.tabs.setTabBar(EmphasisTabBar(self.tabs))     # the Nexus's tab dividers
        self.views: dict[str, QPlainTextEdit] = {}
        for key, title in LOG_TABS:
            view = QPlainTextEdit(readOnly=True)
            view.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
            self.views[key] = view
            self.tabs.addTab(view, title)
        self.tabs.currentChanged.connect(lambda _i: self.refresh_log())
        lv.addWidget(self.tabs, 1)
        tools = QHBoxLayout()
        self.follow = QCheckBox("Follow")
        self.follow.setChecked(True)
        self.follow.setToolTip("Keep the end of the log in view as it grows")
        self.wrap = QCheckBox("Wrap lines")
        self.wrap.toggled.connect(self._set_wrap)
        tools.addWidget(self.follow)
        tools.addWidget(self.wrap)
        tools.addStretch(1)
        for text, step in (("A−", -1), ("A+", 1)):
            b = QToolButton(text=text)
            b.setToolTip("Smaller text" if step < 0 else "Larger text")
            b.clicked.connect(lambda _c=False, s=step: self._zoom(s))
            tools.addWidget(b)
        lv.addLayout(tools)
        self.split.addWidget(logs)
        self.split.setStretchFactor(0, 1)
        self.split.setStretchFactor(1, 4)
        v.addWidget(self.split, 1)
        self._apply_font()

        self.timer = QTimer(self)
        self.timer.setInterval(10_000)
        self.timer.timeout.connect(self.refresh_log)
        self._log_read = None                 # (job id, tab key) of the read in flight
        self._log_again = False               # selection or tab changed during that read
        from .run_page import Relay
        self._relay = Relay(self)
        self._relay.done.connect(self._log_done, Qt.ConnectionType.QueuedConnection)
        # Estimates move between the 30 s polls: refresh the progress once a second
        self.tick = QTimer(self)
        self.tick.setInterval(1000)
        self.tick.timeout.connect(self._tick)
        page.jobs_changed.connect(self.show_jobs)
        self.show_jobs()

    # side tabs call refresh() when shown
    def refresh(self) -> None:
        self.show_jobs()

    # ── table ───────────────────────────────────────────────────────────
    def show_jobs(self) -> None:
        jobs = self.page.jobs
        selected = {it.text(0) for it in self.tree.selectedItems()}
        self.tree.blockSignals(True)
        self.tree.clear()
        for j in jobs:
            item = QTreeWidgetItem(job_row(j))
            item.setToolTip(6, j.work_dir and f"Logs in {j.work_dir}")
            if j.time_warning():
                for col in (4, 5):
                    item.setForeground(col, QBrush(QColor(style.active()["red"])))
                item.setToolTip(5, "Little time left: the job may be stopped (TIMEOUT). "
                                   "Raise hpc.timelimit for the next run.")
            if j.finished and not j.state.startswith("COMPLETED"):
                item.setForeground(2, QBrush(QColor(style.active()["red"])))
            if "exceeds time limit" in item.text(PROGRESS_COL):
                item.setForeground(PROGRESS_COL, QBrush(QColor(style.active()["red"])))
            if j.parent:
                item.setToolTip(1, f"Submitted by job {j.parent}")
            self.tree.addTopLevelItem(item)
            item.setSelected(j.job_id in selected)
        for col in range(len(JOB_COLUMNS)):
            self.tree.resizeColumnToContents(col)
        self.tree.blockSignals(False)
        n_open = sum(1 for j in jobs if not j.finished)
        self.label.setText(f"Cluster jobs ({len(jobs)}, {n_open} open)" if jobs
                           else "Cluster jobs (none submitted from this window)")
        has = bool(jobs) and self.page.job_kind is not None
        self.b_refresh.setEnabled(has)
        self.b_cancel.setEnabled(has)
        self.b_clear.setEnabled(any(j.finished for j in jobs))
        self._show_selected_progress()
        if any(j.running and j.progress_log for j in jobs):
            self.tick.start()
        else:
            self.tick.stop()

    # ── progress ────────────────────────────────────────────────────────
    def _tick(self) -> None:
        jobs = {j.job_id: j for j in self.page.jobs}
        for i in range(self.tree.topLevelItemCount()):
            item = self.tree.topLevelItem(i)
            j = jobs.get(item.text(0))
            if j is not None and j.running:
                item.setText(PROGRESS_COL, progress_cell(j))
        self._show_selected_progress()

    def _to_local(self):
        remote = self.page.job_remote
        return remote.to_local_or_none if remote is not None else None

    def _show_selected_progress(self, select_live: bool = False) -> None:
        job = self.selected_job()
        p = hpc_jobs.job_progress(job) if job is not None else None
        self.progress_panel.show_progress(p, time.time(), self._to_local(),
                                          job.expected if job is not None else None)
        if job is not None and p is not None:
            whole = self.progress_panel.overall
            eta = whole[1] if whole is not None else p.estimate(time.time()).eta_s
            self.progress_panel.set_warning(timeout_text(job, eta))
        if job is not None and (self.progress_panel.frames or self.progress_panel.missing):
            self.page.show_job_live(job, p, self._to_local(), select=select_live)
        else:                            # keep following the job the Live view tab shows
            source = self.page.live_page.source
            shown = next((j for j in self.page.jobs if j.job_id == source), None)
            if shown is not None:
                self.page.show_job_live(shown, hpc_jobs.job_progress(shown), self._to_local())

    def selected_id(self) -> str | None:
        items = self.tree.selectedItems()
        return items[-1].text(0) if items else None

    def select(self, job_id: str | None) -> None:
        if job_id is None:
            return
        for i in range(self.tree.topLevelItemCount()):
            item = self.tree.topLevelItem(i)
            item.setSelected(item.text(0) == job_id)
            if item.text(0) == job_id:
                self.tree.scrollToItem(item)

    def selected_job(self) -> hpc_jobs.Job | None:
        job_id = self.selected_id()
        return next((j for j in self.page.jobs if j.job_id == job_id), None)

    def cancel(self) -> None:
        self.page.cancel_jobs([it.text(0) for it in self.tree.selectedItems()])

    # ── logs ────────────────────────────────────────────────────────────
    def _selected(self) -> None:
        for view in self.views.values():
            view.clear()
        self.progress_panel.clear()
        self._show_selected_progress()
        self.refresh_log()

    def refresh_log(self) -> None:
        """Show the selected job's current tab; follow it while the job is open."""
        job = self.selected_job()
        if job is None:
            self.timer.stop()
            self.path_label.setText("Select a job to see its output, errors, and report.")
            return
        if self.page.job_kind != "slurm":
            self.path_label.setText("Job logs are shown for SLURM jobs only.")
            return
        if not job.work_dir:              # submitted, not yet seen by squeue/sacct
            self.page.refresh_jobs()
        key = LOG_TABS[self.tabs.currentIndex()][0]
        if self._log_read is not None:    # one read at a time; read again when it is in
            self._log_again = self._log_read != (job.job_id, key)
            return
        self._log_read = (job.job_id, key)
        from .run_page import submit_background
        self._relay.deliver((job.job_id, key),
                            submit_background(hpc_jobs.read_log, copy.copy(job), key,
                                              self.page.job_remote))

    def log_busy(self) -> bool:
        """A log read is in flight."""
        return self._log_read is not None

    def _log_done(self, tag, fut) -> None:
        self._log_read = None
        again, self._log_again = self._log_again, False
        job = self.selected_job()
        key = LOG_TABS[self.tabs.currentIndex()][0]
        if again or job is None or (job.job_id, key) != tag:   # superseded: show what is selected now
            self.refresh_log()
            return
        try:
            path, text = fut.result()
        except Exception as exc:
            path, text = "", ""
            self.path_label.setText(f"Job {job.job_id}: cannot read the {key} log ({exc}).")
            self.timer.start()
            return
        self._show_log(job, key, path, text)

    def _show_log(self, job: hpc_jobs.Job, key: str, path: str, text: str) -> None:
        view = self.views[key]
        if text != view.toPlainText():
            bar = view.verticalScrollBar()
            pos = bar.value()
            view.setPlainText(text)
            bar.setValue(bar.maximum() if self.follow.isChecked() else pos)
        if path:
            self.path_label.setText(f"Job {job.job_id} — <tt>{html.escape(path)}</tt>")
        elif key == "report":
            self.path_label.setText(f"Job {job.job_id}: no report yet (written by the report "
                                    "job after the job ends).")
        else:
            self.path_label.setText(f"Job {job.job_id}: no {key} log yet"
                                    + ("" if job.work_dir else " (folder not known yet)") + ".")
        if job.finished and (path or key != "report"):    # the report comes after the job
            self.timer.stop()
        else:
            self.timer.start()

    def _set_wrap(self, on: bool) -> None:
        mode = QPlainTextEdit.LineWrapMode.WidgetWidth if on else QPlainTextEdit.LineWrapMode.NoWrap
        for view in self.views.values():
            view.setLineWrapMode(mode)

    def _zoom(self, step: int) -> None:
        session = self.page.window.session
        session.set("job_log_font_pt", max(7, min(24, int(session.get("job_log_font_pt", FONT_PT)) + step)))
        self._apply_font()

    def _apply_font(self) -> None:
        font = style.mono_font()
        font.setPointSize(int(self.page.window.session.get("job_log_font_pt", FONT_PT)))
        for view in self.views.values():
            view.setFont(font)


class JobsWindow(QWidget):
    """The jobs and their logs in a window of their own (RunPage.pop_out_jobs)."""

    def __init__(self, page):
        super().__init__(page, Qt.WindowType.Window)
        self.setWindowTitle("Cluster jobs — PRESTUS Nexus")
        self.view = JobsView(page, popout=True)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.view)
        self.resize(1100, 760)
