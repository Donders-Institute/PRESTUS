"""Cluster jobs: the job table and the selected job's logs (output, errors, report).

The jobs themselves (list, scheduler, SSH host, 30 s polling, cancel) live on the
RunPage; a JobsView only shows them. One view is the Run section's "Jobs" side tab;
"Pop out" opens another in its own window (JobsWindow), which can be maximised or
moved to a second screen. Each view follows its selected job's log every 10 s.
"""

from __future__ import annotations

import html

from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QBrush, QColor, QFont
from PyQt6.QtWidgets import (QAbstractItemView, QCheckBox, QHBoxLayout, QLabel, QPlainTextEdit,
                             QPushButton, QSplitter, QTabWidget, QToolButton, QTreeWidget,
                             QTreeWidgetItem, QVBoxLayout, QWidget)

from .. import jobs as hpc_jobs

JOB_COLUMNS = ["Job ID", "Name", "State", "Elapsed / limit", "Left", "Where", "Peak memory", "Exit"]
LOG_TABS = (("output", "Output"), ("error", "Errors"), ("report", "Report"))
FONT_PT = 11                      # default log font size; A-/A+ change it (session: job_log_font_pt)


def job_row(j: hpc_jobs.Job) -> list[str]:
    timing = f"{j.elapsed} / {j.time_limit}" if j.time_limit else j.elapsed
    where = f"{j.partition}: {j.where}" if j.partition and j.where else (j.where or j.partition)
    return [j.job_id, j.name, j.state, timing, j.time_left, where, j.max_rss, j.exit_code]


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
        head.addWidget(self.b_refresh)
        head.addWidget(self.b_cancel)
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
        self.tree.setMinimumHeight(70)
        self.tree.itemSelectionChanged.connect(self._selected)
        if not popout:
            self.tree.itemDoubleClicked.connect(lambda item, _c: page.pop_out_jobs(item.text(0)))
        self.split.addWidget(self.tree)

        logs = QWidget()
        lv = QVBoxLayout(logs)
        lv.setContentsMargins(0, 4, 0, 0)
        self.path_label = QLabel("Select a job to see its output, errors and report.")
        self.path_label.setWordWrap(True)
        self.path_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        lv.addWidget(self.path_label)
        self.tabs = QTabWidget()
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
            item.setToolTip(5, j.work_dir and f"Logs in {j.work_dir}")
            if j.time_warning():
                for col in (3, 4):
                    item.setForeground(col, QBrush(QColor("#c9302c")))
                item.setToolTip(4, "Little time left: the job may be stopped (TIMEOUT). "
                                   "Raise hpc.timelimit for the next run.")
            if j.finished and not j.state.startswith("COMPLETED"):
                item.setForeground(2, QBrush(QColor("#c9302c")))
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
        self.refresh_log()

    def refresh_log(self) -> None:
        """Show the selected job's current tab; follow it while the job is open."""
        job = self.selected_job()
        if job is None:
            self.timer.stop()
            self.path_label.setText("Select a job to see its output, errors and report.")
            return
        if self.page.job_kind != "slurm":
            self.path_label.setText("Job logs are shown for SLURM jobs only.")
            return
        if not job.work_dir:              # submitted, not yet seen by squeue/sacct
            self.page.refresh_jobs()
        key = LOG_TABS[self.tabs.currentIndex()][0]
        view = self.views[key]
        path, text = hpc_jobs.read_log(job, key, self.page.job_remote)
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
        font = QFont("Menlo")
        font.setStyleHint(QFont.StyleHint.Monospace)
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
