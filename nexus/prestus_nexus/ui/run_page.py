"""Run page: write the launcher, validate in MATLAB, run locally or submit to the cluster
with a live log, and follow submitted jobs."""

from __future__ import annotations

import html
from pathlib import Path

from PyQt6.QtCore import QProcess, QTimer, pyqtSignal
from PyQt6.QtGui import QFont, QTextCursor
from PyQt6.QtWidgets import (QAbstractItemView, QCheckBox, QFileDialog, QHBoxLayout, QLabel, QLineEdit,
                             QMessageBox, QPlainTextEdit, QPushButton, QTextEdit, QToolButton,
                             QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWidget)

from .. import jobs as hpc_jobs
from .. import launcher
from ..model import ConfigModel
from ..predict_outputs import predict, save_summary


class RunPage(QWidget):
    finished = pyqtSignal(str, bool)          # task, success

    def __init__(self, window, parent: QWidget | None = None):
        super().__init__(parent)
        self.window = window
        self.proc: QProcess | None = None
        self.task = ""
        self._ok = False

        v = QVBoxLayout(self)
        title = QLabel("Run in MATLAB")
        title.setObjectName("CardTitle")
        v.addWidget(title)
        intro = QLabel("Runs the saved study config in MATLAB on this computer, or submits it to the "
                       "cluster (on a login node with sbatch/qsub) and follows the jobs. "
                       "Validate first: it only loads the config with <tt>load_parameters</tt>.")
        intro.setWordWrap(True)
        v.addWidget(intro)

        self.deploy_toggle = QCheckBox("Show deployment & save details")
        self.deploy_toggle.setToolTip(
            "What prestus_pipeline_start would submit (single job, or one job per transducer/"
            "target/uncertainty variant) and what io.save_* keeps besides the CSV/report.")
        self.deploy_toggle.setChecked(bool(window.session.get("run_show_deployment", False)))
        self.deploy_toggle.toggled.connect(self._toggle_deployment)
        v.addWidget(self.deploy_toggle)
        self.deploy_text = QPlainTextEdit(readOnly=True)
        self.deploy_text.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        deploy_mono = QFont("Menlo")
        deploy_mono.setStyleHint(QFont.StyleHint.Monospace)
        self.deploy_text.setFont(deploy_mono)
        self.deploy_text.setFixedHeight(180)
        self.deploy_text.setVisible(self.deploy_toggle.isChecked())
        v.addWidget(self.deploy_text)

        row = QHBoxLayout()
        row.addWidget(QLabel("MATLAB"))
        self.matlab = QLineEdit()
        self.matlab.setPlaceholderText("path to the matlab binary")
        browse = QToolButton(text="…")
        browse.clicked.connect(self._browse)
        row.addWidget(self.matlab, 1)
        row.addWidget(browse)
        v.addLayout(row)

        buttons = QHBoxLayout()
        self.b_write = QPushButton("Write launcher")
        self.b_write.setToolTip("Write run_<study>.m next to the config")
        self.b_write.clicked.connect(self.write_launcher)
        self.b_validate = QPushButton("Validate in MATLAB")
        self.b_validate.clicked.connect(self.validate)
        self.b_run = QPushButton("Run locally")
        self.b_run.clicked.connect(self.run)
        self.b_submit = QPushButton("Submit to cluster")
        self.b_submit.clicked.connect(self.submit)
        self.scheduler = hpc_jobs.scheduler()
        self.b_submit.setToolTip(
            "PRESTUS submits the run to the scheduler with the config's hpc.* settings "
            "(prestus_pipeline_start); the jobs are listed below." if self.scheduler else
            "No sbatch or qsub on this computer: run the Nexus on a cluster login node to submit.")
        self.b_stop = QPushButton("Stop")
        self.b_stop.setEnabled(False)
        self.b_stop.clicked.connect(self.stop)
        self.b_save_log = QPushButton("Save log…")
        self.b_save_log.clicked.connect(self._save_log)
        for b in (self.b_write, self.b_validate, self.b_run, self.b_submit, self.b_stop):
            buttons.addWidget(b)
        buttons.addStretch(1)
        buttons.addWidget(self.b_save_log)
        v.addLayout(buttons)

        self.status = QLabel()
        self.status.setWordWrap(True)
        v.addWidget(self.status)
        self.log = QTextEdit(readOnly=True)
        mono = QFont("Menlo")
        mono.setStyleHint(QFont.StyleHint.Monospace)
        self.log.setFont(mono)
        v.addWidget(self.log, 1)

        # Submitted jobs (see jobs.py): polled every 30 s while any is open
        self.jobs: list[hpc_jobs.Job] = []
        jobs_row = QHBoxLayout()
        self.jobs_label = QLabel("Cluster jobs")
        jobs_row.addWidget(self.jobs_label)
        jobs_row.addStretch(1)
        self.b_jobs_refresh = QPushButton("Refresh")
        self.b_jobs_refresh.clicked.connect(self.refresh_jobs)
        self.b_jobs_cancel = QPushButton("Cancel selected")
        self.b_jobs_cancel.clicked.connect(self.cancel_jobs)
        jobs_row.addWidget(self.b_jobs_refresh)
        jobs_row.addWidget(self.b_jobs_cancel)
        v.addLayout(jobs_row)
        self.jobs_tree = QTreeWidget()
        self.jobs_tree.setHeaderLabels(["Job ID", "Name", "State"])
        self.jobs_tree.setRootIsDecorated(False)
        self.jobs_tree.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.jobs_tree.setFixedHeight(110)
        v.addWidget(self.jobs_tree)
        self.jobs_timer = QTimer(self)
        self.jobs_timer.setInterval(30_000)
        self.jobs_timer.timeout.connect(self.refresh_jobs)
        self._show_jobs()

    @property
    def model(self) -> ConfigModel:
        return self.window.model

    def refresh(self) -> None:
        if not self.matlab.text():
            hint = self.model.get("startup.matlab_bin_path")
            found = launcher.find_matlab(self.window.session.get("matlab_path"),
                                         hint if isinstance(hint, str) else None)
            self.matlab.setText(found or "")
        running = self.proc is not None
        for b in (self.b_write, self.b_validate, self.b_run):
            b.setEnabled(not running)
        self.b_submit.setEnabled(not running and self.scheduler is not None)
        self.b_stop.setEnabled(running)
        self.update_deployment_summary()

    def update_deployment_summary(self) -> None:
        """Cheap, no disk I/O — safe to call on every model edit, unlike refresh()."""
        if not self.deploy_toggle.isChecked():
            return
        plan = predict(self.model)
        lines = [f"DEPLOYMENT  ({plan.mode.replace('_', '-')} mode)"]
        lines += [f"  {j}" for j in plan.jobs]
        for w in plan.warnings:
            lines.append(f"  ! {w}")
        lines += ["", "SAVE CHOICES"]
        lines += [f"  {s}" for s in save_summary(self.model)]
        self.deploy_text.setPlainText("\n".join(lines))

    def _toggle_deployment(self, on: bool) -> None:
        self.window.session.set("run_show_deployment", on)
        self.deploy_text.setVisible(on)
        if on:
            self.update_deployment_summary()

    # ── helpers ─────────────────────────────────────────────────────────
    def _browse(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "MATLAB binary", self.matlab.text())
        if path:
            self.matlab.setText(path)

    def _ready(self, need_matlab: bool = True) -> bool:
        """Config saved and (optionally) MATLAB found; returns False after telling the user."""
        if self.window.dirty or self.model.path is None:
            if not self.window.save():
                return False
        if need_matlab and not Path(self.matlab.text()).is_file():
            QMessageBox.warning(self, "MATLAB", "Set the path to the matlab binary first.")
            return False
        if need_matlab:
            self.window.session.set("matlab_path", self.matlab.text())
        return True

    def _append(self, text: str, colour: str | None = None) -> None:
        for line in text.splitlines():
            c = colour
            if c is None:
                if line.startswith(launcher.FAIL_SENTINEL) or "Error" in line[:12] or line.startswith("{Error"):
                    c = "#c9302c"
                elif "Warning" in line[:10] or line.startswith("[Warning"):
                    c = "#b8860b"
                elif line.startswith(launcher.OK_SENTINEL):
                    c = "#3c763d"
            if launcher.OK_SENTINEL in line:
                self._ok = True
            safe = html.escape(line) or "&nbsp;"
            self.log.append(f"<span style='color:{c}'>{safe}</span>" if c else safe)
        self.log.moveCursor(QTextCursor.MoveOperation.End)

    def _start(self, task: str, argv: list[str], cwd: Path) -> None:
        self.task, self._ok = task, False
        self.log.clear()
        self._append("$ " + " ".join(argv[:2]) + " …", "#777")
        self.proc = QProcess(self)
        self.proc.setWorkingDirectory(str(cwd))
        self.proc.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self.proc.readyReadStandardOutput.connect(self._read)
        self.proc.finished.connect(self._done)
        self.proc.errorOccurred.connect(self._error)
        self.status.setText(f"<b>{task}</b> running…")
        self.proc.start(argv[0], argv[1:])
        self.refresh()

    def _read(self) -> None:
        if self.proc is not None:
            self._append(bytes(self.proc.readAllStandardOutput()).decode(errors="replace"))

    def _error(self, err) -> None:
        if err == QProcess.ProcessError.FailedToStart:
            self._append(f"Could not start {self.matlab.text()}", "#c9302c")
            self._done(-1, None)

    def _done(self, code: int, _status) -> None:
        if self.proc is None:
            return
        self._read()
        ok = code == 0 and self._ok
        colour = "#3c763d" if ok else "#c9302c"
        self.status.setText(f"<b>{self.task}</b>: <span style='color:{colour}'>"
                            f"{'finished' if ok else f'failed (exit code {code})'}</span>")
        self.proc.deleteLater()
        self.proc = None
        if self.task == "Submit":
            self._collect_jobs()
        self.refresh()
        self.finished.emit(self.task, ok)

    # ── actions ─────────────────────────────────────────────────────────
    def write_launcher(self) -> Path | None:
        if not self._ready(need_matlab=False):
            return None
        path = launcher.write_launcher(self.model.path, self.window.prestus_root)
        self.status.setText(f"Launcher written: <tt>{path}</tt> — run it in MATLAB with "
                            f"<tt>run('{path}')</tt>")
        return path

    def validate(self) -> None:
        if not self._ready():
            return
        argv = launcher.validate_command(self.matlab.text(), self.model.path, self.window.prestus_root)
        self._start("Validate", argv, self.model.path.parent)

    def run(self) -> None:
        errors = [i for i in self.window.issues if i.level == "error"]
        if errors and QMessageBox.question(
                self, "Run", f"{len(errors)} error(s) in the Issues panel. Run anyway?") \
                != QMessageBox.StandardButton.Yes:
            return
        if self.model.get("simulation.interactive") == 1:
            QMessageBox.warning(self, "Run", "simulation.interactive needs the MATLAB desktop; "
                                "set it to 0 or run the launcher from MATLAB.")
            return
        path = self.write_launcher()
        if path is None or not self._ready():
            return
        self._start("Run", launcher.run_command(self.matlab.text(), path), path.parent)

    def submit(self) -> None:
        """Submit with PRESTUS' own HPC path; the Nexus only follows the jobs."""
        if self.scheduler is None:
            return
        platform = self.model.get("platform")
        if platform == "matlab":
            QMessageBox.warning(self, "Submit", "platform is 'matlab' in this config, which runs "
                                "locally. Set it to 'auto', 'slurm' or 'qsub' (Setup tab) to submit.")
            return
        errors = [i for i in self.window.issues if i.level == "error"]
        if errors and QMessageBox.question(
                self, "Submit", f"{len(errors)} error(s) in the Issues panel. Submit anyway?") \
                != QMessageBox.StandardButton.Yes:
            return
        if not self._ready():
            return
        path = launcher.write_submit_launcher(self.model.path, self.window.prestus_root)
        self._start("Submit", launcher.run_command(self.matlab.text(), path), path.parent)

    def _collect_jobs(self) -> None:
        known = {j.job_id for j in self.jobs}
        new = [j for j in hpc_jobs.parse_submissions(self.log.toPlainText()) if j.job_id not in known]
        if not new:
            self._append("No submitted job was reported (see the log above).", "#b8860b")
        self.jobs.extend(new)
        self.refresh_jobs()

    def refresh_jobs(self) -> None:
        if self.jobs:
            hpc_jobs.query(self.jobs, self.scheduler)
        self._show_jobs()
        if any(not j.finished for j in self.jobs):
            self.jobs_timer.start()
        else:
            self.jobs_timer.stop()

    def _show_jobs(self) -> None:
        self.jobs_tree.clear()
        for j in self.jobs:
            item = QTreeWidgetItem([j.job_id, j.name, j.state])
            self.jobs_tree.addTopLevelItem(item)
        n_open = sum(1 for j in self.jobs if not j.finished)
        self.jobs_label.setText(f"Cluster jobs ({len(self.jobs)}, {n_open} open)" if self.jobs
                                else "Cluster jobs (none submitted from this window)")
        has = bool(self.jobs) and self.scheduler is not None
        self.b_jobs_refresh.setEnabled(has)
        self.b_jobs_cancel.setEnabled(has)

    def cancel_jobs(self) -> None:
        ids = [it.text(0) for it in self.jobs_tree.selectedItems()]
        if not ids:
            return
        if QMessageBox.question(self, "Cancel jobs", f"Cancel job(s) {', '.join(ids)}?") \
                != QMessageBox.StandardButton.Yes:
            return
        out = hpc_jobs.cancel(ids, self.scheduler)
        if out.strip():
            self._append(out, "#777")
        self.refresh_jobs()

    def stop(self) -> None:
        if self.proc is not None:
            self._append("Stopping…", "#777")
            self.proc.kill()

    def _save_log(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, "Save log", "planner_run.log", "Log (*.log *.txt)")
        if path:
            Path(path).write_text(self.log.toPlainText(), encoding="utf-8")
