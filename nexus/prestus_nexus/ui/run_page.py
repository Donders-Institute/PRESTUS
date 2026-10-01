"""Run page: write the launcher, validate in MATLAB, run locally with a live log."""

from __future__ import annotations

import html
from pathlib import Path

from PyQt6.QtCore import QProcess, pyqtSignal
from PyQt6.QtGui import QFont, QTextCursor
from PyQt6.QtWidgets import (QCheckBox, QFileDialog, QHBoxLayout, QLabel, QLineEdit, QMessageBox,
                             QPlainTextEdit, QPushButton, QTextEdit, QToolButton, QVBoxLayout, QWidget)

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
        intro = QLabel("Runs the saved study config in MATLAB on this computer. "
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
        self.b_stop = QPushButton("Stop")
        self.b_stop.setEnabled(False)
        self.b_stop.clicked.connect(self.stop)
        self.b_save_log = QPushButton("Save log…")
        self.b_save_log.clicked.connect(self._save_log)
        for b in (self.b_write, self.b_validate, self.b_run, self.b_stop):
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

    def stop(self) -> None:
        if self.proc is not None:
            self._append("Stopping…", "#777")
            self.proc.kill()

    def _save_log(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, "Save log", "planner_run.log", "Log (*.log *.txt)")
        if path:
            Path(path).write_text(self.log.toPlainText(), encoding="utf-8")
