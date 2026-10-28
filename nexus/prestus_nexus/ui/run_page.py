"""Run page: write the launcher, validate in MATLAB, run locally or submit to the cluster
(here, or over SSH from a computer without a scheduler) with a live log, and follow
submitted jobs."""

from __future__ import annotations

import html
from pathlib import Path

from PyQt6.QtCore import QProcess, QProcessEnvironment, QTimer, pyqtSignal
from PyQt6.QtGui import QTextCursor
from PyQt6.QtWidgets import (QCheckBox, QFileDialog, QFormLayout, QHBoxLayout, QLabel, QLineEdit,
                             QMenu, QMessageBox, QPlainTextEdit, QPushButton, QTextEdit, QToolButton,
                             QVBoxLayout, QWidget)

from .. import jobs as hpc_jobs
from .. import launcher
from .. import remote as rmt
from ..model import ConfigModel
from ..estimate import estimate
from ..predict_outputs import predict, save_summary
from . import style
from .jobs_view import JobsView, JobsWindow


class RunPage(QWidget):
    finished = pyqtSignal(str, bool)          # task, success
    jobs_changed = pyqtSignal()

    def __init__(self, window, parent: QWidget | None = None):
        super().__init__(parent)
        self.window = window
        self.proc: QProcess | None = None
        self.task = ""
        self._ok = False
        self._program = ""
        self._remote_host: rmt.RemoteHost | None = None

        v = QVBoxLayout(self)
        title = QLabel("Run in MATLAB")
        title.setObjectName("CardTitle")
        v.addWidget(title)
        intro = QLabel("Runs the saved study config in MATLAB on this computer, or submits it to the "
                       "cluster (on a login node with sbatch/qsub, or over SSH) and follows the jobs. "
                       "Validate first: it only loads the config with <tt>load_parameters</tt>.")
        intro.setWordWrap(True)
        v.addWidget(intro)

        self.deploy_toggle = QCheckBox("Show deployment, save && resource details")
        self.deploy_toggle.setToolTip(
            "What prestus_pipeline_start would submit (single job, or one job per transducer/"
            "target/uncertainty variant) and what io.save_* keeps besides the CSV/report.")
        self.deploy_toggle.setChecked(bool(window.session.get("run_show_deployment", False)))
        self.deploy_toggle.toggled.connect(self._toggle_deployment)
        v.addWidget(self.deploy_toggle)
        self.deploy_text = QPlainTextEdit(readOnly=True)
        self.deploy_text.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        deploy_mono = style.mono_font()
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

        # Submit over SSH (remote.py): for a computer without sbatch/qsub that sees the
        # study folder on a shared drive
        self.remote_toggle = QCheckBox("Submit over SSH")
        self.remote_toggle.setToolTip(
            "Submit from this computer by running PRESTUS on a cluster login node over SSH. "
            "The study folder must be on a drive the cluster also sees; paths are translated "
            "by the folder mapping. Log in with Connect… (password, once) or an SSH key.")
        self.remote_toggle.setChecked(bool(window.session.get("remote_enabled", False)))
        self.remote_toggle.toggled.connect(self._toggle_remote)
        v.addWidget(self.remote_toggle)
        self.remote_box = QWidget()
        form = QFormLayout(self.remote_box)
        form.setContentsMargins(18, 0, 0, 0)
        sess = window.session
        self.r_host = QLineEdit(sess.get("remote_host", ""))
        self.r_host.setPlaceholderText("user@mentat005.dccn.nl, or an alias from ~/.ssh/config")
        self.r_host.textChanged.connect(lambda _t: self.refresh())
        self.r_local = QLineEdit(sess.get("remote_local_prefix", ""))
        self.r_local.setPlaceholderText("folder on this computer, e.g. /Volumes/2425122.01")
        self.r_remote = QLineEdit(sess.get("remote_remote_prefix", ""))
        self.r_remote.setPlaceholderText("the same folder on the cluster, e.g. /project/2425122.01")
        self.r_setup = QLineEdit(sess.get("remote_setup", ""))
        self.r_setup.setPlaceholderText("e.g. module load slurm matlab/R2024a")
        self.r_matlab = QLineEdit(sess.get("remote_matlab", ""))
        self.r_matlab.setPlaceholderText("matlab")
        self.r_root = QLineEdit(sess.get("remote_prestus_root", ""))
        self.r_root.setPlaceholderText("empty: this PRESTUS checkout, through the folder mapping")
        form.addRow("Host", self.r_host)
        form.addRow("Local folder", self.r_local)
        form.addRow("is on the cluster", self.r_remote)
        form.addRow("Setup commands", self.r_setup)
        form.addRow("Remote MATLAB", self.r_matlab)
        form.addRow("Remote PRESTUS", self.r_root)
        self.b_connect = QPushButton("Connect…")
        self.b_connect.setToolTip(
            "Log in once with your password; submits, status checks and cancels reuse this "
            f"connection until it has been idle for {rmt.CONTROL_PERSIST}. Not needed with an SSH key."
            if rmt.sharing_supported() else "Connection sharing needs OpenSSH on macOS/Linux; use an SSH key.")
        self.b_connect.clicked.connect(self.connect_ssh)
        self.b_disconnect = QPushButton("Disconnect")
        self.b_disconnect.clicked.connect(self.disconnect_ssh)
        self.b_probe = QPushButton("Test connection")
        self.b_probe.setToolTip("Report sbatch/qsub, MATLAB and the PRESTUS folder on the host")
        self.b_probe.clicked.connect(self.test_connection)
        self.r_status = QLabel()
        conn_row = QHBoxLayout()
        for b in (self.b_connect, self.b_disconnect, self.b_probe):
            conn_row.addWidget(b)
        conn_row.addWidget(self.r_status, 1)
        form.addRow("", conn_row)
        self.remote_box.setVisible(self.remote_toggle.isChecked())
        v.addWidget(self.remote_box)

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
        # One primary action, chosen by 'Run on' (refresh); the other way to run, the launcher
        # and the log export sit in 'More' so the row reads as: run · check · stop.
        self.b_more = QToolButton(text="More ▾")
        self.b_more.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self.b_more.setStyleSheet("QToolButton::menu-indicator { image: none; width: 0; }")
        more = QMenu(self.b_more)
        self.act_alt = more.addAction("Submit to cluster")
        self.act_alt.triggered.connect(lambda: (self.submit if self.b_run.property("primary") else self.run)())
        more.addSeparator()
        more.addAction("Write launcher (run_<study>.m)", self.write_launcher)
        more.addAction("Save log…", self._save_log)
        self.b_more.setMenu(more)
        for b in (self.b_run, self.b_submit, self.b_validate, self.b_stop):
            buttons.addWidget(b)
        buttons.addStretch(1)
        buttons.addWidget(self.b_more)
        for b in (self.b_write, self.b_save_log):    # reachable from 'More'; kept for scripting
            b.setVisible(False)
            buttons.addWidget(b)
        v.addLayout(buttons)

        self.status = QLabel()
        self.status.setWordWrap(True)
        v.addWidget(self.status)
        self.log = QTextEdit(readOnly=True)
        mono = style.mono_font()
        self.log.setFont(mono)
        v.addWidget(self.log, 1)

        # Submitted jobs (see jobs.py): polled every 30 s while any is open, with the
        # scheduler and (for SSH submissions) the host they were submitted on
        self.jobs: list[hpc_jobs.Job] = []
        self.job_kind: str | None = self.scheduler
        self.job_remote: rmt.RemoteHost | None = None
        self.jobs_timer = QTimer(self)
        self.jobs_timer.setInterval(30_000)
        self.jobs_timer.timeout.connect(self.refresh_jobs)
        # Shown in the Run section's "Jobs" side tab (main_window) and, on demand, a window
        self.jobs_view = JobsView(self)
        self.jobs_window: JobsWindow | None = None

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
        self.b_submit.setEnabled(not running and (self.scheduler is not None or self._remote_on()))
        for b in (self.b_probe, self.b_disconnect):
            b.setEnabled(not running)
        self.b_connect.setEnabled(not running and rmt.sharing_supported())
        self.b_stop.setEnabled(running)
        self._arrange_run_buttons(running)
        if self._remote_on():
            self.b_submit.setToolTip("Runs PRESTUS on the SSH host, which submits with the config's "
                                     "hpc.* settings; the jobs are listed below.")
        if self.model.path is not None and not self.r_local.text():
            self._guess_mapping(self.model.path)
        self.update_deployment_summary()

    def _cluster_first(self) -> bool:
        """Whether 'Run on' sends this config to a scheduler (else it runs in MATLAB here)."""
        platform = self.model.get("platform")
        if platform in ("slurm", "qsub"):
            return True
        return platform == "auto" and (self.scheduler is not None or self._remote_on())

    def _arrange_run_buttons(self, running: bool) -> None:
        cluster = self._cluster_first()
        primary, other = (self.b_submit, self.b_run) if cluster else (self.b_run, self.b_submit)
        primary.setVisible(True)
        other.setVisible(False)
        for b, on in ((primary, True), (other, False)):
            if b.property("primary") != on:
                b.setProperty("primary", on)
                b.style().unpolish(b)
                b.style().polish(b)
        self.act_alt.setText(other.text())
        self.act_alt.setEnabled(other.isEnabled())
        self.act_alt.setToolTip(other.toolTip())

    def update_deployment_summary(self) -> None:
        """Cheap — safe to call on every model edit, unlike refresh(). The only disk
        access is a stat of run_stats.jsonl (re-read only when it changed)."""
        if not self.deploy_toggle.isChecked():
            return
        plan = predict(self.model)
        lines = [f"DEPLOYMENT  ({plan.mode.replace('_', '-')} mode)"]
        lines += [f"  {j}" for j in plan.jobs]
        for w in plan.warnings:
            lines.append(f"  ! {w}")
        lines += ["", "SAVE CHOICES"]
        lines += [f"  {s}" for s in save_summary(self.model)]
        lines += ["", "RESOURCES  (from earlier runs of this project)"]
        lines += [f"  {s}" for s in estimate(self.model).lines()]
        self.deploy_text.setPlainText("\n".join(lines))

    def _toggle_deployment(self, on: bool) -> None:
        self.window.session.set("run_show_deployment", on)
        self.deploy_text.setVisible(on)
        if on:
            self.update_deployment_summary()

    # ── SSH submission ─────────────────────────────────────────────────
    def _toggle_remote(self, on: bool) -> None:
        self.window.session.set("remote_enabled", on)
        self.remote_box.setVisible(on)
        self.refresh()

    def _remote_on(self) -> bool:
        return self.remote_toggle.isChecked() and bool(self.r_host.text().strip())

    def _guess_mapping(self, study: Path) -> None:
        """A study on /Volumes/<share> is usually /project/<share> on the Donders cluster."""
        parts = Path(study).resolve().parts
        if len(parts) > 2 and parts[1] == "Volumes":
            self.r_local.setText(f"/Volumes/{parts[2]}")
            if not self.r_remote.text():
                self.r_remote.setText(f"/project/{parts[2]}")

    def _remote(self, need_mapping: bool = True) -> rmt.RemoteHost | None:
        """The SSH host from the fields (saved to the session), or None after telling the user."""
        host = self.r_host.text().strip()
        local, remote = self.r_local.text().strip(), self.r_remote.text().strip()
        if not host or (need_mapping and (not local or not remote)):
            QMessageBox.warning(self, "SSH", "Set the host" + (" and the folder mapping" if need_mapping
                                                               else "") + " first.")
            return None
        values = {"remote_host": host, "remote_local_prefix": local, "remote_remote_prefix": remote,
                  "remote_setup": self.r_setup.text().strip(),
                  "remote_matlab": self.r_matlab.text().strip(),
                  "remote_prestus_root": self.r_root.text().strip()}
        for k, val in values.items():
            self.window.session.set(k, val)
        return rmt.RemoteHost(host=host, path_map=[(local, remote)] if local and remote else [],
                              setup=values["remote_setup"], matlab=values["remote_matlab"] or "matlab",
                              prestus_root=values["remote_prestus_root"])

    def connect_ssh(self) -> None:
        """Log in once (password in a Nexus dialog) and keep a shared connection open."""
        r = self._remote(need_mapping=False)
        if r is None:
            return
        self._remote_host = r
        self._start("Connect", rmt.connect_command(r), Path.home(), env=rmt.askpass_env())

    def disconnect_ssh(self) -> None:
        r = self._remote(need_mapping=False)
        if r is None:
            return
        self._remote_host = r
        self._start("Disconnect", rmt.disconnect_command(r), Path.home())

    def _connection_done(self, code: int) -> None:
        host = self._remote_host.host if self._remote_host else ""
        if self.task == "Connect":
            self._ok = code == 0
            self.r_status.setText(f"Connected to {host}" if self._ok else "Not connected")
        elif self.task == "Disconnect":
            self._ok = True                       # 'no master running' is fine too
            self.r_status.setText("Not connected")

    def test_connection(self) -> None:
        r = self._remote()
        if r is None:
            return
        try:
            argv = rmt.probe_command(r, self.window.prestus_root)
        except rmt.PathMapError as e:
            QMessageBox.warning(self, "SSH", f"PRESTUS: {e}. Set the remote PRESTUS folder.")
            return
        self._start("Test connection", argv, Path.home())

    def _probe_summary(self) -> None:
        found = rmt.parse_probe(self.log.toPlainText())
        if "host" not in found:
            return
        missing = [k for k in ("matlab", "prestus") if found.get(k, "-").startswith("-")]
        if found.get("sbatch", "-") == "-" and found.get("qsub", "-") == "-":
            missing.insert(0, "sbatch/qsub")
        self._ok = not missing
        if missing:
            self._append(f"Not found on {found['host']}: {', '.join(missing)} "
                         "(check the setup commands and the remote PRESTUS folder).", style.active()["red"])
        else:
            self._append(f"Ready to submit from {found['host']}.", style.active()["green"])

    def _submit_remote(self) -> None:
        r = self._remote()
        if r is None or not self._ready(need_matlab=False):
            return
        mapped, unmapped = rmt.path_overrides(self.model, r)
        if unmapped and QMessageBox.question(
                self, "Submit over SSH",
                "These config values are paths on this computer outside the mapped folder, so "
                "they cannot be translated for the cluster:\n\n" +
                "\n".join(unmapped[:8]) +
                "\n\nSubmit anyway?") != QMessageBox.StandardButton.Yes:
            return
        try:
            argv = rmt.submit_command(r, self.model.path, self.window.prestus_root, mapped)
        except rmt.PathMapError as e:
            QMessageBox.warning(self, "SSH", f"{e}.\nThe study folder (and PRESTUS, unless a remote "
                                "PRESTUS folder is set) must be under the local folder.")
            return
        self.job_remote = r
        self._start("Submit", argv, self.model.path.parent)
        for key, value in mapped:
            self._append(f"  {key} → {value}  (on the cluster)", "#777")

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
                    c = style.active()["red"]
                elif "Warning" in line[:10] or line.startswith("[Warning"):
                    c = style.active()["amber"]
                elif line.startswith(launcher.OK_SENTINEL):
                    c = style.active()["green"]
            if launcher.OK_SENTINEL in line:
                self._ok = True
            safe = html.escape(line) or "&nbsp;"
            # Always a span: a bare escaped line is taken as plain text and keeps its
            # entities (&quot;), which breaks parse_submissions on the log
            self.log.append(f"<span style='color:{c}'>{safe}</span>" if c else f"<span>{safe}</span>")
        self.log.moveCursor(QTextCursor.MoveOperation.End)

    def _start(self, task: str, argv: list[str], cwd: Path, env: dict[str, str] | None = None) -> None:
        self.task, self._ok, self._program = task, False, argv[0]
        self.log.clear()
        shown = argv[:2] if argv[0] != "ssh" else ["ssh", rmt.host_of(argv)]
        self._append("$ " + " ".join(shown) + " …", "#777")
        self.proc = QProcess(self)
        self.proc.setWorkingDirectory(str(cwd))
        if env:
            penv = QProcessEnvironment.systemEnvironment()
            for k, val in env.items():
                penv.insert(k, val)
            self.proc.setProcessEnvironment(penv)
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
            self._append(f"Could not start {self._program}", style.active()["red"])
            self._done(-1, None)

    def _done(self, code: int, _status) -> None:
        if self.proc is None:
            return
        self._read()
        if self.task == "Test connection":
            self._ok = code == 0
            self._probe_summary()
        self._connection_done(code)
        # ssh's own refusal reads "Permission denied (publickey,password)."; a bare
        # "Permission denied" comes from the remote side (e.g. MATLAB's mkdir)
        if self._program == "ssh" and code != 0 and "Permission denied (" in self.log.toPlainText():
            self._append("Not logged in: click Connect… to log in with your password, "
                         "or set up an SSH key.", style.active()["amber"])
        if self._program == "ssh" and "DISPLAY" in self.log.toPlainText():
            # site 'matlab' wrappers (e.g. Donders' module) refuse to start without X
            self._append("The host's 'matlab' is a wrapper that needs a display. Set Remote MATLAB "
                         "to the MATLAB binary itself (e.g. <matlabroot>/bin/matlab; on the host: "
                         "module load matlab/<version>; type -a matlab).", style.active()["amber"])
        ok = code == 0 and self._ok
        colour = style.active()["green"] if ok else style.active()["red"]
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

    def run_script(self, task: str, script: Path) -> bool:
        """Run a generated MATLAB script (Localite batch ingestion, k-Wave drive check)
        with the live log; the config is saved first."""
        if self.proc is not None:
            QMessageBox.warning(self, task, f"MATLAB is still running “{self.task}” on the Run tab.")
            return False
        if not self._ready():
            return False
        self._start(task, launcher.run_command(self.matlab.text(), Path(script)), Path(script).parent)
        return True

    def submit(self) -> None:
        """Submit with PRESTUS' own HPC path; the Nexus only follows the jobs."""
        if self.scheduler is None and not self._remote_on():
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
        if self._remote_on():
            self._submit_remote()
            return
        if not self._ready():
            return
        self.job_remote = None
        path = launcher.write_submit_launcher(self.model.path, self.window.prestus_root)
        self._start("Submit", launcher.run_command(self.matlab.text(), path), path.parent)

    def _collect_jobs(self) -> None:
        known = {j.job_id for j in self.jobs}
        new = [j for j in hpc_jobs.parse_submissions(self.log.toPlainText()) if j.job_id not in known]
        if not new:
            self._append("No submitted job was reported (see the log above).", style.active()["amber"])
        elif self.job_remote is not None:
            self.job_kind = rmt.scheduler_from_log(self.log.toPlainText()) or "slurm"
        else:
            self.job_kind = self.scheduler
        self.jobs.extend(new)
        self.refresh_jobs()
        if new:
            self.show_jobs_tab(new[0].job_id)

    def refresh_jobs(self) -> None:
        if self.jobs:
            hpc_jobs.query(self.jobs, self.job_kind, self.job_remote)
        self.jobs_changed.emit()
        if any(not j.finished for j in self.jobs):
            self.jobs_timer.start()
        else:
            self.jobs_timer.stop()

    def show_jobs_tab(self, job_id: str | None = None) -> None:
        """Bring the Jobs side tab forward, with job_id selected."""
        side = getattr(self.window, "run_side", None)
        if side is not None:
            side.setCurrentWidget(self.jobs_view)
        self.jobs_view.select(job_id)

    def pop_out_jobs(self, job_id: str | None = None) -> None:
        """The jobs and their logs in a separate, resizable window."""
        if self.jobs_window is None:
            self.jobs_window = JobsWindow(self)
        self.jobs_window.view.select(job_id or self.jobs_view.selected_id())
        self.jobs_window.show()
        self.jobs_window.raise_()
        self.jobs_window.activateWindow()

    def cancel_jobs(self, ids: list[str]) -> None:
        if not ids:
            return
        if QMessageBox.question(self, "Cancel jobs", f"Cancel job(s) {', '.join(ids)}?") \
                != QMessageBox.StandardButton.Yes:
            return
        out = hpc_jobs.cancel(ids, self.job_kind, self.job_remote)
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
