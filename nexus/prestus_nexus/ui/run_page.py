"""Run page: write the launcher, validate in MATLAB, run locally or submit to the cluster
(here, or over SSH from a computer without a scheduler) with a live log, and follow
submitted jobs."""

from __future__ import annotations

import copy
import html
import shlex
import time
from concurrent.futures import Future, ThreadPoolExecutor
from pathlib import Path, PurePosixPath

from PyQt6.QtCore import QObject, QProcess, QProcessEnvironment, Qt, QTimer, pyqtSignal
from PyQt6.QtGui import QTextCursor
from PyQt6.QtWidgets import (QCheckBox, QFileDialog, QFormLayout, QHBoxLayout, QLabel, QLineEdit,
                             QMenu, QMessageBox, QPushButton, QTextEdit, QToolButton,
                             QVBoxLayout, QWidget)

from .. import jobs as hpc_jobs
from .. import launcher
from .. import progress as prog
from .. import proctree
from .. import remote as rmt
from .. import sequence as sq
from ..labels import plural
from ..model import ConfigModel
from . import icons, style
from .jobs_view import JobsView, JobsWindow
from .notify import notify
from .progress_panel import LiveViewPage, ProgressPanel
from .run_plan_card import RunPlanCard


# Scheduler queries, cancels and log reads run here: over SSH they take seconds, and the
# GUI thread must not wait for them. Results come back through a Relay signal.
_POOL = ThreadPoolExecutor(max_workers=3, thread_name_prefix="nexus-jobs")


def submit_background(fn, *args) -> Future:
    return _POOL.submit(fn, *args)


class Relay(QObject):
    """Hands a finished Future from a worker thread to the GUI thread."""
    done = pyqtSignal(object, object)         # tag, Future

    def deliver(self, tag, fut: Future) -> None:
        def emit(f, t=tag):
            try:
                self.done.emit(t, f)
            except RuntimeError:  # best-effort: the receiver was closed meanwhile
                pass
        fut.add_done_callback(emit)


class RunPage(QWidget):
    finished = pyqtSignal(str, bool)          # task, success
    jobs_changed = pyqtSignal()
    job_completed = pyqtSignal(str)           # job id of a cluster job that ended COMPLETED

    def __init__(self, window, parent: QWidget | None = None):
        super().__init__(parent)
        self.window = window
        self.proc: QProcess | None = None
        self.task = ""
        self._ok = False
        self._program = ""
        self._remote_host: rmt.RemoteHost | None = None

        v = QVBoxLayout(self)
        # 'What this run does' (run_plan_card): its own Run sub-tab, 'Plan' (main_window)
        self.plan_card = RunPlanCard()
        self.plan_card.refresh = lambda: self.plan_card.update_plan(self.model)   # TabPage._sub_tab_changed
        title = QLabel("Run in MATLAB")
        title.setObjectName("CardTitle")
        v.addWidget(title)
        intro = QLabel("Runs the saved study config in MATLAB on this computer, or submits it to the "
                       "cluster (on a login node with sbatch/qsub, or over SSH) and follows the jobs. "
                       "Validate first: it only loads the config with <tt>load_parameters</tt>.")
        intro.setWordWrap(True)
        v.addWidget(intro)

        row = QHBoxLayout()
        row.addWidget(QLabel("MATLAB"))
        self.matlab = QLineEdit()
        self.matlab.setPlaceholderText("path to the matlab binary")
        browse = icons.PathButton("file", "Choose the MATLAB binary")
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
            "Log in once with your password; submits, status checks, and cancels reuse this "
            f"connection until it has been idle for {rmt.CONTROL_PERSIST}. Not needed with an SSH key."
            if rmt.sharing_supported() else "Connection sharing needs OpenSSH on macOS/Linux; use an SSH key.")
        self.b_connect.clicked.connect(self.connect_ssh)
        self.b_disconnect = QPushButton("Disconnect")
        self.b_disconnect.clicked.connect(self.disconnect_ssh)
        self.b_probe = QPushButton("Test connection")
        self.b_probe.setToolTip("Report sbatch/qsub, MATLAB, and the PRESTUS folder on the host")
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
        # Progress of the run in MATLAB here (progress.py), from its output as it arrives
        self.progress = prog.RunProgress()
        self.progress_panel = ProgressPanel()
        self.progress_panel.live_requested.connect(self.open_live_view)
        v.addWidget(self.progress_panel)
        # simulation.live_view frames: the Run section's 'Live view' sub-tab (main_window adds
        # it, hidden until a run has frames), for this local run or a cluster job (jobs view)
        self.live_page = LiveViewPage()
        self._expected: dict | None = None        # predicted seconds per stage of the local run
        self._submit_expected: dict | None = None # ... of a whole-pipeline submission
        self._live_switched = False           # switched to 'Live view' for this local run
        self.progress_timer = QTimer(self)
        self.progress_timer.setInterval(1000)  # time-based estimates move between log lines
        self.progress_timer.timeout.connect(self._show_progress)
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
        self._jobs_poll: Future | None = None     # scheduler query in flight
        self._relay = Relay(self)
        self._relay.done.connect(self._background_done, Qt.ConnectionType.QueuedConnection)
        # Shown in the Run section's "Jobs" side tab (main_window) and, on demand, a window
        self.jobs_view = JobsView(self)
        self.jobs_window: JobsWindow | None = None
        # Jobs submitted earlier (also from an earlier session) are followed again
        self._restore_jobs()

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
        """Cheap — safe to call on every model edit, unlike refresh(). Disk access: a stat of
        run_stats.jsonl (re-read only when it changed) and the existence checks of
        run_plan (segmentation, pseudo-CT, measured CT)."""
        self.plan_card.update_plan(self.model)

    # ── SSH submission ─────────────────────────────────────────────────
    def _toggle_remote(self, on: bool) -> None:
        self.window.session.set("remote_enabled", on)
        self.remote_box.setVisible(on)
        self.window.update_context()          # scheduler settings show while submitting over SSH
        self.refresh()
        tab = self.window.tabs.currentWidget()
        if hasattr(tab, "refresh"):
            tab.refresh()

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
        target = self._launch_target()
        if target is None:
            return
        study, options = target
        run_model = self.model if study == self.model.path else \
            ConfigModel.from_file(self.model.schema, study, self.model.library)
        mapped, unmapped = rmt.path_overrides(run_model, r)
        seq_over: dict = {}
        for cfg in options.sequential_configs:
            p = Path(cfg).resolve()
            try:
                child = ConfigModel.from_file(self.model.schema, p, self.model.library)
            except Exception as e:  # reported: in a dialog; the user decides
                QMessageBox.warning(self, "Submit over SSH", f"Cannot read the sequential config {p}:\n{e}")
                return
            seq_over[str(p)], child_unmapped = rmt.path_overrides(child, r)
            unmapped += [f"{p.name}: {k}" for k in child_unmapped]
        if unmapped and QMessageBox.question(
                self, "Submit over SSH",
                "These config values are paths on this computer outside the mapped folder, so "
                "they cannot be translated for the cluster:\n\n" +
                "\n".join(unmapped[:8]) +
                "\n\nSubmit anyway?") != QMessageBox.StandardButton.Yes:
            return
        try:
            argv = rmt.submit_command(r, study, self.window.prestus_root, mapped,
                                      options=options, sequential_overrides=seq_over)
        except rmt.PathMapError as e:
            QMessageBox.warning(self, "SSH", f"{e}.\nThe study folder (and PRESTUS, unless a remote "
                                "PRESTUS folder is set) must be under the local folder.")
            return
        self.job_remote = r
        self._start("Submit", argv, self.model.path.parent)
        for key, value in mapped:
            self._append(f"  {key} → {value}  (on the cluster)", "#777")

    # ── helpers ─────────────────────────────────────────────────────────
    def _launch_target(self) -> tuple[Path, launcher.PipelineOptions] | None:
        """The config to run and its pipeline options: a sonication sequence simulated as a
        chain runs its generated step 1 with the other steps as its sequential chain
        (sequence.launch_target). None after telling the user why it cannot run."""
        try:
            study, options, notes = sq.launch_target(
                self.model, self._pipeline_options(),
                sq.load_backend(self.window.session, self.model.path))
        except sq.SequenceError as e:
            QMessageBox.warning(self, "Sonication sequence", str(e))
            return None
        for note in notes:
            self._append(note, "#777")
        return study, options

    def _pipeline_options(self):
        """Pipeline options stored for this study (Run tab ▸ Pipeline options)."""
        return launcher.load_options(self.window.session, self.model.path)

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
            # pre-wrap keeps runs of spaces (k-Wave's progress table stays aligned)
            css = f"white-space:pre-wrap;color:{c}" if c else "white-space:pre-wrap"
            self.log.append(f"<span style='{css}'>{safe}</span>")
        self.log.moveCursor(QTextCursor.MoveOperation.End)

    def _start(self, task: str, argv: list[str], cwd: Path, env: dict[str, str] | None = None) -> None:
        self.task, self._ok, self._program = task, False, argv[0]
        self._expected = None
        self.log.clear()
        self.progress = prog.RunProgress()
        self.progress_panel.clear()
        self._live_switched = False
        if self.live_page.source == "local":     # the last run's frames are gone with it
            self.live_page.source = None
            self._set_live_tab(False)
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
        proctree.new_session(self.proc)       # Stop ends the program and what it started
        self.status.setText(f"<b>{task}</b> running…")
        self.proc.start(argv[0], argv[1:])
        self.progress_timer.start()
        self.refresh()

    def _read(self) -> None:
        if self.proc is not None:
            text = bytes(self.proc.readAllStandardOutput()).decode(errors="replace")
            self.progress.feed(text, time.time())
            self._append(text)

    def _show_progress(self) -> None:
        """The local run's progress bar (Submit only submits: its jobs show theirs in Jobs)."""
        if self.task == "Submit":
            self.progress_panel.clear()
            return
        self.progress_panel.show_progress(self.progress, time.time(), expected=self._expected)
        has_frames = bool(self.progress_panel.frames or self.progress_panel.missing)
        # a running local run takes the Live view tab; else it keeps the job it shows
        if has_frames and (self.proc is not None or self.live_page.source in (None, "local")):
            self._show_local_live()
            if not self._live_switched:            # once per run: switch to it
                self._live_switched = True
                self._set_live_tab(True, select=True)

    def _show_local_live(self) -> None:
        name = self.model.path.name if self.model.path is not None else "run"
        self.live_page.show_run("local", self.progress, time.time(), f"Local run · {name}")
        self._set_live_tab(True)

    def open_live_view(self) -> None:
        """The card's Live view button: this local run's fields, in the Live view tab."""
        self._show_local_live()
        self._set_live_tab(True, select=True)

    def show_job_live(self, job, p, to_local=None, select: bool = False) -> None:
        """A cluster job's fields in the Live view tab (jobs view). A running local run keeps
        the tab unless the user asked for the job (select)."""
        if self.live_page.source == "local" and self.proc is not None and not select:
            return
        self.live_page.show_run(job.job_id, p, time.time(), f"Job {job.job_id} · {job.name}", to_local)
        self._set_live_tab(True, select)

    def _set_live_tab(self, visible: bool, select: bool = False) -> None:
        side = getattr(self.window, "run_side", None)
        i = side.indexOf(self.live_page) if side is not None else -1
        if i < 0:
            return
        side.setTabVisible(i, visible)
        if visible and select:
            side.setCurrentWidget(self.live_page)

    def _error(self, err) -> None:
        if err == QProcess.ProcessError.FailedToStart:
            self._append(f"Could not start {self._program}", style.active()["red"])
            self._done(-1, None)

    def _done(self, code: int, _status) -> None:
        if self.proc is None:
            return
        self._read()
        self.progress.flush(time.time())
        self.progress_timer.stop()
        if self.progress.has_news() and not (self.progress.finished or self.progress.failed):
            if code == 0:
                self.progress.finished, self.progress.open_stages = True, []
            else:
                self.progress.failed = True
        self._show_progress()
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
        if self.task not in ("Validate", "Submit", "Connect", "Disconnect", "Test connection"):
            notify(self.window, f"{self.task} {'finished' if ok else 'failed'}",
                   self.model.path.name if self.model.path is not None else "")
        self.finished.emit(self.task, ok)

    # ── actions ─────────────────────────────────────────────────────────
    def write_launcher(self) -> Path | None:
        if not self._ready(need_matlab=False):
            return None
        target = self._launch_target()
        if target is None:
            return None
        study, options = target
        path = launcher.write_launcher(study, self.window.prestus_root, options=options)
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
                self, "Run", f"{plural(len(errors), 'error')} in the Issues panel. Run anyway?") \
                != QMessageBox.StandardButton.Yes:
            return
        if self.model.get("simulation.interactive") == 1:
            QMessageBox.warning(self, "Run", "simulation.interactive needs the MATLAB desktop; "
                                "set it to 0 or run the launcher from MATLAB. To choose a PlanTUS "
                                "position in Workbench, set placement.plantus.interactive instead.")
            return
        path = self.write_launcher()
        if path is None or not self._ready():
            return
        self._start("Run", launcher.run_command(self.matlab.text(), path), path.parent)
        self._expected = self._predict_expected()

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
                                "locally. Set it to 'auto', 'slurm', or 'qsub' (Setup tab) to submit.")
            return
        errors = [i for i in self.window.issues if i.level == "error"]
        if errors and QMessageBox.question(
                self, "Submit", f"{plural(len(errors), 'error')} in the Issues panel. Submit anyway?") \
                != QMessageBox.StandardButton.Yes:
            return
        self._submit_expected = self._predict_expected()      # the whole pipeline
        if self._remote_on():
            self._submit_remote()
            return
        if not self._ready():
            return
        self.job_remote = None
        target = self._launch_target()
        if target is None:
            return
        study, options = target
        path = launcher.write_submit_launcher(study, self.window.prestus_root, options=options)
        self._start("Submit", launcher.run_command(self.matlab.text(), path), path.parent)

    def remote_host(self) -> rmt.RemoteHost | None:
        """The SSH host when Submit goes over SSH (asks for missing fields), else None."""
        return self._remote() if self._remote_on() else None

    def submit_kwave_drive(self, script: Path, fit_mat: Path, sim_config: Path | None,
                           remote: rmt.RemoteHost | None = None) -> bool:
        """Submit a k-Wave drive-check script as a scheduler job (kwave_drive_submit): to this
        computer's scheduler, or over SSH when remote is given (script written with the host's
        paths). The job submits one simulation per depth and waits for them on a compute node.
        Returns whether it started."""
        self._submit_expected = None           # not a whole pipeline: no whole-run prediction
        if self.proc is not None:
            QMessageBox.warning(self, "Submit", f"“{self.task}” is still running on the Run tab.")
            return False
        script = Path(script)
        if remote is not None:
            if not self._ready(need_matlab=False):
                return False
            try:
                r_script, r_mat = remote.to_remote(script), remote.to_remote(fit_mat)
                r_cfg = remote.to_remote(sim_config) if sim_config else ""
            except rmt.PathMapError as e:
                QMessageBox.warning(self, "SSH", f"{e}.")
                return False
            body = launcher.kwave_drive_submit_body(r_script, r_mat, r_cfg,
                                                    remote.remote_root(self.window.prestus_root))
            r_dir = str(PurePosixPath(r_script).parent)
            argv = remote.ssh(f"cd {shlex.quote(r_dir)} && {remote.matlab} -batch "
                              f"{shlex.quote(launcher._wrapped(body))}")
            self.job_remote = remote
        else:
            if self.scheduler is None:
                QMessageBox.warning(self, "Submit", "No SLURM or qsub scheduler on this computer. "
                                    "Turn on Submit over SSH (Run tab), or run the script on the cluster.")
                return False
            if not self._ready():
                return False
            path = launcher.write_kwave_drive_submit(script, fit_mat, sim_config, self.window.prestus_root)
            argv = launcher.run_command(self.matlab.text(), path)
            self.job_remote = None
        self._start("Submit", argv, script.parent)
        self._append("k-Wave drive check as a job: it submits one free-water simulation per depth "
                     "and waits for them, then saves the corrected entry.", "#777")
        return True

    def submits_to_cluster(self) -> bool:
        """Whether a run from this page goes to a scheduler (here or over SSH)."""
        return self._cluster_first()

    def submit_placement(self, mode: str = "plantus") -> bool:
        """Submit only the placement stage of the config (PlanTUS, TUS_entry, heuristic) as a
        job, no simulation. Returns whether it started."""
        return self.submit_partial(mode, f"the {mode} placement", launcher.placement_overrides(mode),
                                   launcher.PLACEMENT_ONLY)

    def submit_partial(self, name: str, what: str, overrides: list, extra: str = "") -> bool:
        """Submit part of the pipeline (launcher.partial_script: ``overrides`` on the config,
        ``extra`` statements) as a job, the way Submit does: to this computer's scheduler, or
        over SSH. Returns whether it started."""
        self._submit_expected = None           # part of the pipeline: no whole-run prediction
        if self.proc is not None:
            QMessageBox.warning(self, "Submit", f"“{self.task}” is still running on the Run tab.")
            return False
        if self.model.get("platform") == "matlab":
            QMessageBox.warning(self, "Submit", "platform is 'matlab' in this config, which runs "
                                "locally. Set it to 'auto', 'slurm', or 'qsub' (Setup tab) to submit.")
            return False
        if self._remote_on():
            r = self._remote()
            if r is None or not self._ready(need_matlab=False):
                return False
            mapped, unmapped = rmt.path_overrides(self.model, r)
            if unmapped and QMessageBox.question(
                    self, "Submit over SSH",
                    "These config values are paths on this computer outside the mapped folder:\n\n"
                    + "\n".join(unmapped[:8]) + "\n\nSubmit anyway?") != QMessageBox.StandardButton.Yes:
                return False
            try:
                argv = rmt.submit_command(r, self.model.path, self.window.prestus_root,
                                          mapped + list(overrides), extra=extra)
            except rmt.PathMapError as e:
                QMessageBox.warning(self, "SSH", f"{e}.")
                return False
            self.job_remote = r
        else:
            if self.scheduler is None or not self._ready():
                return False
            self.job_remote = None
            path = launcher.write_partial_script(self.model.path, self.window.prestus_root, name, what,
                                                 list(overrides), extra, local=False)
            argv = launcher.run_command(self.matlab.text(), path)
        self._start("Submit", argv, self.model.path.parent)
        self._append(f"Only {what} (no simulation).", "#777")
        return True

    def _collect_jobs(self) -> None:
        known = {j.job_id for j in self.jobs}
        new = [j for j in hpc_jobs.parse_submissions(self.log.toPlainText()) if j.job_id not in known]
        if not new:
            self._append("No submitted job was reported (see the log above).", style.active()["amber"])
        elif self.job_remote is not None:
            self.job_kind = rmt.scheduler_from_log(self.log.toPlainText()) or "slurm"
        else:
            self.job_kind = self.scheduler
        for j in new:
            j.submitted = time.time()
        if len(new) == 1 and self._submit_expected:     # one job runs the whole pipeline
            new[0].expected = dict(self._submit_expected)
        self._submit_expected = None
        self.jobs.extend(new)
        self._store_jobs()
        self.refresh_jobs()
        if new:
            self.show_jobs_tab(new[0].job_id)

    def _predict_expected(self) -> dict | None:
        """Predicted seconds per stage of a pipeline run of this config (workload.py), for
        the whole-run progress; None when the prediction is incomplete."""
        try:
            from .. import workload
            pred = workload.predict(self.model)
            if pred.total() is None:
                return None
            return workload.expected_seconds(pred) or None
        except Exception:  # best-effort: the prediction is a convenience, never a blocker
            return None

    # ── the job list across sessions ────────────────────────────────────
    def _store_jobs(self) -> None:
        self.window.session.set(hpc_jobs.STORE_KEY, hpc_jobs.to_store(
            self.jobs, self.job_kind, self.job_remote is not None))
        self.window.session.save()

    def _restore_jobs(self) -> None:
        jobs, kind, remote = hpc_jobs.from_store(self.window.session.get(hpc_jobs.STORE_KEY))
        if not jobs:
            return
        self.jobs, self.job_kind = jobs, kind or self.job_kind
        self.job_remote = self._remote_from_session() if remote else None
        if any(not j.finished for j in self.jobs):
            QTimer.singleShot(0, self.refresh_jobs)

    def _remote_from_session(self) -> rmt.RemoteHost | None:
        """The SSH host of the session's settings, without asking (None when not set)."""
        sess = self.window.session
        host = (sess.get("remote_host") or "").strip()
        if not host:
            return None
        local, remote = sess.get("remote_local_prefix") or "", sess.get("remote_remote_prefix") or ""
        return rmt.RemoteHost(host=host, path_map=[(local, remote)] if local and remote else [],
                              setup=sess.get("remote_setup") or "", matlab=sess.get("remote_matlab") or "matlab",
                              prestus_root=sess.get("remote_prestus_root") or "")

    def clear_finished_jobs(self) -> None:
        """Forget the jobs that have ended (the Jobs tab's 'Clear finished')."""
        self.jobs = [j for j in self.jobs if not j.finished]
        self._store_jobs()
        self.jobs_changed.emit()

    def refresh_jobs(self) -> None:
        """Ask the scheduler for the jobs' state in the background (one query at a time);
        jobs_changed follows when the answer is in."""
        if self._jobs_poll is not None:
            return
        if not self.jobs:
            self._jobs_polled([])
            return
        snapshot = copy.deepcopy(self.jobs)       # the worker updates copies, never the shown jobs
        self._jobs_poll = submit_background(hpc_jobs.poll, snapshot, self.job_kind, self.job_remote)
        self._relay.deliver("query", self._jobs_poll)

    def jobs_busy(self) -> bool:
        """A scheduler query or cancel is in flight."""
        return self._jobs_poll is not None

    def _background_done(self, tag, fut: Future) -> None:
        if tag == "query":
            self._jobs_poll = None
            try:
                updated = fut.result()
            except Exception as exc:  # reported: in the run log; the last known state is kept
                self._append(f"Could not query the scheduler: {exc}", style.active()["amber"])
                updated = []
            self._jobs_polled(updated)
        elif tag == "cancel":
            self._jobs_poll = None
            try:
                out = fut.result()
            except Exception as exc:  # reported: in the run log
                out = f"Cancel failed: {exc}"
            if out.strip():
                self._append(out, "#777")
            self.refresh_jobs()

    def _jobs_polled(self, updated: list[hpc_jobs.Job]) -> None:
        fresh = {j.job_id: j for j in updated}
        completed, ended = [], []
        for job in self.jobs:                     # jobs added meanwhile keep their state
            if job.job_id in fresh:
                was_open = not job.finished
                job.__dict__.update(fresh[job.job_id].__dict__)
                if was_open and job.finished:
                    ended.append(job)
                    if job.state.startswith("COMPLETED"):
                        completed.append(job.job_id)
        # jobs that the followed jobs submitted themselves (drive check, pipeline stages)
        children = hpc_jobs.child_jobs(self.jobs)
        self.jobs.extend(children)
        if updated or children:
            self._store_jobs()
        self.jobs_changed.emit()
        for job_id in completed:                  # Results tab: read the new outputs
            self.job_completed.emit(job_id)
        if ended:
            failed = [j for j in ended if not j.state.startswith("COMPLETED")]
            text = ", ".join(f"{j.job_id} {j.state.split()[0].lower()}" for j in ended[:6])
            notify(self.window, f"{plural(len(ended), 'cluster job')} ended"
                   + (f" ({len(failed)} not completed)" if failed else ""), text)
        if children:
            QTimer.singleShot(0, self.refresh_jobs)   # their state and logs now
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
        if QMessageBox.question(self, "Cancel jobs", f"Cancel {'job' if len(ids) == 1 else 'jobs'} {', '.join(ids)}?") \
                != QMessageBox.StandardButton.Yes:
            return
        if self._jobs_poll is not None:           # after the query in flight
            QTimer.singleShot(300, lambda: self._cancel_now(ids))
            return
        self._cancel_now(ids)

    def _cancel_now(self, ids: list[str]) -> None:
        if self._jobs_poll is not None:
            QTimer.singleShot(300, lambda: self._cancel_now(ids))
            return
        self._append(f"Cancelling {'job' if len(ids) == 1 else 'jobs'} {', '.join(ids)} …", "#777")
        self._jobs_poll = submit_background(hpc_jobs.cancel, ids, self.job_kind, self.job_remote)
        self._relay.deliver("cancel", self._jobs_poll)

    def stop(self) -> None:
        """End the running program and everything it started (MATLAB runs under a launcher)."""
        if self.proc is not None:
            self._append("Stopping…", "#777")
            proctree.kill_tree(int(self.proc.processId() or 0))
            if self.proc is not None:
                self.proc.kill()

    def _save_log(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, "Save log", "planner_run.log", "Log (*.log *.txt)")
        if path:
            Path(path).write_text(self.log.toPlainText(), encoding="utf-8")
