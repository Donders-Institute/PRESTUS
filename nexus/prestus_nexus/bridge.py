"""MATLAB bridge: one persistent MATLAB session (functions/gui/nexus_bridge.m) that serves
the Nexus's interactive MATLAB steps -- Localite detect/write, transducer calibration
predict/fit/save -- without a MATLAB start-up per click.

Protocol (all in one bridge folder; see nexus_bridge.m):

  req/<id>.json       {"id", "fn", "args"}; written to a temporary name, then renamed
  resp/<id>.json      {"id", "ok", "result" | "error", "identifier", "report", "log"}
  progress/<id>.txt   last progress message of a running request
  ready.json          the session can take requests
  heartbeat           touched by the Nexus; MATLAB exits when it goes stale
  stop                ask the session to end

``BridgeClient`` is the file protocol alone (no Qt, used by the tests); ``MatlabBridge``
starts MATLAB with QProcess, polls on a timer and calls back on the GUI thread.
Long work (batch ingestion, k-Wave) does not use the bridge; it runs on the Run tab.
"""

from __future__ import annotations

import itertools
import json
import shutil
import tempfile
import threading
import time
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

from .fsio import write_text_atomic

ALLOWED = ("nexus_calibration", "nexus_neuronav", "nexus_pct", "nexus_ping")
IDLE_TIMEOUT_S = 1800          # MATLAB ends after this long without requests (it holds a licence)
HEARTBEAT_TIMEOUT_S = 60       # ... or when the Nexus has not touched the heartbeat for this long
HEARTBEAT_EVERY_S = 10
STARTUP_TIMEOUT_S = 180        # MATLAB must report ready within this (licence checks can be slow)
REQUEST_TIMEOUT_S = 600        # default time a request may take; MatlabBridge.call(timeout=...) overrides
STOP_GRACE_S = 10              # after a stop request, time MATLAB gets before it is ended
STALE_RESPONSE_S = 5           # an unreadable response this old is reported, not retried
TICK_BUSY_MS, TICK_IDLE_MS = 50, 500


class BridgeError(Exception):
    """A request that failed in MATLAB (or never reached it)."""

    def __init__(self, message: str, identifier: str = "", report: str = "", log: str = ""):
        super().__init__(message)
        self.message = message
        self.identifier = identifier
        self.report = report
        self.log = log


def _q(text: str | Path) -> str:
    return "'" + str(text).replace("'", "''") + "'"


def bridge_command(matlab: str, prestus_root: Path, bridge_dir: Path,
                   idle_s: int = IDLE_TIMEOUT_S, heartbeat_s: int = HEARTBEAT_TIMEOUT_S) -> list[str]:
    """matlab -batch command that starts nexus_bridge on bridge_dir."""
    fn = Path(prestus_root) / "functions"
    body = (f"addpath({_q(fn / 'helper')}, {_q(fn / 'gui')}); "
            f"nexus_bridge({_q(bridge_dir)}, IdleTimeout={idle_s}, HeartbeatTimeout={heartbeat_s})")
    return [matlab, "-batch", body]


def _write_atomic(path: Path, text: str) -> None:
    write_text_atomic(path, text)


class BridgeClient:
    """The bridge folder seen from the Nexus side."""

    _counter = itertools.count(1)

    def __init__(self, bridge_dir: str | Path):
        self.dir = Path(bridge_dir)
        for sub in ("req", "resp", "progress"):
            (self.dir / sub).mkdir(parents=True, exist_ok=True)

    def ready(self) -> dict | None:
        f = self.dir / "ready.json"
        try:
            return json.loads(f.read_text(encoding="utf-8")) if f.is_file() else None
        except (OSError, ValueError):
            return None

    def heartbeat(self) -> None:
        (self.dir / "heartbeat").touch()

    def stop(self) -> None:
        (self.dir / "stop").touch()

    def submit(self, fn: str, args: dict | None = None) -> str:
        if fn not in ALLOWED:
            raise BridgeError(f"{fn} cannot be called through the bridge.", "nexus_bridge:notAllowed")
        rid = f"{next(self._counter):06d}_{uuid.uuid4().hex[:8]}"
        _write_atomic(self.dir / "req" / f"{rid}.json",
                      json.dumps({"id": rid, "fn": fn, "args": args or {}}))
        return rid

    def poll(self) -> list[dict]:
        """Finished responses (each read once)."""
        out = []
        for f in sorted((self.dir / "resp").glob("*.json")):
            try:
                data = json.loads(f.read_text(encoding="utf-8"))
                if not isinstance(data, dict):
                    raise ValueError("not a JSON object")
            except (OSError, ValueError) as exc:
                try:                           # still being renamed into place, or broken
                    age = time.time() - f.stat().st_mtime
                except OSError:
                    continue
                if age < STALE_RESPONSE_S:
                    continue
                data = {"id": f.stem, "ok": False, "identifier": "nexus_bridge:badResponse",
                        "error": f"MATLAB's answer could not be read ({exc}).", "report": "", "log": ""}
            try:
                f.unlink(missing_ok=True)
            except OSError:
                pass
            out.append(data)
        return out

    def progress(self, rid: str) -> str | None:
        f = self.dir / "progress" / f"{rid}.txt"
        try:
            return f.read_text(encoding="utf-8") if f.is_file() else None
        except OSError:
            return None

    def call_blocking(self, fn: str, args: dict | None = None, timeout: float = 600) -> Any:
        """Submit and wait (scripts and tests); raises BridgeError on failure."""
        rid = self.submit(fn, args)
        t0 = time.monotonic()
        pending: dict[str, dict] = {}
        while time.monotonic() - t0 < timeout:
            self.heartbeat()
            for resp in self.poll():
                pending[resp.get("id", "")] = resp
            if rid in pending:
                return result_of(pending[rid])
            time.sleep(0.05)
        raise BridgeError(f"No answer from MATLAB within {timeout:g} s.", "nexus_bridge:timeout")


def result_of(resp: dict) -> Any:
    if resp.get("ok"):
        return resp.get("result")
    raise BridgeError(resp.get("error") or "MATLAB error", resp.get("identifier", ""),
                      resp.get("report", ""), resp.get("log", ""))


# ── JSON helpers: jsonencode writes a 1-element vector as a scalar and [] as [] ──
def as_list(v: Any) -> list:
    if v is None:
        return []
    return v if isinstance(v, list) else [v]


def as_floats(v: Any) -> list[float]:
    return [float("nan") if x is None else float(x) for x in as_list(v)]


def as_columns(v: Any) -> list[list[float]]:
    """A MATLAB matrix sent as num2cell(M, 1): one list per column."""
    return [as_floats(c) for c in as_list(v)]


# ── Qt side ──────────────────────────────────────────────────────────────
try:
    from PyQt6.QtCore import QObject, QProcess, QTimer, pyqtSignal
except ImportError:                      # pragma: no cover - the GUI always has PyQt6
    QObject = object                     # type: ignore[misc,assignment]


@dataclass
class _Pending:
    fn: str
    args: dict
    on_done: Callable[[Any], None] | None
    on_error: Callable[[BridgeError], None] | None
    on_progress: Callable[[str], None] | None
    timeout: float | None = None             # seconds once sent; None: the bridge's default
    rid: str = ""
    deadline: float = 0.0                     # time.monotonic() by which the answer is due
    last_progress: str = field(default="")


class _Heartbeat:
    """Touches the heartbeat file from a thread of its own, so a busy GUI thread (a long
    plot, a modal dialog) never lets MATLAB think the Nexus has gone."""

    def __init__(self, client: BridgeClient, every_s: float = HEARTBEAT_EVERY_S):
        self._client = client
        self._every = every_s
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, name="nexus-bridge-heartbeat", daemon=True)
        self._thread.start()

    def _run(self) -> None:
        while not self._stop.is_set():
            try:
                self._client.heartbeat()
            except OSError:
                pass
            self._stop.wait(self._every)

    def stop(self) -> None:
        self._stop.set()


def _end_in_background(pid: int, bridge_dir: Path | None, grace_s: float = STOP_GRACE_S,
                       pgid: int = 0) -> None:
    """After a stop request: give the session grace_s to end by itself, then end its
    process tree (and process group: MATLAB outlives its launcher); remove the bridge
    folder once it is gone."""
    from . import proctree

    def running() -> bool:
        return proctree.alive(pid) or proctree.group_alive(pgid)

    def run() -> None:
        t0 = time.monotonic()
        while running() and time.monotonic() - t0 < grace_s:
            time.sleep(0.2)
        if running():
            proctree.kill_tree(pid)
            proctree.kill_group(pgid)
        if bridge_dir is not None:
            time.sleep(0.5)
            shutil.rmtree(bridge_dir, ignore_errors=True)

    threading.Thread(target=run, name="nexus-bridge-stop", daemon=False).start()


class MatlabBridge(QObject):
    """One MATLAB session for the whole Nexus: started on the first call, kept while used.

    Deadlines: MATLAB must report ready within ``startup_timeout_s`` and answer each request
    within its timeout (``request_timeout_s`` unless ``call`` gives one); otherwise the
    session is stopped and the waiting requests fail with a message saying so."""

    state_changed = pyqtSignal(str)          # 'off' | 'starting' | 'ready' | 'busy' | 'failed'
    log = pyqtSignal(str)                    # MATLAB console output and request logs

    def __init__(self, prestus_root: Path, matlab: Callable[[], str | None],
                 command: Callable[[Path], list[str]] | None = None, parent=None,
                 startup_timeout_s: float = STARTUP_TIMEOUT_S,
                 request_timeout_s: float = REQUEST_TIMEOUT_S):
        super().__init__(parent)
        self.prestus_root = Path(prestus_root)
        self._matlab = matlab
        self._command = command
        self.startup_timeout_s = startup_timeout_s
        self.request_timeout_s = request_timeout_s
        self.proc: QProcess | None = None
        self.client: BridgeClient | None = None
        self.state = "off"
        self.output = ""                     # console output of the session (diagnostics)
        self.error = ""
        self._queue: list[_Pending] = []     # waiting for the session to be ready
        self._running: dict[str, _Pending] = {}
        self._timer = QTimer(self, interval=TICK_IDLE_MS)
        self._timer.timeout.connect(self._tick)
        self._heartbeat: _Heartbeat | None = None
        self._started_at = 0.0
        self._ending: tuple[str, str] | None = None     # (message, identifier) of a deliberate end
        self._culprit: dict[str, BridgeError] = {}      # rid -> its own error (timeout)

    # ── lifecycle ───────────────────────────────────────────────────────
    def _set_state(self, state: str) -> None:
        if state != self.state:
            self.state = state
            self.state_changed.emit(state)
        self._timer.setInterval(TICK_BUSY_MS if self.busy() or state == "starting" else TICK_IDLE_MS)

    def start(self) -> bool:
        if self.proc is not None:
            return True
        matlab = self._matlab() if self._command is None else "custom"
        if not matlab:
            self.error = "MATLAB not found; set its path on the Run tab."
            self._set_state("failed")
            return False
        bridge_dir = Path(tempfile.mkdtemp(prefix="prestus_nexus_bridge_"))
        self.client = BridgeClient(bridge_dir)
        self._heartbeat = _Heartbeat(self.client)
        argv = (self._command(bridge_dir) if self._command is not None
                else bridge_command(matlab, self.prestus_root, bridge_dir))
        self.output, self.error = "", ""
        self._ending = None
        self._culprit = {}
        self.proc = QProcess(self)
        from . import proctree
        proctree.new_session(self.proc)       # stop() can end MATLAB with its launcher
        self.proc.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self.proc.readyReadStandardOutput.connect(self._read)
        self.proc.finished.connect(self._finished)
        self.proc.errorOccurred.connect(self._proc_error)
        self._started_at = time.monotonic()
        self._set_state("starting")
        self.proc.start(argv[0], argv[1:])
        self._timer.start()
        return True

    def stop(self, wait_ms: int = 1000) -> None:
        """Ask the session to end, waiting at most wait_ms (the GUI is blocked meanwhile).
        A session that has not ended by then is ended in the background."""
        if self.client is not None:
            try:
                self.client.stop()
            except OSError:
                pass
        if self.proc is None:
            return
        self._ending = self._ending or ("The MATLAB session was stopped.", "nexus_bridge:stopped")
        if self.proc.waitForFinished(wait_ms) or self.proc is None:
            return
        from . import proctree
        pid = int(self.proc.processId() or 0)
        pgid = proctree.group_of(pid)
        bridge_dir = self.client.dir if self.client is not None else None
        self._teardown(-1, remove_dir=False)
        _end_in_background(pid, bridge_dir, pgid=pgid)

    def cancel(self, message: str = "Cancelled.") -> None:
        """Stop the session now and fail every waiting request with message."""
        if self.proc is None:
            failed = self._queue[:]
            self._queue.clear()
            for p in failed:
                if p.on_error is not None:
                    p.on_error(BridgeError(message, "nexus_bridge:cancelled"))
            return
        self._abort(message, "nexus_bridge:cancelled")

    def _abort(self, message: str, identifier: str) -> None:
        """End the session at once (stop file, then its process tree) with this error."""
        self._ending = (message, identifier)
        if self.client is not None:
            try:
                self.client.stop()
            except OSError:
                pass
        if self.proc is not None:
            pid = int(self.proc.processId() or 0)
            bridge_dir = self.client.dir if self.client is not None else None
            from . import proctree
            pgid = proctree.group_of(pid)
            proctree.kill_tree(pid)
            if self.proc is not None and not self.proc.waitForFinished(500):
                self._teardown(-1, remove_dir=False)
                _end_in_background(pid, bridge_dir, grace_s=0, pgid=pgid)

    def _read(self) -> None:
        if self.proc is None:
            return
        text = bytes(self.proc.readAllStandardOutput()).decode(errors="replace")
        self.output += text
        self.log.emit(text)

    def _proc_error(self, err) -> None:
        if err == QProcess.ProcessError.FailedToStart:
            self.error = f"Could not start MATLAB ({self.proc.program() if self.proc else ''})."
            self._finished(-1, None)

    def _finished(self, code: int, _status) -> None:
        if self.proc is None:
            return
        self._read()
        self._teardown(code)

    def _teardown(self, code: int, remove_dir: bool = True) -> None:
        """The session has ended (or is being left to end): fail what waits, reset."""
        proc, self.proc = self.proc, None
        was_ready = self.state in ("ready", "busy")
        tail = "\n".join(self.output.strip().splitlines()[-12:])
        identifier = "nexus_bridge:ended"
        if self._ending is not None:
            self.error, identifier = self._ending
        elif not self.error:
            self.error = ("The MATLAB session ended." if was_ready
                          else f"MATLAB could not start the bridge (exit code {code}).")
        if tail and self._ending is None:
            self.error += "\n\n" + tail
        if proc is not None:
            for sig in (proc.readyReadStandardOutput, proc.finished, proc.errorOccurred):
                try:
                    sig.disconnect()
                except (TypeError, RuntimeError):
                    pass
            proc.deleteLater()
        self._timer.stop()
        if self._heartbeat is not None:
            self._heartbeat.stop()
            self._heartbeat = None
        failed = list(self._running.values()) + self._queue
        self._running.clear()
        self._queue.clear()
        if self.client is not None:
            if remove_dir:
                shutil.rmtree(self.client.dir, ignore_errors=True)
            self.client = None
        deliberate = self._ending is not None and self._ending[1] == "nexus_bridge:stopped"
        self._set_state("off" if (was_ready or deliberate) and not failed else "failed")
        culprit, self._culprit = self._culprit, {}
        for p in failed:
            if p.on_error is not None:
                p.on_error(culprit.get(p.rid) or BridgeError(self.error, identifier))

    # ── requests ────────────────────────────────────────────────────────
    def call(self, fn: str, args: dict | None = None, on_done: Callable[[Any], None] | None = None,
             on_error: Callable[[BridgeError], None] | None = None,
             on_progress: Callable[[str], None] | None = None,
             timeout: float | None = None) -> None:
        """Run fn(args) in MATLAB; callbacks run on the GUI thread. timeout: seconds MATLAB
        has for this request once it is sent (default: request_timeout_s)."""
        p = _Pending(fn, args or {}, on_done, on_error, on_progress, timeout)
        if self.state in ("ready", "busy"):
            self._send(p)
            return
        self._queue.append(p)
        if self.proc is None and not self.start():
            self._queue.remove(p)
            if on_error is not None:
                on_error(BridgeError(self.error, "nexus_bridge:noMatlab"))

    def _send(self, p: _Pending) -> None:
        try:
            p.rid = self.client.submit(p.fn, p.args)
        except (BridgeError, OSError) as e:
            if p.on_error is not None:
                p.on_error(e if isinstance(e, BridgeError) else
                           BridgeError(f"Cannot send the request to MATLAB: {e}", "nexus_bridge:io"))
            return
        p.deadline = time.monotonic() + (self.request_timeout_s if p.timeout is None else p.timeout)
        self._running[p.rid] = p
        self._set_state("busy")

    def busy(self) -> bool:
        return bool(self._running or self._queue)

    def _tick(self) -> None:
        c = self.client
        if c is None:
            return
        now = time.monotonic()
        if self.state == "starting":
            if c.ready() is not None:
                self._set_state("ready")
                queued, self._queue = self._queue, []
                for p in queued:
                    self._send(p)
            elif now - self._started_at > self.startup_timeout_s:
                self._abort(f"MATLAB did not start the bridge within {self.startup_timeout_s:g} s; "
                            "the session was stopped. Check the MATLAB path and licence, then try "
                            "again.", "nexus_bridge:startTimeout")
                return
        for p in list(self._running.values()):
            if p.on_progress is not None:
                msg = c.progress(p.rid)
                if msg and msg != p.last_progress:
                    p.last_progress = msg
                    p.on_progress(msg)
        late = [p for p in list(self._running.values()) if now > p.deadline]
        if late:
            for p in late:
                limit = self.request_timeout_s if p.timeout is None else p.timeout
                self._culprit[p.rid] = BridgeError(
                    f"No answer from MATLAB to {p.fn} within {limit:g} s; the MATLAB session was "
                    "stopped (it starts again with the next request).", "nexus_bridge:timeout")
            self._abort("The MATLAB session was stopped because another request took too long.",
                        "nexus_bridge:ended")
            return
        for resp in c.poll():
            p = self._running.pop(resp.get("id", ""), None)
            if resp.get("log"):
                self.log.emit(resp["log"])
            if p is None:
                continue
            try:
                result = result_of(resp)
            except BridgeError as e:
                if p.on_error is not None:
                    p.on_error(e)
            else:
                if p.on_done is not None:
                    p.on_done(result)
            if self.client is not c:          # a callback stopped the session
                return
        if self.state == "busy" and not self._running:
            self._set_state("ready")


_SHARED: MatlabBridge | None = None


def shared_bridge(prestus_root: Path, matlab: Callable[[], str | None]) -> MatlabBridge:
    """The app-wide session: switching study windows keeps the same MATLAB."""
    global _SHARED
    if _SHARED is None:
        _SHARED = MatlabBridge(prestus_root, matlab)
    else:
        _SHARED._matlab = matlab
    return _SHARED


def shutdown_shared() -> None:
    """End the app-wide session (at quit): blocks the GUI for at most about a second."""
    global _SHARED
    if _SHARED is not None:
        _SHARED.stop()
        _SHARED = None
