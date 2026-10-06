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
import os
import shutil
import tempfile
import time
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

ALLOWED = ("nexus_calibration", "nexus_neuronav", "nexus_ping")
IDLE_TIMEOUT_S = 1800          # MATLAB ends after this long without requests (it holds a licence)
HEARTBEAT_TIMEOUT_S = 60       # ... or when the Nexus has not touched the heartbeat for this long
HEARTBEAT_EVERY_S = 10


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
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


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
            except (OSError, ValueError):
                continue                       # still being renamed into place
            f.unlink(missing_ok=True)
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
    rid: str = ""
    last_progress: str = field(default="")


class MatlabBridge(QObject):
    """One MATLAB session for the whole Nexus: started on the first call, kept while used."""

    state_changed = pyqtSignal(str)          # 'off' | 'starting' | 'ready' | 'busy' | 'failed'
    log = pyqtSignal(str)                    # MATLAB console output and request logs

    def __init__(self, prestus_root: Path, matlab: Callable[[], str | None],
                 command: Callable[[Path], list[str]] | None = None, parent=None):
        super().__init__(parent)
        self.prestus_root = Path(prestus_root)
        self._matlab = matlab
        self._command = command
        self.proc: QProcess | None = None
        self.client: BridgeClient | None = None
        self.state = "off"
        self.output = ""                     # console output of the session (diagnostics)
        self.error = ""
        self._queue: list[_Pending] = []     # waiting for the session to be ready
        self._running: dict[str, _Pending] = {}
        self._timer = QTimer(self, interval=50)
        self._timer.timeout.connect(self._tick)
        self._last_heartbeat = 0.0

    # ── lifecycle ───────────────────────────────────────────────────────
    def _set_state(self, state: str) -> None:
        if state != self.state:
            self.state = state
            self.state_changed.emit(state)

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
        self.client.heartbeat()
        argv = (self._command(bridge_dir) if self._command is not None
                else bridge_command(matlab, self.prestus_root, bridge_dir))
        self.output, self.error = "", ""
        self.proc = QProcess(self)
        self.proc.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self.proc.readyReadStandardOutput.connect(self._read)
        self.proc.finished.connect(self._finished)
        self.proc.errorOccurred.connect(self._proc_error)
        self._set_state("starting")
        self.proc.start(argv[0], argv[1:])
        self._timer.start()
        return True

    def stop(self) -> None:
        """Ask the session to end (it also ends by itself when idle or orphaned)."""
        if self.client is not None:
            self.client.stop()
        if self.proc is not None and not self.proc.waitForFinished(5000):
            self.proc.kill()
            self.proc.waitForFinished(2000)

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
        was_ready = self.state in ("ready", "busy")
        tail = "\n".join(self.output.strip().splitlines()[-12:])
        if not self.error:
            self.error = ("The MATLAB session ended." if was_ready
                          else f"MATLAB could not start the bridge (exit code {code}).")
        if tail:
            self.error += "\n\n" + tail
        self.proc.deleteLater()
        self.proc = None
        self._timer.stop()
        failed = list(self._running.values()) + self._queue
        self._running.clear()
        self._queue.clear()
        if self.client is not None:
            shutil.rmtree(self.client.dir, ignore_errors=True)
            self.client = None
        self._set_state("off" if was_ready and not failed else "failed")
        for p in failed:
            if p.on_error is not None:
                p.on_error(BridgeError(self.error, "nexus_bridge:ended"))

    # ── requests ────────────────────────────────────────────────────────
    def call(self, fn: str, args: dict | None = None, on_done: Callable[[Any], None] | None = None,
             on_error: Callable[[BridgeError], None] | None = None,
             on_progress: Callable[[str], None] | None = None) -> None:
        """Run fn(args) in MATLAB; callbacks run on the GUI thread."""
        p = _Pending(fn, args or {}, on_done, on_error, on_progress)
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
        except BridgeError as e:
            if p.on_error is not None:
                p.on_error(e)
            return
        self._running[p.rid] = p
        self._set_state("busy")

    def busy(self) -> bool:
        return bool(self._running or self._queue)

    def _tick(self) -> None:
        c = self.client
        if c is None:
            return
        now = time.monotonic()
        if now - self._last_heartbeat > HEARTBEAT_EVERY_S:
            c.heartbeat()
            self._last_heartbeat = now
        if self.state == "starting" and c.ready() is not None:
            self._set_state("ready")
            queued, self._queue = self._queue, []
            for p in queued:
                self._send(p)
        for p in self._running.values():
            if p.on_progress is not None:
                msg = c.progress(p.rid)
                if msg and msg != p.last_progress:
                    p.last_progress = msg
                    p.on_progress(msg)
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
    global _SHARED
    if _SHARED is not None:
        _SHARED.stop()
        _SHARED = None
