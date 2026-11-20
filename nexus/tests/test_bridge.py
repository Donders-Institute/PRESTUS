"""MATLAB bridge: the file protocol, the Qt client (against a fake server), and -- with
MATLAB -- the real nexus_bridge serving nexus_ping and nexus_calibration."""

import glob
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from prestus_nexus import bridge as br  # noqa: E402

FAKE = Path(__file__).with_name("fake_bridge.py")


def test_bridge_command_starts_nexus_bridge(tmp_path):
    argv = br.bridge_command("/m/matlab", Path("/opt/PRESTUS"), tmp_path / "b")
    assert argv[:2] == ["/m/matlab", "-batch"]
    fn = Path("/opt/PRESTUS") / "functions"          # the platform's separators
    assert f"addpath('{fn / 'helper'}', '{fn / 'gui'}')" in argv[2]
    assert f"nexus_bridge('{tmp_path / 'b'}', IdleTimeout=" in argv[2]


def test_json_helpers_undo_matlab_scalar_collapse():
    assert br.as_list(None) == [] and br.as_list(3) == [3] and br.as_list([1, 2]) == [1, 2]
    assert br.as_floats([1, None]) == [1.0, pytest.approx(float("nan"), nan_ok=True)]
    assert br.as_columns([[1, 2], 3]) == [[1.0, 2.0], [3.0]]


def test_client_refuses_functions_outside_the_allow_list(tmp_path):
    c = br.BridgeClient(tmp_path)
    with pytest.raises(br.BridgeError):
        c.submit("system", {"cmd": "rm -rf /"})
    assert not list((tmp_path / "req").iterdir())


def test_client_round_trip_with_fake_server(tmp_path):
    proc = subprocess.Popen([sys.executable, str(FAKE), str(tmp_path)])
    try:
        c = br.BridgeClient(tmp_path)
        t0 = time.monotonic()
        while c.ready() is None and time.monotonic() - t0 < 10:
            time.sleep(0.02)
        assert c.ready()["matlab"] == "fake"
        assert c.call_blocking("nexus_ping", {"x": 1}, timeout=10)["echo"] == {"x": 1}
        with pytest.raises(br.BridgeError) as err:
            c.call_blocking("nexus_ping", {"action": "fail"}, timeout=10)
        assert err.value.identifier == "fake:fail" and err.value.log == "printed"
    finally:
        br.BridgeClient(tmp_path).stop()
        proc.wait(10)


@pytest.fixture(scope="module")
def app():
    pytest.importorskip("PyQt6")
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _wait(app, cond, timeout=10.0):
    t0 = time.monotonic()
    while not cond() and time.monotonic() - t0 < timeout:
        app.processEvents()
        time.sleep(0.01)
    return cond()


def test_qt_bridge_queues_until_ready_then_calls_back(app, tmp_path):
    b = br.MatlabBridge(tmp_path, lambda: "unused",
                        command=lambda d: [sys.executable, str(FAKE), str(d)])
    states, done, errors, progress = [], [], [], []
    b.state_changed.connect(states.append)
    try:
        b.call("nexus_ping", {"a": 1}, on_done=done.append, on_error=errors.append)
        b.call("nexus_ping", {"progress": "Depth 1/2"}, on_done=done.append, on_progress=progress.append)
        b.call("nexus_ping", {"action": "fail"}, on_done=done.append, on_error=errors.append)
        assert _wait(app, lambda: len(done) == 2 and len(errors) == 1)
        assert done[0]["echo"] == {"a": 1} and errors[0].message == "it failed"
        assert progress == ["Depth 1/2"]
        assert states[:2] == ["starting", "ready"] and b.state == "ready"
        bridge_dir = b.client.dir
    finally:
        b.stop()
    assert _wait(app, lambda: b.proc is None) and b.state == "off"
    assert not bridge_dir.exists()


def test_qt_bridge_reports_a_failed_start(app, tmp_path):
    b = br.MatlabBridge(tmp_path, lambda: "unused",
                        command=lambda d: [sys.executable, str(FAKE), str(d), "--fail-start"])
    errors = []
    b.call("nexus_ping", {}, on_error=errors.append)
    assert _wait(app, lambda: bool(errors))
    assert b.state == "failed" and "no licence" in errors[0].message


def test_qt_bridge_without_matlab_fails_at_once(app, tmp_path):
    b = br.MatlabBridge(tmp_path, lambda: None)
    errors = []
    b.call("nexus_ping", {}, on_error=errors.append)
    assert errors and "MATLAB not found" in errors[0].message and b.state == "failed"


def test_qt_bridge_startup_deadline(app, tmp_path):
    b = br.MatlabBridge(tmp_path, lambda: "unused", startup_timeout_s=0.5,
                        command=lambda d: [sys.executable, str(FAKE), str(d), "--never-ready"])
    errors = []
    b.call("nexus_ping", {}, on_error=errors.append)
    assert _wait(app, lambda: bool(errors))
    assert errors[0].identifier == "nexus_bridge:startTimeout" and "within 0.5 s" in errors[0].message
    assert b.proc is None and b.state == "failed"


def test_qt_bridge_request_deadline_stops_the_session(app, tmp_path):
    b = br.MatlabBridge(tmp_path, lambda: "unused", request_timeout_s=30,
                        command=lambda d: [sys.executable, str(FAKE), str(d)])
    done, errors = [], []
    b.call("nexus_ping", {"sleep": 20}, on_done=done.append, on_error=errors.append, timeout=0.5)
    assert _wait(app, lambda: bool(errors))
    assert not done and errors[0].identifier == "nexus_bridge:timeout"
    assert "nexus_ping within 0.5 s" in errors[0].message
    assert b.proc is None and b.client is None
    # the next call starts a fresh session
    b.call("nexus_ping", {"a": 2}, on_done=done.append, on_error=errors.append)
    try:
        assert _wait(app, lambda: bool(done)) and done[0]["echo"] == {"a": 2}
    finally:
        b.stop()


def test_qt_bridge_cancel_fails_pending_requests(app, tmp_path):
    b = br.MatlabBridge(tmp_path, lambda: "unused",
                        command=lambda d: [sys.executable, str(FAKE), str(d)])
    errors = []
    b.call("nexus_ping", {"sleep": 20}, on_error=errors.append)
    assert _wait(app, lambda: b.state == "busy")
    t0 = time.monotonic()
    b.cancel()
    assert time.monotonic() - t0 < 2
    assert errors and errors[0].identifier == "nexus_bridge:cancelled"
    assert b.proc is None and not b.busy()


def test_qt_bridge_stop_returns_quickly_when_matlab_hangs(app, tmp_path):
    b = br.MatlabBridge(tmp_path, lambda: "unused",
                        command=lambda d: [sys.executable, str(FAKE), str(d), "--never-ready"])
    b.start()
    assert _wait(app, lambda: b.proc is not None and b.proc.processId() > 0)
    pid = b.proc.processId()
    t0 = time.monotonic()
    b.stop(wait_ms=300)
    assert time.monotonic() - t0 < 1.5 and b.proc is None
    from prestus_nexus import proctree
    t1 = time.monotonic()
    while proctree.alive(pid) and time.monotonic() - t1 < br.STOP_GRACE_S + 5:
        time.sleep(0.2)
    assert not proctree.alive(pid)


def test_heartbeat_runs_without_the_event_loop(tmp_path):
    c = br.BridgeClient(tmp_path)
    hb = br._Heartbeat(c, every_s=0.05)
    try:
        time.sleep(0.2)                      # no processEvents: a blocked GUI thread
        assert (tmp_path / "heartbeat").is_file()
        first = (tmp_path / "heartbeat").stat().st_mtime_ns
        time.sleep(0.3)
        assert (tmp_path / "heartbeat").stat().st_mtime_ns >= first
    finally:
        hb.stop()


def test_poll_reports_unreadable_old_responses(tmp_path):
    c = br.BridgeClient(tmp_path)
    bad = tmp_path / "resp" / "000001_abc.json"
    bad.write_text("{not json", encoding="utf-8")
    assert c.poll() == [] and bad.exists()          # may still be arriving
    old = time.time() - br.STALE_RESPONSE_S - 1
    os.utime(bad, (old, old))
    out = c.poll()
    assert len(out) == 1 and out[0]["id"] == "000001_abc" and out[0]["ok"] is False
    assert out[0]["identifier"] == "nexus_bridge:badResponse" and not bad.exists()
    with pytest.raises(br.BridgeError):
        br.result_of(out[0])


# ── with MATLAB ──────────────────────────────────────────────────────────
def _find_matlab():
    cand = [os.environ.get("PRESTUS_MATLAB"), shutil.which("matlab")]
    cand += sorted(glob.glob("/Applications/MATLAB_R20*.app/bin/matlab"), reverse=True)
    return next((c for c in cand if c and Path(c).exists()), None)


MATLAB = _find_matlab()


@pytest.fixture(scope="module")
def matlab_bridge(prestus_root, tmp_path_factory):
    if MATLAB is None:
        pytest.skip("MATLAB not available")
    d = tmp_path_factory.mktemp("bridge")
    c = br.BridgeClient(d)
    c.heartbeat()
    proc = subprocess.Popen(br.bridge_command(MATLAB, prestus_root, d),
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    t0 = time.monotonic()
    while c.ready() is None and proc.poll() is None and time.monotonic() - t0 < 240:
        c.heartbeat()
        time.sleep(0.2)
    if c.ready() is None:
        proc.kill()
        pytest.fail("nexus_bridge did not start:\n" + (proc.stdout.read() if proc.stdout else ""))
    yield c
    c.stop()
    proc.wait(60)


@pytest.mark.matlab
def test_matlab_bridge_ping_and_errors(matlab_bridge):
    out = matlab_bridge.call_blocking("nexus_ping", {"x": [1, 2]}, timeout=60)
    assert out["echo"]["x"] == [1, 2] and out["prestus_root"]
    with pytest.raises(br.BridgeError) as err:
        matlab_bridge.call_blocking("nexus_calibration", {"action": "view", "session": "none"}, timeout=60)
    assert err.value.identifier == "nexus_calibration:noSession"


@pytest.mark.matlab
def test_matlab_bridge_calibration_predict(matlab_bridge, prestus_root):
    out = matlab_bridge.call_blocking("nexus_calibration", {
        "action": "open", "session": "t", "serial": "CTX_250_001", "driving_system": "SC_105_010",
        "library_path": ""}, timeout=120)
    assert out["combo"] == "CTX_250_001_SC_105_010" and out["calibration_type"]
    assert br.as_floats(out["default_depths"])
    view = matlab_bridge.call_blocking("nexus_calibration", {
        "action": "view", "session": "t", "predicted": True, "depths": [30, 45],
        "settings": {"n_starts": 2}}, timeout=120)
    res = view["results"]
    assert {"geometric", "library"} <= set(res)
    assert len(br.as_floats(res["library"]["focus"])) == 2
    assert len(br.as_columns(res["library"]["Imod"])) == 2
    assert view["messages"][0].startswith("Predicted 2 depths")
