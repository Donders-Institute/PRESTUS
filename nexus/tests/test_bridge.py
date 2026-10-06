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
    assert "addpath('/opt/PRESTUS/functions/helper', '/opt/PRESTUS/functions/gui')" in argv[2]
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
