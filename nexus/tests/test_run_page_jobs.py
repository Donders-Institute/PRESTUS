"""Run page: scheduler queries run off the GUI thread, one at a time; Stop ends the
process tree."""

import os
import sys
import threading
import time

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")

from prestus_nexus import jobs, proctree  # noqa: E402
from prestus_nexus.model import ConfigModel  # noqa: E402


def _wait(app, cond, timeout=10.0):
    t0 = time.monotonic()
    while not cond() and time.monotonic() - t0 < timeout:
        app.processEvents()
        time.sleep(0.01)
    return cond()


@pytest.fixture
def page(schema, prestus_root, tmp_path):
    from PyQt6.QtWidgets import QApplication

    from prestus_nexus.session import Session
    from prestus_nexus.ui.main_window import RUN, MainWindow
    app = QApplication.instance() or QApplication([])
    m = ConfigModel(schema)
    m.save(tmp_path / "config_j.yaml")
    window = MainWindow(m, Session(tmp_path / "s.yaml"), prestus_root)
    p = window._page(RUN)
    p.app = app
    yield p
    p.jobs_timer.stop()


def test_poll_runs_in_background_one_at_a_time(page, monkeypatch):
    release = threading.Event()
    calls = []

    def slow_query(js, kind=None, remote=None):
        calls.append(threading.current_thread().name)
        release.wait(5)
        for j in js:
            j.state = "RUNNING"
        return js
    monkeypatch.setattr(jobs, "query", slow_query)
    page.job_kind = "slurm"
    page.jobs = [jobs.Job("11", "PS")]
    changed = []
    page.jobs_changed.connect(lambda: changed.append(1))
    t0 = time.monotonic()
    page.refresh_jobs()
    page.refresh_jobs()                         # in flight: not started again
    assert time.monotonic() - t0 < 1 and page.jobs_busy()
    assert page.jobs[0].state == "SUBMITTED"    # the shown job is not touched by the worker
    release.set()
    assert _wait(page.app, lambda: not page.jobs_busy())
    assert len(calls) == 1 and calls[0] != threading.current_thread().name
    assert page.jobs[0].state == "RUNNING" and changed
    assert page.jobs_timer.isActive()


def test_failed_query_keeps_last_state(page, monkeypatch):
    def broken(js, kind=None, remote=None):
        raise RuntimeError("ssh died")
    monkeypatch.setattr(jobs, "query", broken)
    page.job_kind = "slurm"
    page.jobs = [jobs.Job("12", state="PENDING")]
    page.refresh_jobs()
    assert _wait(page.app, lambda: not page.jobs_busy())
    assert page.jobs[0].state == "PENDING" and "ssh died" in page.log.toPlainText()


@pytest.mark.skipif(os.name == "nt", reason="POSIX shell")
def test_stop_kills_the_child_of_a_launcher(page, tmp_path):
    pidfile = tmp_path / "child.pid"
    script = (f"{sys.executable} -c 'import time; time.sleep(60)' & echo $! > {pidfile}; wait")
    page._start("Run", ["/bin/sh", "-c", script], tmp_path)
    assert _wait(page.app, lambda: pidfile.is_file() and pidfile.read_text().strip())
    child = int(pidfile.read_text())
    assert proctree.alive(child)
    page.stop()
    assert _wait(page.app, lambda: page.proc is None)
    t0 = time.monotonic()
    while proctree.alive(child) and time.monotonic() - t0 < 5:
        time.sleep(0.05)
    assert not proctree.alive(child)
