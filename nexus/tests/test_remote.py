"""Submitting over SSH: path mapping, the ssh command and its quoting, jobs on the host."""

import os
import stat
import subprocess

import pytest

from prestus_nexus import jobs, launcher
from prestus_nexus import remote as rmt
from prestus_nexus.model import ConfigModel


def _wait_until(app, cond, timeout=10.0):
    import time
    t0 = time.monotonic()
    while not cond() and time.monotonic() - t0 < timeout:
        app.processEvents()
        time.sleep(0.01)
    return cond()


def _exe(path, body):
    path.write_text("#!/bin/sh\n" + body)
    path.chmod(path.stat().st_mode | stat.S_IEXEC)


def test_to_remote_uses_longest_prefix(tmp_path):
    (tmp_path / "proj" / "sub").mkdir(parents=True)
    r = rmt.RemoteHost("h", path_map=[(str(tmp_path), "/a"), (str(tmp_path / "proj"), "/project/x")])
    assert r.to_remote(tmp_path / "proj" / "sub" / "c.yaml") == "/project/x/sub/c.yaml"
    assert r.to_remote(tmp_path / "other") == "/a/other"
    assert r.to_remote(tmp_path / "proj") == "/project/x"
    with pytest.raises(rmt.PathMapError):
        r.to_remote("/elsewhere/c.yaml")
    # sibling folder sharing a name prefix is not under the mapping
    with pytest.raises(rmt.PathMapError):
        rmt.RemoteHost("h", path_map=[(str(tmp_path / "pro"), "/p")]).to_remote(tmp_path / "proj")


def test_remote_root(tmp_path):
    r = rmt.RemoteHost("h", path_map=[(str(tmp_path), "/project/x")])
    assert r.remote_root(tmp_path / "PRESTUS") == "/project/x/PRESTUS"
    r.prestus_root = "/home/me/PRESTUS/"
    assert r.remote_root(tmp_path / "PRESTUS") == "/home/me/PRESTUS"


def test_ssh_argv():
    r = rmt.RemoteHost("me@login", setup="module load slurm")
    argv = r.wrap(["squeue", "-o", "%i %T"])
    assert argv[:-2] == ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=15",
                         "-o", f"ControlPath={rmt.control_path()}"]
    assert rmt.host_of(argv) == "me@login"
    assert argv[-2] == "me@login"
    assert argv[-1] == "bash -lc 'module load slurm && squeue -o '\"'\"'%i %T'\"'\"''"


def test_shared_connection_commands():
    r = rmt.RemoteHost("me@login")
    cp = f"ControlPath={rmt.control_path()}"
    connect = rmt.connect_command(r)
    assert "ControlMaster=yes" in connect and cp in connect and "BatchMode=yes" not in connect
    assert connect[-2:] == ["me@login", "true"] and rmt.host_of(connect) == "me@login"
    assert rmt.check_command(r)[-3:] == ["-O", "check", "me@login"]
    assert rmt.disconnect_command(r)[-3:] == ["-O", "exit", "me@login"]
    assert rmt.host_of(rmt.disconnect_command(r)) == "me@login"


def test_askpass_wrapper_starts_the_nexus_askpass(tmp_path, monkeypatch):
    env = rmt.askpass_env(tmp_path)
    assert env["SSH_ASKPASS_REQUIRE"] == "force"
    script = tmp_path / "askpass.sh"
    assert env["SSH_ASKPASS"] == str(script) and os.access(script, os.X_OK)
    assert "-m prestus_nexus --askpass" in script.read_text()
    # main() hands --askpass to askpass.run, without starting the app
    from prestus_nexus import app, askpass
    seen = []
    monkeypatch.setattr(askpass, "run", lambda prompt: seen.append(prompt) or 0)
    assert app.main(["--askpass", "me@login's", "password:"]) == 0
    assert seen == ["me@login's password:"]


def test_submit_command_survives_both_shells(tmp_path, prestus_root):
    """Run the remote command string as sshd would (the login shell parses it once) with a
    fake matlab that prints its -batch argument: the MATLAB statement must arrive intact."""
    study_dir = tmp_path / "my study's"           # space and quote in the path
    study_dir.mkdir()
    study = study_dir / "config_a.yaml"
    study.write_text("{}\n")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    _exe(bin_dir / "matlab", 'pwd\nprintf "%s\\n" "$2"\n')
    r = rmt.RemoteHost("h", path_map=[(str(tmp_path), str(tmp_path))],
                       setup=f"export PATH={bin_dir}:$PATH", prestus_root="/opt/PRESTUS")
    argv = rmt.submit_command(r, study, prestus_root, [("path.sim", "/p/it's")])
    assert argv[0] == "ssh"
    out = subprocess.run(["sh", "-c", argv[-1]], capture_output=True, text=True, timeout=60).stdout
    lines = out.strip().splitlines()
    assert lines[-2] == str(study_dir.resolve())
    rdir = str(study_dir.resolve()).replace("'", "''")
    assert lines[-1] == launcher._wrapped(
        "addpath(genpath(fullfile('/opt/PRESTUS', 'functions'))); "
        f"parameters = load_parameters('config_a.yaml', '{rdir}', {{'path.sim', '/p/it''s'}}); "
        "prestus_pipeline_start(parameters)")


def test_parse_probe_and_scheduler_from_log():
    found = rmt.parse_probe("Welcome!\nhost: mentat005\nsbatch: /usr/bin/sbatch\nqsub: -\n"
                            "matlab: /opt/matlab\nprestus: - /x/PRESTUS\n")
    assert found == {"host": "mentat005", "sbatch": "/usr/bin/sbatch", "qsub": "-",
                     "matlab": "/opt/matlab", "prestus": "- /x/PRESTUS"}
    assert rmt.scheduler_from_log("x\n➤ auto-detected: SLURM\n") == "slurm"
    assert rmt.scheduler_from_log("➤ deploying: QSUB\n") == "qsub"
    assert rmt.scheduler_from_log("➤ deploying: MATLAB\n") is None


def test_map_value():
    r = rmt.RemoteHost("h", path_map=[("/Volumes/2425122.01", "/project/2425122.01")])
    assert r.map_value("/Volumes/2425122.01/data/sub-%03d/T1.nii") == "/project/2425122.01/data/sub-%03d/T1.nii"
    assert r.map_value("/Volumes/2425122.01") == "/project/2425122.01"
    assert r.map_value("/Volumes/2425122.010/x") is None
    assert r.map_value("data/simnibs") is None


def test_map_value_windows_paths():
    r = rmt.RemoteHost("h", path_map=[("Z:\\2425122.01", "/project/2425122.01"),
                                      ("\\\\fileserver\\project", "/project")])
    assert r.map_value("Z:\\2425122.01\\data\\sub-%03d\\T1.nii") == \
        "/project/2425122.01/data/sub-%03d/T1.nii"
    assert r.map_value("z:\\2425122.01\\data") == "/project/2425122.01/data"      # drive letter case
    assert r.map_value("Z:/2425122.01/Data") == "/project/2425122.01/Data"          # forward slashes
    assert r.map_value("Z:\\2425122.01") == "/project/2425122.01"
    assert r.map_value("Z:\\2425122.010\\x") is None
    assert r.map_value("Y:\\2425122.01\\x") is None
    assert r.map_value("\\\\fileserver\\project\\2425122.02\\sim") == "/project/2425122.02/sim"
    assert r.map_value("\\\\FILESERVER\\Project\\x") == "/project/x"            # UNC, any case
    assert r.map_value("//fileserver/project/x") == "/project/x"
    assert r.map_value("\\\\otherserver\\project\\x") is None
    lower = rmt.RemoteHost("h", path_map=[("z:/2425122.01", "/project/2425122.01")])
    assert lower.map_value("Z:\\2425122.01\\a") == "/project/2425122.01/a"
    posix = rmt.RemoteHost("h", path_map=[("/Volumes/P", "/project/p")])
    assert posix.map_value("/volumes/P/x") is None                                # POSIX stays case-sensitive


def test_path_overrides(schema):
    m = ConfigModel(schema)
    m.set("path.anat", "/Volumes/2425122.01/v05/data/simnibs")
    m.set("path.sim", "/Volumes/2425122.01/v05/data/tussim/demo")
    m.set("path.seg", "/Users/me/seg")                  # local, outside the mapping
    m.set("path.data", "/project/2425122.01/raw")       # already a cluster path
    r = rmt.RemoteHost("h", path_map=[("/Volumes/2425122.01", "/project/2425122.01")])
    mapped, unmapped = rmt.path_overrides(m, r)
    assert ("path.anat", "/project/2425122.01/v05/data/simnibs") in mapped
    assert ("path.sim", "/project/2425122.01/v05/data/tussim/demo") in mapped
    assert all(k != "path.data" for k, _ in mapped)
    assert unmapped == ["path.seg"]


def test_path_overrides_keys_transducers(schema):
    m = ConfigModel(schema)
    m.set("transducer.library_path", "/Volumes/2425122.01/lib")
    r = rmt.RemoteHost("h", path_map=[("/Volumes/2425122.01", "/project/2425122.01")])
    assert ("transducer.library_path", "/project/2425122.01/lib") in rmt.path_overrides(m, r)[0]
    m.add_transducer()
    keys = [k for k, _ in rmt.path_overrides(m, r)[0]]
    assert "transducer(1).library_path" in keys and "transducer(2).library_path" in keys


def test_query_and_cancel_over_ssh(monkeypatch):
    calls = []

    def fake_run(argv, timeout=15):
        calls.append(argv)
        if "squeue" in argv[-1]:
            return "1|RUNNING|0:10|4:00:00|3:59:50|gpu|node1|/w\n"
        return "" if "sstat" in argv[-1] else "2|COMPLETED|00:05:00|04:00:00|0:0||n|gpu|/w\n"
    monkeypatch.setattr(jobs, "_run", fake_run)
    monkeypatch.setattr(jobs.shutil, "which", lambda name: None)   # nothing local
    r = rmt.RemoteHost("me@login")
    tracked = [jobs.Job("1"), jobs.Job("2")]
    jobs.query(tracked, "slurm", r)
    assert [j.state for j in tracked] == ["RUNNING", "COMPLETED"]
    assert all(c[0] == "ssh" for c in calls) and "sacct" in calls[-1][-1]
    jobs.cancel(["1"], "slurm", r)
    assert calls[-1][-1] == "bash -lc 'scancel 1'"


# ── Run page, with a fake ssh that runs the remote command here ─────────

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def test_run_page_submits_over_ssh(schema, prestus_root, tmp_path, monkeypatch):
    pytest.importorskip("PyQt6")
    from PyQt6.QtCore import QEventLoop, QTimer
    from PyQt6.QtWidgets import QApplication, QMessageBox

    from prestus_nexus.session import Session
    from prestus_nexus.ui.main_window import RUN, MainWindow
    app = QApplication.instance() or QApplication([])  # noqa: F841

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    _exe(bin_dir / "ssh", 'for a; do last=$a; done\nexec sh -c "$last"\n')   # as sshd: run the last argument
    _exe(bin_dir / "matlab", "echo '➤ auto-detected: SLURM'\n"
                             "echo 'Job \"PS_sub-001\" (ID: 4242) submitted successfully'\n"
                             f"echo '{launcher.OK_SENTINEL}'\n")
    _exe(bin_dir / "squeue", "echo '4242|PENDING|0:00|8:00:00|8:00:00|gpu40g|(Priority)|/w'\n")
    monkeypatch.setenv("PATH", f"{bin_dir}:{os.environ['PATH']}")
    monkeypatch.setattr(QMessageBox, "question",
                        lambda *a, **k: QMessageBox.StandardButton.Yes)   # Issues-panel errors

    proj = tmp_path / "proj"
    proj.mkdir()
    m = ConfigModel(schema)
    m.save(proj / "config_r.yaml")
    window = MainWindow(m, Session(tmp_path / "s.yaml"), prestus_root)
    page = window._page(RUN)
    page.remote_toggle.setChecked(True)
    page.r_host.setText("me@login")
    page.r_local.setText(str(proj))
    page.r_remote.setText(str(proj))
    page.r_setup.setText(f"export PATH={bin_dir}:$PATH")
    page.r_root.setText(str(prestus_root))
    page.refresh()
    assert page.b_submit.isEnabled()

    results = []
    loop = QEventLoop()
    page.finished.connect(lambda task, ok: (results.append((task, ok)), loop.quit()))
    QTimer.singleShot(30000, loop.quit)
    page.submit()
    loop.exec()
    assert results == [("Submit", True)], page.log.toPlainText()
    _wait_until(app, lambda: not page.jobs_busy())         # the scheduler query runs in the background
    assert [j.job_id for j in page.jobs] == ["4242"], page.log.toPlainText()
    assert page.job_kind == "slurm" and page.job_remote is not None
    assert page.jobs[0].state == "PENDING"            # squeue ran through the fake ssh
    assert window.session.get("remote_host") == "me@login"
    page.jobs_timer.stop()


def test_run_page_connects_with_askpass(schema, prestus_root, tmp_path, monkeypatch):
    """Connect… runs the ControlMaster login with SSH_ASKPASS pointing at the Nexus helper."""
    pytest.importorskip("PyQt6")
    from PyQt6.QtCore import QEventLoop, QTimer
    from PyQt6.QtWidgets import QApplication

    from prestus_nexus.session import Session
    from prestus_nexus.ui.main_window import RUN, MainWindow
    app = QApplication.instance() or QApplication([])  # noqa: F841

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    _exe(bin_dir / "ssh", 'echo "askpass=$SSH_ASKPASS require=$SSH_ASKPASS_REQUIRE args=$*"\n'
                          'test -x "$SSH_ASKPASS"\n')
    monkeypatch.setenv("PATH", f"{bin_dir}:{os.environ['PATH']}")
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "cfg"))
    m = ConfigModel(schema)
    m.save(tmp_path / "config_c.yaml")
    window = MainWindow(m, Session(tmp_path / "s.yaml"), prestus_root)
    page = window._page(RUN)
    page.remote_toggle.setChecked(True)
    page.r_host.setText("me@login")

    results = []
    loop = QEventLoop()
    page.finished.connect(lambda task, ok: (results.append((task, ok)), loop.quit()))
    QTimer.singleShot(20000, loop.quit)
    page.connect_ssh()
    loop.exec()
    log = page.log.toPlainText()
    assert results == [("Connect", True)], log
    assert "require=force" in log and "ControlMaster=yes" in log and "askpass.sh" in log
    assert page.r_status.text() == "Connected to me@login"


def test_run_page_shows_selected_job_log(schema, prestus_root, tmp_path, monkeypatch):
    """Clicking a job shows its output log, read through the folder mapping, and the job row
    carries timing, place and memory."""
    pytest.importorskip("PyQt6")
    from PyQt6.QtWidgets import QApplication

    from prestus_nexus.session import Session
    from prestus_nexus.ui.main_window import RUN, MainWindow
    app = QApplication.instance() or QApplication([])  # noqa: F841

    logs = tmp_path / "proj" / "sim" / "sub-001" / "log_hpc"
    logs.mkdir(parents=True)
    (logs / "sub-001_slurm_output_4242.log").write_text("Running k-Wave simulation...\n")
    m = ConfigModel(schema)
    m.save(tmp_path / "config_v.yaml")
    window = MainWindow(m, Session(tmp_path / "s.yaml"), prestus_root)
    app = QApplication.instance()
    page = window._page(RUN)
    page.job_kind = "slurm"
    page.job_remote = rmt.RemoteHost("h", path_map=[(str(tmp_path / "proj"), "/project/x")])
    monkeypatch.setattr(jobs, "_run", lambda argv, timeout=15: "")     # no scheduler calls
    page.jobs = [jobs.Job("4242", "PS_sub-001", state="RUNNING", elapsed="3:50:00",
                          time_limit="4:00:00", time_left="10:00", partition="gpu40g",
                          where="dccn-c078", max_rss="12.0 GB",
                          work_dir="/project/x/sim/sub-001/log_hpc")]
    jv = page.jobs_view
    page.jobs_changed.emit()
    item = jv.tree.topLevelItem(0)
    assert [item.text(c) for c in range(3, 7)] == ["3:50:00 / 4:00:00", "10:00",
                                                    "gpu40g: dccn-c078", "12.0 GB"]
    from prestus_nexus.ui import style
    assert item.foreground(4).color().name() == style.active()["red"]      # time running out
    item.setSelected(True)
    _wait_until(app, lambda: jv.views["output"].toPlainText() != "")    # read on a worker thread
    assert jv.views["output"].toPlainText() == "Running k-Wave simulation...\n"
    assert "sub-001_slurm_output_4242.log" in jv.path_label.text()
    assert jv.timer.isActive()                               # following a running job
    jv.timer.stop()
    # the Jobs side tab holds this view; Pop out opens a second one on the same job
    assert window.run_side.indexOf(jv) >= 0
    page.pop_out_jobs()
    pv = page.jobs_window.view
    assert pv.selected_id() == "4242"
    _wait_until(app, lambda: pv.views["output"].toPlainText() != "")
    assert pv.views["output"].toPlainText() == "Running k-Wave simulation...\n"
    pv.timer.stop()
    page.jobs_window.close()


