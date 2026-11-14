"""HPC job tracking: IDs from PRESTUS logs, scheduler output parsing, cancel commands."""

from prestus_nexus import jobs

LOG = """Running ...
SLURM command: sbatch --parsable /x/temp_slurm.sh (in /x)
Job "PS_sub-001" (ID: 4815162) submitted successfully
Monitor: squeue --user $USER ...
Job "PS_sub-001_rep" (ID: 4815163) submitted successfully
Job "PS_sub-001" (ID: 4815162) submitted successfully
qsub: Job "PS_sub-002" (ID: 77.torque01) submitted successfully
"""


def test_parse_submissions_dedupes_and_keeps_order():
    found = jobs.parse_submissions(LOG)
    assert [j.job_id for j in found] == ["4815162", "4815163", "77.torque01"]
    assert found[0].name == "PS_sub-001"
    assert all(j.state == "SUBMITTED" and not j.finished for j in found)


def test_parse_squeue():
    text = ("4815162|RUNNING|1:02:03|4:00:00|2:57:57|gpu40g|dccn-c078|/project/x/sub-001/log_hpc\n"
            "4815163|PENDING|0:00|8:00:00|8:00:00|gpu40g|(Priority)|/project/x/sub-002/log_hpc\n")
    out = jobs.parse_squeue(text)
    assert out["4815162"] == {"state": "RUNNING", "elapsed": "1:02:03", "time_limit": "4:00:00",
                              "time_left": "2:57:57", "partition": "gpu40g", "where": "dccn-c078",
                              "work_dir": "/project/x/sub-001/log_hpc"}
    assert out["4815163"]["where"] == "(Priority)"


def test_parse_sacct_takes_peak_memory_from_steps():
    text = ("4815162|COMPLETED|01:10:00|04:00:00|0:0||dccn-c078|gpu40g|/project/x/log_hpc\n"
            "4815162.batch|COMPLETED|01:10:00||0:0|31457280K|dccn-c078||\n"
            "4815162.extern|COMPLETED|01:10:00||0:0|1024K|dccn-c078||\n"
            "4815163|CANCELLED by 123|00:00:00|08:00:00|0:0||None assigned|gpu40g|/p\n")
    out = jobs.parse_sacct(text)
    assert out["4815162"]["state"] == "COMPLETED" and out["4815162"]["max_rss"] == "30.0 GB"
    assert out["4815162"]["exit_code"] == "0:0" and out["4815162"]["work_dir"] == "/project/x/log_hpc"
    assert out["4815163"]["state"] == "CANCELLED" and "max_rss" not in out["4815163"]


def test_parse_sstat_and_units():
    assert jobs.parse_sstat("4815162.batch|2097152K\n4815162.extern|512K\n") == {"4815162": "2.0 GB"}
    assert jobs.format_rss("1.5G") == "1.5 GB" and jobs.format_rss("") == ""
    assert jobs.slurm_seconds("1-02:03:04") == 93784 and jobs.slurm_seconds("05:00") == 300
    assert jobs.slurm_seconds("UNLIMITED") is None


def test_time_warning():
    assert jobs.Job("1", state="RUNNING", time_left="10:00", time_limit="4:00:00").time_warning()
    assert jobs.Job("1", state="RUNNING", time_left="50:00", time_limit="1-00:00:00").time_warning()
    assert not jobs.Job("1", state="RUNNING", time_left="2:00:00", time_limit="1-00:00:00").time_warning()
    assert jobs.Job("1", state="RUNNING", time_left="40:00", time_limit="8:00:00").time_warning()  # < 10 %
    assert not jobs.Job("1", state="PENDING", time_left="5:00", time_limit="4:00:00").time_warning()


def test_parse_qstat():
    text = ("Job ID          Name       User   Time Use S Queue\n"
            "--------------- ---------- ------ -------- - -----\n"
            "77.torque01     PS_sub-002 me     00:01:02 R batch\n")
    assert jobs.parse_qstat(text) == {"77.torque01": "RUNNING"}


def test_finished_states():
    assert jobs.Job("1", state="COMPLETED").finished
    assert jobs.Job("1", state="CANCELLED").finished
    assert not jobs.Job("1", state="PENDING").finished


def test_query_updates_open_jobs_only(monkeypatch):
    calls = []

    def fake_run(argv, timeout=15):
        calls.append(argv[0])
        if argv[0] == "squeue":
            return "1|RUNNING|0:10|4:00:00|3:59:50|gpu|node1|/w\n"
        if argv[0] == "sstat":
            return "1.batch|1048576K\n"
        return "2|COMPLETED|00:05:00|04:00:00|0:0||node2|gpu|/w\n2.batch|COMPLETED|||0:0|2048K|||\n"
    monkeypatch.setattr(jobs, "_run", fake_run)
    monkeypatch.setattr(jobs.shutil, "which", lambda name: "/usr/bin/" + name)
    tracked = [jobs.Job("1"), jobs.Job("2"), jobs.Job("3", state="FAILED")]
    jobs.query(tracked, "slurm")
    assert [j.state for j in tracked] == ["RUNNING", "COMPLETED", "FAILED"]
    assert calls == ["squeue", "sstat", "sacct"]
    assert tracked[0].max_rss == "1.0 GB" and tracked[0].where == "node1" and tracked[0].work_dir == "/w"
    assert tracked[1].max_rss == "2.0 MB" and tracked[1].exit_code == "0:0"


def test_read_log_from_folder(tmp_path):
    (tmp_path / "sub-001_slurm_output_77.log").write_text("line 1\nline 2\n")
    job = jobs.Job("77", work_dir=str(tmp_path))
    path, text = jobs.read_log(job, "output")
    assert path.endswith("sub-001_slurm_output_77.log") and text == "line 1\nline 2\n"
    assert jobs.read_log(job, "report") == ("", "")
    assert jobs.read_log(jobs.Job("77"), "output") == ("", "")


def test_read_log_through_folder_mapping(tmp_path):
    from prestus_nexus.remote import RemoteHost
    (tmp_path / "log_hpc").mkdir()
    (tmp_path / "log_hpc" / "sub-001_job77_report.txt").write_text("Wall time: 1 h\n")
    r = RemoteHost("h", path_map=[(str(tmp_path), "/project/x")])
    path, text = jobs.read_log(jobs.Job("77", work_dir="/project/x/log_hpc"), "report", r)
    assert text == "Wall time: 1 h\n" and path.startswith(str(tmp_path))


def test_cancel_command():
    assert jobs.cancel_command(["1", "2"], "slurm") == ["scancel", "1", "2"]
    assert jobs.cancel_command(["77.torque01"], "qsub") == ["qdel", "77.torque01"]


def test_run_replaces_undecodable_output():
    import sys
    out = jobs._run([sys.executable, "-c",
                     "import sys; sys.stdout.buffer.write(b'ok \\xff\\xfe end\\n')"])
    assert out.startswith("ok ") and out.rstrip().endswith("end") and "�" in out


def test_run_timeout_and_missing_program_return_empty():
    import sys
    assert jobs._run([sys.executable, "-c", "import time; time.sleep(5)"], timeout=0.3) == ""
    assert jobs._run(["prestus-no-such-program-xyz"]) == ""


def test_run_swallows_unexpected_errors(monkeypatch):
    def boom(*a, **k):
        raise ValueError("bad argv")
    monkeypatch.setattr(jobs.subprocess, "run", boom)
    assert jobs._run(["squeue"]) == ""
