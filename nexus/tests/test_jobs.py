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


def test_parse_squeue_and_sacct():
    assert jobs.parse_squeue("4815162 RUNNING\n4815163 PENDING\n") == {
        "4815162": "RUNNING", "4815163": "PENDING"}
    assert jobs.parse_sacct("4815162|COMPLETED\n4815163|CANCELLED by 123\n") == {
        "4815162": "COMPLETED", "4815163": "CANCELLED"}


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
            return "1 RUNNING\n"
        return "2|COMPLETED\n"
    monkeypatch.setattr(jobs, "_run", fake_run)
    monkeypatch.setattr(jobs.shutil, "which", lambda name: "/usr/bin/" + name)
    tracked = [jobs.Job("1"), jobs.Job("2"), jobs.Job("3", state="FAILED")]
    jobs.query(tracked, "slurm")
    assert [j.state for j in tracked] == ["RUNNING", "COMPLETED", "FAILED"]
    assert calls == ["squeue", "sacct"]


def test_cancel_command():
    assert jobs.cancel_command(["1", "2"], "slurm") == ["scancel", "1", "2"]
    assert jobs.cancel_command(["77.torque01"], "qsub") == ["qdel", "77.torque01"]
