"""HPC jobs submitted through PRESTUS: find their IDs, query and cancel them.

The Nexus does not write job scripts. "Submit to cluster" runs the study's
launcher with the configured platform, so PRESTUS submits as usual
(prestus_pipeline_start -> hpc_submit_job), and every submission prints
    Job "<name>" (ID: <id>) submitted successfully
which parse_submissions reads from the log. Chained pipelines (async,
multi-ISPPA, uncertainty) submit several jobs; each prints that line.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass

SUBMITTED = re.compile(r'Job "([^"]*)" \(ID: ([^)\s]+)\) submitted successfully')

# Scheduler states after which a job no longer changes
TERMINAL = {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL",
            "PREEMPTED", "BOOT_FAIL", "DEADLINE", "DONE", "FINISHED"}


@dataclass
class Job:
    job_id: str
    name: str = ""
    state: str = "SUBMITTED"

    @property
    def finished(self) -> bool:
        return self.state.split()[0].rstrip("+") in TERMINAL


def scheduler() -> str | None:
    """'slurm' or 'qsub' when its submit command is on PATH, else None."""
    if shutil.which("sbatch"):
        return "slurm"
    if shutil.which("qsub"):
        return "qsub"
    return None


def parse_submissions(text: str) -> list[Job]:
    """Jobs announced in a PRESTUS log, in order, without duplicates."""
    jobs: list[Job] = []
    seen: set[str] = set()
    for name, job_id in SUBMITTED.findall(text):
        if job_id not in seen:
            seen.add(job_id)
            jobs.append(Job(job_id=job_id, name=name))
    return jobs


def _run(argv: list[str], timeout: float = 15) -> str:
    try:
        out = subprocess.run(argv, capture_output=True, text=True, timeout=timeout, check=False)
    except (OSError, subprocess.TimeoutExpired):
        return ""
    return out.stdout


def parse_squeue(text: str) -> dict[str, str]:
    """'<id> <STATE>' lines (squeue -h -o '%i %T') -> {id: state}."""
    out: dict[str, str] = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 2:
            out[parts[0]] = parts[1]
    return out


def parse_sacct(text: str) -> dict[str, str]:
    """'<id>|<STATE>' lines (sacct -n -X -P -o JobID,State) -> {id: state}."""
    out: dict[str, str] = {}
    for line in text.splitlines():
        parts = line.strip().split("|")
        if len(parts) >= 2 and parts[0]:
            out[parts[0]] = parts[1].split()[0]
    return out


def parse_qstat(text: str) -> dict[str, str]:
    """qstat table rows -> {id: state}; state letters as qstat prints them (Q, R, C, ...)."""
    out: dict[str, str] = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 5 and re.match(r"^\d+", parts[0]):
            letter = parts[-2] if len(parts) >= 6 else parts[4]
            out[parts[0]] = {"Q": "PENDING", "R": "RUNNING", "C": "COMPLETED", "E": "RUNNING",
                             "H": "HELD", "F": "FINISHED"}.get(letter, letter)
    return out


def query(jobs: list[Job], kind: str | None = None) -> list[Job]:
    """Update the jobs' states from the scheduler (jobs it no longer lists keep their state)."""
    kind = kind or scheduler()
    open_jobs = [j for j in jobs if not j.finished]
    if not open_jobs or kind is None:
        return jobs
    ids = [j.job_id for j in open_jobs]
    states: dict[str, str] = {}
    if kind == "slurm":
        states.update(parse_squeue(_run(["squeue", "-h", "-j", ",".join(ids), "-o", "%i %T"])))
        missing = [i for i in ids if i not in states]
        if missing and shutil.which("sacct"):
            states.update(parse_sacct(_run(["sacct", "-n", "-X", "-P", "-j", ",".join(missing),
                                            "-o", "JobID,State"])))
    else:
        states.update(parse_qstat(_run(["qstat"] + ids)))
    for job in open_jobs:
        if job.job_id in states:
            job.state = states[job.job_id]
    return jobs


def cancel_command(job_ids: list[str], kind: str | None = None) -> list[str]:
    kind = kind or scheduler()
    return (["scancel"] if kind == "slurm" else ["qdel"]) + list(job_ids)


def cancel(job_ids: list[str], kind: str | None = None) -> str:
    """Cancel jobs; returns the scheduler's output."""
    if not job_ids:
        return ""
    return _run(cancel_command(job_ids, kind))
