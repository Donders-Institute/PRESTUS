"""HPC jobs submitted through PRESTUS: find their IDs, query and cancel them, read their logs.

The Nexus does not write job scripts. "Submit to cluster" runs the study's
launcher with the configured platform, so PRESTUS submits as usual
(prestus_pipeline_start -> hpc_submit_job), and every submission prints
    Job "<name>" (ID: <id>) submitted successfully
which parse_submissions reads from the log. Chained pipelines (async,
multi-ISPPA, uncertainty) submit several jobs; each prints that line.

With a remote host (remote.RemoteHost), the scheduler commands run over SSH.

SLURM jobs also report their timing (elapsed, limit, time left), where they
run (node, or why they are pending), peak memory (sstat while running, sacct
after) and exit code. Their logs are in the job's working directory, PRESTUS'
<output>/log_hpc: <subj>_slurm_output_<id>.log, <subj>_slurm_error_<id>.log and,
once the report job has run, <subj>_job<id>_report.txt (hpc_job_report).
"""

from __future__ import annotations

import re
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

SUBMITTED = re.compile(r'Job "([^"]*)" \(ID: ([^)\s]+)\) submitted successfully')

# Scheduler states after which a job no longer changes
TERMINAL = {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL",
            "PREEMPTED", "BOOT_FAIL", "DEADLINE", "DONE", "FINISHED"}

# Job log files in the job's working directory (see hpc_submit_job, hpc_job_report)
LOG_PATTERNS = {"output": "*_slurm_output_{id}.log",
                "error": "*_slurm_error_{id}.log",
                "report": "*_job{id}_report.txt"}
TAIL_BYTES = 200_000


@dataclass
class Job:
    job_id: str
    name: str = ""
    state: str = "SUBMITTED"
    elapsed: str = ""          # [D-]HH:MM:SS as SLURM prints it
    time_limit: str = ""
    time_left: str = ""
    partition: str = ""
    where: str = ""            # node list while running, (Reason) while pending
    max_rss: str = ""          # peak memory, e.g. '12.3 GB'
    exit_code: str = ""
    work_dir: str = ""         # on the scheduler's side (remote path for SSH submissions)

    @property
    def finished(self) -> bool:
        return self.state.split()[0].rstrip("+") in TERMINAL

    @property
    def running(self) -> bool:
        return self.state.startswith("RUNNING")

    def time_warning(self) -> bool:
        """Running with little time left: under 10 % of the limit, but at least 15 min
        and at most 1 h (so a 24 h job warns at 1 h, a 1 h job at 15 min)."""
        left, limit = slurm_seconds(self.time_left), slurm_seconds(self.time_limit)
        if not self.running or left is None:
            return False
        return left < min(max(15 * 60, 0.1 * (limit or 0)), 3600)


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


def _on(argv: list[str], remote=None) -> str:
    """Run argv here, or on the remote host over SSH (slower, so a longer timeout)."""
    if remote is None:
        return _run(argv)
    return _run(remote.wrap(argv), timeout=45)


# ── SLURM output ────────────────────────────────────────────────────────

def slurm_seconds(text: str) -> int | None:
    """[D-]HH:MM:SS, MM:SS or MM -> seconds; None for UNLIMITED, INVALID, empty."""
    m = re.fullmatch(r"(?:(\d+)-)?(\d+)(?::(\d+))?(?::(\d+))?", (text or "").strip())
    if not m:
        return None
    days = int(m.group(1) or 0)
    parts = [int(p) for p in m.group(2, 3, 4) if p is not None]
    if len(parts) == 3:
        h, mi, s = parts
    elif len(parts) == 2:
        h, mi, s = (0, *parts) if not m.group(1) else (parts[0], parts[1], 0)
    else:
        h, mi, s = (parts[0], 0, 0) if m.group(1) else (0, parts[0], 0)
    return ((days * 24 + h) * 60 + mi) * 60 + s


def format_rss(text: str) -> str:
    """SLURM memory ('2048K', '1.5G'; a bare number is KiB) -> '2.0 MB' style; '' if empty."""
    value = _rss_bytes(text or "")
    if value is None:
        return ""
    for unit in ("B", "KB", "MB", "GB"):
        if value < 1024:
            return f"{value:.0f} B" if unit == "B" else f"{value:.1f} {unit}"
        value /= 1024
    return f"{value:.1f} TB"


SQUEUE_FORMAT = "%i|%T|%M|%l|%L|%P|%R|%Z"


def parse_squeue(text: str) -> dict[str, dict]:
    """squeue -h -o SQUEUE_FORMAT lines -> {id: fields}."""
    out: dict[str, dict] = {}
    for line in text.splitlines():
        parts = line.strip().split("|")
        if len(parts) >= 2 and parts[0]:
            parts += [""] * (8 - len(parts))
            out[parts[0]] = {"state": parts[1], "elapsed": parts[2], "time_limit": parts[3],
                             "time_left": parts[4], "partition": parts[5], "where": parts[6],
                             "work_dir": parts[7]}
    return out


SACCT_FIELDS = "JobID,State,Elapsed,Timelimit,ExitCode,MaxRSS,NodeList,Partition,WorkDir"


def parse_sacct(text: str) -> dict[str, dict]:
    """sacct -n -P -o SACCT_FIELDS lines (job and step rows) -> {id: fields}.

    The job row gives state and timing; peak memory is the largest MaxRSS of its steps.
    """
    out: dict[str, dict] = {}
    rss: dict[str, float] = {}
    for line in text.splitlines():
        parts = line.strip().split("|")
        if len(parts) < 2 or not parts[0]:
            continue
        parts += [""] * (9 - len(parts))
        job, _, step = parts[0].partition(".")
        if not step:
            out[job] = {"state": parts[1].split()[0], "elapsed": parts[2], "time_limit": parts[3],
                        "exit_code": parts[4], "where": parts[6], "partition": parts[7],
                        "work_dir": parts[8]}
        if parts[5]:
            b = _rss_bytes(parts[5])
            if b is not None and b >= rss.get(job, -1):
                rss[job] = b
                out.setdefault(job, {})["max_rss_raw"] = parts[5]
    for job, fields in out.items():
        if "max_rss_raw" in fields:
            fields["max_rss"] = format_rss(fields.pop("max_rss_raw"))
    return out


def parse_sstat(text: str) -> dict[str, str]:
    """sstat -n -P -o JobID,MaxRSS lines -> {job id: peak memory so far}."""
    best: dict[str, tuple[float, str]] = {}
    for line in text.splitlines():
        parts = line.strip().split("|")
        if len(parts) >= 2 and parts[0] and parts[1]:
            job = parts[0].split(".")[0]
            b = _rss_bytes(parts[1])
            if b is not None and b >= best.get(job, (-1, ""))[0]:
                best[job] = (b, parts[1])
    return {job: format_rss(raw) for job, (_, raw) in best.items()}


def _rss_bytes(text: str) -> float | None:
    m = re.fullmatch(r"([\d.]+)([KMGTP]?)", text.strip())
    if not m:
        return None
    return float(m.group(1)) * 1024 ** ("KMGTP".find(m.group(2) or "K") + 1)


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


def _apply(job: Job, fields: dict) -> None:
    for key, value in fields.items():
        if value or key == "where":
            setattr(job, key, value)


def query(jobs: list[Job], kind: str | None = None, remote=None) -> list[Job]:
    """Update the jobs from the scheduler (jobs it no longer lists keep their last values).

    SLURM: squeue for queued/running jobs (plus sstat for their peak memory so far),
    sacct for jobs squeue no longer lists.
    """
    kind = kind or (None if remote else scheduler())
    open_jobs = [j for j in jobs if not j.finished]
    if not open_jobs or kind is None:
        return jobs
    ids = [j.job_id for j in open_jobs]
    by_id = {j.job_id: j for j in open_jobs}
    if kind == "slurm":
        queued = parse_squeue(_on(["squeue", "-h", "-j", ",".join(ids), "-o", SQUEUE_FORMAT], remote))
        for job_id, fields in queued.items():
            if job_id in by_id:
                _apply(by_id[job_id], fields)
        running = [i for i, f in queued.items() if f["state"] == "RUNNING"]
        if running and (remote is not None or shutil.which("sstat")):
            peaks = parse_sstat(_on(["sstat", "-n", "-P", "-a", "-j", ",".join(running),
                                     "-o", "JobID,MaxRSS"], remote))
            for job_id, peak in peaks.items():
                if job_id in by_id and peak:
                    by_id[job_id].max_rss = peak
        missing = [i for i in ids if i not in queued]
        if missing and (remote is not None or shutil.which("sacct")):
            done = parse_sacct(_on(["sacct", "-n", "-P", "-j", ",".join(missing),
                                    "-o", SACCT_FIELDS], remote))
            for job_id, fields in done.items():
                if job_id in by_id:
                    _apply(by_id[job_id], fields)
                    by_id[job_id].time_left = ""
    else:
        states = parse_qstat(_on(["qstat"] + ids, remote))
        for job in open_jobs:
            if job.job_id in states:
                job.state = states[job.job_id]
    return jobs


def cancel_command(job_ids: list[str], kind: str | None = None) -> list[str]:
    kind = kind or scheduler()
    return (["scancel"] if kind == "slurm" else ["qdel"]) + list(job_ids)


def cancel(job_ids: list[str], kind: str | None = None, remote=None) -> str:
    """Cancel jobs; returns the scheduler's output."""
    if not job_ids:
        return ""
    return _on(cancel_command(job_ids, kind), remote)


# ── job logs ────────────────────────────────────────────────────────────

def _tail_file(path: Path, nbytes: int = TAIL_BYTES) -> str:
    with open(path, "rb") as f:
        f.seek(0, 2)
        size = f.tell()
        f.seek(max(0, size - nbytes))
        data = f.read()
    text = data.decode(errors="replace")
    return text if size <= nbytes else "…\n" + text.split("\n", 1)[-1]


def read_log(job: Job, which: str, remote=None) -> tuple[str, str]:
    """(path, last TAIL_BYTES of the text) of one of the job's logs (LOG_PATTERNS key).

    Read from the file system when the job's folder is visible here (on the cluster,
    or the project drive mounted locally via the SSH folder mapping), else over SSH.
    ('', '') when the folder is unknown or the file does not exist (yet).
    """
    if not job.work_dir:
        return "", ""
    pattern = LOG_PATTERNS[which].format(id=job.job_id)
    local_dir = job.work_dir if remote is None else remote.to_local_or_none(job.work_dir)
    if local_dir and Path(local_dir).is_dir():
        hits = sorted(Path(local_dir).glob(pattern))
        if not hits:
            return "", ""
        try:
            return str(hits[-1]), _tail_file(hits[-1])
        except OSError:
            return str(hits[-1]), ""
    if remote is None:
        return "", ""
    script = (f"cd {shlex.quote(job.work_dir)} && f=$(ls -1 {pattern} 2>/dev/null | tail -n 1) "
              f"&& [ -n \"$f\" ] && echo \"$PWD/$f\" && tail -c {TAIL_BYTES} \"$f\"")
    out = _run(remote.ssh(script), timeout=45)
    path, _, text = out.partition("\n")
    return path.strip(), text
