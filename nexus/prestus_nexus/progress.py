"""Progress of a PRESTUS run, read from its log: the MATLAB output of a local run or the
SLURM output log of a cluster job.

What the log provides, and what is read from it:

- stages: log_timer prints ``▶ <label> [#n]`` when a stage starts and ``⏱ <label> …``
  when it stops (prestus_pipeline: segmentation, placement, preproc, medium, source,
  acoustic, thermal, …).
- solver start: ``[progress] <stage> 0/1 @<epoch> solver <name>`` (prestus_progress,
  printed by acoustic_simulation, jwave_wrapper, babelbrain_wrapper, thermal_simulation).
- steps: ``[progress] <stage> k/n @<epoch> <detail>`` (the thermal segment loop; detail
  carries the simulated time and the peak temperature so far).
- k-Wave C++ binaries (cpp_cpu, cpp_gpu, cpp_mlx): their table rows
  ``|   45% |   12.345s |   15.000s | <finish> |`` (progress, elapsed, time to go).
- k-Wave MATLAB solvers (matlab_cpu/gpu, axisymmetric, elastic): ``precomputation completed
  in <t>`` and, after 50 time steps, ``estimated simulation time <t>...``. They print no
  further progress, so the fraction is estimated from the clock.
- jwave time-domain solvers: ``[run_jwave_sim]   step n/Nt`` and ``estimated remaining time``.
  The Helmholtz solver prints no progress: only the elapsed time is known.
- BabelBrain: functions/babelbrain/bb_progress.py samples its solvers' time loops every 5 s
  and prints ``[progress]`` lines (thermal detail ``Tnow=<°C>``: the current maximum).
- live view (simulation.live_view): ``[live] <stage> <png path>``.

The parser keeps no clock of its own: estimate(now) turns what was read into a fraction
and the time left at a given time.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

STAGE_START = re.compile(r"^▶ (\S+) \[#\d+\]")
STAGE_STOP = re.compile(r"^⏱ (\S+)\s")
PROGRESS = re.compile(r"^\[progress\] (\S+) (\d+)/(\d+) @([\d.]+)\s*(.*)$")
LIVE = re.compile(r"^\[live\] (\S+) (.+?)\s*$")
CPP_ROW = re.compile(r"^\|\s*(\d+(?:\.\d+)?)%\s*\|\s*([\d.]+)s\s*\|\s*([\d.]+)s\s*\|")
PRECOMP = re.compile(r"^\s*precomputation completed in (.+?)\s*$")
MATLAB_EST = re.compile(r"^\s*estimated simulation time (.+?)\.\.\.\s*$")
JWAVE_STEP = re.compile(r"^\[run_jwave_sim\]\s+step\s+(\d+)/(\d+)")
JWAVE_ETA = re.compile(r"^\[run_jwave_sim\]\s+estimated remaining time: ([\d.]+) s")
CPP_PHASES = {"Reading simulation configuration": "reading the configuration",
              "Reading input file": "reading the input file", "Data loading": "loading the input",
              "FFT plans creation": "creating FFT plans", "Pre-processing phase": "pre-processing",
              "Sampled data post-processing": "post-processing", "Creating output file": "writing the output",
              "Reading output file": "reading the output", "Creating checkpoint": "writing a checkpoint"}
CPP_PHASE = re.compile(r"^\|?\s*\+?\s*(" + "|".join(CPP_PHASES) + r")\b")
SUBMITTED = re.compile(r'Job "([^"]*)" \(ID: ([^)\s]+)\) submitted successfully')   # as jobs.SUBMITTED
TEMP = re.compile(r"\b(Tmax|Tnow)=([-\d.]+)C")
SIM_TIME = re.compile(r"\bt=([\d.]+)s")
# Progress over a series of runs (not a solver): drive-check depths, study rows, variants
OUTER_STAGES = {"calibration": "Depth", "batch": "Run"}
DONE_OK = ("[PRESTUS Nexus] OK", "Pipeline finished successfully")
DONE_FAIL = "[PRESTUS Nexus] FAILED"

# Lines worth keeping from a job log (remote grep -E and local filtering)
RELEVANT = re.compile(
    r"^(▶ |⏱ |\[progress\] |\[live\] |\|\s*\d+(\.\d+)?%|\s*precomputation completed in"
    r"|\s*estimated simulation time|\[run_jwave_sim\]\s+(step|estimated remaining)"
    r"|Pipeline finished successfully|\[PRESTUS Nexus\]|Job \".*\" \(ID: .*submitted successfully"
    r"|\|?\s*\+?\s*(" + "|".join(CPP_PHASES) + r"))")
RELEVANT_ERE = (r"^(▶ |⏱ |\[progress\] |\[live\] |\|[[:space:]]*[0-9.]+%"
                r"|[[:space:]]*precomputation completed in|[[:space:]]*estimated simulation time"
                r"|\[run_jwave_sim\][[:space:]]+(step|estimated remaining)"
                r"|Pipeline finished successfully|\[PRESTUS Nexus\]|Job \".*\" \(ID: .*submitted successfully"
                r"|\|?[[:space:]]*\+?[[:space:]]*(" + "|".join(CPP_PHASES) + r"))")

STAGE_NAMES = {
    "prestus_pipeline": "Pipeline", "segmentation": "Segmentation", "placement": "Placement",
    "preproc": "Head model", "medium": "Medium", "nifti_medium": "Medium maps", "source": "Source setup",
    "freefield_baseline": "Free-water baseline", "time_reversal": "Time reversal",
    "acoustic": "Acoustic simulation", "acoustic_analysis": "Acoustic analysis",
    "nifti_acoustic": "Acoustic maps", "thermal": "Thermal simulation",
    "thermal_analysis": "Thermal analysis", "nifti_thermal": "Thermal maps", "report": "Report",
}
SOLVER_NAMES = {
    "matlab_cpu": "k-Wave MATLAB, CPU", "matlab_gpu": "k-Wave MATLAB, GPU",
    "cpp_cpu": "k-Wave C++, CPU", "cpp_gpu": "k-Wave C++, CUDA", "cpp_mlx": "k-Wave MLX",
    "kwave_elastic": "k-Wave elastic", "jwave_helmholtz": "jwave Helmholtz",
    "jwave_pseudo_spectral": "jwave pseudo-spectral", "jwave_finite_differences": "jwave finite differences",
    "babelbrain": "BabelBrain", "kWaveDiffusion": "k-Wave diffusion",
}
# Solvers that print nothing between start and end: elapsed time only
SILENT_SOLVERS = {"jwave_helmholtz"}


def scale_time_seconds(text: str) -> float | None:
    """k-Wave's scaleTime text ('1min 23.4s', '2hours 3min 4s', '1 days, 2 hours, and 3 min',
    '12.5s') -> seconds; None when it is not a time."""
    units = {"year": 365 * 86400, "week": 7 * 86400, "day": 86400, "hour": 3600,
             "min": 60, "s": 1}
    total, found = 0.0, False
    for value, unit in re.findall(r"([\d.]+)\s*(years?|weeks?|days?|hours?|min|s)\b", text):
        try:
            total += float(value) * units[unit.rstrip("s") if unit != "s" else "s"]
        except (KeyError, ValueError):
            return None
        found = True
    return total if found else None


@dataclass
class Solver:
    stage: str                       # 'acoustic' | 'thermal'
    name: str                        # code_type / backend ('cpp_gpu', 'jwave_helmholtz', …)
    start: float                     # epoch [s] at the solver call
    owner: str = ""                  # innermost open stage when it started
    k: int = 0                       # [progress] steps
    n: int = 0
    t_first: float | None = None     # epoch of the first step line after the start
    k_first: int = 0
    t_last: float | None = None
    detail: str = ""
    pct: float | None = None         # C++ table or jwave steps [0-100]
    pct_seen: float | None = None    # when that row was read (None: unknown)
    time_to_go: float | None = None  # [s] at pct_seen
    precomp_s: float | None = None   # k-Wave MATLAB precomputation [s]
    loop_est_s: float | None = None  # k-Wave MATLAB estimated time loop [s]
    phase: str = ""                  # k-Wave C++ phase (CPP_PHASES), before and after the time loop


@dataclass
class Estimate:
    fraction: float | None           # 0-1; None: unknown (busy indicator)
    eta_s: float | None              # time left [s]; None: unknown
    estimated: bool                  # fraction from the clock, not from the solver
    text: str                        # one line for the user (title, percentage, and detail)
    title: str = ""                  # the stage ('Thermal simulation')
    detail: str = ""                 # solver, time left, and solver detail, ' · '-joined
    outer: str = ""                  # progress over a series of runs ('Depth 2 of 5 · 45.0 mm')

    @property
    def percent(self) -> str:
        """'45 %', '~63 %' (estimated), or '' when unknown."""
        if self.fraction is None:
            return ""
        return f"{'~' if self.estimated else ''}{100 * self.fraction:.0f} %"


@dataclass
class RunProgress:
    open_stages: list[str] = field(default_factory=list)    # started, not stopped (oldest first)
    done_stages: list[str] = field(default_factory=list)
    solver: Solver | None = None
    live: dict[str, str] = field(default_factory=dict)       # stage -> frame path (.png or .json)
    # peak temperature over the thermal run: (x, °C, x is simulated time [s] (else % done), kind)
    temperature: list[tuple[float, float, bool, str]] = field(default_factory=list)
    outer: tuple[str, int, int, str] | None = None            # (stage, k, n, detail): OUTER_STAGES
    children: list[tuple[str, str]] = field(default_factory=list)   # (name, id) submitted by the run
    finished: bool = False
    failed: bool = False
    _partial: str = ""

    @classmethod
    def from_text(cls, text: str, now: float | None = None) -> "RunProgress":
        p = cls()
        p.feed(text, now)
        p.flush(now)
        return p

    @property
    def stage(self) -> str:
        """The current stage: the innermost open one below the pipeline itself."""
        inner = [s for s in self.open_stages if s != "prestus_pipeline"]
        return inner[-1] if inner else ""

    def has_news(self) -> bool:
        return bool(self.open_stages or self.done_stages or self.solver or self.finished or self.failed)

    def feed(self, text: str, now: float | None = None) -> None:
        """Read more log text (it may end in the middle of a line). now: when it arrived
        (epoch s), for live logs; None for a log read later (cluster jobs)."""
        data = self._partial + text
        lines = data.split("\n")
        self._partial = lines.pop()
        for line in lines:
            self._line(line.rstrip("\r"), now)

    def flush(self, now: float | None = None) -> None:
        if self._partial:
            self._line(self._partial, now)
            self._partial = ""

    def _line(self, line: str, now: float | None) -> None:
        if not line:
            return
        m = STAGE_START.match(line)
        if m:
            label = m.group(1)
            if label == "prestus_pipeline":            # a new (nested) pipeline run
                self.finished = self.failed = False
            # a solver that ran outside this new stage is over (e.g. a baseline between stages)
            if self.solver is not None and self.solver.owner not in self.open_stages:
                self.solver = None
            if label in self.open_stages:
                self.open_stages.remove(label)
            self.open_stages.append(label)
            return
        m = STAGE_STOP.match(line)
        if m:
            label = m.group(1)
            if label in self.open_stages:
                self.open_stages.remove(label)
            if label not in self.done_stages and label != "prestus_pipeline":
                self.done_stages.append(label)
            if self.solver is not None and label == self.solver.owner:
                self.solver = None
            return
        m = PROGRESS.match(line)
        if m:
            stage, k, n, t, detail = m.group(1), int(m.group(2)), int(m.group(3)), float(m.group(4)), m.group(5)
            if stage in OUTER_STAGES:
                self.outer = (stage, k, n, detail)
                return
            if k == 0:
                if stage == "thermal":
                    self.temperature = []
                name = detail.split(" ", 1)[1] if detail.startswith("solver ") else detail
                self.solver = Solver(stage, name.strip(), t, owner=self.stage)
                return
            s = self.solver
            if s is None or s.stage != stage:
                s = self.solver = Solver(stage, "", t, owner=self.stage)
            if s.t_first is None or k < s.k:
                s.t_first, s.k_first = t, k
            s.k, s.n, s.t_last, s.detail = k, n, t, detail
            mt = TEMP.search(detail)
            if stage == "thermal" and mt and n > 0:
                ts = SIM_TIME.search(detail)
                x = float(ts.group(1)) if ts else 100.0 * k / n
                self.temperature.append((x, float(mt.group(2)), bool(ts), mt.group(1)))
            return
        m = LIVE.match(line)
        if m:
            self.live[m.group(1)] = m.group(2)
            return
        if line.startswith(DONE_OK):
            self.finished = True
            return
        if line.startswith(DONE_FAIL):
            self.failed = True
            return
        m = SUBMITTED.search(line)
        if m:
            if (m.group(1), m.group(2)) not in self.children:
                self.children.append((m.group(1), m.group(2)))
            return
        s = self.solver
        if s is None:
            return
        m = CPP_ROW.match(line)
        if m:
            s.pct, s.time_to_go, s.pct_seen = float(m.group(1)), float(m.group(3)), now
            s.phase = ""
            return
        m = CPP_PHASE.match(line)
        if m:
            s.phase = CPP_PHASES[m.group(1)]
            return
        m = JWAVE_STEP.match(line)
        if m and int(m.group(2)) > 0:
            s.pct, s.pct_seen = 100.0 * int(m.group(1)) / int(m.group(2)), now
            return
        m = JWAVE_ETA.match(line)
        if m:
            s.time_to_go, s.pct_seen = float(m.group(1)), now
            return
        m = PRECOMP.match(line)
        if m:
            s.precomp_s = scale_time_seconds(m.group(1))
            s.loop_est_s = None
            return
        m = MATLAB_EST.match(line)
        if m:
            s.loop_est_s = scale_time_seconds(m.group(1))

    # ── estimate ──────────────────────────────────────────────────────────
    def overall(self, expected: dict[str, float] | None, now: float) -> tuple[float, float] | None:
        """(fraction, seconds left) of the whole run from predicted seconds per stage label
        (workload.expected_seconds): finished stages count in full, the current stage by its
        solver's fraction (and time left), later stages by their prediction. None without a
        prediction."""
        if not expected or sum(expected.values()) <= 0:
            return None
        if self.finished and not self.open_stages:
            return 1.0, 0.0
        total = sum(expected.values())
        done = sum(v for k, v in expected.items() if k in self.done_stages)
        cur = self.solver.owner if self.solver is not None and self.solver.owner in expected else self.stage
        left_cur = 0.0
        if cur in expected and cur not in self.done_stages:
            e = self._estimate(now)
            f = e.fraction if self.solver is not None and e.fraction is not None else 0.0
            done += f * expected[cur]
            left_cur = e.eta_s if self.solver is not None and e.eta_s is not None else (1 - f) * expected[cur]
        later = sum(v for k, v in expected.items() if k not in self.done_stages and k != cur)
        return min(0.99, done / total), left_cur + later

    def estimate(self, now: float) -> Estimate:
        e = self._estimate(now)
        if self.outer is not None:
            stage, k, n, detail = self.outer
            e.outer = f"{OUTER_STAGES[stage]} {k} of {n}" + (f" · {detail}" if detail else "")
        return e

    def _estimate(self, now: float) -> Estimate:
        if self.failed:
            return Estimate(None, None, False, "Failed: see the log", "Failed", "see the log")
        if self.finished and not self.open_stages:
            return Estimate(1.0, 0.0, False, "Finished", "Finished")
        stage = STAGE_NAMES.get(self.stage, self.stage.replace("_", " ").capitalize()) if self.stage else ""
        s = self.solver
        if s is None:
            text = stage or ("Starting" if not self.done_stages else "Working")
            return Estimate(None, None, False, text, text)
        frac, eta, est, note = self._solver_fraction(s, now)
        title = stage or STAGE_NAMES.get(s.stage, s.stage)
        solver = SOLVER_NAMES.get(s.name, s.name) if s.name else ""
        detail = [solver] if solver else []
        if s.phase and (frac is None or (s.pct or 0) >= 100):   # before or after the time loop
            detail.append(s.phase)
        if eta is not None:
            detail.append(f"about {format_duration(eta)} left")
        elif note:
            detail.append(note)
        extra = describe_detail(s.detail)
        if extra:
            detail.append(extra)
        e = Estimate(frac, eta, est, "", title, " · ".join(detail))
        head = title + (f" ({solver})" if solver else "")
        e.text = " · ".join([head] + ([e.percent] if e.percent else []) + detail[1 if solver else 0:])
        return e

    @staticmethod
    def _solver_fraction(s: Solver, now: float):
        """(fraction, eta, estimated, note) from the best source the solver gave."""
        elapsed = max(0.0, now - s.start)
        if s.n > 1 and s.k > 0:
            frac = min(1.0, s.k / s.n)
            eta = None
            if s.t_first is not None and s.t_last is not None and s.k > s.k_first and s.t_last > s.t_first:
                rate = (s.t_last - s.t_first) / (s.k - s.k_first)
                eta = max(0.0, rate * (s.n - s.k) - max(0.0, now - s.t_last))
            elif s.t_last is not None and s.k > 0 and s.t_last > s.start:
                eta = max(0.0, (s.t_last - s.start) / s.k * (s.n - s.k) - max(0.0, now - s.t_last))
            return frac, eta, False, ""
        if s.pct is not None:
            eta = None
            if s.time_to_go is not None:
                since = max(0.0, now - s.pct_seen) if s.pct_seen is not None else 0.0
                eta = max(0.0, s.time_to_go - since)
            return min(1.0, s.pct / 100.0), eta, False, ""
        if s.loop_est_s:
            loop_start = s.start + (s.precomp_s or 0.0)
            done = max(0.0, now - loop_start)
            frac = min(0.99, done / s.loop_est_s)
            return frac, max(0.0, s.loop_est_s - done), True, ""
        if s.name in SILENT_SOLVERS:
            return None, None, False, f"{format_duration(elapsed)} elapsed (this solver reports no progress)"
        return None, None, False, f"{format_duration(elapsed)} elapsed"


def describe_detail(detail: str) -> str:
    """Thermal step detail 't=12.00s Tmax=37.85C' -> 'peak 37.85 °C at t = 12.0 s' (k-Wave:
    the maximum so far); 'Tnow=38.10C' -> 'now 38.10 °C' (BabelBrain: the current maximum)."""
    if not detail:
        return ""
    tnow = re.search(r"Tnow=([-\d.]+)C", detail)
    if tnow:
        return f"now {float(tnow.group(1)):.2f} °C"
    tmax = re.search(r"Tmax=([-\d.]+)C", detail)
    t = re.search(r"\bt=([\d.]+)s", detail)
    if tmax:
        return f"peak {float(tmax.group(1)):.2f} °C" + (f" at t = {float(t.group(1)):.1f} s" if t else "")
    return detail


def format_duration(seconds: float) -> str:
    seconds = max(0, int(round(seconds)))
    if seconds < 60:
        return f"{seconds} s"
    minutes, s = divmod(seconds, 60)
    if minutes < 60:
        return f"{minutes} min {s:02d} s" if minutes < 10 else f"{minutes} min"
    hours, minutes = divmod(minutes, 60)
    if hours < 48:
        return f"{hours} h {minutes:02d} min"
    return f"{hours // 24} d {hours % 24} h"


def relevant_lines(text: str) -> str:
    """The lines of a log that the parser reads (keeps job-log transfers small)."""
    return "\n".join(line for line in text.splitlines() if RELEVANT.match(line))
