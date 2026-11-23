"""Sonication sequences: stances fired one after another, with a pause between them.

A sequence is stored with PRESTUS' own keys: every stance is a transducer of an async run
with the 'custom' schedule (simulation.transducer_coupling: async, timing.async_schedule:
custom). transducer(i).timing_offset_s is when its step starts, transducer(i).timing.* its
pulse timing and duration where they differ from the timing block, which holds the protocol
of the step that contains transducer 1 (thermal_schedule.m). Transducers that start at the
same time form one step and fire together.

Two ways to simulate a plan (plan() picks one and says why):

- joint: the config runs as it is, one async run. Each stance gets its own acoustic
  simulation; one thermal simulation covers the whole sequence (async_transducer_pipeline.m).
  The transducers of a step fire together as a firing group (firing_groups.m): their complex
  fields are summed, so they interfere (on per-transducer grids after resampling them with
  the carrier removed, remap_complex_field.m).
- chain: one config per step, run as sequential follow-ups (sequential_pipeline.m). A step
  with several transducers is one coherent run. Each step starts from the previous step's
  end temperature and CEM43 maps, written after its post_ptri_dur (thermal_simulation.m),
  so the pause before step k+1 becomes step k's timing.post_ptri_dur. Each step gets its
  own grid (io.preproc_affix); sequential_pipeline.m would otherwise reuse step 1's grid.

Qt-free.
"""

from __future__ import annotations

import copy
import datetime as _dt
import math
from dataclasses import dataclass, field
from pathlib import Path

from .fsio import write_text_atomic
from .merge import mergestruct
from .model import TX, ConfigModel, Issue
from .study import sanitize_affix
from .timing import TIMING_OVERRIDES, AsyncSchedule, Protocol

TOL = 1e-6                                   # s: offsets closer than this start together
BACKENDS = ("auto", "joint", "chain")
BACKEND_LABELS = {"auto": "Automatic", "joint": "One joint run", "chain": "Chain of runs"}
SESSION_KEY = "sequence_backend"             # session: resolved study path -> backend
STEP_KEYS = ("timing_offset_s", "timing")   # per-transducer keys a step config drops


class SequenceError(ValueError):
    """The plan cannot be run with the chosen backend (its errors, joined)."""


def _num(v) -> float | None:
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
        return None
    return float(v)


def _close(a, b) -> bool:
    a, b = _num(a), _num(b)
    if a is None or b is None:
        return a is None and b is None
    return abs(a - b) <= TOL * max(1.0, abs(a), abs(b))


def fmt_s(v: float) -> str:
    """'90 s', '2.5 min', '1.5 h', '250 ms'."""
    if v is None or not math.isfinite(v):
        return "–"
    a = abs(v)
    if a == 0:
        return "0 s"
    if a < 1:
        return f"{v * 1e3:.3g} ms"
    if a < 120:
        return f"{v:.3g} s"
    if a < 7200:
        return f"{v / 60:.3g} min"
    return f"{v / 3600:.3g} h"


# ── reading ────────────────────────────────────────────────────────────────────

def is_sequence(model: ConfigModel) -> bool:
    """Several transducers on the async 'custom' schedule (thermal_schedule.m)."""
    return (model.n_transducers > 1 and model.get("simulation.transducer_coupling") == "async"
            and str(model.get("timing.async_schedule") or "") == "custom")


@dataclass
class Step:
    transducers: list[int]                   # 0-based transducer indices, ascending
    start_s: float
    protocol: Protocol
    freqs: list[float] = field(default_factory=list)
    mixed: bool = False                      # its transducers have different timing

    @property
    def duration_s(self) -> float:
        """How long the step sonicates: its PTRD, as one run spans timing.ptrd."""
        return self.protocol.total_duration_s

    @property
    def end_s(self) -> float:
        return self.start_s + self.duration_s

    @property
    def together(self) -> bool:
        return len(self.transducers) > 1

    def timing(self) -> dict[str, float]:
        return {f: float(getattr(self.protocol, f)) for f in TIMING_OVERRIDES}


def _same_protocol(a: Protocol, b: Protocol) -> bool:
    return all(_close(getattr(a, f), getattr(b, f)) for f in TIMING_OVERRIDES)


def steps(model: ConfigModel) -> list[Step]:
    """The steps of a sequence in start order ([] when the config is not one)."""
    if not is_sequence(model):
        return []
    sched = AsyncSchedule.from_model(model)
    order = sorted(range(sched.n), key=lambda i: (sched.sources[i].offset, i))
    out: list[Step] = []
    for i in order:
        src = sched.sources[i]
        if out and abs(src.offset - out[-1].start_s) <= TOL:
            step = out[-1]
            step.transducers.append(i)
            step.freqs.append(src.freq_hz)
            step.mixed = step.mixed or not _same_protocol(step.protocol, src.protocol)
        else:
            out.append(Step([i], src.offset, src.protocol, [src.freq_hz]))
    for s in out:
        s.transducers.sort()
    return out


def pauses(stps: list[Step]) -> list[float]:
    """Pause after each step but the last [s]; negative when the next starts early."""
    return [b.start_s - a.end_s for a, b in zip(stps, stps[1:])]


# ── planning ───────────────────────────────────────────────────────────────────

@dataclass
class Plan:
    steps: list[Step]
    pauses: list[float]
    post_s: float                            # timing.post_ptri_dur, after the last step
    preference: str                          # 'auto' | 'joint' | 'chain'
    backend: str                             # 'joint' | 'chain' | '' (cannot run)
    reason: str
    issues: list[Issue] = field(default_factory=list)

    @property
    def n(self) -> int:
        return len(self.steps)

    @property
    def total_s(self) -> float:
        return (max(s.end_s for s in self.steps) if self.steps else 0.0) + max(0.0, self.post_s)

    @property
    def errors(self) -> list[Issue]:
        return [i for i in self.issues if i.level == "error"]


def _issue(level: str, path: str, message: str, ti: int | None = None) -> Issue:
    return Issue(level, path, message, ti)


def plan(model: ConfigModel, preference: str = "auto") -> Plan | None:
    """How the sequence of ``model`` is simulated, with what stops it (None: no sequence)."""
    stps = steps(model)
    if not stps:
        return None
    preference = preference if preference in BACKENDS else "auto"
    g = model.get
    gaps = pauses(stps)
    post = _num(g("timing.post_ptri_dur")) or 0.0
    issues: list[Issue] = []
    heating = g("modules.run_heating_sims") in (1, True)

    together = [k for k, s in enumerate(stps) if s.together]
    overlap = [k for k, p in enumerate(gaps) if p < -TOL]
    per_tx_grid = g("grid.async_layout") == "per_transducer"

    # what each backend needs: a step firing together is summed coherently in one joint
    # run (firing_groups.m); a chain runs it as one coherent run
    if preference == "chain":
        backend = "chain"
        reason = ("Chosen here: the steps run one after another, each on its own grid, and each "
                  "starts from the previous step's end temperature and thermal dose")
    else:
        backend = "joint"
        reason = ("Each transducer gets its own acoustic simulation, transducers that start "
                  "together are summed coherently, and one thermal simulation covers the whole "
                  "sequence" if together else
                  "Every step is a single stance: each gets its own acoustic simulation and one "
                  "thermal simulation covers the whole sequence")

    for s in stps:
        if s.mixed:
            names = ", ".join(f"T{i + 1}" for i in s.transducers)
            issues.append(_issue("error", "transducer.timing", f"{names} start together but have "
                                 "different pulse timing. Transducers that fire together share one "
                                 "timing; give them the same timing or different start times",
                                 s.transducers[-1]))
        p = s.protocol
        if heating and p.has_repetition and abs(p.ptrd / p.ptri - round(p.ptrd / p.ptri)) > 1e-6:
            issues.append(_issue("error", "transducer.timing.ptrd",
                                 f"Step {stps.index(s) + 1} sonicates {fmt_s(p.ptrd)}, which is not a "
                                 f"whole number of pulse train repetitions ({fmt_s(p.ptri)}); "
                                 "thermal_parameters.m stops the run", s.transducers[0]))

    if backend == "chain":
        for k in overlap:
            issues.append(_issue("error", "transducer.timing_offset_s",
                                 f"Step {k + 2} starts {fmt_s(-gaps[k])} before step {k + 1} ends. "
                                 "A chain needs each step to finish first", stps[k + 1].transducers[0]))
        dt = _num(g("timing.post_pt_timestep")) or 0.0
        for k, p in enumerate(gaps):
            if p > TOL and (dt <= 0 or abs(p / dt - round(p / dt)) > 1e-6):
                issues.append(_issue("error", "transducer.timing_offset_s",
                                     f"The pause after step {k + 1} ({fmt_s(p)}) must be a whole number "
                                     f"of thermal time steps (timing.post_pt_timestep = {fmt_s(dt)}) "
                                     "in a chain, where it is that step's post_ptri_dur",
                                     stps[k + 1].transducers[0]))
        if g("simulation.uncertainty") is True:
            issues.append(_issue("error", "simulation.uncertainty",
                                 "Uncertainty mode cannot run a sequence as a chain: its variants "
                                 "share one preprocessing cache (uncertainty_pipeline.m), and the "
                                 "steps' stances need their own"))
        multi = [i + 1 for i, t in enumerate(model.effective_transducers())
                 if isinstance(t.get("target_isppa_wcm2"), list)
                 and len([x for x in t["target_isppa_wcm2"] if _num(x)]) > 1]
        if multi:
            issues.append(_issue("error", "transducer.target_isppa_wcm2",
                                 "A chain passes the end state of plain runs only; a list of target "
                                 f"ISPPAs (transducer {multi[0]}) makes a step a multi-ISPPA run",
                                 multi[0] - 1))
    elif backend == "joint":
        if per_tx_grid:
            issues.append(_issue("warning", "grid.async_layout",
                                 "With one grid per transducer, the heat of later stances is "
                                 "resampled onto transducer 1's grid and dropped outside it. The "
                                 "shared grid, or a chain of runs, keeps it"))
        if per_tx_grid and together:
            issues.append(_issue("warning", "grid.async_layout",
                                 f"Step {together[0] + 1} fires several transducers together: their "
                                 "complex fields are resampled onto transducer 1's grid with the "
                                 "carrier removed, which is accurate along the beam but less so in "
                                 "the skull. The shared grid needs no resampling"))

    if not heating:
        issues.append(_issue("warning", "modules.run_heating_sims",
                             "Heating is off: the pauses and durations of a sequence only matter "
                             "for the thermal simulation"))
    order = {"error": 0, "warning": 1, "info": 2}
    issues.sort(key=lambda i: order.get(i.level, 3))
    return Plan(stps, gaps, post, preference, backend, reason, issues)


# ── editing ────────────────────────────────────────────────────────────────────

@dataclass
class StepSpec:
    """A step as the card edits it: who fires, with which timing, and the pause after it."""
    transducers: list[int]
    timing: dict[str, float]
    pause_after: float = 0.0


def specs(model: ConfigModel) -> list[StepSpec]:
    stps = steps(model)
    gaps = pauses(stps) + [0.0]
    return [StepSpec(list(s.transducers), s.timing(), max(0.0, p)) for s, p in zip(stps, gaps)]


def write(model: ConfigModel, spec: list[StepSpec], label: str = "Edit sequence") -> None:
    """Lay the steps out back to back with their pauses: start times into
    transducer.timing_offset_s, the timing of the step with transducer 1 into the timing
    block, and every other step's differences from it into transducer.timing.*."""
    with model.edit(label):
        model.set("simulation.transducer_coupling", "async")
        model.set("timing.async_schedule", "custom")
        assigned: dict[int, dict[str, float]] = {}
        start = 0.0
        for s in spec:
            for ti in s.transducers:
                assigned[ti] = {"timing_offset_s": start,
                                **{f"timing.{f}": float(s.timing[f]) for f in TIMING_OVERRIDES}}
            start += float(s.timing.get("ptrd", 0.0)) + max(0.0, float(s.pause_after))
        if 0 in assigned:
            for f in TIMING_OVERRIDES:
                model.set(f"timing.{f}", assigned[0][f"timing.{f}"])
        for ti in sorted(assigned):
            for key, value in assigned[ti].items():
                if ti == 0:
                    inherited = 0.0 if key == "timing_offset_s" else _num(model.get(key))
                else:
                    inherited = assigned.get(0, {}).get(key)
                path = f"{TX}.{key}"
                if _close(value, inherited):
                    if model.is_overridden(path, ti):
                        model.reset(path, ti)
                else:
                    model.set(path, value, ti)


def make_sequence(model: ConfigModel) -> None:
    """Fire the existing transducers one after another (each its own step, no pause),
    keeping a sequence that is already set up."""
    if is_sequence(model):
        return
    sched_timing = Protocol.from_model(model.get)
    base = {f: float(getattr(sched_timing, f)) for f in TIMING_OVERRIDES}
    spec = []
    for ti in range(model.n_transducers):
        t = dict(base)
        for f in TIMING_OVERRIDES:
            v = _num(model.get(f"{TX}.timing.{f}", ti)) if ti else None
            if v:
                t[f] = v
        spec.append(StepSpec([ti], t, 0.0))
    write(model, spec, "Fire the transducers one after another")


def add_step(model: ConfigModel, pause_s: float = 0.0) -> int:
    """Append a stance as a new last step with the timing of the current last step;
    return its transducer index (placed where transducer 1 is, until moved)."""
    with model.edit("Add a step"):
        if model.n_transducers == 1:
            ti = model.add_transducer()
            last = StepSpec([0], {f: float(getattr(Protocol.from_model(model.get), f))
                                  for f in TIMING_OVERRIDES})
            spec = [last]
        else:
            make_sequence(model)
            spec = specs(model)
            ti = model.add_transducer()
        spec[-1].pause_after = max(0.0, float(pause_s))
        spec.append(StepSpec([ti], dict(spec[-1].timing)))
        write(model, spec, "Add a step")
    return ti


def set_pause(model: ConfigModel, k: int, seconds: float) -> None:
    spec = specs(model)
    spec[k].pause_after = max(0.0, float(seconds))
    write(model, spec, f"Pause after step {k + 1}")


def set_timing(model: ConfigModel, k: int, values: dict[str, float]) -> None:
    """Change step k's timing (pd, pri, ptd, ptri, ptrd); later steps move, pauses stay."""
    spec = specs(model)
    spec[k].timing.update({f: float(v) for f, v in values.items() if f in TIMING_OVERRIDES})
    write(model, spec, f"Timing of step {k + 1}")


def whole_trains(duration: float, ptri: float) -> float:
    """``duration`` rounded to a whole number (>= 1) of pulse train repetitions."""
    if not ptri or ptri <= 0:
        return duration
    return max(1, math.floor(duration / ptri + 0.5)) * ptri     # half up, not to even


def join_previous(model: ConfigModel, k: int) -> None:
    """Fire step k together with step k-1 (with step k-1's timing)."""
    spec = specs(model)
    if k <= 0 or k >= len(spec):
        return
    spec[k - 1].transducers += spec[k].transducers
    spec[k - 1].pause_after = spec[k].pause_after
    del spec[k]
    write(model, spec, f"Fire step {k + 1} with step {k}")


def separate(model: ConfigModel, k: int) -> None:
    """Give each transducer of step k its own step, one after another."""
    spec = specs(model)
    s = spec[k]
    if len(s.transducers) < 2:
        return
    parts = [StepSpec([ti], dict(s.timing)) for ti in s.transducers]
    parts[-1].pause_after = s.pause_after
    spec[k:k + 1] = parts
    write(model, spec, f"Separate step {k + 1}")


def remove_step(model: ConfigModel, k: int) -> None:
    """Remove step k and its transducers; not the step with transducer 1."""
    spec = specs(model)
    if not 0 <= k < len(spec) or 0 in spec[k].transducers:
        return
    gone = sorted(spec[k].transducers)
    del spec[k]                    # the pause before it stays; later steps move earlier
    with model.edit(f"Remove step {k + 1}"):
        for ti in reversed(gone):
            model.remove_transducer(ti)
        for s in spec:
            s.transducers = [ti - sum(1 for g in gone if g < ti) for ti in s.transducers]
        if model.n_transducers > 1:
            write(model, spec, f"Remove step {k + 1}")
        else:
            for key in ("timing_offset_s",) + tuple(f"timing.{f}" for f in TIMING_OVERRIDES):
                if model.is_overridden(f"{TX}.{key}", 0):
                    model.reset(f"{TX}.{key}", 0)


# ── chain of runs ──────────────────────────────────────────────────────────────

def _raw_affix(affix: str) -> str:
    return affix[len("_desc-"):] if affix.startswith("_desc-") else affix.lstrip("_")


def step_affix(model: ConfigModel, k: int) -> str:
    """Raw io.output_affix of step k (1-based): the plan's label plus 'step<k>'."""
    base = _raw_affix(sanitize_affix(model.get("io.output_affix") or ""))
    return f"{base}_step{k}" if base else f"step{k}"


def step_path(plan_path: Path, k: int) -> Path:
    plan_path = Path(plan_path)
    return plan_path.with_name(f"{plan_path.stem}_step-{k}{plan_path.suffix or '.yaml'}")


def chain_configs(model: ConfigModel, pl: Plan) -> list[ConfigModel]:
    """One plain config per step (step 1 runs first; the rest are its sequential
    follow-ups). Raises SequenceError when the plan has errors or is not a chain."""
    if pl.backend != "chain" or pl.errors:
        raise SequenceError("; ".join(i.message for i in pl.errors) or "The sequence is not a chain")
    out = []
    created = _dt.datetime.now().replace(microsecond=0).isoformat()
    for k, (s, after) in enumerate(zip(pl.steps, pl.pauses + [pl.post_s]), start=1):
        child = ConfigModel(model.schema, model.library)
        child.load_dict(model.override_dict())
        txs = []
        for ti in s.transducers:
            t = copy.deepcopy(model.tx[0]) if ti == 0 else mergestruct(model.tx[0], model.tx[ti])
            for key in STEP_KEYS:
                t.pop(key, None)
            txs.append(t)
        child.tx = txs or [{}]
        child.tx_as_list = len(txs) > 1
        for f in TIMING_OVERRIDES:
            child.set(f"timing.{f}", float(getattr(s.protocol, f)))
        child.set("timing.post_ptri_dur", max(0.0, float(after)))
        child.reset("timing.async_schedule")
        child.reset("grid.async_layout")
        if s.together:
            child.set("simulation.transducer_coupling", "coherent")
        else:
            child.reset("simulation.transducer_coupling")
        child.set("modules.run_heating_sims", 1)
        child.set("modules.run_nifti_creation", 1)
        for key in ("io.acoustic_cache_affix", "io.thermal_cache_affix", "io.preproc_affix"):
            child.reset(key)
        affix = step_affix(model, k)
        child.set("io.output_affix", affix)
        if k > 1:
            # own grid: sequential_pipeline.m would reuse step 1's (another stance)
            child.set("io.preproc_affix", sanitize_affix(affix))
            for key in ("io.adopted_heatmap", "io.adopted_cem43", "io.adopted_cem43_iso"):
                child.reset(key)                       # set by sequential_pipeline.m
        parent = str(model.path.resolve()) if model.path is not None else ""
        child.lineage = {"parent_config": parent, "type": "sequence_step", "step": k,
                         "steps": pl.n, "created": created}
        out.append(child)
    return out


HEADER = ("Step {k} of {n} of the sonication sequence in {plan}. Written by the PRESTUS Nexus "
          "when the sequence is run; edit the sequence there, not this file.")


def write_chain(model: ConfigModel, pl: Plan) -> list[Path]:
    """Write the step configs next to the plan config; return their paths in order."""
    if model.path is None:
        raise SequenceError("Save the config first")
    paths = []
    for k, child in enumerate(chain_configs(model, pl), start=1):
        p = step_path(model.path, k)
        write_text_atomic(p, child.to_yaml(HEADER.format(k=k, n=pl.n, plan=model.path.name)))
        child.path = p
        paths.append(p)
    return paths


def launch_target(model: ConfigModel, options, preference: str = "auto"):
    """(config to run, pipeline options, notes) for the Run tab. A chained sequence runs
    its step 1 with steps 2.. (and then the study's own follow-ups) as its sequential
    chain; everything else runs the config itself. Raises SequenceError for a sequence
    that cannot run."""
    pl = plan(model, preference)
    if pl is None:
        return model.path, options, []
    if pl.errors:
        raise SequenceError("The sonication sequence cannot run: " +
                            "; ".join(i.message for i in pl.errors))
    if pl.backend != "chain":
        return model.path, options, []
    paths = write_chain(model, pl)
    opts = copy.deepcopy(options)
    opts.sequential_configs = [str(p.resolve()) for p in paths[1:]] + \
        [c for c in options.sequential_configs if c not in {str(p.resolve()) for p in paths}]
    notes = [f"Sequence: {pl.n} steps as a chain of runs, starting with {paths[0].name}"]
    return paths[0], opts, notes


# ── backend preference (Nexus setting, per study, in the session) ───────────────

def load_backend(session, study: Path | None) -> str:
    if session is None or study is None:
        return "auto"
    store = session.get(SESSION_KEY) or {}
    v = store.get(str(Path(study).resolve())) if isinstance(store, dict) else None
    return v if v in BACKENDS else "auto"


def save_backend(session, study: Path, backend: str) -> None:
    store = dict(session.get(SESSION_KEY) or {})
    key = str(Path(study).resolve())
    if backend == "auto" or backend not in BACKENDS:
        store.pop(key, None)
    else:
        store[key] = backend
    session.set(SESSION_KEY, store)


def sequence_issues(model: ConfigModel) -> list[Issue]:
    """Issues of the sequence under the automatic backend (checks.all_issues)."""
    pl = plan(model)
    return list(pl.issues) if pl is not None else []


def describe(model: ConfigModel, preference: str = "auto") -> str | None:
    """One line for the Run tab: how the sequence runs ('3 steps as a chain of 3 runs; step 2
    fires T2 + T3 together'); None without a sequence."""
    pl = plan(model, preference)
    if pl is None:
        return None
    together = [f"step {k + 1} fires " + " + ".join(f"T{i + 1}" for i in s.transducers) + " together"
                for k, s in enumerate(pl.steps) if s.together]
    how = {"joint": "as one joint run", "chain": f"as a chain of {pl.n} runs"}.get(pl.backend, "(cannot run)")
    order = ("overlapping in time" if any(p < -TOL for p in pl.pauses) else "one after another")
    return f"{pl.n} steps, {order}, {how}" + ("; " + ", ".join(together) if together else "")


def has_together_steps(model: ConfigModel) -> bool:
    return any(s.together for s in steps(model))


__all__ = ["Plan", "Step", "StepSpec", "SequenceError", "is_sequence", "steps", "plan", "specs",
           "write", "make_sequence", "add_step", "set_pause", "set_timing", "join_previous",
           "separate", "remove_step", "chain_configs", "write_chain", "launch_target",
           "load_backend", "save_backend", "sequence_issues", "whole_trains", "fmt_s"]
