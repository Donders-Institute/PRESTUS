"""Sonication protocol (timing.*): derived quantities and TUS calculator exchange.

PRESTUS describes a protocol in seconds (timing.pd, pri, ptd, ptri, ptrd and
the pulse ramp timing.pulse_ramp_dur / pulse_ramp_shape). The TUS calculator
(https://www.socsci.ru.nl/fusinitiative/tuscalculator/, Radboud University,
MIT licence) imports and exports its settings as JSON/YAML keyed by its field
labels, in milliseconds. Both use the same ramp convention: the ramp duration is
the rise (= fall) time at each end of a pulse.

PRESTUS models ramps on pulses only. Train and train-repetition ramps of the
calculator have no PRESTUS equivalent; they are exported as 0 and reported
when an imported file sets them.
"""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .fsio import write_text_atomic

TUS_CALCULATOR_URL = "https://www.socsci.ru.nl/fusinitiative/tuscalculator/"

# PRESTUS key (seconds) -> calculator label (milliseconds)
_MS_KEYS = {
    "timing.pd": "Pulse Duration [ms]",
    "timing.pulse_ramp_dur": "Ramp Duration [ms]",
    "timing.pri": "Pulse Repetition Interval [ms]",
    "timing.ptd": "Pulse Train Duration [ms]",
    "timing.ptri": "Pulse Train Repetition Interval [ms]",
    "timing.ptrd": "Pulse Train Repetition Duration [ms]",
}
_SHAPE = "Pulse Ramp Shape"
_TRAIN_ON = "Pulse Train is enabled"
_REP_ON = "Pulse Train Repetition is enabled"
_UNSUPPORTED = ("Pulse Train Ramp Duration [ms]", "Pulse Train Repetition Ramp Duration [ms]")
_FREQ = "Fundamental Frequency [MHz]"
_ISPPA = "Isppa [W/cm2]"
_DIAMETER = "Diameter of transducer [mm]"
_DEPTH = "Focal depth [mm]"


def _num(v: Any) -> float:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return 0.0
    return f if math.isfinite(f) else 0.0


@dataclass
class Protocol:
    pd: float = 0.0
    pri: float = 0.0
    ptd: float = 0.0
    ptri: float = 0.0
    ptrd: float = 0.0
    ramp: float = 0.0
    ramp_shape: str = "linear"

    @classmethod
    def from_model(cls, get) -> "Protocol":
        """get(path) -> value, e.g. ConfigModel.get."""
        shape = get("timing.pulse_ramp_shape") or "linear"
        return cls(pd=_num(get("timing.pd")), pri=_num(get("timing.pri")), ptd=_num(get("timing.ptd")),
                   ptri=_num(get("timing.ptri")), ptrd=_num(get("timing.ptrd")),
                   ramp=_num(get("timing.pulse_ramp_dur")),
                   ramp_shape=str(shape).lower())

    # ── levels ──────────────────────────────────────────────────────────
    @property
    def has_pulse(self) -> bool:
        return self.pd > 0

    @property
    def has_train(self) -> bool:
        return self.has_pulse and self.pri > 0 and self.ptd > 0

    @property
    def has_repetition(self) -> bool:
        return self.has_train and self.ptri > 0 and self.ptrd > 0

    # ── derived ─────────────────────────────────────────────────────────
    @property
    def prf_hz(self) -> float:
        return 1.0 / self.pri if self.pri > 0 else math.nan

    @property
    def pulses_per_train(self) -> int:
        return int(math.floor(self.ptd / self.pri + 1e-9)) if self.has_train else 1

    @property
    def trains(self) -> int:
        return int(math.floor(self.ptrd / self.ptri + 1e-9)) if self.has_repetition else 1

    @property
    def duty_cycle_pulse_train(self) -> float:
        return self.pd / self.pri if self.has_train else 1.0

    @property
    def duty_cycle_repetition(self) -> float:
        return self.ptd / self.ptri if self.has_repetition else 1.0

    @property
    def ramp_factor(self) -> float:
        """Isppa / Isptp for the ramped pulse (compute_ramp_factor.m)."""
        if self.ramp <= 0 or self.pd <= 0:
            return 1.0
        k = 3 / 8 if self.ramp_shape == "tukey" else 1 / 3
        return (2 * k * self.ramp + max(0.0, self.pd - 2 * self.ramp)) / self.pd

    @property
    def overall_duty_cycle(self) -> float:
        return self.duty_cycle_pulse_train * self.duty_cycle_repetition

    @property
    def total_duration_s(self) -> float:
        if self.has_repetition:
            return self.ptrd
        if self.has_train:
            return self.ptd
        return self.pd

    @property
    def sonication_time_s(self) -> float:
        """Time the transducer is on (sum of pulse durations)."""
        return self.pd * self.pulses_per_train * self.trains

    def problems(self) -> list[str]:
        """Inconsistent settings. Comparisons allow a relative 1e-9, so values that are equal
        on screen (e.g. 700 ms typed with a unit = 0.7000000000000001 s) are not flagged."""
        def gt(a: float, b: float) -> bool:
            return a > b * (1 + 1e-9) + 1e-15

        out = []
        if not self.has_pulse:
            out.append("No pulse duration: thermal simulations need timing.pd (and the train settings).")
            return out
        if self.pri and gt(self.pd, self.pri):
            out.append("Pulse duration exceeds the pulse repetition interval.")
        if self.ptd and self.pri and gt(self.pri, self.ptd):
            out.append("Pulse train is shorter than one pulse repetition interval.")
        if self.ptri and gt(self.ptd, self.ptri):
            out.append("Pulse train duration exceeds the train repetition interval.")
        if self.ptrd and self.ptri and gt(self.ptri, self.ptrd):
            out.append("Repetition duration is shorter than one train repetition interval.")
        if gt(2 * self.ramp, self.pd):
            out.append("Rise + fall ramp is longer than the pulse.")
        if self.ramp_shape not in ("linear", "tukey"):
            out.append("Pulse ramp shape must be 'linear' or 'tukey'.")
        return out


# ── TUS calculator exchange ─────────────────────────────────────────────

def to_tus_calculator(protocol: Protocol, freq_hz: float | None = None, isppa_wcm2: float | None = None,
                      diameter_mm: float | None = None, focal_depth_mm: float | None = None) -> dict:
    """Settings dict for the calculator's "Upload settings" (JSON)."""
    out: dict[str, Any] = {}
    if freq_hz:
        out[_FREQ] = round(freq_hz / 1e6, 6)
    if isppa_wcm2:
        out[_ISPPA] = float(isppa_wcm2)
    if diameter_mm:
        out[_DIAMETER] = float(diameter_mm)
    if focal_depth_mm:
        out[_DEPTH] = float(focal_depth_mm)
    values = {"timing.pd": protocol.pd, "timing.pulse_ramp_dur": protocol.ramp, "timing.pri": protocol.pri,
              "timing.ptd": protocol.ptd, "timing.ptri": protocol.ptri, "timing.ptrd": protocol.ptrd}
    for key, lab in _MS_KEYS.items():
        out[lab] = round(values[key] * 1000, 9)
    out[_SHAPE] = "Tukey" if protocol.ramp_shape == "tukey" else "Linear"
    out[_TRAIN_ON] = protocol.has_train
    out[_REP_ON] = protocol.has_repetition
    for lab in _UNSUPPORTED:
        out[lab] = 0
    return out


@dataclass
class Imported:
    timing: dict[str, Any] = field(default_factory=dict)      # dotted PRESTUS key -> value
    optional: dict[str, Any] = field(default_factory=dict)    # frequency / Isppa, applied on request
    notes: list[str] = field(default_factory=list)


def from_tus_calculator(settings: dict) -> Imported:
    """PRESTUS values from calculator settings (labels or its internal state keys)."""
    alias = {"pulse_duration_ms": "Pulse Duration [ms]", "pulse_ramp_duration_ms": "Ramp Duration [ms]",
             "pulse_repetition_interval_ms": "Pulse Repetition Interval [ms]",
             "pulse_train_duration_ms": "Pulse Train Duration [ms]",
             "pulse_train_repetition_interval_ms": "Pulse Train Repetition Interval [ms]",
             "pulse_train_repetition_duration_ms": "Pulse Train Repetition Duration [ms]",
             "pulse_ramp_shape": _SHAPE, "isSecondLevelTimingParameters": _TRAIN_ON,
             "isThirdLevelTimingParameters": _REP_ON, "frequency_MHz": _FREQ, "isppa": _ISPPA,
             "pulse_train_ramp_duration_ms": _UNSUPPORTED[0],
             "pulse_train_repetition_ramp_duration_ms": _UNSUPPORTED[1]}
    s = {alias.get(k, k): v for k, v in settings.items()}
    res = Imported()
    if not any(lab in s for lab in _MS_KEYS.values()):
        res.notes.append("No TUS calculator timing fields found in this file.")
        return res
    on = lambda k: str(s.get(k, True)).lower() in ("true", "1")  # noqa: E731
    train, rep = on(_TRAIN_ON), on(_REP_ON) and on(_TRAIN_ON)
    for key, lab in _MS_KEYS.items():
        if lab not in s:
            continue
        v = _num(s[lab]) / 1000
        if key in ("timing.pri", "timing.ptd") and not train:
            v = 0.0
        if key in ("timing.ptri", "timing.ptrd") and not rep:
            v = 0.0
        res.timing[key] = v
    if _SHAPE in s:
        res.timing["timing.pulse_ramp_shape"] = "tukey" if str(s[_SHAPE]).lower() == "tukey" else "linear"
    if not train:
        res.notes.append("Pulse train disabled in the calculator: PRI and PTD set to 0 (single pulse).")
    elif not rep:
        res.notes.append("Train repetition disabled in the calculator: PTRI and PTRD set to 0.")
    for lab in _UNSUPPORTED:
        if _num(s.get(lab, 0)) > 0:
            res.notes.append(f"'{lab}' = {s[lab]} is not modelled by PRESTUS (pulse ramps only) and was ignored.")
    if _FREQ in s and _num(s[_FREQ]) > 0:
        res.optional["transducer.freq_hz"] = _num(s[_FREQ]) * 1e6
    if _ISPPA in s and _num(s[_ISPPA]) > 0:
        res.optional["transducer.target_isppa_wcm2"] = _num(s[_ISPPA])
    return res


def read_settings(path: str | Path) -> dict:
    """The calculator's JSON, or its simple 'label: value' YAML export."""
    text = Path(path).read_text(encoding="utf-8")
    try:
        data = json.loads(text)
        if isinstance(data, dict):
            return data
    except json.JSONDecodeError:
        pass
    out: dict[str, Any] = {}
    for line in text.splitlines():
        if ":" not in line or line.lstrip().startswith("#"):
            continue
        k, _, v = line.rpartition(":")          # labels contain no ':'; values are plain
        k, v = k.strip().strip("'\""), v.strip().strip("'\"")
        out[k] = v
    return out


def write_settings(path: str | Path, settings: dict) -> Path:
    path = Path(path)
    write_text_atomic(path, json.dumps(settings, indent=2) + "\n")
    return path


# ── Editing from the diagram ────────────────────────────────────────────

# Protocol field -> config key
FIELD_KEYS = {"pd": "timing.pd", "ramp": "timing.pulse_ramp_dur", "pri": "timing.pri",
              "ptd": "timing.ptd", "ptri": "timing.ptri", "ptrd": "timing.ptrd"}
FIELD_NAMES = {"pd": "pulse duration", "ramp": "ramp duration", "pri": "pulse repetition interval",
               "ptd": "pulse train duration", "ptri": "pulse train repetition interval",
               "ptrd": "pulse train repetition duration"}
_UNITS = {"s": 1.0, "ms": 1e-3, "us": 1e-6, "µs": 1e-6, "μs": 1e-6, "min": 60.0}


def snap(v: float, digits: int = 2) -> float:
    """`v` rounded to `digits` significant digits (drag resolution)."""
    if v <= 0 or not math.isfinite(v):
        return 0.0
    return round(v, digits - 1 - int(math.floor(math.log10(v))))


def limits(p: Protocol, name: str) -> tuple[float, float]:
    """Range a dragged value may take so the protocol stays valid (see Protocol.problems)."""
    inf = math.inf
    if name == "pd":
        return 2 * p.ramp, (p.pri if p.has_train else inf)
    if name == "ramp":
        return 0.0, p.pd / 2
    if name == "pri":
        return p.pd, p.ptd
    if name == "ptd":
        return p.pri, (p.ptri if p.has_repetition else inf)
    if name == "ptri":
        return p.ptd, p.ptrd
    if name == "ptrd":
        return p.ptri, inf
    raise KeyError(name)


def dragged(p: Protocol, name: str, value: float) -> float:
    """A dragged value: snapped, then kept inside `limits`."""
    lo, hi = limits(p, name)
    v = snap(value)
    return min(max(v, lo), hi)


def parse_duration(text: str, hz: bool = False) -> float | None:
    """Seconds from '20 ms', '0.5', '500 us', '2 min' (no unit = seconds); with hz=True also
    '5 Hz' / '1 kHz' (a repetition frequency, returned as its interval). None if unreadable."""
    t = text.strip().replace(",", ".").lower()
    m = re.fullmatch(r"([0-9]*\.?[0-9]+(?:e[-+]?[0-9]+)?)\s*([a-zµμ]*)", t)
    if not m:
        return None
    v, unit = float(m.group(1)), m.group(2)
    if hz and unit in ("hz", "khz"):
        f = v * (1e3 if unit == "khz" else 1.0)
        return 1.0 / f if f > 0 else None
    scale = _UNITS.get(unit or "s")
    return None if scale is None else float(f"{v * scale:.12g}")     # 700 ms -> 0.7, not 0.7000000000000001


# ── Async firing schedule (timing.async_schedule, thermal_schedule.m) ─

ASYNC_MODES = ("pulses", "trains", "custom")
TIMING_OVERRIDES = ("pd", "pri", "ptd", "ptri", "ptrd")  # transducer.timing.<field>, read in 'custom'
_MAX_WINDOWS = 200_000                               # per transducer, for the overlap check


@dataclass
class SourceTiming:
    """One transducer's protocol and start delay [s]."""
    protocol: Protocol
    offset: float = 0.0
    freq_hz: float = 0.0
    overridden: tuple[str, ...] = ()

    def pulse_windows(self, t_max: float = math.inf, limit: int = _MAX_WINDOWS) -> list[tuple[float, float]]:
        """(start, end) of each pulse starting before t_max, as thermal_schedule.m lays them out."""
        p = self.protocol
        if not p.has_pulse:
            return []
        n_pulses = p.pulses_per_train if p.has_train else 1
        n_trains = p.trains if p.has_repetition else 1
        out = []
        for r in range(n_trains):
            for k in range(n_pulses):
                s = self.offset + r * (p.ptri if p.has_repetition else 0.0) + k * (p.pri if p.has_train else 0.0)
                if s >= t_max or len(out) >= limit:
                    return out
                out.append((s, s + p.pd))
        return out

    def train_windows(self, t_max: float = math.inf) -> list[tuple[float, float]]:
        p = self.protocol
        if not p.has_train:
            return [(self.offset, self.offset + p.pd)] if p.has_pulse else []
        n_trains = p.trains if p.has_repetition else 1
        out = []
        for r in range(n_trains):
            s = self.offset + r * (p.ptri if p.has_repetition else 0.0)
            if s >= t_max:
                break
            out.append((s, s + p.ptd))
        return out

    @property
    def overall_duty_cycle(self) -> float:
        return self.protocol.overall_duty_cycle


@dataclass
class AsyncSchedule:
    """How several transducers take turns firing (async coupling)."""
    mode: str
    sources: list[SourceTiming]

    @classmethod
    def build(cls, base: Protocol, mode: str, n: int, offsets: list | None = None,
              overrides: list[dict] | None = None, freqs: list | None = None) -> "AsyncSchedule":
        """Offsets as thermal_schedule.m: 'pulses' (i-1)·PRI/N, 'trains' (i-1)·PTRI/N,
        'custom' the given offsets and per-transducer overrides of pd/pri/ptd/ptri/ptrd."""
        mode = mode if mode in ASYNC_MODES else "pulses"
        offsets = list(offsets or [])
        overrides = list(overrides or [])
        freqs = list(freqs or [])
        sources = []
        for i in range(n):
            ov = {k: v for k, v in (overrides[i] if i < len(overrides) and overrides[i] else {}).items()
                  if k in TIMING_OVERRIDES and _num(v) > 0}
            p = base
            if mode == "custom" and ov:
                p = Protocol(**{**base.__dict__, **{k: _num(v) for k, v in ov.items()}})
            if mode == "pulses":
                off = i * base.pri / n if base.has_train else 0.0
            elif mode == "trains":
                off = i * base.ptri / n if base.has_repetition else 0.0
            else:
                off = _num(offsets[i]) if i < len(offsets) else 0.0
            sources.append(SourceTiming(p, max(0.0, off), _num(freqs[i]) if i < len(freqs) else 0.0,
                                        tuple(sorted(ov))))
        return cls(mode, sources)

    @classmethod
    def from_model(cls, model) -> "AsyncSchedule | None":
        """The schedule of a config model with async coupling and several transducers, else None."""
        n = model.n_transducers
        if n < 2 or model.get("simulation.transducer_coupling") != "async":
            return None
        get = model.get
        base = Protocol.from_model(get)
        mode = str(get("timing.async_schedule") or "pulses")
        offsets = [get("transducer.timing_offset_s", i) for i in range(n)]
        overrides = [{k: get(f"transducer.timing.{k}", i) for k in TIMING_OVERRIDES} for i in range(n)]
        freqs = [get("transducer.freq_hz", i) for i in range(n)]
        return cls.build(base, mode, n, offsets, overrides, freqs)

    @property
    def n(self) -> int:
        return len(self.sources)

    def cycle(self) -> float:
        """Time span the diagram shows: two PRIs ('pulses'), one PTRI ('trains'), or until
        every transducer has started a full cycle ('custom')."""
        p = self.sources[0].protocol
        if self.mode == "pulses":
            span = 2 * p.pri if p.has_train else p.pd
        elif self.mode == "trains":
            span = p.ptri if p.has_repetition else p.ptd
        else:
            span = max(s.offset + (s.protocol.ptri if s.protocol.has_repetition else
                                   s.protocol.pri if s.protocol.has_train else s.protocol.pd)
                       for s in self.sources)
        total = max((s.offset + s.protocol.total_duration_s for s in self.sources), default=span)
        return min(span, total) if total > 0 else span

    def blocks(self, i: int, t_win: float, max_blocks: int = 60) -> tuple[list[tuple[float, float]], str]:
        """On-blocks of transducer i within [0, t_win]: its pulses, or its trains when more
        pulses than max_blocks would be drawn."""
        s = self.sources[i]
        w = s.pulse_windows(t_win, max_blocks + 1)
        if len(w) <= max_blocks:
            return [(a, min(b, t_win)) for a, b in w], "pulses"
        return [(a, min(b, t_win)) for a, b in s.train_windows(t_win)], "trains"

    def problems(self, chained: bool = False) -> list[str]:
        """The checks of thermal_schedule.m that would stop the run. ``chained``: the
        sequence is run as a chain of runs (sequence.py), where transducers that start
        together form one coherent run, so their overlap is not an error."""
        out = []
        base = self.sources[0].protocol
        if not base.has_pulse:
            return out
        n = self.n
        tight = lambda a, b: a > b * (1 + 1e-9) + 1e-15        # noqa: E731
        if self.mode != "custom":
            over = [i + 1 for i, s in enumerate(self.sources) if s.overridden]
            if over:
                out.append(f"Transducer{'s' if len(over) > 1 else ''} {', '.join(map(str, over))} "
                           f"override{'' if len(over) > 1 else 's'} the timing, which only the "
                           "'custom' schedule reads.")
        if self.mode == "pulses" and base.has_train and tight(n * base.pd, base.pri):
            out.append(f"'pulses' fires the {n} transducers in turn within each PRI: needs "
                       f"{n} × PD ≤ PRI ({n} × {base.pd * 1e3:.4g} ms > {base.pri * 1e3:.4g} ms). "
                       "Use 'trains', a shorter PD, or 'custom'.")
        if self.mode == "trains":
            if not base.has_repetition:
                out.append("'trains' needs a train repetition interval (PTRI and PTRD).")
            elif tight(n * base.ptd, base.ptri):
                out.append(f"'trains' fires the {n} transducers in turn within each PTRI: needs "
                           f"{n} × PTD ≤ PTRI ({n} × {base.ptd:.4g} s > {base.ptri:.4g} s). "
                           "Use 'pulses', a shorter PTD, or 'custom'.")
        same, _ = self._overlap_pairs()
        if same and not chained:
            out.append(f"{', '.join(f'T{i + 1}/T{j + 1}' for i, j in same)} overlap at the same "
                       "frequency: their fields interfere, which async coupling does not model. "
                       "Give them the same start and timing to fire them together, change the "
                       "offsets, or use coherent coupling.")
        return out

    def notes(self) -> list[str]:
        """What thermal_schedule.m warns about without stopping the run."""
        _, other = self._overlap_pairs()
        if not other:
            return []
        return [f"{', '.join(f'T{i + 1}/T{j + 1}' for i, j in other)} overlap at different "
                "frequencies: their intensities are summed (valid when a pulse is long against "
                "1/|f1 − f2|)."]

    def groups(self) -> list[int]:
        """Firing group of each transducer (0-based), as firing_groups.m: on the 'custom'
        schedule, transducers with the same offset and the same timing fire together (their
        fields are summed coherently); otherwise every transducer is its own group."""
        if self.mode != "custom":
            return list(range(self.n))
        out: list[int] = []
        for i, s in enumerate(self.sources):
            for j in range(i):
                t = self.sources[j]
                if abs(s.offset - t.offset) <= 1e-6 and all(
                        abs(getattr(s.protocol, f) - getattr(t.protocol, f))
                        <= 1e-9 * max(1.0, abs(getattr(t.protocol, f))) for f in TIMING_OVERRIDES):
                    out.append(out[j])
                    break
            else:
                out.append(max(out, default=-1) + 1)
        return out

    def _overlap_pairs(self) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
        """Overlapping pairs ('custom' only) at the same frequency, and at different ones. As
        thermal_schedule.m: equal known frequencies, or none known, count as the same.
        Transducers of one firing group (groups()) fire together by design: not a pair."""
        if self.mode != "custom" or not self.sources[0].protocol.has_pulse:
            return [], []
        grp = self.groups()
        clash = [(i, j) for i, j in self.overlaps() if grp[i] != grp[j]]
        same = [(i, j) for i, j in clash if self.sources[i].freq_hz == self.sources[j].freq_hz]
        return same, [pq for pq in clash if pq not in same]

    def overlaps(self) -> list[tuple[int, int]]:
        """Pairs of transducers whose pulses overlap in time."""
        events = sorted((a, b, i) for i, s in enumerate(self.sources) for a, b in s.pulse_windows())
        tol = 1e-9 * max((b for _, b, _ in events), default=1.0)
        pairs: set[tuple[int, int]] = set()
        active: list[tuple[float, int]] = []          # (end, source) of pulses still on
        for a, b, i in events:
            active = [(e, j) for e, j in active if e > a + tol]
            for _, j in active:
                if j != i:
                    pairs.add((min(i, j), max(i, j)))
            active.append((b, i))
        return sorted(pairs)
