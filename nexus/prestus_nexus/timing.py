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
        out = []
        if not self.has_pulse:
            out.append("No pulse duration: thermal simulations need timing.pd (and the train settings).")
            return out
        if self.pri and self.pd > self.pri:
            out.append("Pulse duration exceeds the pulse repetition interval.")
        if self.ptd and self.pri and self.ptd < self.pri:
            out.append("Pulse train is shorter than one pulse repetition interval.")
        if self.ptri and self.ptd > self.ptri:
            out.append("Pulse train duration exceeds the train repetition interval.")
        if self.ptrd and self.ptri and self.ptrd < self.ptri:
            out.append("Repetition duration is shorter than one train repetition interval.")
        if 2 * self.ramp > self.pd:
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
