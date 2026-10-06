"""Calibration page logic without Qt: depth lists, the study's own focal depths, fit settings
to and from the library entry's meta.fit (calibration_save_entry.m)."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Any

from . import placement as pl
from .model import TX, ConfigModel, is_empty

# Advanced FIT_TRANSDUCER_MODEL options and their fitter defaults, as advanced_defaults()
# in calibration_gui.m; only values that differ are sent.
ADVANCED_DEFAULTS: dict[str, Any] = {
    "PhaseConstraint": "none", "PhaseStepBounds": [], "AmplitudeFit": "shape",
    "SharedVelocity": False, "Optimizer": "lsqnonlin", "PopSize": 5000, "MaxFunEvals": 200000,
    "Polish": True, "FocalWeight": 0, "FitWindowMm": [], "Lambda": 1e-3, "LambdaAmp": 1e-2,
    "ZShiftBounds": [-3, 3],
}
_MANIFEST_EP = re.compile(r"^transducer(?:\((\d+)\))?\.focal_distance_ep$")


def parse_depths(text: str) -> list[float] | None:
    """'a:b:c', 'a:c' (step 1) or a comma/space-separated list, as calibration_gui;
    rounded to 0.01 mm, unique, sorted. None when the text is not a depth list."""
    text = (text or "").strip()
    if not text:
        return None
    parts = text.split(":")
    try:
        if len(parts) == 1:
            vals = [float(x) for x in re.split(r"[,\s;]+", text) if x]
        elif len(parts) in (2, 3):
            x = [float(p) for p in parts]
            a, step, b = (x[0], 1.0, x[1]) if len(x) == 2 else x
            if step <= 0 or b < a or (b - a) / step > 200:
                return None
            n = int(math.floor((b - a) / step + 1e-9)) + 1
            vals = [a + k * step for k in range(n)]
        else:
            return None
    except ValueError:
        return None
    if not vals or not all(math.isfinite(v) for v in vals):
        return None
    return sorted({round(v, 2) for v in vals})


def format_depths(depths: list[float]) -> str:
    return ", ".join(f"{d:g}" for d in depths)


def _uses_device(model: ConfigModel, serial: str, ds: str) -> list[bool]:
    """Per transducer of the config: whether it is serial driven by ds."""
    out = []
    for ti in range(model.n_transducers):
        t_serial = model.get(f"{TX}.serial", ti)
        t_ds = model.get(f"{TX}.combo.ds_serial", ti)
        out.append(str(t_serial) == serial and (str(t_ds) if not is_empty(t_ds) else "") == ds)
    return out


def study_depths(model: ConfigModel, serial: str, ds: str, manifest=None) -> list[float]:
    """Focal depth settings (from the exit plane) this study uses with serial/ds: every
    transducer of the config with that device, and manifest columns
    transducer[(i)].focal_distance_ep."""
    out: set[float] = set()
    txs = model.effective_transducers()
    match = _uses_device(model, serial, ds)
    for ti, same in enumerate(match):
        if same and ti < len(txs):
            bowl = pl.Bowl.from_transducer(txs[ti])
            ep = pl.resolve_focal_distance_ep(txs[ti], bowl)
            if ep is not None:
                out.add(round(float(ep), 2))
    if manifest is not None:
        for col in manifest.columns:
            m = _MANIFEST_EP.match(col)
            if not m:
                continue
            ti = int(m.group(1)) - 1 if m.group(1) else 0
            if ti >= len(match) or not match[ti]:
                continue
            for row in manifest.rows:
                try:
                    out.add(round(float(row.get(col, "")), 2))
                except (TypeError, ValueError):
                    pass
    return sorted(out)


def target_depths(model: ConfigModel, serial: str, ds: str, depth_for) -> dict[int, float]:
    """Planned target depths (exit plane -> focus_pos, 0.1 mm) of the transducers that use
    serial/ds, by transducer index; depth_for(ti) gives the depth or None (no placement)."""
    out = {}
    for ti, same in enumerate(_uses_device(model, serial, ds)):
        d = depth_for(ti) if same else None
        if d is not None and math.isfinite(d) and d > 0:
            out[ti] = round(float(d), 1)
    return out


def desired_depths(settings: list[float], targets: dict[int, float]) -> list[float]:
    """The depths to calibrate for: the focal depth settings and the planned target depths."""
    return sorted({round(float(d), 2) for d in [*settings, *targets.values()]})


def describe_depths(settings: list[float], targets: dict[int, float]) -> str:
    """Where desired_depths() came from, for the page log."""
    parts = []
    if settings:
        parts.append(f"focal depth setting {format_depths(settings)} mm")
    for ti, d in targets.items():
        same = any(abs(s - d) < 0.05 for s in settings)
        parts.append(f"transducer {ti + 1} target {d:g} mm" + (" (= setting)" if same else ""))
    return "; ".join(parts)


@dataclass
class FitSettings:
    """The Calibration page's controls, in nexus_calibration.m's 'settings' form."""
    mode: str = "per_depth"
    free_params: list[str] = field(default_factory=lambda: ["phase", "amp"])
    forward_model: str = "rayleigh"
    near_field_skip_mm: float = 10.0
    domain: str = "intensity"
    huber_delta: float = 0.8
    n_starts: int = 16
    isppa_wcm2: float = 30.0
    align: str = "label"
    max_extrapolation_mm: float = 0.0
    calibrate_at: str = "measured"
    advanced: dict[str, Any] = field(default_factory=dict)

    def bridge_settings(self) -> dict:
        return {"forward_model": self.forward_model, "near_field_skip_mm": self.near_field_skip_mm,
                "huber_delta": self.huber_delta, "n_starts": self.n_starts, "domain": self.domain,
                "isppa_wcm2": self.isppa_wcm2, "align": self.align,
                "max_extrapolation_mm": self.max_extrapolation_mm, "advanced": dict(self.advanced)}


def _same(a: Any, b: Any) -> bool:
    if isinstance(a, (list, tuple)) or isinstance(b, (list, tuple)):
        a, b = list(a) if isinstance(a, (list, tuple)) else [a], list(b) if isinstance(b, (list, tuple)) else [b]
        return len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b))
    if isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        return math.isclose(float(a), float(b), rel_tol=1e-9, abs_tol=1e-12)
    return a == b


def settings_from_record(rec: dict | None, current: FitSettings) -> tuple[FitSettings, list[float] | None]:
    """Settings a saved entry was fitted with (meta.fit), on top of the current ones, and the
    desired depths when it was calibrated at desired depths -- as apply_fit_record in
    calibration_gui.m."""
    if not isinstance(rec, dict):
        return current, None
    fo = rec.get("fit_options") or {}
    s = FitSettings(**{k: (list(v) if isinstance(v, list) else dict(v) if isinstance(v, dict) else v)
                       for k, v in current.__dict__.items()})
    s.mode = "static" if rec.get("mode") == "static" else "per_depth"
    free = fo.get("FreeParams", ["phase", "amp"])
    s.free_params = [free] if isinstance(free, str) else list(free)
    s.forward_model = str(fo.get("ForwardModel", s.forward_model))
    s.domain = str(fo.get("Domain", s.domain))
    for key, attr in (("NearFieldSkipMm", "near_field_skip_mm"), ("HuberDelta", "huber_delta")):
        if key in fo:
            setattr(s, attr, float(fo[key]))
    if "NStarts" in fo:
        s.n_starts = int(fo["NStarts"])
    if rec.get("isppa_wcm2") is not None:
        s.isppa_wcm2 = float(rec["isppa_wcm2"])
    po = rec.get("profile_options") or {}
    if isinstance(po, dict):
        s.align = str(po.get("Align", s.align))
        if "MaxExtrapolationMm" in po:
            s.max_extrapolation_mm = float(po["MaxExtrapolationMm"])
    s.calibrate_at = "desired" if rec.get("calibrate_at") == "desired" else "measured"
    s.advanced = {k: fo[k] for k, d in ADVANCED_DEFAULTS.items() if k in fo and not _same(fo[k], d)}
    depths = None
    if s.calibrate_at == "desired":
        d = rec.get("depths_ep_mm")
        depths = [float(x) for x in (d if isinstance(d, list) else [d])]
    return s, depths
