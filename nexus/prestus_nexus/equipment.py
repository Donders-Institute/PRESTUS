"""Transducer equipment library (config/equipment) and calibration coverage.

Mirrors load_equipment_config.m and the file lookup in
resolve_transducer_from_serial.m: calibration entries are
<serial>_<ds_serial>.yaml, falling back to <serial>.yaml, in
transducer.library_path (if it holds an entry) or config/calibration.
"""

from __future__ import annotations

import hashlib
import os
import pickle
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from . import yamlio
from .localcopy import cache_root


@dataclass
class DrivingSystem:
    serial: str
    name: str = ""
    manufact: str = ""
    available_ch: int | None = None


@dataclass
class Device:
    serial: str
    name: str
    manufact: str
    n_elem: int | None
    min_foc: float | None
    max_foc: float | None
    geometry: dict                      # the 'transducer' block: {annular|matrix: {...}, freq_hz}
    combos: dict[str, DrivingSystem] = field(default_factory=dict)

    @property
    def type(self) -> str | None:
        for kind in ("annular", "matrix"):
            if self.geometry.get(kind):
                return kind
        return None

    @property
    def freq_hz(self) -> float | None:
        return self.geometry.get("freq_hz")

    def label(self) -> str:
        return f"{self.serial} — {self.name}" if self.name and self.name != self.serial else self.serial


@dataclass
class DepthStatus:
    """Calibration.at_depth() result: where a focal_distance_ep falls relative to the
    library entry's calibrated (refined) depths, mirroring the blend rule
    library_model_at_depth.m applies at run time -- pure arithmetic, no forward model."""
    depth_mm: float
    in_range: bool                      # within [min(depths_mm), max(depths_mm)]
    bracket: tuple[float, float] | None = None       # the two depths either side, if in range
    nearest_mm: float | None = None                  # nearest calibrated depth, if out of range
    extrapolation_mm: float = 0.0                     # distance beyond the calibrated range (0 if in range)
    nearest_r2: float | None = None                   # fit_r2 at the nearest fit_depths_ep_mm entry
    phase_model: str = "static"                        # 'refined' | 'blend' | 'static'


@dataclass
class Calibration:
    """Calibration library coverage for a transducer / driving-system pair."""
    path: Path | None                   # file PRESTUS would use, or None
    key: str                            # library key tried first (<serial>_<ds> or <serial>)
    fallback: bool = False              # DS-agnostic <serial>.yaml used
    calibration_type: str = ""
    intensities: list[float] = field(default_factory=list)   # W/cm² with amplitude scaling
    depths_mm: list[float] = field(default_factory=list)     # calibrated (refined) focal depths (ep)
    updated: str = ""
    blend_mm: float = 5.0                               # refinement.blend_mm (run-time default)
    fit_depths_mm: list[float] = field(default_factory=list)   # global_model.fit_depths_ep_mm
    fit_r2: list[float] = field(default_factory=list)          # paired with fit_depths_mm
    amplitude_depths_mm: list[float] = field(default_factory=list)   # amplitude_by_depth.depths_ep_mm
    amplitude_isppa_wcm2: float | None = None

    @property
    def found(self) -> bool:
        return self.path is not None

    def summary(self) -> str:
        if not self.found:
            return f"no calibration entry ({self.key}.yaml)"
        parts = [self.path.name]
        if self.calibration_type:
            parts.append(self.calibration_type)
        if self.depths_mm:
            parts.append(f"depths {min(self.depths_mm):g}–{max(self.depths_mm):g} mm")
        if self.intensities:
            parts.append("amplitude at " + ", ".join(f"{i:g}" for i in self.intensities) + " W/cm²")
        if self.updated:
            parts.append(f"updated {self.updated}")
        return " · ".join(parts)

    def at_depth(self, focal_distance_ep_mm: float) -> DepthStatus | None:
        """Coverage of this entry at a given depth: in/out of the calibrated range, the
        bracketing or nearest calibrated depth(s), nearby fit quality, and the qualitative
        phase_model state library_model_at_depth.m would resolve (static/blend/refined,
        by distance to the nearest refined depth vs. blend_mm) -- no physics is re-run."""
        if not self.found or not self.depths_mm:
            return None
        d = float(focal_distance_ep_mm)
        depths = sorted(self.depths_mm)
        lo, hi = depths[0], depths[-1]
        in_range = lo <= d <= hi
        bracket = nearest = None
        extrap = 0.0
        if in_range:
            below = [x for x in depths if x <= d]
            above = [x for x in depths if x >= d]
            bracket = (below[-1], above[0])
        else:
            nearest = lo if d < lo else hi
            extrap = lo - d if d < lo else d - hi
        nearest_fit_r2 = None
        if self.fit_depths_mm and self.fit_r2:
            j = min(range(len(self.fit_depths_mm)), key=lambda k: abs(self.fit_depths_mm[k] - d))
            nearest_fit_r2 = self.fit_r2[j]
        # library_model_at_depth.m: weight w=1 (fully 'refined') throughout the interpolated
        # span; outside it, w fades linearly to 0 ('static') over blend_mm past the nearest
        # endpoint -- not a function of distance to the *nearest* refined depth in general.
        blend = self.blend_mm or 5.0
        phase_model = "refined" if in_range else ("blend" if extrap <= blend else "static")
        return DepthStatus(depth_mm=d, in_range=in_range, bracket=bracket, nearest_mm=nearest,
                           extrapolation_mm=extrap, nearest_r2=nearest_fit_r2, phase_model=phase_model)


class EquipmentLibrary:
    def __init__(self, prestus_root: str | Path):
        self.root = Path(prestus_root)
        self.equipment_dir = self.root / "config" / "equipment"
        self.calibration_dir = self.root / "config" / "calibration"
        self.devices: dict[str, Device] = {}
        self._load()

    def _load(self) -> None:
        for name, data in sorted(self._read_entries().items()):
            if not isinstance(data, dict) or data.get("type") != "transducer" or "serial" not in data:
                continue
            combos = {}
            for ds_key, entry in (data.get("combos") or {}).items():
                ds = (entry or {}).get("ds", {}) or {}
                combos[str(ds_key)] = DrivingSystem(str(ds.get("serial", ds_key)), ds.get("name", ""),
                                                    ds.get("manufact", ""), ds.get("available_ch"))
            self.devices[str(data["serial"])] = Device(
                serial=str(data["serial"]), name=str(data.get("name", "")),
                manufact=str(data.get("manufact", "")), n_elem=data.get("n_elem"),
                min_foc=data.get("min_foc"), max_foc=data.get("max_foc"),
                geometry=data.get("transducer", {}) or {}, combos=combos)

    def _read_entries(self) -> dict[str, Any]:
        """file name -> parsed YAML of config/equipment/*.yaml. The parsed files are kept in
        the user's cache folder and reused while a file's size and modification time are
        unchanged: on a network drive, reading the ~20 YAML files took ~8 s per start."""
        try:
            with os.scandir(self.equipment_dir) as it:
                stamps = {e.name: (e.stat().st_size, e.stat().st_mtime_ns) for e in it
                          if e.is_file() and e.name.endswith(".yaml") and e.name != "equipment_info.yaml"}
        except OSError:
            return {}
        cache = cache_root() / ("equipment-" + hashlib.sha1(
            str(self.equipment_dir.resolve()).encode()).hexdigest()[:10] + ".pkl")
        try:
            old = pickle.loads(cache.read_bytes())
        except Exception:
            old = {}
        entries, fresh = {}, {}
        for name, stamp in stamps.items():
            hit = old.get(name)
            if hit is not None and hit[0] == stamp:
                fresh[name] = hit
            else:
                try:
                    fresh[name] = (stamp, yamlio.load_file(self.equipment_dir / name) or {})
                except Exception:
                    continue
            entries[name] = fresh[name][1]
        if fresh != old:
            try:
                cache.parent.mkdir(parents=True, exist_ok=True)
                cache.write_bytes(pickle.dumps(fresh))
            except Exception:
                pass
        return entries

    def get(self, serial: Any) -> Device | None:
        return self.devices.get(str(serial)) if serial else None

    def calibration(self, serial: str, ds_serial: str | None = None,
                    library_path: str | None = None) -> Calibration:
        key = f"{serial}_{ds_serial}" if ds_serial else serial
        folder = self.calibration_dir
        if library_path:
            user = Path(library_path)
            if (user / f"{key}.yaml").is_file() or (user / f"{serial}.yaml").is_file():
                folder = user
        for name, fallback in ((key, False), (serial, key != serial)):
            path = folder / f"{name}.yaml"
            if path.is_file():
                return self._describe(path, key, fallback)
        return Calibration(path=None, key=key)

    @staticmethod
    def _describe(path: Path, key: str, fallback: bool) -> Calibration:
        cal = Calibration(path=path, key=key, fallback=fallback)
        try:
            data = yamlio.load_file(path) or {}
        except Exception:
            return cal
        meta = data.get("meta", {}) or {}
        model = data.get("global_model", {}) or {}
        refinement = data.get("refinement", {}) or {}
        cal.calibration_type = str(meta.get("calibration_type", ""))
        cal.updated = str(meta.get("last_updated", ""))
        for k in (model.get("amplitude_scaling") or {}):
            m = re.fullmatch(r"i(\d+(?:p\d+)?)", str(k))
            if m:
                cal.intensities.append(float(m.group(1).replace("p", ".")))
        # The calibrated/refined depths: refinement.depths_ep_mm in the current schema
        # (library_model_at_depth.m), global_model.depths_ep_mm for the legacy direct form
        # (its docstring), else no refinement -- fit_depths_mm below is still shown.
        depths = refinement.get("depths_ep_mm") or model.get("depths_ep_mm") or []
        cal.depths_mm = _floats(depths)
        cal.blend_mm = _float(refinement.get("blend_mm"), 5.0)
        cal.fit_depths_mm = _floats(model.get("fit_depths_ep_mm"))
        cal.fit_r2 = _floats(model.get("fit_r2"))
        amp_by_depth = model.get("amplitude_by_depth") or {}
        cal.amplitude_depths_mm = _floats(amp_by_depth.get("depths_ep_mm"))
        isppa = amp_by_depth.get("isppa_wcm2")
        cal.amplitude_isppa_wcm2 = _float(isppa, None)
        return cal


def _floats(v: Any) -> list[float]:
    if v is None:
        return []
    return [float(x) for x in v] if isinstance(v, list) else [float(v)]


def _float(v: Any, default: float | None) -> float | None:
    return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else default
