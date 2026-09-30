"""Transducer equipment library (config/equipment) and calibration coverage.

Mirrors load_equipment_config.m and the file lookup in
resolve_transducer_from_serial.m: calibration entries are
<serial>_<ds_serial>.yaml, falling back to <serial>.yaml, in
transducer.library_path (if it holds an entry) or config/calibration.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from . import yamlio


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
class Calibration:
    """Calibration library coverage for a transducer / driving-system pair."""
    path: Path | None                   # file PRESTUS would use, or None
    key: str                            # library key tried first (<serial>_<ds> or <serial>)
    fallback: bool = False              # DS-agnostic <serial>.yaml used
    calibration_type: str = ""
    intensities: list[float] = field(default_factory=list)   # W/cm² with amplitude scaling
    depths_mm: list[float] = field(default_factory=list)     # calibrated focal depths (ep)
    updated: str = ""

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


class EquipmentLibrary:
    def __init__(self, prestus_root: str | Path):
        self.root = Path(prestus_root)
        self.equipment_dir = self.root / "config" / "equipment"
        self.calibration_dir = self.root / "config" / "calibration"
        self.devices: dict[str, Device] = {}
        self._load()

    def _load(self) -> None:
        if not self.equipment_dir.is_dir():
            return
        for path in sorted(self.equipment_dir.glob("*.yaml")):
            if path.name == "equipment_info.yaml":
                continue
            try:
                data = yamlio.load_file(path) or {}
            except Exception:
                continue
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
        cal.calibration_type = str(meta.get("calibration_type", ""))
        cal.updated = str(meta.get("last_updated", ""))
        for k in (model.get("amplitude_scaling") or {}):
            m = re.fullmatch(r"i(\d+(?:p\d+)?)", str(k))
            if m:
                cal.intensities.append(float(m.group(1).replace("p", ".")))
        depths = model.get("depths_ep_mm") or []
        cal.depths_mm = [float(d) for d in depths] if isinstance(depths, list) else [float(depths)]
        return cal
