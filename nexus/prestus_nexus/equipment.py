"""Transducer equipment library (config/equipment) and calibration coverage.

Mirrors load_equipment_config.m (entries in subfolders, units taking what they omit
from their model, driving systems in ds/) and the lookup in
resolve_transducer_from_serial.m: the entry for <serial>_<ds_serial>, falling back to
<serial>, flat (<key>.yaml) or nested (model-*/tran-*_calib.yaml, library_entry_file.m),
in transducer.library_path (if it holds an entry) or config/calibration.
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
    model: str = ""                     # serial of the device model (config/equipment/README.md)
    housing: dict = field(default_factory=dict)   # display only: diameter_mm, depth_mm
    status: str = ""                    # 'deprecated' | 'excluded': kept for earlier studies, not for new ones
    superseded_by: str = ""
    user: bool = False                  # from one of your own libraries
    file: str = ""                      # the YAML the entry came from

    @property
    def type(self) -> str | None:
        for kind in ("annular", "matrix"):
            if self.geometry.get(kind):
                return kind
        return None

    @property
    def freq_hz(self) -> float | None:
        return self.geometry.get("freq_hz")

    @property
    def deprecated(self) -> bool:
        return self.status == "deprecated"

    @property
    def excluded(self) -> bool:
        """Withheld from the library offering (not listed); still resolves for earlier studies."""
        return self.status == "excluded"

    def label(self) -> str:
        name = self.name.replace(self.serial, "").strip()     # equipment names often end in the serial
        text = f"{self.serial} — {name}" if name else self.serial
        if self.user:
            text += " (own library)"
        return f"{text} ({self.status})" if self.deprecated or self.excluded else text


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
        self.user_folders: tuple[Path, ...] = ()
        self._shipped_entries: dict[str, Any] | None = None
        self._load()

    def set_user_folders(self, folders) -> None:
        """Also read the equipment entries (*_transducer.yaml, *_driver.yaml) in these
        folders, the study's transducer.library_path folders, as load_equipment_config.m
        does; they add serials or replace shipped ones. A no-op when unchanged."""
        folders = tuple(Path(f) for f in folders or ())
        if folders == self.user_folders:
            return
        self.user_folders = folders
        self.devices = {}
        self._load()

    def _user_entries(self) -> dict[str, Any]:
        """Own equipment entries, keyed so that sorting puts the most specific folder (the
        first) last: it replaces entries of the same serial from later folders."""
        entries: dict[str, Any] = {}
        n = len(self.user_folders)
        for rank, folder in enumerate(self.user_folders):
            try:
                if not folder.is_dir() or folder.resolve() == self.equipment_dir.resolve():
                    continue
                files = sorted([*folder.rglob("*_transducer.yaml"), *folder.rglob("*_driver.yaml")])
            except OSError:
                continue
            for f in files:
                try:
                    entries[f"user:{n - 1 - rank:03d}:{f}"] = yamlio.load_file_cached(f) or {}
                except Exception:
                    continue
        return entries

    def _load(self) -> None:
        if self._shipped_entries is None:
            self._shipped_entries = self._read_entries()
        entries = {**self._shipped_entries, **self._user_entries()}
        drivers = {str(d["serial"]): d for d in entries.values()
                   if isinstance(d, dict) and d.get("type") == "driving_system" and "serial" in d}
        for name, data in sorted(entries.items()):
            if not isinstance(data, dict) or data.get("type") != "transducer" or "serial" not in data:
                continue
            combos = {}
            for ds_key, entry in (data.get("combos") or {}).items():
                ds = (entry or {}).get("ds") or drivers.get(str(ds_key)) or {}
                combos[str(ds_key)] = DrivingSystem(str(ds.get("serial", ds_key)), ds.get("name", ""),
                                                    ds.get("manufact", ""), ds.get("available_ch"))
            self.devices[str(data["serial"])] = Device(
                serial=str(data["serial"]), name=str(data.get("name", "")),
                manufact=str(data.get("manufact", "")), n_elem=data.get("n_elem"),
                min_foc=data.get("min_foc"), max_foc=data.get("max_foc"),
                geometry=data.get("transducer", {}) or {}, combos=combos,
                model=str(data.get("model") or ""),
                housing=data.get("housing") if isinstance(data.get("housing"), dict) else {},
                status=str(data.get("status") or ""), superseded_by=str(data.get("superseded_by") or ""),
                user=name.startswith("user:"),
                file=name.split(":", 2)[2] if name.startswith("user:") else str(self.equipment_dir / name))
        # A unit takes what it omits from its model entry, as load_equipment_config.m
        for dev in self.devices.values():
            base = self.devices.get(dev.model) if dev.model != dev.serial else None
            if base is None:
                continue
            dev.geometry = dev.geometry or dict(base.geometry)
            dev.housing = dev.housing or dict(base.housing)
            dev.min_foc = base.min_foc if dev.min_foc is None else dev.min_foc
            dev.max_foc = base.max_foc if dev.max_foc is None else dev.max_foc
            dev.n_elem = base.n_elem if dev.n_elem is None else dev.n_elem
            dev.manufact = dev.manufact or base.manufact

    def _read_entries(self) -> dict[str, Any]:
        """relative path -> parsed YAML of every config/equipment YAML, subfolders included
        (model-<MODEL>/tran-<SERIAL>_transducer.yaml, ds/ds-<DS>_driver.yaml, or a flat
        <SERIAL>.yaml). The parsed files are kept in the user's cache folder and reused while
        a file's size and modification time are unchanged: on a network drive, reading the
        ~30 YAML files took ~8 s per start."""
        stamps = {}
        try:
            for dirpath, dirnames, filenames in os.walk(self.equipment_dir):
                dirnames[:] = [d for d in dirnames if not d.startswith(".")]
                for fn in filenames:
                    if fn.endswith(".yaml") and fn != "equipment_info.yaml":
                        full = Path(dirpath) / fn
                        st = full.stat()
                        stamps[full.relative_to(self.equipment_dir).as_posix()] = (st.st_size, st.st_mtime_ns)
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

    def model_for(self, serial: Any) -> str | None:
        """The model a unit serial '<MODEL>_<unit id>[_<variant>]' that is not in the library
        falls back to, as resolve_transducer_from_serial.m: the longest '<MODEL>_' prefix
        whose entry is a model entry (model: its own serial); None if there is none."""
        serial = str(serial) if serial else ""
        if not serial or serial in self.devices:
            return None
        models = [s for s, d in self.devices.items() if d.model == s and serial.startswith(s + "_")]
        return max(models, key=len) if models else None

    def lookup(self, serial: Any) -> Device | None:
        """The device PRESTUS uses for a serial: its own entry, else its model's."""
        model = self.model_for(serial)
        return self.devices.get(model) if model else self.get(serial)

    def housing(self, serial: Any) -> dict:
        """The housing block (display only) of a device: its entry's 'housing', else its
        model's (a unit entry usually has none); {} if neither has one."""
        dev = self.lookup(serial) if serial else None
        if dev is None:
            return {}
        h = dev.housing
        if not h and dev.model and self.get(dev.model) is not None:
            h = self.get(dev.model).housing
        return h if isinstance(h, dict) else {}

    def calibration(self, serial: str, ds_serial: str | None = None,
                    library_path: str | None = None) -> Calibration:
        key = f"{serial}_{ds_serial}" if ds_serial else serial
        # Your own libraries first (this study's library_path, then the folders set by
        # set_user_folders: study, project, PRESTUS_TRANSDUCER_LIBRARY), as
        # resolve_transducer_from_serial.m: the combo in any of them, else the serial alone
        own = ([Path(library_path)] if library_path else []) + list(self.user_folders)
        for name, fallback in ((key, False), (serial, key != serial)):
            for folder in own:
                path = library_entry_file(folder, name)
                if path is not None:
                    return self._describe(path, key, fallback)
        for name, fallback in ((key, False), (serial, key != serial)):
            path = library_entry_file(self.calibration_dir, name)
            if path is not None:
                return self._describe(path, key, fallback)
        return Calibration(path=None, key=key)

    def profile_file(self, rel: str | None, serial: str | None = None) -> Path | None:
        """A char_data_path resolved as resolve_profile_file.m: absolute as is, else in
        profiles/ next to an own equipment entry, then each own library's profiles/, then
        <PRESTUS>/local/profiles."""
        if not rel:
            return None
        p = Path(str(rel))
        if p.is_absolute():
            return p if p.is_file() else None
        cands = []
        dev = self.devices.get(str(serial)) if serial else None
        if dev is not None and dev.user and dev.file:
            cands += [Path(dev.file).parent / "profiles" / p, Path(dev.file).parent / p]
        cands += [f / "profiles" / p for f in self.user_folders]
        cands.append(self.root / "local" / "profiles" / p)
        return next((c for c in cands if c.is_file()), None)

    @staticmethod
    def _describe(path: Path, key: str, fallback: bool) -> Calibration:
        cal = Calibration(path=path, key=key, fallback=fallback)
        try:
            data = yamlio.load_file_cached(path) or {}
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


def library_entry_file(folder: str | Path, key: str) -> Path | None:
    """The calibration-library entry for key, as library_entry_file.m: flat <KEY>.yaml, else
    nested model-<MODEL>/tran-<TRAN>[_ds-<DS>]_calib.yaml (key = <TRAN>[_<DS>]); None if absent."""
    folder = Path(folder)
    flat = folder / f"{key}.yaml"
    if flat.is_file():
        return flat
    try:
        nested = sorted(folder.glob("model-*/tran-*_calib.yaml"))
    except OSError:
        return None
    for path in nested:
        stem = path.name[len("tran-"):-len("_calib.yaml")]
        tran, _, ds = stem.partition("_ds-")
        if (f"{tran}_{ds}" if ds else tran) == key:
            return path
    return None


def _floats(v: Any) -> list[float]:
    if v is None:
        return []
    return [float(x) for x in v] if isinstance(v, list) else [float(v)]


def _float(v: Any, default: float | None) -> float | None:
    return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else default
