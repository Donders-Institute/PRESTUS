"""Transducer entries (annular or matrix) and driving-system combinations for the equipment
library, as calibration_gui's 'New…' form writes them (on_new_transducer / collect_new in
functions/gui/calibration_gui.m), without Qt so the rules are testable.

A form gives the transducer fields and its combinations, each a driving-system serial,
optionally a new driving system (written as ds-<DS>_driver.yaml next to the transducer)
and an axial-profile CSV (annular only; stored relative to local/profiles when under it).
Entries are saved into your own library folders; an own entry with a shipped serial
replaces the shipped one (load_equipment_config.m). Maintainers may save into the shipped
config/equipment (calibration_gui's AllowShipped), with the shipped layout.

A matrix entry's element count comes from its element layout, which MATLAB computes
(convert_to_element_pos, through nexus_equipment 'layout'): build_entry takes it as
matrix_n_elem.
"""
from __future__ import annotations

import copy
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from . import yamlio

SERIAL_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")     # MATLAB isvarname, as calibration_gui
HEADER = "Written by the PRESTUS Nexus equipment editor (Transducer tab > Equipment)."


@dataclass
class NewDrivingSystem:
    serial: str
    name: str = ""
    manufact: str = ""
    available_ch: int = 4


@dataclass
class Combo:
    ds_serial: str
    csv: str = ""                                  # axial-profile CSV ('' = none yet; annular only)
    new_ds: NewDrivingSystem | None = None         # define this driving system (not in the library)
    keep: dict = field(default_factory=dict)       # further fields of an existing combo (drive_limit, ...)


@dataclass
class MatrixFields:
    steering: str = "3D"                           # '3D' | '1D'
    is_curved: bool = True
    outer_diameter_mm: float = 60.0
    elem_shape: str = "rect"                       # 'rect' | 'disc' | 'bowl'
    elem_width_mm: float = 3.0
    elem_height_mm: float = 3.0
    grid: str = "rect"                             # 'rect' | 'fibonacci' | 'fermat'
    n_row: int = 16
    n_col: int = 16
    spacing_height_mm: float = 0.5
    spacing_width_mm: float = 0.5
    sparsity: float = 1.0
    fib_n: int = 256
    fib_kerf_mm: float = 0.5
    fermat_n: int = 256


@dataclass
class TransducerForm:
    kind: str = "annular"                          # 'annular' | 'matrix'
    serial: str = ""
    name: str = ""
    manufact: str = ""
    model: str = ""
    freq_khz: float = 250.0
    curv_radius_mm: float = 63.0
    dist_geom_ep_mm: float = 52.0
    min_foc: float = 10.0
    max_foc: float = 60.0
    n_channels: int = 4                            # annular: physical channels
    elem_id_mm: str = ""                           # annular: one value per ring, comma- or space-separated
    elem_od_mm: str = ""
    matrix: MatrixFields = field(default_factory=MatrixFields)
    combos: list[Combo] = field(default_factory=list)


def parse_list(text: str) -> list[float] | None:
    """Numbers separated by commas, spaces, semicolons or brackets; None if any is not a number."""
    tokens = [t for t in re.split(r"[\[\],;\s]+", str(text).strip()) if t]
    try:
        return [float(t) for t in tokens]
    except ValueError:
        return None


def _fmt(values) -> str:
    return ", ".join(f"{float(v):g}" for v in values)


def profile_path(csv: str, profile_root: Path) -> str:
    """The CSV relative to local/profiles when it lies there (as the shipped entries store it)."""
    csv = csv.strip()
    if not csv:
        return ""
    try:
        return Path(csv).resolve().relative_to(Path(profile_root).resolve()).as_posix()
    except (ValueError, OSError):
        return csv


def matrix_block(form: TransducerForm) -> dict:
    """transducer.matrix of the form (the fields calibration_gui writes; config_default names)."""
    m = form.matrix
    if m.grid == "rect":
        gs = {"type": "rect", "rect": {"elem_n_row": int(m.n_row), "elem_n_col": int(m.n_col),
                                       "elem_spacing_height_mm": float(m.spacing_height_mm),
                                       "elem_spacing_width_mm": float(m.spacing_width_mm),
                                       "sparsity_factor": float(m.sparsity)}}
    elif m.grid == "fibonacci":
        gs = {"type": "fibonacci", "fibonacci": {"elem_n": int(m.fib_n), "kerf_mm": float(m.fib_kerf_mm)}}
    else:
        gs = {"type": "fermat", "fermat": {"elem_n": int(m.fermat_n)}}
    return {"steering": m.steering, "elem_shape": m.elem_shape, "elem_width_mm": float(m.elem_width_mm),
            "elem_height_mm": float(m.elem_height_mm), "outer_diameter_mm": float(m.outer_diameter_mm),
            "is_curved": bool(m.is_curved), "curv_radius_mm": float(form.curv_radius_mm),
            "dist_geom_ep_mm": float(form.dist_geom_ep_mm),
            "matrix_shape": {"type": "define_here", "define_here": {"grid_shape": gs}}}


def build_entry(form: TransducerForm, profile_root: Path, base: dict | None = None,
                matrix_n_elem: int | None = None) -> tuple[dict, list[str]]:
    """Transducer entry (dict for YAML) and the issues that block saving. base: the raw entry
    being edited (same serial), whose further fields (housing, status, corrections, combos
    not in the form) are kept. matrix_n_elem: element count of the matrix layout (MATLAB)."""
    issues: list[str] = []
    serial = form.serial.strip()
    if not SERIAL_RE.match(serial):
        issues.append("Serial must start with a letter and contain only letters, digits and underscores.")
    if form.min_foc >= form.max_foc:
        issues.append("Focal range: minimum must be below maximum.")
    entry: dict[str, Any] = copy.deepcopy(base) if base else {}
    entry["type"] = "transducer"
    entry["serial"] = serial
    if form.model.strip():
        entry["model"] = form.model.strip()
    entry["name"] = form.name.strip() or serial
    if form.manufact.strip():
        entry["manufact"] = form.manufact.strip()
    else:
        entry.pop("manufact", None)
    entry["min_foc"] = float(form.min_foc)
    entry["max_foc"] = float(form.max_foc)
    tr = entry.get("transducer") if isinstance(entry.get("transducer"), dict) else {}
    if form.kind == "matrix":
        tr.pop("annular", None)
        m = matrix_block(form)
        if m["is_curved"] and m["curv_radius_mm"] <= m["outer_diameter_mm"] / 2:
            issues.append("Curvature radius must exceed half the aperture diameter.")
        if form.matrix.grid == "rect" and (form.matrix.spacing_height_mm <= 0 or form.matrix.spacing_width_mm <= 0):
            issues.append("Element spacing must be greater than 0.")
        if matrix_n_elem is None:
            issues.append("Element layout not computed yet (MATLAB).")
        elif matrix_n_elem < 1:
            issues.append("No element of the grid lies inside the aperture.")
        entry["n_elem"] = int(matrix_n_elem or 0)
        tr["matrix"] = m
    else:
        tr.pop("matrix", None)
        idd, od = parse_list(form.elem_id_mm), parse_list(form.elem_od_mm)
        if not idd:
            issues.append("Ring inner Ø: enter one number per ring.")
        if not od:
            issues.append("Ring outer Ø: enter one number per ring.")
        if idd and od:
            if len(idd) != len(od):
                issues.append(f"{len(idd)} inner but {len(od)} outer diameters.")
            elif any(o <= i for i, o in zip(idd, od)):
                issues.append("Each outer diameter must exceed its inner diameter.")
            if form.n_channels < 1 or len(idd) % form.n_channels:
                issues.append(f"Ring count ({len(idd)}) must be a multiple of the physical channel count "
                              f"({form.n_channels}); rings are grouped into channels in order.")
        entry["n_elem"] = int(form.n_channels)
        tr["annular"] = {"elem_n": len(idd or []), "elem_id_mm": idd or [], "elem_od_mm": od or [],
                         "curv_radius_mm": float(form.curv_radius_mm),
                         "dist_geom_ep_mm": float(form.dist_geom_ep_mm)}
    tr["freq_hz"] = float(form.freq_khz) * 1e3
    entry["transducer"] = tr

    combos: dict[str, Any] = {}
    seen: set[str] = set()
    for c in form.combos:
        ds = c.ds_serial.strip()
        if not ds:
            if c.csv.strip():
                issues.append("A profile CSV belongs to a driving system: enter a DS serial.")
            continue
        if not SERIAL_RE.match(ds):
            issues.append(f"DS serial '{ds}' must start with a letter and contain only letters, "
                          "digits and underscores.")
        if ds in seen:
            issues.append(f"Driving system {ds} is listed twice.")
        seen.add(ds)
        if c.new_ds is not None and c.new_ds.available_ch < entry["n_elem"]:
            issues.append(f"Driving system {ds} has {c.new_ds.available_ch} channels, fewer than the "
                          f"transducer's {entry['n_elem']}.")
        item = dict(c.keep)
        if c.new_ds is not None:
            item.pop("ds", None)                   # taken from its ds-<DS>_driver.yaml
        csv = profile_path(c.csv, profile_root) if form.kind == "annular" else ""
        if csv:
            item["char_data_path"] = csv
        elif form.kind == "annular":
            item.pop("char_data_path", None)
        combos[ds] = item
    if combos:
        entry["combos"] = combos
    else:
        entry.pop("combos", None)
    return entry, issues


def ds_entry(ds: NewDrivingSystem) -> dict:
    out = {"type": "driving_system", "serial": ds.serial.strip(), "name": ds.name.strip() or ds.serial.strip()}
    if ds.manufact.strip():
        out["manufact"] = ds.manufact.strip()
    out["available_ch"] = int(ds.available_ch)
    return out


def form_from_device(dev, raw: dict | None, keep_serial: bool) -> tuple[TransducerForm, dict | None]:
    """Form filled from a library device. keep_serial: edit that entry (its combinations
    listed, its other fields kept); else copy the geometry under a new serial, as
    calibration_gui's 'Copy from'. Returns the form and the base entry for build_entry."""
    geo = dev.geometry or {}
    form = TransducerForm(
        serial=dev.serial if keep_serial else "", name=dev.name if keep_serial else "",
        manufact=dev.manufact or "", model=(dev.model or "") if keep_serial else "",
        freq_khz=(dev.freq_hz or 250e3) / 1e3,
        min_foc=float(dev.min_foc if dev.min_foc is not None else 10.0),
        max_foc=float(dev.max_foc if dev.max_foc is not None else 60.0))
    if geo.get("matrix"):
        m = geo["matrix"]
        form.kind = "matrix"
        form.curv_radius_mm = float(m.get("curv_radius_mm") or form.curv_radius_mm)
        form.dist_geom_ep_mm = float(m.get("dist_geom_ep_mm") or form.dist_geom_ep_mm)
        mf = form.matrix
        mf.steering = str(m.get("steering") or mf.steering)
        mf.is_curved = bool(m.get("is_curved", True))
        mf.outer_diameter_mm = float(m.get("outer_diameter_mm") or mf.outer_diameter_mm)
        mf.elem_shape = str(m.get("elem_shape") or mf.elem_shape)
        mf.elem_width_mm = float(m.get("elem_width_mm") or mf.elem_width_mm)
        mf.elem_height_mm = float(m.get("elem_height_mm") or mf.elem_height_mm)
        gs = (((m.get("matrix_shape") or {}).get("define_here") or {}).get("grid_shape") or {})
        mf.grid = str(gs.get("type") or mf.grid)
        g = gs.get(mf.grid) or {}
        if mf.grid == "rect":
            mf.n_row, mf.n_col = int(g.get("elem_n_row", mf.n_row)), int(g.get("elem_n_col", mf.n_col))
            mf.spacing_height_mm = float(g.get("elem_spacing_height_mm", mf.spacing_height_mm))
            mf.spacing_width_mm = float(g.get("elem_spacing_width_mm", mf.spacing_width_mm))
            mf.sparsity = float(g.get("sparsity_factor") or 1.0)
        elif mf.grid == "fibonacci":
            mf.fib_n, mf.fib_kerf_mm = int(g.get("elem_n", mf.fib_n)), float(g.get("kerf_mm") or 0.0)
        elif mf.grid == "fermat":
            mf.fermat_n = int(g.get("elem_n", mf.fermat_n))
    else:
        a = geo.get("annular") or {}
        form.curv_radius_mm = float(a.get("curv_radius_mm", 63.0))
        form.dist_geom_ep_mm = float(a.get("dist_geom_ep_mm", 52.0))
        form.n_channels = int(dev.n_elem or a.get("elem_n") or 1)
        form.elem_id_mm = _fmt(_as_list(a.get("elem_id_mm")))
        form.elem_od_mm = _fmt(_as_list(a.get("elem_od_mm")))
    if not keep_serial:
        return form, None
    base = copy.deepcopy(raw) if isinstance(raw, dict) else {}
    for ds, item in ((base.get("combos") or {}).items()):
        item = item if isinstance(item, dict) else {}
        form.combos.append(Combo(ds_serial=str(ds), csv=str(item.get("char_data_path", "")),
                                 keep={k: v for k, v in item.items() if k != "char_data_path"}))
    return form, base


def _as_list(v) -> list:
    if v is None:
        return []
    return list(v) if isinstance(v, (list, tuple)) else [v]


def target_files(entry: dict, folder: Path, dev, new_ds: list[NewDrivingSystem],
                 shipped: Path | None = None) -> dict[Path, dict]:
    """Files to write. Own library (shipped None): the transducer entry into its own file when
    it already lives in one of your libraries, else <folder>/tran-<serial>_transducer.yaml;
    each new driving system as <folder>/ds-<DS>_driver.yaml. Shipped library (maintainers,
    shipped = config/equipment): into the entry's existing file, else
    model-<MODEL>/tran-<serial>_transducer.yaml when that model folder exists, else
    <serial>.yaml (as calibration_gui); new driving systems as ds/ds-<DS>_driver.yaml."""
    if shipped is not None:
        shipped = Path(shipped)
        model = str(entry.get("model") or "")
        if dev is not None and not getattr(dev, "user", False) and getattr(dev, "file", ""):
            tfile = Path(dev.file)
        elif model and (shipped / f"model-{model}").is_dir():
            tfile = shipped / f"model-{model}" / f"tran-{entry['serial']}_transducer.yaml"
        else:
            tfile = shipped / f"{entry['serial']}.yaml"
        out = {tfile: entry}
        for ds in new_ds:
            out[shipped / "ds" / f"ds-{ds.serial.strip()}_driver.yaml"] = ds_entry(ds)
        return out
    folder = Path(folder)
    if dev is not None and getattr(dev, "user", False) and getattr(dev, "file", ""):
        tfile = Path(dev.file)
    else:
        tfile = folder / f"tran-{entry['serial']}_transducer.yaml"
    out = {tfile: entry}
    for ds in new_ds:
        out[folder / f"ds-{ds.serial.strip()}_driver.yaml"] = ds_entry(ds)
    return out


def write(files: dict[Path, dict]) -> list[Path]:
    written = []
    for path, data in files.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        yamlio.dump_file(path, data, header=HEADER)
        written.append(path)
    return written
