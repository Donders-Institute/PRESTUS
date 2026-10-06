"""Finds the runs whose results sit in one subject's output folder, for the Results tab.

One folder (<path.sim>/sub-XXX) usually holds several runs, told apart by their output
affix: a study's conditions, an uncertainty variant, the per-transducer sub-runs and the
joint result of async multi-transducer mode, the steps of a sequential run, ... Each run
writes ``<tag>_<medium><affix>.csv``, ``<tag>_<medium>_report<affix>.html`` and NIfTIs
``nii/<tag>_<medium>_T1w<affix>_<map>.nii.gz`` (nifti_acoustic.m / nifti_thermal.m), and
the placement stage ``<tag>_placement<affix>.json`` (write_placement_positions.m).
Read-only and Qt-free.
"""

from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass, field
from pathlib import Path

from .model import ConfigModel
from .outputs import _subject_id, dir_output
from .predict_outputs import predict
from .study import sanitize_affix

# Overlay maps: label -> NIfTI map names to try, newest name first.
MAPS: dict[str, tuple[str, ...]] = {
    "Intensity (ISPPA)": ("Itp", "intensity"),
    "Pulse-average intensity": ("Ipa",),
    "Pressure (peak)": ("pressure",),
    "Mechanical index": ("MI",),
    "Focal volume (-6 dB)": ("focal_mask_6dB",),
    "Max temperature": ("heating",),
    "Max temperature rise": ("heatrise",),
    "Temperature at end": ("heating_end",),
    "CEM43": ("CEM43_end", "CEM43"),
}
THERMAL_MAPS = ("Max temperature", "Max temperature rise", "Temperature at end", "CEM43")

# CSV columns shown when runs are compared side by side.
COMPARE_COLUMNS = (
    ("Isppa_brain", "ISPPA brain [W/cm²]"),
    ("Isppa", "ISPPA peak [W/cm²]"),
    ("MI_brain", "MI brain"),
    ("real_focal_distance_mm", "Focal distance [mm]"),
    ("halfmax_ISPPA_volume_brain_mm3", "Half-max volume [mm³]"),
    ("maxT_brain", "Max T brain [°C]"),
    ("riseT_brain", "ΔT brain [°C]"),
    ("riseT_skull", "ΔT skull [°C]"),
    ("CEM43_brain", "CEM43 brain [min]"),
)

_TRANS = re.compile(r"_trans-(\d+)$")
_JOINT = re.compile(r"_desc-asyncJoint(?:S(\d+))?$")


@dataclass
class Run:
    """The outputs of one run (one output affix) in a subject folder."""
    affix: str
    medium: str
    base: Path
    tag: str
    label: str = ""
    kind: str = "run"                 # 'run' | 'joint' | 'transducer'
    transducer: int | None = None     # 1-based, for kind == 'transducer'
    csv: Path | None = None
    report: Path | None = None
    maps: dict[str, Path] = field(default_factory=dict)

    @property
    def stem(self) -> str:
        return f"{self.tag}_{self.medium}{self.affix}"

    @property
    def heat(self) -> Path:
        """Thermal time courses of the run (thermal_plot_sim.m), read by read_heat."""
        return self.base / "cache" / f"{self.tag}_{self.medium}_HEAT{self.affix}.mat"


def describe_affix(affix: str) -> tuple[str, str, int | None]:
    """(label, kind, transducer) for an output affix, e.g. '_desc-async_trans-02' ->
    ('async · transducer 2 alone', 'transducer', 2)."""
    kind, ti = "run", None
    rest = affix
    m = _TRANS.search(rest)
    if m:
        kind, ti = "transducer", int(m.group(1))
        rest = rest[:m.start()]
    j = _JOINT.search(rest)
    if j:
        kind = "joint"
        rest = rest[:j.start()]
    parts = [p for p in re.split(r"_desc-|_", rest) if p]
    label = " · ".join(parts) if parts else "(no affix)"
    if kind == "transducer":
        label += f" · transducer {ti} alone"
    elif kind == "joint":
        label += " · both transducers (async)" + (f", sweep {int(j.group(1))}" if j.group(1) else "")
    return label, kind, ti


def _maps_for(base: Path, tag: str, medium: str, affix: str) -> dict[str, Path]:
    nii = base / "nii"
    found: dict[str, Path] = {}
    for label, names in MAPS.items():
        for name in names:
            p = nii / f"{tag}_{medium}_T1w{affix}_{name}.nii.gz"
            if p.is_file():
                found[label] = p
                break
    return found


def discover_runs(base: Path, tag: str, media: tuple[str, ...] = ("layered", "water", "phantom")) -> list[Run]:
    """Every run with a results CSV or a result NIfTI in ``base``, sorted by affix."""
    affixes: dict[tuple[str, str], None] = {}
    for medium in media:
        prefix = f"{tag}_{medium}"
        for p in base.glob(f"{prefix}*.csv"):
            affix = p.stem[len(prefix):]
            if affix and not affix.startswith("_"):
                continue
            if affix.endswith("_targeting"):
                continue
            affixes[(medium, affix)] = None
        nii = base / "nii"
        if nii.is_dir():
            t1w = f"{prefix}_T1w"
            for p in nii.glob(f"{t1w}*.nii.gz"):
                rest = p.name[len(t1w):-len(".nii.gz")]
                for names in MAPS.values():
                    hit = next((n for n in names if rest.endswith("_" + n)), None)
                    if hit:
                        affixes[(medium, rest[:-len(hit) - 1])] = None
                        break
    runs = []
    for medium, affix in sorted(affixes, key=lambda k: (k[0] != "layered", k[1])):
        label, kind, ti = describe_affix(affix)
        csv_path = base / f"{tag}_{medium}{affix}.csv"
        report = base / f"{tag}_{medium}_report{affix}.html"
        maps = _maps_for(base, tag, medium, affix)
        if not csv_path.is_file() and not maps:
            continue
        if medium != "layered":
            label = f"{medium} · {label}"
        runs.append(Run(affix=affix, medium=medium, base=base, tag=tag, label=label, kind=kind,
                        transducer=ti, csv=csv_path if csv_path.is_file() else None,
                        report=report if report.is_file() else None, maps=maps))
    return runs


def runs_for(model: ConfigModel) -> list[Run]:
    base, sid = dir_output(model), _subject_id(model)
    if base is None or sid is None or not base.is_dir():
        return []
    return discover_runs(base, f"sub-{sid:03d}")


def default_run(runs: list[Run], model: ConfigModel) -> int:
    """Index of the run the open config produces: its primary CSV as predict_outputs
    lists it (the joint result for async coupling), else its plain output affix."""
    if not runs:
        return -1
    names = [a.name for a in predict(model).artifacts if a.name.endswith(".csv")]
    stems = {r.stem + ".csv": i for i, r in enumerate(runs)}
    for name in names:
        if not name.endswith("_targeting.csv") and name in stems:
            return stems[name]
    affix = sanitize_affix(model.get("io.output_affix") or "")
    for i, r in enumerate(runs):
        if r.affix == affix:
            return i
    return 0


def related_runs(runs: list[Run], run: Run) -> list[Run]:
    """For an async joint result: the per-transducer runs it combines (same base affix)."""
    if run.kind != "joint":
        return []
    base = _JOINT.sub("", run.affix)
    return [r for r in runs if r.kind == "transducer" and _TRANS.sub("", r.affix) == base]


def read_metrics(path: Path | None) -> dict[str, str]:
    """Last row of a PRESTUS results CSV as {column: value}."""
    if path is None or not path.is_file():
        return {}
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    if len(rows) < 2:
        return {}
    return dict(zip(rows[0], rows[-1]))


def placement_positions(run: Run) -> list[tuple[list[float] | None, list[float] | None]]:
    """(trans_pos, focus_pos) per transducer, 1-based T1 voxels, as the run's placement
    stage recorded them. Async sub-runs and joint results fall back to the record of their
    base affix and of the async preprocessing stage; [] when nothing is recorded."""
    stem = _JOINT.sub("", _TRANS.sub("", run.affix))
    candidates = [run.affix, _TRANS.sub("", run.affix), stem, stem + "_desc-asyncPreproc"]
    for affix in dict.fromkeys(candidates):
        p = run.base / f"{run.tag}_placement{affix}.json"
        try:
            data = json.loads(p.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        txs = data.get("transducers") if isinstance(data, dict) else None
        if isinstance(txs, dict):
            txs = [txs]
        if not isinstance(txs, list):
            continue
        out = []
        for t in txs:
            t = t if isinstance(t, dict) else {}
            out.append(tuple(v if isinstance(v, list) and len(v) == 3 else None
                             for v in (t.get("trans_pos"), t.get("focus_pos"))))
        if run.kind == "transducer" and run.transducer and len(out) >= run.transducer:
            return [out[run.transducer - 1]]
        return out
    return []


def other_reports(base: Path) -> list[Path]:
    """Every HTML report in the folder (per-run, sequential and multi-ISPPA summaries)."""
    return sorted(base.glob("*report*.html"))


# ── Thermal time courses (cache/<tag>_<medium>_HEAT<affix>.mat) ──────────────────────
# label -> (struct, field): 'timeseries' holds the max over each tissue layer
# (thermal_update_timeseries.m); 'HEAT' the max over each tissue on the focal axis plus
# the TARGET voxel (thermal_plot_sim.m).
HEAT_METRICS: dict[str, tuple[str, str, str]] = {   # label -> (struct, field, unit)
    "Temperature (tissue max)": ("timeseries", "T", "°C"),
    "Temperature rise (tissue max)": ("timeseries", "Tdiff", "°C"),
    "CEM43 k-Wave (tissue max)": ("timeseries", "CEM43", "min"),
    "CEM43 ISO (tissue max)": ("timeseries", "CEM43_iso", "min"),
    "Temperature (focal axis)": ("HEAT", "temp", "°C"),
    "Temperature rise (focal axis)": ("HEAT", "tempRise", "°C"),
    "CEM43 (focal axis)": ("HEAT", "CEM43", "min"),
}


@dataclass
class HeatSeries:
    """One metric's time course per tissue."""
    time: "np.ndarray"
    tissues: dict[str, "np.ndarray"]
    unit: str


@dataclass
class HeatData:
    metrics: dict[str, HeatSeries]                  # HEAT_METRICS label -> series
    on: list[tuple[float, float]]                   # sonication-on intervals [s]


class HeatReadError(RuntimeError):
    """The HEAT file is missing, unreadable, or needs a package that is not installed."""


def read_heat(path: Path) -> HeatData:
    """Thermal time courses from a run's HEAT .mat file. The time axis follows
    thermal_plot_sim.m: TimeProfileTarget when it has one entry per sample (BabelBrain),
    else HEAT.recordedtime without its t=0 entry. Sonication-on intervals come from
    time_status_seq, which PRESTUS saves there since 2026-10-05 (else none)."""
    import numpy as np
    if not path.is_file():
        raise HeatReadError(f"No thermal time courses for this run ({path.name} not found).")
    try:
        raw = _load_mat(path)
    except HeatReadError:
        raise
    except Exception as exc:       # corrupt or unexpected layout
        raise HeatReadError(f"Cannot read {path.name}: {exc}") from exc
    ts, heat = raw.get("timeseries") or {}, raw.get("HEAT") or {}
    rec = np.ravel(np.asarray(heat.get("recordedtime", []), dtype=float))
    metrics: dict[str, HeatSeries] = {}
    for label, (struct, fld, unit) in HEAT_METRICS.items():
        src = (ts if struct == "timeseries" else heat).get(fld)
        if not isinstance(src, dict):
            continue
        tissues = {k: np.ravel(np.asarray(v, dtype=float)) for k, v in src.items()
                   if v is not None and not isinstance(v, (dict, list)) and np.size(v) > 0}
        if not tissues:
            continue
        n = min(len(v) for v in tissues.values())
        if struct == "timeseries":
            tp = np.ravel(np.asarray(ts.get("TimeProfileTarget", []), dtype=float))
            if tp.size == n:
                t = tp
            elif rec.size == n + 1:
                t = rec[1:]
            else:
                t = rec
        else:
            t = rec
        if t.size == 0:
            t = np.arange(n, dtype=float)
        n = min(n, t.size)
        metrics[label] = HeatSeries(t[:n], {k: v[:n] for k, v in tissues.items()}, unit)
    return HeatData(metrics, _on_intervals(raw.get("time_status_seq")))


def _on_intervals(seq) -> list[tuple[float, float]]:
    """(start, end) of every 'on' step of time_status_seq, as thermal_plot_sim.m shades
    them (step i spans the times of steps i-1 and i; first and last step skipped)."""
    if not isinstance(seq, list) or len(seq) < 3:
        return []
    out = []
    for i in range(1, len(seq) - 1):
        if str(seq[i].get("status", "")) == "on":
            out.append((float(seq[i - 1].get("time", 0.0)), float(seq[i].get("time", 0.0))))
    return out


def _load_mat(path: Path) -> dict:
    """A MATLAB -v7 (scipy) or -v7.3 (h5py) file as nested dicts/lists of arrays."""
    try:
        import scipy.io as sio
    except ImportError:
        sio = None
    if sio is not None:
        try:
            d = sio.loadmat(str(path), squeeze_me=True, struct_as_record=False)
            return {k: _from_scipy(v) for k, v in d.items() if not k.startswith("__")}
        except NotImplementedError:    # -v7.3
            pass
    try:
        import h5py
    except ImportError:
        raise HeatReadError(
            "Reading thermal time courses needs scipy (MATLAB -v7 files) or h5py (-v7.3): "
            "conda install scipy h5py") from None
    with h5py.File(str(path), "r") as h:
        return {k: _from_h5(h[k]) for k in ("timeseries", "HEAT") if k in h}


def _from_scipy(v):
    import numpy as np
    if hasattr(v, "_fieldnames"):
        return {f: _from_scipy(getattr(v, f)) for f in v._fieldnames}
    if isinstance(v, np.ndarray) and v.dtype == object:
        return [_from_scipy(x) for x in v.ravel()]
    return v


def _from_h5(node):
    """Numeric fields of a -v7.3 struct (struct arrays such as time_status_seq skipped)."""
    import h5py
    import numpy as np
    if isinstance(node, h5py.Group):
        return {k: _from_h5(node[k]) for k in node}
    data = node[()]
    return np.ravel(data) if np.issubdtype(np.asarray(data).dtype, np.number) else None
