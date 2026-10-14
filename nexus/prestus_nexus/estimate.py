"""Resource estimate for a study from earlier runs of the same project.

Each successful PRESTUS run appends one line to <path.sim>/run_stats.jsonl
(functions/helper/run_stats.m): what it simulated (grid voxels, time steps,
thermal steps), how long the acoustic and thermal simulations took, and the
job's peak RAM. This module picks the earlier runs that match the current
config (solver, medium, resolution, partition) and reports their acoustic
time, thermal time (scaled by the current protocol's step count) and peak
RAM as median and maximum.

Only the acoustic and thermal simulations are estimated. Segmentation,
pseudo-CT, BabelBrain preprocessing and queue wait are not.

Pure and fast: reads the config and one small file, never MATLAB.

Headless:
    python -m prestus_nexus.estimate study.yaml [--stats run_stats.jsonl]
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
from dataclasses import dataclass, field
from pathlib import Path

from .model import ConfigModel, is_empty

STATS_FILE = "run_stats.jsonl"
NOT_ESTIMATED = "segmentation, pseudo-CT, BabelBrain preprocessing, queue wait"


# ── inputs ───────────────────────────────────────────────────────────────
def stats_path(model: ConfigModel) -> Path | None:
    """<path.sim>/run_stats.jsonl, or None when path.sim is unset.

    path.sim may be a folder name inside path.anat (see config_default.yaml).
    """
    sim = model.get("path.sim")
    if is_empty(sim):
        return None
    sim = Path(str(sim)).expanduser()
    if not sim.is_absolute():
        anat = model.get("path.anat")
        if is_empty(anat):
            return None
        sim = Path(str(anat)).expanduser() / sim
    return sim / STATS_FILE


def load_records(path: str | Path | None) -> list[dict]:
    """Records from a run_stats.jsonl file; unreadable lines are skipped."""
    if path is None:
        return []
    try:
        text = Path(path).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return []
    out = []
    for line in text.splitlines():
        try:
            rec = json.loads(line)
        except ValueError:
            continue
        if isinstance(rec, dict):
            out.append(rec)
    return out


_CACHE: dict[Path, tuple[tuple[float, int], list[dict]]] = {}


def load_records_cached(path: Path | None) -> list[dict]:
    """load_records, re-read only when the file's mtime or size changed."""
    if path is None:
        return []
    try:
        st = path.stat()
    except OSError:
        return []
    stamp = (st.st_mtime, st.st_size)
    hit = _CACHE.get(path)
    if hit is None or hit[0] != stamp:
        hit = (stamp, load_records(path))
        _CACHE[path] = hit
    return hit[1]


def thermal_steps(model: ConfigModel) -> int | None:
    """Thermal time steps of the configured protocol.

    Mirrors functions/thermal/thermal_parameters.m (and the step count
    thermal_stage.m records); None when the protocol is incomplete.
    """
    try:
        t = {k: float(model.get(f"timing.{k}")) for k in
             ("pd", "pri", "ptd", "pt_timestep", "ptri", "ptrd", "post_ptri_dur", "post_pt_timestep")}
        equal = bool(model.get("timing.equal_step_duration"))
        if t["pri"] <= 0 or t["ptri"] <= 0 or t["post_pt_timestep"] <= 0:
            return None
        if equal:
            if t["pt_timestep"] <= 0:
                return None
            on_n = round(t["pd"] / t["pt_timestep"])
            off_n = round((t["pri"] - t["pd"]) / t["pt_timestep"])
        else:
            on_n = off_n = 1
        pulses = round(t["ptd"] / t["pri"])
        reps = round(t["ptrd"] / t["ptri"])
        ptri_off = t["ptri"] - t["ptd"]
        ptri_off_n = round(ptri_off / t["post_pt_timestep"]) if ptri_off > 0 else 0
        post_n = round(t["post_ptri_dur"] / t["post_pt_timestep"]) if t["post_ptri_dur"] > 0 else 0
    except (TypeError, ValueError):
        return None
    return reps * (pulses * (on_n + off_n) + ptri_off_n) + post_n


# ── matching ─────────────────────────────────────────────────────────────
def _num(v) -> float | None:
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
        return None
    return float(v)


def _str(v) -> str:
    return "" if is_empty(v) else str(v)


def _keys(model: ConfigModel) -> dict:
    return {
        "code_type": _str(model.get("simulation.code_type")),
        "acoustic_backend": _str(model.get("acoustic_backend")) or "kwave",
        "medium": _str(model.get("simulation.medium")),
        "resolution_mm": _num(model.get("grid.resolution_mm")),
        "hpc_partition": _str(model.get("hpc.partition")),
    }


def _matches(rec: dict, keys: dict, fields: tuple[str, ...]) -> bool:
    for f in fields:
        want, have = keys[f], rec.get(f)
        if f == "resolution_mm":
            have = _num(have)
            if want is None or have is None or abs(want - have) > 1e-9:
                return False
        elif _str(have) != want:
            return False
    return True


# Tried in order; the partition is dropped when no run on it matches.
_ACOUSTIC_FIELDS = (("code_type", "acoustic_backend", "medium", "resolution_mm", "hpc_partition"),
                    ("code_type", "acoustic_backend", "medium", "resolution_mm"))
_THERMAL_FIELDS = (("code_type", "medium", "resolution_mm", "hpc_partition"),
                   ("code_type", "medium", "resolution_mm"))


def _select(records: list[dict], keys: dict, tiers, value) -> tuple[list[float], tuple[str, ...]]:
    """Values from the first tier with any matching record that has one."""
    for fields in tiers:
        vals = [v for r in records if _matches(r, keys, fields) and (v := value(r)) is not None]
        if vals:
            return vals, fields
    return [], ()


# ── estimate ─────────────────────────────────────────────────────────────
@dataclass
class Figure:
    median: float
    maximum: float
    n: int
    basis: tuple[str, ...]          # config fields the runs matched on


@dataclass
class Estimate:
    source: Path | None
    n_records: int
    acoustic_s: Figure | None = None
    thermal_s: Figure | None = None
    peak_ram_gb: Figure | None = None
    notes: list[str] = field(default_factory=list)

    def lines(self) -> list[str]:
        if self.n_records == 0:
            where = self.source or "<path.sim>/" + STATS_FILE
            return [f"No earlier runs recorded ({where}).",
                    "Estimates appear once runs of this project have finished."]
        out = []
        if self.acoustic_s:
            out.append(_line("Acoustic sim", self.acoustic_s, _hms))
        if self.thermal_s:
            out.append(_line("Thermal sim", self.thermal_s, _hms))
        if self.peak_ram_gb:
            out.append(_line("Peak RAM (job)", self.peak_ram_gb, lambda g: f"{g:.1f} GB"))
        out += self.notes
        out.append(f"Not estimated: {NOT_ESTIMATED}.")
        return out


def _line(label: str, f: Figure, fmt) -> str:
    runs = f"{f.n} run" + ("s" if f.n != 1 else "")
    partition = "" if "hpc_partition" in f.basis else ", any partition"
    return f"{label:<15} ~{fmt(f.median)} (max {fmt(f.maximum)}; {runs}{partition})"


def _hms(s: float) -> str:
    m = int(round(s / 60))
    return f"{m // 60}h{m % 60:02d}" if m >= 60 else f"{max(m, 1)} min"


def _figure(vals: list[float], basis, scale: float = 1.0) -> Figure | None:
    if not vals:
        return None
    return Figure(statistics.median(vals) * scale, max(vals) * scale, len(vals), basis)


def estimate(model: ConfigModel, records: list[dict] | None = None,
             source: Path | None = None) -> Estimate:
    """Estimate from the given records, or from <path.sim>/run_stats.jsonl."""
    if records is None:
        source = source or stats_path(model)
        records = load_records_cached(source)
    est = Estimate(source=source, n_records=len(records))
    if not records:
        return est
    keys = _keys(model)

    if model.get("modules.run_acoustic_sims"):
        vals, basis = _select(records, keys, _ACOUSTIC_FIELDS, lambda r: _num(r.get("acoustic_s")))
        est.acoustic_s = _figure(vals, basis)
        if est.acoustic_s is None:
            est.notes.append("Acoustic sim: no earlier run with this solver, medium and resolution.")

    if model.get("modules.run_heating_sims"):
        steps = thermal_steps(model)

        def rate(r):
            s, n = _num(r.get("thermal_s")), _num(r.get("thermal_steps"))
            return s / n if s is not None and n else None
        vals, basis = _select(records, keys, _THERMAL_FIELDS, rate)
        if steps is None:
            est.notes.append("Thermal sim: timing protocol incomplete.")
        elif not vals:
            est.notes.append("Thermal sim: no earlier run with this solver, medium and resolution.")
        else:
            est.thermal_s = _figure(vals, basis, scale=steps)

    vals, basis = _select(records, keys, _ACOUSTIC_FIELDS, lambda r: _num(r.get("peak_ram_gb")))
    est.peak_ram_gb = _figure(vals, basis)
    return est


# ── headless ─────────────────────────────────────────────────────────────
def main(argv: list[str] | None = None) -> int:
    from .paths import base_default_for, find_prestus_root, toolbox_default
    from .schema import Schema

    ap = argparse.ArgumentParser(prog="python -m prestus_nexus.estimate",
                                 description="Estimate simulation time and peak RAM from earlier runs.")
    ap.add_argument("study", help="study YAML")
    ap.add_argument("--stats", help=f"run statistics file (default: <path.sim>/{STATS_FILE})")
    ap.add_argument("--prestus-root", help="PRESTUS checkout (default: this one)")
    args = ap.parse_args(argv)

    root = find_prestus_root(args.prestus_root)
    if root is None:
        print("PRESTUS root not found; pass --prestus-root.", file=sys.stderr)
        return 2
    schema = Schema(base_default_for(args.study, root), toolbox_path=toolbox_default(root))
    model = ConfigModel.from_file(schema, args.study)
    est = estimate(model, source=Path(args.stats) if args.stats else None)
    print("\n".join(est.lines()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
