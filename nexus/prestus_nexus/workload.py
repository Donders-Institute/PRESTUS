"""Predicted simulation workload, memory and time for the current config, before any run.

Workload follows the PRESTUS code:
  - acoustic grid: water_grid_dims.m (water/phantom, grid.default_dims 'auto'), the
    acoustic FOV box (resolve_acoustic_fov_params.m), or for a head the grid of this
    project's earlier runs (run_stats.jsonl acoustic_dims), else a reference head;
  - PML: getOptimalPMLSize (k-Wave) for pml_size 'auto';
  - time steps: source_sensor_setup.m (dt from PPW and CFL or source_ppp, t_end from
    steady_state_transits or sim_time_us, the jwave margin);
  - thermal steps: estimate.thermal_steps (thermal_parameters.m);
  - GPU memory: the acoustic_wrapper.m estimate (40 arrays, 70 for elastic).

Time, per stage of the run plan:
  - acoustic simulation and free-water baseline: t0 + grid points (PML included) x time
    steps x ns per point-step, by device and solver;
  - thermal simulation: thermal points x steps x ns;
  - stages that do not simulate (preprocessing, analysis, outputs; CPU work): seconds by
    medium and resolution.
Rates come from this project's runs when they recorded the size (or stage time), else from
throughput.yaml, which `fit` refits from run_stats.jsonl files or telemetry exports.

Headless:
    python -m prestus_nexus.workload study.yaml
    python -m prestus_nexus.workload fit [--update] <run_stats.jsonl | telemetry export> [...]
"""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import statistics
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

from . import yamlio
from .estimate import load_records_cached, stats_path, thermal_steps
from .model import ConfigModel, is_empty
from .placement import Bowl, resolve_focal_distance_ep
from .predict_outputs import is_async_mode, is_multi_isppa_mode

THROUGHPUT_FILE = Path(__file__).parent / "throughput.yaml"

# Reference head for a head model without earlier runs: the SimNIBS example head
# (sim_ernie m2m_sub-001 final_tissues, all labels incl. neck), axis-aligned extent
# 168 x 225 x 237 mm; its bounding box over 300 random orientations holds 9.1-14.4 L
# (median 11.8 L). Measured 2026-11-16; the transducer stand-off is not included.
REF_HEAD_MM = (168.0, 225.0, 237.0)
REF_HEAD_L = (9.07, 11.84, 14.42)          # low, median, high


def _num(v) -> float | None:
    if isinstance(v, list):
        v = v[0] if len(v) == 1 else None
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
        return None
    return float(v)


# ── k-Wave helpers ───────────────────────────────────────────────────────
def _max_prime(n: int) -> int:
    p, f, out = n, 2, 1
    while f * f <= p:
        while p % f == 0:
            out, p = f, p // f
        f += 1
    return max(out, p) if p > 1 else out


def optimal_pml(n: int, lo: int = 10, hi: int = 40) -> int:
    """getOptimalPMLSize for one dimension: the PML in [lo, hi] whose n + 2*pml has the
    smallest largest prime factor (first one on ties)."""
    return min(range(lo, hi + 1), key=lambda p: (_max_prime(n + 2 * p), p))


def pml_sizes(model: ConfigModel, dims: tuple[int, ...]) -> tuple[int, ...]:
    v = model.get("grid.pml_size")
    fixed = _num(v)
    if fixed is not None:
        return tuple(int(fixed) for _ in dims)
    rng = model.get("grid.pml_size_range")
    lo, hi = (int(rng[0]), int(rng[1])) if isinstance(rng, list) and len(rng) == 2 else (10, 40)
    return tuple(optimal_pml(n, lo, hi) for n in dims)


# ── grids ────────────────────────────────────────────────────────────────
@dataclass
class Grid:
    dims: tuple[int, ...]          # simulated grid, PML not included
    basis: str                     # where the size comes from
    approx: bool = False           # True: a reference size, not this subject's
    points_range: tuple[float, float] | None = None   # low/high grid points when approx

    @property
    def points(self) -> int:
        return math.prod(self.dims)


def _tx(model: ConfigModel) -> dict:
    txs = model.effective_transducers()
    return txs[0] if txs else {}


def _aperture_and_bowl_depth(model: ConfigModel) -> tuple[float | None, float | None]:
    """Outer element diameter and bowl apex -> focus distance [mm] (water_grid_dims.m)."""
    tx = _tx(model)
    kind = tx.get("type")
    sub = (tx.get(kind) or {}) if kind in ("annular", "matrix") else {}
    ap = None
    od = sub.get("elem_od_mm") if kind == "annular" else None
    od = od if od is not None else tx.get("elem_od_mm")
    if isinstance(od, list) and od and all(_num(x) is not None for x in od):
        ap = max(float(x) for x in od)
    elif _num(od) is not None:
        ap = _num(od)
    elif kind == "matrix":
        ap = _num(sub.get("outer_diameter_mm"))
    fb = _num(tx.get("focal_distance_bowl"))
    if fb is None:
        bowl = Bowl.from_transducer(tx)
        ep = resolve_focal_distance_ep(tx, bowl)
        if ep is not None and bowl is not None:
            fb = ep + bowl.offset
    return ap, fb


def water_dims(model: ConfigModel, dx: float) -> tuple[int, ...] | None:
    """water_grid_dims.m: lateral aperture + 2 x margin, axial apex->focus + beyond."""
    ap, fb = _aperture_and_bowl_depth(model)
    if ap is None or fb is None:
        return None
    margin = _num(model.get("grid.water_margin_mm")) or 5.0
    beyond = _num(model.get("grid.z_beyond_focus_mm"))
    beyond = 40.0 if beyond is None else beyond
    n_lat = 2 * math.ceil((ap + 2 * margin) / dx / 2) + 1
    n_ax = math.ceil((fb + beyond) / dx) + 2
    return (n_lat, n_lat, n_ax)


def _fov_dims(model: ConfigModel, dx: float) -> tuple[int, ...] | None:
    d, ln = model.get("grid.acoustic_fov_diameter_mm"), model.get("grid.acoustic_fov_length_mm")
    if is_empty(d) or is_empty(ln):
        return None
    ap, fb = _aperture_and_bowl_depth(model)
    beyond = _num(model.get("grid.z_beyond_focus_mm"))
    beyond = 40.0 if beyond is None else beyond
    d = ap if d == "auto" else _num(d)
    ln = (fb + beyond if fb is not None else None) if ln == "auto" else _num(ln)
    if d is None or ln is None:
        return None
    n = math.ceil(d / dx)
    return (n, n, math.ceil(ln / dx))


def _project_head_dims(records: list[dict], dx: float, medium: str) -> tuple[tuple[int, ...], int] | None:
    """Median 3D acoustic dims of this project's earlier full-grid runs at this resolution."""
    dims = [tuple(int(x) for x in r["acoustic_dims"]) for r in records
            if r.get("medium") == medium and _num(r.get("resolution_mm")) == dx
            and isinstance(r.get("acoustic_dims"), list) and len(r["acoustic_dims"]) == 3
            and not r.get("axisymmetric")]
    if not dims:
        return None
    by_points = sorted(dims, key=math.prod)
    return by_points[len(by_points) // 2], len(dims)


def acoustic_grid(model: ConfigModel, records: list[dict]) -> Grid | None:
    dx = _num(model.get("grid.resolution_mm"))
    if dx is None or dx <= 0:
        return None
    medium = model.get("simulation.medium") or "layered"
    fov = _fov_dims(model, dx)
    if fov is not None:
        return Grid(fov, "acoustic FOV")
    if medium in ("water", "phantom"):
        dd = model.get("grid.default_dims")
        if isinstance(dd, list) and len(dd) in (2, 3) and all(_num(x) for x in dd):
            dims = tuple(int(x) for x in dd)
            basis = "grid.default_dims"
        else:
            dims = water_dims(model, dx)
            basis = "grid sized to the transducer"
            if dims is None:
                return None
        if medium == "phantom" and basis != "grid.default_dims":
            return None        # phantom grids follow the phantom geometry, not the transducer
        return Grid(dims, basis)
    hit = _project_head_dims(records, dx, medium)
    if hit is not None:
        dims, n = hit
        return Grid(dims, f"median of {n} earlier head run" + ("s" if n != 1 else ""))
    lo, mid, hi = (v * 1e6 / dx ** 3 for v in REF_HEAD_L)
    scale = (REF_HEAD_L[1] * 1e6 / math.prod(REF_HEAD_MM)) ** (1 / 3)
    dims = tuple(math.ceil(e * scale / dx) for e in REF_HEAD_MM)
    return Grid(dims, "reference head (SimNIBS example)", approx=True, points_range=(lo, hi))


def axisymmetric(model: ConfigModel) -> bool:
    return _num(model.get("grid.axisymmetric")) == 1


def simulated_dims(model: ConfigModel, grid: Grid) -> tuple[int, ...]:
    """The plane k-Wave runs for an axisymmetric grid (axial x lateral // 2 + 1, as the
    ITRUSST logs show: 140 lateral -> 71), else the grid."""
    if axisymmetric(model) and len(grid.dims) == 3:
        return (grid.dims[2], grid.dims[0] // 2 + 1)
    if len(grid.dims) == 2 and axisymmetric(model):
        return (grid.dims[1], grid.dims[0] // 2 + 1)
    return grid.dims


# ── time steps ───────────────────────────────────────────────────────────
def sound_speeds(model: ConfigModel) -> tuple[float, float] | None:
    medium = model.get("simulation.medium") or "layered"
    props = model.get("medium_properties") or {}
    layers = model.get("layers")
    used = set(layers) | {"water"} if isinstance(layers, dict) and layers else set(props)
    tissues = ["water"] if medium == "water" else [t for t in props if t in used]
    c = [_num((props.get(t) or {}).get("sound_speed")) for t in tissues if isinstance(props.get(t), dict)]
    c = [x for x in c if x]
    return (min(c), max(c)) if c else None


def acoustic_steps(model: ConfigModel, dims: tuple[int, ...]) -> int | None:
    """Nt as source_sensor_setup.m sets it."""
    dx = _num(model.get("grid.resolution_mm"))
    freqs = [_num(t.get("freq_hz")) for t in model.effective_transducers()]
    freqs = [f for f in freqs if f]
    speeds = sound_speeds(model)
    if not dx or not freqs or speeds is None:
        return None
    c_min, c_max = speeds
    f = max(freqs)
    period = 1 / f
    backend = str(model.get("acoustic_backend") or "kwave")
    ppp = _num(model.get("grid.source_ppp"))
    if ppp is None:
        ppw = _num(model.get("grid.source_ppw")) or c_max / (f * dx / 1e3)
        cfl = _num(model.get("grid.source_cfl"))
        if cfl is None:
            cfl = 0.10 if backend.startswith("jwave") else 0.15
        ppp = math.ceil(ppw / cfl)
    dt = period / ppp
    sim_us = _num(model.get("grid.sim_time_us"))
    if sim_us is not None:
        t_end = sim_us * 1e-6
    else:
        transits = _num(model.get("grid.steady_state_transits")) or 1.0
        t_end = transits * math.sqrt(sum((n * dx / 1e3) ** 2 for n in dims)) / c_min
        if backend == "jwave":
            pml = _num(model.get("grid.pml_size")) or 40
            t_end += 2 * pml * dx / 1e3 / c_min
    return round(t_end / dt)


# ── prediction ───────────────────────────────────────────────────────────
@dataclass
class Rate:
    ns: float                       # per grid point per step (acoustic: PML included)
    lo: float
    hi: float
    basis: str
    t0: float = 0.0                 # fixed seconds per simulation (solver start, file I/O)


@dataclass
class Part:
    """A simulation: seconds = runs x (t0 + points x steps x ns)."""
    label: str
    runs: int                       # simulations of this kind in the run
    points: float                   # grid points per simulation (with PML for acoustic)
    steps: int
    rate: Rate | None = None

    @property
    def point_steps(self) -> float:
        return self.points * self.steps * self.runs

    def seconds(self) -> tuple[float, float, float] | None:
        if self.rate is None:
            return None
        w = self.points * self.steps * 1e-9
        return tuple(self.runs * (self.rate.t0 + w * ns) for ns in (self.rate.ns, self.rate.lo, self.rate.hi))


@dataclass
class Fixed:
    """A stage that does not simulate (CPU work, outputs): seconds from earlier runs."""
    label: str
    s: float
    lo: float
    hi: float
    basis: str

    def seconds(self) -> tuple[float, float, float]:
        return self.s, self.lo, self.hi


@dataclass
class Prediction:
    grid: Grid | None = None
    parts: list[Part] = field(default_factory=list)
    stages: list[Fixed] = field(default_factory=list)
    untimed: list[str] = field(default_factory=list)      # stages that run but have no timing data
    gpu_mem_gb: float | None = None
    notes: list[str] = field(default_factory=list)

    def total(self) -> tuple[float, float, float] | None:
        """Sum over simulations and stages; None while a simulation is untimed."""
        secs = [p.seconds() for p in self.parts] + [s.seconds() for s in self.stages]
        if not secs or any(s is None for s in secs):
            return None
        return tuple(sum(s[i] for s in secs) for i in range(3))

    def lines(self) -> list[str]:
        out = []
        if self.grid is not None:
            g = self.grid
            size = " x ".join(str(n) for n in g.dims)
            if g.approx and g.points_range:
                out.append(f"Grid: {size} ({_m(g.points)} points; {_m(g.points_range[0])}-"
                           f"{_m(g.points_range[1])} for heads of that size), from the {g.basis}")
            else:
                out.append(f"Grid: {size} ({_m(g.points)} points), {g.basis}")
        for p in self.parts:
            runs = f"{p.runs} x " if p.runs > 1 else ""
            work = f"{runs}{_m(p.points)} points x {p.steps:,} steps"
            s = p.seconds()
            if s is None:
                out.append(f"{p.label}: {work}; time unknown: no timing data yet for this solver and device")
            else:
                out.append(f"{p.label}: ~{_span(s)}; {work}; {p.rate.basis}")
        if self.stages:
            st = [s.seconds() for s in self.stages]
            out.append(f"Other stages: ~{_span(tuple(sum(x[i] for x in st) for i in range(3)))} ("
                       + ", ".join(f"{s.label.lower()} {_hms(s.s)}" for s in self.stages)
                       + "); CPU, the same on any GPU")
        if self.untimed:
            out.append("Not timed: " + ", ".join(self.untimed))
        t = self.total()
        if t is not None and (len(self.parts) + len(self.stages)) > 1:
            out.append(f"Total: ~{_span(t)}" + (" plus the stages not timed" if self.untimed else ""))
        if self.gpu_mem_gb is not None:
            mem = self.gpu_mem_gb
            out.append(f"GPU memory (acoustic): ~{mem:.0f} GB" if mem >= 10 else
                       f"GPU memory (acoustic): ~{mem:.1f} GB" if mem >= 0.1 else "GPU memory (acoustic): <0.1 GB")
        return out + self.notes

    def summary(self) -> str:
        """One short line for a page footer: the total and GPU memory, or which simulation has
        no timing yet. The breakdown is in lines() (the footer's tooltip and Details)."""
        if self.grid is None:
            return "Set the resolution and the transducer to predict the load."
        if not self.parts:
            return "No acoustic or thermal simulation"
        if any(n.startswith(NO_LOCAL_GPU) for n in self.notes):
            return f"{NO_LOCAL_GPU} on this computer: time unknown"
        bits = []
        t = self.total()
        if t is not None:
            bits.append(f"~{_span(t)}" + (" plus untimed stages" if self.untimed else ""))
        else:
            unknown = [p.label.replace(" sim", "").lower() for p in self.parts if p.seconds() is None]
            more = f" and {len(unknown) - 1} more" if len(unknown) > 1 else ""
            bits.append(f"{unknown[0]}{more}: time unknown")
        if self.gpu_mem_gb is not None:
            mem = self.gpu_mem_gb
            bits.append(f"GPU ~{mem:.0f} GB" if mem >= 10 else f"GPU ~{mem:.1f} GB" if mem >= 0.1
                        else "GPU <0.1 GB")
        return " · ".join(bits)

def _m(n: float) -> str:
    return f"{n / 1e6:.1f} M" if n >= 1e6 else f"{n / 1e3:.0f} k"


def _hms(s: float) -> str:
    m = int(round(s / 60))
    if m >= 60:
        return f"{m // 60}h{m % 60:02d}"
    return f"{m} min" if s >= 90 else f"{max(round(s), 1)} s"


def _span(s: tuple[float, float, float]) -> str:
    """'7 min (5-11 min)'; no range when both ends round alike."""
    lo, hi = _hms(s[1]), _hms(s[2])
    return _hms(s[0]) if lo == hi else f"{_hms(s[0])} ({lo}-{hi})"


def load_throughput(path: Path | None = None) -> dict:
    try:
        return yamlio.load_file_cached(path or THROUGHPUT_FILE) or {}
    except OSError:
        return {}


def _device_kind(model: ConfigModel) -> str:
    return "gpu" if "gpu" in str(model.get("simulation.code_type") or "") else "cpu"


def _norm(text: str) -> str:
    """GPU names compared across sources: 'NVIDIA A100-PCIE-40GB' (gpuDevice) and
    'gpu:nvidia_a100-pcie-40gb:1' (SLURM gres) both contain 'a100-pcie-40gb'."""
    t = re.sub(r"[\s_]+", "-", str(text).lower())
    return t.replace("nvidia-", "").replace("tesla-", "")


def _aliases(device: str) -> list[str]:
    """Match keys of a GPU name: the name, and for a MIG slice also the SLURM gres
    spelling ('A100 80GB PCIe MIG 3g.40gb' -> 'a100-80gb-pcie-3g', as in
    'nvidia_a100_80gb_pcie_3g.39gb')."""
    n = _norm(device)
    out = [n]
    mig = re.sub(r"-mig-(\d+g)\.\d+gb$", r"-\1", n)
    if mig != n:
        out.append(mig)
    return out


def _best_match(items: list[dict], hint: str) -> dict | None:
    """The item whose longest match key occurs in hint ('a100-80gb-pcie-3g' beats 'a100-80gb-pcie')."""
    best, length = None, 0
    for e in items:
        for m in e.get("match") or []:
            k = _norm(m)
            if k and k in hint and len(k) > length:
                best, length = e, len(k)
    return best


def _solver(model: ConfigModel) -> dict:
    return {"code_type": str(model.get("simulation.code_type") or ""),
            "backend": str(model.get("acoustic_backend") or "kwave"),
            "axisymmetric": int(axisymmetric(model)),
            "precision": str(model.get("simulation.precision") or "single")}


def _entry_fits(e: dict, want: dict) -> bool:
    """Optional solver keys of an entry (acoustic entries from `fit`) must match."""
    return all(k not in e or e[k] == v for k, v in want.items())


def _entry_rate(e: dict, factor: float = 1.0, note: str = "") -> Rate:
    lo, hi = (e.get("range") or [e["ns"], e["ns"]])[:2]
    return Rate(float(e["ns"]) * factor, float(lo) * factor, float(hi) * factor,
                f"{e['device']} (n={e.get('n', '?')}){note}", float(e.get("t0_s") or 0.0))


def _fallbacks(pool: list[dict], want: dict, kind: str, factors: dict) -> list[tuple[list[dict], float, str]]:
    """Candidate entries in order of preference: this solver (factor 1), then the reference
    solver (cpp_gpu / cpp_cpu, single precision) times the code-type and precision factors."""
    out = [([e for e in pool if _entry_fits(e, want)], 1.0, "")]
    over, notes, factor = {}, [], 1.0
    ref = "cpp_gpu" if kind == "gpu" else "cpp_cpu"
    fx = (factors.get("code_type") or {}).get(want["code_type"])
    if want["code_type"] != ref and isinstance(fx, dict) and _num(fx.get("x")):
        over["code_type"], factor = ref, factor * float(fx["x"])
        notes.append(f"x{float(fx['x']):.1f} for {want['code_type']}")
    fp = (factors.get("precision") or {}).get(want["precision"])
    if want["precision"] != "single" and isinstance(fp, dict) and _num(fp.get("x")):
        over["precision"], factor = "single", factor * float(fp["x"])
        notes.append(f"x{float(fp['x']):.1f} for {want['precision']} precision")
    if over:
        out.append(([e for e in pool if _entry_fits(e, {**want, **over})], factor, "; " + ", ".join(notes)))
    return [c for c in out if c[0]]


_LOCAL_GPU: list = []          # cache: [name or None] once nvidia-smi has been asked


def local_gpu() -> str | None:
    """Name of this computer's first NVIDIA GPU (nvidia-smi), asked once; None without one."""
    if _LOCAL_GPU:
        return _LOCAL_GPU[0]
    name = None
    exe = shutil.which("nvidia-smi")
    if exe:
        try:
            out = subprocess.run([exe, "--query-gpu=name", "--format=csv,noheader"], capture_output=True,
                                 text=True, timeout=5, creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            lines = [ln.strip() for ln in out.stdout.splitlines() if ln.strip()]
            name = lines[0] if out.returncode == 0 and lines else None
        except (OSError, subprocess.SubprocessError):
            name = None
    _LOCAL_GPU.append(name)
    return name


def device_hint(model: ConfigModel) -> tuple[str, str | None]:
    """(match text, detected GPU): hpc.gpu and hpc.partition; for a GPU run on this computer
    (platform: matlab) without hpc.gpu, the local GPU's name from nvidia-smi."""
    hint = _norm(f"{model.get('hpc.gpu') or ''} {model.get('hpc.partition') or ''}")
    if (not is_empty(model.get("hpc.gpu")) or _device_kind(model) != "gpu"
            or str(model.get("platform") or "") != "matlab"):
        return hint, None
    name = local_gpu()
    return (_norm(name), name) if name else (hint, None)


NO_LOCAL_GPU = "No NVIDIA GPU"


def no_local_gpu(model: ConfigModel) -> bool:
    """A GPU run on this computer (platform: matlab, no hpc.gpu) that has no NVIDIA GPU: the
    cluster's GPU timings do not apply, and MATLAB's GPU code cannot run here."""
    return (_device_kind(model) == "gpu" and str(model.get("platform") or "") == "matlab"
            and is_empty(model.get("hpc.gpu")) and local_gpu() is None)


def table_rate(entries: list, model: ConfigModel, table: dict | None = None) -> Rate | None:
    """table_rate_for on the device of device_hint; a detected local GPU is named in the basis.
    None for a GPU run on this computer without a GPU (no_local_gpu)."""
    if no_local_gpu(model):
        return None
    hint, detected = device_hint(model)
    rate = table_rate_for(entries, model, table, hint)
    if rate is not None and detected:
        rate.basis += f"; this computer: {detected}"
    return rate


def table_rate_for(entries: list, model: ConfigModel, table: dict | None, hint: str) -> Rate | None:
    """Rate for the config's solver on the device named in hint, from throughput.yaml entries,
    first that holds:
    1. the entry for this GPU and solver;
    2. this GPU's entry for the reference solver (cpp_gpu / cpp_cpu, single) x the factors;
    3. (GPU) entries of the same generation scaled by memory bandwidth, this solver first;
    4. the span of all entries of the kind, this solver first."""
    table = table or {}
    kind = _device_kind(model)
    pool = [e for e in entries or [] if isinstance(e, dict) and e.get("kind") == kind and _num(e.get("ns"))]
    options = _fallbacks(pool, _solver(model), kind, table.get("factors") or {})
    if not options:
        return None
    for cands, factor, note in options:
        hit = _best_match(cands, hint)
        if hit is not None:
            return _entry_rate(hit, factor, note)
    if kind == "gpu":
        for cands, factor, note in options:
            bw = _bandwidth_rate(cands, hint, table.get("bandwidth") or [])
            if bw is not None:
                return Rate(bw.ns * factor, bw.lo * factor, bw.hi * factor, bw.basis + note, bw.t0)
    cands, factor, note = options[0]
    if len(cands) == 1:
        return _entry_rate(cands[0], factor, note)
    ns = [float(e["ns"]) for e in cands]
    t0 = statistics.median(float(e.get("t0_s") or 0.0) for e in cands)
    return Rate(statistics.median(ns) * factor, min(ns) * factor, max(ns) * factor,
                f"any {kind.upper()} ({', '.join(e['device'] for e in cands)}); set hpc.gpu to narrow" + note, t0)


def _bandwidth_rate(cands: list[dict], hint: str, bandwidth: list[dict]) -> Rate | None:
    """A GPU without a measured entry: k-Wave is memory-bound, so ns x GB/s is about
    constant within a GPU generation. Scale the measured entries of the same generation."""
    target = _best_match(bandwidth, hint)
    if target is None or not _num(target.get("gbs")):
        return None
    k, devs, t0 = [], [], []
    for e in cands:
        spec = _best_match(bandwidth, " ".join(_norm(m) for m in e.get("match") or [e["device"]]))
        if spec is not None and spec.get("arch") == target.get("arch") and _num(spec.get("gbs")):
            k.append(float(e["ns"]) * float(spec["gbs"]))
            devs.append(e["device"])
            t0.append(float(e.get("t0_s") or 0.0))
    if not k:
        return None
    ns = statistics.median(k) / float(target["gbs"])
    lo, hi = min(k) / float(target["gbs"]) * 0.8, max(k) / float(target["gbs"]) * 1.25
    name = target.get("name") or (target.get("match") or ["this GPU"])[0]
    return Rate(ns, lo, hi, f"{name}: scaled by memory bandwidth from {', '.join(devs)}", statistics.median(t0))


# ── records: run_stats.jsonl lines and telemetry events ─────────────────
def flatten(r: dict) -> dict | None:
    """A run_stats.jsonl line as is; a telemetry run_end event (its payload and sim_stats
    merged to the top level, duration_s as total_s). Other events and failed runs: None."""
    if not isinstance(r, dict):
        return None
    if isinstance(r.get("payload"), dict):
        r = {**r, **r["payload"]}
    if r.get("event") is not None and r.get("event") != "run_end":
        return None
    if r.get("status") not in (None, "", "success", "ok"):
        return None
    if isinstance(r.get("sim_stats"), dict):
        r = {**r, **r["sim_stats"]}
    if _num(r.get("total_s")) is None and _num(r.get("duration_s")) is not None:
        r = {**r, "total_s": r["duration_s"]}
    return r


def _device(r: dict) -> tuple[str, str]:
    code = str(r.get("code_type") or "")
    kind = "gpu" if "gpu" in code else "cpu"
    if kind == "gpu":
        return kind, str(r.get("gpu_model") or "")
    n = _num(r.get("n_cpu")) or _num(r.get("n_cpu_alloc"))
    return kind, f"CPU, {int(n)} cores" if n else "CPU"


def _record_solver(r: dict) -> dict:
    return {"code_type": str(r.get("code_type") or ""),
            "backend": str(r.get("acoustic_backend") or "kwave"),
            "axisymmetric": int(bool(_num(r.get("axisymmetric")))),
            "precision": str(r.get("precision") or "single")}


def acoustic_work(r: dict) -> tuple[float, bool] | None:
    """(grid points with PML x time steps, approximate) of a record's acoustic simulation.
    Records from before 2026-11-16 hold no dims or PML: a cube of the recorded voxels with
    k-Wave's optimal PML stands in (approximate). Nt 'auto' (axisymmetric runs then): None."""
    nt = _num(r.get("acoustic_nt"))
    if not nt or not _num(r.get("acoustic_s")):
        return None
    dims, pml = r.get("acoustic_dims"), _num(r.get("acoustic_pml"))
    if isinstance(dims, list) and dims and pml is not None:
        return math.prod(int(d) + 2 * int(pml) for d in dims) * nt, False
    vox = _num(r.get("acoustic_voxels"))
    if not vox or _num(r.get("axisymmetric")):
        return None
    n = round(vox ** (1 / 3))
    return (n + 2 * optimal_pml(n)) ** 3 * nt, True


def _medium_group(medium) -> str:
    m = str(medium or "")
    return m if m in ("water", "phantom") else "head"


def stage_seconds(r: dict) -> dict[str, float]:
    """stage_<label>_s fields of a record; a segmentation under a minute (pseudo-CT time
    taken out) reused charm output and is left out."""
    out = {}
    for k, v in r.items():
        m = re.fullmatch(r"stage_(\w+)_s", k)
        if m and _num(v) is not None and _num(v) >= 0:
            out[m.group(1)] = _num(v)
    # the segmentation stage holds pseudo-CT generation (preproc_segmentation.m): charm alone
    if "segmentation" in out:
        out["segmentation"] -= out.get("pct", 0.0)
        if out["segmentation"] < 60:
            del out["segmentation"]
    # head preprocessing has two paths (preproc_head.m): a shared head model placed per
    # run (headmodel.reuse_across_placements; times preproc_placement), else rotation,
    # scaling, and skull smoothing per run (preproc_rotscale, preproc_skullsmooth). Without
    # either sub-step the run loaded its preprocessed head from the cache. The three differ
    # in cost, so they are fitted apart.
    if "preproc" in out:
        if "preproc_placement" in out:
            out["preproc_shared"] = out.pop("preproc")
        elif not {"preproc_rotscale", "preproc_skullsmooth"} & set(out):
            out["preproc_cached"] = out.pop("preproc")
    return out


def _preproc_label(model: ConfigModel) -> str:
    return "preproc_shared" if model.get("headmodel.reuse_across_placements") in (True, 1) else "preproc"


def _record_rates(records: list[dict], model: ConfigModel, stage: str, t0: float = 0.0) -> list[float]:
    """ns per point-step of this project's runs with the same solver (and partition when any)."""
    want = _solver(model)
    part = str(model.get("hpc.partition") or "")

    def rate(r: dict) -> float | None:
        if r.get("code_type") != want["code_type"]:
            return None
        if stage == "acoustic":
            got = _record_solver(r)
            if any(got[k] != want[k] for k in ("backend", "axisymmetric")):
                return None
            work = acoustic_work(r)
            if work is None or work[1]:
                return None
            return max(_num(r["acoustic_s"]) - t0, 0.1 * _num(r["acoustic_s"])) / work[0] * 1e9
        n, s, st = _num(r.get("thermal_voxels")), _num(r.get("thermal_s")), _num(r.get("thermal_steps"))
        return s / (n * st) * 1e9 if n and s and st else None

    same = [v for r in records if str(r.get("hpc_partition") or "") == part and (v := rate(r))]
    return same or [v for r in records if (v := rate(r))]


def _rate(records, model, stage, table) -> Rate | None:
    tab = table_rate(table.get(stage) or [], model, table)
    t0 = tab.t0 if tab is not None else (DEFAULT_T0_S if stage == "acoustic" else 0.0)
    vals = _record_rates(records, model, stage, t0)
    if vals:
        q = statistics.quantiles(vals, n=10) if len(vals) >= 10 else (min(vals), max(vals))
        return Rate(statistics.median(vals), q[0], q[-1], f"this project's runs (n={len(vals)})", t0)
    return tab


DEFAULT_T0_S = 5.0      # acoustic start-up per simulation when no entry gives one (k-Wave, 2026-11-17 fit)

# Run plan stage -> log_timer labels (stage_<label>_s). Stages that simulate are Parts.
STAGE_TIMERS = {
    "Segmentation": ("segmentation",),
    "Placement": ("placement",),
    "Grid and head model": ("preproc",),
    "Medium": ("medium", "nifti_medium"),
    "Source": ("source",),
    "Time reversal": ("time_reversal",),
    "Acoustic analysis": ("acoustic_analysis",),
    "NIfTI maps": ("nifti_acoustic",),
    "Thermal analysis": ("thermal_analysis",),
    "Skull CT": ("pct",),
    "Report": ("report",),
    "Group report": ("group_report",),
}
SIM_TIMERS = ("acoustic", "thermal", "freefield_baseline")
UNTIMED = {"Targeting analysis": None}     # None: short (inside the NIfTI stage), not listed
UNTIMED_NAMES = {"Segmentation": "segmentation (SimNIBS charm)", "Skull CT": "pseudo-CT",
                 "Report": "HTML report", "Group report": "group report"}


def _stage_estimate(label: str, model: ConfigModel, records: list[dict], table: dict
                    ) -> tuple[float, float, float, str] | None:
    """Seconds for one log_timer label: this project's runs at this resolution and medium
    (3 or more), else the table (s_1mm x (1 mm / dx)^(3 x exponent))."""
    dx = _num(model.get("grid.resolution_mm"))
    group = _medium_group(model.get("simulation.medium"))
    own = [v for r in records if _medium_group(r.get("medium")) == group
           and _num(r.get("resolution_mm")) == dx and (v := stage_seconds(r).get(label)) is not None]
    if len(own) >= 3:
        q = statistics.quantiles(own, n=10) if len(own) >= 10 else (min(own), max(own))
        return statistics.median(own), q[0], q[-1], f"this project's runs (n={len(own)})"
    for e in table.get("stages") or []:
        if isinstance(e, dict) and e.get("stage") == label and e.get("medium") == group and _num(e.get("s_1mm")):
            s = float(e["s_1mm"]) * ((1.0 / dx) ** (3 * float(e.get("exponent") or 0)) if dx else 1.0)
            lo, hi = (e.get("ratio") or [1, 1])[:2]
            return s, s * float(lo), s * float(hi), f"n={e.get('n', '?')}"
    return None


def baseline_dims(model: ConfigModel, dims: tuple[int, ...]) -> tuple[int, ...]:
    """The free-water baseline's grid (water_baseline_box.m): a box around the source and
    the beam to 40 mm past the focus (10 mm margin, sides rounded up to 8), else the
    full grid (axisymmetric, positions unknown, or a box over 80% of the grid)."""
    dx = _num(model.get("grid.resolution_mm"))
    ap, fb = _aperture_and_bowl_depth(model)
    if axisymmetric(model) or len(dims) != 3 or not dx or ap is None or fb is None:
        return dims
    tx = _tx(model)
    curv = _num((tx.get(tx.get("type")) or {}).get("curv_radius_mm")) if tx.get("type") else None
    reach = max(fb, curv or 0.0)
    lat = math.ceil((ap / dx + 2 * math.ceil(10 / dx)) / 8) * 8
    ax = math.ceil(((reach + 40) / dx + 2 * math.ceil(10 / dx)) / 8) * 8
    box = (lat, lat, ax)
    if math.prod(box) > 0.8 * math.prod(dims):
        return dims
    return box


def run_counts(model: ConfigModel) -> tuple[int, int]:
    """(acoustic simulations, thermal simulations) of one pipeline run, by mode."""
    n_tx = len(model.tx)
    targets = max((len([x for x in (v if isinstance(v, list) else [v]) if _num(x) is not None])
                   for v in (model.get("transducer.target_isppa_wcm2", ti) for ti in range(n_tx))), default=1)
    n_ac = n_tx if is_async_mode(model) else 1
    n_th = max(targets, 1) if (is_multi_isppa_mode(model) or is_async_mode(model)) else 1
    if model.get("simulation.uncertainty") and str(model.get("simulation.medium") or "") != "water":
        return 5 * n_ac, 3 * n_th        # default + Temp-H/L heat; Trans-H/L acoustic only
    return n_ac, n_th


def predict(model: ConfigModel, records: list[dict] | None = None,
            throughput: dict | None = None) -> Prediction:
    from . import run_plan as rp
    if records is None:
        records = load_records_cached(stats_path(model))
    table = throughput if throughput is not None else load_throughput()
    pred = Prediction()
    stages = {s.name: s.state for s in rp.plan(model).stages}
    grid = acoustic_grid(model, records)
    pred.grid = grid
    if grid is None:
        pred.notes.append("Grid size unknown: set the resolution, transducer aperture, and focal distance"
                          " (phantoms: grid.default_dims).")
        return pred
    n_ac, n_th = run_counts(model)
    elastic = str(model.get("acoustic_backend") or "") == "kwave_elastic"
    precision = str(model.get("simulation.precision") or "single")
    heating = stages.get("Thermal simulation") == rp.RUN

    if stages.get("Acoustic simulation") == rp.RUN:
        if str(model.get("acoustic_backend") or "") == "babelbrain":
            pred.notes.append("BabelBrain runs on its own grid: acoustic time not predicted.")
        else:
            sim = simulated_dims(model, grid)
            pml = pml_sizes(model, sim)
            comp = math.prod(n + 2 * p for n, p in zip(sim, pml))
            nt = acoustic_steps(model, sim)
            if nt is not None:
                rate = _rate(records, model, "acoustic", table)
                pred.parts.append(Part("Acoustic sim", n_ac, comp, nt, rate))
                water = {"Free-water baseline": 1, "Post-hoc free-water run": 1}
                if str(model.get("simulation.medium")) != "water":
                    for name, runs in water.items():
                        if stages.get(name) == rp.RUN:
                            box = baseline_dims(model, sim)
                            pts = math.prod(n + 2 * p for n, p in zip(box, pml_sizes(model, box)))
                            pred.parts.append(Part(name, runs, pts, nt, rate))
            if _device_kind(model) == "gpu":
                pred.gpu_mem_gb = (70 if elastic else 40) * comp * (8 if precision == "double" else 4) / 1e9
    if heating:
        steps = thermal_steps(model)
        if steps is None:
            pred.notes.append("Thermal sim: timing protocol incomplete.")
        else:
            th = grid.points
            if axisymmetric(model) and len(grid.dims) == 2:
                th = grid.dims[0] ** 2 * grid.dims[1]      # expanded to 3D for heating
            res_th = _num(model.get("grid.thermal_resolution_mm"))
            dx = _num(model.get("grid.resolution_mm"))
            if res_th and dx:
                th *= (dx / res_th) ** 3
            pred.parts.append(Part("Thermal sim", n_th, th, steps, _rate(records, model, "thermal", table)))

    repeats = {"acoustic_analysis": n_ac, "nifti_acoustic": n_ac, "thermal_analysis": n_th, "nifti_thermal": n_th}
    for name, state in stages.items():
        if state != rp.RUN:
            continue
        labels = STAGE_TIMERS.get(name)
        if labels is None:
            if UNTIMED.get(name):
                pred.untimed.append(UNTIMED[name])
            continue
        labels = labels + (("nifti_thermal",) if name == "NIfTI maps" and heating else ())
        labels = tuple(_preproc_label(model) if lb == "preproc" else lb for lb in labels)
        est = [(lb, _stage_estimate(lb, model, records, table)) for lb in labels]
        if any(e is None for _, e in est):
            pred.untimed.append(UNTIMED_NAMES.get(name, name.lower()))
            continue
        k = [repeats.get(lb, 1) for lb, _ in est]
        s = tuple(sum(n * e[i] for n, (_, e) in zip(k, est)) for i in range(3))
        pred.stages.append(Fixed(name, *s, est[0][1][3]))
    if grid.approx:
        pred.notes.append("Head size assumed: the subject's own grid replaces it after its first run.")
    if no_local_gpu(model):
        pred.notes.append(f"{NO_LOCAL_GPU} found on this computer: a GPU code type cannot run here. "
                          "Submit to a cluster, or choose a CPU code type (simulation.code_type).")
    return pred


# ── refit throughput.yaml ────────────────────────────────────────────────
def _quantiles(vals: list[float]) -> tuple[float, float]:
    if len(vals) >= 10:
        q = statistics.quantiles(vals, n=10)
        return q[0], q[-1]
    return min(vals), max(vals)


def _fit_two_term(w: list[float], t: list[float], t0: float | None) -> tuple[float, float]:
    """t = t0 + c1 x w, least squares on relative error. t0 given: fit c1 only. A free fit
    needs 6 records over a 3-fold range of work; a negative t0 falls back to the given one."""
    if t0 is None and len(w) >= 6 and max(w) >= 3 * min(w):
        a = [1 / x for x in t]
        s11 = sum(x * x for x in a)
        s12 = sum(x * x * y for x, y in zip(a, w))
        s22 = sum((x * y) ** 2 for x, y in zip(a, w))
        b1, b2 = sum(x for x in a), sum(x * y for x, y in zip(a, w))     # sum a^2 t = sum a
        det = s11 * s22 - s12 ** 2
        if det > 0:
            f0 = (b1 * s22 - b2 * s12) / det
            c1 = (s11 * b2 - s12 * b1) / det
            if 0 <= f0 <= 120 and c1 > 0:
                return f0, c1
    t0 = DEFAULT_T0_S if t0 is None else t0
    return t0, statistics.median(max(y - t0, 0.1 * y) / x for x, y in zip(w, t))


def fit(records: list[dict]) -> dict:
    """throughput.yaml sections from run records (run_stats.jsonl lines or telemetry run_end
    events): acoustic t0 + ns per point-step by device and solver, thermal ns per point-step,
    non-simulating stages by medium and resolution, and code-type and precision factors.
    Records with a 'study' key also yield the number of studies per entry."""
    recs = [f for r in records if (f := flatten(r)) is not None]
    out: dict[str, list | dict] = {"acoustic": [], "thermal": [], "stages": [], "factors": {}}

    # acoustic: per device and solver; free t0 where the work varies enough, else pooled
    groups: dict[tuple, list[tuple[float, float, str]]] = {}
    for r in recs:
        work = acoustic_work(r)
        if work is None:
            continue
        kind, dev = _device(r)
        sv = _record_solver(r)
        key = (kind, dev, sv["code_type"], sv["backend"], sv["axisymmetric"], sv["precision"])
        groups.setdefault(key, []).append((work[0], _num(r["acoustic_s"]), str(r.get("study") or "")))
    fits: dict[tuple, tuple[float, float, list]] = {}
    for key, rows in groups.items():
        for _ in range(2):                                   # once more without outliers
            t0, c1 = _fit_two_term([x[0] for x in rows], [x[1] for x in rows], None)
            ratios = [s / (t0 + c1 * w) for w, s, _ in rows]
            keep = [row for row, q in zip(rows, ratios) if 0.4 <= q <= 2.5]
            if len(keep) == len(rows) or len(keep) < 2:
                break
            rows = keep
        fits[key] = (t0, c1, rows)
    for key, (t0, c1, rows) in sorted(fits.items(), key=lambda kv: str(kv[0])):
        kind, dev, code, backend, axi, prec = key
        ratios = [s / (t0 + c1 * w) for w, s, _ in rows]
        lo, hi = _quantiles(ratios)
        e = {"device": dev, "kind": kind, "match": _aliases(dev) if kind == "gpu" else [],
             "code_type": code, "backend": backend, "axisymmetric": axi, "precision": prec,
             "ns": round(c1 * 1e9, 3), "t0_s": round(t0, 1),
             "range": [round(c1 * 1e9 * lo, 3), round(c1 * 1e9 * hi, 3)], "n": len(rows)}
        studies = {st for _, _, st in rows if st}
        if studies:
            e["studies"] = len(studies)
        out["acoustic"].append(e)

    # factors: same device, backend and geometry; code type vs the reference solver, precision vs single
    def ratio_pairs(vary: str, ref_of) -> dict:
        acc: dict[str, list[float]] = {}
        for k1, (_, c1, _) in fits.items():
            kind, dev, code, backend, axi, prec = k1
            val = {"code": code, "prec": prec}[vary]
            ref = ref_of(kind)
            if val == ref:
                continue
            k0 = (kind, dev, ref if vary == "code" else code, backend, axi, "single" if vary == "prec" else prec)
            if vary == "code" and prec != "single":
                continue
            if k0 in fits:
                acc.setdefault(val, []).append(c1 / fits[k0][1])
        return {v: {"x": round(statistics.median(xs), 2), "n": len(xs)} for v, xs in sorted(acc.items())}
    out["factors"] = {"code_type": ratio_pairs("code", lambda k: "cpp_gpu" if k == "gpu" else "cpp_cpu"),
                      "precision": ratio_pairs("prec", lambda k: "single")}

    # thermal: ns per point-step by device
    tgroups: dict[tuple, list[float]] = {}
    for r in recs:
        n, s, st = _num(r.get("thermal_voxels")), _num(r.get("thermal_s")), _num(r.get("thermal_steps"))
        if n and s and st:
            tgroups.setdefault(_device(r), []).append(s / (n * st) * 1e9)
    for (kind, dev), vals in sorted(tgroups.items()):
        lo, hi = _quantiles(vals)
        out["thermal"].append({"device": dev, "kind": kind, "match": _aliases(dev) if kind == "gpu" else [],
                               "ns": round(statistics.median(vals), 2), "range": [round(lo, 2), round(hi, 2)],
                               "n": len(vals)})

    # stages: seconds = s_1mm x (1 / dx^3)^exponent, by medium group (CPU work, any device)
    sgroups: dict[tuple, list[tuple[float, float]]] = {}
    for r in recs:
        dx = _num(r.get("resolution_mm"))
        if not dx:
            continue
        heated = (_num(r.get("thermal_s")) or 0) > 0 or (_num(r.get("stage_thermal_s")) or 0) >= 1
        for label, s in stage_seconds(r).items():
            # simulations are modelled from their work; preproc_* are parts of preproc;
            # thermal outputs only count from runs that simulated heating
            if (s > 0 and label not in SIM_TIMERS
                    and (label in ("preproc_shared", "preproc_cached") or not label.startswith("preproc_"))
                    and (heated or label not in ("thermal_analysis", "nifti_thermal"))):
                sgroups.setdefault((label, _medium_group(r.get("medium"))), []).append((1.0 / dx ** 3, s))
    for (label, group), rows in sorted(sgroups.items()):
        lx = [math.log(x) for x, _ in rows]
        ly = [math.log(s) for _, s in rows]
        b = 0.0
        xs = [x for x, _ in rows]
        if len(rows) >= 6 and max(xs) >= 1.3 ** 3 * min(xs):
            mx, my = statistics.fmean(lx), statistics.fmean(ly)
            sxx = sum((v - mx) ** 2 for v in lx)
            b = min(max(sum((u - mx) * (v - my) for u, v in zip(lx, ly)) / sxx, 0.0), 1.5) if sxx else 0.0
        resid = [v - b * u for u, v in zip(lx, ly)]
        a = statistics.median(resid)
        lo, hi = _quantiles([math.exp(v - a) for v in resid])
        out["stages"].append({"stage": label, "medium": group, "s_1mm": round(math.exp(a), 1),
                              "exponent": round(b, 2), "ratio": [round(lo, 2), round(hi, 2)], "n": len(rows)})
    return out


_ENTRY_KEYS = {"acoustic": ("device", "kind", "code_type", "backend", "axisymmetric", "precision"),
               "thermal": ("device", "kind"), "stages": ("stage", "medium")}


def merge(old: dict, new: dict) -> dict:
    """The refitted table: a fitted entry replaces the old entry with the same key (device
    and solver; stage and medium), other old entries stay; fitted factors replace the old
    ones by name; the hand-entered bandwidth table is kept."""
    out = dict(old or {})
    for sec, keys in _ENTRY_KEYS.items():
        fresh = [e for e in new.get(sec) or [] if isinstance(e, dict)]
        if not fresh:
            continue
        ids = {tuple(e.get(k) for k in keys) for e in fresh}
        kept = [e for e in out.get(sec) or [] if isinstance(e, dict) and tuple(e.get(k) for k in keys) not in ids]
        out[sec] = kept + fresh
    for kind, vals in (new.get("factors") or {}).items():
        if vals:
            out.setdefault("factors", {})
            out["factors"] = {**out["factors"], kind: {**(out["factors"].get(kind) or {}), **vals}}
    return out


def validate(records: list[dict], by: str = "study") -> dict:
    """Leave-one-group-out check of the acoustic model: fit without a group's records,
    predict them, and return the observed/predicted ratios (median, p10, p90, n)."""
    recs = [f for r in records if (f := flatten(r)) is not None]
    groups = sorted({str(r.get(by) or "") for r in recs} - {""})
    ratios = []
    for g in groups:
        train = [r for r in recs if str(r.get(by) or "") != g]
        table = fit(train)
        for r in recs:
            if str(r.get(by) or "") != g or (work := acoustic_work(r)) is None:
                continue
            kind, dev = _device(r)
            sv = _record_solver(r)
            for e in table["acoustic"]:
                if e["device"] == dev and e["kind"] == kind and all(e[k] == sv[k] for k in sv):
                    ratios.append(_num(r["acoustic_s"]) / (e["t0_s"] + e["ns"] * 1e-9 * work[0]))
                    break
    if not ratios:
        return {"n": 0}
    lo, hi = _quantiles(ratios)
    return {"median": round(statistics.median(ratios), 3), "p10": round(lo, 3), "p90": round(hi, 3),
            "n": len(ratios)}


def read_records(path) -> list[dict]:
    """Records from a run_stats.jsonl file or a telemetry export (JSON lines, or one JSON
    array of events). Unreadable lines are skipped."""
    from .estimate import load_records
    try:
        text = Path(path).read_text(encoding="utf-8", errors="replace").lstrip()
    except OSError:
        return []
    if text.startswith("["):
        try:
            data = json.loads(text)
        except ValueError:
            return []
        return [r for r in data if isinstance(r, dict)]
    return load_records(path)


THROUGHPUT_HEADER = """Simulation and stage times for the Nexus load prediction (workload.py).

acoustic: seconds = t0_s + ns x 1e-9 x grid points (PML included) x time steps, per device
  and solver (code_type, backend, axisymmetric, precision). range: ns at the p10/p90 of
  observed/predicted (min/max when n < 10).
thermal: ns per thermal grid point per step.
stages: stages that do not simulate (CPU work, outputs; the same on any GPU), by medium
  (head, water, phantom): seconds = s_1mm x (1 mm / resolution)^(3 x exponent); ratio: p10/p90.
factors: a solver without its own entry = the reference solver (cpp_gpu / cpp_cpu, single) x.
bandwidth: GPU name, memory bandwidth [GB/s, NVIDIA spec], and generation, to scale a
  measured entry to a GPU without one (k-Wave is memory-bound). Entered by hand; kept by refits.
match: lower-case substrings of hpc.gpu / hpc.partition; the longest match wins.

Refit (records: run_stats.jsonl lines or telemetry run_end events, JSONL or JSON):
  python -m prestus_nexus.workload fit [--update] <records> [...]"""


# ── headless ─────────────────────────────────────────────────────────────
def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "fit":
        fp = argparse.ArgumentParser(prog="python -m prestus_nexus.workload fit",
                                     description="Refit throughput.yaml from run records.")
        fp.add_argument("records", nargs="+", help="run_stats.jsonl files or telemetry exports (JSONL or JSON)")
        fp.add_argument("--update", action="store_true", help=f"write {THROUGHPUT_FILE.name} (else print)")
        args = fp.parse_args(argv[1:])
        recs = [r for p in args.records for r in read_records(p)]
        table = merge(load_throughput(), fit(recs))
        if args.update:
            yamlio.dump_file(THROUGHPUT_FILE, table, THROUGHPUT_HEADER)
            print(f"{THROUGHPUT_FILE}: {len(table.get('acoustic') or [])} acoustic, "
                  f"{len(table.get('stages') or [])} stage entries from {len(recs)} records")
        else:
            print(yamlio.dump(table, THROUGHPUT_HEADER))
        return 0
    from .paths import base_default_for, find_prestus_root, toolbox_default
    from .schema import Schema

    ap = argparse.ArgumentParser(prog="python -m prestus_nexus.workload",
                                 description="Predict simulation workload, memory, and time.")
    ap.add_argument("study", help="study YAML")
    ap.add_argument("--prestus-root", help="PRESTUS checkout (default: this one)")
    args = ap.parse_args(argv)
    root = find_prestus_root(args.prestus_root)
    if root is None:
        print("PRESTUS root not found; pass --prestus-root.", file=sys.stderr)
        return 2
    schema = Schema(base_default_for(args.study, root), toolbox_path=toolbox_default(root))
    model = ConfigModel.from_file(schema, args.study)
    print("\n".join(predict(model).lines()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
