"""−3 dB focal widths (FWHM of the focal intensity) for the Results tab's overview. Qt-free.

Runs since PRESTUS writes them carry ``focal_fwhm_axial_mm``, ``focal_fwhm_lateral_min_mm``
and ``focal_fwhm_lateral_max_mm`` in the metrics CSV (compute_focal_fwhm.m, on the
simulation grid). For earlier runs the Nexus measures the same widths on the run's T1w
intensity map, with the same method:

- the focal peak is the intensity maximum inside the run's −6 dB focal mask
  (brain-restricted for layered runs, beyond the exit plane for free water), else the
  maximum of the whole map;
- the axial width is measured along the transducer axis (bowl apex to planned focus, from
  the run's placement record), the lateral widths every 5° in the plane perpendicular to
  it through the peak; the narrowest and widest are reported;
- profiles are sampled every quarter of the smallest voxel by trilinear interpolation, and
  each crossing of half the peak intensity is placed by linear interpolation between the
  two samples that bracket it. A side without a crossing inside the map gives no width.

The T1w map is resampled from the simulation grid, so these widths can differ from the
grid's by a fraction of a voxel. Results are cached per map (path, size, modification
time and axis) in memory and in the user's cache folder.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

SOURCE_MAP = "T1w intensity map (Nexus)"
STEP_FRACTION = 0.25                 # sample spacing, as a fraction of the smallest voxel
LATERAL_ANGLES_DEG = tuple(range(0, 180, 5))
_VERSION = 1                         # bump when the method changes (invalidates the disk cache)


@dataclass
class FocalWidths:
    axial_mm: float | None
    lateral_min_mm: float | None
    lateral_max_mm: float | None
    source: str = SOURCE_MAP
    note: str = ""


# ── method ─────────────────────────────────────────────────────────────────

def crossing_width(profile: np.ndarray, centre: int, step: float, threshold: float) -> float | None:
    """Distance between the first samples below ``threshold`` on either side of ``centre``,
    each crossing placed by linear interpolation; None when a side runs out of samples or
    into NaN (outside the map) first."""
    edges = []
    for sgn in (1, -1):
        prev = profile[centre]
        i = centre + sgn
        edge = None
        while 0 <= i < len(profile):
            cur = profile[i]
            if np.isnan(cur):
                break
            if cur < threshold:
                frac = (prev - threshold) / (prev - cur)
                edge = (abs(i - centre) - 1 + frac) * step
                break
            prev = cur
            i += sgn
        if edge is None:
            return None
        edges.append(edge)
    return float(sum(edges))


def _sample(data: np.ndarray, inv_affine: np.ndarray, points_mm: np.ndarray) -> np.ndarray:
    from scipy.ndimage import map_coordinates
    vox = (inv_affine @ np.c_[points_mm, np.ones(len(points_mm))].T)[:3]
    return map_coordinates(data, vox, order=1, mode="constant", cval=np.nan, prefilter=False)


def widths_at(data: np.ndarray, affine: np.ndarray, peak_vox, axis_mm) -> tuple[float | None, ...]:
    """(axial, lateral min, lateral max) FWHM [mm] through the voxel ``peak_vox`` (0-based)
    of ``data``; ``axis_mm`` is the transducer axis in world (affine) coordinates."""
    data = np.asarray(data, dtype=float)
    affine = np.asarray(affine, dtype=float)
    inv = np.linalg.inv(affine)
    u = np.asarray(axis_mm, dtype=float)
    if not np.linalg.norm(u):
        return None, None, None
    u = u / np.linalg.norm(u)
    vsize = np.sqrt((affine[:3, :3] ** 2).sum(axis=0))
    step = float(vsize.min()) * STEP_FRACTION
    reach = float(np.linalg.norm(np.array(data.shape) * vsize))
    n = int(np.ceil(reach / step))
    s = np.arange(-n, n + 1) * step
    p0 = (affine @ np.r_[np.asarray(peak_vox, dtype=float), 1.0])[:3]
    peak = float(data[tuple(int(round(c)) for c in peak_vox)])
    if not peak > 0:
        return None, None, None
    thr = peak / 2

    def width(d) -> float | None:
        prof = _sample(data, inv, p0 + s[:, None] * d[None, :])
        return crossing_width(prof, n, step, thr)

    axial = width(u)
    e = np.zeros(3)
    e[int(np.argmin(np.abs(u)))] = 1.0
    v1 = np.cross(u, e)
    v1 /= np.linalg.norm(v1)
    v2 = np.cross(u, v1)
    lat = [width(np.cos(a) * v1 + np.sin(a) * v2) for a in np.deg2rad(LATERAL_ANGLES_DEG)]
    if any(w is None for w in lat):
        return axial, None, None
    return axial, min(lat), max(lat)


def focal_peak(data: np.ndarray, mask: np.ndarray | None) -> tuple[int, int, int]:
    """Voxel of the intensity maximum, inside ``mask`` when it has any voxel."""
    if mask is not None and mask.shape == data.shape and np.any(mask > 0.5):
        masked = np.where(mask > 0.5, data, -np.inf)
        return tuple(int(i) for i in np.unravel_index(np.argmax(masked), data.shape))
    return tuple(int(i) for i in np.unravel_index(np.nanargmax(data), data.shape))


# ── per run ────────────────────────────────────────────────────────────────

_memory: dict[str, FocalWidths] = {}


def _cache_key(map_path: Path, mask_path: Path | None, axis: tuple) -> str | None:
    parts = [str(_VERSION)]
    for p in (map_path, mask_path):
        if p is None:
            parts.append("-")
            continue
        try:
            st = p.stat()
        except OSError:
            return None
        parts.append(f"{p.absolute()}|{st.st_size}|{st.st_mtime_ns}")
    parts.append(repr(axis))
    return hashlib.sha1("\n".join(parts).encode()).hexdigest()[:20]


def _cache_file(key: str) -> Path:
    from .localcopy import cache_root
    return cache_root() / "focal_widths" / f"{key}.json"


def cached(map_path: Path, mask_path: Path | None, axis: tuple) -> FocalWidths | None:
    """The widths of an unchanged map if they were measured before (no map read)."""
    key = _cache_key(map_path, mask_path, axis)
    if key is None:
        return None
    if key in _memory:
        return _memory[key]
    try:
        d = json.loads(_cache_file(key).read_text(encoding="utf-8"))
        fw = FocalWidths(**d)
    except (OSError, ValueError, TypeError):
        return None
    _memory[key] = fw
    return fw


def measure(map_path: Path, mask_path: Path | None, trans_vox1, focus_vox1) -> FocalWidths:
    """Widths on a T1w intensity map; trans/focus are the run's 1-based T1 voxels (placement
    record). Reads the map (slow on a network share): call off the GUI thread."""
    axis = tuple(float(a) for a in (*trans_vox1, *focus_vox1))
    hit = cached(map_path, mask_path, axis)
    if hit is not None:
        return hit
    from .placement import Volume
    vol = Volume(map_path)
    mask = None
    note = ""
    if mask_path is not None:
        try:
            mask = Volume(mask_path).data
        except Exception:  # best-effort: fall back to the map maximum (nibabel raises many types)
            from .ui.errors import log_current
            log_current(f"focal widths: cannot read the focal mask {mask_path}")
            mask = None
    if mask is None:
        note = "peak is the map maximum (no focal mask map)"
    t = vol.to_ras(np.asarray(trans_vox1, dtype=float) - 1)[0]
    f = vol.to_ras(np.asarray(focus_vox1, dtype=float) - 1)[0]
    peak = focal_peak(np.asarray(vol.data, dtype=float), mask)
    ax, lo, hi = widths_at(vol.data, vol.affine, peak, f - t)
    fw = FocalWidths(ax, lo, hi, SOURCE_MAP, note)
    key = _cache_key(map_path, mask_path, axis)
    if key is not None:
        _memory[key] = fw
        try:
            from .fsio import write_text_atomic
            path = _cache_file(key)
            path.parent.mkdir(parents=True, exist_ok=True)
            write_text_atomic(path, json.dumps(asdict(fw)))
        except OSError:  # best-effort: disk cache only; the result is kept in memory
            pass
    return fw


def inputs_for(run) -> tuple[Path, Path | None, list, list] | None:
    """(intensity map, −6 dB mask map or None, trans_pos, focus_pos) for a run, or None when
    it has no T1w intensity map or no recorded transducer axis."""
    from . import results as rs
    path = run.maps.get("Intensity (ISPPA)")
    if path is None:
        return None
    pos = rs.placement_positions(run)
    if not pos or pos[0][0] is None or pos[0][1] is None:
        return None
    return path, run.maps.get("Focal volume (-6 dB)"), pos[0][0], pos[0][1]


def from_csv(metrics: dict) -> FocalWidths | None:
    """The widths PRESTUS wrote into the run's CSV, if it has the columns."""
    keys = ("focal_fwhm_axial_mm", "focal_fwhm_lateral_min_mm", "focal_fwhm_lateral_max_mm")
    if not any(k in metrics for k in keys):
        return None
    vals = []
    for k in keys:
        try:
            x = float(metrics.get(k, "nan"))
        except (TypeError, ValueError):
            x = float("nan")
        vals.append(None if np.isnan(x) else x)
    return FocalWidths(*vals, source="run CSV")
