"""Grid orientation of a multi-transducer run: which direction becomes the grid z axis
under grid.mode, and how far each beam is tilted from it.

Ports of the MATLAB logic (code is ground truth):

- functions/head/mean_beam_axis.m: the beams are averaged as axes (lines): a beam pointing
  against transducer(1)'s is reversed first; below ``min_norm`` the mean axis is ambiguous
  and transducer(1)'s axis is used instead.
- functions/head/preproc_head.m (grid.mode switch, ~lines 137-165): 'focal_axis' (default)
  -> preproc_align_to_mean_focal_axis, 'transducer_axis' -> preproc_align_to_focal_axis,
  'ras_plus' -> preproc_rescale_only (no rotation), 'principal_axis' -> error
  PRESTUS:gridMode:removedPrincipalAxis, anything else -> warning and 'focal_axis'.
- functions/head/preproc_align_to_mean_focal_axis.m (~lines 66-89): min_norm is
  parameters.grid.focal_axis_min_norm when set, else 0.5; transducer(1)'s beam honours
  align_to_focus (false: the beam is the +z voxel axis), transducers 2..N use
  focus_pos - trans_pos.

MATLAB works in the voxel grid of the T1 after ensure_ras_plus (preproc_head.m ~79-93),
i.e. voxel axes along R, A, S. Positions passed here must be in such a space: RAS+
canonical voxel indices (placement.Volume.ras_to_canon) are exact; RAS mm agree for
isotropic voxels. Only directions matter, so the voxel origin does not.
"""

from __future__ import annotations

import math
from typing import Any, Callable, Iterable

import numpy as np

MODES = ("focal_axis", "transducer_axis", "ras_plus")
REMOVED_MODES = {"principal_axis": "grid.mode 'principal_axis' has been removed; use 'focal_axis' "
                                   "(the default), which aligns z with the mean beam axis of all "
                                   "transducers (preproc_head.m)"}
DEFAULT_MIN_NORM = 0.5          # preproc_align_to_mean_focal_axis.m / mean_beam_axis.m default
RAS_Z = np.array([0.0, 0.0, 1.0])

# Guidance for the "tilted beam on a coarse grid" note (not a PRESTUS rule). Evidence:
# doc_multitransducer_async.md ('Grid mode choice and accuracy'): at 2 mm (2.5 points per
# wavelength at 300 kHz) a beam 41 deg off axis lost 16 % of its free-water ISPPA. The
# thresholds are set well inside that case: tilts from 20 deg, grids below 6 points per
# wavelength in water (the grid.min_ppw default of source_sensor_setup.m).
TILT_WARN_DEG = 20.0
PPW_WARN = 6.0
WATER_SOUND_SPEED = 1500.0      # m/s, as the PPW check of checks.grid_issues for water


class GridModeError(ValueError):
    """grid.mode PRESTUS refuses (principal_axis), or positions missing."""


def _vec(v: Any) -> np.ndarray | None:
    if v is None:
        return None
    try:
        a = np.asarray(v, dtype=float).reshape(-1)
    except (TypeError, ValueError):
        return None
    return a if a.size == 3 and np.all(np.isfinite(a)) else None


def _truthy(v: Any, default: bool = True) -> bool:
    if v is None or (isinstance(v, (list, str)) and len(v) == 0):
        return default
    return bool(v) and v not in ("0", "false", "False")


def beam_unit(trans: Any, focus: Any) -> np.ndarray | None:
    """Unit vector transducer -> focus, or None when either is missing or they coincide."""
    t, f = _vec(trans), _vec(focus)
    if t is None or f is None:
        return None
    d = f - t
    n = float(np.linalg.norm(d))
    return d / n if n > 0 else None


def beam_units(transducers: list[dict]) -> list[np.ndarray | None]:
    """Beam unit vector per transducer, as preproc_align_to_mean_focal_axis.m builds them:
    transducer 1 with ``align_to_focus`` false points along +z (its 'natural focus' lies at
    trans_pos + [0, 0, curv_radius / resolution]); the others use focus_pos - trans_pos."""
    out: list[np.ndarray | None] = []
    for i, tx in enumerate(transducers):
        if i == 0 and not _truthy(tx.get("align_to_focus"), True):
            out.append(RAS_Z.copy() if _vec(tx.get("trans_pos")) is not None else None)
        else:
            out.append(beam_unit(tx.get("trans_pos"), tx.get("focus_pos")))
    return out


def mean_beam_axis(units: Iterable[Any], min_norm: float = DEFAULT_MIN_NORM) -> tuple[np.ndarray, float, bool]:
    """Port of mean_beam_axis.m: (axis unit vector, norm of the mean before normalising,
    fell_back). Row 0 is the reference beam; beams against it are flipped first."""
    u = np.asarray([np.asarray(x, float) for x in units], dtype=float)
    if u.ndim != 2 or u.shape[1] != 3 or len(u) == 0:
        raise GridModeError("mean_beam_axis needs at least one 3-D beam vector")
    u = u / np.linalg.norm(u, axis=1, keepdims=True)
    flip = (u @ u[0]) < 0
    u[flip] = -u[flip]
    mean = u.mean(axis=0)
    norm = float(np.linalg.norm(mean))
    if norm < min_norm:
        return u[0].copy(), norm, True
    return mean / norm, norm, False


def grid_z_axis(transducers: list[dict], mode: str | None = "focal_axis",
                min_norm: float | None = DEFAULT_MIN_NORM) -> tuple[np.ndarray, bool]:
    """The direction PRESTUS rotates onto the grid z axis, and whether 'focal_axis' fell
    back to transducer 1's axis (mean axis norm < min_norm; only possible with 3+ beams).

    ``transducers``: dicts with trans_pos / focus_pos (and optionally align_to_focus) in a
    common RAS+-oriented space. Raises GridModeError for grid.mode 'principal_axis' and when
    a needed position is missing."""
    mode = mode or "focal_axis"
    if mode in REMOVED_MODES:
        raise GridModeError(REMOVED_MODES[mode])
    if mode == "ras_plus":
        return RAS_Z.copy(), False
    if mode not in MODES:
        mode = "focal_axis"               # preproc_head.m: unrecognised -> warn, focal_axis
    if not transducers:
        raise GridModeError("no transducer configured")
    units = beam_units(transducers)
    if units[0] is None:
        raise GridModeError("transducer 1 needs trans_pos and focus_pos")
    if mode == "transducer_axis":
        return units[0], False
    missing = [i + 1 for i, u in enumerate(units) if u is None]
    if missing:
        raise GridModeError("transducer " + ", ".join(map(str, missing)) + " needs trans_pos and focus_pos")
    mn = DEFAULT_MIN_NORM if min_norm is None else float(min_norm)
    axis, _norm, fell_back = mean_beam_axis(units, mn)
    return axis, fell_back


def beam_tilts(transducers: list[dict], z: Any) -> list[float | None]:
    """Angle [deg] between each beam's axis (a line) and the grid z axis, 0-90: how
    obliquely the beam crosses the grid. None for a transducer without positions."""
    zv = np.asarray(z, float)
    zv = zv / np.linalg.norm(zv)
    out: list[float | None] = []
    for u in beam_units(transducers):
        if u is None:
            out.append(None)
        else:
            c = min(1.0, abs(float(np.dot(u, zv))))
            out.append(math.degrees(math.acos(c)))
    return out


def min_norm_of(model) -> float:
    """grid.focal_axis_min_norm as preproc_align_to_mean_focal_axis.m reads it."""
    v = model.get("grid.focal_axis_min_norm")
    return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else DEFAULT_MIN_NORM


def points_per_wavelength(model, sound_speed: float = WATER_SOUND_SPEED) -> float | None:
    """Grid points per wavelength at the highest transducer frequency (None if unknown)."""
    dx = model.get("grid.resolution_mm")
    freqs = []
    for t in model.effective_transducers():
        f = t.get("freq_hz")
        f = f[0] if isinstance(f, (list, tuple)) and f else f
        if isinstance(f, (int, float)) and not isinstance(f, bool) and f > 0:
            freqs.append(float(f))
    if not freqs or not isinstance(dx, (int, float)) or isinstance(dx, bool) or dx <= 0:
        return None
    return sound_speed / (max(freqs) * dx * 1e-3)


def model_transducers(model, positions: Callable[[int], tuple] | None = None,
                      to_space: Callable[[np.ndarray], np.ndarray] | None = None) -> list[dict]:
    """The model's transducers with trans_pos / focus_pos in one RAS+-oriented space.

    ``positions(ti)`` -> (trans, focus) in RAS mm (e.g. PlacementPage._positions_for, which
    also resolves Localite / heuristic placements); ``to_space`` maps RAS mm to RAS+ voxel
    indices (Volume.ras_to_canon) so anisotropic T1 voxels are handled as MATLAB does.
    Without ``positions`` the configured values are used as they are: voxel indices are
    taken to be stored RAS+ (true for SimNIBS T1 images), mm as RAS."""
    out = []
    for ti, t in enumerate(model.effective_transducers()):
        d = {"align_to_focus": t.get("align_to_focus")}
        if positions is not None:
            tr, fo = positions(ti)
            if to_space is not None:
                tr = None if tr is None else np.asarray(to_space(np.asarray(tr, float)), float).reshape(-1)
                fo = None if fo is None else np.asarray(to_space(np.asarray(fo, float)), float).reshape(-1)
            d["trans_pos"], d["focus_pos"] = tr, fo
        else:
            d["trans_pos"], d["focus_pos"] = t.get("trans_pos"), t.get("focus_pos")
        out.append(d)
    return out


def orientation(model, transducers: list[dict] | None = None) -> dict:
    """Everything the Transducers page and the viewer show about the grid orientation:
    {'mode', 'z', 'fell_back', 'tilts', 'error', 'min_norm', 'ppw', 'coarse_tilted'}."""
    mode = str(model.get("grid.mode") or "focal_axis")
    txs = transducers if transducers is not None else model_transducers(model)
    mn = min_norm_of(model)
    res = {"mode": mode, "z": None, "fell_back": False, "tilts": [None] * len(txs), "error": "",
           "min_norm": mn, "ppw": points_per_wavelength(model), "coarse_tilted": []}
    try:
        z, fell = grid_z_axis(txs, mode, mn)
    except GridModeError as exc:
        res["error"] = str(exc)
        return res
    res["z"], res["fell_back"] = z, fell
    res["tilts"] = beam_tilts(txs, z)
    ppw = res["ppw"]
    if ppw is not None and ppw < PPW_WARN:
        res["coarse_tilted"] = [i for i, a in enumerate(res["tilts"]) if a is not None and a >= TILT_WARN_DEG]
    return res
