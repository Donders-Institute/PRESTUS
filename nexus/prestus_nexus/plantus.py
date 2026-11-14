"""PlanTUS page logic without Qt: the maps of a PlanTUS run, the skin vertex a click lands
on, the colours a map is drawn with, and the readout of one vertex.

A PlanTUS run (position_transducer_plantus.m with save_metrics) writes plantus_metrics.mat
to its placement run folder: the skin surface (RAS mm, 0-based faces), the five criteria per
skin vertex, the avoidance mask, the composite score (lower is better, inf where not
eligible), its optimum and the vertex placed. The page draws these maps and lets the user
choose a vertex; placement.plantus.vertex then places it without recomputing the maps.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .model import ConfigModel
from .outputs import dir_placement_run

METRICS_FILENAME = "plantus_metrics.mat"
TX_CONFIG_FILENAME = "PlanTUSTxConfig.yaml"

# key in plantus_metrics.mat, label, unit, which end is better: 'low', 'high' or 'optimal'
# (the distance criterion: closeness to optimal_distance)
MAPS = (
    ("score", "Composite score", "", "low"),
    ("skin_target_distances", "Skin–target distance", "mm", "optimal"),
    ("skin_target_angles", "Beam angle at skin", "°", "low"),
    ("skin_target_intersections", "Beam–target overlap", "mm", "high"),
    ("skin_skull_angles", "Skin–skull angle", "°", "low"),
    ("skull_thickness", "Skull thickness", "mm", "low"),
)
MAP_LABELS = {k: label for k, label, _u, _b in MAPS}
EXCLUDED_RGBA = (0.62, 0.62, 0.62, 1.0)


@dataclass
class PlantusMaps:
    path: Path
    coords: np.ndarray                  # (N, 3) skin vertices, RAS mm
    faces: np.ndarray                   # (M, 3) 0-based
    target: np.ndarray                  # (3,) target centre, RAS mm
    maps: dict[str, np.ndarray]         # MAPS keys -> (N,)
    avoidance: np.ndarray               # (N,) 1 = may be used
    best: int | None                    # composite-score optimum
    selected: int | None                # vertex placed by that run (or a later re-placement)
    selection_mode: str
    min_distance: float
    max_distance: float
    optimal_distance: float
    weights: dict[str, float] = field(default_factory=dict)

    @property
    def n(self) -> int:
        return len(self.coords)

    @property
    def score(self) -> np.ndarray:
        return self.maps["score"]

    def eligible(self, v: int) -> bool:
        return bool(np.isfinite(self.score[v]))

    def rank(self, v: int) -> tuple[int, int] | None:
        """(rank, number eligible) under the composite score; None if not eligible."""
        s = self.score
        if not np.isfinite(s[v]):
            return None
        return int((s < s[v]).sum()) + 1, int(np.isfinite(s).sum())

    def readout(self, v: int) -> list[tuple[str, str]]:
        """(label, value) rows describing skin vertex v, as PlanTUS reports them."""
        if not 0 <= v < self.n:
            return [("Vertex", f"{v} is not on the skin surface ({self.n} vertices)")]

        def fmt(x, unit):
            return "n/a" if not np.isfinite(x) else f"{x:.1f} {unit}".rstrip()

        rows = [("Vertex", str(v))]
        for key, label, unit, _better in MAPS[1:]:
            rows.append((label, fmt(float(self.maps[key][v]), unit)))
        rk = self.rank(v)
        if rk is not None:
            rows.append(("Composite score", f"{self.score[v]:.3f} (rank {rk[0]} of {rk[1]})"))
        else:
            reasons = []
            if self.avoidance[v] == 0:
                reasons.append("avoidance area (eyes, ears, air cavities, face)")
            d = float(self.maps["skin_target_distances"][v])
            if d < self.min_distance or d > self.max_distance:
                reasons.append(f"outside the distance window {self.min_distance:.0f}–"
                               f"{self.max_distance:.0f} mm")
            rows.append(("Composite score", "not eligible" + (": " + "; ".join(reasons) if reasons else "")))
        return rows


def metrics_path(model: ConfigModel) -> Path | None:
    run = dir_placement_run(model, "plantus")
    return None if run is None else run / METRICS_FILENAME


def _scalar(d: dict, key: str, default=np.nan) -> float:
    v = d.get(key)
    if v is None:
        return default
    v = np.asarray(v).ravel()
    return float(v[0]) if v.size else default


def _vertex(d: dict, key: str) -> int | None:
    v = _scalar(d, key, -1)
    return None if not np.isfinite(v) or v < 0 else int(v)


def load(path: str | Path) -> PlantusMaps:
    """Read plantus_metrics.mat (one read of the file: it sits on the network share)."""
    import io

    from scipy.io import loadmat
    path = Path(path)
    d = loadmat(io.BytesIO(path.read_bytes()))
    n = len(d["skin_coordinates"])

    def vec(key):
        a = d.get(key)
        return np.full(n, np.nan) if a is None else np.asarray(a, float).ravel()

    maps = {k: vec(k) for k, *_ in MAPS}
    names = [str(np.asarray(x).ravel()[0]).strip() for x in np.asarray(d.get("weight_names", [])).ravel()]
    w = np.asarray(d.get("weights", []), float).ravel()
    mode = d.get("selection_mode")
    return PlantusMaps(
        path=path,
        coords=np.asarray(d["skin_coordinates"], float),
        faces=np.asarray(d["skin_faces"], np.int64),
        target=np.asarray(d.get("target_center", [np.nan] * 3), float).ravel()[:3],
        maps=maps,
        avoidance=vec("avoidance"),
        best=_vertex(d, "best_vertex"),
        selected=_vertex(d, "selected_vertex"),
        selection_mode="" if mode is None else str(np.asarray(mode).ravel()[0]),
        min_distance=_scalar(d, "min_distance"),
        max_distance=_scalar(d, "max_distance"),
        optimal_distance=_scalar(d, "optimal_distance"),
        weights=dict(zip(names, map(float, w))) if len(names) == len(w) else {},
    )


def read_tx_config(run_dir: str | Path) -> dict:
    """The PlanTUS transducer config the run was given (diameter, offsets); {} if absent."""
    import yaml
    p = Path(run_dir) / TX_CONFIG_FILENAME
    try:
        cfg = yaml.safe_load(p.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError):
        return {}
    return cfg if isinstance(cfg, dict) else {}


# ── colours ──────────────────────────────────────────────────────────────

def goodness(maps: PlantusMaps, key: str) -> np.ndarray:
    """Map values scaled to [0, 1] with 1 = best end of the criterion; NaN where the map
    has no value or the vertex lies in the avoidance area (PlanTUS masks its maps there)."""
    better = {k: b for k, _l, _u, b in MAPS}[key]
    x = np.array(maps.maps[key], float)
    if better == "optimal":
        x = np.abs(x - maps.optimal_distance)
        better = "low"
    x[~np.isfinite(x)] = np.nan
    x[maps.avoidance == 0] = np.nan
    ok = np.isfinite(x)
    if not ok.any():
        return x
    lo, hi = np.nanpercentile(x[ok], [2, 98])            # robust: a few outliers don't wash it out
    g = np.clip((x - lo) / (hi - lo), 0, 1) if hi > lo else np.where(ok, 0.5, np.nan)
    return 1 - g if better == "low" else g


def colours(good: np.ndarray, cmap: str = "viridis") -> np.ndarray:
    """(N, 4) float32 RGBA: the colour map over ``good`` (1 = bright), grey where NaN."""
    from matplotlib import colormaps
    rgba = np.asarray(colormaps[cmap](np.nan_to_num(good, nan=0.0)), np.float32)
    rgba[~np.isfinite(good)] = EXCLUDED_RGBA
    return rgba


# ── picking ──────────────────────────────────────────────────────────────

def ray_hit(origin, direction, verts: np.ndarray, faces: np.ndarray) -> tuple[int, np.ndarray] | None:
    """(face, point) where the ray first meets the mesh (Möller–Trumbore), or None."""
    o = np.asarray(origin, float)
    d = np.asarray(direction, float)
    d = d / np.linalg.norm(d)
    a, b, c = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    e1, e2 = b - a, c - a
    p = np.cross(d, e2)
    det = np.einsum("ij,ij->i", e1, p)
    ok = np.abs(det) > 1e-12
    inv = np.zeros_like(det)
    inv[ok] = 1.0 / det[ok]
    t0 = o - a
    u = np.einsum("ij,ij->i", t0, p) * inv
    q = np.cross(t0, e1)
    v = (q @ d) * inv
    t = np.einsum("ij,ij->i", e2, q) * inv
    hit = ok & (u >= 0) & (v >= 0) & (u + v <= 1) & (t > 0)
    if not hit.any():
        return None
    idx = np.flatnonzero(hit)
    k = idx[np.argmin(t[idx])]
    return int(k), o + t[k] * d


def pick_vertex(origin, direction, verts: np.ndarray, faces: np.ndarray) -> int | None:
    """The skin vertex nearest to where the ray first meets the surface, or None."""
    h = ray_hit(origin, direction, verts, faces)
    if h is None:
        return None
    face, point = h
    corners = faces[face]
    return int(corners[np.argmin(np.linalg.norm(verts[corners] - point, axis=1))])


# ── the transducer at a vertex ───────────────────────────────────────────

def preview(maps: PlantusMaps, v: int, tx: dict) -> dict:
    """Where PlanTUS puts the transducer for vertex v (its default 'target_centered'
    orientation): the beam axis from the skin vertex towards the target centre, the
    transducer centre (plane_offset + additional_offset) behind the skin."""
    skin = maps.coords[v]
    axis = maps.target - skin
    dist = float(np.linalg.norm(axis))
    axis = axis / dist if dist > 0 else -maps.coords[v] / max(np.linalg.norm(maps.coords[v]), 1e-9)
    back = float(tx.get("plane_offset", 0.0) or 0.0) + float(tx.get("additional_offset", 0.0) or 0.0)
    return {"skin": skin, "axis": axis, "centre": skin - back * axis, "distance": dist,
            "diameter": float(tx.get("transducer_diameter", 0.0) or 0.0),
            "skin_gap": float(tx.get("additional_offset", 0.0) or 0.0)}
