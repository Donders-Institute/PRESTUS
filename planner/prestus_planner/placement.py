"""Placement geometry: T1/segmentation access, PRESTUS voxel conventions,
transducer bowl geometry, scalp snapping and trajectory import.

Coordinates are subject RAS+ mm unless noted. PRESTUS voxel indices are
1-based MATLAB indices into the T1 array (niftiread order).

Conventions (checked against MATLAB R2023b): niftiinfo's Transform.T is the
raw sform/qform, i.e. it maps 0-based NIfTI indices to RAS, so MATLAB array
element p (1-based) sits at A @ (p - 1). PRESTUS' ras_to_grid / grid_to_ras
use the same rule.
"""

from __future__ import annotations

import csv
import math
import os
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# ── Volumes ──────────────────────────────────────────────────────────────

class Volume:
    """A NIfTI volume with RAS <-> voxel helpers (needs nibabel)."""

    def __init__(self, path: str | Path):
        import nibabel as nib
        self.path = Path(path)
        img = nib.load(str(self.path))
        self.affine = np.asarray(img.affine, dtype=float)       # 0-based voxel -> RAS
        self.inv_affine = np.linalg.inv(self.affine)
        self.shape = img.shape[:3]
        self._img = img
        self._data = None
        canon = nib.as_closest_canonical(img)
        self.canon_affine = np.asarray(canon.affine, dtype=float)
        self.canon_inv = np.linalg.inv(self.canon_affine)
        self._canon_img = canon
        self._canon = None

    @property
    def data(self) -> np.ndarray:
        if self._data is None:
            self._data = np.asanyarray(self._img.dataobj)
            if self._data.ndim > 3:
                self._data = self._data[..., 0]
        return self._data

    @property
    def canon(self) -> np.ndarray:
        """Data reoriented to the closest RAS axes (for display)."""
        if self._canon is None:
            d = np.asanyarray(self._canon_img.dataobj)
            self._canon = d[..., 0] if d.ndim > 3 else d
        return self._canon

    @property
    def voxel_size(self) -> np.ndarray:
        return np.sqrt((self.affine[:3, :3] ** 2).sum(axis=0))

    # 0-based voxel coordinates (float) <-> RAS mm
    def to_voxel(self, ras) -> np.ndarray:
        ras = np.atleast_2d(np.asarray(ras, dtype=float))
        return (self.inv_affine @ np.c_[ras, np.ones(len(ras))].T).T[:, :3]

    def to_ras(self, vox) -> np.ndarray:
        vox = np.atleast_2d(np.asarray(vox, dtype=float))
        return (self.affine @ np.c_[vox, np.ones(len(vox))].T).T[:, :3]

    def canon_to_ras(self, vox) -> np.ndarray:
        vox = np.atleast_2d(np.asarray(vox, dtype=float))
        return (self.canon_affine @ np.c_[vox, np.ones(len(vox))].T).T[:, :3]

    def ras_to_canon(self, ras) -> np.ndarray:
        ras = np.atleast_2d(np.asarray(ras, dtype=float))
        return (self.canon_inv @ np.c_[ras, np.ones(len(ras))].T).T[:, :3]

    def sample(self, ras) -> np.ndarray:
        """Nearest-neighbour values at RAS points (0 outside the volume)."""
        idx = np.rint(self.to_voxel(ras)).astype(int)
        inside = np.all((idx >= 0) & (idx < np.array(self.shape)), axis=1)
        out = np.zeros(len(idx), dtype=self.data.dtype)
        i = idx[inside]
        out[inside] = self.data[i[:, 0], i[:, 1], i[:, 2]]
        return out


# ── PRESTUS voxel conventions ────────────────────────────────────────────

def prestus_voxel_to_ras(vol: Volume, voxel) -> np.ndarray:
    """RAS of MATLAB array element ``voxel`` (1-based), i.e. where PRESTUS simulates."""
    return vol.to_ras(np.asarray(voxel, dtype=float) - 1)[0]


def ras_to_prestus_voxel(vol: Volume, ras) -> list[int]:
    """1-based array element containing ``ras`` (as ras_to_grid.m)."""
    return [int(v) for v in np.rint(vol.to_voxel(ras)[0]) + 1]


# ── Transducer geometry ──────────────────────────────────────────────────

@dataclass
class Bowl:
    curv_radius: float          # R [mm]
    aperture: float             # outer diameter D [mm]
    dist_geom_ep: float         # exit plane -> geometric focus [mm]

    @property
    def offset(self) -> float:
        """Bowl apex -> exit plane distance (focal_distance_offset)."""
        return self.curv_radius - self.dist_geom_ep

    @property
    def sagitta(self) -> float:
        r = self.aperture / 2
        return self.curv_radius - math.sqrt(max(self.curv_radius ** 2 - r ** 2, 0.0))

    @classmethod
    def from_transducer(cls, tx: dict) -> "Bowl | None":
        """Bowl from an effective transducer dict (annular or curved matrix)."""
        kind = tx.get("type")
        sub = tx.get(kind) or {} if kind in ("annular", "matrix") else {}
        R = _num(sub.get("curv_radius_mm"))
        if kind == "annular":
            od = sub.get("elem_od_mm")
            D = max(od) if isinstance(od, list) and od else _num(od)
        else:
            D = _num(sub.get("outer_diameter_mm"))
        if R is None or D is None or not math.isfinite(R):
            return None
        ep = _num(sub.get("dist_geom_ep_mm"))
        if ep is None:
            ep = math.sqrt(max(R ** 2 - (D / 2) ** 2, 0.0))     # validate_annular_transducer default
        return cls(R, D, ep)

    def outline(self, apex, axis, n: int = 72) -> dict[str, np.ndarray]:
        """3-D points for drawing: rim circle, exit plane centre, geometric focus."""
        apex = np.asarray(apex, float)
        axis = _unit(axis)
        u, v = _orthonormal(axis)
        t = np.linspace(0, 2 * np.pi, n)
        rim_c = apex + axis * self.sagitta
        rim = rim_c + (self.aperture / 2) * (np.outer(np.cos(t), u) + np.outer(np.sin(t), v))
        return {"rim": rim, "exit": apex + axis * self.offset,
                "geom_focus": apex + axis * self.curv_radius}


def _num(v) -> float | None:
    return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else None


def _unit(v) -> np.ndarray:
    v = np.asarray(v, float)
    n = np.linalg.norm(v)
    if n == 0:
        raise ValueError("zero-length direction")
    return v / n


def _orthonormal(axis: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    helper = np.array([1.0, 0, 0]) if abs(axis[0]) < 0.9 else np.array([0, 1.0, 0])
    u = np.cross(axis, helper)
    u /= np.linalg.norm(u)
    return u, np.cross(axis, u)


def place_on_axis(focus, towards_transducer, focal_distance_bowl: float) -> np.ndarray:
    """Bowl apex on the line from the focus, focal_distance_bowl away."""
    return np.asarray(focus, float) + _unit(towards_transducer) * focal_distance_bowl


def scalp_exit(seg: Volume, focus, towards_transducer, max_mm: float = 200, step: float = 0.25):
    """Last head voxel (label > 0) on the ray from the focus outwards, or None."""
    d = _unit(towards_transducer)
    s = np.arange(0, max_mm, step)
    pts = np.asarray(focus, float) + np.outer(s, d)
    inside = seg.sample(pts) > 0
    if not inside.any():
        return None
    return pts[np.nonzero(inside)[0][-1]]


def snap_to_scalp(seg: Volume, focus, towards_transducer, bowl: Bowl, gap_mm: float):
    """Apex such that the exit plane sits gap_mm outside the scalp along the axis."""
    exit_pt = scalp_exit(seg, focus, towards_transducer)
    if exit_pt is None:
        return None
    d = _unit(towards_transducer)
    return exit_pt + d * (gap_mm + bowl.offset)


# ── Trajectory files (Brainsight, 3D Slicer) ─────────────────────────────

@dataclass
class Trajectory:
    name: str
    target: np.ndarray            # RAS mm (focus)
    towards_transducer: np.ndarray  # unit vector target -> transducer
    source: str


def read_trajectory(path: str | Path) -> Trajectory:
    """Brainsight target export (.txt) or 3D Slicer ITK transform (.tfm/.txt).

    Follows BabelBrain (ConvMatTransform.py / BabelDatasetPreps.py): the
    4x4 location is the target and the transducer lies along the matrix's
    +z column. Slicer transforms are converted to that convention first.
    """
    path = Path(path)
    text = path.read_text(errors="ignore")
    if "Insight Transform File" in text or "Transform: AffineTransform" in text:
        return _read_slicer(text, path)
    if "brainsight" in text.lower() or "Loc. X" in text or _has_target_row(text):
        return _read_brainsight(text, path)
    raise ValueError(f"{path.name}: not a Brainsight export or Slicer transform")


def _has_target_row(text: str) -> bool:
    """A Brainsight-style data row: name followed by 12 numbers, tab separated."""
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        cols = line.split("\t")
        try:
            [float(c) for c in cols[1:13]]
        except ValueError:
            return False
        return len(cols) >= 13
    return False


def _to_trajectory(mat: np.ndarray, name: str, source: str) -> Trajectory:
    return Trajectory(name, mat[:3, 3].copy(), _unit(mat[:3, 2]), source)


def _read_brainsight(text: str, path: Path) -> Trajectory:
    header = [ln for ln in text.splitlines() if ln.startswith("#")]
    cs = next((ln.split(":", 1)[1].strip() for ln in header if "Coordinate system" in ln), "")
    if cs and "NIfTI" not in cs:
        raise ValueError(f"Brainsight coordinate system '{cs}' is not NIfTI; export in NIfTI space")
    rows = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    if not rows:
        raise ValueError("no target rows")
    cols = rows[0].split("\t")
    vals = [float(c) for c in cols[1:13]]
    mat = np.eye(4)
    mat[:3, 3] = vals[0:3]
    mat[:3, 0] = vals[3:6]
    mat[:3, 1] = vals[6:9]
    mat[:3, 2] = vals[9:12]
    return _to_trajectory(mat, cols[0], "Brainsight")


def _read_slicer(text: str, path: Path) -> Trajectory:
    params = next(ln for ln in text.splitlines() if ln.startswith("Parameters:"))
    fixed = next((ln for ln in text.splitlines() if ln.startswith("FixedParameters:")), "FixedParameters: 0 0 0")
    p = [float(x) for x in params.split(":", 1)[1].split()]
    c = np.array([float(x) for x in fixed.split(":", 1)[1].split()][:3])
    M = np.array(p[:9]).reshape(3, 3)
    itk = np.eye(4)
    itk[:3, :3] = M
    itk[:3, 3] = M @ (-c) + c + np.array(p[9:12])
    ras2lps = np.diag([-1.0, -1, 1, 1])
    bs = np.linalg.inv(ras2lps @ itk @ ras2lps)      # itk_to_BSight
    bs[:3, :3] = np.diag([-1.0, -1, -1]) @ bs[:3, :3]
    return _to_trajectory(bs, path.stem, "3D Slicer")


def find_t1(anat: str, pattern: str, subject_id: int) -> Path | None:
    from .paths import resolve_pattern
    hits = resolve_pattern(anat, pattern, subject_id)
    return hits[0] if hits else None


def find_segmentation(seg_root: str, subject_id: int) -> Path | None:
    p = Path(seg_root) / f"m2m_sub-{subject_id:03d}" / "final_tissues.nii.gz"
    return p if p.is_file() else None


# ── MNI <-> subject space (SimNIBS CLI) ──────────────────────────────────

class SimnibsCoordError(RuntimeError):
    """mni2subject_coords failed, or the tool or the m2m folder is missing."""


def mni_to_subject(mni_pt, m2m: Path, simnibs_bin: str, ld_library_path: str | None = None,
                   transformation_type: str = "nonl") -> np.ndarray:
    """Subject RAS+ mm for an MNI mm point.

    Shells out to SimNIBS's mni2subject_coords CLI tool via a CSV round-trip, exactly as
    mni2subject_coords_LDfix.m does -- no simnibs Python package needed, only
    startup.simnibs_bin_path and the subject's m2m folder. Has a large per-call overhead
    (a few seconds); call once per user action, never from a redraw/scroll handler.
    """
    exe = Path(simnibs_bin) / "mni2subject_coords"
    if not exe.is_file():
        raise SimnibsCoordError(f"mni2subject_coords not found in {simnibs_bin}")
    if not m2m.is_dir():
        raise SimnibsCoordError(f"{m2m} not found")
    env = os.environ.copy()
    if ld_library_path:
        env["LD_LIBRARY_PATH"] = ld_library_path
    with tempfile.TemporaryDirectory() as tmp:
        fin, fout = Path(tmp) / "in.csv", Path(tmp) / "out.csv"
        with fin.open("w", newline="") as f:
            csv.writer(f).writerow(["Generic", *(f"{c:.6f}" for c in mni_pt)])
        result = subprocess.run([str(exe), "-m", str(m2m), "-s", str(fin), "-o", str(fout),
                                 "-t", transformation_type], capture_output=True, text=True, env=env)
        if result.returncode != 0 or not fout.is_file():
            raise SimnibsCoordError((result.stderr or result.stdout).strip() or "mni2subject_coords failed")
        with fout.open(newline="") as f:
            row = next(csv.reader(f))
    return np.array([float(v) for v in row[1:4]])

