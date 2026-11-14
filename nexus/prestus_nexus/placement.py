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
    """A NIfTI volume with RAS <-> voxel helpers (needs nibabel).

    Only the header is read up front, from one read of the file's start. The voxels are
    read on first use in one go: nibabel's own .nii.gz reader opens the file several times
    and fetches it in small chunks, which on a network drive costs a round trip each
    (several seconds per volume)."""

    def __init__(self, path: str | Path):
        import nibabel as nib
        from nibabel.orientations import inv_ornt_aff, io_orientation
        self.path = Path(path)
        self._img = None
        hdr = _read_header(self.path)
        if hdr is None:                                         # not a plain NIfTI file
            self._img = nib.load(str(self.path))
            affine, shape = self._img.affine, self._img.shape
        else:
            affine, shape = hdr.get_best_affine(), hdr.get_data_shape()
        self.affine = np.asarray(affine, dtype=float)           # 0-based voxel -> RAS
        self.inv_affine = np.linalg.inv(self.affine)
        self.shape = tuple(int(n) for n in shape[:3])
        self._data = None
        # Closest-canonical (RAS) orientation, as nibabel.as_closest_canonical
        ornt = io_orientation(affine)
        self._ornt = None if np.array_equal(ornt, [[0, 1], [1, 1], [2, 1]]) else ornt
        canon_affine = affine if self._ornt is None else np.asarray(affine).dot(inv_ornt_aff(ornt, shape[:3]))
        self.canon_affine = np.asarray(canon_affine, dtype=float)
        self.canon_inv = np.linalg.inv(self.canon_affine)
        self._canon = None

    @property
    def data(self) -> np.ndarray:
        if self._data is None:
            self._data = _read_voxels(self.path, self._nib_image)
            if self._data.ndim > 3:
                self._data = self._data[..., 0]
        return self._data

    def _nib_image(self):
        """The nibabel image, opened only when the voxels are read from the original."""
        if self._img is None:
            import nibabel as nib
            self._img = nib.load(str(self.path))
        return self._img

    @property
    def canon(self) -> np.ndarray:
        """Data reoriented to the closest RAS axes (for display)."""
        if self._canon is None:
            if self._ornt is None:
                self._canon = self.data
            else:
                from nibabel.orientations import apply_orientation
                self._canon = apply_orientation(self.data, self._ornt)
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


VOXEL_CACHE_KEEP = 8          # decompressed volumes kept in the user's cache folder


def _voxel_cache_file(path: Path) -> Path | None:
    """Where the decompressed voxels of ``path`` are kept: named by its absolute path,
    size and modification time, so an edited or replaced image is read again."""
    import hashlib

    from .localcopy import cache_root
    try:
        st = path.stat()
    except OSError:
        return None
    tag = f"{path.absolute()}|{st.st_size}|{st.st_mtime_ns}"
    return cache_root() / "volumes" / (hashlib.sha1(tag.encode()).hexdigest()[:16] + ".npy")


def _read_header(path: Path):
    """The NIfTI-1/2 header of a single-file .nii or .nii.gz from one read of its first
    64 kB; None for anything else or when that is not enough (nibabel then opens it)."""
    import io
    import zlib

    import nibabel as nib
    name = path.name.lower()
    if not (name.endswith(".nii") or name.endswith(".nii.gz")):
        return None
    try:
        with open(path, "rb") as f:
            head = f.read(1 << 16)
        if name.endswith(".gz"):
            head = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(head, 1 << 16)
        for cls in (nib.Nifti1Header, nib.Nifti2Header):
            if len(head) < cls.sizeof_hdr:
                continue
            size = (int.from_bytes(head[:4], "little"), int.from_bytes(head[:4], "big"))
            if cls.sizeof_hdr not in size:
                continue
            hdr = cls.from_fileobj(io.BytesIO(head), check=False)
            if hdr["magic"].item() not in (b"n+1", b"n+2"):
                return None
            return hdr
    except Exception:
        return None
    return None


def _read_voxels(path: Path, load_img) -> np.ndarray:
    """The image's (scaled) voxel array. The project volume is slow for whole images
    (several seconds each), so the array is also kept in the user's cache folder and
    read from there while the image is unchanged (the original stays the only source).
    load_img() returns the nibabel image; it is called only when the cache misses."""
    cached = _voxel_cache_file(path)
    if cached is not None and cached.is_file():
        try:
            data = np.load(cached, allow_pickle=False)
            os.utime(cached)                      # most recently used: kept longest
            return data
        except Exception:
            pass
    data = _read_voxels_from(path, load_img())
    if cached is not None:
        try:
            cached.parent.mkdir(parents=True, exist_ok=True)
            tmp = cached.with_suffix(".tmp.npy")
            np.save(tmp, np.asarray(data), allow_pickle=False)
            os.replace(tmp, cached)
            old = sorted(cached.parent.glob("*.npy"), key=lambda p: p.stat().st_mtime, reverse=True)
            for stale in old[VOXEL_CACHE_KEEP:]:
                stale.unlink(missing_ok=True)
        except Exception:
            pass
    return data


def _read_voxels_from(path: Path, img) -> np.ndarray:
    """The image's (scaled) voxel array; a .nii.gz is read whole and decompressed in
    memory, falling back to nibabel's reader for anything that fails."""
    if path.name.endswith(".gz"):
        try:
            import gzip
            whole = type(img).from_bytes(gzip.decompress(path.read_bytes()))
            return np.asanyarray(whole.dataobj)
        except Exception:
            pass
    return np.asanyarray(img.dataobj)


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
            if sub.get("is_curved") in (0, False):       # flat (validate_matrix_transducer)
                return None
            D = _num(sub.get("outer_diameter_mm"))
        if R is None or D is None or not math.isfinite(R):
            return None
        ep = _num(sub.get("dist_geom_ep_mm"))
        if ep is None:
            # exit plane at the bowl rim, as validate_annular/matrix_transducer.m
            ep = math.sqrt(max(R ** 2 - (D / 2) ** 2, 0.0))
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


def fov_box(trans, focus, diameter_mm: float, length_mm: float) -> dict[str, np.ndarray]:
    """Beam-aligned box for the acoustic FOV crop (doc_limited_fov.md): axial extent
    length_mm starting at the transducer centre, lateral extent diameter_mm (full
    width) centred on the focus. Returns the 8 corners and the 12 edges (each a (2,3)
    point pair) in RAS mm, for drawing -- a projection of the 3-D box, not a true
    cross-section, same convention as Bowl.outline."""
    trans = np.asarray(trans, float)
    focus = np.asarray(focus, float)
    axis = _unit(focus - trans)
    u, v = _orthonormal(axis)
    r = diameter_mm / 2
    corners = np.array([trans + axis * s + u * r * su + v * r * sv
                        for s in (0.0, length_mm) for su, sv in ((1, 1), (1, -1), (-1, -1), (-1, 1))])
    near, far = corners[:4], corners[4:]
    edges = [near[[i, (i + 1) % 4]] for i in range(4)]
    edges += [far[[i, (i + 1) % 4]] for i in range(4)]
    edges += [np.vstack([near[i], far[i]]) for i in range(4)]
    return {"corners": corners, "edges": edges}


def resolve_focal_distance_ep(tx: dict, bowl: "Bowl | None") -> float | None:
    """focal_distance_ep (exit plane -> focus) for an effective transducer dict: the
    stored value, or derived from focal_distance_bowl (apex -> focus) via bowl.offset
    when only that is set. The inverse direction of the bowl-apex placement used by
    PlacementPage._expected_bowl_distance (which resolves the apex distance instead)."""
    ep = _num(tx.get("focal_distance_ep"))
    if ep is not None:
        return ep
    fb = _num(tx.get("focal_distance_bowl"))
    if fb is not None and bowl is not None:
        return fb - bowl.offset
    return None


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


# ── Beam QC (same measures as the placement comparison) ──────────────────

BONE_LABELS = (4, 7, 8)       # SimNIBS charm: bone, compact bone, spongy bone
EYE_LABEL = 6
INTRACRANIAL_LABELS = (1, 2, 3)  # WM, GM, CSF: inside the skull


def _inward_normal(seg: Volume, p: np.ndarray, labels, radius: float = 5.0) -> np.ndarray | None:
    """Unit normal of the tissue surface at p, pointing into the tissue: from the
    centroid of the surrounding non-tissue samples to that of the tissue samples."""
    g = np.arange(-radius, radius + 1e-9, 1.0)
    off = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3)
    off = off[np.linalg.norm(off, axis=1) <= radius]
    inside = np.isin(seg.sample(p + off), labels)
    if inside.all() or not inside.any():
        return None
    n = off[inside].mean(axis=0) - off[~inside].mean(axis=0)
    return n / np.linalg.norm(n) if np.linalg.norm(n) > 0 else None


def beam_qc(seg: Volume, trans, focus, bowl: Bowl, step: float = 0.25) -> dict[str, float]:
    """Geometry of a placement along its beam, from the tissue labels
    (code/report/02_placement/placement_geometry_metrics.m uses the same measures):

    skull_path_mm               bone on the central ray, exit plane -> focus
    skull_incidence_deg         beam vs skull normal where the central ray enters bone
    beam_skull_path_mean/sd_mm  bone on 37 rays from the aperture (centre + 3 rings x 12)
    beam_skull_incidence_mean_deg   to the focus, and their mean incidence angle
    rim_gap_mean/sd/min_mm      exit-plane rim -> scalp along the beam (36 points);
                                SD = coupling evenness, min < 0 = the rim cuts into the head
    eye_in_beam                 fraction of the 37 rays that pass an eye
    Values are NaN where they cannot be measured. Incidence uses the outer skull
    normal from a 5 mm neighbourhood; it can differ by a few degrees from the
    MATLAB evaluator, which takes the gradient of a smoothed bone mask."""
    trans, focus = np.asarray(trans, float), np.asarray(focus, float)
    axis = _unit(focus - trans)
    ep = trans + axis * bowl.offset
    u, v = _orthonormal(axis)
    nan = float("nan")
    out: dict[str, float] = {}

    def ray(origin, target):
        d = target - origin
        L = float(np.linalg.norm(d))
        d = d / L
        s = np.arange(0, L, step)
        pts = origin + np.outer(s, d)
        lab = seg.sample(pts)
        bone = np.isin(lab, BONE_LABELS)
        inc = nan
        if bone.any():
            # outer skull surface: bone and everything inside it count as inside,
            # so a thin skull does not see the brain side as 'outside'
            n = _inward_normal(seg, pts[np.argmax(bone)], BONE_LABELS + INTRACRANIAL_LABELS)
            if n is not None:
                inc = float(np.degrees(np.arccos(min(1.0, abs(float(n @ d))))))
        return step * float(bone.sum()), inc, bool((lab == EYE_LABEL).any())

    out["skull_path_mm"], out["skull_incidence_deg"], _ = ray(ep, focus)

    R = bowl.aperture / 2
    origins = [ep]
    for r in (R / 3, 2 * R / 3, R):
        for a in np.radians(np.arange(0, 360, 30)):
            origins.append(ep + r * (np.cos(a) * u + np.sin(a) * v))
    res = [ray(o, focus) for o in origins]
    paths = np.array([r[0] for r in res])
    incs = np.array([r[1] for r in res])
    out["beam_skull_path_mean_mm"] = float(paths.mean())
    out["beam_skull_path_sd_mm"] = float(paths.std(ddof=1))
    out["beam_skull_incidence_mean_deg"] = float(np.nanmean(incs)) if np.isfinite(incs).any() else nan
    out["eye_in_beam"] = float(np.mean([r[2] for r in res]))

    back = 30.0
    t = np.arange(-back, 60, step)
    gaps = []
    for a in np.radians(np.arange(0, 360, 10)):
        rim = ep + R * (np.cos(a) * u + np.sin(a) * v)
        head = seg.sample(rim + np.outer(t, axis)) > 0
        gaps.append(float(t[np.argmax(head)]) if head.any() else nan)
    gaps = np.array(gaps)
    ok = np.isfinite(gaps)
    out["rim_gap_mean_mm"] = float(gaps[ok].mean()) if ok.any() else nan
    out["rim_gap_sd_mm"] = float(gaps[ok].std(ddof=1)) if ok.sum() > 1 else nan
    out["rim_gap_min_mm"] = float(gaps[ok].min()) if ok.any() else nan
    return out


# ── Beam plane: a 2-D section that contains the beam axis ────────────────

def beam_plane_basis(trans, focus, angle_deg: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """(axis, lateral) unit vectors of the plane through the beam axis. At 0° the lateral
    direction points as close to superior as the axis allows (anterior for a near-vertical
    beam); angle_deg turns the plane about the axis."""
    axis = _unit(np.asarray(focus, float) - np.asarray(trans, float))
    ref = np.array([0, 0, 1.0]) if abs(axis[2]) < 0.9 else np.array([0, 1.0, 0])
    lat0 = _unit(ref - axis * float(ref @ axis))
    other = np.cross(axis, lat0)
    a = math.radians(angle_deg)
    return axis, math.cos(a) * lat0 + math.sin(a) * other


def beam_plane_points(trans, axis, lateral, depth: np.ndarray, offsets: np.ndarray) -> np.ndarray:
    """RAS points of the plane, shape (len(offsets), len(depth), 3): row i is the line
    offsets[i] mm to the side of the axis, column j depth[j] mm from the bowl apex."""
    d, l = np.meshgrid(depth, offsets)
    return np.asarray(trans, float) + d[..., None] * axis + l[..., None] * lateral


def _label_list(v) -> list[int]:
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        return [int(v)]
    return [int(x) for x in v if isinstance(x, (int, float)) and not isinstance(x, bool)] \
        if isinstance(v, list) else []


def segmentation_to_medium(labels: np.ndarray, layers: dict, medium_names: list[str]) -> np.ndarray:
    """SimNIBS labels -> medium index (1-based position in medium_names), the mapping of
    preproc_medium_mask.m without its smoothing: layers in config order, later ones win;
    with both skull_cortical and skull_trabecular, every skull label is cortical and
    trabecular overlays it. Unmapped voxels are water, as medium_setup.m does."""
    out = np.zeros(np.shape(labels), dtype=int)
    names = [n for n in layers if n in medium_names]
    multi = "skull_cortical" in names and "skull_trabecular" in names
    skulls = ("skull", "skull_cortical", "skull_trabecular")
    for n in names:
        if not (multi and n in skulls):
            out[np.isin(labels, _label_list(layers[n]))] = medium_names.index(n) + 1
    if multi:
        shell = sum((_label_list(layers.get(n)) for n in skulls), [])
        out[np.isin(labels, shell)] = medium_names.index("skull_cortical") + 1
        out[np.isin(labels, _label_list(layers["skull_trabecular"]))] = \
            medium_names.index("skull_trabecular") + 1
    if "water" in medium_names:
        out[out == 0] = medium_names.index("water") + 1
    return out


def medium_clearance(medium: np.ndarray, depth: np.ndarray, offsets: np.ndarray, bowl: Bowl,
                     water_index: int) -> tuple[float, float]:
    """(on the axis, smallest across the aperture) distance from the exit plane to the first
    non-water medium voxel in a beam-plane section [mm], searched from the bowl surface
    forwards. Negative: tissue sits behind the exit plane, inside the transducer's space.
    0 (outside the simulation grid) counts as no tissue. NaN where no tissue is met."""
    R, a = bowl.curv_radius, bowl.aperture / 2
    tissue = (medium != water_index) & (medium != 0)
    nan = float("nan")

    def gap(i: int) -> float:
        y = offsets[i]
        sag = R - math.sqrt(max(R * R - y * y, 0.0))
        hit = tissue[i] & (depth >= sag)
        return float(depth[np.argmax(hit)]) - bowl.offset if hit.any() else nan

    rows = [i for i, y in enumerate(offsets) if abs(y) <= a]
    gaps = np.array([gap(i) for i in rows])
    centre = gap(int(np.argmin(np.abs(offsets))))
    return centre, float(np.nanmin(gaps)) if np.isfinite(gaps).any() else nan


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

