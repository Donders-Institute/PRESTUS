"""Triangle meshes for the 3D placement view (numpy and scipy only, no Qt).

Every function returns ``(vertices, faces)``: vertices as an (N, 3) float array in subject
RAS mm, faces as an (M, 3) int array of vertex indices, wound counter-clockwise seen from
outside (normals from ``vertex_normals`` point outwards).

Transducer geometry follows PRESTUS (doc_transducer.md): the bowl apex is
``transducer.trans_pos``, the beam axis points from the apex to the focus, ring ``i`` covers
the bowl between the inner and outer diameters ``elem_id_mm[i]`` and ``elem_od_mm[i]``
(measured across the aperture, as makeAnnularArray does), and the exit plane lies
``curv_radius_mm - dist_geom_ep_mm`` in front of the apex.
"""

from __future__ import annotations

import math

import numpy as np

from .placement import Bowl, Volume, _orthonormal, _unit

Mesh = tuple[np.ndarray, np.ndarray]


# ── Surfaces from label volumes ──────────────────────────────────────────

def voxel_surface(mask: np.ndarray, affine: np.ndarray, smooth_iter: int = 12) -> Mesh:
    """Outer boundary of a binary volume as a smoothed triangle mesh in world coordinates.

    The boundary faces between inside and outside voxels form a closed, consistently
    oriented surface; Taubin smoothing (which does not shrink it) then rounds off the
    voxel steps. ``affine`` maps the mask's 0-based voxel indices to RAS mm."""
    m = np.pad(np.asarray(mask, bool), 1)
    quads = []
    for a in range(3):
        b, c = (a + 1) % 3, (a + 2) % 3
        d = np.diff(m.astype(np.int8), axis=a)       # +1: outside -> inside along +a
        for sign in (-1, 1):
            idx = np.argwhere(d == sign)            # voxel i (and i+1 along a)
            if not len(idx):
                continue
            base = idx.copy()
            base[:, a] += 1                         # the face lies on corner plane i+1
            corners = []
            for db, dc in ((0, 0), (1, 0), (1, 1), (0, 1)):    # normal along +a
                q = base.copy()
                q[:, b] += db
                q[:, c] += dc
                corners.append(q)
            if sign == 1:                           # inside lies on the +a side: flip
                corners = corners[::-1]
            quads.append(np.stack(corners, axis=1))
    if not quads:
        return np.empty((0, 3)), np.empty((0, 3), int)
    quads = np.concatenate(quads)                   # (n, 4, 3) corner grid indices
    grid = np.array(m.shape) + 1
    flat = np.ravel_multi_index(quads.reshape(-1, 3).T, grid)
    uniq, inverse = np.unique(flat, return_inverse=True)
    q = inverse.reshape(-1, 4)
    faces = np.concatenate([q[:, [0, 1, 2]], q[:, [0, 2, 3]]])
    corner = np.stack(np.unravel_index(uniq, grid), axis=1).astype(float)
    verts = corner - 1.5                            # corner c sits at voxel c - 1.5 (1 pad)
    verts = taubin_smooth(verts, faces, smooth_iter)
    world = (np.c_[verts, np.ones(len(verts))] @ np.asarray(affine, float).T)[:, :3]
    if np.linalg.det(np.asarray(affine, float)[:3, :3]) < 0:
        faces = faces[:, ::-1]                      # a mirroring affine flips the winding
    return world, faces


def mask_surface(data: np.ndarray, affine: np.ndarray, smooth_iter: int = 4) -> Mesh:
    """Surface of a binary mask volume (> 0.5) in RAS mm, built on its bounding box only, so
    a small target mask in a whole-head volume costs little. Few smoothing passes: a mask of
    a few voxels would otherwise shrink to a blob."""
    m = np.asarray(data) > 0.5
    idx = np.argwhere(m)
    if not len(idx):
        return np.empty((0, 3)), np.empty((0, 3), int)
    lo, hi = idx.min(axis=0), idx.max(axis=0) + 1
    sub = m[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]
    shift = np.eye(4)
    shift[:3, 3] = lo
    return voxel_surface(sub, np.asarray(affine, float) @ shift, smooth_iter)


def taubin_smooth(verts: np.ndarray, faces: np.ndarray, iterations: int,
                  lam: float = 0.5, mu: float = -0.53) -> np.ndarray:
    """Taubin lambda/mu smoothing with uniform neighbour weights."""
    if iterations <= 0 or not len(faces):
        return verts
    from scipy import sparse
    n = len(verts)
    i = np.concatenate([faces[:, 0], faces[:, 1], faces[:, 2], faces[:, 1], faces[:, 2], faces[:, 0]])
    j = np.concatenate([faces[:, 1], faces[:, 2], faces[:, 0], faces[:, 0], faces[:, 1], faces[:, 2]])
    adj = sparse.coo_matrix((np.ones(len(i)), (i, j)), shape=(n, n)).tocsr()
    adj.data[:] = 1.0                               # duplicate edges count once
    deg = np.asarray(adj.sum(axis=1)).ravel()
    deg[deg == 0] = 1
    v = verts.copy()
    for _ in range(iterations):
        for f in (lam, mu):
            v = v + f * (adj @ v / deg[:, None] - v)
    return v


def head_surfaces(seg: Volume, target_mm: float = 2.0) -> dict[str, Mesh]:
    """Scalp (every head label, holes filled) and brain (white + grey matter) surfaces of a
    SimNIBS charm segmentation, on a grid of about ``target_mm``."""
    from scipy import ndimage
    labels = seg.data
    step = max(1, int(round(target_mm / float(np.min(seg.voxel_size)))))
    sub = labels[::step, ::step, ::step]
    affine = seg.affine @ np.diag([step, step, step, 1.0])
    out = {}
    head = ndimage.binary_fill_holes(sub > 0)
    head = _largest_component(head)
    out["scalp"] = voxel_surface(head, affine)
    brain = ndimage.binary_fill_holes(np.isin(sub, (1, 2)))
    out["brain"] = voxel_surface(_largest_component(brain), affine) if brain.any() \
        else (np.empty((0, 3)), np.empty((0, 3), int))
    face = face_region(sub, affine)
    out["scalp_defaced"] = scramble_region(*out["scalp"], *face) if face is not None else out["scalp"]
    return out


FACE_MARGIN_MM = 15.0       # the face region reaches this far behind and above the eyes' centre
SCRAMBLE_ITERATIONS = 40    # smoothing passes over the face: removes nose, lips and brows
SCRAMBLE_NOISE_MM = 1.5     # then random displacement (SD), seeded: the same head looks the same


def face_region(labels: np.ndarray, affine: np.ndarray) -> tuple[float, float] | None:
    """(A from, S up to) [RAS mm] of the face, for defacing the scalp: in front of a plane
    FACE_MARGIN_MM behind the eyes' centre and below one FACE_MARGIN_MM above it (eyes,
    brows, nose, mouth and chin). Without an eye label (SimNIBS 6), from the brain (labels
    1, 2): in front of its anterior pole less the margin, below its centre. Assumes the head
    roughly upright in RAS, as a T1 is acquired; None without eyes or brain."""
    from .placement import EYE_LABEL
    for sel, anchor in ((labels == EYE_LABEL, "eyes"), (np.isin(labels, (1, 2)), "brain")):
        idx = np.argwhere(sel)
        if len(idx):
            break
    else:
        return None
    ras = (affine @ np.c_[idx, np.ones(len(idx))].T).T[:, :3]
    if anchor == "eyes":
        c = ras.mean(axis=0)
        return float(c[1] - FACE_MARGIN_MM), float(c[2] + FACE_MARGIN_MM)
    return float(ras[:, 1].max() - FACE_MARGIN_MM), float(ras[:, 2].mean())


def scramble_region(verts: np.ndarray, faces: np.ndarray, a_from: float, s_to: float) -> Mesh:
    """The mesh with the face scrambled rather than cut out: vertices in the box A > a_from,
    S < s_to are smoothed towards their neighbours (SCRAMBLE_ITERATIONS passes; the box edge
    stays put, so the surface stays closed) and then displaced at random (SCRAMBLE_NOISE_MM).
    The head keeps its outline; the facial features are not recoverable from the view."""
    from scipy import sparse
    verts = np.array(verts, float)
    inside = (verts[:, 1] > a_from) & (verts[:, 2] < s_to)
    if not inside.any() or not len(faces):
        return verts, faces
    n = len(verts)
    i = np.concatenate([faces[:, 0], faces[:, 1], faces[:, 2], faces[:, 1], faces[:, 2], faces[:, 0]])
    j = np.concatenate([faces[:, 1], faces[:, 2], faces[:, 0], faces[:, 0], faces[:, 1], faces[:, 2]])
    adj = sparse.coo_matrix((np.ones(len(i)), (i, j)), shape=(n, n)).tocsr()
    adj.data[:] = 1.0
    deg = np.asarray(adj.sum(axis=1)).ravel()
    deg[deg == 0] = 1.0
    for _ in range(SCRAMBLE_ITERATIONS):
        mean = adj @ verts / deg[:, None]
        verts[inside] = mean[inside]
    rng = np.random.default_rng(int(n))
    verts[inside] += rng.normal(0.0, SCRAMBLE_NOISE_MM, (int(inside.sum()), 3))
    return verts, faces


def without_region(verts: np.ndarray, faces: np.ndarray, a_from: float, s_to: float) -> Mesh:
    """The mesh without the triangles that touch the box A > a_from and S < s_to."""
    inside = (verts[:, 1] > a_from) & (verts[:, 2] < s_to)
    keep = ~inside[faces].any(axis=1)
    return verts, faces[keep]


def _largest_component(mask: np.ndarray) -> np.ndarray:
    from scipy import ndimage
    lab, n = ndimage.label(mask)
    if n <= 1:
        return mask
    sizes = ndimage.sum(mask, lab, range(1, n + 1))
    return lab == (int(np.argmax(sizes)) + 1)


# ── Transducer ───────────────────────────────────────────────────────────

def ring_meshes(bowl: Bowl, elem_id_mm, elem_od_mm, apex, axis,
                n_phi: int = 96, n_rad: int = 6) -> list[Mesh]:
    """One mesh per element ring on the bowl (a flat transducer when the radius of
    curvature is infinite). Without ring diameters, one ring covers the aperture."""
    ids = list(elem_id_mm) if elem_id_mm is not None else []
    ods = list(elem_od_mm) if elem_od_mm is not None else []
    if not ods or len(ids) != len(ods):
        ids, ods = [0.0], [bowl.aperture]
    out = []
    for d_in, d_out in zip(ids, ods):
        r = np.linspace(float(d_in) / 2, float(d_out) / 2, n_rad + 1)
        out.append(_world(*_annulus(r, bowl.curv_radius, n_phi), apex, axis))
    return out


def _annulus(radii: np.ndarray, curv_radius: float, n_phi: int) -> Mesh:
    """Annulus on the bowl in the local frame: apex at the origin, beam along +z.
    ``radii`` are distances from the axis (half the ring diameters)."""
    phi = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
    rr, pp = np.meshgrid(radii, phi, indexing="ij")
    if curv_radius is None or not math.isfinite(curv_radius):
        z = np.zeros_like(rr)
    else:
        R = float(curv_radius)
        z = R - np.sqrt(np.maximum(R ** 2 - rr ** 2, 0.0))       # sagitta at radius r
    verts = np.stack([rr * np.cos(pp), rr * np.sin(pp), z], -1).reshape(-1, 3)
    n_r = len(radii)
    faces = []
    for i in range(n_r - 1):
        for k in range(n_phi):
            a, b = i * n_phi + k, i * n_phi + (k + 1) % n_phi
            c, d = (i + 1) * n_phi + (k + 1) % n_phi, (i + 1) * n_phi + k
            faces += [(a, d, c), (a, c, b)]          # normal along +z: the radiating side
    return verts, np.array(faces, int).reshape(-1, 3)


def exit_plane_disc(bowl: Bowl, apex, axis, n_phi: int = 96) -> Mesh:
    """The exit plane as a disc of the aperture's diameter."""
    phi = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
    r = bowl.aperture / 2
    rim = np.stack([r * np.cos(phi), r * np.sin(phi), np.full(n_phi, bowl.offset)], -1)
    verts = np.vstack([[0, 0, bowl.offset], rim])
    faces = np.array([(0, 1 + (k + 1) % n_phi, 1 + k) for k in range(n_phi)], int)
    return _world(verts, faces, apex, axis)


GENERIC_HOUSING_DEPTH_MM = 15.0     # depth of the generic housing behind the bowl apex
HOUSING_BACK_MARGIN_MM = 2.0        # a housing reaches at least this far behind the bowl apex


def housing_size(aperture: float, offset: float,
                 housing: dict | None) -> tuple[float, float, float, bool]:
    """(diameter, depth [mm], front face [mm from the bowl apex along the beam], known) of
    the housing of a transducer with this aperture [mm] and apex-to-exit-plane ``offset`` [mm].
    The front face is the exit plane, so the housing holds the elements. With a size in the
    equipment entry's 'housing' block, that size, made at least as wide as the aperture and deep
    enough to reach behind the bowl apex. Without one, a generic cylinder: the aperture's
    diameter, from the exit plane to GENERIC_HOUSING_DEPTH_MM behind the apex."""
    h = housing if isinstance(housing, dict) else {}
    d, depth = h.get("diameter_mm"), h.get("depth_mm")
    front = float(offset)
    if isinstance(d, (int, float)) and isinstance(depth, (int, float)) and d > 0 and depth > 0:
        return (max(float(d), float(aperture)), max(float(depth), front + HOUSING_BACK_MARGIN_MM),
                front, True)
    return float(aperture), front + GENERIC_HOUSING_DEPTH_MM, front, False


def housing_puck(diameter: float, depth: float, front: float, apex, axis, n_phi: int = 96) -> Mesh:
    """The housing as a closed cylinder: front face ``front`` mm
    from the bowl apex along the beam, back face ``depth`` mm behind it."""
    phi = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
    ring = np.stack([np.cos(phi), np.sin(phi)], -1) * diameter / 2
    zf, zb = front, front - depth
    verts = np.vstack([np.c_[ring, np.full(n_phi, zf)], np.c_[ring, np.full(n_phi, zb)],
                       [[0.0, 0.0, zf], [0.0, 0.0, zb]]])
    cf, cb = 2 * n_phi, 2 * n_phi + 1
    faces = []
    for k in range(n_phi):
        a, b = k, (k + 1) % n_phi
        c, d = n_phi + b, n_phi + a
        faces += [(a, d, c), (a, c, b), (cf, a, b), (cb, c, d)]   # wall, front, back: outward
    return _world(verts, np.array(faces, int), apex, axis)


def _world(verts: np.ndarray, faces: np.ndarray, apex, axis) -> Mesh:
    """Local transducer frame (apex at the origin, beam along +z) -> RAS mm."""
    axis = _unit(axis)
    u, v = _orthonormal(axis)
    rot = np.stack([u, v, axis], axis=1)           # columns: local x, y, z in world
    if np.linalg.det(rot) < 0:
        rot[:, 0] = -rot[:, 0]                      # keep the frame right-handed
    return np.asarray(apex, float) + verts @ rot.T, faces


# ── Markers ──────────────────────────────────────────────────────────────

def sphere(centre, radius: float, n: int = 16) -> Mesh:
    theta = np.linspace(0, np.pi, n + 1)
    phi = np.linspace(0, 2 * np.pi, 2 * n, endpoint=False)
    tt, pp = np.meshgrid(theta, phi, indexing="ij")
    verts = radius * np.stack([np.sin(tt) * np.cos(pp), np.sin(tt) * np.sin(pp), np.cos(tt)], -1).reshape(-1, 3)
    m = 2 * n
    faces = []
    for i in range(n):
        for k in range(m):
            a, b = i * m + k, i * m + (k + 1) % m
            c, d = (i + 1) * m + (k + 1) % m, (i + 1) * m + k
            faces += [(a, d, c), (a, c, b)]
    return np.asarray(centre, float) + verts, np.array(faces, int)


def cylinder(p0, p1, radius: float, n: int = 12) -> Mesh:
    """Open tube from p0 to p1 (for the beam axis)."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    length = float(np.linalg.norm(p1 - p0))
    if length <= 0:
        return np.empty((0, 3)), np.empty((0, 3), int)
    phi = np.linspace(0, 2 * np.pi, n, endpoint=False)
    ring = np.stack([radius * np.cos(phi), radius * np.sin(phi)], -1)
    verts = np.vstack([np.c_[ring, np.zeros(n)], np.c_[ring, np.full(n, length)]])
    faces = []
    for k in range(n):
        a, b, c, d = k, (k + 1) % n, n + (k + 1) % n, n + k
        faces += [(a, b, c), (a, c, d)]
    return _world(verts, np.array(faces, int), p0, p1 - p0)


def vertex_normals(verts: np.ndarray, faces: np.ndarray) -> np.ndarray:
    """Area-weighted vertex normals (unit length; zero for unused vertices)."""
    n = np.zeros_like(verts, dtype=float)
    if len(faces):
        fn = np.cross(verts[faces[:, 1]] - verts[faces[:, 0]], verts[faces[:, 2]] - verts[faces[:, 0]])
        for k in range(3):
            np.add.at(n, faces[:, k], fn)
    norm = np.linalg.norm(n, axis=1, keepdims=True)
    return np.divide(n, norm, out=np.zeros_like(n), where=norm > 0)
