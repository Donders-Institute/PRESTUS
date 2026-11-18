"""Placement geometry, PRESTUS voxel conventions and trajectory import."""

import math

import numpy as np
import pytest

nib = pytest.importorskip("nibabel")

from prestus_nexus import placement as pl  # noqa: E402

AFFINE = np.array([[-1.1, 0.05, 0, 90.3], [0.02, 1.2, 0.1, -126.7],
                   [0, -0.08, 0.9, -72.2], [0, 0, 0, 1]])


@pytest.fixture
def vol(tmp_path):
    data = np.zeros((40, 50, 60), np.float32)
    data[10, 20, 30] = 1
    path = tmp_path / "t1.nii.gz"
    nib.save(nib.Nifti1Image(data, AFFINE), str(path))
    return pl.Volume(path)


def test_voxel_conventions(vol):
    ras = (AFFINE @ [10, 20, 30, 1])[:3]           # centre of the bright voxel
    # MATLAB element (1-based) holding it is [11 21 31] (checked in R2023b)
    assert pl.ras_to_prestus_voxel(vol, ras) == [11, 21, 31]
    np.testing.assert_allclose(pl.prestus_voxel_to_ras(vol, [11, 21, 31]), ras)
    assert vol.sample([ras])[0] == 1


def test_canonical_round_trip(vol):
    ras = np.array([[12.3, -40.1, 5.5]])
    np.testing.assert_allclose(vol.canon_to_ras(vol.ras_to_canon(ras)), ras, atol=1e-9)


def test_bowl_geometry():
    tx = {"type": "annular", "annular": {"curv_radius_mm": 63.2, "dist_geom_ep_mm": 52.38,
                                         "elem_od_mm": [32.3, 45.4, 55.4, 64]}}
    bowl = pl.Bowl.from_transducer(tx)
    assert bowl.aperture == 64 and math.isclose(bowl.offset, 63.2 - 52.38)
    out = bowl.outline([0, 0, 0], [0, 0, 1])
    np.testing.assert_allclose(out["geom_focus"], [0, 0, 63.2])
    radii = np.linalg.norm(out["rim"][:, :2], axis=1)
    np.testing.assert_allclose(radii, 32, atol=1e-9)
    np.testing.assert_allclose(out["rim"][:, 2], bowl.sagitta, atol=1e-9)
    assert pl.Bowl.from_transducer({"type": "annular", "annular": {}}) is None


def test_bowl_without_dist_geom_ep_uses_rim():
    # as validate_annular/matrix_transducer.m: exit plane at the rim of the aperture
    ann = pl.Bowl.from_transducer({"type": "annular", "annular": {
        "curv_radius_mm": 63.2, "elem_od_mm": [32.3, 45.4, 55.4, 64]}})
    assert ann.offset == pytest.approx(ann.sagitta)
    assert ann.dist_geom_ep == pytest.approx(math.sqrt(63.2 ** 2 - 32 ** 2))
    mat = pl.Bowl.from_transducer({"type": "matrix", "matrix": {
        "is_curved": 1, "curv_radius_mm": 72, "outer_diameter_mm": 60}})
    assert mat.offset == pytest.approx(mat.sagitta) and mat.aperture == 60
    assert pl.Bowl.from_transducer({"type": "matrix", "matrix": {
        "is_curved": 0, "curv_radius_mm": 72, "outer_diameter_mm": 60}}) is None


def test_fov_box_geometry():
    trans, focus = np.array([0.0, 0, 0]), np.array([0.0, 0, 50])
    box = pl.fov_box(trans, focus, diameter_mm=20, length_mm=120)
    assert box["corners"].shape == (8, 3) and len(box["edges"]) == 12
    # axial: starts at the transducer (z=0), ends length_mm along the beam axis
    np.testing.assert_allclose(sorted(box["corners"][:, 2].tolist()), [0.0] * 4 + [120.0] * 4,
                               atol=1e-9)
    # lateral: full width is diameter_mm regardless of where the focus sits along the axis
    assert np.ptp(box["corners"][:, 0]) == pytest.approx(20) and \
        np.ptp(box["corners"][:, 1]) == pytest.approx(20)
    for edge in box["edges"]:
        assert edge.shape == (2, 3)


def test_resolve_focal_distance_ep():
    bowl = pl.Bowl.from_transducer({"type": "annular",
                                    "annular": {"curv_radius_mm": 100, "dist_geom_ep_mm": 92.7,
                                               "elem_od_mm": [65]}})
    assert pl.resolve_focal_distance_ep({"focal_distance_ep": 85.0}, bowl) == 85.0
    # derived from focal_distance_bowl via bowl.offset when ep itself isn't set
    derived = pl.resolve_focal_distance_ep({"focal_distance_bowl": 92.3}, bowl)
    assert derived == pytest.approx(92.3 - bowl.offset)
    assert pl.resolve_focal_distance_ep({}, bowl) is None
    assert pl.resolve_focal_distance_ep({"focal_distance_bowl": 92.3}, None) is None


def test_snap_to_scalp(tmp_path):
    # 1 mm grid, a 'head' sphere of radius 50 mm around the origin
    shape = (121, 121, 121)
    affine = np.eye(4)
    affine[:3, 3] = -60
    g = np.indices(shape).reshape(3, -1).T - 60
    head = (np.linalg.norm(g, axis=1) <= 50).reshape(shape).astype(np.uint8)
    path = tmp_path / "seg.nii.gz"
    nib.save(nib.Nifti1Image(head, affine), str(path))
    seg = pl.Volume(path)
    bowl = pl.Bowl(63.2, 64, 52.38)
    apex = pl.snap_to_scalp(seg, [0, 0, 0], [0, 0, 1], bowl, gap_mm=5)
    # exit plane 5 mm outside the scalp (z = 50), apex a further offset behind
    assert apex[2] == pytest.approx(50 + 5 + bowl.offset, abs=0.6)


def test_read_brainsight(tmp_path):
    text = ("# Version: 13\n# Coordinate system: NIfTI:Aligned\n"
            "# Target name\tLoc. X\tLoc. Y\tLoc. Z\tm0n0\tm0n1\tm0n2\tm1n0\tm1n1\tm1n2\tm2n0\tm2n1\tm2n2\n"
            "left_amy\t-21.5\t-3.0\t-18.0\t1\t0\t0\t0\t1\t0\t0\t0\t1\n")
    path = tmp_path / "traj.txt"
    path.write_text(text)
    t = pl.read_trajectory(path)
    assert t.name == "left_amy" and t.source == "Brainsight"
    np.testing.assert_allclose(t.target, [-21.5, -3, -18])
    np.testing.assert_allclose(t.towards_transducer, [0, 0, 1])   # BabelBrain: transducer along +z


def test_brainsight_rejects_non_nifti(tmp_path):
    path = tmp_path / "traj.txt"
    path.write_text("# Coordinate system: Scanner\n# Target name\tLoc. X\n"
                    "t\t0\t0\t0\t1\t0\t0\t0\t1\t0\t0\t0\t1\n")
    with pytest.raises(ValueError, match="NIfTI"):
        pl.read_trajectory(path)


def test_slicer_matches_babelbrain_conversion(tmp_path):
    # Identity rotation about a fixed centre, translated in LPS
    text = ("#Insight Transform File V1.0\n#Transform 0\n"
            "Transform: AffineTransform_double_3_3\n"
            "Parameters: 1 0 0 0 1 0 0 0 1 10 20 30\nFixedParameters: 0 0 0\n")
    path = tmp_path / "traj.tfm"
    path.write_text(text)
    t = pl.read_trajectory(path)
    # itk_to_BSight: inv(ras2lps @ M @ ras2lps), rotation negated
    itk = np.eye(4)
    itk[:3, 3] = [10, 20, 30]
    ras2lps = np.diag([-1.0, -1, 1, 1])
    ref = np.linalg.inv(ras2lps @ itk @ ras2lps)
    ref[:3, :3] = -ref[:3, :3]
    np.testing.assert_allclose(t.target, ref[:3, 3])
    np.testing.assert_allclose(t.towards_transducer, ref[:3, 2])


def _kplan_file(path, apex_vox, axes, target_vox=None, vox=None, kind="kdp"):
    """A k-Plan file in the layout of k-Plan TPM v1.2.0: world [m] = 0-based voxel x voxel size."""
    h5py = pytest.importorskip("h5py")
    vox = np.asarray(vox, float)
    T = np.zeros((4, 4))
    T[:3, :3] = axes                                   # rows: transducer x, y, z axes (beam along -z)
    T[3, :3] = np.asarray(apex_vox, float) * vox / 1e3
    T[3, 3] = 1
    with h5py.File(str(path), "w") as h:
        if kind == "kps":
            h.attrs["file_type"] = np.bytes_("k-Plan Transducer Position")
            g = h.create_group("1")
            g.attrs["transform_label"] = np.bytes_("pose")
            g["position_transform"] = T.reshape(1, 4, 4).astype(np.float32)
        else:
            h.attrs["file_type"] = np.bytes_("k-Plan Device Parameters")
            h.attrs["plan_label"] = np.bytes_("plan")
            sp = h.create_group("sonications/1/sonication_parameters")
            sp["position_transform"] = T.reshape(1, 4, 4).astype(np.float32)
            sp["target_position"] = (np.asarray(target_vox, float) * vox / 1e3).reshape(1, 1, 3).astype(np.float32)
            sp["focal_distance"] = np.full((1, 1, 1), 0.05, np.float32)
            h.create_group("transducer").attrs["transducer_id"] = np.bytes_("TX-1")
            h["transducer/properties/focal_distance_surface_offset"] = np.full((1, 1, 1), 0.0108, np.float32)
    return path


def test_read_kplan_maps_world_to_t1_ras(vol, tmp_path):
    # Oblique, flipped, anisotropic T1 (AFFINE): k-Plan received the copy with the
    # rotation removed, so world = 0-based voxel x voxel size on the same voxels.
    vox = vol.voxel_size
    apex_v, target_v = np.array([2.0, 25, 30]), np.array([30.0, 26, 31])
    d = (apex_v - target_v) * vox
    z = d / np.linalg.norm(d)                         # k-Plan world direction target -> apex
    x = np.cross([0, 0, 1], z); x /= np.linalg.norm(x); y = np.cross(z, x)
    t = pl.read_trajectory(_kplan_file(tmp_path / "p.kdp", apex_v, np.vstack([x, y, z]), target_v, vox), t1=vol)
    assert t.source == "k-Plan" and t.name == "plan"
    np.testing.assert_allclose(t.transducer, vol.to_ras(apex_v)[0], atol=1e-3)
    np.testing.assert_allclose(t.target, vol.to_ras(target_v)[0], atol=1e-3)
    ref = vol.to_ras(apex_v)[0] - vol.to_ras(target_v)[0]
    np.testing.assert_allclose(t.towards_transducer, ref / np.linalg.norm(ref), atol=1e-5)
    assert any("rotated or flipped" in n for n in t.notes)
    assert any("50.00 mm" in n and "10.80 mm" in n for n in t.notes)


def test_read_kplan_positions_have_no_target(vol, tmp_path):
    t = pl.read_trajectory(_kplan_file(tmp_path / "p.kps", [2, 25, 30], np.eye(3), vox=vol.voxel_size,
                                       kind="kps"), t1=vol)
    assert t.target is None and t.name == "pose"
    np.testing.assert_allclose(t.transducer, vol.to_ras([2, 25, 30])[0], atol=1e-3)


def test_read_kplan_checks(vol, tmp_path):
    f = _kplan_file(tmp_path / "p.kdp", [2, 25, 30], np.eye(3), [300, 25, 30], vol.voxel_size)
    with pytest.raises(ValueError, match="outside this T1"):
        pl.read_trajectory(f, t1=vol)
    with pytest.raises(ValueError, match="T1"):
        pl.read_trajectory(f, t1=None)
    h5py = pytest.importorskip("h5py")
    other = tmp_path / "other.h5"
    with h5py.File(str(other), "w") as h:
        h["x"] = 1
    with pytest.raises(ValueError, match="not a k-Plan"):
        pl.read_trajectory(other, t1=vol)


def _fake_mni2subject_coords(bin_dir, offset=(1.0, 2.0, 3.0)):
    """A stand-in for SimNIBS's CLI tool: reads the input CSV, adds `offset`, writes it back
    in the same 'Generic, x, y, z' format -- enough to test the subprocess/CSV plumbing without
    a real SimNIBS install or m2m folder."""
    import stat
    exe = bin_dir / "mni2subject_coords"
    exe.write_text(
        "#!/usr/bin/env python3\n"
        "import csv, sys\n"
        "args = dict(zip(sys.argv[1::2], sys.argv[2::2]))\n"
        f"offset = {tuple(offset)!r}\n"
        "with open(args['-s']) as f: rows = list(csv.reader(f))\n"
        "with open(args['-o'], 'w', newline='') as f:\n"
        "    w = csv.writer(f)\n"
        "    for row in rows:\n"
        "        w.writerow([row[0]] + [float(v) + o for v, o in zip(row[1:4], offset)])\n")
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    return exe


def test_mni_to_subject(tmp_path):
    bin_dir, m2m = tmp_path / "bin", tmp_path / "m2m_sub-001"
    bin_dir.mkdir()
    m2m.mkdir()
    _fake_mni2subject_coords(bin_dir)
    out = pl.mni_to_subject([10.0, -20.0, 30.0], m2m, str(bin_dir))
    np.testing.assert_allclose(out, [11.0, -18.0, 33.0])


def test_mni_to_subject_missing_tool_or_m2m(tmp_path):
    bin_dir, m2m = tmp_path / "bin", tmp_path / "m2m_sub-001"
    bin_dir.mkdir()
    m2m.mkdir()
    with pytest.raises(pl.SimnibsCoordError, match="not found"):
        pl.mni_to_subject([0, 0, 0], m2m, str(bin_dir))          # tool missing
    _fake_mni2subject_coords(bin_dir)
    with pytest.raises(pl.SimnibsCoordError, match="not found"):
        pl.mni_to_subject([0, 0, 0], tmp_path / "no_m2m", str(bin_dir))


def test_beam_qc_layered_sphere(tmp_path):
    # 1 mm grid: brain r <= 40 (label 2), skull 40-46 (4), scalp 46-50 (5)
    shape = (121, 121, 121)
    affine = np.eye(4)
    affine[:3, 3] = -60
    g = np.indices(shape).reshape(3, -1).T - 60
    r = np.linalg.norm(g, axis=1)
    lab = np.zeros(len(r), np.uint8)
    lab[r <= 50] = 5
    lab[r <= 46] = 4
    lab[r <= 40] = 2
    path = tmp_path / "seg.nii.gz"
    nib.save(nib.Nifti1Image(lab.reshape(shape), affine), str(path))
    seg = pl.Volume(path)
    bowl = pl.Bowl(63.2, 64, 52.38)
    apex = pl.snap_to_scalp(seg, [0, 0, 0], [0, 0, 1], bowl, gap_mm=5)
    q = pl.beam_qc(seg, apex, [0, 0, 0], bowl)
    assert q["skull_path_mm"] == pytest.approx(6, abs=1)
    assert q["skull_incidence_deg"] < 8                    # radial beam: normal incidence
    assert q["beam_skull_path_mean_mm"] >= q["skull_path_mm"] - 0.5
    # rim (radius 32) above a sphere of radius 50: scalp at z = sqrt(50^2 - 32^2)
    exit_z = apex[2] - bowl.offset
    assert q["rim_gap_mean_mm"] == pytest.approx(exit_z - math.sqrt(50 ** 2 - 32 ** 2), abs=1)
    assert q["rim_gap_sd_mm"] < 1                          # symmetric: even contact (voxel staircase only)
    assert q["eye_in_beam"] == 0
    # transducer pushed into the head: the rim cuts into the scalp
    q = pl.beam_qc(seg, np.asarray(apex) - [0, 0, 20], [0, 0, 0], bowl)
    assert q["rim_gap_min_mm"] < 0


def test_beam_plane_basis():
    axis, lat = pl.beam_plane_basis([0, 0, 50], [0, 30, 0])       # oblique beam
    assert axis @ lat == pytest.approx(0, abs=1e-12)
    assert lat[2] > 0                                              # 0°: towards superior
    _, lat90 = pl.beam_plane_basis([0, 0, 50], [0, 30, 0], 90)
    assert lat @ lat90 == pytest.approx(0, abs=1e-12)
    _, lat_v = pl.beam_plane_basis([0, 0, 50], [0, 0, 0])         # vertical beam: anterior
    np.testing.assert_allclose(lat_v, [0, 1, 0], atol=1e-12)


def test_segmentation_to_medium():
    names = ["water", "brain", "skin", "skull", "skull_trabecular", "skull_cortical"]
    layers = {"water": [0, 3, 6, 9, 10], "brain": [1, 2], "skin": [5], "skull": [4],
              "skull_cortical": [7], "skull_trabecular": [8]}
    labels = np.array([0, 1, 2, 3, 4, 5, 7, 8, 11])
    med = pl.segmentation_to_medium(labels, layers, names)
    # multi-layer skull: bone (4) and compact bone (7) are cortical, spongy (8) trabecular;
    # an unmapped label (11) is water
    assert [names[i - 1] for i in med] == ["water", "brain", "brain", "water", "skull_cortical",
                                           "skin", "skull_cortical", "skull_trabecular", "water"]
    single = {"water": [0], "brain": [1, 2], "skin": [5], "skull": [4, 7, 8]}
    assert names[pl.segmentation_to_medium(np.array([7]), single, names)[0] - 1] == "skull"


def test_medium_clearance():
    bowl = pl.Bowl(63.2, 64, 52.38)                 # exit plane 10.82 mm from the apex
    depth = np.arange(-10, 100, 0.5)
    offsets = np.arange(-40, 40.25, 0.5)
    med = np.ones((len(offsets), len(depth)), int)  # all water (index 1)
    med[:, depth >= 20] = 3                          # flat skin 20 mm from the apex
    on_axis, worst = pl.medium_clearance(med, depth, offsets, bowl, water_index=1)
    assert on_axis == pytest.approx(20 - bowl.offset)
    assert worst == pytest.approx(on_axis)
    med[np.abs(offsets) > 25, :] = 3                 # tissue beside the axis, behind the exit plane
    _, worst = pl.medium_clearance(med, depth, offsets, bowl, water_index=1)
    assert worst < 0
    med[:] = 0                                       # outside the simulation grid: no tissue
    assert all(math.isnan(v) for v in pl.medium_clearance(med, depth, offsets, bowl, 1))


def test_voxels_are_cached_locally_until_the_image_changes(tmp_path, monkeypatch):
    path = tmp_path / "t1.nii.gz"
    data = np.arange(24, dtype=np.float32).reshape(2, 3, 4)
    nib.save(nib.Nifti1Image(data, AFFINE), str(path))
    np.testing.assert_array_equal(pl.Volume(path).data, data)
    cached = pl._voxel_cache_file(path)
    assert cached is not None and cached.is_file()

    reads = []
    real = pl._read_voxels_from
    monkeypatch.setattr(pl, "_read_voxels_from", lambda p, img: reads.append(p) or real(p, img))
    np.testing.assert_array_equal(pl.Volume(path).data, data)
    assert reads == []                                   # served from the cache

    nib.save(nib.Nifti1Image(data + 1, AFFINE), str(path))
    import os
    st = path.stat()
    os.utime(path, ns=(st.st_atime_ns, st.st_mtime_ns + 10**9))   # a newer modification time
    np.testing.assert_array_equal(pl.Volume(path).data, data + 1)
    assert reads == [path]                               # changed image: read again


def test_voxel_cache_keeps_the_most_recent_volumes(tmp_path, monkeypatch):
    monkeypatch.setattr(pl, "VOXEL_CACHE_KEEP", 2)
    for i in range(4):
        path = tmp_path / f"v{i}.nii.gz"
        nib.save(nib.Nifti1Image(np.full((2, 2, 2), i, np.float32), AFFINE), str(path))
        pl.Volume(path).data
    folder = pl._voxel_cache_file(path).parent
    assert len(list(folder.glob("*.npy"))) == 2
