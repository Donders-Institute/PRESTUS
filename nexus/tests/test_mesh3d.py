"""Meshes of the 3D placement view: head surface, bowl rings, camera rotation."""

import math

import numpy as np
import pytest

from prestus_nexus import mesh3d as m3
from prestus_nexus.placement import Bowl


def _closed(faces: np.ndarray) -> bool:
    edges = np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    return bool(np.all(counts == 2))


def test_voxel_surface_of_a_ball_is_closed_outward_and_in_world_mm():
    g = np.arange(20) - 9.5
    ball = (g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2) < 7 ** 2
    affine = np.diag([2.0, 2.0, 2.0, 1.0])
    affine[:3, 3] = [-19.0, -19.0, -19.0]                  # voxel 9.5 -> 0 mm
    v, f = m3.voxel_surface(ball, affine, smooth_iter=8)
    assert _closed(f)
    n = m3.vertex_normals(v, f)
    assert np.mean(np.sum(v * n, axis=1) > 0) > 0.99      # normals point outwards
    r = np.linalg.norm(v, axis=1)
    assert 11 < np.median(r) < 16                          # radius 7 voxels x 2 mm, smoothed


def test_mirroring_affine_keeps_normals_outward():
    ball = np.zeros((9, 9, 9), bool)
    ball[2:7, 2:7, 2:7] = True
    affine = np.diag([-1.0, 1.0, 1.0, 1.0])
    affine[:3, 3] = [4.0, -4.0, -4.0]
    v, f = m3.voxel_surface(ball, affine, smooth_iter=0)
    n = m3.vertex_normals(v, f)
    assert np.mean(np.sum(v * n, axis=1) > 0) > 0.99


def test_rings_lie_on_the_bowl_between_their_diameters():
    bowl = Bowl(curv_radius=100.0, aperture=65.0, dist_geom_ep=92.7)
    apex, axis = np.array([80.0, 0.0, 0.0]), np.array([-1.0, 0.0, 0.0])
    rings = m3.ring_meshes(bowl, [10.0, 22.1], [21.1, 28.8], apex, axis)
    assert len(rings) == 2
    centre = apex + axis * bowl.curv_radius                # every point is R from the centre
    for (v, _), (d_in, d_out) in zip(rings, [(10.0, 21.1), (22.1, 28.8)]):
        np.testing.assert_allclose(np.linalg.norm(v - centre, axis=1), 100.0, atol=1e-9)
        lateral = np.linalg.norm(np.cross(v - apex, axis), axis=1)
        assert lateral.min() == pytest.approx(d_in / 2) and lateral.max() == pytest.approx(d_out / 2)


def test_flat_rings_and_single_aperture_default():
    flat = Bowl(curv_radius=math.inf, aperture=33.6, dist_geom_ep=math.inf)
    (v, _), = m3.ring_meshes(flat, None, None, [0, 0, 0], [0, 0, 1])
    assert np.allclose(v[:, 2], 0) and np.linalg.norm(v[:, :2], axis=1).max() == pytest.approx(16.8)


def test_exit_plane_sits_at_the_offset():
    bowl = Bowl(curv_radius=62.94, aperture=64.0, dist_geom_ep=52.38)
    v, _ = m3.exit_plane_disc(bowl, [0, 0, 0], [0, 1, 0])
    np.testing.assert_allclose(v[:, 1], bowl.offset)


def test_look_rotation_turns_plus_z_onto_the_direction():
    pytest.importorskip("PyQt6.QtGui")
    from prestus_nexus.ui.view3d import look_rotation, to_scene
    for d in ([1, 0, 0], [0, 0, -1], [0.3, 0.9, 0.2], [0, -1, 0.01], [0, 0, 1]):
        d = np.asarray(d, float) / np.linalg.norm(d)
        w, x, y, z = look_rotation(d)
        q = np.array([x, y, z])
        k = np.array([0.0, 0.0, 1.0])
        rotated = k + 2 * w * np.cross(q, k) + 2 * np.cross(q, np.cross(q, k))
        np.testing.assert_allclose(rotated, d, atol=1e-9)
    s = to_scene([[1.0, 2.0, 3.0]])[0]                   # R, A, S -> x = R, y = S, z = -A
    np.testing.assert_allclose(s, [1.0, 3.0, -2.0])


def test_housing_cylinder_holds_the_elements():
    bowl = Bowl(curv_radius=62.94, aperture=60.83, dist_geom_ep=52.38)
    # a size on record (BabelBrain's BSonix case): front face in the exit plane, back past the apex
    d, depth, front, known = m3.housing_size(bowl.aperture, bowl.offset, {"diameter_mm": 70, "depth_mm": 34})
    assert (d, depth, front, known) == (70.0, 34.0, pytest.approx(bowl.offset), True)
    d, depth, _front, _known = m3.housing_size(bowl.aperture, bowl.offset, {"diameter_mm": 70, "depth_mm": 5})
    assert depth == pytest.approx(bowl.offset + m3.HOUSING_BACK_MARGIN_MM)   # too shallow: reaches past the apex
    # none: a generic cylinder of the aperture, exit plane to 15 mm behind the apex
    d, depth, front, known = m3.housing_size(bowl.aperture, bowl.offset, {})
    assert (d, depth, front, known) == (pytest.approx(60.83), pytest.approx(bowl.offset + 15.0),
                                        pytest.approx(bowl.offset), False)
    v, f = m3.housing_puck(d, depth, front, [0, 0, 0], [0, 0, 1])
    assert _closed(f)
    assert v[:, 2].max() == pytest.approx(bowl.offset) and v[:, 2].min() == pytest.approx(-15.0)
    assert np.hypot(v[:, 0], v[:, 1]).max() == pytest.approx(60.83 / 2)
    # the element rings (bowl from apex to exit plane) lie inside it
    rings = np.vstack([rv for rv, _rf in m3.ring_meshes(bowl, None, None, [0, 0, 0], [0, 0, 1])])
    assert rings[:, 2].min() >= v[:, 2].min() and rings[:, 2].max() <= v[:, 2].max() + 1e-6
    assert np.hypot(rings[:, 0], rings[:, 1]).max() <= 60.83 / 2 + 1e-6
    n = m3.vertex_normals(v, f)
    centred = v - np.array([0, 0, (v[:, 2].max() + v[:, 2].min()) / 2])
    assert np.mean(np.sum(centred * n, axis=1) > 0) > 0.99         # normals point outwards


def test_face_region_and_scrambling():
    from prestus_nexus.placement import EYE_LABEL
    labels = np.zeros((40, 40, 40), np.uint8)
    labels[18:22, 30:33, 18:22] = EYE_LABEL                    # eyes centred at A = 31, S = 19.5 (1 mm)
    a_from, s_to = m3.face_region(labels, np.eye(4))
    assert a_from == pytest.approx(31 - m3.FACE_MARGIN_MM) and s_to == pytest.approx(19.5 + m3.FACE_MARGIN_MM)
    assert m3.face_region(np.zeros((5, 5, 5), np.uint8), np.eye(4)) is None
    # a head-sized box surface with a 'nose' bump on the face
    head = np.zeros((40, 40, 40), bool)
    head[5:35, 5:35, 5:35] = True
    head[18:22, 35:39, 12:18] = True                           # nose: in front, below the eyes
    v, f = m3.voxel_surface(head, np.eye(4))
    sv, sf = m3.scramble_region(v, f, a_from, s_to)
    assert np.array_equal(sf, f)                               # no triangles removed: still closed
    face = (v[:, 1] > a_from) & (v[:, 2] < s_to)
    assert np.allclose(sv[~face], v[~face])                    # outside the face: unchanged
    nose = v[:, 1] > 35.5
    assert sv[nose, 1].mean() < v[nose, 1].mean() - 1.5       # the nose is flattened
    again, _ = m3.scramble_region(v, f, a_from, s_to)
    assert np.allclose(again, sv)                              # seeded: the same head looks the same


def test_equipment_housing_reaches_units_through_their_model(prestus_root):
    from prestus_nexus.equipment import EquipmentLibrary
    lib = EquipmentLibrary(prestus_root)
    model = lib.get("CTX_500")
    assert model.housing == {"diameter_mm": 70, "depth_mm": 20}
    assert lib.housing("CTX_500_026") == model.housing         # a unit gets its model's
    assert lib.housing("H246") == {}                            # none: the generic cylinder


def test_3d_page_toggles_housing_and_face(schema, prestus_root, tmp_path):
    """Housing and Face are checkboxes like Scalp and Brain; the face is off by default."""
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PyQt6.QtQuick3D")
    nib = pytest.importorskip("nibabel")
    from PyQt6.QtWidgets import QApplication

    from prestus_nexus.equipment import EquipmentLibrary
    from prestus_nexus.model import ConfigModel
    from prestus_nexus.placement import EYE_LABEL
    from prestus_nexus.ui.placement_viewer import PlacementPage
    from prestus_nexus.ui.view3d import HOUSING, SCALP, Placement3DPage
    app = QApplication.instance() or QApplication([])  # noqa: F841 (kept alive for the widgets)
    g = np.indices((81, 81, 81)).reshape(3, -1).T - 40
    r = np.linalg.norm(g, axis=1).reshape(81, 81, 81)
    affine = np.eye(4)
    affine[:3, 3] = -40
    seg = np.where(r <= 25, 2, np.where(r <= 35, 5, 0)).astype(np.uint8)
    seg[34:38, 70:73, 34:38] = EYE_LABEL                               # eyes in front (A ~ 31)
    anat = tmp_path / "anat"
    (anat / "m2m_sub-004").mkdir(parents=True)
    nib.save(nib.Nifti1Image((r <= 35).astype(np.float32), affine), str(anat / "sub-004_T1w.nii.gz"))
    nib.save(nib.Nifti1Image(seg, affine), str(anat / "m2m_sub-004" / "final_tissues.nii.gz"))
    m = ConfigModel(schema, EquipmentLibrary(prestus_root))
    for k, v in {"subject_id": 4, "path.anat": str(anat), "transducer.serial": "CTX_500_026",
                 "transducer.focal_distance_ep": 40, "placement.manual.coordinate_space": "ras",
                 "transducer.trans_pos": [0.0, 0.0, 50.0], "transducer.focus_pos": [0.0, 0.0, 0.0]}.items():
        m.set(k, v)
    viewer = PlacementPage(m)
    viewer.refresh()
    page = Placement3DPage(viewer)
    assert not page.show_face.isChecked() and page.show_housing.isChecked()
    colours = lambda: [c for _v, _f, c, _a in page.meshes()]          # noqa: E731
    scalp = lambda: next(f for _v, f, c, _a in page.meshes() if c == SCALP[0])   # noqa: E731
    assert HOUSING in colours()
    assert page._housing_known                                         # CTX_500's size, via the unit
    opaque = next(a for _v, _f, c, a in page.meshes() if c == HOUSING)
    assert opaque == 1.0                                               # opaque by default
    page.housing_opacity.setValue(40)
    assert next(a for _v, _f, c, a in page.meshes() if c == HOUSING) == pytest.approx(0.4)
    scrambled = next(v for v, _f, c, _a in page.meshes() if c == SCALP[0])
    page.show_face.setChecked(True)
    assert not np.allclose(next(v for v, _f, c, _a in page.meshes() if c == SCALP[0]), scrambled)  # the face is back
    page.show_housing.setChecked(False)
    assert HOUSING not in colours()


def test_mask_surface_crops_to_the_mask_and_keeps_world_mm():
    vol = np.zeros((60, 60, 60), np.float32)
    g = np.indices(vol.shape).reshape(3, -1).T
    vol.reshape(-1)[np.linalg.norm(g - [40, 20, 30], axis=1) <= 3] = 1.0
    affine = np.diag([2.0, 2.0, 2.0, 1.0])
    affine[:3, 3] = [-60.0, -60.0, -60.0]
    v, f = m3.mask_surface(vol, affine)
    w, wf = m3.voxel_surface(vol > 0.5, affine, smooth_iter=4)
    assert _closed(f) and len(f) == len(wf)
    assert np.allclose(np.sort(v, axis=0), np.sort(w, axis=0))           # same surface as uncropped
    assert np.allclose(v.mean(axis=0), affine[:3, :3] @ [40, 20, 30] + affine[:3, 3], atol=0.5)
    e, ef = m3.mask_surface(np.zeros((5, 5, 5)), np.eye(4))
    assert e.shape == (0, 3) and ef.shape == (0, 3)


def test_target_mask_path_follows_the_mode(schema, tmp_path):
    from prestus_nexus.model import ConfigModel
    from prestus_nexus.ui.view3d import target_mask_path
    m = ConfigModel(schema)
    m.set("path.sim", str(tmp_path))
    m.set("subject_id", 4)
    assert target_mask_path(m) is None                                   # manual: no mask
    m.set("placement.mode", "plantus")
    m.set("placement.plantus.target_name", "M1")
    assert target_mask_path(m) == tmp_path / "placement/sub-004/plantus/M1/sub-004_M1_mask.nii.gz"
    m.set("placement.mode", "tus_entry")
    assert target_mask_path(m).name == "sub-004_target_target_mask.nii.gz"
    m.set("placement.tus_entry.target_mask", str(tmp_path / "roi.nii.gz"))
    assert target_mask_path(m) == tmp_path / "roi.nii.gz"


def test_3d_page_opacity_sliders_and_target_mask(schema, prestus_root, tmp_path):
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PyQt6.QtQuick3D")
    nib = pytest.importorskip("nibabel")
    from PyQt6.QtWidgets import QApplication

    from prestus_nexus.equipment import EquipmentLibrary
    from prestus_nexus.model import ConfigModel
    from prestus_nexus.ui.placement_viewer import PlacementPage
    from prestus_nexus.ui.view3d import BRAIN, MASK, SCALP, TARGET, Placement3DPage, target_mask_path
    app = QApplication.instance() or QApplication([])  # noqa: F841
    g = np.indices((81, 81, 81)).reshape(3, -1).T - 40
    r = np.linalg.norm(g, axis=1).reshape(81, 81, 81)
    affine = np.eye(4)
    affine[:3, 3] = -40
    seg = np.where(r <= 25, 2, np.where(r <= 35, 5, 0)).astype(np.uint8)
    anat = tmp_path / "anat"
    (anat / "m2m_sub-004").mkdir(parents=True)
    nib.save(nib.Nifti1Image((r <= 35).astype(np.float32), affine), str(anat / "sub-004_T1w.nii.gz"))
    nib.save(nib.Nifti1Image(seg, affine), str(anat / "m2m_sub-004" / "final_tissues.nii.gz"))
    m = ConfigModel(schema, EquipmentLibrary(prestus_root))
    for k, v in {"subject_id": 4, "path.anat": str(anat), "path.sim": str(tmp_path / "sim"),
                 "transducer.serial": "CTX_500_026", "transducer.focal_distance_ep": 40,
                 "placement.manual.coordinate_space": "ras",
                 "transducer.trans_pos": [0.0, 0.0, 50.0], "transducer.focus_pos": [0.0, 0.0, 0.0]}.items():
        m.set(k, v)
    viewer = PlacementPage(m)
    viewer.refresh()
    page = Placement3DPage(viewer)
    alpha = lambda colour: next(a for _v, _f, c, a in page.meshes() if c == colour)   # noqa: E731
    assert alpha(SCALP[0]) == pytest.approx(SCALP[1]) and alpha(BRAIN[0]) == pytest.approx(BRAIN[1])
    page.scalp_opacity.setValue(70)
    page.brain_opacity.setValue(100)
    assert alpha(SCALP[0]) == pytest.approx(0.7) and alpha(BRAIN[0]) == pytest.approx(1.0)
    page.show_brain.setChecked(False)
    assert not page.brain_opacity.isEnabled()
    assert sum(c == TARGET for _v, _f, c, _a in page.meshes()) == 4      # sphere + 3-axis crosshair
    # a PlanTUS run's target mask is drawn once found
    page._load_mask()
    assert page._mask is None and page.show_mask.isHidden()
    m.set("placement.mode", "plantus")
    m.set("placement.plantus.target_name", "M1")
    path = target_mask_path(m)
    path.parent.mkdir(parents=True)
    nib.save(nib.Nifti1Image((r <= 4).astype(np.float32), affine), str(path))
    page._load_mask()
    assert not page.show_mask.isHidden()
    assert MASK[0] in [c for _v, _f, c, _a in page.meshes()]
    page.show_mask.setChecked(False)
    assert MASK[0] not in [c for _v, _f, c, _a in page.meshes()]
