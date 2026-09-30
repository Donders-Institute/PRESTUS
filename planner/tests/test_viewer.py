"""Placement viewer on a synthetic head, offscreen."""

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")
pytest.importorskip("matplotlib")
nib = pytest.importorskip("nibabel")

from PyQt6.QtWidgets import QApplication  # noqa: E402

from prestus_planner import placement as pl  # noqa: E402
from prestus_planner.equipment import EquipmentLibrary  # noqa: E402
from prestus_planner.model import ConfigModel  # noqa: E402
from prestus_planner.ui.placement_viewer import PlacementPage  # noqa: E402


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def subject(tmp_path):
    """1 mm T1 with a 50 mm 'head' sphere and matching segmentation."""
    shape = (121, 121, 121)
    affine = np.eye(4)
    affine[:3, 3] = -60
    g = np.indices(shape).reshape(3, -1).T - 60
    r = np.linalg.norm(g, axis=1).reshape(shape)
    anat = tmp_path / "anat"
    (anat / "m2m_sub-004").mkdir(parents=True)
    nib.save(nib.Nifti1Image((r <= 50).astype(np.float32) * 100, affine), str(anat / "sub-004_T1w.nii.gz"))
    seg = np.where(r <= 40, 2, np.where(r <= 45, 4, np.where(r <= 50, 5, 0))).astype(np.uint8)
    nib.save(nib.Nifti1Image(seg, affine), str(anat / "m2m_sub-004" / "final_tissues.nii.gz"))
    return anat


def _page(app, schema, prestus_root, anat, space):
    m = ConfigModel(schema, EquipmentLibrary(prestus_root))
    for k, v in {"subject_id": 4, "path.anat": str(anat), "transducer.serial": "CTX_500_026",
                 "transducer.combo.ds_serial": "SC_203_035", "transducer.focal_distance_ep": 60,
                 "placement.manual.coordinate_space": space}.items():
        m.set(k, v)
    page = PlacementPage(m)
    page.refresh()
    return m, page


@pytest.mark.parametrize("space", ["ras", "voxel"])
def test_place_snap_and_match(app, schema, prestus_root, subject, space):
    m, page = _page(app, schema, prestus_root, subject, space)
    assert page.t1 is not None and page.seg is not None
    page._write("focus_pos", np.array([0.0, 0.0, 10.0]))
    page.direction = np.array([0.0, 0.0, 1.0])
    page._snap()
    trans, focus = page.positions()
    bowl = page._bowl()
    # exit plane 5 mm (default gap) outside the scalp at z = 50
    assert trans[2] - bowl.offset == pytest.approx(55, abs=1.1)
    assert page._tissue_at(focus) == "grey matter"
    if space == "voxel":
        assert m.get("transducer.focus_pos") == pl.ras_to_prestus_voxel(page.t1, [0, 0, 10])
    page._match_focal_distance()
    ep = m.get("transducer.focal_distance_ep")
    assert ep == pytest.approx(np.linalg.norm(trans - focus) - bowl.offset, abs=0.05)
    page._draw()                                   # renders without error
    assert "Exit plane vs scalp" in page.info.text()


def _fake_mni2subject_coords(bin_dir, fail_for: str | None = None, offset=5.0):
    """A stand-in SimNIBS CLI: adds `offset` to every coordinate, or exits non-zero when called
    with `-t {fail_for}` (to simulate a missing/failed nonlinear warp)."""
    import stat
    exe = bin_dir / "mni2subject_coords"
    exe.write_text(
        "#!/usr/bin/env python3\n"
        "import csv, sys\n"
        "args = dict(zip(sys.argv[1::2], sys.argv[2::2]))\n"
        f"if args['-t'] == {fail_for!r}:\n"
        "    print('nonlinear transform not found', file=sys.stderr); sys.exit(1)\n"
        "with open(args['-s']) as f: rows = list(csv.reader(f))\n"
        "with open(args['-o'], 'w', newline='') as f:\n"
        "    w = csv.writer(f)\n"
        "    for row in rows:\n"
        f"        w.writerow([row[0]] + [float(v) + {offset} for v in row[1:4]])\n")
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    return exe


def test_coordinate_entry_ras_and_voxel(app, schema, prestus_root, subject):
    m, page = _page(app, schema, prestus_root, subject, "ras")
    page.coord_space.setCurrentIndex(0)                     # RAS mm
    page.cx.setValue(1.0)
    page.cy.setValue(2.0)
    page.cz.setValue(3.0)
    page._set_coords("focus_pos")
    assert m.get("transducer.focus_pos") == [1.0, 2.0, 3.0]

    page.coord_space.setCurrentIndex(1)                     # Voxel
    page.cx.setValue(61)
    page.cy.setValue(61)
    page.cz.setValue(61)
    page._goto_coords()
    np.testing.assert_allclose(page.cursor, pl.prestus_voxel_to_ras(page.t1, [61, 61, 61]))


def test_mni_entry_uses_nonl_for_focus_and_affine_for_transducer(app, schema, prestus_root,
                                                                 subject, tmp_path, monkeypatch):
    from PyQt6.QtWidgets import QMessageBox
    monkeypatch.setattr(QMessageBox, "warning", lambda *a, **k: pytest.fail(str(a[2:])))
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    _fake_mni2subject_coords(bin_dir)                       # both transforms succeed, +5 offset
    m, page = _page(app, schema, prestus_root, subject, "ras")
    m.set("startup.simnibs_bin_path", str(bin_dir))
    page.coord_space.setCurrentIndex(2)                     # MNI mm
    page.cx.setValue(0.0)
    page.cy.setValue(0.0)
    page.cz.setValue(0.0)

    assert page._transform_for("focus_pos") == "nonl"
    assert page._transform_for("trans_pos") == "12dof"
    page._set_coords("focus_pos")
    np.testing.assert_allclose(m.get("transducer.focus_pos"), [5.0, 5.0, 5.0])
    page._set_coords("trans_pos")
    np.testing.assert_allclose(m.get("transducer.trans_pos"), [5.0, 5.0, 5.0])


def test_mni_entry_falls_back_to_affine_when_nonlinear_fails(app, schema, prestus_root, subject,
                                                              tmp_path, monkeypatch):
    from PyQt6.QtWidgets import QMessageBox
    warnings = []
    monkeypatch.setattr(QMessageBox, "warning", lambda *a, **k: warnings.append(a[2]))
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    _fake_mni2subject_coords(bin_dir, fail_for="nonl")      # nonlinear warp "missing"
    m, page = _page(app, schema, prestus_root, subject, "ras")
    m.set("startup.simnibs_bin_path", str(bin_dir))
    page.coord_space.setCurrentIndex(2)                     # MNI mm
    page.cx.setValue(0.0)
    page.cy.setValue(0.0)
    page.cz.setValue(0.0)
    page._set_coords("focus_pos")                           # nonl fails -> falls back to 12dof
    np.testing.assert_allclose(m.get("transducer.focus_pos"), [5.0, 5.0, 5.0])
    assert warnings and "affine" in warnings[0]


def test_trajectory_import(app, schema, prestus_root, subject, tmp_path, monkeypatch):
    from PyQt6.QtWidgets import QFileDialog, QMessageBox
    traj = tmp_path / "t.txt"
    traj.write_text("# Coordinate system: NIfTI:Aligned\n"
                    "tgt\t0\t0\t10\t1\t0\t0\t0\t1\t0\t0\t0\t1\n")
    monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **k: (str(traj), ""))
    monkeypatch.setattr(QMessageBox, "information", lambda *a, **k: None)
    monkeypatch.setattr(QMessageBox, "warning", lambda *a, **k: pytest.fail(str(a[2:])))
    m, page = _page(app, schema, prestus_root, subject, "ras")
    page._import_trajectory()
    trans, focus = page.positions()
    np.testing.assert_allclose(focus, [0, 0, 10])
    # transducer on +z at focal_distance_bowl = ep + offset
    assert trans[2] - focus[2] == pytest.approx(60 + page._bowl().offset, abs=0.1)
