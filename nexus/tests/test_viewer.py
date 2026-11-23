"""Placement viewer on a synthetic head, offscreen."""

import math
import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")
pytest.importorskip("matplotlib")
nib = pytest.importorskip("nibabel")

from PyQt6.QtWidgets import QApplication  # noqa: E402

from prestus_nexus import placement as pl  # noqa: E402
from prestus_nexus.equipment import EquipmentLibrary  # noqa: E402
from prestus_nexus.model import ConfigModel  # noqa: E402
from prestus_nexus.ui.placement_viewer import PlacementPage  # noqa: E402


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
    assert "Exit plane → scalp" in page.info.text()
    assert "Rim → scalp" in page.info.text()           # beam QC (placement.beam_qc)


def test_qc_issues_flag_depth_mismatch_and_exit_plane_in_scalp(app, schema, prestus_root, subject):
    """The viewer's red flags reach the Issues panel (MainWindow adds page.qc_issues())."""
    m, page = _page(app, schema, prestus_root, subject, "ras")
    page._write("focus_pos", np.array([0.0, 0.0, 10.0]))
    page.direction = np.array([0.0, 0.0, 1.0])
    page._snap()
    page._match_focal_distance()
    assert page.qc_issues() == []                  # on target, 5 mm outside the scalp
    m.set("transducer.focal_distance_ep", m.get("transducer.focal_distance_ep") + 10)
    msgs = [i.message for i in page.qc_issues()]
    assert len(msgs) == 1 and "short of the focus" in msgs[0]
    page.gap.setValue(-8)                          # exit plane 8 mm inside the scalp
    page._snap()
    page._match_focal_distance()
    msgs = [i.message for i in page.qc_issues()]
    assert len(msgs) == 1 and "inside the scalp" in msgs[0]
    m.set("placement.mode", "localite")            # resolved at run time: not checked here
    assert page.qc_issues() == []


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


@pytest.mark.posix
def test_mni_entry_uses_nonl_for_focus_and_transducer(app, schema, prestus_root,
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
    assert page._transform_for("trans_pos") == "nonl"          # SimNIBS 4: no affine
    page._set_coords("focus_pos")
    np.testing.assert_allclose(m.get("transducer.focus_pos"), [5.0, 5.0, 5.0])
    page._set_coords("trans_pos")
    np.testing.assert_allclose(m.get("transducer.trans_pos"), [5.0, 5.0, 5.0])


@pytest.mark.posix
def test_mni_entry_reports_a_failed_nonlinear_warp(app, schema, prestus_root, subject,
                                                   tmp_path, monkeypatch):
    """No affine fallback (SimNIBS 4 writes none): a failed warp is reported, nothing is set."""
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
    page._set_coords("focus_pos")
    assert not m.get("transducer.focus_pos")                # nothing written
    assert len(warnings) == 1 and "affine" not in warnings[0]


def test_mni_entry_without_simnibs_having_run_names_charm(
        app, schema, prestus_root, tmp_path, monkeypatch):
    """The nonlinear warp comes from charm's output: without an m2m folder the message must say
    so, distinctly from just 'startup.simnibs_bin_path is unset'."""
    from PyQt6.QtWidgets import QMessageBox
    infos = []
    monkeypatch.setattr(QMessageBox, "information", lambda *a, **k: infos.append(a[2]))
    anat = tmp_path / "anat"                                # no m2m_sub-XXX under it
    anat.mkdir()
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    _fake_mni2subject_coords(bin_dir)
    m, page = _page(app, schema, prestus_root, anat, "ras")
    m.set("startup.simnibs_bin_path", str(bin_dir))
    page.coord_space.setCurrentIndex(2)                     # MNI mm
    page._set_coords("focus_pos")
    assert not m.get("transducer.focus_pos")                  # nothing written
    assert infos and "charm" in infos[-1] and "m2m" in infos[-1]


def test_mni_entry_distinguishes_missing_bin_path_from_missing_m2m(app, schema, prestus_root,
                                                                    subject, monkeypatch):
    from PyQt6.QtWidgets import QMessageBox
    infos = []
    monkeypatch.setattr(QMessageBox, "information", lambda *a, **k: infos.append(a[2]))
    m, page = _page(app, schema, prestus_root, subject, "ras")      # m2m exists, bin path unset
    page.coord_space.setCurrentIndex(2)                     # MNI mm
    page._goto_coords()
    assert infos and "simnibs_bin_path" in infos[-1] and "charm" not in infos[-1]


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


class _FakeBridge:
    """Answers nexus_neuronav 'file_series' with a recorded MATLAB answer."""
    def __init__(self, out):
        self.out, self.calls = out, []

    def call(self, fn, args, on_done=None, on_error=None, **_):
        self.calls.append((fn, args))
        on_done(self.out)


def _localite_answer(xml):
    pose = {"is_pose": True, "tracked": True, "error": ""}
    return {"file": str(xml), "markertype": "TriggerMarkers", "note": "", "series": [
        {"index": 1, "description": "", "n_markers": 1, "pos_dev_mm": 0, "is_pose": False,
         "tracked": True, "trans_pos": [], "focus_pos": [], "error": ""},       # planned point
        dict(pose, index=2, description="L_M1", n_markers=80, pos_dev_mm=1.2,
             trans_pos=[61, 61, 111], focus_pos=[61, 61, 71]),
        dict(pose, index=3, description="", n_markers=80, pos_dev_mm=0.8,
             trans_pos=[61, 21, 61], focus_pos=[61, 61, 61])]}


@pytest.mark.parametrize("pick, expect_target, expect_position", [
    (0, "L_M1", None),          # unique description: chosen by name
    (1, None, 3)])              # no description: chosen by series number
def test_localite_import_sets_file_mode(app, schema, prestus_root, subject, tmp_path, monkeypatch,
                                        pick, expect_target, expect_position):
    """A Localite XML in Import trajectory: localite mode with placement.localite.file and
    the chosen series; the viewer shows the positions MATLAB resolved for that series."""
    from PyQt6.QtWidgets import QFileDialog, QInputDialog, QMessageBox
    xml = tmp_path / "TriggerMarkers_Coil0_1.xml"
    xml.write_text("<TriggerMarkerList/>")
    shown = {}
    monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **k: (str(xml), ""))
    monkeypatch.setattr(QInputDialog, "getItem",
                        lambda _p, _t, _l, items, *a, **k: (shown.setdefault("items", items)[pick], True))
    monkeypatch.setattr(QMessageBox, "information", lambda *a, **k: shown.setdefault("info", a[2]))
    monkeypatch.setattr(QMessageBox, "warning", lambda *a, **k: pytest.fail(str(a[2:])))
    m, page = _page(app, schema, prestus_root, subject, "voxel")
    m.set("transducer.trans_pos", [1, 2, 3])
    bridge = _FakeBridge(_localite_answer(xml))
    page.bridge_getter = lambda: bridge
    page.saved_config_getter = lambda: str(tmp_path / "study.yaml")
    page._import_trajectory()                        # manual mode not needed for a Localite file
    assert bridge.calls[0][1]["action"] == "file_series" and bridge.calls[0][1]["file"] == str(xml)
    assert len(shown["items"]) == 2 and "L_M1" in shown["items"][0]   # planned point not offered
    assert m.get("placement.mode") == "localite"
    assert m.get("placement.localite.file") == str(xml)
    assert (m.get("placement.localite.target") or None) == expect_target
    if expect_position is not None:
        assert m.get("placement.localite.position") == expect_position
    assert m.get("transducer.trans_pos") in (None, [], "")       # resolved from the file
    chosen = _localite_answer(xml)["series"][1 + pick]
    trans, focus = page.positions()
    np.testing.assert_allclose(trans, pl.prestus_voxel_to_ras(page.t1, chosen["trans_pos"]))
    np.testing.assert_allclose(focus, pl.prestus_voxel_to_ras(page.t1, chosen["focus_pos"]))
    assert "Localite" in shown["info"]
    # another series configured by hand: the imported positions no longer apply
    m.set("placement.localite.target", "")
    m.set("placement.localite.position", 1 if pick else 2)
    page._forget()
    assert page.positions() == (None, None)


def test_placement_undo_redo(app, schema, prestus_root, tmp_path):
    """The viewer's buttons step through the window's one edit history (Edit ▸ Undo)."""
    import numpy as np
    from prestus_nexus.model import ConfigModel
    from prestus_nexus.session import Session
    from prestus_nexus.ui.main_window import VIEWER, MainWindow
    m = ConfigModel(schema)
    m.set("placement.manual.coordinate_space", "ras")
    w = MainWindow(m, Session(tmp_path / "s.yaml"), prestus_root)
    v = w.pages[VIEWER]
    assert v.history is w.history
    assert not v.undo_btn.isEnabled() and not v.redo_btn.isEnabled()
    v._write("focus_pos", np.array([1.0, 2.0, 3.0]))
    app.processEvents()                                   # one undo step per action
    v._write("focus_pos", np.array([4.0, 5.0, 6.0]))
    v._write("trans_pos", np.array([0.0, 0.0, -60.0]))   # same action: no extra step
    app.processEvents()
    assert len(w.history._undo) == 2 and v.undo_btn.isEnabled()
    assert w.undo_action.text() == "Undo move focus"
    m.set("transducer.freq_hz", 250e3)                   # an edit on another tab: the same history
    v.undo()
    assert m.get("transducer.freq_hz") != 250e3 and m.get("transducer.focus_pos") == [4.0, 5.0, 6.0]
    v.undo()
    assert m.get("transducer.focus_pos") == [1.0, 2.0, 3.0]
    assert not m.is_overridden("transducer.trans_pos")
    w.undo_action.trigger()                               # the menu and the buttons are one history
    assert not m.is_overridden("transducer.focus_pos") and not v.undo_btn.isEnabled()
    v.redo()
    v.redo()
    assert m.get("transducer.focus_pos") == [4.0, 5.0, 6.0]
    assert m.get("transducer.trans_pos") == [0.0, 0.0, -60.0] and v.redo_btn.isEnabled()
    v._add_transducer()                                   # adding a transducer undoes too
    app.processEvents()
    assert m.n_transducers == 2 and not v.redo_btn.isEnabled()
    v.undo()
    assert m.n_transducers == 1 and v.ti == 0 and m.get("transducer.focus_pos") == [4.0, 5.0, 6.0]
    w.dirty = False
    w.close()


def test_standalone_viewer_keeps_its_own_history(app, schema):
    import numpy as np
    from prestus_nexus.model import ConfigModel
    m = ConfigModel(schema)
    m.set("placement.manual.coordinate_space", "ras")
    v = PlacementPage(m)
    v._write("focus_pos", np.array([1.0, 2.0, 3.0]))
    app.processEvents()
    assert v.undo_btn.isEnabled()
    v.undo()
    assert not m.is_overridden("transducer.focus_pos") and v.redo_btn.isEnabled()


def test_loading_card_goes_when_a_waiting_caller_applied_the_read(app, schema, prestus_root, subject):
    """A caller that waits for the read (the Issues check asks for the target depth) applies
    it before the background relay arrives; the page must still swap 'Loading' for the slices."""
    import time
    m = ConfigModel(schema, EquipmentLibrary(prestus_root))
    for k, v in {"subject_id": 4, "path.anat": str(subject), "transducer.serial": "CTX_500_026",
                 "placement.manual.coordinate_space": "ras",
                 "transducer.trans_pos": [-30.0, 0.0, 40.0], "transducer.focus_pos": [0.0, 0.0, 0.0]}.items():
        m.set(k, v)
    page = PlacementPage(m)
    page.refresh()                                   # shows 'Loading', read runs in the background
    assert page.images_loading() and page.stack.is_empty()
    page.target_depth_ep(0)                          # waits for and applies the read, no redraw
    assert not page.images_loading()
    deadline = time.time() + 20
    while page.stack.is_empty() and time.time() < deadline:
        app.processEvents()                          # the relay's late signal redraws
        time.sleep(0.01)
    assert not page.stack.is_empty()


def test_images_load_in_the_background_and_grid_z_shows_with_two_transducers(
        app, schema, prestus_root, subject):
    import time
    m = ConfigModel(schema, EquipmentLibrary(prestus_root))
    for k, v in {"subject_id": 4, "path.anat": str(subject), "transducer.serial": "CTX_500_026",
                 "placement.manual.coordinate_space": "ras",
                 "transducer.trans_pos": [-30.0, 0.0, 40.0], "transducer.focus_pos": [0.0, 0.0, 0.0]}.items():
        m.set(k, v)
    page = PlacementPage(m)
    page.refresh()                                   # does not wait for the read
    deadline = time.time() + 20
    while page.images_loading() and time.time() < deadline:
        app.processEvents()                          # the relay draws when the read is done
        time.sleep(0.01)
    assert not page.images_loading() and page._t1 is not None
    assert not page.stack.is_empty()
    assert not page.grid_z_action.isChecked()        # one transducer: off by default
    m.add_transducer()
    m.set("transducer.trans_pos", [30.0, 0.0, 40.0], 1)
    m.set("transducer.focus_pos", [0.0, 0.0, 0.0], 1)
    page.refresh()
    assert page.grid_z_action.isChecked()
    line = page._grid_z_line()
    d = line[1] - line[0]
    np.testing.assert_allclose(np.abs(d / np.linalg.norm(d)), [0, 0, 1], atol=1e-6)   # bisector
    page.grid_z_action.setChecked(False)             # the user's choice sticks
    page.refresh()
    assert not page.grid_z_action.isChecked()
