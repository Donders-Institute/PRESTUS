"""Start tab: inputs derived from one T1 file, and the page writing them into the model."""

import os

import pytest

from prestus_nexus import quickstart as qs
from prestus_nexus.paths import resolve_pattern


def _touch(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"")
    return path


def test_bids_name_gives_subject_pattern_t2_and_sibling_segmentation(tmp_path):
    anat = tmp_path / "sub-012" / "anat"
    t1 = _touch(anat / "sub-012_T1w.nii.gz")
    _touch(anat / "sub-012_T2w.nii.gz")
    _touch(anat / "sub-013_T1w.nii.gz")                    # another subject in the same folder
    (tmp_path / "sub-012" / "simnibs" / "m2m_sub-012").mkdir(parents=True)
    got = qs.from_t1(t1)
    assert got.subject_id == 12
    assert got.t1_pattern == "sub-%1$03d_T1w.nii.gz"
    assert resolve_pattern(got.anat, got.t1_pattern, 12) == [t1.resolve()]
    assert got.t2_pattern == "sub-%1$03d_T2w.nii.gz"
    assert got.m2m_found and got.seg == str((tmp_path / "sub-012" / "simnibs").resolve())


def test_plain_name_is_kept_literally_and_keeps_the_subject_id(tmp_path):
    t1 = _touch(tmp_path / "ernie" / "ernie_T1.nii.gz")
    (tmp_path / "other" / "m2m_sub-001").mkdir(parents=True)       # a neighbour's segmentation
    got = qs.from_t1(t1, current_subject=1)
    assert got.subject_id is None and "subject_id" not in got.as_settings()
    assert got.t1_pattern == "ernie_T1.nii.gz"
    assert got.t2_pattern == "" and not got.m2m_found and got.seg == ""


def test_percent_in_name_is_escaped_for_sprintf(tmp_path):
    t1 = _touch(tmp_path / "50%_T1.nii")
    got = qs.from_t1(t1)
    assert got.t1_pattern == "50%%_T1.nii"
    assert resolve_pattern(got.anat, got.t1_pattern, 0) == [t1.resolve()]


def test_segmentation_in_the_anatomy_folder_leaves_path_seg_empty(tmp_path):
    t1 = _touch(tmp_path / "sub-003_T1w.nii")
    (tmp_path / "m2m_sub-003").mkdir()
    got = qs.from_t1(t1)
    assert got.m2m_found and got.seg == ""


def test_default_output_is_beside_not_inside_the_anatomy(tmp_path):
    anat = tmp_path / "data" / "anat"
    assert qs.default_output(str(anat), None) == tmp_path / "data" / "prestus_sim"
    cfg = tmp_path / "study" / "config.yaml"
    assert qs.default_output(str(anat), cfg) == (tmp_path / "study").resolve() / "prestus_sim"


def test_simnibs_bin_is_recognised_by_charm(tmp_path):
    name = "charm.cmd" if os.name == "nt" else "charm"
    _touch(tmp_path / "bin" / name)
    assert qs.is_simnibs_bin(tmp_path / "bin")
    assert not qs.is_simnibs_bin(tmp_path)


os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def test_start_page_writes_the_model(tmp_path, schema, prestus_root, monkeypatch):
    pytest.importorskip("PyQt6")
    from PyQt6.QtWidgets import QApplication

    from prestus_nexus.equipment import EquipmentLibrary
    from prestus_nexus.model import ConfigModel
    from prestus_nexus.session import Session
    from prestus_nexus.ui.main_window import START, MainWindow

    app = QApplication.instance() or QApplication([])
    monkeypatch.setattr(qs, "find_simnibs_bin", lambda: None)
    w = MainWindow(ConfigModel(schema, EquipmentLibrary(prestus_root)), Session(tmp_path / "s.yaml"),
                   prestus_root)
    assert w.tabs.widget(0) is w.pages[START] and w.tabs.currentIndex() == 0
    start = w.pages[START]
    t1 = _touch(tmp_path / "sub-004" / "anat" / "sub-004_T1w.nii.gz")
    start.apply_t1(t1)
    m = w.model
    assert m.get("subject_id") == 4 and m.get("path.t1_pattern") == "sub-%1$03d_T1w.nii.gz"
    assert m.get("path.sim") == str(tmp_path / "sub-004" / "prestus_sim")
    assert w.dirty
    start.rb_auto.click()
    assert m.get("placement.mode") == "heuristic"
    for sp, v in zip(start.focus, (-29, -6, 7)):
        sp.setValue(v)
    assert m.get("placement.heuristic.mni_target_mm") == [-29.0, -6.0, 7.0]
    start.rb_entry.click()                       # the target moves with the mode
    assert m.get("placement.mode") == "mni" and m.get("placement.mni.focus_pos_mm") == [-29.0, -6.0, 7.0]
    start.ep.setValue(50)
    assert m.get("transducer.focal_distance_ep") == 50.0
    start.ep.setValue(0)                         # 'Not set'
    assert m.get("transducer.focal_distance_ep") in ([], None)
    app.processEvents()
