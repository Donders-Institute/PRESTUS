"""Reads what transducer_positioning.m / tp_select_heuristic_position.m write to
the heuristic placement run folder (doc_placement_heuristic.md)."""

import numpy as np

from prestus_nexus import heuristic_import as hi


def test_candidates_reads_tpars_csv(tmp_path):
    (tmp_path / "tpars_sub-009_M1.csv").write_text(
        "idx,trans_x,trans_y,trans_z,targ_x,targ_y,targ_z,dist_to_target,prop_intersect,"
        "mean_dist_skin,var_dist_skin,mean_dist_skull,var_dist_skull\n"
        "1,161,138,246,117,139,161,85.2,0.02,3.1,0.4,5.2,0.6\n"
        "2,160,137,245,117,139,161,84.9,0.31,2.9,0.5,5.0,0.7\n")
    cands = hi.candidates(tmp_path, 9, "M1")
    assert len(cands) == 2
    assert cands[0].idx == 1
    np.testing.assert_allclose(cands[0].trans_vox, [161, 138, 246])
    assert cands[1].prop_intersect == 0.31
    assert hi.candidates(tmp_path, 9, "NOPE") == []


def test_selected_reads_key_value_tsv(tmp_path):
    (tmp_path / "sub-009_M1.txt").write_text(
        "trans_x\t161\n"
        "trans_y\t138\n"
        "trans_z\t246\n"
        "targ_x\t117\n"
        "targ_y\t139\n"
        "targ_z\t161\n"
        "LOCALITE_Transducer_mm\tInf\n"
        "loc_trans_x_mm\t144.9\n"
        "RAS_Transducer_mm\tInf\n"
        "ras_trans_x_mm\t43.35\n"
        "ras_trans_y_mm\t-16.17\n"
        "ras_trans_z_mm\t82.3\n"
        "RAS_Target_mm\tInf\n"
        "ras_targ_x_mm\t5.45\n"
        "ras_targ_y_mm\t-20.71\n"
        "ras_targ_z_mm\t5.19\n")
    sel = hi.selected(tmp_path, 9, "M1")
    assert sel is not None
    np.testing.assert_allclose(sel.trans_ras, [43.35, -16.17, 82.3])
    np.testing.assert_allclose(sel.focus_ras, [5.45, -20.71, 5.19])
    np.testing.assert_allclose(sel.trans_vox, [161, 138, 246])
    assert hi.selected(tmp_path, 9, "NOPE") is None


def test_selected_without_ras_columns_falls_back_to_voxel(tmp_path):
    (tmp_path / "sub-001_X.txt").write_text("trans_x\t1\ntrans_y\t2\ntrans_z\t3\n"
                                            "targ_x\t4\ntarg_y\t5\ntarg_z\t6\n")
    sel = hi.selected(tmp_path, 1, "X")
    assert sel is not None and sel.trans_ras is None
    np.testing.assert_allclose(sel.trans_vox, [1, 2, 3])


def test_dir_placement_run_mirrors_placement_run_dir(schema, tmp_path):
    from prestus_nexus.model import ConfigModel
    from prestus_nexus.outputs import dir_placement_run

    m = ConfigModel(schema)
    m.set("subject_id", 9)
    assert dir_placement_run(m, "heuristic") is None              # no path.sim / path.placement
    m.set("path.sim", str(tmp_path / "sims"))
    m.set("placement.heuristic.target_name", "M1")
    assert dir_placement_run(m, "heuristic") == tmp_path / "sims" / "placement" / "sub-009" / "heuristic" / "M1"
    m.set("path.placement", str(tmp_path / "shared"))
    assert dir_placement_run(m, "plantus") == tmp_path / "shared" / "sub-009" / "plantus" / "target"


def test_dir_pct_is_next_to_seg_and_finds_legacy(schema, tmp_path):
    from prestus_nexus.model import ConfigModel
    from prestus_nexus.outputs import dir_pct, pseudoct_file

    m = ConfigModel(schema)
    m.set("subject_id", 9)
    m.set("path.seg", str(tmp_path / "simnibs"))
    assert dir_pct(m) == tmp_path / "pct" / "sub-009"              # mirrors get_pct_dir.m
    assert pseudoct_file(m) == (None, False)
    legacy = tmp_path / "simnibs" / "m2m_sub-009" / "pseudoCT.nii.gz"
    legacy.parent.mkdir(parents=True)
    legacy.write_bytes(b"")
    assert pseudoct_file(m) == (legacy, True)
    new = tmp_path / "pct" / "sub-009" / "pseudoCT.nii.gz"
    new.parent.mkdir(parents=True)
    new.write_bytes(b"")
    assert pseudoct_file(m) == (new, False)
    m.set("path.pct", str(tmp_path / "pct_k"))
    assert dir_pct(m) == tmp_path / "pct_k" / "sub-009"
