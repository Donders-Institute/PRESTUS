"""Study runs from neuronav_gui / neuronav_ingest_batch outputs."""

import json

import numpy as np

from prestus_nexus import neuronav_import as ni


def _post(tmp_path):
    post = tmp_path / "localite_post"
    (post / "sub-001" / "ses-02").mkdir(parents=True, exist_ok=True)
    (post / "neuronav_ingest_summary.csv").write_text(
        "sub,ses,target,status\n"
        "sub-001,ses-02,R_MD,2 target(s) written\n"
        "sub-001,ses-02,L_MD,2 target(s) written\n"
        "sub-002,ses-02,L_MD,1 target(s) written\n"
        "sub-003,ses-02,,no usable markers\n")
    # a session written from the single-session GUI only (no summary row)
    d = post / "sub-004" / "ses-03"
    d.mkdir(parents=True, exist_ok=True)
    (d / "InstrumentMarker_sub-004_ses-03.json").write_text(json.dumps(
        {"L_MD": {"series_index": 1}, "R_MD": {"series_index": 2}}))
    return post


def test_scan_reads_summary_and_json(tmp_path):
    entries = ni.scan(_post(tmp_path))
    assert [(e.sub, e.ses, e.target, e.source) for e in entries] == [
        ("sub-001", "ses-02", "L_MD", "summary"), ("sub-001", "ses-02", "R_MD", "summary"),
        ("sub-002", "ses-02", "L_MD", "summary"),
        ("sub-004", "ses-03", "L_MD", "json"), ("sub-004", "ses-03", "R_MD", "json")]


def test_one_run_per_target(tmp_path):
    m, notes = ni.build_manifest(ni.scan(_post(tmp_path)), "target")
    assert [(r["subject_id"], r["condition"], r["placement.localite.target"]) for r in m.rows] == [
        ("1", "L_MD", "L_MD"), ("1", "R_MD", "R_MD"), ("2", "L_MD", "L_MD"), ("4", "L_MD", "L_MD"),
        ("4", "R_MD", "R_MD")]
    assert all(r["placement.mode"] == "localite" for r in m.rows) and not notes
    assert not m.issues()


def test_one_run_per_session_bilateral(tmp_path):
    m, notes = ni.build_manifest(ni.scan(_post(tmp_path)), "session", "lr", n_transducers=2)
    assert [(r["subject_id"], r["condition"], r["placement.localite.target"]) for r in m.rows] == [
        ("1", "ses-02", "[L_MD, R_MD]"), ("4", "ses-03", "[L_MD, R_MD]")]
    assert len(notes) == 1 and "sub-002" in notes[0]          # only one target: skipped
    m2, _ = ni.build_manifest(ni.scan(_post(tmp_path)), "session", ["R_MD", "L_MD"], n_transducers=2)
    assert m2.rows[0]["placement.localite.target"] == "[R_MD, L_MD]"


def test_scan_reads_qc_and_coordinates(tmp_path):
    post = tmp_path / "localite_post"
    post.mkdir()
    (post / "neuronav_ingest_summary.csv").write_text(
        "sub,ses,target,pos_dev_mm,rot_dev_rad,axis_dist_mm,depth_mm,n_pulses,"
        "trans_ras_x,trans_ras_y,trans_ras_z,focus_ras_x,focus_ras_y,focus_ras_z\n"
        "sub-001,ses-02,L_MD,0.42,0.01,3.1,68.5,80,43.35,-16.17,82.3,5.45,-20.71,5.19\n"
        "sub-001,ses-02,R_MD,0.5,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN\n")
    entries = {e.target: e for e in ni.scan(post)}
    l_md = entries["L_MD"]
    assert l_md.pos_dev_mm == 0.42 and l_md.rot_dev_rad == 0.01 and l_md.n_pulses == 80
    np.testing.assert_allclose(l_md.trans_ras, [43.35, -16.17, 82.3])
    np.testing.assert_allclose(l_md.focus_ras, [5.45, -20.71, 5.19])
    # a NaN cell (MATLAB writetable's literal 'NaN') is missing, not a parsed float
    r_md = entries["R_MD"]
    assert r_md.pos_dev_mm == 0.5 and r_md.rot_dev_rad is None and r_md.trans_ras is None

    found = ni.coordinates_for(post, "ses-02", "L_MD", subject_id=1)
    assert found is not None and found.target == "L_MD"
    np.testing.assert_allclose(found.trans_ras, [43.35, -16.17, 82.3])
    by_int = ni.coordinates_for(post, 2, "L_MD", subject_id=1)          # session as plain int
    assert by_int is not None and by_int.sub == found.sub and by_int.ses == found.ses
    assert ni.coordinates_for(post, "ses-02", "L_MD", subject_id=2) is None   # wrong subject
    assert ni.coordinates_for(post, "ses-02", "NOPE", subject_id=1) is None


def test_order_targets():
    assert ni.order_targets(["R_PUL", "L_PUL", "M1"], "lr") == ["L_PUL", "R_PUL", "M1"]
    assert ni.order_targets(["b", "a"], "alphabetical") == ["a", "b"]
    assert ni.order_targets(["x", "y", "z"], ["z", "x"]) == ["z", "x", "y"]
