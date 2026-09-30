"""Study runs from neuronav_gui / neuronav_ingest_batch outputs."""

import json

from prestus_planner import neuronav_import as ni


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


def test_order_targets():
    assert ni.order_targets(["R_PUL", "L_PUL", "M1"], "lr") == ["L_PUL", "R_PUL", "M1"]
    assert ni.order_targets(["b", "a"], "alphabetical") == ["a", "b"]
    assert ni.order_targets(["x", "y", "z"], ["z", "x"]) == ["z", "x", "y"]
