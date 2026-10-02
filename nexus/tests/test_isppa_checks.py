"""Issues for transducer.target_isppa_wcm2 (needs_water_baseline.m / isppa_scale_factor.m),
the affix rule (sanitize_output_affix.m) and the placement record (write_placement_positions.m)."""

import json

from prestus_nexus import checks
from prestus_nexus.model import ConfigModel
from prestus_nexus.outputs import placement_record
from prestus_nexus.study import sanitize_affix


def _custom(schema, mode, target):
    m = ConfigModel(schema)
    m.set("modules.run_water_baseline", mode)
    m.set("transducer.type", "annular")
    m.set("transducer.annular.elem_p_amp", [50000])
    m.set("transducer.target_isppa_wcm2", target)
    return m


def _levels(issues):
    return [i.level for i in issues]


def test_uncalibrated_target_with_baseline_off_is_an_error(schema):
    issues = checks.isppa_issues(_custom(schema, 0, 30))
    assert _levels(issues) == ["error"]
    assert "noBaseline" in issues[0].message


def test_uncalibrated_target_with_auto_is_info(schema):
    issues = checks.isppa_issues(_custom(schema, "auto", 30))
    assert _levels(issues) == ["info"] and "extra free-water simulation" in issues[0].message


def test_baseline_always_on_or_no_target_is_silent(schema):
    assert checks.isppa_issues(_custom(schema, 1, 30)) == []
    assert checks.isppa_issues(_custom(schema, 0, [])) == []
    assert checks.isppa_issues(_custom(schema, 0, [10, 20])) == []     # multi-ISPPA runs its own


def test_library_transducer_is_calibrated(schema):
    m = ConfigModel(schema)
    m.set("modules.run_water_baseline", 0)
    m.set("transducer.serial", "CTX-500-026")
    assert checks.isppa_issues(m) == []


def test_sanitize_affix_matches_matlab():
    assert sanitize_affix("") == ""
    assert sanitize_affix("run 1") == "_desc-run_1"
    assert sanitize_affix("_myRun") == "_desc-myRun"
    assert sanitize_affix("desc-A") == "_desc-A"            # leading desc- is not repeated
    assert sanitize_affix("_desc-A") == "_desc-A"
    assert sanitize_affix("___") == ""


def test_placement_record(schema, tmp_path):
    m = ConfigModel(schema)
    m.set("subject_id", 7)
    m.set("path.sim", str(tmp_path))
    m.set("io.output_affix", "A")
    assert placement_record(m) is None
    out = tmp_path / "sub-007"
    out.mkdir()
    (out / "sub-007_placement_desc-A.json").write_text(json.dumps({
        "subject_id": 7, "mode": "heuristic", "space": "t1_voxel_1based",
        "transducers": [{"index": 1, "trans_pos": [1, 2, 3], "focus_pos": [4, 5, 6]}]}))
    rec = placement_record(m)
    assert rec["mode"] == "heuristic" and rec["transducers"][0]["focus_pos"] == [4, 5, 6]
