from prestus_nexus import checks
from prestus_nexus.model import ConfigModel
from prestus_nexus.paths import matlab_sprintf


def _msgs(issues, level=None):
    return {(i.path, i.level): i.message for i in issues if level in (None, i.level)}


def test_matlab_sprintf():
    assert matlab_sprintf("sub-%1$03d_T1w.nii*", 7) == "sub-007_T1w.nii*"
    assert matlab_sprintf("sub-%03d/%03d", 7, 8) == "sub-007/008"
    assert matlab_sprintf("100%%_%d", 3) == "100%_3"


def test_subject_inputs(schema, tmp_path):
    anat = tmp_path / "anat"
    anat.mkdir()
    (anat / "sub-003_T1w.nii.gz").write_bytes(b"")
    (tmp_path / "anat" / "m2m_sub-003").mkdir()
    m = ConfigModel(schema)
    for key, value in {"subject_id": 3, "path.anat": str(anat), "path.sim": str(tmp_path / "sim")}.items():
        m.set(key, value)
    found = _msgs(checks.subject_issues(m))
    assert "sub-003_T1w.nii.gz" in found[("path.t1_pattern", "info")]
    assert ("path.seg", "info") in found                    # m2m found (seg defaults to anat)
    assert not [i for i in checks.subject_issues(m) if i.level == "error"]

    m.set("subject_id", 4)
    found = _msgs(checks.subject_issues(m), "error")
    assert ("path.t1_pattern", "error") in found
    assert ("startup.simnibs_bin_path", "error") in found   # no m2m and no SimNIBS


def test_relative_sim_path_warns(schema):
    m = ConfigModel(schema)
    m.set("path.sim", "sims")
    assert ("path.sim", "warning") in _msgs(checks.subject_issues(m))


def test_coordinate_space_requires_placement_stage(schema):
    m = ConfigModel(schema)
    m.set("placement.manual.coordinate_space", "ras")
    m.set("modules.run_transducer_placement", 0)
    assert ("modules.run_transducer_placement", "error") in _msgs(checks.rule_issues(m))


def test_transducer_rules(schema):
    m = ConfigModel(schema)
    m.set("transducer.freq_hz", 250e3)
    m.add_transducer()
    m.set("transducer.freq_hz", 500e3, 1)
    m.set("transducer.serial", "CTX_500_026", 1)           # not on transducer 1
    found = _msgs(checks.rule_issues(m), "error")
    assert ("transducer.freq_hz", "error") in found
    assert ("transducer.serial", "error") in found


def test_babelbrain_geometry_rules(schema):
    m = ConfigModel(schema)
    m.set("acoustic_backend", "babelbrain")
    assert ("acoustic_backend", "error") not in _msgs(checks.rule_issues(m))
    m.add_transducer()
    assert ("acoustic_backend", "error") in _msgs(checks.rule_issues(m))
    m = ConfigModel(schema)
    m.set("acoustic_backend", "babelbrain")
    m.set("grid.mode", "ras_plus")
    assert ("grid.mode", "error") in _msgs(checks.rule_issues(m))
    m.set("acoustic_backend", "kwave")
    assert ("grid.mode", "error") not in _msgs(checks.rule_issues(m))


def test_coherent_sweep_rule(schema):
    m = ConfigModel(schema)
    m.set("transducer.target_isppa_wcm2", [10, 20])
    assert ("transducer.target_isppa_wcm2", "error") not in _msgs(checks.rule_issues(m))   # one transducer
    m.add_transducer()
    assert ("transducer.target_isppa_wcm2", "error") in _msgs(checks.rule_issues(m))
    m.set("simulation.transducer_coupling", "async")
    assert ("transducer.target_isppa_wcm2", "error") not in _msgs(checks.rule_issues(m))


def test_value_checks(schema):
    m = ConfigModel(schema)
    m.load_dict({"transducer": {"trans_pos": [1, 2]}, "simulation": {"medium": "vacuum"},
                 "subject_id": 12.0})
    found = _msgs(checks.value_issues(m), "error")
    assert ("transducer.trans_pos", "error") in found
    assert ("simulation.medium", "error") in found
    assert ("subject_id", "error") not in found            # 12.0 is an integer to MATLAB
