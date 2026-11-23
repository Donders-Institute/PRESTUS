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
    assert ("transducer.freq_hz", "error") not in found      # mixed frequencies are allowed
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


def _subject(schema, tmp_path, *, m2m=(), simnibs=False):
    """Model for sub-003 with a T1 and the given m2m_sub-003 files."""
    anat = tmp_path / "anat"
    anat.mkdir(parents=True)
    (anat / "sub-003_T1w.nii.gz").write_bytes(b"")
    for rel in m2m:
        f = anat / "m2m_sub-003" / rel
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_bytes(b"")
    m = ConfigModel(schema)
    for key, value in {"subject_id": 3, "path.anat": str(anat), "path.sim": str(tmp_path / "sim")}.items():
        m.set(key, value)
    if simnibs:
        (tmp_path / "simnibs").mkdir()
        m.set("startup.simnibs_bin_path", str(tmp_path / "simnibs"))
    return m


SEG = ("final_tissues.nii.gz", "toMNI/MNI2Conform_nonl.nii.gz")


def test_mode_status_without_segmentation(schema, tmp_path):
    st = checks.placement_mode_status(_subject(schema, tmp_path))
    assert st["manual"].state == "ready"
    for mode in ("heuristic", "mni", "plantus", "tus_entry"):
        assert st[mode].state == "unavailable"
        assert "SimNIBS bin folder" in st[mode].reason
    assert st["localite"].state == "pending"          # no pose location set yet: selectable to set it
    assert "Localite exports folder" in st["localite"].reason


def test_localite_unavailable_names_the_folders(schema, tmp_path):
    model = _subject(schema, tmp_path)
    raw = tmp_path / "localite_raw"
    raw.mkdir()
    model.set("path.localite_raw", str(raw))
    st = checks.placement_mode_status(model)
    assert st["localite"].state == "unavailable"
    assert str(raw) in st["localite"].reason


def test_mode_status_segmentation_runs_first(schema, tmp_path):
    st = checks.placement_mode_status(_subject(schema, tmp_path, simnibs=True))
    assert st["heuristic"].state == "pending"
    assert "charm runs first" in st["heuristic"].reason and "subject" in st["heuristic"].reason


def test_mode_status_plantus_needs_mesh_and_paths(schema, tmp_path):
    m = _subject(schema, tmp_path, m2m=SEG)
    st = checks.placement_mode_status(m)
    assert st["heuristic"].state == "ready"
    assert st["plantus"].state == "unavailable"            # no sub-003.msh, no SimNIBS
    assert "sub-003.msh" in st["plantus"].reason
    m = _subject(schema, tmp_path / "b", m2m=SEG + ("sub-003.msh",))
    st = checks.placement_mode_status(m)
    assert st["plantus"].state == "pending"                # tool paths not set: still selectable
    assert "PlanTUS script" in st["plantus"].reason


def test_mode_status_localite(schema, tmp_path):
    m = _subject(schema, tmp_path)
    post = tmp_path / "post"
    (post / "sub-003").mkdir(parents=True)
    m.set("path.localite_post", str(post))
    assert checks.placement_mode_status(m)["localite"].state == "ready"


def test_placement_issues(schema, tmp_path):
    m = _subject(schema, tmp_path)
    m.set("placement.mode", "heuristic")
    found = _msgs(checks.placement_issues(m), "error")
    assert ("placement.mode", "error") in found                       # no segmentation
    m = _subject(schema, tmp_path / "b", m2m=SEG + ("sub-003.msh",))
    m.set("placement.mode", "plantus")
    found = _msgs(checks.placement_issues(m), "error")
    assert ("placement.plantus.script_path", "error") in found
    assert ("placement.plantus.focal_distance_list", "error") in found
    assert ("placement.mode", "error") not in found
    m.set("placement.mode", "tus_entry")
    m.set("placement.tus_entry.script_path", str(tmp_path / "missing.R"))
    found = _msgs(checks.placement_issues(m), "error")
    assert found[("placement.tus_entry.script_path", "error")] == "File not found"


def test_grid_resolution_issues(schema):
    from prestus_nexus.checks import STAIRCASE_MM, grid_issues
    from prestus_nexus.model import ConfigModel
    m = ConfigModel(schema)
    m.set("transducer.freq_hz", 500e3)
    m.set("simulation.medium", "layered")
    m.set("grid.resolution_mm", 0.5)                      # 6 PPW in water at 500 kHz: fine
    assert not grid_issues(m)
    m.set("grid.resolution_mm", 1.0)                      # 3 PPW, and staircasing
    msgs = [(i.level, i.message) for i in grid_issues(m)]
    assert any(lvl == "warning" and "points per wavelength" in t and "0.50 mm" in t for lvl, t in msgs)
    assert any("staircasing" in t for _lvl, t in msgs)
    m.set("acoustic_backend", "jwave")
    assert any(i.level == "error" for i in grid_issues(m))   # jwave diverges below min PPW
    m.reset("acoustic_backend")
    m.set("grid.resolution_mm", 0.6)
    m.set("simulation.medium", "water")                    # no head: no staircasing warning
    assert not any("staircasing" in i.message for i in grid_issues(m))
    m.set("simulation.medium", "phantom")
    assert 0.6 <= STAIRCASE_MM and not any("staircasing" in i.message for i in grid_issues(m))


def test_measured_ct_pattern(schema, tmp_path):
    from prestus_nexus.outputs import measured_ct_file
    anat = tmp_path / "anat"
    anat.mkdir()
    (anat / "sub-003_T1w.nii.gz").write_bytes(b"")
    m = ConfigModel(schema)
    for key, value in {"subject_id": 3, "path.anat": str(anat), "path.sim": str(tmp_path / "sim"),
                       "simulation.medium": "layered", "pct.enabled": 1,
                       "path.ct_pattern": "sub-%1$03d_CT.nii.gz"}.items():
        m.set(key, value)
    assert measured_ct_file(m) is None
    assert ("path.ct_pattern", "error") in _msgs(checks.subject_issues(m))
    (anat / "sub-003_CT.nii.gz").write_bytes(b"")
    checks.clear_status_cache()           # file probes are cached while the settings are unchanged
    assert measured_ct_file(m) == anat / "sub-003_CT.nii.gz"
    assert ("path.ct_pattern", "info") in _msgs(checks.subject_issues(m))
    # pseudo-CT generation keys apply only without a measured CT
    assert not m.is_visible("pct.skull_mapping") and m.is_visible("pct.mapping_density")
    m.set("path.ct_pattern", "")
    assert measured_ct_file(m) is None and m.is_visible("pct.skull_mapping")


def test_io_preview_ct_rows(schema, tmp_path):
    from PyQt6.QtWidgets import QApplication
    from prestus_nexus.ui.io_preview import IOPreviewPanel
    app = QApplication.instance() or QApplication([])  # noqa: F841
    m = ConfigModel(schema)
    for key, value in {"subject_id": 3, "path.anat": str(tmp_path / "anat"), "path.sim": str(tmp_path / "sim"),
                       "simulation.medium": "layered", "pct.enabled": 1}.items():
        m.set(key, value)
    panel = IOPreviewPanel(m)
    text = panel.plain_text()
    assert f"pseudo-CT  {tmp_path / 'pct' / 'sub-003'}/" in text      # next to the segmentation, not in m2m
    m.set("path.ct_pattern", "sub-%1$03d_CT.nii.gz")
    assert f"measured CT  {tmp_path / 'anat' / 'sub-003_CT.nii.gz'}" in panel.plain_text()
    panel.refresh()
    assert "nexus-icon:folder" in panel._shown


def test_path_button_icons():
    from PyQt6.QtWidgets import QApplication
    from prestus_nexus.ui.icons import PathButton
    app = QApplication.instance() or QApplication([])  # noqa: F841
    assert PathButton("file").toolTip() == "Choose a file"
    assert PathButton("dir").toolTip() == "Choose a folder" and not PathButton("dir").icon().isNull()


def test_transducers_at_different_frequencies_allowed(schema):
    """load_transducer_parameters.m runs each frequency as its own simulation: no error, for
    async and coherent coupling alike."""
    from prestus_nexus import checks
    from prestus_nexus.model import ConfigModel
    m = ConfigModel(schema)
    m.set("transducer.freq_hz", 500000)
    m.add_transducer()
    m.set("transducer.freq_hz", 300000, 1)
    for coupling in ("async", "coherent"):
        m.set("simulation.transducer_coupling", coupling)
        assert not [i for i in checks.all_issues(m) if i.path == "transducer.freq_hz" and i.level == "error"]


def test_air_cavity_needs_charm_labelling(schema, tmp_path):
    """source 'simnibs' (the default) warns when an existing m2m folder lacks charm's
    labelling, and not when it has it, charm runs first, the check is off, or the medium
    is not a head."""
    key = ("headmodel.air_cavities.source", "warning")
    m = _subject(schema, tmp_path / "a", m2m=SEG)
    m.set("simulation.medium", "layered")
    assert m.get("headmodel.air_cavities.source") == "simnibs"
    assert key in _msgs(checks.subject_issues(m))
    checks.clear_status_cache()
    m.set("io.overwrite_simnibs", 1)                    # charm re-runs and writes it
    assert key not in _msgs(checks.subject_issues(m))
    m.set("io.overwrite_simnibs", 0)
    for setting, value in (("headmodel.air_cavities.source", "hu"),
                           ("headmodel.air_cavities.source", "none"),
                           ("simulation.medium", "water")):
        m.set(setting, value)
        assert key not in _msgs(checks.subject_issues(m)), (setting, value)
        m.set("simulation.medium", "layered")
        m.set("headmodel.air_cavities.source", "simnibs")
    ok = _subject(schema, tmp_path / "b", m2m=SEG + ("segmentation/labeling.nii.gz",))
    ok.set("simulation.medium", "layered")
    assert key not in _msgs(checks.subject_issues(ok))
    none = _subject(schema, tmp_path / "c")             # no segmentation yet: charm runs first
    none.set("simulation.medium", "layered")
    assert key not in _msgs(checks.subject_issues(none))
