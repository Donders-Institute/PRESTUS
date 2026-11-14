"""Follow-up runs (followup.py): keys per type, validation, the lineage header, pose tables,
lineage discovery in the study folder, and the launcher's pipeline options."""

import json
from pathlib import Path

import pytest

from prestus_nexus import checks, followup as fu, launcher
from prestus_nexus import results as rs
from prestus_nexus.model import ConfigModel, parse_lineage
from prestus_nexus.session import Session

TAG, MED = "sub-001", "layered"


def _parent(tmp_path, schema, affix="", sidecar=True, level=30.0, extra=""):
    sims = tmp_path / "sims"
    base = sims / TAG
    (base / "nii").mkdir(parents=True)
    raw = affix[len("_desc-"):] if affix else ""
    study = tmp_path / "config_pilot.yaml"
    study.write_text(f"subject_id: 1\npath:\n  sim: {sims}\nio:\n  output_affix: '{raw}'\n"
                     "transducer:\n  target_isppa_wcm2: 30\n" + extra)
    model = ConfigModel.from_file(schema, study)
    (base / f"{TAG}_{MED}{affix}.csv").write_text("Isppa\n1\n")
    nii = base / "nii" / f"{TAG}_{MED}_T1w{affix}_pressure.nii.gz"
    nii.write_bytes(b"")
    if sidecar:
        sig = {"signature_version": 1, "medium": "layered", "subject_id": 1,
               "grid": {"resolution_mm": model.get("grid.resolution_mm")}}
        (base / "nii" / f"{TAG}_{MED}_T1w{affix}_pressure.json").write_text(json.dumps(
            {"sidecar_version": 1, "kind": "acoustic_pressure",
             "isppa": {"level_wcm2": level}, "signature": sig}))
    run = rs.discover_runs(base, TAG)[0]
    return model, run


def test_pressure_paths_follow_matlab_names(tmp_path, schema):
    _, run = _parent(tmp_path, schema, "_desc-pilot")
    assert fu.pressure_nifti(run).name == "sub-001_layered_T1w_desc-pilot_pressure.nii.gz"
    assert fu.sidecar_of(fu.pressure_nifti(run)).name == "sub-001_layered_T1w_desc-pilot_pressure.json"
    assert [p.name for p in fu.thermal_end_maps(run)] == [
        "sub-001_layered_T1w_desc-pilot_heating_end.nii.gz",
        "sub-001_layered_T1w_desc-pilot_CEM43_end.nii.gz",
        "sub-001_layered_T1w_desc-pilot_CEM43_iso_end.nii.gz"]


def test_new_intensity_keys(tmp_path, schema):
    parent, run = _parent(tmp_path, schema, "_desc-pilot")
    child = fu.create_followup(parent, run, "new_intensity", {"target_isppa_wcm2": 12.5})
    assert child.get("transducer.target_isppa_wcm2") == 12.5
    assert child.get("path.acoustic_input") == str(fu.pressure_nifti(run))
    assert child.get("io.output_affix") == "pilot_isppa12p5"
    assert child.get("modules.run_posthoc_water_sims") == 0
    assert child.get("subject_id") == 1                      # parent overrides copied
    assert child.lineage["parent_affix"] == "_desc-pilot" and child.lineage["type"] == "new_intensity"
    assert parent.get("path.acoustic_input") in ("", None)  # parent untouched
    issues = fu.validate_followup(parent, run, "new_intensity", child)
    assert not [i for i in issues if i.level == "error"], issues


def test_new_timing_and_placement_and_sequential_keys(tmp_path, schema):
    parent, run = _parent(tmp_path, schema)
    t = fu.create_followup(parent, run, "new_timing", {"timing": {"pd": 0.02, "timing.pri": 0.1}})
    assert (t.get("timing.pd"), t.get("timing.pri"), t.get("modules.run_heating_sims")) == (0.02, 0.1, 1)
    assert t.get("path.acoustic_input").endswith("_T1w_pressure.nii.gz")
    assert t.get("io.output_affix") == "timing"
    pl = fu.create_followup(parent, run, "new_placement_same_head",
                            {"changes": {"transducer.trans_pos": [1, 2, 3]}})
    assert pl.get("headmodel.reuse_across_placements") is True
    assert pl.get("transducer.trans_pos") == [1, 2, 3]
    assert pl.get("path.acoustic_input") in ("", None)
    seq = fu.create_followup(parent, run, "sequential", {"seq_index": 3})
    assert seq.get("io.output_affix") == "seq3"
    assert seq.get("modules.run_heating_sims") == 1 and seq.get("modules.run_nifti_creation") == 1


def test_validation_errors(tmp_path, schema):
    parent, run = _parent(tmp_path, schema, sidecar=False)
    child = fu.create_followup(parent, run, "new_intensity", {"target_isppa_wcm2": 10, "affix": ""})
    child.set("io.output_affix", "")                         # same affix as the parent
    errs = {i.path for i in fu.validate_followup(parent, run, "new_intensity", child) if i.level == "error"}
    assert "io.output_affix" in errs
    assert "path.acoustic_input" in errs                      # no sidecar: cannot rescale
    child.set("simulation.medium", "water")
    child.set("grid.axisymmetric", 1)
    child.set("thermal_backend", "babelbrain")
    errs = {i.path for i in fu.validate_followup(parent, run, "new_intensity", child) if i.level == "error"}
    assert {"simulation.medium", "grid.axisymmetric", "thermal_backend"} <= errs
    # sequential needs a parent with a thermal stage
    seq = fu.create_followup(parent, run, "sequential")
    errs = {i.path for i in fu.validate_followup(parent, run, "sequential", seq) if i.level == "error"}
    assert "modules.run_heating_sims" in errs


def test_multi_transducer_target_refused(tmp_path, schema):
    parent, run = _parent(tmp_path, schema)
    parent.add_transducer()
    child = fu.create_followup(parent, run, "new_intensity", {"target_isppa_wcm2": 10})
    msgs = [i.message for i in fu.validate_followup(parent, run, "new_intensity", child) if i.level == "error"]
    assert any("combined pressure map" in m for m in msgs)


def test_signature_diffs(tmp_path, schema):
    parent, _ = _parent(tmp_path, schema)
    eff = parent.effective()
    stored = {"medium": "layered", "subject_id": 1, "grid": {"resolution_mm": eff["grid"]["resolution_mm"],
                                                             "thermal_x": 9},
              "transducer": {"trans_pos": [9, 9, 9], "freq_hz": eff["transducer"]["freq_hz"]}}
    assert fu.signature_diffs(stored, eff) == []
    stored["grid"]["resolution_mm"] = 99
    stored["subject_id"] = 2
    assert fu.signature_diffs(stored, eff) == ["grid.resolution_mm", "subject_id"]


def test_adopted_heatmap_missing_is_an_error(tmp_path, schema):
    parent, _ = _parent(tmp_path, schema)
    parent.set("io.adopted_heatmap", str(tmp_path / "nope.nii.gz"))
    assert [i.path for i in fu.adopted_issues(parent)] == ["io.adopted_heatmap"]
    assert any(i.path == "io.adopted_heatmap" and i.level == "error" for i in checks.all_issues(parent))


def test_acoustic_input_check(tmp_path, schema):
    parent, run = _parent(tmp_path, schema)
    parent.set("path.acoustic_input", str(tmp_path / "missing_pressure.nii.gz"))
    assert any(i.path == "path.acoustic_input" and i.level == "error" for i in checks.all_issues(parent))
    checks.clear_status_cache()
    parent.set("path.acoustic_input", str(fu.pressure_nifti(run)))
    issues = [i for i in checks.all_issues(parent) if i.path == "path.acoustic_input"]
    assert not [i for i in issues if i.level == "error"], issues


def test_lineage_round_trip(tmp_path, schema):
    parent, run = _parent(tmp_path, schema)
    child = fu.create_followup(parent, run, "new_intensity", {"target_isppa_wcm2": 10})
    before = child.override_dict()
    path = fu.save_followup(child, tmp_path / fu.followup_filename(parent.path, "new_intensity",
                                                                    child.get("io.output_affix")),
                            header="PRESTUS study config")
    assert path.name == "config_pilot_new-intensity_isppa10.yaml"
    text = path.read_text()
    assert text.startswith("# PRESTUS study config\n# nexus-followup: {")
    loaded = ConfigModel.from_file(schema, path)
    assert loaded.lineage == child.lineage
    assert loaded.override_dict() == before                  # the header holds no config key
    loaded.save(header="again")                              # survives a later save
    assert parse_lineage(path.read_text()) == child.lineage
    assert "nexus-followup" not in ConfigModel.from_file(schema, parent.path).to_yaml("x")
    # a missing parent is reported
    loaded.lineage = dict(loaded.lineage, parent_config=str(tmp_path / "gone.yaml"))
    assert fu.lineage_issues(loaded)


def test_lineage_discovery(tmp_path, schema):
    parent, run = _parent(tmp_path, schema, "_desc-pilot")
    child = fu.create_followup(parent, run, "new_intensity", {"target_isppa_wcm2": 10})
    fu.save_followup(child, tmp_path / "child.yaml")
    links = rs.lineage_links(tmp_path)
    assert len(links) == 1 and links[0].child_affix == "_desc-pilot_isppa10"
    up, down = rs.lineage_for(run, links, parent.path)
    assert up is None and [d.child_config.name for d in down] == ["child.yaml"]
    child_run = rs.Run(affix="_desc-pilot_isppa10", medium=MED, base=run.base, tag=TAG)
    up, down = rs.lineage_for(child_run, links)
    assert up is not None and up.parent_affix == "_desc-pilot" and down == []


def test_pose_table_and_pose_followup(tmp_path, schema):
    parent, run = _parent(tmp_path, schema, extra=(
        "  type: annular\n  annular:\n    curv_radius_mm: 63\n    dist_geom_ep_mm: 50\n"
        "    elem_od_mm: [10, 20]\n    elem_id_mm: [0, 12]\n"))
    head = ["pose", "kind", "polar_deg", "p_corrected_rel", "complete"] + list(fu.POSE_RAS_COLUMNS)
    row = ["7", "orbit", "5", "1.2", "true", "10", "0", "0", "40", "0", "0", "41", "0", "0"]
    csv_path = run.base / f"{TAG}_{MED}_timereversal_poses_T1.csv"
    csv_path.write_text(",".join(head) + "\n" + ",".join(row) + "\n")
    assert fu.pose_tables(run) == [csv_path]
    pose = fu.read_pose_table(csv_path)[0]
    assert pose["pose"] == 7 and pose["kind"] == "orbit" and pose["complete"] is True
    child = fu.create_followup(parent, run, "pose_from_time_reversal", {"pose": pose})
    assert child.get("transducer.trans_pos") == [10.0, 0.0, 0.0]
    assert child.get("transducer.focus_pos") == [40.0, 0.0, 0.0]        # corrected aim point
    assert child.get("transducer.focal_distance_ep") == pytest.approx(30 - 13)   # minus bowl depth
    assert child.get("placement.mode") == "manual"
    assert child.get("placement.manual.coordinate_space") == "ras"
    assert child.get("time_reversal.enabled") is False
    assert child.get("io.output_affix") == "pose007"
    geo = fu.create_followup(parent, run, "pose_from_time_reversal", {"pose": pose, "aim": "geometric"})
    assert geo.get("transducer.focus_pos") == [41.0, 0.0, 0.0]
    ph = fu.create_followup(parent, run, "pose_from_time_reversal", {"aim": "apply_phases"})
    assert ph.get("time_reversal.apply_phases") is True and ph.get("time_reversal.enabled") is True


# ── launcher options ───────────────────────────────────────────────────────────

def test_launcher_unchanged_without_options(tmp_path, prestus_root):
    study = tmp_path / "config_pilot.yaml"
    study.write_text("subject_id: 1\n")
    plain = launcher.launcher_script(study, prestus_root)
    assert launcher.launcher_script(study, prestus_root, options=launcher.PipelineOptions()) == plain
    assert plain.splitlines()[-1] == "prestus_pipeline_start(parameters);"


def test_launcher_with_sequential_and_stage_options(tmp_path, prestus_root):
    study = tmp_path / "config_pilot.yaml"
    study.write_text("subject_id: 1\n")
    a, b = tmp_path / "c2.yaml", tmp_path / "sub" / "c3.yaml"
    opts = launcher.PipelineOptions([str(a), str(b)], {"thermal_timelimit": "03:00:00",
                                                       "sim_memorylimit": 40}, False)
    text = launcher.launcher_script(study, prestus_root, options=opts)
    assert "options = struct();" in text
    assert "options.thermal_timelimit = '03:00:00';" in text
    assert "options.sim_memorylimit = 40;" in text
    assert "options.generate_summary_report = false;" in text
    assert f"options.sequential_configs.config_2 = load_parameters('c2.yaml', '{a.parent.resolve()}');" in text
    assert "options.sequential_configs.config_3 = load_parameters('c3.yaml'" in text
    assert "options.sequential_configs.config_2.platform = 'matlab';" in text
    assert text.rstrip().endswith("prestus_pipeline_start(parameters, options);")
    hpc = launcher.launcher_script(study, prestus_root, local=False, options=opts)
    assert ".platform = 'matlab'" not in hpc


def test_options_stored_in_session(tmp_path):
    s = Session(tmp_path / "session.yaml")
    study = tmp_path / "config_pilot.yaml"
    assert launcher.load_options(s, study).is_empty()
    launcher.add_to_chain(s, study, tmp_path / "child.yaml")
    launcher.add_to_chain(s, study, tmp_path / "child.yaml")
    opts = launcher.load_options(s, study)
    assert opts.sequential_configs == [str((tmp_path / "child.yaml").resolve())]
    launcher.save_options(s, study, launcher.PipelineOptions())
    assert launcher.load_options(s, study).is_empty()


# ── UI ─────────────────────────────────────────────────────────────────────────

@pytest.fixture
def qapp():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_dialog_creates_followup(tmp_path, schema, qapp):
    from types import SimpleNamespace
    from PyQt6.QtWidgets import QWidget
    from prestus_nexus.ui.followup_dialog import FollowupDialog
    parent, run = _parent(tmp_path, schema, "_desc-pilot")
    host = QWidget()
    win = SimpleNamespace(model=parent, session=Session(tmp_path / "s.yaml"))
    dlg = FollowupDialog(win, run, "new_intensity", parent=host)
    dlg.isppa.setValue(15)
    dlg._update()
    assert dlg.affix.text() == "pilot_isppa15"
    assert dlg.create_btn.isEnabled(), [i.message for i in dlg.issues]
    got = []
    dlg.created.connect(got.append)
    path = dlg.create()
    assert path == tmp_path / "config_pilot_new-intensity_pilot-isppa15.yaml" and got == [path]
    assert ConfigModel.from_file(schema, path).lineage["type"] == "new_intensity"
    # sequential: the child joins the parent's chain
    seq = FollowupDialog(win, run, "sequential", parent=host)
    seq._update()
    assert not seq.create_btn.isEnabled()                     # parent has no thermal stage
    parent.set("modules.run_heating_sims", 1)
    seq._update()
    p2 = seq.create()
    assert launcher.load_options(win.session, parent.path).sequential_configs == [str(p2.resolve())]


def test_results_page_followup_and_poses(tmp_path, schema, qapp):
    from types import SimpleNamespace
    from prestus_nexus.ui.results_page import ResultsPage
    parent, run = _parent(tmp_path, schema)
    head = ["pose", "kind", "p_corrected_rel"] + list(fu.POSE_RAS_COLUMNS)
    (run.base / f"{TAG}_{MED}_timereversal_poses_T1.csv").write_text(
        ",".join(head) + "\n3,orbit,1.1," + ",".join(["1"] * 9) + "\n")
    child = fu.create_followup(parent, run, "new_intensity", {"target_isppa_wcm2": 10})
    fu.save_followup(child, tmp_path / "child.yaml")
    page = ResultsPage(SimpleNamespace(model=parent, session=Session(tmp_path / "s.yaml")))
    page.refresh(force=True)
    assert page.followup_btn.isEnabled()
    assert "Follow-ups of this run" in page.lineage.text() and "child.yaml" in page.lineage.text()
    assert page.subtabs.isTabVisible(page.subtabs.indexOf(page.poses))
    assert page.subtabs.tabText(page.subtabs.indexOf(page.table)) == "Data table"
    page.subtabs.setCurrentWidget(page.poses)
    assert page.poses.rows[0]["pose"] == 3 and page.poses.use.isEnabled()
    dlg = page.open_followup("pose_from_time_reversal", page.poses.rows[0])
    dlg._update()
    assert dlg.child.get("transducer.trans_pos") == [1.0, 1.0, 1.0]
    dlg.close()


def test_pipeline_options_panel(tmp_path, schema, qapp):
    from types import SimpleNamespace
    from prestus_nexus.ui.pipeline_options import PipelineOptionsPanel
    parent, _ = _parent(tmp_path, schema)
    win = SimpleNamespace(model=parent, session=Session(tmp_path / "s.yaml"))
    launcher.add_to_chain(win.session, parent.path, tmp_path / "a.yaml")
    panel = PipelineOptionsPanel(win)
    assert panel.chain.count() == 1 and not panel.multi_box.isVisibleTo(panel)
    parent.set("simulation.uncertainty", True)
    panel.refresh()
    assert panel.unc_box.isVisibleTo(panel.body)
    panel.limits["sim"][0].setText("08:00:00")
    panel._store()
    assert launcher.load_options(win.session, parent.path).stages == {"sim_timelimit": "08:00:00"}


def test_sweep_problems():
    from prestus_nexus.ui.isppa_sweep import sweep_problems
    assert sweep_problems([10, 20]) == []
    assert sweep_problems([10, 10])
    assert sweep_problems([10, -1])
    assert sweep_problems([10, 20], 2, "coherent") and not sweep_problems([10, 20], 2, "async")
