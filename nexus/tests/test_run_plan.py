"""run_plan.plan: the stages prestus_pipeline.m runs for a config."""

from prestus_nexus import run_plan as rp
from prestus_nexus.model import ConfigModel


def _stages(m):
    return {s.name: s for s in rp.plan(m).stages}


def _layered(schema, tmp_path):
    anat = tmp_path / "anat"
    (anat / "m2m_sub-003").mkdir(parents=True)
    m = ConfigModel(schema)
    for k, v in {"subject_id": 3, "path.anat": str(anat), "simulation.medium": "layered"}.items():
        m.set(k, v)
    return m, anat


def test_default_layered(schema, tmp_path):
    m, anat = _layered(schema, tmp_path)
    st = _stages(m)
    assert st["Segmentation"].state == rp.RUN                 # m2m folder without final_tissues
    (anat / "m2m_sub-003" / "final_tissues.nii.gz").write_bytes(b"")
    st = _stages(m)
    assert st["Segmentation"].state == rp.REUSE
    assert st["Skull CT"].state == rp.SKIP
    assert st["Thermal simulation"].state == rp.SKIP           # run_heating_sims: 0 by default
    assert st["Acoustic simulation"].state == rp.RUN
    m.set("modules.run_heating_sims", 1)
    assert _stages(m)["Thermal simulation"].state == rp.RUN


def test_skull_ct_sources(schema, tmp_path):
    m, anat = _layered(schema, tmp_path)
    m.set("pct.enabled", 1)
    assert _stages(m)["Skull CT"].state == rp.RUN              # pseudo-CT generated
    m.set("path.ct_pattern", "sub-%1$03d_CT.nii.gz")
    assert _stages(m)["Skull CT"].state == rp.STOP             # pattern matches nothing
    (anat / "sub-003_CT.nii.gz").write_bytes(b"")
    s = _stages(m)["Skull CT"]
    assert s.state == rp.REUSE and "measured CT" in s.detail


def test_segmentation_only_and_water(schema, tmp_path):
    m, _ = _layered(schema, tmp_path)
    m.set("modules.segmentation_only", 1)
    names = [s.name for s in rp.plan(m).stages]
    assert names == ["Segmentation", "Everything after segmentation"]
    w = ConfigModel(schema)
    w.set("simulation.medium", "water")
    st = _stages(w)
    assert st["Segmentation"].state == rp.SKIP and "Skull CT" not in st and "Targeting analysis" not in st


def test_plan_card_html(schema, tmp_path):
    from PyQt6.QtWidgets import QApplication
    from prestus_nexus.ui.run_plan_card import RunPlanCard, plan_html
    app = QApplication.instance() or QApplication([])  # noqa: F841  (kept alive for the widget)
    m, _ = _layered(schema, tmp_path)
    card = RunPlanCard()
    text = plan_html(m, card.view.document())
    assert "Mode: single" in text and "Acoustic simulation" in text and "Always written" in text
    assert "nexus-icon:folder" in text and "nexus-icon:file" in text
    card.update_plan(m)
    plain = card.view.toPlainText()
    assert "Stages" in plain and "Saved" in plain and "Not saved" in plain


def test_save_items(schema):
    from prestus_nexus.predict_outputs import save_items
    m = ConfigModel(schema)
    m.set("subject_id", 3)
    items = {i.label: i for i in save_items(m)}
    src = items["source matrices"]
    assert src.saved and src.where == "cache/sub-003_layered_kwave_source.mat" and "default on" in src.reason
    assert not items["acoustic matrices"].saved
    assert not items["thermal matrices"].applies                       # no heating by default
    assert items["MNI-space maps"].saved                                # save_MNI: 1
    m.set("io.save_MNI", [])
    assert items["MNI-space maps"].saved == {i.label: i for i in save_items(m)}["MNI-space maps"].saved
    m.set("io.save_matrices", 1)
    acoustic = {i.label: i for i in save_items(m)}["acoustic matrices"]
    assert acoustic.saved and "follows io.save_matrices: 1" in acoustic.reason


def test_existing_outputs(schema, tmp_path):
    from prestus_nexus.outputs import existing_outputs
    m = ConfigModel(schema)
    assert existing_outputs(m) is None                                # no path.sim / subject yet
    for k, v in {"subject_id": 3, "path.sim": str(tmp_path / "sim")}.items():
        m.set(k, v)
    found = existing_outputs(m)
    assert found and not any(found.values())
    out = tmp_path / "sim" / "sub-003"
    (out / "cache").mkdir(parents=True)
    (out / "cache" / "sub-003_layered_kwave_source.mat").write_bytes(b"")
    (out / "nii").mkdir()
    (out / "nii" / "sub-003_layered_MNI_Itp.nii.gz").write_bytes(b"")
    found = existing_outputs(m)
    assert found["source matrices"] and found["MNI-space maps"]
    assert not found["thermal matrices"] and not found["debug diagnostics"]
    m.set("io.output_affix", "_run2")                                # another affix: not the same files
    assert not existing_outputs(m)["source matrices"]


def test_ute_lines_from_script(schema, prestus_root, tmp_path):
    from prestus_nexus.pct_lines import read_script, ute_line
    r = read_script(prestus_root)
    assert r["lines"]["wiesinger"] == (-2000.0, 2042.0) and r["soft"] == 42.0 and r["air"] == -1000.0
    m = ConfigModel(schema)
    m.set("pct.skull_mapping", "kosciessa")
    assert ute_line(m, prestus_root, tmp_path).slope is None          # no fit yet
    (tmp_path / "debug_pct").mkdir()
    (tmp_path / "debug_pct" / "pCT_skull_mapping.txt").write_text("-2100.5\n2300\n")
    line = ute_line(m, prestus_root, tmp_path)
    assert (line.slope, line.intercept) == (-2100.5, 2300.0)
    m.set("pct.skull_mapping", "custom")
    m.set("pct.custom_slope", -1500)
    m.set("pct.custom_intercept", 1800)
    assert (ute_line(m, prestus_root, None).slope, ute_line(m, prestus_root, None).intercept) == (-1500.0, 1800.0)


def test_stage_icons_exist(schema, tmp_path):
    from PyQt6.QtWidgets import QApplication
    from prestus_nexus.ui import icons
    app = QApplication.instance() or QApplication([])  # noqa: F841
    m, _ = _layered(schema, tmp_path)
    for k, v in {"modules.run_heating_sims": 1, "time_reversal.enabled": True, "placement.mode": "localite"}.items():
        m.set(k, v)
    stages = rp.plan(m).stages
    assert all(s.icon in icons.names() for s in stages), [s.name for s in stages if s.icon not in icons.names()]
    assert {s.name: s.icon for s in stages}["Placement"] == "place_localite"


def test_predicted_load_table():
    """Total first, one table row per simulation and one for the other stages, and the GPU
    list of stages timed on 'any GPU' said once."""
    import html
    from prestus_nexus import workload as w
    from prestus_nexus.ui import style
    from prestus_nexus.ui.run_plan_card import _merge_bases, prediction_html
    rate_a = w.Rate(1.0, 0.5, 2.0, "any GPU (A100, L40S); set hpc.gpu to narrow", 10.0)
    rate_t = w.Rate(1.0, 0.5, 2.0, "any GPU (A100, P100); set hpc.gpu to narrow", 5.0)
    pred = w.Prediction(parts=[w.Part("Acoustic sim", 3, 1e8, 10_000, rate_a),
                               w.Part("Thermal sim", 1, 9e7, 260, rate_t)],
                        stages=[w.Fixed("Placement", 30, 20, 60, "runs")],
                        untimed=["segmentation (SimNIBS charm)"], gpu_mem_gb=20)
    out = prediction_html(pred, style.active(), html.escape)
    assert out.index("Total") < out.index("Acoustic sim") < out.index("Thermal sim") < out.index("Other stages")
    assert out.count("Timing basis") == 1 and "A100, L40S, P100" in out
    assert "GPU memory ~20 GB" in out and "Not timed, add: segmentation" in out
    assert _merge_bases(["measured on X"]) == ["measured on X"]          # other bases unchanged
    assert prediction_html(w.Prediction(), style.active(), html.escape) == ""
