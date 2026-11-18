"""predict_outputs.py mirrors functions/helper/predict_outputs.m (and the mode
predicates is_async_mode.m / is_multi_isppa_mode.m it calls). The nested
uncertainty cases have MATLAB counterparts in
tests/test_uncertainty_nesting.m (predictOutputsNested) — if one side changes,
the other should too."""

from prestus_nexus.model import ConfigModel
from prestus_nexus.predict_outputs import predict, save_summary


def _base(schema) -> ConfigModel:
    m = ConfigModel(schema)
    m.set("subject_id", 3)
    m.set("path.sim", "/tmp/prestus_test_helper")
    m.set("simulation.medium", "layered")
    m.set("simulation.transducer_coupling", "coherent")
    m.set("simulation.uncertainty", False)
    m.set("simulation.debug", 0)
    m.set("transducer.target_isppa_wcm2", 30)
    return m


def test_single_mode(schema):
    plan = predict(_base(schema))
    assert plan.mode == "single"
    assert plan.nested == ""
    names = [a.name for a in plan.artifacts]
    assert names == ["sub-003_layered.csv", "sub-003_layered_report.html",
                      "sub-003_layered_targeting.csv"]
    assert len(plan.jobs) == 1


def test_single_mode_water_has_no_targeting(schema):
    m = _base(schema)
    m.set("simulation.medium", "water")
    plan = predict(m)
    names = [a.name for a in plan.artifacts]
    assert len(names) == 2
    assert not any("targeting" in n for n in names)


def test_async_mode(schema):
    m = _base(schema)
    m.set("simulation.transducer_coupling", "async")
    m.add_transducer()
    m.set("transducer.trans_pos", [52, -8, 58], 1)       # a second stance
    m.set("modules.run_posthoc_water_sims", 1)
    plan = predict(m)
    assert plan.mode == "async"
    names = [a.name for a in plan.artifacts]
    assert "sub-003_water_trans-01.csv" in names
    assert "sub-003_water_trans-02.csv" in names
    assert "sub-003_layered_desc-asyncJoint.csv" in names
    assert "sub-003_layered_report_desc-asyncJoint.html" in names
    assert any("2 acoustic jobs" in j for j in plan.jobs)
    assert any("2 free-water validation jobs" in j for j in plan.jobs)


def test_async_repeated_stance_is_simulated_once(schema):
    # T3 returns to T1's stance with other timing and target: one field for both
    # (acoustic_store_key.m); other placement modes resolve every entry on its own
    m = _base(schema)
    m.set("simulation.transducer_coupling", "async")
    m.add_transducer()
    m.add_transducer()
    m.set("transducer.trans_pos", [52, -8, 58], 1)
    m.set("transducer.timing_offset_s", 100.0, 2)
    m.set("transducer.target_isppa_wcm2", 20.0, 2)
    m.set("modules.run_posthoc_water_sims", 1)
    plan = predict(m)
    names = [a.name for a in plan.artifacts]
    assert "sub-003_water_trans-03.csv" not in names
    assert any(j.startswith("2 acoustic jobs — one per distinct field (3 transducers") for j in plan.jobs)
    assert any("2 free-water validation jobs" in j for j in plan.jobs)
    m.set("placement.mode", "heuristic")                 # negative control
    assert any(j.startswith("3 acoustic jobs — one per transducer") for j in predict(m).jobs)


def test_async_sweep(schema):
    m = _base(schema)
    m.set("simulation.transducer_coupling", "async")
    m.add_transducer()
    m.set("transducer.target_isppa_wcm2", [10, 20], 0)
    m.set("transducer.target_isppa_wcm2", [10, 20], 1)
    m.set("modules.run_posthoc_water_sims", 0)
    plan = predict(m)
    names = [a.name for a in plan.artifacts]
    assert len(names) == 4   # 2 sweep points x (csv + html)
    assert "sub-003_layered_desc-asyncJointS01.csv" in names
    assert "sub-003_layered_desc-asyncJointS02.csv" in names


def test_multi_isppa_mode(schema):
    m = _base(schema)
    m.set("transducer.target_isppa_wcm2", [10, 20, 30])
    plan = predict(m)
    assert plan.mode == "multi_isppa"
    names = [a.name for a in plan.artifacts]
    assert len(names) == 7   # 3 targets x (csv + html) + summary
    assert "sub-003_layered_desc-isppa010.csv" in names
    assert "sub-003_layered_desc-multiIsppa_report.html" in names
    assert any("3 thermal jobs" in j for j in plan.jobs)


def test_uncertainty_mode(schema):
    m = _base(schema)
    m.set("simulation.uncertainty", True)
    plan = predict(m)
    assert plan.mode == "uncertainty"
    names = [a.name for a in plan.artifacts]
    assert len(names) == 11  # 5 variants x (csv + html) + combined report
    assert "sub-003_layered_desc-temphigh.csv" in names
    assert "sub-003_layered_desc-translow.csv" in names
    assert any("5 variant simulation jobs" in j and "(stages 2-6)" in j for j in plan.jobs)
    assert any("1 uncertainty report job (stage 7)" in j for j in plan.jobs)
    assert "sub-003_layered_desc-uncertainty_report.html" in names


def test_uncertainty_is_outer_mode_on_every_platform(schema):
    # Uncertainty nests an async run inside each variant (uncertainty_pipeline.m),
    # with the same precedence on every platform (prestus_pipeline_start.m).
    m = _base(schema)
    m.set("simulation.transducer_coupling", "async")
    m.add_transducer()
    m.set("simulation.uncertainty", True)
    for platform in ("matlab", "slurm", "auto"):
        m.set("platform", platform)
        plan = predict(m)
        assert plan.mode == "uncertainty"
        assert plan.nested == "async"


def test_uncertainty_nests_multi_isppa_sweep(schema):
    m = _base(schema)
    m.set("simulation.uncertainty", True)
    m.set("transducer.target_isppa_wcm2", [10, 30])
    plan = predict(m)
    assert plan.mode == "uncertainty"
    assert plan.nested == "multi_isppa"
    names = [a.name for a in plan.artifacts]
    # 5 variants x 2 targets x (csv + html) + 5 summaries + 2 uncertainty reports
    assert len(names) == 27
    assert "sub-003_layered_desc-temphigh_desc-isppa030.csv" in names
    assert "sub-003_layered_desc-multiIsppa_report_desc-transhigh.html" in names
    assert "sub-003_layered_desc-uncertainty_report_desc-isppa010.html" in names
    assert any("2 uncertainty report job(s)" in j for j in plan.jobs)


def test_uncertainty_nests_async_sweep(schema):
    m = _base(schema)
    m.set("simulation.transducer_coupling", "async")
    m.add_transducer()
    m.set("transducer.target_isppa_wcm2", [10, 20], 0)
    m.set("transducer.target_isppa_wcm2", [10, 20], 1)
    m.set("simulation.uncertainty", True)
    plan = predict(m)
    assert plan.nested == "async"
    names = [a.name for a in plan.artifacts]
    assert "sub-003_layered_desc-templow_desc-asyncJointS02.csv" in names
    assert "sub-003_layered_desc-uncertainty_report_desc-asyncJointS01.html" in names


def test_async_multi_isppa_sweep_is_async(schema):
    # A vector target under async coupling is async's own internal sweep
    # (see async_transducer_pipeline.m), not a conflict with multi_isppa_pipeline.
    m = _base(schema)
    m.set("simulation.transducer_coupling", "async")
    m.add_transducer()
    m.set("transducer.target_isppa_wcm2", [10, 20], 0)
    m.set("transducer.target_isppa_wcm2", [10, 20], 1)
    plan = predict(m)
    assert plan.mode == "async"
    assert plan.nested == ""


def test_common_dirs_include_debug_when_enabled(schema):
    m = _base(schema)
    plan = predict(m)
    assert "debug/" not in plan.common_dirs
    m.set("simulation.debug", 1)
    plan = predict(m)
    assert "debug/" in plan.common_dirs


def test_save_summary_reflects_flags(schema):
    m = _base(schema)
    lines = save_summary(m)
    text = "\n".join(lines)
    assert "Global default (io.save_matrices): unset" in text
    # config_default.yaml leaves the matrix flags empty: built-in defaults
    # save_source_matrices=1, save_acoustic_matrices=0 apply
    assert "✓ source matrices" in text
    assert "✗ acoustic matrices" in text
    assert "Debug diagnostics (simulation.debug): OFF" in text

    # empty per-step flags inherit io.save_matrices
    m.set("io.save_matrices", 1)
    assert "✓ acoustic matrices" in "\n".join(save_summary(m))
    m.set("io.save_matrices", 0)
    assert "✗ source matrices" in "\n".join(save_summary(m))

    m.set("io.save_acoustic_matrices", 1)
    m.set("simulation.debug", 1)
    text = "\n".join(save_summary(m))
    assert "✓ acoustic matrices" in text
    assert "Debug diagnostics (simulation.debug): ON" in text
