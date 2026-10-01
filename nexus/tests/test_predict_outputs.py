"""predict_outputs.py mirrors functions/helper/predict_outputs.m (and the mode
predicates is_async_mode.m / is_multi_isppa_mode.m it calls). These cases mirror
tests/test_helper.m's predict_outputs block one for one — if one side changes,
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
    assert not plan.ambiguous
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
    assert len(names) == 7   # 3 variants x (csv + html) + combined report
    assert "sub-003_layered_desc-liberal.csv" in names
    assert "sub-003_layered_uncertainty_report.html" in names


def test_ambiguous_precedence_by_platform(schema):
    m = _base(schema)
    m.set("simulation.transducer_coupling", "async")
    m.add_transducer()
    m.set("simulation.uncertainty", True)

    m.set("platform", "matlab")
    plan_matlab = predict(m)
    assert plan_matlab.mode == "async"
    assert plan_matlab.ambiguous

    m.set("platform", "slurm")
    plan_hpc = predict(m)
    assert plan_hpc.mode == "uncertainty"
    assert plan_hpc.ambiguous


def test_ambiguous_with_platform_auto_is_disclosed_not_guessed(schema):
    # platform='auto' (config_default.yaml's own default) is only resolved to a
    # concrete platform by hpc_detect_system() inside prestus_pipeline_start.m at
    # submission time — a prediction made beforehand cannot know which way it goes.
    m = _base(schema)
    m.set("simulation.transducer_coupling", "async")
    m.add_transducer()
    m.set("simulation.uncertainty", True)
    m.set("platform", "auto")

    plan = predict(m)
    assert plan.platform_used == "auto"
    assert plan.ambiguous
    assert any("platform is 'auto'" in w for w in plan.warnings)


def test_async_multi_isppa_sweep_is_not_ambiguous(schema):
    # A vector target under async coupling is async's own internal sweep
    # (see async_transducer_pipeline.m), not a conflict with multi_isppa_pipeline.
    m = _base(schema)
    m.set("simulation.transducer_coupling", "async")
    m.add_transducer()
    m.set("transducer.target_isppa_wcm2", [10, 20], 0)
    m.set("transducer.target_isppa_wcm2", [10, 20], 1)
    plan = predict(m)
    assert plan.mode == "async"
    assert not plan.ambiguous


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
