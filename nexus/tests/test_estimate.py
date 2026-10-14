"""estimate.py: resource estimate from <path.sim>/run_stats.jsonl
(written by functions/helper/run_stats.m)."""

import json

from prestus_nexus.estimate import estimate, load_records, main, stats_path, thermal_steps
from prestus_nexus.model import ConfigModel


def _model(schema, tmp_path) -> ConfigModel:
    m = ConfigModel(schema)
    m.set("path.sim", str(tmp_path))
    m.set("simulation.medium", "layered")
    m.set("simulation.code_type", "cpp_gpu")
    m.set("grid.resolution_mm", 0.5)
    m.set("hpc.partition", "gpu")
    m.set("modules.run_acoustic_sims", 1)
    m.set("modules.run_heating_sims", 0)
    return m


def _rec(**kw) -> dict:
    rec = {"code_type": "cpp_gpu", "acoustic_backend": "kwave", "medium": "layered",
           "resolution_mm": 0.5, "hpc_partition": "gpu", "acoustic_s": 1200.0, "peak_ram_gb": 20.0}
    rec.update(kw)
    return rec


def _write(tmp_path, recs) -> None:
    (tmp_path / "run_stats.jsonl").write_text("\n".join(json.dumps(r) for r in recs) + "\n")


def test_no_records(schema, tmp_path):
    est = estimate(_model(schema, tmp_path))
    assert est.n_records == 0
    assert "No earlier runs" in est.lines()[0]


def test_acoustic_median_and_max(schema, tmp_path):
    _write(tmp_path, [_rec(acoustic_s=600), _rec(acoustic_s=1200), _rec(acoustic_s=3000),
                      _rec(acoustic_s=9, medium="water")])         # different medium: ignored
    est = estimate(_model(schema, tmp_path))
    assert est.acoustic_s.median == 1200 and est.acoustic_s.maximum == 3000
    assert est.acoustic_s.n == 3
    assert est.peak_ram_gb.median == 20
    assert est.thermal_s is None                                   # heating off


def test_partition_falls_back(schema, tmp_path):
    _write(tmp_path, [_rec(hpc_partition="gpu40g", acoustic_s=900)])
    est = estimate(_model(schema, tmp_path))
    assert est.acoustic_s.median == 900
    assert "any partition" in est.lines()[0]


def test_unmatched_resolution(schema, tmp_path):
    _write(tmp_path, [_rec(resolution_mm=1.0)])
    est = estimate(_model(schema, tmp_path))
    assert est.acoustic_s is None
    assert any("no earlier run" in s for s in est.lines())


def test_thermal_scales_with_steps(schema, tmp_path):
    m = _model(schema, tmp_path)
    m.set("modules.run_heating_sims", 1)
    for k, v in {"pd": 0.02, "pri": 0.2, "ptd": 20, "ptri": 40, "ptrd": 80,
                 "post_ptri_dur": 10, "post_pt_timestep": 1, "pt_timestep": 0.02,
                 "equal_step_duration": 0}.items():
        m.set(f"timing.{k}", v)
    # 2 reps x (100 pulses x 2 + 20 off steps) + 10 post steps
    assert thermal_steps(m) == 450
    _write(tmp_path, [_rec(thermal_s=100.0, thermal_steps=1000)])  # 0.1 s/step
    est = estimate(m)
    assert abs(est.thermal_s.median - 45.0) < 1e-9


def test_thermal_equal_steps(schema, tmp_path):
    m = _model(schema, tmp_path)
    for k, v in {"pd": 0.02, "pri": 0.1, "ptd": 1, "ptri": 1, "ptrd": 2,
                 "post_ptri_dur": 0, "post_pt_timestep": 1, "pt_timestep": 0.02,
                 "equal_step_duration": 1}.items():
        m.set(f"timing.{k}", v)
    # 2 reps x 10 pulses x (1 on + 4 off) steps, no PTRI-off or post steps
    assert thermal_steps(m) == 100


def test_bad_lines_skipped(tmp_path):
    p = tmp_path / "run_stats.jsonl"
    p.write_text('{"acoustic_s": 1}\nnot json\n[1, 2]\n{"acoustic_s": NaN}\n')
    assert len(load_records(p)) == 2       # Python's json accepts NaN; _num drops it later


def test_relative_path_sim(schema, tmp_path):
    m = ConfigModel(schema)
    m.set("path.anat", str(tmp_path))
    m.set("path.sim", "sims")
    assert stats_path(m) == tmp_path / "sims" / "run_stats.jsonl"


def test_headless(schema, tmp_path, capsys, prestus_root):
    _write(tmp_path, [_rec()])
    study = tmp_path / "study.yaml"
    _model(schema, tmp_path).save(study)
    assert main([str(study), "--prestus-root", str(prestus_root)]) == 0
    assert "Acoustic sim" in capsys.readouterr().out
