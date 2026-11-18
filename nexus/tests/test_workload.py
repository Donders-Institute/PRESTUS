"""workload.py: config -> grid, time steps, memory and time, before any run."""

import json

from prestus_nexus import workload as w
from prestus_nexus.model import ConfigModel


def _model(schema, tmp_path) -> ConfigModel:
    m = ConfigModel(schema)
    m.set("path.sim", str(tmp_path))
    m.set("simulation.medium", "phantom")
    m.set("simulation.code_type", "cpp_gpu")
    m.set("grid.resolution_mm", 0.5)
    m.set("grid.pml_size", 40)
    m.set("grid.axisymmetric", 1)
    m.set("grid.default_dims", [140, 240])
    m.set("transducer.freq_hz", 500000)
    m.set("layers", {"water": [0], "brain": [1], "skull": [4]})
    m.set("medium_properties.water.sound_speed", 1482)
    m.set("medium_properties.brain.sound_speed", 1560)
    m.set("medium_properties.skull.sound_speed", 2550)
    return m


def test_nt_matches_itrusst_log(schema, tmp_path):
    # PRESTUS_itrusst sub-361 log: grid [240 71], dt 2.9412e-08 s, Nt 2871
    m = _model(schema, tmp_path)
    assert w.sound_speeds(m) == (1482.0, 2550.0)       # unused layers (skull_cortical) ignored
    grid = w.acoustic_grid(m, [])
    assert w.simulated_dims(m, grid) == (240, 71)
    assert w.acoustic_steps(m, (240, 71)) == 2871


def test_optimal_pml_minimises_prime_factor():
    p = w.optimal_pml(240)
    assert 10 <= p <= 40
    assert w._max_prime(240 + 2 * p) == min(w._max_prime(240 + 2 * q) for q in range(10, 41))


def test_water_dims_follow_water_grid_dims():
    # aperture 64 mm, apex->focus 60 mm, 0.5 mm: lateral 2*ceil(74/0.5/2)+1, axial ceil(100/0.5)+2
    class M:
        def get(self, k, ti=0):
            return {"grid.water_margin_mm": None, "grid.z_beyond_focus_mm": None}.get(k)
    w_ap = w._aperture_and_bowl_depth
    try:
        w._aperture_and_bowl_depth = lambda model: (64.0, 60.0)
        assert w.water_dims(M(), 0.5) == (149, 149, 202)
    finally:
        w._aperture_and_bowl_depth = w_ap


def test_table_rate_matches_gpu_or_spans_kind(schema, tmp_path):
    table = [{"device": "A", "kind": "gpu", "match": ["l40s"], "ns": 16, "range": [15, 17], "n": 4},
             {"device": "B", "kind": "gpu", "match": ["p100"], "ns": 36, "range": [35, 50], "n": 4},
             {"device": "C", "kind": "cpu", "match": [], "ns": 700, "range": [550, 1350], "n": 9}]
    m = _model(schema, tmp_path)
    m.set("hpc.gpu", "nvidia_l40s:1")
    assert w.table_rate(table, m).ns == 16
    m.set("hpc.gpu", "gpu:1")
    r = w.table_rate(table, m)
    assert (r.ns, r.lo, r.hi) == (26, 16, 36)
    m.set("simulation.code_type", "cpp_cpu")
    assert w.table_rate(table, m).ns == 700


def test_project_records_override_table(schema, tmp_path):
    m = _model(schema, tmp_path)
    rec = {"code_type": "cpp_gpu", "thermal_voxels": 1e6, "thermal_steps": 100, "thermal_s": 2.0}
    assert w._rate([rec], m, "thermal", {}).ns == 20.0


def test_thermal_part_uses_expanded_grid(schema, tmp_path):
    m = _model(schema, tmp_path)
    m.set("modules.run_heating_sims", 1)
    for k, v in {"pd": 0.02, "pri": 0.2, "ptd": 1, "pt_timestep": 0.01, "ptri": 1, "ptrd": 1,
                 "post_ptri_dur": 0, "post_pt_timestep": 1, "equal_step_duration": 0}.items():
        m.set(f"timing.{k}", v)
    pred = w.predict(m, records=[], throughput={"thermal": [
        {"device": "G", "kind": "gpu", "match": [], "ns": 20, "range": [20, 20], "n": 1}]})
    th = [p for p in pred.parts if p.label == "Thermal sim"]
    assert th and th[0].points == 140 * 140 * 240 and th[0].steps == 10
    assert "(" not in pred.lines()[2].split(";")[0]          # no range when low == high


def _rec(dev="NVIDIA L40S", code="cpp_gpu", n=100, nt=1000, t0=4.0, ns=0.6, **kw):
    """A run_stats line whose acoustic time follows t0 + ns x points x steps exactly."""
    s = t0 + ns * 1e-9 * n ** 3 * nt
    return {"code_type": code, "gpu_model": dev, "acoustic_backend": "kwave", "precision": "single",
            "acoustic_dims": [n, n, n], "acoustic_pml": 0, "acoustic_nt": nt, "acoustic_s": s, **kw}


def test_fit_recovers_t0_and_rate():
    recs = [_rec(n=n, nt=nt) for n in (80, 120, 200, 260) for nt in (2000, 5000)]
    ac = w.fit(recs)["acoustic"][0]
    assert abs(ac["t0_s"] - 4.0) < 0.15 and abs(ac["ns"] - 0.6) < 0.01
    assert (ac["code_type"], ac["backend"], ac["axisymmetric"], ac["precision"]) == ("cpp_gpu", "kwave", 0, "single")
    one = w.fit([_rec(n=200, nt=2000)])["acoustic"][0]          # too few records: pooled t0
    assert one["t0_s"] == w.DEFAULT_T0_S


def test_fit_thermal_and_roundtrip():
    rec = {"code_type": "cpp_gpu", "gpu_model": "NVIDIA L40S", "thermal_voxels": 1e6, "thermal_steps": 100,
           "thermal_s": 1.6}
    out = w.fit([rec])
    assert out["thermal"][0]["ns"] == 16.0 and out["thermal"][0]["match"] == ["l40s"]
    json.dumps(out)


def test_fit_drops_outliers():
    recs = [_rec(n=n, nt=3000) for n in (80, 100, 120, 160, 200, 240, 260)]
    recs.append(dict(recs[3], acoustic_s=recs[3]["acoustic_s"] * 6))    # e.g. a sum over transducers
    ac = w.fit(recs)["acoustic"][0]
    assert ac["n"] == 7 and abs(ac["ns"] - 0.6) < 0.01


def test_telemetry_events_fit_like_run_stats():
    line = _rec(n=150, nt=3000)
    sim = {k: line[k] for k in ("acoustic_dims", "acoustic_pml", "acoustic_nt", "acoustic_s")}
    event = {"event": "run_end", "status": "success", "code_type": "cpp_gpu", "gpu_model": "NVIDIA L40S",
             "precision": "single", "acoustic_backend": "kwave", "sim_stats": sim, "duration_s": 99}
    wrapped = {"id": 7, "payload": event}
    other = dict(event, event="run_start")
    failed = dict(event, status="error")
    assert w.fit([line])["acoustic"] == w.fit([event])["acoustic"] == w.fit([wrapped, other, failed])["acoustic"]
    assert w.flatten(event)["total_s"] == 99 and w.flatten(other) is None and w.flatten(failed) is None


def test_acoustic_work_exact_approx_and_auto():
    exact = w.acoustic_work({"acoustic_dims": [100, 80], "acoustic_pml": 10, "acoustic_nt": 10, "acoustic_s": 1})
    assert exact == (120 * 100 * 10, False)
    old = w.acoustic_work({"acoustic_voxels": 1e6, "acoustic_nt": 10, "acoustic_s": 1})
    assert old[1] and old[0] == (100 + 2 * w.optimal_pml(100)) ** 3 * 10
    assert w.acoustic_work({"acoustic_voxels": 1e6, "acoustic_nt": "auto", "acoustic_s": 1}) is None


def test_factors_from_paired_solvers(schema, tmp_path):
    t0 = w.DEFAULT_T0_S                                        # two records per solver: pooled t0
    recs = [_rec(n=n, nt=3000, t0=t0) for n in (100, 200)]
    recs += [_rec(code="matlab_gpu", n=n, nt=3000, ns=1.8, t0=t0) for n in (100, 200)]
    recs += [_rec(n=n, nt=3000, ns=0.96, t0=t0, precision="double") for n in (100, 200)]
    table = w.fit(recs)
    assert table["factors"]["code_type"]["matlab_gpu"]["n"] == 1
    assert abs(table["factors"]["code_type"]["matlab_gpu"]["x"] - 3.0) < 0.1
    assert abs(table["factors"]["precision"]["double"]["x"] - 1.6) < 0.1
    # a solver without an entry on a device: the reference entry x the factor
    only_cpp = {"acoustic": [e for e in table["acoustic"] if e["code_type"] == "cpp_gpu"
                             and e["precision"] == "single"], "factors": table["factors"]}
    m = _model(schema, tmp_path)
    m.set("grid.axisymmetric", 0)
    m.set("hpc.gpu", "gpu:nvidia_l40s:1")
    m.set("simulation.code_type", "matlab_gpu")
    r = w.table_rate(only_cpp["acoustic"], m, only_cpp)
    assert abs(r.ns - 1.8) < 0.1 and "matlab_gpu" in r.basis


def test_gpu_matching_longest_key_and_bandwidth(schema, tmp_path):
    full = {"device": "NVIDIA A100 80GB PCIe", "kind": "gpu", "match": w._aliases("NVIDIA A100 80GB PCIe"),
            "ns": 0.6, "t0_s": 5, "n": 7}
    mig = {"device": "NVIDIA A100 80GB PCIe MIG 3g.40gb", "kind": "gpu",
           "match": w._aliases("NVIDIA A100 80GB PCIe MIG 3g.40gb"), "ns": 0.9, "t0_s": 5, "n": 9}
    l40s = {"device": "NVIDIA L40S", "kind": "gpu", "match": ["l40s"], "ns": 0.8, "t0_s": 1, "n": 4}
    bw = [{"match": ["a100-80gb-pcie"], "gbs": 1935, "arch": "ampere"},
          {"match": ["a100-80gb-pcie-3g"], "gbs": 968, "arch": "ampere"},
          {"match": ["rtx-a6000"], "gbs": 768, "arch": "ampere"},
          {"match": ["l40s"], "gbs": 864, "arch": "ada"},
          {"match": ["h100-pcie"], "gbs": 2000, "arch": "hopper"}]
    m = _model(schema, tmp_path)
    m.set("grid.axisymmetric", 0)
    for gres, want in (("gpu:nvidia_a100_80gb_pcie_3g.39gb:1", 0.9), ("gpu:nvidia_a100_80gb_pcie:1", 0.6)):
        m.set("hpc.gpu", gres)
        assert w.table_rate([full, mig, l40s], m, {"bandwidth": bw}).ns == want
    m.set("hpc.gpu", "gpu:nvidia_rtx_a6000:1")                 # unmeasured, same generation as the A100s
    r = w.table_rate([full, mig, l40s], m, {"bandwidth": bw})
    k = sorted([0.6 * 1935, 0.9 * 968])
    assert abs(r.ns - (k[0] + k[1]) / 2 / 768) < 1e-9 and "bandwidth" in r.basis and r.lo < r.ns < r.hi
    m.set("hpc.gpu", "gpu:nvidia_h100_pcie:1")                 # no measured generation: span of all
    assert w.table_rate([full, mig, l40s], m, {"bandwidth": bw}).basis.startswith("any GPU")


def test_validate_leave_one_study_out():
    recs = [_rec(n=n, nt=nt, study=s) for s in ("a", "b", "c") for n in (80, 160, 240) for nt in (2000, 6000)]
    v = w.validate(recs)
    assert v["n"] == len(recs) and abs(v["median"] - 1) < 0.01


def test_stage_fit_scales_with_resolution():
    recs = [{"medium": "layered", "resolution_mm": dx, "stage_preproc_s": 60 / dx ** 3, "stage_preproc_skullsmooth_s": 40 / dx ** 3,
             "stage_placement_s": 30,
             "stage_segmentation_s": 0.2} for dx in (0.5, 0.75, 1.0, 1.0, 1.25, 1.5)]
    st = {e["stage"]: e for e in w.fit(recs)["stages"]}
    assert st["preproc"]["medium"] == "head" and abs(st["preproc"]["exponent"] - 1) < 0.01
    assert abs(st["preproc"]["s_1mm"] - 60) < 0.5 and st["placement"]["exponent"] == 0
    assert "segmentation" not in st                            # under a minute: charm output reused


def _heating(m):
    m.set("modules.run_heating_sims", 1)
    for k, v in {"pd": 0.02, "pri": 0.2, "ptd": 1, "pt_timestep": 0.01, "ptri": 1, "ptrd": 1,
                 "post_ptri_dur": 0, "post_pt_timestep": 1, "equal_step_duration": 0}.items():
        m.set(f"timing.{k}", v)


def test_predict_covers_every_running_stage(schema, tmp_path):
    from prestus_nexus import run_plan as rp
    m = _model(schema, tmp_path)
    _heating(m)
    labels = {lb for v in w.STAGE_TIMERS.values() for lb in v} | {"nifti_thermal"}
    table = {"stages": [{"stage": lb, "medium": "phantom", "s_1mm": 10, "exponent": 0, "ratio": [0.5, 2], "n": 3}
                        for lb in labels],
             "thermal": [{"device": "G", "kind": "gpu", "match": [], "ns": 20, "range": [20, 20], "n": 1}],
             "acoustic": [{"device": "G", "kind": "gpu", "match": [], "ns": 1.0, "t0_s": 5, "range": [1, 1], "n": 1}]}
    pred = w.predict(m, records=[], throughput=table)
    running = {s.name for s in rp.plan(m).stages if s.state == rp.RUN}
    covered = {p.label for p in pred.parts} | {s.label for s in pred.stages}
    sims = {"Acoustic simulation": "Acoustic sim", "Thermal simulation": "Thermal sim"}
    for name in running:
        assert sims.get(name, name) in covered or any(name.lower() in u or (w.UNTIMED.get(name) or "") == u
                                                      for u in pred.untimed) or name in w.UNTIMED, name
    nifti = next(s for s in pred.stages if s.label == "NIfTI maps")
    assert nifti.s == 20                                       # acoustic + thermal maps
    t = pred.total()
    assert t is not None and t[0] == sum(p.seconds()[0] for p in pred.parts) + sum(s.s for s in pred.stages)
    line = pred.summary()
    assert line.startswith("~") and "other stages" not in line        # total only; breakdown in lines()
    assert any(x.startswith("Other stages:") for x in pred.lines())
    # without stage entries the stages are listed as not timed
    bare = w.predict(m, records=[], throughput={k: v for k, v in table.items() if k != "stages"})
    assert not bare.stages and bare.untimed and "Not timed" in "\n".join(bare.lines())


def test_project_stage_times_override_table(schema, tmp_path):
    m = _model(schema, tmp_path)
    recs = [{"medium": "phantom", "resolution_mm": 0.5, "stage_source_s": s} for s in (7, 8, 9)]
    est = w._stage_estimate("source", m, recs, {"stages": [
        {"stage": "source", "medium": "phantom", "s_1mm": 100, "exponent": 0, "ratio": [1, 1], "n": 1}]})
    assert est[0] == 8 and "project" in est[3]


def test_baseline_box(schema, tmp_path):
    m = _model(schema, tmp_path)
    m.set("grid.axisymmetric", 0)
    w_ap = w._aperture_and_bowl_depth
    try:
        w._aperture_and_bowl_depth = lambda model: (64.0, 100.0)
        # 1 mm: lateral ceil((64 + 20) / 8) * 8 = 88, axial ceil((140 + 20) / 8) * 8 = 160
        m.set("grid.resolution_mm", 1.0)
        assert w.baseline_dims(m, (256, 216, 192)) == (88, 88, 160)
        assert w.baseline_dims(m, (90, 90, 170)) == (90, 90, 170)      # box over 80%: full grid
        m.set("grid.axisymmetric", 1)
        assert w.baseline_dims(m, (240, 71)) == (240, 71)
    finally:
        w._aperture_and_bowl_depth = w_ap


def test_merge_replaces_by_key_and_keeps_the_rest():
    old = {"bandwidth": [{"match": ["l40s"], "gbs": 864, "arch": "ada"}],
           "thermal": [{"device": "x", "kind": "cpu", "ns": 700}, {"device": "y", "kind": "gpu", "ns": 20}],
           "factors": {"code_type": {"matlab_cpu": {"x": 4.7, "n": 1}}}}
    new = {"acoustic": [{"device": "y", "kind": "gpu", "ns": 1}], "thermal": [{"device": "y", "kind": "gpu", "ns": 16}],
           "stages": [], "factors": {"code_type": {"matlab_gpu": {"x": 3.5, "n": 5}}, "precision": {}}}
    out = w.merge(old, new)
    assert out["bandwidth"] == old["bandwidth"] and out["acoustic"] == new["acoustic"]
    assert {(e["device"], e["ns"]) for e in out["thermal"]} == {("x", 700), ("y", 16)}
    assert set(out["factors"]["code_type"]) == {"matlab_cpu", "matlab_gpu"}


def test_stage_fit_skips_thermal_outputs_without_heating():
    recs = [{"medium": "layered", "resolution_mm": 1.0, "stage_thermal_analysis_s": 0.1, "stage_thermal_s": 0.1},
            {"medium": "layered", "resolution_mm": 1.0, "stage_thermal_analysis_s": 12, "thermal_s": 40}]
    st = [e for e in w.fit(recs)["stages"] if e["stage"] == "thermal_analysis"]
    assert st and st[0]["n"] == 1 and st[0]["s_1mm"] == 12


def test_fitted_entry_matches_gres_and_solver(schema, tmp_path):
    entry = {"device": "NVIDIA A100-PCIE-40GB", "kind": "gpu", "match": ["a100-pcie-40gb"],
             "code_type": "cpp_gpu", "backend": "kwave", "axisymmetric": 0, "ns": 3.0, "n": 2}
    other = dict(entry, device="other", match=["l40s"], ns=1.0)
    m = _model(schema, tmp_path)
    m.set("grid.axisymmetric", 0)
    m.set("hpc.gpu", "gpu:nvidia_a100-pcie-40gb:1")          # SLURM gres spelling
    assert w.table_rate([other, entry], m).ns == 3.0
    m.set("simulation.code_type", "matlab_gpu")               # entries hold for cpp_gpu only
    assert w.table_rate([other, entry], m) is None


def test_shipped_table_loads():
    t = w.load_throughput()
    assert any(e["kind"] == "gpu" for e in t["thermal"])
    assert t["acoustic"] and all(w._num(e.get("ns")) and w._num(e.get("t0_s")) is not None and e.get("n")
                                 for e in t["acoustic"])
    assert all(w._num(b.get("gbs")) and b.get("arch") and b.get("match") for b in t["bandwidth"])
    assert {"preproc", "placement", "acoustic_analysis", "nifti_acoustic"} <= {
        e["stage"] for e in t["stages"] if e["medium"] == "head"}


def test_summary_line(schema, tmp_path):
    m = _model(schema, tmp_path)
    m.set("modules.run_heating_sims", 1)
    for k, v in {"pd": 0.02, "pri": 0.2, "ptd": 1, "pt_timestep": 0.01, "ptri": 1, "ptrd": 1,
                 "post_ptri_dur": 0, "post_pt_timestep": 1, "equal_step_duration": 0}.items():
        m.set(f"timing.{k}", v)
    pred = w.predict(m, records=[], throughput={"thermal": [
        {"device": "G", "kind": "gpu", "match": [], "ns": 20, "range": [20, 20], "n": 1}]})
    line = pred.summary()
    assert line.startswith("acoustic and 2 more: time unknown")   # no total while a simulation is untimed
    assert "thermal" not in line
    assert w.Prediction().summary().startswith("Set the resolution")


def test_preproc_paths_fitted_and_chosen_apart(schema, tmp_path):
    shared = {"medium": "layered", "resolution_mm": 0.75, "stage_preproc_s": 550, "stage_preproc_placement_s": 520}
    classic = {"medium": "layered", "resolution_mm": 0.75, "stage_preproc_s": 215, "stage_preproc_skullsmooth_s": 150}
    cached = {"medium": "layered", "resolution_mm": 0.75, "stage_preproc_s": 2}
    st = {e["stage"]: e for e in w.fit([shared, classic, cached])["stages"]}
    assert st["preproc_shared"]["s_1mm"] == 550 and st["preproc"]["s_1mm"] == 215
    assert st["preproc"]["n"] == 1 and st["preproc_cached"]["s_1mm"] == 2      # cache loads kept apart
    table = {"stages": [st["preproc"], st["preproc_shared"]]}
    m = _model(schema, tmp_path)
    m.set("simulation.medium", "layered")
    m.set("grid.resolution_mm", 1.0)
    for reuse, want in ((False, 215), (True, 550)):
        m.set("headmodel.reuse_across_placements", reuse)
        assert w._stage_estimate(w._preproc_label(m), m, [], table)[0] == want


def test_doc_calculator_data_matches_table_rate(schema, tmp_path):
    """The docs calculator (tests/docs/gen_load_data.py) resolves rates as table_rate does."""
    import importlib.util
    from pathlib import Path
    gen_path = Path(__file__).resolve().parents[2] / "tests" / "docs" / "gen_load_data.py"
    spec = importlib.util.spec_from_file_location("gen_load_data", gen_path)
    gen = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gen)
    table = w.load_throughput()
    devices = gen.devices(table)
    assert any(not d["measured"] for d in devices) and any(d["kind"] == "cpu" for d in devices)
    m = _model(schema, tmp_path)
    m.set("grid.axisymmetric", 0)
    m.set("acoustic_backend", "kwave")
    for d in devices:
        m.set("hpc.gpu", f"gpu:{d['hint']}:1" if d["kind"] == "gpu" else "")
        for key, r in d["rates"].items():
            solver, precision = key.split("|")
            m.set("simulation.code_type", solver)
            m.set("simulation.precision", precision)
            want = w.table_rate(table["acoustic"], m, table)
            assert want is not None, (d["id"], key)
            assert abs(want.ns - r["ns"]) < 1e-3 and abs(want.t0 - r["t0"]) < 1e-6, (d["id"], key, want, r)
        th = w.table_rate(table["thermal"], m, table) if d["thermal"] else None
        if d["thermal"]:
            assert th is not None and abs(th.ns - d["thermal"]["ns"]) < 1e-3, d["id"]


def test_local_gpu_detected_only_for_local_gpu_runs_without_hpc_gpu(schema, tmp_path, monkeypatch):
    l40s = {"device": "NVIDIA L40S", "kind": "gpu", "match": ["l40s"], "ns": 0.8, "t0_s": 1, "n": 4}
    a100 = {"device": "NVIDIA A100-SXM4-40GB", "kind": "gpu", "match": ["a100-sxm4-40gb"], "ns": 0.55, "t0_s": 5, "n": 9}
    cpu = {"device": "CPU", "kind": "cpu", "match": [], "ns": 90, "t0_s": 5, "n": 1}
    entries = [l40s, a100, cpu]
    monkeypatch.setattr(w, "local_gpu", lambda: "NVIDIA L40S")
    m = _model(schema, tmp_path)
    m.set("grid.axisymmetric", 0)
    m.set("hpc.gpu", "")
    cases = [  # platform, hpc.gpu, code_type -> ns, detected in basis
        ("matlab", "", "cpp_gpu", 0.8, True),
        ("matlab", "gpu:nvidia_a100-sxm4-40gb:1", "cpp_gpu", 0.55, False),   # hpc.gpu wins
        ("slurm", "", "cpp_gpu", None, False),                                # cluster: span
        ("auto", "", "cpp_gpu", None, False),
        ("matlab", "", "cpp_cpu", 90, False),                                 # CPU run: no GPU lookup
    ]
    for platform, gpu, code, ns, detected in cases:
        m.set("platform", platform)
        m.set("hpc.gpu", gpu)
        m.set("simulation.code_type", code)
        r = w.table_rate(entries, m, {})
        assert ("this computer" in r.basis) == detected, (platform, gpu, code, r.basis)
        if ns is None:
            assert r.basis.startswith("any GPU"), r.basis
        else:
            assert r.ns == ns, (platform, gpu, code)
    monkeypatch.setattr(w, "local_gpu", lambda: None)                        # no NVIDIA GPU found
    m.set("platform", "matlab")
    m.set("hpc.gpu", "")
    m.set("simulation.code_type", "cpp_gpu")
    assert w.table_rate(entries, m, {}) is None          # GPU code cannot run here: no cluster timings


def test_local_gpu_reads_nvidia_smi_once(monkeypatch):
    calls = []

    class Done:
        returncode, stdout = 0, "NVIDIA GeForce RTX 4090\nNVIDIA GeForce RTX 4090\n"
    monkeypatch.setattr(w, "_LOCAL_GPU", [])
    monkeypatch.setattr(w.shutil, "which", lambda exe: "/usr/bin/nvidia-smi")
    monkeypatch.setattr(w.subprocess, "run", lambda *a, **k: calls.append(a) or Done())
    assert w.local_gpu() == "NVIDIA GeForce RTX 4090" and w.local_gpu() == "NVIDIA GeForce RTX 4090"
    assert len(calls) == 1
    monkeypatch.setattr(w, "_LOCAL_GPU", [])
    monkeypatch.setattr(w.shutil, "which", lambda exe: None)                # no driver: no GPU
    assert w.local_gpu() is None


def test_local_run_without_gpu(schema, tmp_path, monkeypatch):
    """A GPU run on this computer (platform matlab) with no NVIDIA GPU: no cluster timings,
    a note why, and the footer says so; with hpc.gpu set or on the cluster, timings as before."""
    m = _model(schema, tmp_path)
    m.set("simulation.code_type", "matlab_gpu")
    table = {"acoustic": [{"device": "A100", "kind": "gpu", "match": ["a100"], "ns": 1, "range": [1, 1], "n": 1}]}
    monkeypatch.setattr(w, "local_gpu", lambda: None)
    m.set("platform", "slurm")
    assert w.predict(m, records=[], throughput=table).parts[0].seconds() is not None
    m.set("platform", "matlab")
    pred = w.predict(m, records=[], throughput=table)
    assert pred.parts[0].seconds() is None and any(n.startswith(w.NO_LOCAL_GPU) for n in pred.notes)
    assert pred.summary().startswith("No NVIDIA GPU on this computer")
    monkeypatch.setattr(w, "local_gpu", lambda: "NVIDIA A100-PCIE-40GB")
    assert w.predict(m, records=[], throughput=table).parts[0].seconds() is not None


def test_segmentation_time_excludes_pseudo_ct_and_reuse():
    built = w.stage_seconds({"stage_segmentation_s": 3900, "stage_pct_s": 600})
    assert built["segmentation"] == 3300 and built["pct"] == 600
    reused = w.stage_seconds({"stage_segmentation_s": 620, "stage_pct_s": 600})   # charm reused, pCT built
    assert "segmentation" not in reused and reused["pct"] == 600
    assert "segmentation" not in w.stage_seconds({"stage_segmentation_s": 0.2})


def test_new_stage_timers_are_mapped_and_queue_listed(schema, tmp_path, monkeypatch):
    monkeypatch.setattr(w, "local_gpu", lambda: None)
    for name, label in (("Skull CT", "pct"), ("Report", "report"), ("Group report", "group_report"),
                        ("Segmentation", "segmentation")):
        assert w.STAGE_TIMERS[name] == (label,)
    m = _model(schema, tmp_path)
    for platform in ("slurm", "auto", "matlab"):     # queue wait is not a stage: never listed as untimed
        m.set("platform", platform)
        pred = w.predict(m, records=[], throughput={})
        assert not any("queue" in u for u in pred.untimed), platform
