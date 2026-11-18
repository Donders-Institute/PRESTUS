"""Sonication sequences (sequence.py): steps read from the async 'custom' keys, the backend
choice for every case, the edit operations, and the chain of step configs."""

from __future__ import annotations

from pathlib import Path

import pytest

from prestus_nexus import sequence as sq
from prestus_nexus import timing as tm
from prestus_nexus.launcher import PipelineOptions
from prestus_nexus.model import ConfigModel


def _model(schema, n_steps: int = 2, pause: float = 5.0) -> ConfigModel:
    """Two stances: 10 s each (5 trains of 1 s every 2 s), a pause between them."""
    m = ConfigModel(schema)
    for key, v in {"timing.pd": 0.02, "timing.pri": 0.1, "timing.ptd": 1.0, "timing.ptri": 2.0,
                   "timing.ptrd": 10.0, "timing.post_pt_timestep": 1.0,
                   "timing.post_ptri_dur": 0.0, "modules.run_heating_sims": 1}.items():
        m.set(key, v)
    m.set("transducer.freq_hz", 5e5)
    m.set("transducer.trans_pos", [10, 20, 30])
    m.set("transducer.focus_pos", [10, 60, 30])
    for _ in range(n_steps - 1):
        ti = sq.add_step(m, pause)
        m.set("transducer.trans_pos", [40 + ti, 20, 30], ti)
    return m


# ── reading ───────────────────────────────────────────────────────────────────

def test_not_a_sequence(schema):
    m = ConfigModel(schema)
    assert not sq.is_sequence(m) and sq.plan(m) is None and sq.steps(m) == []
    m.add_transducer()
    for coupling, schedule in (("coherent", "custom"), ("async", "pulses"), ("async", "trains")):
        m.set("simulation.transducer_coupling", coupling)
        m.set("timing.async_schedule", schedule)
        assert sq.plan(m) is None


def test_add_step_lays_out_back_to_back(schema):
    m = _model(schema, 3, pause=5.0)
    assert m.get("simulation.transducer_coupling") == "async"
    assert m.get("timing.async_schedule") == "custom"
    st = sq.steps(m)
    assert [s.transducers for s in st] == [[0], [1], [2]]
    assert [s.start_s for s in st] == pytest.approx([0, 15, 30])
    assert sq.pauses(st) == pytest.approx([5, 5])
    # identical timing: nothing but the start delay is written per transducer
    assert "timing" not in m.tx[1] and m.tx[1]["timing_offset_s"] == pytest.approx(15)


def test_step_timing_goes_to_overrides_and_later_steps_move(schema):
    m = _model(schema, 3)
    sq.set_timing(m, 1, {"ptrd": 4.0})
    assert m.tx[1]["timing"] == {"ptrd": 4.0}
    assert m.get("transducer.timing.ptrd", 2) in ([], None) or "timing" not in m.tx[2]
    st = sq.steps(m)
    assert st[1].duration_s == pytest.approx(4)
    assert st[2].start_s == pytest.approx(15 + 4 + 5)             # pause kept
    # step 1's timing is the timing block; the other steps keep their own values
    sq.set_timing(m, 0, {"ptrd": 6.0})
    assert m.get("timing.ptrd") == pytest.approx(6) and "timing" not in m.tx[0]
    assert m.tx[2]["timing"]["ptrd"] == pytest.approx(10)
    assert [s.duration_s for s in sq.steps(m)] == pytest.approx([6, 4, 10])


def test_schedule_matches_steps(schema):
    """The async schedule the timing diagram draws has the same starts and durations."""
    m = _model(schema, 3)
    sq.set_timing(m, 1, {"ptrd": 4.0})
    sched = tm.AsyncSchedule.from_model(m)
    for s in sq.steps(m):
        src = sched.sources[s.transducers[0]]
        assert src.offset == pytest.approx(s.start_s)
        assert src.protocol.trains * src.protocol.ptri == pytest.approx(s.duration_s)


# ── backend choice: every case ────────────────────────────────────────────────

def test_single_stances_run_jointly(schema):
    pl = sq.plan(_model(schema))
    assert pl.backend == "joint" and not pl.errors


def test_chain_when_chosen(schema):
    pl = sq.plan(_model(schema), "chain")
    assert pl.backend == "chain" and not pl.errors


def test_together_runs_jointly_on_the_shared_grid(schema):
    # a step firing together is a firing group of one joint run (firing_groups.m)
    m = _model(schema, 3)
    sq.join_previous(m, 2)
    st = sq.steps(m)
    assert [s.transducers for s in st] == [[0], [1, 2]]
    assert sq.plan(m).backend == "joint" and not sq.plan(m).errors
    assert sq.plan(m, "chain").backend == "chain" and not sq.plan(m, "chain").errors
    # uncertainty mode wraps the joint run (the chain still refuses it)
    m.set("simulation.uncertainty", True)
    assert not sq.plan(m).errors
    assert any(i.path == "simulation.uncertainty" for i in sq.plan(m, "chain").errors)
    # the group's overlap is not the timing panel's same-frequency error
    sched = tm.AsyncSchedule.from_model(m)
    assert sched.groups() == [0, 1, 1]
    assert not any("same frequency" in p for p in sched.problems())


def test_together_on_per_transducer_grids_needs_a_chain(schema):
    # negative control: complex fields are not interpolated between grids
    m = _model(schema, 3)
    sq.join_previous(m, 2)
    m.set("grid.async_layout", "per_transducer")
    assert sq.plan(m).backend == "chain" and not sq.plan(m).errors
    joint = sq.plan(m, "joint")
    assert joint.backend == "" and any(i.path == "grid.async_layout" for i in joint.errors)


def test_overlap_without_shared_timing_is_still_an_error(schema):
    # negative control: same start, different timing is no group; overlapping at one
    # frequency stays the timing panel's error (thermal_schedule.m, coherentOverlap)
    m = _model(schema, 3)
    sq.join_previous(m, 2)
    m.set("transducer.timing.pd", 0.04, 2)
    sched = tm.AsyncSchedule.from_model(m)
    assert sched.groups() == [0, 1, 2]
    assert any("same frequency" in p for p in sched.problems())


def test_separate_undoes_join(schema):
    m = _model(schema, 3)
    sq.join_previous(m, 2)
    sq.separate(m, 1)
    assert [s.transducers for s in sq.steps(m)] == [[0], [1], [2]]
    assert sq.plan(m).backend == "joint"


def test_mixed_timing_within_a_step_is_an_error(schema):
    m = _model(schema, 3)
    sq.join_previous(m, 2)
    m.set("transducer.timing.pd", 0.04, 2)
    for pref in sq.BACKENDS:
        assert any("different pulse timing" in i.message for i in sq.plan(m, pref).errors)


def test_partial_trains_are_an_error_for_both_backends(schema):
    m = _model(schema)
    m.set("transducer.timing.ptrd", 3.0, 1)             # 3 s / 2 s PTRI
    for pref in ("joint", "chain"):
        assert any("whole number of pulse train" in i.message for i in sq.plan(m, pref).errors)
    assert sq.whole_trains(3.0, 2.0) == pytest.approx(4.0)
    assert sq.whole_trains(0.5, 2.0) == pytest.approx(2.0)
    assert sq.whole_trains(5.0, 2.0) == pytest.approx(6.0)      # half up


def test_chain_refuses_overlap_and_off_grid_pauses(schema):
    m = _model(schema)
    m.set("transducer.timing_offset_s", 5.0, 1)          # starts before step 1 ends
    assert sq.plan(m).backend == "joint"                 # joint: interleaved firing is fine
    assert any("before step 1 ends" in i.message for i in sq.plan(m, "chain").errors)
    m = _model(schema, pause=2.5)                        # 2.5 s, time step 1 s
    assert not sq.plan(m).errors
    assert any("whole number of thermal time steps" in i.message
               for i in sq.plan(m, "chain").errors)


def test_chain_refuses_uncertainty_and_isppa_sweeps(schema):
    m = _model(schema)
    m.set("simulation.uncertainty", True)
    assert not sq.plan(m).errors
    assert any(i.path == "simulation.uncertainty" for i in sq.plan(m, "chain").errors)
    m = _model(schema)
    m.set("transducer.target_isppa_wcm2", [10, 20], 1)
    assert any(i.path == "transducer.target_isppa_wcm2" for i in sq.plan(m, "chain").errors)


def test_joint_warns_about_per_transducer_grids(schema):
    m = _model(schema)
    m.set("grid.async_layout", "per_transducer")
    assert any(i.path == "grid.async_layout" for i in sq.plan(m).issues)
    assert not any(i.path == "grid.async_layout" for i in sq.plan(m, "chain").issues)


def test_heating_off_warns(schema):
    m = _model(schema)
    m.set("modules.run_heating_sims", 0)
    assert any(i.path == "modules.run_heating_sims" and i.level == "warning"
               for i in sq.plan(m).issues)


# ── editing ───────────────────────────────────────────────────────────────────

def test_make_sequence_from_existing_transducers(schema):
    m = ConfigModel(schema)
    m.set("timing.ptrd", 10.0)
    m.set("timing.ptri", 2.0)
    m.set("timing.ptd", 1.0)
    m.set("timing.pri", 0.1)
    m.set("timing.pd", 0.02)
    m.add_transducer()
    sq.make_sequence(m)
    assert [s.start_s for s in sq.steps(m)] == pytest.approx([0, 10])


def test_remove_step_keeps_the_rest(schema):
    m = _model(schema, 3)
    sq.remove_step(m, 0)                                 # never the step with transducer 1
    assert m.n_transducers == 3
    sq.remove_step(m, 1)
    assert m.n_transducers == 2
    assert [s.start_s for s in sq.steps(m)] == pytest.approx([0, 15])
    assert m.get("transducer.trans_pos", 1) == [42, 20, 30]   # the third stance moved up
    sq.remove_step(m, 1)
    assert m.n_transducers == 1 and "timing_offset_s" not in m.tx[0]


def test_pause_edit(schema):
    m = _model(schema, 3)
    sq.set_pause(m, 0, 60.0)
    assert [s.start_s for s in sq.steps(m)] == pytest.approx([0, 70, 85])


def test_config_round_trip(schema, tmp_path):
    m = _model(schema, 3)
    sq.set_timing(m, 1, {"ptrd": 4.0})
    p = m.save(tmp_path / "config_seq.yaml")
    txs = m.override_dict()["transducer"]          # MATLAB concatenates these: same keys
    assert all(set(t) == set(txs[0]) for t in txs)
    back = ConfigModel.from_file(schema, p)
    assert [(s.transducers, s.start_s, s.duration_s) for s in sq.steps(back)] == \
        [(s.transducers, s.start_s, s.duration_s) for s in sq.steps(m)]


# ── chain of step configs ─────────────────────────────────────────────────────

def test_chain_configs(schema, tmp_path):
    m = _model(schema, 3)
    m.set("io.output_affix", "pilot")
    m.set("timing.post_ptri_dur", 30.0)
    sq.set_timing(m, 1, {"ptrd": 4.0})
    sq.join_previous(m, 2)                               # step 2: T2 + T3 together
    m.path = tmp_path / "config_seq.yaml"
    pl = sq.plan(m, "chain")
    cfgs = sq.chain_configs(m, pl)
    assert len(cfgs) == 2
    one, two = cfgs
    # step 1: one transducer, its own timing, the pause as its post_ptri_dur, base affix + step1
    assert one.n_transducers == 1 and one.get("timing.ptrd") == pytest.approx(10)
    assert one.get("timing.post_ptri_dur") == pytest.approx(5)
    assert one.get("io.output_affix") == "pilot_step1"
    assert one.is_overridden("io.preproc_affix") is False
    # step 2: both transducers as one coherent run, step 2's timing, the final steady-state
    assert two.n_transducers == 2
    assert two.get("simulation.transducer_coupling") == "coherent"
    assert two.get("timing.ptrd") == pytest.approx(4)
    assert two.get("timing.post_ptri_dur") == pytest.approx(30)
    assert two.get("io.output_affix") == "pilot_step2"
    assert two.get("io.preproc_affix") == "_desc-pilot_step2"
    assert two.get("transducer.trans_pos", 0) == [41, 20, 30]
    assert two.get("transducer.freq_hz", 1) == pytest.approx(5e5)   # inherited from T1
    for c in cfgs:
        assert c.get("modules.run_nifti_creation") == 1
        assert not c.is_overridden("timing.async_schedule")
        for ti in range(c.n_transducers):
            assert "timing_offset_s" not in c.tx[ti] and "timing" not in c.tx[ti]
        assert c.lineage["type"] == "sequence_step"


def test_chain_configs_refuse_errors(schema):
    m = _model(schema, pause=2.5)
    with pytest.raises(sq.SequenceError):
        sq.chain_configs(m, sq.plan(m, "chain"))
    with pytest.raises(sq.SequenceError):
        sq.chain_configs(m, sq.plan(m, "joint"))


def test_launch_target(schema, tmp_path):
    m = _model(schema, 3)
    p = m.save(tmp_path / "config_seq.yaml")
    opts = PipelineOptions(sequential_configs=[str(tmp_path / "later.yaml")])
    # joint: the config itself, options unchanged
    path, o, notes = sq.launch_target(m, opts)
    assert path == p and o is opts and not notes
    # chain: step 1 runs, steps 2.. first in the chain, then the study's own follow-ups
    path, o, notes = sq.launch_target(m, opts, "chain")
    assert path == sq.step_path(p, 1) and path.is_file()
    assert o.sequential_configs == [str(sq.step_path(p, 2).resolve()),
                                    str(sq.step_path(p, 3).resolve()),
                                    str(tmp_path / "later.yaml")]
    assert opts.sequential_configs == [str(tmp_path / "later.yaml")]   # not changed in place
    back = ConfigModel.from_file(schema, sq.step_path(p, 2))
    assert back.lineage["step"] == 2 and back.n_transducers == 1
    # no sequence: the config itself
    plain = ConfigModel(schema)
    plain.path = p
    assert sq.launch_target(plain, opts)[0] == p


def test_launch_target_refuses_errors(schema, tmp_path):
    m = _model(schema, 3)
    sq.join_previous(m, 2)
    m.set("grid.async_layout", "per_transducer")     # a group needs the shared grid
    m.save(tmp_path / "config_seq.yaml")
    with pytest.raises(sq.SequenceError):
        sq.launch_target(m, PipelineOptions(), "joint")


def test_backend_preference_in_session(tmp_path):
    class S(dict):
        def set(self, k, v):
            self[k] = v
    s = S()
    study = tmp_path / "config_seq.yaml"
    assert sq.load_backend(s, study) == "auto"
    sq.save_backend(s, study, "chain")
    assert sq.load_backend(s, study) == "chain"
    sq.save_backend(s, study, "auto")
    assert sq.load_backend(s, study) == "auto" and s[sq.SESSION_KEY] == {}


def test_checks_report_sequence_errors(schema):
    from prestus_nexus import checks
    m = _model(schema, 3)
    sq.join_previous(m, 2)
    m.set("transducer.timing.pd", 0.04, 2)
    assert any("different pulse timing" in i.message for i in checks.all_issues(m))
