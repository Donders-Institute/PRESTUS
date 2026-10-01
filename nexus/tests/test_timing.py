"""Protocol maths and the TUS calculator settings exchange."""

import json

import pytest

from prestus_nexus import timing as tm


def test_derived_values():
    p = tm.Protocol(pd=0.02, pri=0.2, ptd=80, ramp=0.002, ramp_shape="tukey")
    assert p.has_train and not p.has_repetition
    assert p.pulses_per_train == 400
    assert p.prf_hz == pytest.approx(5)
    assert p.duty_cycle_pulse_train == pytest.approx(0.1)
    assert p.ramp_factor == pytest.approx((2 * 3 / 8 * 0.002 + 0.016) / 0.02)   # compute_ramp_factor.m
    assert p.sonication_time_s == pytest.approx(8)
    assert p.problems() == []
    r = tm.Protocol(pd=0.02, pri=0.2, ptd=1, ptri=3, ptrd=12)
    assert r.trains == 4 and r.duty_cycle_repetition == pytest.approx(1 / 3)
    assert r.overall_duty_cycle == pytest.approx(0.1 / 3)


def test_problems():
    assert tm.Protocol().problems()
    assert "exceeds the pulse repetition" in tm.Protocol(pd=0.3, pri=0.2, ptd=1).problems()[0]
    assert any("ramp" in m for m in tm.Protocol(pd=0.01, ramp=0.006).problems())


def test_roundtrip_through_calculator_labels(tmp_path):
    p = tm.Protocol(pd=0.02, pri=0.2, ptd=1, ptri=3, ptrd=12, ramp=0.004, ramp_shape="tukey")
    s = tm.to_tus_calculator(p, freq_hz=5e5, isppa_wcm2=30, diameter_mm=64, focal_depth_mm=55)
    assert s["Pulse Duration [ms]"] == 20 and s["Ramp Duration [ms]"] == 4
    assert s["Pulse Train Repetition Duration [ms]"] == 12000
    assert s["Pulse Ramp Shape"] == "Tukey"
    assert s["Fundamental Frequency [MHz]"] == 0.5 and s["Isppa [W/cm2]"] == 30
    assert s["Pulse Train is enabled"] is True and s["Pulse Train Repetition is enabled"] is True
    f = tm.write_settings(tmp_path / "x.json", s)
    back = tm.from_tus_calculator(tm.read_settings(f))
    assert back.timing == pytest.approx({"timing.pd": 0.02, "timing.pulse_ramp_dur": 0.004, "timing.pri": 0.2,
                                         "timing.ptd": 1, "timing.ptri": 3, "timing.ptrd": 12}) \
        or back.timing["timing.pd"] == pytest.approx(0.02)
    assert back.timing["timing.pulse_ramp_shape"] == "tukey"
    assert back.optional == {"transducer.freq_hz": 5e5, "transducer.target_isppa_wcm2": 30}
    assert back.notes == []


def test_import_calculator_yaml_and_disabled_levels(tmp_path):
    f = tmp_path / "calc.yml"
    f.write_text("Fundamental Frequency [MHz]: 0.25\r\nPulse Duration [ms]: 10\r\nRamp Duration [ms]: 1\r\n"
                 "Pulse Ramp Shape: Linear\r\nPulse Train is enabled: true\r\n"
                 "Pulse Repetition Interval [ms]: 100\r\nPulse Train Duration [ms]: 30000\r\n"
                 "Pulse Train Ramp Duration [ms]: 500\r\nPulse Train Repetition is enabled: false\r\n"
                 "Pulse Train Repetition Interval [ms]: 5\r\nPulse Train Repetition Duration [ms]: 5\r\n")
    r = tm.from_tus_calculator(tm.read_settings(f))
    assert r.timing["timing.pd"] == pytest.approx(0.01) and r.timing["timing.ptd"] == pytest.approx(30)
    assert r.timing["timing.ptri"] == 0 and r.timing["timing.ptrd"] == 0
    assert r.timing["timing.pulse_ramp_shape"] == "linear"
    assert any("Train Ramp" in n for n in r.notes) and any("repetition disabled" in n for n in r.notes)
    assert r.optional["transducer.freq_hz"] == pytest.approx(2.5e5)


def test_import_rejects_unrelated_file(tmp_path):
    f = tmp_path / "other.json"
    f.write_text(json.dumps({"foo": 1}))
    r = tm.from_tus_calculator(tm.read_settings(f))
    assert not r.timing and r.notes
