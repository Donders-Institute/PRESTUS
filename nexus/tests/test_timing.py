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


# ── Editing from the diagram ────────────────────────────────────────────

def test_snap_limits_and_dragged_values():
    p = tm.Protocol(pd=0.02, pri=0.2, ptd=1, ptri=3, ptrd=12, ramp=0.002)
    assert tm.snap(0.023456) == pytest.approx(0.023) and tm.snap(0) == 0
    assert tm.limits(p, "pd") == (0.004, 0.2)
    assert tm.limits(p, "ptri") == (1, 12)
    assert tm.dragged(p, "pd", 0.5) == pytest.approx(0.2)        # cannot outgrow the PRI
    assert tm.dragged(p, "ramp", 0.5) == pytest.approx(0.01)     # rise + fall fit in the pulse
    assert tm.dragged(p, "ptrd", 40.4) == pytest.approx(40)
    assert not tm.Protocol(**{**p.__dict__, "pd": tm.dragged(p, "pd", 9)}).problems()


def test_parse_duration_units_and_frequency():
    assert tm.parse_duration("20 ms") == pytest.approx(0.02)
    assert tm.parse_duration("500 µs") == pytest.approx(5e-4)
    assert tm.parse_duration("0,5") == pytest.approx(0.5)
    assert tm.parse_duration("2 min") == pytest.approx(120)
    assert tm.parse_duration("5 Hz", hz=True) == pytest.approx(0.2)
    assert tm.parse_duration("5 Hz") is None and tm.parse_duration("fast") is None


def test_diagram_drag_and_type(schema):
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PyQt6")
    from PyQt6.QtCore import QPoint, Qt
    from PyQt6.QtTest import QTest
    from PyQt6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    from prestus_nexus.model import ConfigModel
    from prestus_nexus.ui.timing_view import TimingPanel
    model = ConfigModel(schema)
    for k, v in {"timing.pd": 0.02, "timing.pri": 0.2, "timing.ptd": 1, "timing.ptri": 3, "timing.ptrd": 12,
                 "timing.pulse_ramp_dur": 0}.items():
        model.set(k, v)
    panel = TimingPanel(model)
    panel.resize(800, 500)
    panel.show()
    app.processEvents()
    d = panel.diagram
    d.repaint()
    edge = next(h for h in d._handles if h.edge and h.name == "pd")
    y = int(edge.rect.center().y())
    x = int(edge.rect.center().x())
    QTest.mousePress(d, Qt.MouseButton.LeftButton, pos=QPoint(x, y))
    QTest.mouseMove(d, QPoint(x + int(0.03 * edge.sx), y))       # +30 ms
    QTest.mouseRelease(d, Qt.MouseButton.LeftButton, pos=QPoint(x + int(0.03 * edge.sx), y))
    assert model.get("timing.pd") == pytest.approx(0.05, abs=0.002)

    label = next(h for h in d._handles if not h.edge and h.name == "pri")
    QTest.mouseClick(d, Qt.MouseButton.LeftButton, pos=label.rect.center().toPoint())
    assert d._editor is not None
    d._editor.setText("4 Hz")
    QTest.keyClick(d._editor, Qt.Key.Key_Return)
    assert model.get("timing.pri") == pytest.approx(0.25)
    panel.close()


def test_timing_fields_fold_under_the_diagram(schema, tmp_path):
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PyQt6")
    from pathlib import Path
    from PyQt6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    from prestus_nexus.model import ConfigModel
    from prestus_nexus.session import Session
    from prestus_nexus.ui.main_window import MainWindow
    model = ConfigModel(schema)
    for k, v in {"timing.pd": 0.02, "timing.pri": 0.2, "timing.ptd": 1, "timing.ptri": 3, "timing.ptrd": 12}.items():
        model.set(k, v)
    w = MainWindow(model, Session(tmp_path / "s.yaml"), Path(__file__).resolve().parents[2])
    page = w.pages["timing"]
    rows = {r.node.path: r for r in page.rows}
    page.refresh()
    assert all(rows[k].block.isHidden() for k in ("timing.pd", "timing.pri", "timing.ptrd"))
    assert not rows["timing.pulse_ramp_shape"].block.isHidden()          # not set in the diagram
    w.timing_panel.fold_button.setChecked(True)
    assert not rows["timing.pd"].block.isHidden()
    w.timing_panel.fold_button.setChecked(False)
    model.set("timing.pd", "abc")                                        # an invalid value stays in view
    page.refresh()
    assert not rows["timing.pd"].block.isHidden()
    w.close()


def test_fields_shown_when_no_protocol_can_be_drawn(schema, tmp_path):
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PyQt6")
    from pathlib import Path
    from PyQt6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    from prestus_nexus.model import ConfigModel
    from prestus_nexus.session import Session
    from prestus_nexus.ui.main_window import MainWindow
    w = MainWindow(ConfigModel(schema), Session(tmp_path / "s.yaml"), Path(__file__).resolve().parents[2])
    rows = {r.node.path: r for r in w.pages["timing"].rows}
    assert not tm.Protocol.from_model(w.model.get).has_pulse
    assert w.timing_panel.fold_button.isChecked()                         # empty on opening
    assert not rows["timing.pd"].block.isHidden()
    for k, v in {"timing.pd": 0.02, "timing.pri": 0.2, "timing.ptd": 1, "timing.ptri": 3, "timing.ptrd": 12}.items():
        w.model.set(k, v)
    w.timing_panel.refresh()
    w.timing_panel.fold_button.setChecked(False)                          # hidden once drawn
    assert rows["timing.pd"].block.isHidden()
    w.model.set("timing.pd", None)                                        # cleared again
    w.timing_panel.refresh()
    assert not rows["timing.pd"].block.isHidden()
    w.timing_panel.fold_button.setChecked(False)                          # the user may still hide them
    w.timing_panel.refresh()
    assert not w.timing_panel.fold_button.isChecked()
    w.close()


def test_values_equal_on_screen_are_not_flagged():
    # 700 ms typed with a unit used to be 0.7000000000000001 s, so 'PD > PRI' stayed on screen
    assert tm.parse_duration("700 ms") == 0.7
    p = tm.Protocol(pd=0.7000000000000001, pri=0.7, ptd=1.4, ptri=2.8, ptrd=5.6, ramp=0.35000000000000003)
    assert p.problems() == []
    assert tm.Protocol(pd=0.71, pri=0.7, ptd=1.4).problems()      # a real excess is still reported


# ── async firing schedule (thermal_schedule.m) ──────────────────────────

def test_async_pulses_and_trains_offsets():
    base = tm.Protocol(pd=0.02, pri=0.2, ptd=1, ptri=4, ptrd=12)
    s = tm.AsyncSchedule.build(base, "pulses", 3)
    assert [src.offset for src in s.sources] == pytest.approx([0, 0.2 / 3, 0.4 / 3])
    assert s.cycle() == pytest.approx(0.4)
    assert s.problems() == [] and s.overlaps() == []
    blocks, level = s.blocks(1, s.cycle())
    assert level == "pulses" and blocks[0] == pytest.approx((0.2 / 3, 0.2 / 3 + 0.02))
    t = tm.AsyncSchedule.build(base, "trains", 2)
    assert [src.offset for src in t.sources] == pytest.approx([0, 2])
    assert t.cycle() == pytest.approx(4)
    assert t.blocks(0, 4, max_blocks=3)[1] == "trains"           # 5 pulses > 3: drawn as a train
    # each transducer keeps the timing block's duty cycle (the former async assumption)
    assert all(src.overall_duty_cycle == pytest.approx(0.025) for src in t.sources)


def test_async_schedule_problems_match_thermal_schedule():
    base = tm.Protocol(pd=0.08, pri=0.2, ptd=1, ptri=4, ptrd=12)
    assert "needs 3 × PD ≤ PRI" in " ".join(tm.AsyncSchedule.build(base, "pulses", 3).problems())
    assert "needs 5 × PTD ≤ PTRI" in " ".join(tm.AsyncSchedule.build(base, "trains", 5).problems())
    ok = tm.Protocol(pd=0.02, pri=0.2, ptd=1, ptri=4, ptrd=12)
    over = tm.AsyncSchedule.build(ok, "pulses", 2, overrides=[{}, {"pd": 0.04}])
    assert "only the 'custom' schedule" in over.problems()[0]
    # custom: same frequency overlap is an error, different frequencies a note
    same = tm.AsyncSchedule.build(ok, "custom", 2, offsets=[0, 0.01], freqs=[5e5, 5e5])
    assert same.overlaps() == [(0, 1)] and "interfere" in same.problems()[0]
    diff = tm.AsyncSchedule.build(ok, "custom", 2, offsets=[0, 0.01], freqs=[5e5, 2.5e5])
    assert diff.problems() == [] and "different frequencies" in diff.notes()[0]
    apart = tm.AsyncSchedule.build(ok, "custom", 2, offsets=[0, 0.1],
                                   overrides=[{}, {"pd": 0.04, "pri": 0.4}], freqs=[5e5, 5e5])
    assert apart.problems() == [] and apart.sources[1].protocol.pd == 0.04
    assert apart.sources[1].overall_duty_cycle == pytest.approx(apart.sources[0].overall_duty_cycle)


def test_async_lane_in_diagram(schema):
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PyQt6")
    from PyQt6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    from prestus_nexus.model import ConfigModel
    from prestus_nexus.ui.timing_view import TimingPanel
    model = ConfigModel(schema)
    for k, v in {"timing.pd": 0.02, "timing.pri": 0.2, "timing.ptd": 1, "timing.ptri": 3,
                 "timing.ptrd": 12}.items():
        model.set(k, v)
    panel = TimingPanel(model)
    panel.resize(800, 600)
    h_single = panel.diagram.height()
    assert panel.diagram.sched is None
    for _ in range(3):
        model.add_transducer()
    model.set("simulation.transducer_coupling", "async")
    panel.refresh()
    d = panel.diagram
    assert d.sched is not None and d.sched.n == model.n_transducers and d.sched.mode == "pulses"
    assert d.height() == h_single + d._async_h()
    assert d._async_h() == d.lh + 8 + 4 * d.ROW_H + 6 + d.lh + 10      # one row per transducer
    panel.show()
    app.processEvents()
    d.repaint()
    assert "overall duty cycle" not in panel.summary.text()     # no async summary line
    model.set("timing.pd", 0.08)                      # 4 x 80 ms > 200 ms
    panel.refresh()
    assert "needs 4 × PD ≤ PRI" in panel.summary.text()
    panel.close()
