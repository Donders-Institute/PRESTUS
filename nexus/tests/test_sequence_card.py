"""Sonication › Timing › Sonication sequence card, offscreen; and its main-window wiring."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")

from PyQt6.QtWidgets import QApplication  # noqa: E402

from prestus_nexus import sequence as sq  # noqa: E402
from prestus_nexus.model import ConfigModel  # noqa: E402
from prestus_nexus.session import Session  # noqa: E402
from prestus_nexus.ui.sequence_card import COL_DURATION, COL_PAUSE, SequenceCard  # noqa: E402


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def _model(schema):
    m = ConfigModel(schema)
    for key, v in {"timing.pd": 0.02, "timing.pri": 0.1, "timing.ptd": 1.0, "timing.ptri": 2.0,
                   "timing.ptrd": 10.0, "timing.post_pt_timestep": 1.0,
                   "modules.run_heating_sims": 1}.items():
        m.set(key, v)
    m.set("transducer.trans_pos", [10, 20, 30])
    return m


def test_one_stance_offers_a_later_sonication(app, schema):
    m = _model(schema)
    card = SequenceCard(m)
    assert card.body.isHidden() and not card.empty.isHidden()
    placed = []
    card.edit_in_viewer.connect(placed.append)
    card.start_btn.click()
    assert sq.is_sequence(m) and m.n_transducers == 2
    assert not card.body.isHidden() and card.table.rowCount() == 2
    assert placed == [1]
    assert "T2 has transducer 1's position" in card.place_note.text()
    m.set("transducer.trans_pos", [50, 20, 30], 1)
    card.refresh()
    assert card.place_note.isHidden()


def test_several_transducers_offer_one_after_another(app, schema):
    m = _model(schema)
    m.add_transducer()
    card = SequenceCard(m)
    assert "fire together" in card.empty_text.text()
    card.start_btn.click()
    assert [s.start_s for s in sq.steps(m)] == pytest.approx([0, 10])


def test_table_edits_duration_and_pause(app, schema):
    m = _model(schema)
    sq.add_step(m)
    card = SequenceCard(m)
    card.table.item(0, COL_PAUSE).setText("30 s")
    assert sq.steps(m)[1].start_s == pytest.approx(40)
    card.table.item(1, COL_DURATION).setText("5 s")          # rounded to whole 2 s trains
    assert sq.steps(m)[1].duration_s == pytest.approx(6)
    assert card.table.item(1, COL_PAUSE).text() == "–"       # nothing after the last step
    card.table.item(0, COL_PAUSE).setText("soon")              # unreadable: unchanged
    assert sq.steps(m)[1].start_s == pytest.approx(40)


def test_buttons_follow_the_selection(app, schema):
    m = _model(schema)
    sq.add_step(m)
    sq.add_step(m)
    card = SequenceCard(m)
    card.table.selectRow(0)
    assert not card.join_btn.isEnabled() and not card.remove_btn.isEnabled()
    card.table.selectRow(2)
    assert card.join_btn.isEnabled() and card.remove_btn.isEnabled()
    card.join_btn.click()
    assert [s.transducers for s in sq.steps(m)] == [[0], [1, 2]]
    assert "summed coherently" in card.status.text()     # a firing group of one joint run
    card.table.selectRow(1)
    assert card.split_btn.isEnabled()
    card.split_btn.click()
    assert "One joint run" in card.status.text()


def test_backend_choice_is_stored_per_study(app, schema, tmp_path):
    m = _model(schema)
    sq.add_step(m)
    m.path = tmp_path / "config_seq.yaml"
    session = Session(tmp_path / "session.yaml")
    card = SequenceCard(m, lambda: (session, m.path))
    assert card.backend.isEnabled()
    card.backend.setCurrentIndex(card.backend.findData("chain"))
    assert sq.load_backend(session, m.path) == "chain"
    assert "A chain of 2 runs" in card.status.text()
    unsaved = SequenceCard(_model(schema))
    sq.add_step(unsaved.model)
    unsaved.refresh()
    assert not unsaved.backend.isEnabled()


def test_errors_show_in_the_card(app, schema):
    m = _model(schema)
    sq.add_step(m, 2.5)                       # off the 1 s thermal time step
    card = SequenceCard(m)
    assert "whole number of thermal time steps" not in card.status.text()   # joint: fine
    m.path = None
    session = Session.__new__(Session)
    session.data = {}
    session.path = None
    card.session_getter = lambda: (session, "/tmp/config_seq.yaml")
    sq.save_backend(session, "/tmp/config_seq.yaml", "chain")
    card.refresh()
    assert "whole number of thermal time steps" in card.status.text()


def test_main_window_has_the_card(app, schema, prestus_root, tmp_path):
    from prestus_nexus.ui.main_window import MainWindow
    w = MainWindow(ConfigModel(schema), Session(tmp_path / "session.yaml"), prestus_root)
    try:
        if not hasattr(w, "sequence_card"):
            w.goto("timing.pd")                 # timing page built on first show
        assert hasattr(w, "sequence_card")
        w.sequence_card.start_btn.click()
        assert w.model.n_transducers == 2 and w.dirty
    finally:
        w.dirty = False
        w.close()


def test_timing_diagram_draws_heads_beside_async_lane(app, schema):
    """The async lane gets the head sketch (taller lane, hover text) when transducers have
    positions, and none without them or without async coupling."""
    from PyQt6.QtCore import QPointF, Qt
    from PyQt6.QtGui import QMouseEvent
    from prestus_nexus.ui.timing_view import TimingPanel
    m = _model(schema)
    m.set("transducer.focus_pos", [12, 60, 30])
    ti = sq.add_step(m)
    m.set("transducer.trans_pos", [60, 20, 30], ti)
    m.set("transducer.focus_pos", [40, 60, 30], ti)
    panel = TimingPanel(m)
    d = panel.diagram
    assert len(d.heads) == 2 and all(h is not None for h in d.heads)
    with_heads = d.height()
    d.resize(900, with_heads)
    img = d.grab()                                # paints: sets the sketch's box
    assert not img.isNull() and d._head_rect is not None
    c = d._head_rect.center()
    d.mouseMoveEvent(QMouseEvent(QMouseEvent.Type.MouseMove, QPointF(c), Qt.MouseButton.NoButton,
                                 Qt.MouseButton.NoButton, Qt.KeyboardModifier.NoModifier))
    assert d.toolTip().startswith("Seen from above")
    # no positions: no sketch
    for i in range(m.n_transducers):
        m.reset("transducer.trans_pos", i)
    panel.refresh()
    assert not d._has_heads()
    # no async coupling: no lane, no sketch
    m.set("simulation.transducer_coupling", "coherent")
    panel.refresh()
    assert d.heads == [] and d.sched is None


def test_coupling_prediction_hidden_for_a_sequence(app, schema):
    from prestus_nexus.ui.coupling_view import CouplingPrediction
    m = _model(schema)
    m.add_transducer()
    w = CouplingPrediction(m, lambda: None)
    assert not w.isHidden()                       # two transducers firing together
    sq.make_sequence(m)
    w.refresh()
    assert w.isHidden()
    w.toggle.setChecked(True)
    assert w.text.text() == ""


def test_card_toggles_open_and_closed(app, schema):
    """Closed by default, also with a sequence; the toggle and adding a step open it."""
    m = _model(schema)
    card = SequenceCard(m)
    assert not card.is_open() and card.content.isHidden()
    card.toggle.click()
    assert card.is_open() and not card.content.isHidden()
    card.toggle.click()
    assert not card.is_open()
    card.start_btn.click()                       # adding a step opens it
    assert card.is_open() and not card.content.isHidden()
    seq = _model(schema)
    sq.add_step(seq)
    assert not SequenceCard(seq).is_open()       # closed by default, also with a sequence


def test_strip_colours_follow_the_transducers(app, schema):
    """Each step block in the colour of its transducer (as the async lanes and head sketch);
    a step that fires several transducers in one band per transducer."""
    from PyQt6.QtGui import QColor
    from prestus_nexus.ui import style
    m = _model(schema)
    for _ in range(3):
        sq.add_step(m)
    sq.join_previous(m, 2)                       # steps: [T1], [T2, T3], [T4]
    card = SequenceCard(m)
    strip = card.strip
    strip.resize(900, strip.height())
    img = strip.grab().toImage()
    pl = strip.plan
    sx = (strip.width() - 12.0) / pl.total_s

    def colour_at(t, frac):                      # time [s], height fraction within the block
        return QColor(img.pixel(int(6 + t * sx), int(4 + 26 * frac))).name()

    for st in pl.steps:
        mid = st.start_s + 0.1 * st.duration_s     # away from the centred label
        n = len(st.transducers)
        for j, ti in enumerate(st.transducers):
            assert colour_at(mid, (j + 0.5) / n) == QColor(style.series(ti)).name()


def test_protocol_diagram_follows_the_selected_step(app, schema):
    """The protocol lanes show the step selected in the sequence table; edits in the diagram
    change that step and move the later steps, keeping the pauses."""
    from prestus_nexus.ui.timing_view import TimingPanel
    m = _model(schema)
    sq.add_step(m, 30.0)
    sq.add_step(m, 30.0)
    sq.set_timing(m, 2, {"pd": 0.05, "ptrd": 4.0})
    panel = TimingPanel(m)
    card = SequenceCard(m)
    card.step_selected.connect(panel.show_step)
    assert panel.diagram.p.pd == pytest.approx(0.02) and "step 1" in panel.summary.text()
    card.table.selectRow(2)
    assert panel.diagram.p.pd == pytest.approx(0.05) and panel.diagram.p.ptrd == pytest.approx(4)
    assert "step 3" in panel.summary.text()
    # editing step 1 in the diagram: its duration changes, later steps move, pauses stay
    card.table.selectRow(0)
    panel._apply({"timing.ptrd": 6.0})
    st = sq.steps(m)
    assert st[0].duration_s == pytest.approx(6) and sq.pauses(st) == pytest.approx([30, 30])
    # editing step 3: its own override, step 1 unchanged
    panel.show_step(2)
    panel._apply({"timing.pd": 0.04})
    assert sq.steps(m)[2].protocol.pd == pytest.approx(0.04) and m.get("timing.pd") == pytest.approx(0.02)
    # the ramp is one for every step: the timing block
    panel._apply({"timing.pulse_ramp_dur": 0.002})
    assert m.get("timing.pulse_ramp_dur") == pytest.approx(0.002)
    # no sequence: the timing block as before
    plain = _model(schema)
    p2 = TimingPanel(plain)
    p2._apply({"timing.ptrd": 12.0})
    assert plain.get("timing.ptrd") == pytest.approx(12) and "Showing step" not in p2.summary.text()
