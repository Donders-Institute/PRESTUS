"""Undo/redo of study edits (undo.History)."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")

from prestus_nexus import project_defaults  # noqa: E402
from prestus_nexus.model import ConfigModel  # noqa: E402
from prestus_nexus.undo import History  # noqa: E402


class Clock:
    def __init__(self):
        self.t = 0.0

    def __call__(self):
        return self.t


def _hist(schema):
    m = ConfigModel(schema)
    clock = Clock()
    return m, History(m, clock=clock), clock


def test_undo_redo_set_and_reset(schema):
    m, h, clock = _hist(schema)
    assert not h.can_undo() and not h.can_redo() and h.undo_text() == "Undo"
    m.set("subject_id", 5)
    clock.t = 10
    m.set("simulation.medium", "water")
    clock.t = 20
    m.reset("subject_id")
    assert h.undo_text() == "Undo reset subject_id"
    assert h.undo() == "Reset" and m.get("subject_id") == 5
    assert h.undo() == "Change" and m.get("simulation.medium") != "water"
    assert h.can_redo() and h.redo_text() == "Redo change simulation.medium"
    h.redo()
    assert m.get("simulation.medium") == "water"
    h.undo()
    h.undo()
    assert not m.overrides and not h.can_undo()
    h.redo()
    h.redo()
    h.redo()
    assert "subject_id" not in m.overrides and m.get("simulation.medium") == "water"


def test_typing_in_one_field_is_one_step(schema):
    m, h, clock = _hist(schema)
    for i, v in enumerate((1, 12, 123)):
        clock.t = i * 0.3
        m.set("subject_id", v)
    clock.t = 5                                  # a pause: a new step
    m.set("subject_id", 124)
    h.undo()
    assert m.get("subject_id") == 123
    h.undo()
    assert "subject_id" not in m.overrides and not h.can_undo()


def test_seal_and_other_keys_break_coalescing(schema):
    m, h, clock = _hist(schema)
    m.set("subject_id", 1)
    h.seal()
    m.set("subject_id", 2)
    m.set("simulation.medium", "water")
    m.set("subject_id", 3)
    assert len(h._undo) == 4


def test_new_edit_clears_redo_and_no_op_is_not_recorded(schema):
    m, h, _ = _hist(schema)
    m.set("subject_id", 1)
    h.undo()
    m.set("subject_id", m.baseline("subject_id"))   # no change
    assert h.can_redo() and not h.can_undo()
    m.set("subject_id", 2)
    assert not h.can_redo()


def test_transducer_list_edits(schema):
    m, h, clock = _hist(schema)
    ti = m.add_transducer()
    clock.t = 5
    m.set("transducer.serial", "X", ti)
    clock.t = 10
    m.remove_transducer(0)
    assert m.n_transducers == 1
    h.undo()
    assert m.n_transducers == 2 and m.tx[1].get("serial") == "X"
    h.undo()
    h.undo()
    assert m.n_transducers == 1


def test_grouped_edit_and_prune_are_one_step(schema):
    m, h, _ = _hist(schema)
    with m.edit("Apply preset"):
        m.set("subject_id", 4)
        m.set("simulation.medium", "water")
    assert len(h._undo) == 1 and h.undo_text() == "Undo apply preset"
    h.undo()
    assert not m.overrides
    m.overrides["subject_id"] = m.baseline("subject_id")    # a default written by hand
    m.prune_defaults()
    h.undo()
    assert "subject_id" in m.overrides


def test_carry_pinned_values_is_one_step(schema):
    m, h, _ = _hist(schema)
    written = project_defaults.carry_pinned_values(m, {"subject_id": (9, 1), "simulation.medium": ("water", "layered")})
    assert sorted(written) == ["simulation.medium", "subject_id"]
    h.undo()
    assert not m.overrides


def test_signals_and_attach(schema):
    m, h, _ = _hist(schema)
    changed, restored = [], []
    h.changed.connect(lambda: changed.append(1))
    h.restored.connect(restored.append)
    m.set("subject_id", 2)
    h.undo()
    assert changed and restored == ["Change"]
    other = ConfigModel(schema)
    h.attach(other)
    assert not h.can_undo() and not h.can_redo()
    m.set("subject_id", 3)                                # the old model is no longer followed
    assert not h.can_undo()
    other.set("subject_id", 4)
    assert h.can_undo()


def test_model_without_history_is_unchanged(schema):
    m = ConfigModel(schema)
    m.set("subject_id", 6)
    assert m.get("subject_id") == 6 and "_edit_hooks" not in m.__dict__
