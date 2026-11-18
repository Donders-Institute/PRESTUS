"""Undo and redo for the whole Nexus window (Edit ▸ Undo / Redo).

``History`` records every change made through ConfigModel.set / reset / add_transducer /
remove_transducer / prune_defaults, and any block wrapped in ``model.edit(label)`` or
``history.group(label)``. Other editable state (the Study tab's runs table) joins with
``history.track(target)``: the target's snapshot() / restore(snap) are saved and put back
with the config's, so one linear history covers every tab. A step holds the state before
and after the change. Changes to the same key within COALESCE_S of each other (typing in
a field) form one step.

    history = History(model)
    history.track(study_page)                   # snapshot() / restore(snap)
    history.changed.connect(update_menu)        # can_undo / can_redo / undo_text changed
    history.restored.connect(refresh_pages)     # the model's values were replaced
    with history.group("Add runs"): ...         # one step
    history.group_pass("Move focus")            # one step: everything until the event loop runs
    history.undo(); history.redo()

Opening another study: ``history.attach(new_model)`` (clears the history). Saving does
not clear it; ``seal()`` ends coalescing so the next edit starts a new step.
"""
from __future__ import annotations

import contextlib
import time
from dataclasses import dataclass

from PyQt6.QtCore import QObject, QTimer, pyqtSignal

COALESCE_S = 1.5
LIMIT = 200


@dataclass
class Step:
    label: str
    key: tuple
    before: tuple
    after: tuple
    stamp: float
    sealed: bool = False


class History(QObject):
    changed = pyqtSignal()            # undo/redo availability or texts changed
    restored = pyqtSignal(str)        # undo/redo replaced the model's values (step label)

    def __init__(self, model=None, limit: int = LIMIT, coalesce_s: float = COALESCE_S,
                 clock=time.monotonic, parent=None):
        super().__init__(parent)
        self.model = None
        self.limit = limit
        self.coalesce_s = coalesce_s
        self._clock = clock
        self._undo: list[Step] = []
        self._redo: list[Step] = []
        self._depth = 0
        self._open: tuple[str, tuple, tuple] | None = None     # (label, key, before)
        self._applying = False
        self._pass_open = False
        self._targets: list = []
        if model is not None:
            self.attach(model)

    def track(self, target) -> None:
        """Record target's state with every step (target.snapshot() / target.restore(snap)).
        Track before the first edit: steps recorded earlier do not hold its state."""
        if target not in self._targets:
            self._targets.append(target)

    def _snapshot(self) -> tuple:
        return (self.model.snapshot(), tuple(t.snapshot() for t in self._targets))

    # ── model ───────────────────────────────────────────────────────────
    def attach(self, model) -> None:
        """Follow model (another study opened); the history starts empty."""
        if self.model is not None:
            self.model.remove_edit_hook(self)
        self.model = model
        self._depth, self._open, self._pass_open = 0, None, False
        if model is not None:
            model.add_edit_hook(self)
        self.clear()

    def clear(self) -> None:
        self._undo.clear()
        self._redo.clear()
        self.changed.emit()

    def seal(self) -> None:
        """The next edit starts a new step, even on the same key."""
        if self._undo:
            self._undo[-1].sealed = True

    @contextlib.contextmanager
    def group(self, label: str = "Edit", key: tuple = ()):
        """The changes made inside (config and tracked targets) form one step."""
        self.begin(label, key)
        try:
            yield
        finally:
            self.end()

    def group_pass(self, label: str = "Edit") -> None:
        """The changes made until control returns to the event loop form one step: a click
        or button whose handler writes several values, however it is split into calls."""
        if self._pass_open or self._applying or self.model is None:
            return
        self._pass_open = True
        self.begin(label, ())
        QTimer.singleShot(0, self._end_pass)

    def _end_pass(self) -> None:
        if self._pass_open:
            self._pass_open = False
            self.end()

    # hooks called by ConfigModel.edit
    def begin(self, label: str, key: tuple) -> None:
        if self._applying:
            return
        if self._depth == 0:
            self._open = (label, tuple(key), self._snapshot())
        self._depth += 1

    def end(self) -> None:
        if self._applying or self._depth == 0:
            return
        self._depth -= 1
        if self._depth or self._open is None:
            return
        label, key, before = self._open
        self._open = None
        after = self._snapshot()
        if after == before:
            return
        now = self._clock()
        top = self._undo[-1] if self._undo else None
        if top is not None and not top.sealed and key and top.key == key and \
                now - top.stamp < self.coalesce_s:
            top.after, top.stamp = after, now
        else:
            self._undo.append(Step(label, key, before, after, now))
            del self._undo[:-self.limit]
        self._redo.clear()
        self.changed.emit()

    # ── undo / redo ─────────────────────────────────────────────────────
    def can_undo(self) -> bool:
        return bool(self._undo)

    def can_redo(self) -> bool:
        return bool(self._redo)

    def undo_text(self) -> str:
        """'Undo <what>' for the menu ('Undo' when there is nothing to undo)."""
        return f"Undo {_describe(self._undo[-1])}" if self._undo else "Undo"

    def redo_text(self) -> str:
        return f"Redo {_describe(self._redo[-1])}" if self._redo else "Redo"

    def undo(self) -> str | None:
        """Undo the last step; returns its label (None if there was nothing to undo)."""
        self._end_pass()                    # a pending click action is a step first
        if not self._undo:
            return None
        step = self._undo.pop()
        self._apply(step.before)
        step.sealed = True
        self._redo.append(step)
        self.changed.emit()
        self.restored.emit(step.label)
        return step.label

    def redo(self) -> str | None:
        self._end_pass()
        if not self._redo:
            return None
        step = self._redo.pop()
        self._apply(step.after)
        self._undo.append(step)
        self.changed.emit()
        self.restored.emit(step.label)
        return step.label

    def _apply(self, snap: tuple) -> None:
        self._applying = True
        try:
            model_snap, target_snaps = snap
            self.model.restore(model_snap)
            for target, s in zip(self._targets, target_snaps):
                target.restore(s)
        finally:
            self._applying = False


def _describe(step: Step) -> str:
    """'change transducer.focus_depth' style text for a step."""
    path = step.key[1] if len(step.key) > 1 and isinstance(step.key[1], str) else ""
    return f"{step.label.lower()} {path}".strip() if path else step.label.lower()
