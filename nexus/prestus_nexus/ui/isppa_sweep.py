"""transducer.target_isppa_wcm2 editor: the usual line edit plus a "Sweep…" button that edits
a list of targets. More than one value runs the multi-ISPPA pipeline
(is_multi_isppa_mode.m): one acoustic simulation, one thermal run per target
(doc_multiintensity.md)."""

from __future__ import annotations

import math

from PyQt6.QtWidgets import (QDialog, QDialogButtonBox, QDoubleSpinBox, QHBoxLayout, QLabel,
                             QListWidget, QPushButton, QToolButton, QVBoxLayout, QWidget)

from ..schema import Node
from .fields import TextEditor, _num_text


def sweep_problems(values: list, n_transducers: int = 1, coupling: str | None = None) -> list[str]:
    """Problems with a list of targets: positive finite numbers; a sweep with several coherent
    transducers is refused (checks.rule_issues: coherent_drive_calibration.m /
    multi_isppa_pipeline.m set each drive from one target per transducer)."""
    out = []
    nums = [v for v in values if isinstance(v, (int, float)) and not isinstance(v, bool)]
    if len(nums) != len(values) or any(not math.isfinite(v) or v <= 0 for v in nums):
        out.append("Every target must be a positive number [W/cm²]")
    if len(set(nums)) != len(nums):
        out.append("Duplicate targets: each value runs once")
    if len(nums) > 1 and n_transducers > 1 and coupling != "async":
        out.append("Coherent multi-transducer runs take one target per transducer; use async "
                   "coupling for a sweep")
    return out


class SweepDialog(QDialog):
    def __init__(self, values: list, n_transducers: int = 1, coupling: str | None = None,
                 parent: QWidget | None = None):
        super().__init__(parent)
        self.setWindowTitle("ISPPA sweep")
        self._ctx = (n_transducers, coupling)
        v = QVBoxLayout(self)
        intro = QLabel("Targets of a multi-ISPPA sweep: one acoustic simulation, then one thermal "
                       "run and report per target, plus a summary report.")
        intro.setWordWrap(True)
        v.addWidget(intro)
        self.list = QListWidget()
        v.addWidget(self.list, 1)
        row = QHBoxLayout()
        self.spin = QDoubleSpinBox()
        self.spin.setRange(0.01, 1000.0)
        self.spin.setDecimals(2)
        self.spin.setSuffix(" W/cm²")
        self.spin.setValue(10.0)
        row.addWidget(self.spin)
        add = QPushButton("Add")
        add.clicked.connect(lambda: self._add(self.spin.value()))
        row.addWidget(add)
        rem = QPushButton("Remove")
        rem.clicked.connect(self._remove)
        row.addWidget(rem)
        srt = QPushButton("Sort")
        srt.clicked.connect(self._sort)
        row.addWidget(srt)
        row.addStretch(1)
        v.addLayout(row)
        self.problems = QLabel()
        self.problems.setWordWrap(True)
        v.addWidget(self.problems)
        self.buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok |
                                        QDialogButtonBox.StandardButton.Cancel)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        v.addWidget(self.buttons)
        for x in values:
            self._add(x, check=False)
        self._check()

    def values(self) -> list[float]:
        return [float(self.list.item(i).text()) for i in range(self.list.count())]

    def _add(self, x, check: bool = True) -> None:
        if isinstance(x, (int, float)) and not isinstance(x, bool):
            self.list.addItem(_num_text(float(x)))
        if check:
            self._check()

    def _remove(self) -> None:
        r = self.list.currentRow()
        if r >= 0:
            self.list.takeItem(r)
            self._check()

    def _sort(self) -> None:
        vals = sorted(self.values())
        self.list.clear()
        for x in vals:
            self._add(x, check=False)
        self._check()

    def _check(self) -> None:
        from . import style
        probs = sweep_problems(self.values(), *self._ctx)
        if not self.values():
            probs = ["Add at least one target"]
        self.problems.setText("<br>".join(f"<span style='color:{style.active()['red']}'>{p}</span>"
                                          for p in probs))
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(not probs)


class IsppaSweepEditor(TextEditor):
    """TextEditor for target_isppa_wcm2 with a "Sweep…" button."""

    def __init__(self, node: Node, parent: QWidget | None = None):
        super().__init__(node, parent)
        self.sweep = QToolButton(text="Sweep…")
        self.sweep.setToolTip("Edit several targets (multi-ISPPA sweep)")
        self.sweep.clicked.connect(self._edit_sweep)
        self._layout.addWidget(self.sweep)

    def _context(self) -> tuple[int, str | None]:
        model = getattr(self.window(), "model", None)
        if model is None:
            return 1, None
        return model.n_transducers, model.get("simulation.transducer_coupling")

    def _edit_sweep(self) -> None:
        try:
            current = self.parse(self.edit.text())
        except Exception:  # best-effort: unparsable text: the sweep starts empty
            current = []
        values = current if isinstance(current, list) else [current]
        dlg = SweepDialog([v for v in values if isinstance(v, (int, float))], *self._context(), parent=self)
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return
        vals = dlg.values()
        value = vals[0] if len(vals) == 1 else vals
        self.edit.setText(", ".join(_num_text(v) for v in vals))
        self._shown = self.edit.text()
        self.set_error(None)
        self.edited.emit(value)
