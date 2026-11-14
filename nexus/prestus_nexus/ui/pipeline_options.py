"""Run tab: pipeline options, the second argument of prestus_pipeline_start that the run
script passes (launcher.PipelineOptions). Nexus-side settings, not config keys: they are
stored per study in the session and written into run_<study>.m / submit_<study>.m.

Collapsed by default. Shows the sequential chain (follow-up configs run after this one
with heat carry-over, doc_sequential.md), and, when the config uses them, the per-stage time
and memory limits of the multi-ISPPA and uncertainty pipelines and the multi-ISPPA summary
report switch."""

from __future__ import annotations

import re
from pathlib import Path

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (QCheckBox, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout, QLabel,
                             QLineEdit, QListWidget, QPushButton, QSpinBox, QToolButton,
                             QVBoxLayout, QWidget)

from .. import launcher

_TIME = re.compile(r"^\d{1,3}:[0-5]\d:[0-5]\d$")

# (option prefix, label) per pipeline: <prefix>_timelimit / <prefix>_memorylimit
MULTI_STAGES = (("acoustic", "Acoustic (stage 1)"), ("thermal", "Thermal, per ISPPA"),
                ("summary", "Summary report"))
UNCERTAINTY_STAGES = (("stage1", "Preprocessing (stage 1)"), ("sim", "Variant simulations"),
                      ("report", "Uncertainty report"))
# PRESTUS defaults (multi_isppa_pipeline.m:111-113, uncertainty_pipeline.m:122-124)
DEFAULT_TIME = {"acoustic": "06:00:00", "thermal": "02:00:00", "summary": "00:30:00",
                "stage1": "12:00:00", "sim": "06:00:00", "report": "00:30:00"}


def is_multi_isppa(model) -> bool:
    """is_multi_isppa_mode.m: a transducer with more than one finite target."""
    for t in model.effective_transducers():
        v = t.get("target_isppa_wcm2")
        if isinstance(v, list) and len(v) > 1 and any(
                isinstance(x, (int, float)) and not isinstance(x, bool) for x in v):
            return True
    return False


class PipelineOptionsPanel(QWidget):
    """Compact form bound to the open study's stored PipelineOptions (window.session)."""

    changed = pyqtSignal()

    def __init__(self, window, parent: QWidget | None = None):
        super().__init__(parent)
        self.window = window
        self._setting = False
        v = QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 0)
        self.toggle = QToolButton()
        self.toggle.setText("Pipeline options")
        self.toggle.setCheckable(True)
        self.toggle.setAutoRaise(True)
        self.toggle.setToolTip("Options the run script passes to prestus_pipeline_start: sequential "
                               "follow-ups and per-stage limits. Stored with the Nexus, not in the config.")
        self.toggle.toggled.connect(self._toggled)
        v.addWidget(self.toggle)
        self.body = QWidget()
        v.addWidget(self.body)
        b = QVBoxLayout(self.body)
        b.setContentsMargins(12, 0, 0, 0)

        # sequential chain
        seq = QGroupBox("Sequential follow-ups (heat carry-over)")
        sv = QVBoxLayout(seq)
        self.chain = QListWidget()
        self.chain.setMaximumHeight(110)
        self.chain.setToolTip("Run after this config, in order; each starts from the previous run's end "
                              "temperature and thermal dose (options.sequential_configs.config_2, _3, …)")
        self.chain.currentRowChanged.connect(lambda _r: self._buttons())
        sv.addWidget(self.chain)
        row = QHBoxLayout()
        self.add_btn = QPushButton("Add…")
        self.add_btn.clicked.connect(self._add)
        self.remove_btn = QPushButton("Remove")
        self.remove_btn.clicked.connect(self._remove)
        self.up_btn = QPushButton("Up")
        self.up_btn.clicked.connect(lambda: self._move(-1))
        self.down_btn = QPushButton("Down")
        self.down_btn.clicked.connect(lambda: self._move(1))
        for w in (self.add_btn, self.remove_btn, self.up_btn, self.down_btn):
            row.addWidget(w)
        row.addStretch(1)
        sv.addLayout(row)
        self.chain_note = QLabel()
        self.chain_note.setObjectName("CardSub")
        self.chain_note.setWordWrap(True)
        sv.addWidget(self.chain_note)
        b.addWidget(seq)

        self.limits: dict[str, tuple[QLineEdit, QSpinBox]] = {}
        self.multi_box = self._stage_box("Multi-ISPPA stages (HPC)", MULTI_STAGES)
        self.summary = QCheckBox("Write the multi-ISPPA summary report")
        self.summary.setToolTip("options.generate_summary_report (default on)")
        self.summary.toggled.connect(self._store)
        self.multi_box.layout().addRow(self.summary)
        b.addWidget(self.multi_box)
        self.unc_box = self._stage_box("Uncertainty stages (HPC)", UNCERTAINTY_STAGES)
        b.addWidget(self.unc_box)
        self.toggle.setChecked(False)
        self.body.setVisible(False)
        self.refresh()

    def _stage_box(self, title: str, stages) -> QGroupBox:
        box = QGroupBox(title)
        f = QFormLayout(box)
        for prefix, label in stages:
            t = QLineEdit()
            t.setPlaceholderText(DEFAULT_TIME[prefix])
            t.setToolTip(f"options.{prefix}_timelimit, HH:MM:SS (empty = {DEFAULT_TIME[prefix]})")
            t.setMaximumWidth(110)
            t.editingFinished.connect(self._store)
            m = QSpinBox()
            m.setRange(0, 2048)
            m.setSuffix(" GB")
            m.setSpecialValueText("hpc.memorylimit")
            m.setToolTip(f"options.{prefix}_memorylimit (0 = the config's hpc.memorylimit)")
            m.valueChanged.connect(self._store)
            row = QHBoxLayout()
            row.addWidget(t)
            row.addWidget(m)
            row.addStretch(1)
            f.addRow(label, row)
            self.limits[prefix] = (t, m)
        return box

    # ── binding ─────────────────────────────────────────────────────────
    def _study(self) -> Path | None:
        return self.window.model.path

    def options(self) -> launcher.PipelineOptions:
        return launcher.load_options(getattr(self.window, "session", None), self._study())

    def refresh(self) -> None:
        """Show the open study's stored options."""
        self._setting = True
        opts = self.options()
        study = self._study()
        self.setEnabled(study is not None)
        self.chain.clear()
        for p in opts.sequential_configs:
            self.chain.addItem(Path(p).name)
            item = self.chain.item(self.chain.count() - 1)
            item.setToolTip(p)
            item.setData(Qt.ItemDataRole.UserRole, p)
        for prefix, (t, m) in self.limits.items():
            t.setText(str(opts.stages.get(f"{prefix}_timelimit", "")))
            mem = opts.stages.get(f"{prefix}_memorylimit")
            m.setValue(int(mem) if isinstance(mem, (int, float)) else 0)
        self.summary.setChecked(opts.generate_summary_report is not False)
        model = self.window.model
        self.multi_box.setVisible(is_multi_isppa(model) or any(k.startswith(tuple(p for p, _ in MULTI_STAGES))
                                                             for k in opts.stages))
        self.unc_box.setVisible(model.get("simulation.uncertainty") is True or any(
            k.startswith(tuple(p for p, _ in UNCERTAINTY_STAGES)) for k in opts.stages))
        missing = [p for p in opts.sequential_configs if not Path(p).is_file()]
        self.chain_note.setText(
            ("Missing: " + ", ".join(Path(p).name for p in missing) + ". " if missing else "")
            + ("Chained runs need modules.run_heating_sims: 1 (the end state comes from the thermal stage)."
               if opts.sequential_configs and model.get("modules.run_heating_sims") != 1 else ""))
        self.toggle.setText("Pipeline options" + (" (set)" if not opts.is_empty() else ""))
        self._buttons()
        self._setting = False

    def showEvent(self, event) -> None:          # noqa: N802
        super().showEvent(event)
        self.refresh()

    def _toggled(self, on: bool) -> None:
        self.body.setVisible(on)

    def _buttons(self) -> None:
        r, n = self.chain.currentRow(), self.chain.count()
        self.remove_btn.setEnabled(r >= 0)
        self.up_btn.setEnabled(r > 0)
        self.down_btn.setEnabled(0 <= r < n - 1)

    def _chain_paths(self) -> list[str]:
        return [self.chain.item(i).data(Qt.ItemDataRole.UserRole) for i in range(self.chain.count())]

    def _store(self, *_args) -> None:
        if self._setting or self._study() is None:
            return
        stages: dict = {}
        for prefix, (t, m) in self.limits.items():
            text = t.text().strip()
            if text and _TIME.match(text):
                stages[f"{prefix}_timelimit"] = text
                t.setStyleSheet("")
            elif text:
                from . import style
                t.setStyleSheet(f"color: {style.active()['red']};")
                t.setToolTip("Use HH:MM:SS, e.g. 06:00:00")
            if m.value() > 0:
                stages[f"{prefix}_memorylimit"] = m.value()
        opts = launcher.PipelineOptions(self._chain_paths(), stages,
                                        None if self.summary.isChecked() else False)
        launcher.save_options(self.window.session, self._study(), opts)
        self.window.session.save()
        self.toggle.setText("Pipeline options" + (" (set)" if not opts.is_empty() else ""))
        self.changed.emit()

    def _add(self) -> None:
        start = str(self._study().parent) if self._study() else ""
        files, _ = QFileDialog.getOpenFileNames(self, "Add follow-up configs", start, "YAML (*.yaml *.yml)")
        study = self._study().resolve() if self._study() else None
        have = set(self._chain_paths())
        for f in files:
            p = str(Path(f).resolve())
            if p in have or (study is not None and Path(p) == study):
                continue
            self.chain.addItem(Path(p).name)
            item = self.chain.item(self.chain.count() - 1)
            item.setToolTip(p)
            item.setData(Qt.ItemDataRole.UserRole, p)
        self._store()
        self.refresh()

    def _remove(self) -> None:
        r = self.chain.currentRow()
        if r >= 0:
            self.chain.takeItem(r)
            self._store()
            self.refresh()

    def _move(self, step: int) -> None:
        r = self.chain.currentRow()
        if 0 <= r + step < self.chain.count():
            item = self.chain.takeItem(r)
            self.chain.insertItem(r + step, item)
            self._store()
            self.refresh()
            self.chain.setCurrentRow(r + step)
