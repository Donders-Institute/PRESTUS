"""'Follow-up from this run…': create a new study config that starts from an earlier run's
results (followup.py). The dialog shows the parent run, a type chooser with a one-line
explanation, the fields of the chosen type and a live list of problems; errors block
Create. Create saves the child next to the parent and emits ``created(path)``."""

from __future__ import annotations

from html import escape
from pathlib import Path

from PyQt6.QtCore import Qt, QTimer, pyqtSignal
from PyQt6.QtWidgets import (QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFileDialog,
                             QFormLayout, QHBoxLayout, QLabel, QLineEdit, QMessageBox, QPushButton,
                             QStackedWidget, QVBoxLayout, QWidget)

from .. import followup as fu
from .. import launcher
from ..model import Issue
from . import style
from .pose_table import PoseTable

# timing.* keys of config_default.yaml (#------------ Timing) shown for 'new_timing', in [s]
TIMING_KEYS = (("pd", "Pulse duration"), ("pri", "Pulse repetition interval"),
               ("ptd", "Pulse train duration"), ("ptri", "Pulse train repetition interval"),
               ("ptrd", "Pulse train repetition duration"), ("post_ptri_dur", "Cooling after the trains"))

AIMS = (("corrected", "Aberration-corrected aim point (recommended)"),
        ("geometric", "Geometric aim at the target"),
        ("apply_phases", "Keep the pose, apply the time-reversal phases"))


def _spin(value, maximum: float = 1e5, decimals: int = 4, suffix: str = "") -> QDoubleSpinBox:
    s = QDoubleSpinBox()
    s.setDecimals(decimals)
    s.setRange(0.0, maximum)
    s.setSuffix(suffix)
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        s.setValue(float(value))
    return s


class FollowupDialog(QDialog):
    created = pyqtSignal(object)        # Path of the saved follow-up config

    def __init__(self, window, run, ftype: str | None = None, pose: dict | None = None,
                 parent: QWidget | None = None):
        super().__init__(parent or window)
        self.window = window
        self.parent_model = window.model
        self.run = run
        self.child = None
        self.issues: list = []
        self._affix_edited = False
        self._file_edited = False
        self.setWindowTitle("Follow-up from this run")
        self.setMinimumWidth(620)
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(150)
        self._timer.timeout.connect(self._update)

        v = QVBoxLayout(self)
        cfg = self.parent_model.path
        head = QLabel(f"From run <b>{escape(run.label or run.affix or '(no affix)')}</b> "
                      f"(affix <code>{escape(run.affix or '—')}</code>) of "
                      f"<b>{escape(cfg.name if cfg else 'an unsaved config')}</b>")
        head.setWordWrap(True)
        v.addWidget(head)

        row = QHBoxLayout()
        row.addWidget(QLabel("Follow-up"))
        self.type_combo = QComboBox()
        for t in fu.TYPES.values():
            self.type_combo.addItem(t.label, t.key)
        row.addWidget(self.type_combo, 1)
        v.addLayout(row)
        self.blurb = QLabel()
        self.blurb.setObjectName("CardSub")
        self.blurb.setWordWrap(True)
        v.addWidget(self.blurb)

        self.pages = QStackedWidget()
        self._page = {}
        for key in fu.TYPES:
            w = getattr(self, f"_build_{key}")()
            self._page[key] = self.pages.addWidget(w)
        v.addWidget(self.pages, 1)

        form = QFormLayout()
        self.affix = QLineEdit()
        self.affix.setToolTip("io.output_affix of the follow-up; its outputs are named after it")
        self.affix.textEdited.connect(self._affix_typed)
        form.addRow("Output affix", self.affix)
        file_row = QHBoxLayout()
        self.file = QLineEdit()
        self.file.textEdited.connect(self._file_typed)
        file_row.addWidget(self.file, 1)
        browse = QPushButton("Save follow-up as…")
        browse.clicked.connect(self._browse)
        file_row.addWidget(browse)
        form.addRow("Config file", file_row)
        v.addLayout(form)

        self.problems = QLabel()
        self.problems.setWordWrap(True)
        self.problems.setTextFormat(Qt.TextFormat.RichText)
        v.addWidget(self.problems)

        self.buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Cancel)
        self.create_btn = self.buttons.addButton("Create", QDialogButtonBox.ButtonRole.AcceptRole)
        self.create_btn.setProperty("primary", True)
        self.buttons.rejected.connect(self.reject)
        self.create_btn.clicked.connect(self.create)
        v.addWidget(self.buttons)

        self.type_combo.currentIndexChanged.connect(self._type_changed)
        if pose is not None:
            ftype = "pose_from_time_reversal"
        idx = self.type_combo.findData(ftype) if ftype else 0
        self.type_combo.setCurrentIndex(max(idx, 0))
        self._type_changed()
        if pose is not None:
            self._select_pose(pose)

    # ── pages ───────────────────────────────────────────────────────────
    def _parent_target(self) -> float | None:
        t = self.parent_model.get("transducer.target_isppa_wcm2")
        if isinstance(t, list):
            t = next((x for x in t if isinstance(x, (int, float))), None)
        return t if isinstance(t, (int, float)) and not isinstance(t, bool) else None

    def _build_new_intensity(self) -> QWidget:
        w = QWidget()
        f = QFormLayout(w)
        self.isppa = _spin(self._parent_target() or 10.0, 1000.0, 2, " W/cm²")
        self.isppa.setMinimum(0.01)
        self.isppa.setToolTip("transducer.target_isppa_wcm2: free-water I<sub>SPPA</sub> the pressure map is "
                              "rescaled to")
        self.isppa.valueChanged.connect(self._changed)
        f.addRow("New target I<sub>SPPA</sub>", self.isppa)
        return w

    def _build_new_timing(self) -> QWidget:
        w = QWidget()
        f = QFormLayout(w)
        self.timing = {}
        for key, label in TIMING_KEYS:
            s = _spin(self.parent_model.get(f"timing.{key}"), suffix=" s")
            s.setToolTip(f"timing.{key}")
            s.valueChanged.connect(self._changed)
            self.timing[key] = s
            f.addRow(label, s)
        self.timing_isppa = _spin(self._parent_target() or 10.0, 1000.0, 2, " W/cm²")
        self.timing_isppa.valueChanged.connect(self._changed)
        f.addRow("Target I<sub>SPPA</sub>", self.timing_isppa)
        return w

    def _note_page(self, text: str) -> QWidget:
        w = QWidget()
        lv = QVBoxLayout(w)
        lab = QLabel(text)
        lab.setWordWrap(True)
        lv.addWidget(lab)
        lv.addStretch(1)
        return w

    def _build_new_placement_same_head(self) -> QWidget:
        return self._note_page(
            "The follow-up copies this config with headmodel.reuse_across_placements on. "
            "After Create it opens in the window: set the new position on the Placement tab.")

    def _build_sequential(self) -> QWidget:
        return self._note_page(
            "The follow-up is added to this config's sequential chain (Run tab ▸ Pipeline "
            "options). Running this config then runs the follow-up after it, starting from its end "
            "temperature and thermal dose. Change the follow-up's target or timing after Create.")

    def _build_pose_from_time_reversal(self) -> QWidget:
        w = QWidget()
        lv = QVBoxLayout(w)
        lv.setContentsMargins(0, 0, 0, 0)
        row = QHBoxLayout()
        row.addWidget(QLabel("Aim"))
        self.aim = QComboBox()
        for key, label in AIMS:
            self.aim.addItem(label, key)
        self.aim.currentIndexChanged.connect(self._aim_changed)
        row.addWidget(self.aim, 1)
        lv.addLayout(row)
        self.poses = PoseTable(show_use=False)
        self.poses.selection_changed.connect(lambda _r: self._changed())
        tables = fu.pose_tables(self.run)
        self.poses.load(tables)
        if not tables:
            self.poses.note.setText("This run has no pose table: run it with "
                                    "time_reversal.placement_search.enabled, or apply the phases.")
        lv.addWidget(self.poses, 1)
        return w

    def _aim_changed(self) -> None:
        self.poses.setEnabled(self.aim.currentData() != "apply_phases")
        self._changed()

    def _select_pose(self, pose: dict) -> None:
        for r, row in enumerate(self.poses.rows):
            if row.get("pose") == pose.get("pose"):
                self.poses.table.selectRow(r)
                break

    # ── state ───────────────────────────────────────────────────────────
    def ftype(self) -> str:
        return self.type_combo.currentData()

    def params(self) -> dict:
        t = self.ftype()
        p: dict = {}
        if t == "new_intensity":
            p["target_isppa_wcm2"] = self.isppa.value()
        elif t == "new_timing":
            p["timing"] = {k: s.value() for k, s in self.timing.items()
                           if s.value() != (self.parent_model.get(f"timing.{k}") or 0)}
            if self.timing_isppa.value() != (self._parent_target() or 0):
                p["target_isppa_wcm2"] = self.timing_isppa.value()
        elif t == "pose_from_time_reversal":
            p["aim"] = self.aim.currentData()
            if p["aim"] != "apply_phases":
                p["pose"] = self.poses.selected()
        elif t == "sequential":
            p["seq_index"] = len(self._chain()) + 2
        if self._affix_edited:
            p["affix"] = self.affix.text().strip()
        return p

    def _chain(self) -> list[str]:
        return launcher.load_options(getattr(self.window, "session", None),
                                     self.parent_model.path).sequential_configs

    def _type_changed(self) -> None:
        t = fu.TYPES[self.ftype()]
        self.blurb.setText(t.blurb)
        self.pages.setCurrentIndex(self._page[t.key])
        self._affix_edited = False
        self._changed()

    def _affix_typed(self, _text: str) -> None:
        self._affix_edited = True
        self._changed()

    def _file_typed(self, _text: str) -> None:
        self._file_edited = True

    def _changed(self, *_args) -> None:
        self._timer.start()

    def _default_file(self) -> Path:
        folder = self.parent_model.path.parent if self.parent_model.path else Path.home()
        stem = self.parent_model.path if self.parent_model.path else Path("config.yaml")
        return folder / fu.followup_filename(stem, self.ftype(), self.child.get("io.output_affix") or "")

    def _update(self) -> None:
        """Rebuild the child from the fields and show its problems (errors block Create)."""
        try:
            self.child = fu.create_followup(self.parent_model, self.run, self.ftype(), self.params())
            self.issues = fu.validate_followup(self.parent_model, self.run, self.ftype(), self.child,
                                               self._chain() if self.ftype() == "sequential" else None)
        except (KeyError, ValueError, TypeError) as exc:
            self.child = None
            self.issues = [Issue("error", "", str(exc))]
        if self.child is not None and not self._affix_edited:
            self.affix.setText(str(self.child.get("io.output_affix") or ""))
        if self.child is not None and not self._file_edited:
            self.file.setText(str(self._default_file()))
        self._show_problems()

    def _show_problems(self) -> None:
        c = style.active()
        colour = {"error": c["red"], "warning": c["amber"], "info": c["ink3"]}
        glyph = {"error": "✕", "warning": "!", "info": "i"}
        errors = [i for i in self.issues if i.level == "error"]
        if not self.issues:
            self.problems.setText(f"<span style='color:{c['green']}'>✓ Ready to create.</span>")
        else:
            items = "".join(
                f"<li><span style='color:{colour.get(i.level, c['ink3'])}'>{glyph.get(i.level, '')} "
                f"{escape(i.level)}</span>: {escape(i.message)}"
                + (f" <span style='color:{c['ink3']}'>({escape(i.path)})</span>" if i.path else "")
                + "</li>" for i in self.issues)
            self.problems.setText(f"<ul style='margin:4px 0'>{items}</ul>")
        self.create_btn.setEnabled(self.child is not None and not errors)
        self.create_btn.setToolTip("Fix the errors first" if errors else
                                   "Save the follow-up config and open it")

    def _browse(self) -> None:
        start = self.file.text() or str(self._default_file())
        chosen, _ = QFileDialog.getSaveFileName(self, "Save follow-up as", start, "YAML (*.yaml)")
        if chosen:
            self.file.setText(chosen)
            self._file_edited = True

    # ── create ──────────────────────────────────────────────────────────
    def create(self) -> Path | None:
        self._timer.stop()
        self._update()
        if self.child is None or any(i.level == "error" for i in self.issues):
            return None
        path = Path(self.file.text().strip() or self._default_file())
        if path.suffix.lower() not in (".yaml", ".yml"):
            path = path.with_suffix(".yaml")
        if path.exists() and QMessageBox.question(
                self, "Replace file?", f"{path.name} exists. Replace it?") != QMessageBox.StandardButton.Yes:
            return None
        if self.parent_model.path is not None and path.resolve() == Path(self.parent_model.path).resolve():
            QMessageBox.warning(self, "Same file", "Choose another file name than the parent config's.")
            return None
        from .main_window import study_header
        try:
            fu.save_followup(self.child, path, study_header(path))
        except OSError as exc:
            QMessageBox.warning(self, "Not saved", f"Could not save {path}:\n{exc}")
            return None
        if self.ftype() == "sequential" and getattr(self.window, "session", None) is not None:
            launcher.add_to_chain(self.window.session, self.parent_model.path, path)
            self.window.session.save()
        self.created.emit(path)
        self.accept()
        return path
