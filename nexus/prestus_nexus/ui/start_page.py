"""QuickStart tab: a first simulation in four steps (head image, transducer, target, run).

Everything a first run needs fits on one page; every other setting keeps its PRESTUS default
and stays on the other tabs. The page writes the same config keys as those tabs (it is a
view on the model, not a separate store), so switching between them never loses anything.
"""

from __future__ import annotations

from pathlib import Path

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (QButtonGroup, QComboBox, QDoubleSpinBox, QFileDialog, QFrame, QGridLayout,
                             QHBoxLayout, QLabel, QLineEdit, QPushButton, QRadioButton, QScrollArea,
                             QSpinBox, QVBoxLayout, QWidget)

from .. import docs, labels
from .. import quickstart as qs
from ..model import TX, ConfigModel, is_empty
from ..paths import m2m_dir, resolve_pattern
from . import icons, style

CONTENT_WIDTH = 920
# Which issues belong to which step (prefixes of Issue.path)
STEP_KEYS = {
    "image": ("subject_id", "path.", "startup.simnibs_bin_path"),
    "transducer": (f"{TX}.", "calibration."),
    "target": ("placement.", "placement viewer"),     # the viewer's checks (QC_PATH)
}
# 'Learn more' on each step: the page that explains it (step number -> doc reference)
STEP_DOCS = {1: "doc_nexus_tabs.md#quickstart-tab",
             2: "doc_nexus.md#2-check-the-transducer",
             3: "doc_placement.md#choosing-a-placement-mode",
             4: "doc_nexus.md#5-run-the-simulation"}
TARGET_MODES = {"heuristic": "placement.heuristic.mni_target_mm", "mni": "placement.mni.focus_pos_mm"}


class _Step(QFrame):
    """A numbered card with a title, a one-line status on the right, and a body."""

    def __init__(self, number: int, title: str, subtitle: str):
        super().__init__()
        self.setObjectName("Card")
        v = QVBoxLayout(self)
        v.setContentsMargins(18, 14, 18, 16)
        v.setSpacing(8)
        head = QHBoxLayout()
        self.badge = QLabel(str(number))
        self.badge.setFixedSize(26, 26)
        self.badge.setAlignment(Qt.AlignmentFlag.AlignCenter)
        head.addWidget(self.badge)
        t = QLabel(title)
        t.setObjectName("CardTitle")
        head.addWidget(t)
        head.addStretch(1)
        self.learn = QLabel()                   # set_doc: 'Learn more' link to the docs
        self.learn.setVisible(False)
        self.learn.setToolTip("Open documentation")
        self.learn.linkActivated.connect(docs.open_doc)
        head.addWidget(self.learn)
        self.status = QLabel()
        head.addWidget(self.status)
        v.addLayout(head)
        sub = QLabel(subtitle)
        sub.setObjectName("CardSub")
        sub.setWordWrap(True)
        v.addWidget(sub)
        self.body = QVBoxLayout()
        self.body.setSpacing(8)
        v.addLayout(self.body)
        self.set_state("todo", "To do")

    def set_doc(self, ref: str) -> None:
        """Show a 'Learn more' link that opens doc reference `ref`."""
        self.learn.setText(f"<a href='{ref}'>Learn more</a>")
        self.learn.setVisible(bool(ref))

    def set_state(self, state: str, text: str) -> None:
        c = style.active()
        colour = {"done": c["green"], "todo": c["ink3"], "warn": c["amber"], "error": c["red"]}[state]
        mark = {"done": "✓", "todo": "○", "warn": "⚠", "error": "✗"}[state]
        fill = c["green"] if state == "done" else c["brand"]
        ink = "#ffffff" if state == "done" else c["on_brand"]
        self.badge.setStyleSheet(f"background: {fill}; color: {ink}; border-radius: 13px; font-weight: 600;")
        self.badge.setText("✓" if state == "done" else self.badge.property("n") or self.badge.text())
        self.status.setText(f"<span style='color:{colour}'>{mark} {text}</span>")


def _value_label() -> QLabel:
    lbl = QLabel()
    lbl.setWordWrap(True)
    lbl.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
    return lbl


def _coord_spins() -> list[QDoubleSpinBox]:
    spins = []
    for axis in "xyz":
        s = QDoubleSpinBox()
        s.setRange(-200, 200)
        s.setDecimals(1)
        s.setSuffix(" mm")
        s.setPrefix(f"{axis}  ")
        s.setKeyboardTracking(False)          # commit on Enter / focus out, not per keystroke
        s.setMinimumWidth(110)
        spins.append(s)
    return spins


class StartPage(QWidget):
    """The four steps. `window` is the MainWindow (model, issues, navigation, save)."""

    goto_tab = pyqtSignal(str)
    key = "start"                     # like TabPage.key: the tab's name for navigation and tests

    def __init__(self, window, parent: QWidget | None = None):
        super().__init__(parent)
        self.window = window
        self._setting = False
        self.setAcceptDrops(True)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        outer.addWidget(scroll)
        host = QWidget()
        hl = QHBoxLayout(host)
        hl.setContentsMargins(16, 16, 16, 24)
        col = QWidget()
        col.setMaximumWidth(CONTENT_WIDTH)
        hl.addStretch(1)
        hl.addWidget(col, 100)
        hl.addStretch(1)
        scroll.setWidget(host)
        self.col = QVBoxLayout(col)
        self.col.setSpacing(14)

        title = QLabel("Set up a simulation")
        title.setStyleSheet("font-size: 20px; font-weight: 600;")
        self.col.addWidget(title)
        intro = QLabel("Four steps. Everything else uses the PRESTUS defaults, which suit most "
                       "studies; the other tabs hold every setting if you need to change one.")
        intro.setObjectName("TabIntro")
        intro.setWordWrap(True)
        self.col.addWidget(intro)

        self.step_image = self._build_image_step()
        self.step_tx = self._build_transducer_step()
        self.step_target = self._build_target_step()
        self.step_run = self._build_run_step()
        for i, s in enumerate((self.step_image, self.step_tx, self.step_target, self.step_run), 1):
            s.badge.setProperty("n", str(i))
            s.set_doc(STEP_DOCS.get(i, ""))
            self.col.addWidget(s)
        self.col.addStretch(1)

    @property
    def model(self) -> ConfigModel:
        return self.window.model

    # ── step 1: head image ──────────────────────────────────────────────
    def _build_image_step(self) -> _Step:
        s = _Step(1, "Head image", "The participant's T1-weighted MRI (NIfTI, .nii or .nii.gz). "
                  "PRESTUS finds the T2, the SimNIBS segmentation and the subject number from it. "
                  "You can also drop the file onto this window.")
        row = QHBoxLayout()
        self.t1_edit = QLineEdit(readOnly=True)
        self.t1_edit.setPlaceholderText("No image chosen")
        row.addWidget(self.t1_edit, 1)
        self.t1_btn = QPushButton("Choose T1 image…")
        self.t1_btn.clicked.connect(self.choose_t1)
        row.addWidget(self.t1_btn)
        s.body.addLayout(row)

        self.image_details = QWidget()
        g = QGridLayout(self.image_details)
        g.setContentsMargins(0, 4, 0, 0)
        g.setHorizontalSpacing(14)
        g.setVerticalSpacing(6)
        g.setColumnStretch(1, 1)
        self.sid = QSpinBox()
        self.sid.setRange(0, 999)
        self.sid.setKeyboardTracking(False)
        self.sid.setToolTip("subject_id: names the outputs (sub-XXX) and the segmentation folder (m2m_sub-XXX)")
        self.sid.valueChanged.connect(self._sid_changed)
        self.t2_lbl = _value_label()
        self.seg_lbl = _value_label()
        self.simnibs_btn = QPushButton("Choose SimNIBS folder…")
        self.simnibs_btn.setToolTip("The SimNIBS 'bin' folder, which holds charm (startup.simnibs_bin_path)")
        self.simnibs_btn.clicked.connect(self.choose_simnibs)
        self.out_lbl = _value_label()
        out_btn = QPushButton("Change…")
        out_btn.setToolTip("Where PRESTUS writes the results (path.sim)")
        out_btn.clicked.connect(self.choose_output)
        rows = [("Subject number", self.sid, None), ("T2 image", self.t2_lbl, None),
                ("Segmentation", self.seg_lbl, self.simnibs_btn), ("Results folder", self.out_lbl, out_btn)]
        for r, (name, w, btn) in enumerate(rows):
            lab = QLabel(name)
            lab.setObjectName("CardSub")
            g.addWidget(lab, r, 0, Qt.AlignmentFlag.AlignTop)
            g.addWidget(w, r, 1, Qt.AlignmentFlag.AlignTop | (Qt.AlignmentFlag.AlignLeft if w is self.sid
                                                               else Qt.AlignmentFlag(0)))
            if btn is not None:
                g.addWidget(btn, r, 2, Qt.AlignmentFlag.AlignTop)
        s.body.addWidget(self.image_details)
        return s

    def dragEnterEvent(self, event) -> None:      # noqa: N802 (Qt API)
        if self._dropped_nifti(event) is not None:
            event.acceptProposedAction()

    def dropEvent(self, event) -> None:           # noqa: N802 (Qt API)
        path = self._dropped_nifti(event)
        if path is not None:
            self.apply_t1(path)
            event.acceptProposedAction()

    @staticmethod
    def _dropped_nifti(event) -> Path | None:
        urls = event.mimeData().urls() if event.mimeData().hasUrls() else []
        for u in urls:
            p = Path(u.toLocalFile())
            if p.is_file() and (p.name.endswith(".nii") or p.name.endswith(".nii.gz")):
                return p
        return None

    def choose_t1(self) -> None:
        anat = self.model.get("path.anat")
        start = anat if isinstance(anat, str) and anat and Path(anat).is_dir() else str(Path.home())
        chosen, _ = QFileDialog.getOpenFileName(self, "Choose the T1-weighted image", start,
                                                f"{qs.NIFTI_FILTER};;All files (*)")
        if chosen:
            self.apply_t1(chosen)

    def apply_t1(self, path: str | Path) -> None:
        """Set path.anat, t1/t2 patterns, subject_id and (if found) path.seg from one file,
        then fill in an output folder and SimNIBS if they are still empty."""
        m = self.model
        sid = m.get("subject_id")
        found = qs.from_t1(path, sid if isinstance(sid, int) else None)
        for key, value in found.as_settings().items():
            if key == "path.seg" and not found.m2m_found:
                continue                          # keep a segmentation folder set by hand
            m.set(key, value)
        if is_empty(m.get("path.sim")):
            m.set("path.sim", str(qs.default_output(found.anat, m.path)))
        if not found.m2m_found and is_empty(m.get("startup.simnibs_bin_path")):
            simnibs = qs.find_simnibs_bin()
            if simnibs is not None:
                m.set("startup.simnibs_bin_path", str(simnibs))
        self.window.on_model_edited()

    def _sid_changed(self, value: int) -> None:
        if not self._setting:
            self.model.set("subject_id", value)
            self.window.on_model_edited()

    def choose_simnibs(self) -> None:
        chosen = QFileDialog.getExistingDirectory(self, "Choose the SimNIBS bin folder (holds charm)",
                                                  str(qs.find_simnibs_bin() or Path.home()))
        if chosen:
            self.model.set("startup.simnibs_bin_path", chosen)
            self.window.on_model_edited()

    def choose_output(self) -> None:
        cur = self.model.get("path.sim")
        chosen = QFileDialog.getExistingDirectory(self, "Choose the results folder",
                                                  cur if isinstance(cur, str) and cur else str(Path.home()))
        if chosen:
            self.model.set("path.sim", chosen)
            self.window.on_model_edited()

    def _refresh_image(self) -> None:
        m, c = self.model, style.active()
        anat, pattern, sid = m.get("path.anat"), m.get("path.t1_pattern") or "", m.get("subject_id")
        sid = sid if isinstance(sid, int) else 0
        hits = resolve_pattern(anat, pattern, sid) if isinstance(anat, str) and anat else []
        self.t1_edit.setText(str(hits[0]) if hits else "")
        self.t1_btn.setText("Change…" if hits else "Choose T1 image…")
        self.t1_btn.setProperty("primary", not hits)
        self.t1_btn.style().unpolish(self.t1_btn)
        self.t1_btn.style().polish(self.t1_btn)
        self.image_details.setVisible(bool(hits) or not is_empty(anat))
        self.sid.blockSignals(True)
        self.sid.setValue(sid)
        self.sid.blockSignals(False)
        t2 = m.get("path.t2_pattern")
        t2_hits = resolve_pattern(anat, t2, sid) if not is_empty(t2) and isinstance(anat, str) else []
        self.t2_lbl.setText(t2_hits[0].name if t2_hits else
                            "None (optional: it improves the skull segmentation)")
        seg = m.get("path.seg")
        root = anat if is_empty(seg) else seg
        m2m = m2m_dir(root, sid) if isinstance(root, str) and root else None
        simnibs = m.get("startup.simnibs_bin_path")
        has_simnibs = not is_empty(simnibs) and qs.is_simnibs_bin(simnibs)
        if m2m is not None and m2m.is_dir():
            self.seg_lbl.setText(f"Found: {m2m}")
            self.simnibs_btn.setVisible(False)
        elif has_simnibs:
            self.seg_lbl.setText(f"Not found yet. SimNIBS ({simnibs}) creates <i>{m2m.name if m2m else 'it'}</i> "
                                 "on the first run, which adds a few hours once per participant.")
            self.simnibs_btn.setText("Change SimNIBS…")
            self.simnibs_btn.setVisible(True)
        else:
            self.seg_lbl.setText(f"<span style='color:{c['red']}'>Not found, and SimNIBS is not set.</span> "
                                 "Choose the SimNIBS folder so the segmentation can be created.")
            self.simnibs_btn.setText("Choose SimNIBS folder…")
            self.simnibs_btn.setVisible(True)
        sim = m.get("path.sim")
        self.out_lbl.setText(str(sim) if not is_empty(sim) else "Not set")

    # ── step 2: transducer ──────────────────────────────────────────────
    def _build_transducer_step(self) -> _Step:
        s = _Step(2, "Transducer", "The device and the driving system it is used with, from the equipment "
                  "library. Its calibration sets the phases and amplitude.")
        g = QGridLayout()
        g.setHorizontalSpacing(14)
        g.setColumnStretch(1, 1)
        self.dev_combo = QComboBox()
        self.dev_combo.currentIndexChanged.connect(self._device_chosen)
        self.ds_combo = QComboBox()
        self.ds_combo.currentIndexChanged.connect(self._ds_chosen)
        self.isppa = QDoubleSpinBox()
        self.isppa.setRange(0.1, 200)
        self.isppa.setDecimals(1)
        self.isppa.setSuffix(" W/cm²")
        self.isppa.setKeyboardTracking(False)
        self.isppa.setToolTip("transducer.target_isppa_wcm2: the free-water spatial-peak pulse-average "
                              "intensity the results are scaled to")
        self.isppa.valueChanged.connect(self._isppa_changed)
        self.ep = QDoubleSpinBox()
        self.ep.setRange(0, 200)
        self.ep.setDecimals(1)
        self.ep.setSuffix(" mm")
        self.ep.setSpecialValueText("Not set")
        self.ep.setKeyboardTracking(False)
        self.ep.setToolTip("transducer.focal_distance_ep: the depth the device is set to focus at, "
                           "measured from its exit plane. The calibration is looked up at this depth.")
        self.ep.valueChanged.connect(self._ep_changed)
        self.ep_hint = QLabel()
        self.ep_hint.setObjectName("CardSub")
        ep_row = QHBoxLayout()
        ep_row.setSpacing(10)
        ep_row.addWidget(self.ep)
        ep_row.addWidget(self.ep_hint, 1)
        rows = [("Device", self.dev_combo), ("Driving system", self.ds_combo),
                ("Focal depth (from exit plane)", ep_row), ("Free-water intensity (I_SPPA)", self.isppa)]
        for r, (name, w) in enumerate(rows):
            lab = QLabel(name)
            lab.setObjectName("CardSub")
            g.addWidget(lab, r, 0)
            if isinstance(w, QHBoxLayout):
                g.addLayout(w, r, 1)
            else:
                g.addWidget(w, r, 1, Qt.AlignmentFlag.AlignLeft if w is self.isppa else Qt.AlignmentFlag(0))
        s.body.addLayout(g)
        self.cal_lbl = _value_label()
        s.body.addWidget(self.cal_lbl)
        self.tx_multi = _value_label()
        s.body.addWidget(self.tx_multi)
        return s

    def _device_chosen(self, index: int) -> None:
        if self._setting or index < 0:
            return
        serial = self.dev_combo.itemData(index)
        m = self.model
        m.set(f"{TX}.serial", serial or "")
        dev = m.device(0)
        combos = list(dev.combos) if dev is not None else []
        if m.get(f"{TX}.combo.ds_serial") not in combos:
            m.set(f"{TX}.combo.ds_serial", combos[0] if len(combos) == 1 else "")
        self.window.on_model_edited()

    def _ds_chosen(self, index: int) -> None:
        if self._setting or index < 0:
            return
        self.model.set(f"{TX}.combo.ds_serial", self.ds_combo.itemData(index) or "")
        self.window.on_model_edited()

    def _ep_changed(self, value: float) -> None:
        if not self._setting:
            self.model.set(f"{TX}.focal_distance_ep", round(value, 1) if value > 0 else [])
            self.window.on_model_edited()

    def _isppa_changed(self, value: float) -> None:
        if not self._setting:
            self.model.set(f"{TX}.target_isppa_wcm2", round(value, 1))
            self.window.on_model_edited()

    def _refresh_transducer(self) -> None:
        m = self.model
        lib = m.library
        devices = lib.devices if lib is not None else {}
        serial = m.get(f"{TX}.serial")
        self.dev_combo.clear()
        self.dev_combo.addItem("Choose a device…", "")
        for key, dev in devices.items():
            self.dev_combo.addItem(dev.label(), key)
        idx = self.dev_combo.findData(str(serial)) if not is_empty(serial) else 0
        if idx < 0:
            self.dev_combo.addItem(f"{serial} (not in the library)", str(serial))
            idx = self.dev_combo.count() - 1
        self.dev_combo.setCurrentIndex(idx)
        dev = m.device(0)
        ds = m.get(f"{TX}.combo.ds_serial")
        self.ds_combo.clear()
        self.ds_combo.addItem("None (generic calibration)", "")
        for key, d in (dev.combos.items() if dev is not None else []):
            self.ds_combo.addItem(f"{key} — {d.name}" if d.name else key, key)
        di = self.ds_combo.findData(str(ds)) if not is_empty(ds) else 0
        self.ds_combo.setCurrentIndex(max(di, 0))
        self.ds_combo.setEnabled(dev is not None)
        ep = m.get(f"{TX}.focal_distance_ep")
        self.ep.setValue(float(ep) if isinstance(ep, (int, float)) else 0.0)
        if dev is not None and dev.min_foc is not None and dev.max_foc is not None:
            inside = not isinstance(ep, (int, float)) or dev.min_foc <= ep <= dev.max_foc
            c = style.active()
            rng = f"device range {dev.min_foc:g}–{dev.max_foc:g} mm"
            self.ep_hint.setText(rng if inside else f"<span style='color:{c['amber']}'>outside the {rng}</span>")
        else:
            self.ep_hint.setText("")
        isppa = m.get(f"{TX}.target_isppa_wcm2")
        self.isppa.setEnabled(isinstance(isppa, (int, float)))
        if isinstance(isppa, (int, float)):
            self.isppa.setValue(float(isppa))
        else:
            self.isppa.setToolTip(f"Several intensities are set ({isppa}); edit them on the Transducer tab.")
        if dev is not None and lib is not None:
            cal = lib.calibration(dev.serial, None if is_empty(ds) else str(ds),
                                  m.get(f"{TX}.library_path") or None)
            c = style.active()
            colour = c["green"] if cal.found and not cal.fallback else c["amber"] if cal.found else c["red"]
            self.cal_lbl.setText(f"<span style='color:{colour}'>Calibration:</span> {cal.summary()}")
        else:
            self.cal_lbl.setText("")
        n = m.n_transducers
        self.tx_multi.setVisible(n > 1)
        self.tx_multi.setText(f"This study has {n} transducers; this step edits transducer 1. "
                              "The others are on the Transducer tab.")

    # ── step 3: target ──────────────────────────────────────────────────
    def _build_target_step(self) -> _Step:
        s = _Step(3, "Target", "Where to focus, in MNI coordinates (mm). PRESTUS maps them into the "
                  "participant's head with the SimNIBS segmentation.")
        self.rb_auto = QRadioButton("Target only: PRESTUS finds a transducer position")
        self.rb_auto.setToolTip("placement.mode: heuristic (searches the scalp for a good entry point)")
        self.rb_entry = QRadioButton("Target and transducer position")
        self.rb_entry.setToolTip("placement.mode: mni (the transducer point is snapped to the scalp)")
        node = self.model.schema.get("placement.mode")
        for b, mode in ((self.rb_auto, "heuristic"), (self.rb_entry, "mni")):
            ic = icons.for_choice(node, mode)
            if ic is not None:
                b.setIcon(ic)
        self.mode_group = QButtonGroup(self)
        for b in (self.rb_auto, self.rb_entry):
            self.mode_group.addButton(b)
            s.body.addWidget(b)
        self.mode_group.buttonClicked.connect(self._mode_chosen)
        self.other_mode = _value_label()
        s.body.addWidget(self.other_mode)

        g = QGridLayout()
        g.setHorizontalSpacing(10)
        self.target_name = QLineEdit()
        self.target_name.setPlaceholderText("e.g. left_amygdala")
        self.target_name.setToolTip("placement.heuristic.target_name: used in the output file names")
        self.target_name.setMaximumWidth(320)
        self.target_name.editingFinished.connect(self._name_changed)
        self.focus = _coord_spins()
        self.entry = _coord_spins()
        for sp in self.focus + self.entry:
            sp.valueChanged.connect(self._coords_changed)
        self.name_lbl = QLabel("Target name")
        self.entry_lbl = QLabel("Transducer (scalp)")
        for r, (lab, ws) in enumerate([(self.name_lbl, [self.target_name]), (QLabel("Target"), self.focus),
                                       (self.entry_lbl, self.entry)]):
            lab.setObjectName("CardSub")
            g.addWidget(lab, r, 0)
            for k, w in enumerate(ws):
                g.addWidget(w, r, 1 + k, 1, 3 if len(ws) == 1 else 1)
        g.setColumnStretch(4, 1)
        self.target_grid = QWidget()
        self.target_grid.setLayout(g)
        s.body.addWidget(self.target_grid)
        return s

    def _mode(self) -> str:
        return str(self.model.get("placement.mode") or "manual")

    def _mode_chosen(self, button) -> None:
        if self._setting:
            return
        m = self.model
        new = "heuristic" if button is self.rb_auto else "mni"
        old = self._mode()
        if old in TARGET_MODES and old != new:          # carry the target over
            m.set(TARGET_MODES[new], m.get(TARGET_MODES[old]))
        m.set("placement.mode", new)
        self.window.on_model_edited()

    def _name_changed(self) -> None:
        if not self._setting and self._mode() == "heuristic":
            self.model.set("placement.heuristic.target_name", self.target_name.text().strip())
            self.window.on_model_edited()

    def _coords_changed(self, _v: float = 0.0) -> None:
        if self._setting:
            return
        mode = self._mode()
        if mode not in TARGET_MODES:
            return
        m = self.model
        m.set(TARGET_MODES[mode], [round(s.value(), 1) for s in self.focus])
        if mode == "mni":
            m.set("placement.mni.trans_pos_mm", [round(s.value(), 1) for s in self.entry])
        self.window.on_model_edited()

    def _refresh_target(self) -> None:
        m = self.model
        mode = self._mode()
        known = mode in TARGET_MODES
        self.rb_auto.setChecked(mode == "heuristic")
        self.rb_entry.setChecked(mode == "mni")
        if not known:
            self.mode_group.setExclusive(False)
            self.rb_auto.setChecked(False)
            self.rb_entry.setChecked(False)
            self.mode_group.setExclusive(True)
        self.other_mode.setVisible(not known)
        self.other_mode.setText(f"This study places the transducer with <b>{mode}</b>; that is set on the "
                                "Placement tab. Choose an option above to use MNI coordinates instead.")
        self.target_grid.setVisible(known)
        self.name_lbl.setVisible(mode == "heuristic")
        self.target_name.setVisible(mode == "heuristic")
        self.entry_lbl.setVisible(mode == "mni")
        for sp in self.entry:
            sp.setVisible(mode == "mni")
        if not known:
            return
        if not self.target_name.hasFocus():
            self.target_name.setText(str(m.get("placement.heuristic.target_name") or ""))
        for spins, key in ((self.focus, TARGET_MODES[mode]), (self.entry, "placement.mni.trans_pos_mm")):
            vec = m.get(key)
            if isinstance(vec, list) and len(vec) == 3:
                for sp, val in zip(spins, vec):
                    if not sp.hasFocus():
                        sp.setValue(float(val))

    # ── step 4: run ─────────────────────────────────────────────────────
    def _build_run_step(self) -> _Step:
        s = _Step(4, "Check and run", "Save the study and run it. The first run on a participant also "
                  "creates the segmentation, so it takes longest.")
        where = QHBoxLayout()
        lab = QLabel("Run on")
        lab.setObjectName("CardSub")
        where.addWidget(lab)
        self.where = QComboBox()
        self.where.addItem("Automatic: the cluster if one is found, else this computer", "auto")
        self.where.addItem("This computer (MATLAB)", "matlab")
        self.where.addItem("The cluster, SLURM", "slurm")
        self.where.addItem("The cluster, qsub", "qsub")
        self.where.setToolTip("platform: 'auto' submits to SLURM/qsub when a scheduler is found and "
                              "runs here otherwise; 'matlab' always runs in MATLAB on this computer")
        self.where.currentIndexChanged.connect(self._where_chosen)
        where.addWidget(self.where)
        where.addStretch(1)
        s.body.addLayout(where)
        self.problems = QLabel()
        self.problems.setWordWrap(True)
        self.problems.setTextFormat(Qt.TextFormat.RichText)
        self.problems.linkActivated.connect(self._goto_problem)
        s.body.addWidget(self.problems)
        buttons = QHBoxLayout()
        more = QPushButton("All settings")
        more.setToolTip("Open the Setup tab; every other tab holds the remaining settings")
        more.clicked.connect(lambda: self.goto_tab.emit("setup"))
        buttons.addWidget(more)
        buttons.addStretch(1)
        self.save_btn = QPushButton("Save…")
        self.save_btn.clicked.connect(self.window.save)
        buttons.addWidget(self.save_btn)
        self.run_btn = QPushButton("Continue to Run ›")
        self.run_btn.setProperty("primary", True)
        self.run_btn.setToolTip("Save the study and open the Run tab, where you start the run")
        self.run_btn.clicked.connect(self._continue)
        buttons.addWidget(self.run_btn)
        s.body.addLayout(buttons)
        return s

    def _where_chosen(self, index: int) -> None:
        if not self._setting and index >= 0:
            self.model.set("platform", self.where.itemData(index))
            self.window.on_model_edited()

    def _goto_problem(self, link: str) -> None:
        issue = next((i for i in self.window.issues if f"{i.path}#{i.transducer}" == link), None)
        if issue is not None:
            self.window.goto_issue(issue)

    def _continue(self) -> None:
        if self.window.save():
            self.goto_tab.emit("run")

    def _refresh_run(self) -> None:
        platform = self.model.get("platform")
        self.where.setCurrentIndex(max(self.where.findData(platform), 0))
        errors = [i for i in self.window.issues if i.level == "error"]
        c = style.active()
        if errors:
            items = "".join(f"<li><a href='{i.path}#{i.transducer}' style='color:{c['red']}'>"
                            f"{self._friendly(i)}</a>: {i.message}</li>" for i in errors[:5])
            extra = f"<br>…and {len(errors) - 5} more under Issues in the header." if len(errors) > 5 else ""
            self.problems.setText(f"<b>{len(errors)} thing(s) to fix before running</b> "
                                  f"(click one to go to it):<ul style='margin:4px 0'>{items}</ul>{extra}")
        else:
            self.problems.setText(f"<span style='color:{c['green']}'>✓ Ready to run.</span>")
        self.run_btn.setEnabled(not errors)
        self.run_btn.setToolTip("Fix the items above first" if errors
                                else "Save the study and open the Run tab, where you start the run")
        self.save_btn.setText("Save…" if self.model.path is None else "Save")

    # ── all ─────────────────────────────────────────────────────────────
    def _step_issues(self, step: str) -> list:
        keys = STEP_KEYS[step]
        return [i for i in self.window.issues if i.level in ("error", "warning")
                and any(i.path == k or i.path.startswith(k) for k in keys)]

    def _set_step_state(self, s: _Step, step: str, complete: bool) -> None:
        """'To do' until the user has filled the step in (its errors are expected then);
        after that the count of problems, with the messages in the tooltip."""
        issues = self._step_issues(step)
        errs = [i for i in issues if i.level == "error"]
        warns = [i for i in issues if i.level == "warning"]
        s.status.setToolTip("\n".join(f"{self._friendly(i)}: {i.message}" for i in issues))
        if not complete:
            s.set_state("todo", "To do")
        elif errs:
            s.set_state("error", f"{len(errs)} problem{'s' if len(errs) > 1 else ''}")
        elif warns:
            s.set_state("warn", f"{len(warns)} warning{'s' if len(warns) > 1 else ''}")
        else:
            s.set_state("done", "Done")

    def _friendly(self, issue) -> str:
        """The field's label as the forms show it (the key stays in the link's tooltip)."""
        node = self.model.schema.get(issue.path)
        return labels.label(node) if node is not None else issue.where()

    def refresh(self) -> None:
        self._setting = True
        try:
            self._refresh_image()
            self._refresh_transducer()
            self._refresh_target()
            self._refresh_run()
        finally:
            self._setting = False
        m = self.model
        self._set_step_state(self.step_image, "image", bool(self.t1_edit.text()))
        self._set_step_state(self.step_tx, "transducer", not is_empty(m.get(f"{TX}.serial"))
                             and isinstance(m.get(f"{TX}.focal_distance_ep"), (int, float)))
        mode = self._mode()
        target = m.get(TARGET_MODES.get(mode, "")) if mode in TARGET_MODES else None
        if mode in TARGET_MODES:
            placed = isinstance(target, list) and len(target) == 3 and \
                (mode != "heuristic" or not is_empty(m.get("placement.heuristic.target_name")))
        elif mode == "manual":
            placed = all(not is_empty(m.get(f"{TX}.{k}")) for k in ("trans_pos", "focus_pos"))
        else:
            placed = True                       # localite/plantus/…: their checks are the issues
        self._set_step_state(self.step_target, "target", placed)
        n_err = sum(1 for i in self.window.issues if i.level == "error")
        self.step_run.set_state("done" if not n_err else "todo",
                                "Ready" if not n_err else f"{n_err} error(s) left")
