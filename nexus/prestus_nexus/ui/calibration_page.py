"""Transducer tab ▸ Calibration: compare, fit and save library entries for the study's
transducers (the steps of calibration_gui, computed in MATLAB through the bridge).

The transducer and driving system come from the config, the entry is saved to
transducer.library_path, and the desired depths are prefilled with the depths this study
uses: the focal depth settings and the placement's target depths ("Study depths").
"""

from __future__ import annotations

import datetime as _dt
import math
from pathlib import Path

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (QCheckBox, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox,
                             QFileDialog, QFormLayout, QHBoxLayout, QLabel, QLineEdit,
                             QListWidget, QListWidgetItem, QMessageBox, QPlainTextEdit, QPushButton,
                             QScrollArea, QSpinBox, QSplitter, QTableWidget, QTableWidgetItem,
                             QToolButton, QVBoxLayout, QWidget)

from .. import calibration as cal
from ..bridge import BridgeError, as_columns, as_floats, as_list
from ..model import TX, is_empty
from .matlab_status import MatlabStatus

try:
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
    from matplotlib.figure import Figure
except ImportError:                                    # pragma: no cover
    FigureCanvasQTAgg = None

# result kind, colour, legend name (result_styles in calibration_gui.m)
STYLES = [("geometric", "#999999", "Geometric"), ("library", "#0072b3", "Library"),
          ("fit", "#d95f02", "Fit"), ("kwave", "#33a02b", "k-Wave")]
PROFILE_LINES = [("measured", "-", "Hydrophone"), ("interpolated", "--", "Interpolated profile"),
                 ("extrapolated", ":", "Extrapolated profile")]


class AdvancedDialog(QDialog):
    """Advanced FIT_TRANSDUCER_MODEL options; only values that differ from the defaults are kept."""

    def __init__(self, current: dict, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Advanced fit options")
        cur = dict(cal.ADVANCED_DEFAULTS, **current)
        form = QFormLayout(self)
        self.pc = QComboBox()
        self.pc.addItems(["none", "linear", "monotonic"])
        self.pc.setCurrentText(str(cur["PhaseConstraint"]))
        form.addRow("Phase constraint", self.pc)
        self.sb = QLineEdit(_pair_text(cur["PhaseStepBounds"]))
        self.sb.setPlaceholderText("lo, hi (rad); blank = default")
        form.addRow("Phase step bounds", self.sb)
        self.lam = _spin(cur["Lambda"], 0, 100, 6)
        form.addRow("Phase prior λ", self.lam)
        self.lama = _spin(cur["LambdaAmp"], 0, 100, 6)
        form.addRow("Gain prior λ", self.lama)
        self.af = QComboBox()
        self.af.addItems(["shape", "absolute"])
        self.af.setCurrentText(str(cur["AmplitudeFit"]))
        form.addRow("Amplitude fit", self.af)
        self.sv = QCheckBox("one drive for all depths (absolute)")
        self.sv.setChecked(bool(cur["SharedVelocity"]))
        form.addRow("", self.sv)
        self.fw = _spin(cur["FocalWeight"], 0, 1000, 3)
        form.addRow("Focal weight (0 = uniform)", self.fw)
        self.win = QLineEdit(_pair_text(cur["FitWindowMm"]))
        self.win.setPlaceholderText("lo, hi (mm from EP); blank = all")
        form.addRow("Fit window", self.win)
        self.zb = QLineEdit(_pair_text(cur["ZShiftBounds"]))
        form.addRow("z shift bounds (mm)", self.zb)
        self.opt = QComboBox()
        self.opt.addItems(["lsqnonlin", "fexminimize", "globalsearch"])
        self.opt.setCurrentText(str(cur["Optimizer"]))
        form.addRow("Optimiser", self.opt)
        self.pol = QCheckBox("polish")
        self.pol.setChecked(bool(cur["Polish"]))
        form.addRow("", self.pol)
        self.pop = QSpinBox(minimum=10, maximum=100000, value=int(cur["PopSize"]))
        form.addRow("FEX population", self.pop)
        self.mfe = QSpinBox(minimum=100, maximum=100000000, value=int(cur["MaxFunEvals"]))
        form.addRow("Max. evaluations", self.mfe)
        box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
                               | QDialogButtonBox.StandardButton.Reset)
        box.accepted.connect(self.accept)
        box.rejected.connect(self.reject)
        box.button(QDialogButtonBox.StandardButton.Reset).clicked.connect(self._reset)
        form.addRow(box)
        self._reset_all = False

    def _reset(self) -> None:
        self._reset_all = True
        self.accept()

    def values(self) -> dict:
        if self._reset_all:
            return {}
        v = {"PhaseConstraint": self.pc.currentText(), "PhaseStepBounds": _pair(self.sb.text()),
             "AmplitudeFit": self.af.currentText(), "SharedVelocity": self.sv.isChecked(),
             "Optimizer": self.opt.currentText(), "PopSize": self.pop.value(),
             "MaxFunEvals": self.mfe.value(), "Polish": self.pol.isChecked(),
             "FocalWeight": self.fw.value(), "FitWindowMm": _pair(self.win.text()),
             "Lambda": self.lam.value(), "LambdaAmp": self.lama.value(),
             "ZShiftBounds": _pair(self.zb.text()) or [-3, 3]}
        return {k: x for k, x in v.items() if not cal._same(x, cal.ADVANCED_DEFAULTS[k])}


def _spin(value, lo, hi, decimals) -> QDoubleSpinBox:
    s = QDoubleSpinBox(decimals=decimals, minimum=lo, maximum=hi)
    s.setValue(float(value))
    return s


def _pair_text(v) -> str:
    v = as_list(v)
    return ", ".join(f"{float(x):g}" for x in v) if len(v) == 2 else ""


def _pair(text: str) -> list[float]:
    try:
        v = [float(x) for x in text.replace(";", ",").split(",") if x.strip()]
    except ValueError:
        return []
    return sorted(v) if len(v) == 2 else []


class CalibrationPage(QWidget):
    def __init__(self, window, parent: QWidget | None = None):
        super().__init__(parent)
        self.window = window
        self.session_id = f"nexus-{id(self):x}"
        self.opened: tuple | None = None       # (serial, ds, library_path) open in MATLAB
        self.info: dict = {}                   # last 'open' / 'load' answer
        self.view: dict | None = None          # last 'view' / 'fit' answer
        self.advanced: dict = {}
        self.auto_depths: str | None = None    # desired-depths text this page filled in itself
        self.has_fit = False
        self._busy = False

        outer = QVBoxLayout(self)
        outer.setContentsMargins(12, 8, 12, 8)
        top = QHBoxLayout()
        top.addWidget(QLabel("Transducer"))
        self.tx = QComboBox()
        self.tx.setMinimumWidth(340)
        self.tx.currentIndexChanged.connect(lambda _i: self._open())
        top.addWidget(self.tx)
        self.coverage = QLabel()
        self.coverage.setObjectName("CardSub")
        self.coverage.setWordWrap(True)
        top.addWidget(self.coverage, 1)
        self.status_widget = MatlabStatus(window.bridge())
        top.addWidget(self.status_widget)
        outer.addLayout(top)

        split = QSplitter()
        outer.addWidget(split, 1)
        split.addWidget(self._controls())
        split.addWidget(self._results())
        split.setStretchFactor(0, 0)
        split.setStretchFactor(1, 1)
        split.setSizes([380, 900])

    # ── layout ──────────────────────────────────────────────────────────
    def _controls(self) -> QWidget:
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setMinimumWidth(360)
        body = QWidget()
        form = QFormLayout(body)
        scroll.setWidget(body)

        csv_row = QHBoxLayout()
        self.csv = QLineEdit()
        self.csv.setPlaceholderText("hydrophone axial-profile CSV")
        csv_row.addWidget(self.csv, 1)
        b = QToolButton(text="…")
        b.clicked.connect(self._browse_csv)
        csv_row.addWidget(b)
        self.load_btn = QPushButton("Load")
        self.load_btn.clicked.connect(self._load_csv)
        csv_row.addWidget(self.load_btn)
        form.addRow("Profile CSV", csv_row)

        self.measured = QListWidget()
        self.measured.setMaximumHeight(120)
        self.measured.setToolTip("Measured depths to fit (Calibrate at: measured depths)")
        form.addRow("Fit depths (measured)", self.measured)

        d_row = QHBoxLayout()
        self.depths = QLineEdit()
        self.depths.setPlaceholderText("20:5:60 or 30, 42.5, 50")
        self.depths.setToolTip("Depths from the exit plane")
        d_row.addWidget(self.depths, 1)
        sd = QPushButton("Study depths")
        sd.setToolTip("The depths this study uses with this device: focal depth settings (config and "
                      "manifest) and the target depth of each placement (exit plane → focus_pos)")
        sd.clicked.connect(self.use_study_depths)
        d_row.addWidget(sd)
        self.predict_btn = QPushButton("Predict")
        self.predict_btn.setToolTip("Profiles at these depths (measured, or interpolated from the measured "
                                    "ones) and the geometric, library and fitted models as a study "
                                    "resolves them")
        self.predict_btn.clicked.connect(self.predict)
        d_row.addWidget(self.predict_btn)
        form.addRow("Desired depths (mm)", d_row)

        p_row = QHBoxLayout()
        self.align = QComboBox()
        self.align.addItems(["label", "peak"])
        self.align.setToolTip("Focus alignment when interpolating between measured depths: column "
                              "focal-depth labels (default) or intensity peaks beyond the near field")
        p_row.addWidget(self.align)
        p_row.addWidget(QLabel("extrapolate ≤"))
        self.extrap = _spin(0, 0, 50, 1)
        self.extrap.setSuffix(" mm")
        self.extrap.setToolTip("Desired depths up to this far outside the measured range get the nearest "
                               "measured profile shifted; 0 = no extrapolation")
        p_row.addWidget(self.extrap)
        form.addRow("Desired profiles", p_row)

        self.calat = QComboBox()
        self.calat.addItem("Fit depths (measured)", "measured")
        self.calat.addItem("Desired depths", "desired")
        form.addRow("Calibrate at", self.calat)
        self.mode = QComboBox()
        self.mode.addItem("Per-depth refinement", "per_depth")
        self.mode.addItem("Static model (joint, debug)", "static")
        self.mode.currentIndexChanged.connect(self._mode_changed)
        form.addRow("Fit mode", self.mode)
        fp = QHBoxLayout()
        self.fp_phase = QCheckBox("phase", checked=True)
        self.fp_amp = QCheckBox("gain", checked=True)
        self.fp_c = QCheckBox("c steer")
        for w in (self.fp_phase, self.fp_amp, self.fp_c):
            fp.addWidget(w)
        form.addRow("Free parameters", fp)
        self.fwd = QComboBox()
        self.fwd.addItems(["rayleigh", "oneil"])
        form.addRow("Forward model", self.fwd)
        self.skip = _spin(10, 0, 50, 1)
        self.skip.setSuffix(" mm")
        form.addRow("Near-field skip", self.skip)
        self.domain = QComboBox()
        self.domain.addItems(["intensity", "amplitude"])
        form.addRow("Loss on", self.domain)
        self.huber = _spin(0.8, 0.001, 1, 3)
        form.addRow("Huber δ", self.huber)
        st_row = QHBoxLayout()
        self.nstarts = QSpinBox(minimum=1, maximum=200, value=16)
        st_row.addWidget(self.nstarts)
        self.adv_btn = QPushButton("Advanced…")
        self.adv_btn.setToolTip("Phase constraints, absolute-intensity fit, global optimisers, focal weighting")
        self.adv_btn.clicked.connect(self._advanced)
        st_row.addWidget(self.adv_btn)
        form.addRow("Starts", st_row)
        self.isppa = _spin(30, 0.1, 1000, 1)
        self.isppa.setSuffix(" W/cm²")
        self.isppa.setToolTip("Reference ISPPA of the stored drives; defaults to transducer.target_isppa_wcm2")
        form.addRow("Reference ISPPA", self.isppa)

        a_row = QHBoxLayout()
        self.eval_btn = QPushButton("Evaluate measured")
        self.eval_btn.setToolTip("Geometric, library and fit at the measured depths (R² against the hydrophone)")
        self.eval_btn.clicked.connect(self.evaluate_measured)
        a_row.addWidget(self.eval_btn)
        self.fit_btn = QPushButton("Fit")
        self.fit_btn.setProperty("primary", True)
        self.fit_btn.clicked.connect(self.fit)
        a_row.addWidget(self.fit_btn)
        form.addRow(a_row)

        self.lib_label = QLabel()
        self.lib_label.setWordWrap(True)
        self.lib_label.setObjectName("CardSub")
        form.addRow("Save to", self.lib_label)
        s_row = QHBoxLayout()
        self.save_btn = QPushButton("Save entry")
        self.save_btn.setToolTip("Write the fitted entry to transducer.library_path (with meta.fit)")
        self.save_btn.clicked.connect(self.save)
        s_row.addWidget(self.save_btn)
        ep = QPushButton("Export profiles…")
        ep.clicked.connect(self._export_profiles)
        s_row.addWidget(ep)
        es = QPushButton("Export settings…")
        es.clicked.connect(self._export_settings)
        s_row.addWidget(es)
        form.addRow(s_row)

        k_row = QHBoxLayout()
        self.kw_grid = QComboBox()
        self.kw_grid.addItems(["3D", "axisymmetric"])
        self.kw_grid.setToolTip("Free-water grid, unless the study config sets grid.axisymmetric")
        k_row.addWidget(self.kw_grid)
        self.kw_btn = QPushButton("k-Wave drive check…")
        self.kw_btn.setToolTip("Per-depth fits: simulate each refined depth in free-water k-Wave with the "
                               "study config's grid and platform, correct the drive to the simulated ISPPA "
                               "and save the entry. Runs on the Run tab (minutes per depth).")
        self.kw_btn.clicked.connect(self.kwave_check)
        k_row.addWidget(self.kw_btn)
        form.addRow("k-Wave", k_row)

        self.log = QPlainTextEdit(readOnly=True)
        self.log.setMinimumHeight(140)
        form.addRow(self.log)
        return scroll

    def _results(self) -> QWidget:
        w = QWidget()
        v = QVBoxLayout(w)
        v.setContentsMargins(0, 0, 0, 0)
        h = QHBoxLayout()
        h.addWidget(QLabel("Show depth"))
        self.show_depth = QComboBox()
        self.show_depth.currentIndexChanged.connect(lambda _i: self._draw())
        h.addWidget(self.show_depth)
        h.addStretch(1)
        v.addLayout(h)
        if FigureCanvasQTAgg is not None:
            self.fig = Figure(figsize=(7, 6), layout="constrained")
            self.ax_prof, self.ax_r2 = self.fig.subplots(2, 1, height_ratios=[1.3, 1])
            self.canvas = FigureCanvasQTAgg(self.fig)
            v.addWidget(self.canvas, 3)
        else:
            self.fig = self.canvas = None
            v.addWidget(QLabel("Install matplotlib to see the profiles."), 3)
        self.table = QTableWidget()
        self.table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        v.addWidget(self.table, 1)
        return w

    # ── state ───────────────────────────────────────────────────────────
    @property
    def model(self):
        return self.window.model

    def _say(self, msg: str) -> None:
        self.log.appendPlainText(f"[{_dt.datetime.now():%H:%M:%S}] {msg}")

    def _set_busy(self, busy: bool, msg: str = "") -> None:
        self._busy = busy
        for b in (self.load_btn, self.predict_btn, self.eval_btn, self.fit_btn, self.save_btn, self.kw_btn):
            b.setEnabled(not busy)
        self.tx.setEnabled(not busy)
        if msg:
            self._say(msg)

    def _transducers(self) -> list[tuple[int, str, str]]:
        m = self.model
        out = []
        for ti in range(m.n_transducers):
            serial = m.get(f"{TX}.serial", ti)
            if is_empty(serial):
                continue
            ds = m.get(f"{TX}.combo.ds_serial", ti)
            out.append((ti, str(serial), "" if is_empty(ds) else str(ds)))
        return out

    def current(self) -> tuple[int, str, str] | None:
        data = self.tx.currentData()
        return tuple(data) if data else None

    def library_path(self) -> str:
        cur = self.current()
        lib = self.model.get(f"{TX}.library_path", cur[0] if cur else 0)
        return "" if is_empty(lib) else str(lib)

    def refresh(self) -> None:
        """Called when the sub-tab is shown: follow the config's transducers (opens the
        selected one in MATLAB, once per transducer / driving system / library folder)."""
        items = self._transducers()
        keep = self.current()
        self.tx.blockSignals(True)
        self.tx.clear()
        for ti, serial, ds in items:
            dev = self.model.library.get(serial) if self.model.library else None
            kind = f" ({dev.type})" if dev is not None and dev.type else ""
            self.tx.addItem(f"Transducer {ti + 1}: {serial}{' / ' + ds if ds else ''}{kind}", [ti, serial, ds])
        if keep is not None:
            idx = self.tx.findData(list(keep))
            if idx >= 0:
                self.tx.setCurrentIndex(idx)
        self.tx.blockSignals(False)
        lib = self.library_path()
        self.lib_label.setText(f"<tt>{lib}</tt>" if lib else
                               "<span style='color:#b8860b'>transducer.library_path is not set: Save "
                               "asks for a folder and sets it.</span>")
        if not items:
            self.coverage.setText("No transducer from the equipment library in this config "
                                  "(set transducer.serial on the Settings tab).")
            return
        self._open()
        cur = self.current()
        if not self._busy and cur is not None and self.opened == (cur[1], cur[2], self.library_path()):
            if self._prefill_depths():          # the config or placement changed since
                self.predict()

    def _open(self, then=None) -> None:
        cur = self.current()
        if cur is None:
            return
        key = (cur[1], cur[2], self.library_path())
        if key == self.opened and then is None:
            return
        target = self.model.get(f"{TX}.target_isppa_wcm2", cur[0])
        if isinstance(target, (int, float)) and not isinstance(target, bool) and target > 0:
            self.isppa.setValue(float(target))
        self._set_busy(True, f"Opening {cur[1]}{' / ' + cur[2] if cur[2] else ''} in MATLAB…")
        self.window.bridge().call(
            "nexus_calibration",
            {"action": "open", "session": self.session_id, "serial": cur[1], "driving_system": cur[2],
             "library_path": key[2]},
            on_done=lambda out: self._opened(key, out, then), on_error=self._failed)

    def _opened(self, key, out: dict, then=None) -> None:
        self._set_busy(False)
        if self.opened is None or self.opened[:2] != key[:2]:
            self.auto_depths = None             # another device: the old depths do not apply
            self.depths.clear()
        self.opened = key
        self.info = out
        self.has_fit = False
        self.view = None
        self.csv.setText(out.get("csv") or "")
        self._fill_measured()
        entry = (f"library entry {Path(out['library_file']).name} ({out.get('calibration_type') or '?'}"
                 f"{', ' + out['last_updated'] if out.get('last_updated') else ''})"
                 if out.get("library_file") else "no library entry")
        self.coverage.setText(f"<b>{out['combo']}</b> · {entry}")
        self._say(f"{out['combo']}: {entry}")
        if out.get("matrix"):
            self._say("Matrix transducer: geometry only. Prediction and fitting support annular arrays.")
            return
        rec = out.get("fit_record")
        fallback = as_floats(out.get("default_depths"))
        if isinstance(rec, dict):
            s, depths = cal.settings_from_record(rec, self.settings())
            self.apply_settings(s)
            fallback = depths or fallback
            self._say(f"Fit settings taken from the entry (meta.fit: {rec.get('mode', '')}, "
                      f"{rec.get('tool', '')}).")
        self._prefill_depths(fallback)
        if then is not None:
            then()
        elif out.get("csv"):
            self.evaluate_measured()
        else:
            self.predict()

    def _fill_measured(self) -> None:
        self.measured.clear()
        for d in as_floats(self.info.get("measured_depths")):
            it = QListWidgetItem(f"{d:.2f}")
            it.setData(Qt.ItemDataRole.UserRole, d)
            it.setFlags(it.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            it.setCheckState(Qt.CheckState.Checked)
            self.measured.addItem(it)

    def _manifest(self):
        from .main_window import STUDY
        study = self.window.pages.get(STUDY)
        return getattr(study, "manifest", None)

    def _failed(self, err: BridgeError) -> None:
        self._set_busy(False)
        if err.identifier == "nexus_calibration:noSession":
            self.opened = None              # MATLAB restarted: open again on the next action
        self._say(f"✗ {err.message}")

    # ── settings ────────────────────────────────────────────────────────
    def settings(self) -> cal.FitSettings:
        free = [p for p, w in (("phase", self.fp_phase), ("amp", self.fp_amp)) if w.isChecked()]
        return cal.FitSettings(
            free_params=free, ring_profile=self.ring.currentData(), forward_model=self.fwd.currentText(),
            near_field_skip_mm=self.skip.value(), domain=self.domain.currentText(),
            huber_delta=self.huber.value(), n_starts=self.nstarts.value(), isppa_wcm2=self.isppa.value(),
            align=self.align.currentText(), max_extrapolation_mm=self.extrap.value(),
            calibrate_at=self.calat.currentData(), advanced=dict(self.advanced))

    def apply_settings(self, s: cal.FitSettings) -> None:
        self.ring.setCurrentIndex(max(0, self.ring.findData(s.ring_profile)))
        self.fp_phase.setChecked("phase" in s.free_params)
        self.fp_amp.setChecked("amp" in s.free_params)
        self.fwd.setCurrentText(s.forward_model)
        self.skip.setValue(s.near_field_skip_mm)
        self.domain.setCurrentText(s.domain)
        self.huber.setValue(s.huber_delta)
        self.nstarts.setValue(s.n_starts)
        self.isppa.setValue(s.isppa_wcm2)
        self.align.setCurrentText(s.align)
        self.extrap.setValue(s.max_extrapolation_mm)
        self.calat.setCurrentIndex(self.calat.findData(s.calibrate_at))
        self.advanced = dict(s.advanced)
        self._update_advanced_button()

    def _advanced(self) -> None:
        dlg = AdvancedDialog(self.advanced, self)
        if dlg.exec():
            self.advanced = dlg.values()
            self._update_advanced_button()

    def _update_advanced_button(self) -> None:
        n = len(self.advanced)
        self.adv_btn.setText(f"Advanced ({n})" if n else "Advanced…")

    def _study_depths(self) -> tuple[list[float], str]:
        """Desired depths for the open device from the config, manifest and placement, and
        where they came from."""
        cur = self.current()
        if cur is None:
            return [], ""
        settings = cal.study_depths(self.model, cur[1], cur[2], self._manifest())
        targets = cal.target_depths(self.model, cur[1], cur[2], self._target_depth)
        return cal.desired_depths(settings, targets), cal.describe_depths(settings, targets)

    def _target_depth(self, ti: int) -> float | None:
        from .main_window import VIEWER
        viewer = self.window.pages.get(VIEWER)
        try:
            return viewer.target_depth_ep(ti) if viewer is not None else None
        except Exception:                       # unreadable T1 etc.: no target depth
            return None

    def _prefill_depths(self, fallback: list[float] | None = None) -> bool:
        """Fill the desired depths with the study's depths (else fallback), unless the user
        typed their own. True when the text changed."""
        text = self.depths.text().strip()
        if text and text != self.auto_depths:
            return False
        depths, source = self._study_depths()
        new = cal.format_depths(depths or fallback or [])
        if not new or new == text:
            return False
        self.depths.setText(new)
        self.auto_depths = new
        self._say(f"Desired depths {new} mm from the study ({source})." if depths else
                  f"Desired depths {new} mm: no focal depth setting or placement for this device yet.")
        return True

    def use_study_depths(self) -> None:
        depths, source = self._study_depths()
        if not depths:
            self._say("No focal depth setting (transducer.focal_distance_ep / _bowl, config or "
                      "manifest) or placement target for this device.")
            return
        self.depths.setText(cal.format_depths(depths))
        self.auto_depths = self.depths.text()
        self._say(f"Study depths: {source}.")
        self.predict()

    # ── actions ─────────────────────────────────────────────────────────
    def _call(self, action: str, args: dict, done, busy_msg: str = "") -> None:
        """nexus_calibration action on the open transducer (opened first if needed)."""
        cur = self.current()
        if cur is None:
            return
        if self.opened != (cur[1], cur[2], self.library_path()):
            self._open(then=lambda: self._call(action, args, done, busy_msg))
            return
        self._set_busy(True, busy_msg)
        payload = dict(args, action=action, session=self.session_id)
        self.window.bridge().call("nexus_calibration", payload,
                                  on_done=lambda out: (self._set_busy(False), done(out)),
                                  on_error=self._failed, on_progress=self._say)

    def _browse_csv(self) -> None:
        start = self.csv.text() or str(self.window.prestus_root / "local" / "profiles")
        path, _ = QFileDialog.getOpenFileName(self, "Axial profile CSV", start, "CSV (*.csv)")
        if path:
            self.csv.setText(path)
            self._load_csv()

    def _load_csv(self) -> None:
        csv = self.csv.text().strip()
        if not csv:
            return

        def loaded(out):
            self.info = out
            self.has_fit = False
            self._fill_measured()
            span = as_floats(out.get("measured_span"))
            self._say(f"Loaded {len(as_list(out.get('measured_depths')))} depths "
                      f"({span[0]:.1f}–{span[1]:.1f} mm) from {csv}" if span else f"Loaded {csv}")
            self.evaluate_measured()
        self._call("load", {"csv": csv}, loaded)

    def predict(self) -> None:
        depths = cal.parse_depths(self.depths.text())
        if not depths:
            self._say("Desired depths: enter depths such as 20:5:60 or 30, 42.5, 50.")
            return
        self._call("view", {"predicted": True, "depths": depths,
                            "settings": self.settings().bridge_settings()}, self._show)

    def evaluate_measured(self) -> None:
        if not self.info.get("csv"):
            self.predict()
            return
        self._call("view", {"predicted": False, "settings": self.settings().bridge_settings()}, self._show)

    def fit(self) -> None:
        if not self.info.get("csv"):
            self._say("Fitting needs measured profiles: load a CSV first.")
            return
        s = self.settings()
        args = {"settings": s.bridge_settings(), "mode": s.mode, "free_params": s.free_params,
                "calibrate_at": s.calibrate_at,
                "predicted": bool(self.view and self.view.get("predicted"))}
        depths = cal.parse_depths(self.depths.text())
        if s.calibrate_at == "desired" or args["predicted"]:
            if not depths:
                self._say("Desired depths: enter depths such as 20:5:60 or 30, 42.5, 50.")
                return
            args["depths"] = depths
        if s.calibrate_at == "measured":
            sel = [self.measured.item(i).data(Qt.ItemDataRole.UserRole) for i in range(self.measured.count())
                   if self.measured.item(i).checkState() == Qt.CheckState.Checked]
            if not sel:
                self._say("Select at least one measured depth.")
                return
            args["measured_depths"] = sel

        def fitted(out):
            self.has_fit = True
            self._show(out)
        self._call("fit", args, fitted, "Fitting…")

    def save(self) -> None:
        if not self.has_fit:
            self._say("Nothing to save: run a fit first.")
            return
        cur = self.current()
        out_dir = self.library_path()
        if not out_dir:
            start = str(self.model.path.parent if self.model.path else Path.home())
            out_dir = QFileDialog.getExistingDirectory(self, "Folder for your library entries", start)
            if not out_dir:
                return
            self.model.set(f"{TX}.library_path", out_dir, cur[0])
            self.window._changed()
            self.opened = (cur[1], cur[2], out_dir)     # same MATLAB state; only the folder is new
            self._say(f"transducer.library_path set to {out_dir} (save the config to keep it).")
            self.refresh()
        shipped = self.info.get("shipped_path", "")
        allow = False
        if shipped and Path(out_dir).resolve() == Path(shipped).resolve():
            allow = QMessageBox.question(
                self, "Shipped library", "This overwrites the SHIPPED library entry. Use a folder of your "
                "own and transducer.library_path unless you maintain the library. Overwrite?") \
                == QMessageBox.StandardButton.Yes
            if not allow:
                return

        def saved(out):
            msg = f"Saved {out['file']}"
            if out.get("profiles_csv"):
                msg += f" (+ calibration profiles {Path(out['profiles_csv']).name})"
            self._say(msg)
            self.window._changed()     # the coverage line on the Settings tab reads the new entry
        self._call("save", {"out_dir": out_dir, "allow_shipped": allow}, saved)

    def _export(self, action: str, what: str) -> None:
        name = f"{self.info.get('combo', 'calibration')}_{what}.csv"
        start = str((self.model.path.parent if self.model.path else Path.home()) / name)
        path, _ = QFileDialog.getSaveFileName(self, f"Export {what}", start, "CSV (*.csv)")
        if not path:
            return
        self._call(action, {"path": path, "settings": self.settings().bridge_settings()},
                   lambda out: self._say(f"Exported {out.get('what', '')} {what} to {out['path']}"))

    def _export_profiles(self) -> None:
        self._export("export_profiles", "profiles")

    def _export_settings(self) -> None:
        self._export("export_settings", "settings")

    def kwave_check(self) -> None:
        """Write the unsaved per-depth fit and a script that runs kwave_calibrate_drive and saves
        the corrected entry; run it on the Run tab (or in MATLAB / on the cluster yourself)."""
        from .. import launcher
        from .main_window import RUN
        if not self.has_fit:
            self._say("The k-Wave drive check needs a per-depth fit: run Fit first.")
            return
        if self.model.path is None and not self.window.save():
            return
        lib = self.library_path()
        if not lib:
            self._say("Set transducer.library_path (or Save the analytical fit once) first: the "
                      "corrected entry is saved there.")
            return
        combo = self.info.get("combo", "calibration")
        folder = self.model.path.parent / "calibration"
        mat = folder / f"{combo}_fit.mat"
        script = folder / f"kwave_drive_{combo}.m"

        def written(_out):
            folder.mkdir(parents=True, exist_ok=True)
            launcher.write_kwave_drive_script(
                script, mat, self.window.prestus_root, lib, str(self.model.path),
                self.kw_grid.currentText() == "axisymmetric", self.skip.value())
            box = QMessageBox(self)
            box.setWindowTitle("k-Wave drive check")
            box.setText(f"Script written:\n{script}\n\nIt runs one free-water simulation per refined depth "
                        "(minutes each, with the study config's grid and platform) and saves the corrected "
                        f"entry to {lib}.\n\nRun it now on the Run tab, or run the script yourself "
                        "(e.g. on the cluster)?")
            run = box.addButton("Run on the Run tab", QMessageBox.ButtonRole.AcceptRole)
            box.addButton("Only write the script", QMessageBox.ButtonRole.RejectRole)
            box.exec()
            if box.clickedButton() is run:
                page = self.window.show_special(RUN)
                page.run_script(f"k-Wave drive check ({combo})", script)
            else:
                self._say(f"k-Wave drive check script: {script}")
        self._call("write_fit", {"path": str(mat)}, written)

    # ── drawing ─────────────────────────────────────────────────────────
    def _show(self, out: dict) -> None:
        self.view = out
        for msg in as_list(out.get("messages")):
            self._say(msg)
        depths = as_floats(out.get("depths"))
        prev = self.show_depth.currentText()
        self.show_depth.blockSignals(True)
        self.show_depth.clear()
        items = ([f"{d:.2f}" for d in depths])
        if out.get("predicted"):
            items = ["All depths"] + items
        self.show_depth.addItems(items)
        if prev in items:
            self.show_depth.setCurrentText(prev)
        self.show_depth.blockSignals(False)
        self._draw()
        self._fill_table()

    def _draw(self) -> None:
        if self.canvas is None or self.view is None:
            return
        v = self.view
        ax, ax2 = self.ax_prof, self.ax_r2
        ax.clear()
        ax2.clear()
        res = v.get("results") or {}
        types = as_list(v.get("profile_type"))
        sel = self.show_depth.currentText()
        show_all = sel in ("All depths", "")
        shown_meas = False
        for kind, colour, name in STYLES:
            r = res.get(kind)
            if not r:
                continue
            d = as_floats(r.get("depths"))
            z = as_floats(r.get("z"))
            imod = as_columns(r.get("Imod"))
            imeas = as_columns(r.get("Imeas"))
            r2 = as_floats(r.get("r2"))
            idx = list(range(len(d))) if show_all else [j for j, x in enumerate(d) if abs(x - float(sel)) < 0.01]
            if not idx:
                continue
            if not shown_meas and imeas:
                for ptype, ls, label in PROFILE_LINES:
                    first = True
                    for j in idx:
                        t = types[j] if j < len(types) else "measured"
                        col = imeas[j] if j < len(imeas) else []
                        if t != ptype or not any(math.isfinite(x) for x in col):
                            continue
                        ax.plot(z, col, ls, color="k", lw=1.6, label=label if first else None)
                        first = False
                shown_meas = True
            label = name
            if not show_all and math.isfinite(r2[idx[0]]):
                label += f" (R² {r2[idx[0]]:.2f})"
            settings = as_list(r.get("settings"))
            if not show_all and settings and idx[0] < len(settings):
                label += f" [{settings[idx[0]]['phase_model']}]"
            for n, j in enumerate(idx):
                ax.plot(z, imod[j], "-", color=colour, lw=1.1, label=label if n == 0 else None)
        depths = as_floats(v.get("depths"))
        for x in depths if show_all else [float(sel)]:
            ax.axvline(x, ls=":", color="#bbbbbb", lw=0.8)
        ax.set_xlabel("Distance from exit plane (mm)")
        ax.set_ylabel("Normalised intensity")
        ax.set_title(f"{v.get('combo', '')} — {'all depths' if show_all else sel + ' mm'}"
                     f"{' (prediction, runtime library path)' if v.get('predicted') else ''}", fontsize=10)
        ax.grid(True, alpha=0.3)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=8, loc="upper right")

        predicted = bool(v.get("predicted"))
        for kind, colour, name in STYLES:
            r = res.get(kind)
            if not r:
                continue
            y = as_floats(r.get("focus" if predicted else "r2"))
            ax2.plot(as_floats(r.get("depths")), y, "-o", color=colour, ms=4, lw=1.1, label=name)
        geo = v.get("dist_geom_mm")
        if isinstance(geo, (int, float)):
            ax2.axvline(geo, ls=":", color="#888888", lw=0.8)
            ax2.annotate("natural focus", (geo, 0.02), xycoords=("data", "axes fraction"), fontsize=7,
                         color="#666666", rotation=90)
        if predicted and depths:
            lim = [min(depths), max(depths)]
            ax2.plot(lim, lim, "k--", lw=0.9, label="Target")
            fm = as_floats(v.get("measured_focus"))
            for ptype, face, label in (("measured", "k", "Hydrophone"), ("interpolated", "none", "Interpolated profile"),
                                       ("extrapolated", "none", "Extrapolated profile")):
                pts = [(depths[j], fm[j]) for j in range(min(len(fm), len(depths)))
                       if j < len(types) and types[j] == ptype and math.isfinite(fm[j])]
                if pts:
                    ax2.plot(*zip(*pts), "*", ms=9, mfc=face, mec="k", ls="none", label=label)
            span = as_floats(v.get("measured_span"))
            if len(span) == 2:
                ax2.axvspan(span[0], span[1], color="#888888", alpha=0.08, label="Measured range")
            ax2.set_xlabel("Target depth (mm from EP)")
            ax2.set_ylabel("Focus, FLHM centre (mm)")
        else:
            lo = min([0.0] + [y for k, *_ in STYLES if res.get(k) for y in as_floats(res[k].get("r2"))
                              if math.isfinite(y)])
            ax2.set_ylim(max(-0.2, lo), 1)
            ax2.set_xlabel("Focal depth (mm from EP)")
            ax2.set_ylabel("R² (profile shape)")
        ax2.grid(True, alpha=0.3)
        if ax2.get_legend_handles_labels()[0]:
            ax2.legend(fontsize=8, loc="best")
        self.canvas.draw_idle()

    def _fill_table(self) -> None:
        v = self.view or {}
        res = v.get("results") or {}
        t = self.table
        t.clear()
        if v.get("predicted"):
            src = "fit" if res.get("fit") and as_list(res["fit"].get("settings")) else \
                  "library" if res.get("library") and as_list(res["library"].get("settings")) else None
            if src is None:
                t.setRowCount(0)
                t.setColumnCount(0)
                return
            r = res[src]
            settings = as_list(r.get("settings"))
            n_el = max(len(as_list(s.get("elem_phase_deg"))) for s in settings)
            heads = ["Depth (mm)", "Profile", "Entry", "Phase model", "Drive (kPa)", "Focus (mm)", "R²"] + \
                    [f"Ring {i + 1} (°)" for i in range(n_el)]
            depths, focus, r2 = as_floats(r.get("depths")), as_floats(r.get("focus")), as_floats(r.get("r2"))
            types = as_list(v.get("profile_type"))
            rows = []
            for k, s in enumerate(settings):
                drive = s.get("drive_kpa")
                rows.append([f"{depths[k]:g}", str(types[k]).replace("unavailable", "—") if k < len(types) else "",
                             "Fit (unsaved)" if src == "fit" else "Library", s.get("phase_model", ""),
                             "—" if drive is None else f"{drive:.1f}", f"{focus[k]:.1f}",
                             f"{r2[k]:.3f}" if math.isfinite(r2[k]) else "—"] +
                            [f"{x:.1f}" for x in as_floats(s.get("elem_phase_deg"))])
        else:
            p = v.get("params")
            if not isinstance(p, dict):
                t.setRowCount(0)
                t.setColumnCount(0)
                return
            cols = [as_floats(c) for c in as_list(p.get("phase_deg"))]
            if cols and not isinstance(as_list(p.get("phase_deg"))[0], list):
                cols = [as_floats(p.get("phase_deg"))]
            heads = ["Element"] + [f"{d:.1f} mm (°)" for d in as_floats(p.get("depths"))]
            rows = [[str(i + 1)] + [f"{c[i]:.1f}" for c in cols] for i in range(max(map(len, cols)))]
        t.setColumnCount(len(heads))
        t.setHorizontalHeaderLabels(heads)
        t.setRowCount(len(rows))
        for i, row in enumerate(rows):
            for j, text in enumerate(row):
                t.setItem(i, j, QTableWidgetItem(text))
        t.resizeColumnsToContents()
