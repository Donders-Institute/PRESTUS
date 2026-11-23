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

from PyQt6.QtCore import QSize, Qt, QTimer, QUrl
from PyQt6.QtGui import QDesktopServices, QFontDatabase
from PyQt6.QtWidgets import (QCheckBox, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox,
                             QFileDialog, QFormLayout, QFrame, QHBoxLayout, QLabel, QLineEdit,
                             QListWidget, QListWidgetItem, QMenu, QMessageBox, QPlainTextEdit,
                             QPushButton, QScrollArea, QSizePolicy, QSpinBox, QSplitter,
                             QStackedWidget, QTableWidget, QTableWidgetItem, QToolButton,
                             QVBoxLayout, QWidget)

from .. import calibration as cal
from ..bridge import BridgeError, as_columns, as_floats, as_list
from ..model import TX, is_empty
from . import calibration_plots as cplots
from . import icons, style
from .errors import log_current
from .file_badge import FileBadge, fit_combo_to_items
from .matlab_status import MatlabStatus

try:
    from .mpl_theme import ThemedCanvas as FigureCanvasQTAgg   # drawn in the Nexus theme
    from matplotlib.figure import Figure
except ImportError:                                    # pragma: no cover
    FigureCanvasQTAgg = None

# result kind, colour, legend name (as result_styles in calibration_gui.m, in the colour-blind-
# safe style.SERIES hues; geometric stays a neutral grey)
STYLES = [("geometric", "#999999", "Geometric"), ("library", style.SERIES[0], "Library"),
          ("fit", style.SERIES[1], "Fit"), ("kwave", style.SERIES[2], "k-Wave")]
PROFILE_LINES = [("measured", "-", "Hydrophone"), ("interpolated", "--", "Interpolated profile"),
                 ("extrapolated", ":", "Extrapolated profile")]


class AdvancedDialog(QDialog):
    """Advanced options, as "Advanced…" in calibration_gui.m: FIT_TRANSDUCER_MODEL options, the
    desired-profile tolerance and the k-Wave drive calibration. Only values that differ from the
    defaults are kept. "Start from preset" loads a calibration_fit_options preset; the page then
    also applies its main-panel values (preset())."""

    def __init__(self, current: dict, extra: dict, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Advanced fit options")
        outer = QVBoxLayout(self)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        body = QWidget()
        form = QFormLayout(body)
        scroll.setWidget(body)
        outer.addWidget(scroll)
        self.resize(480, 700)

        self.pre = QComboBox()
        self.pre.addItems(["(current settings)", *cal.PRESETS])
        self.pre.setToolTip("calibration_fit_options preset: sets the main panel (free parameters, forward "
                            "model, loss, near-field skip) and these options to the preset's values")
        self.pre.currentTextChanged.connect(self._load_preset)
        form.addRow("Start from preset", self.pre)
        form.addRow(_header("Phase model"))
        self.pc = QComboBox()
        self.pc.addItems(["none", "linear", "monotonic"])
        form.addRow("Phase constraint", self.pc)
        self.sb = QLineEdit()
        self.sb.setPlaceholderText("lo, hi (rad); blank = default")
        form.addRow("Phase step bounds", self.sb)
        self.lam = _spin(0, 0, 100, 6)
        self.lam.setToolTip("Lambda: pull of the phase offsets towards geometric steering")
        form.addRow("Phase prior λ", self.lam)
        self.lama = _spin(0, 0, 100, 6)
        self.lama.setToolTip("LambdaAmp: pull of the log ring gains towards equal gains")
        form.addRow("Gain prior λ", self.lama)
        form.addRow(_header("Gains"))
        self.ab = QLineEdit()
        self.ab.setToolTip("AmpBounds: bounds on each ring gain relative to the prior")
        form.addRow("Gain bounds (× prior)", self.ab)
        self.pg = QLineEdit()
        self.pg.setToolTip("PlausibleGain: gains outside (× mean) are flagged as a possible hardware fault")
        form.addRow("Plausible gain (flag only)", self.pg)
        form.addRow(_header("Amplitude and loss"))
        self.af = QComboBox()
        self.af.addItems(["shape", "absolute"])
        form.addRow("Amplitude fit", self.af)
        self.sv = QCheckBox("one drive for all depths (absolute)")
        form.addRow("", self.sv)
        self.fw = _spin(0, 0, 1000, 3)
        form.addRow("Focal weight (0 = uniform)", self.fw)
        self.win = QLineEdit()
        self.win.setPlaceholderText("lo, hi (mm from EP); blank = all")
        form.addRow("Fit window", self.win)
        form.addRow(_header("Near field"))
        self.nfp = _spin(0, 0, 1000, 3)
        self.nfp.setToolTip("NearFieldPenalty: penalises model near-field intensity above the measured near "
                            "field (from 'From' to the near-field skip); 0 = off")
        form.addRow("Near-field penalty", self.nfp)
        self.nfm = _spin(0, 0, 10, 3)
        self.nfm.setToolTip("NearFieldMargin: excess (× focal peak) tolerated before the penalty applies")
        form.addRow("Margin (× focal peak)", self.nfm)
        self.nfz = _spin(0, 0, 50, 1)
        self.nfz.setSuffix(" mm")
        self.nfz.setToolTip("NearFieldMinMm: start of the penalised near field (from the exit plane)")
        form.addRow("From", self.nfz)
        form.addRow(_header("Medium and optimiser"))
        self.cw = _spin(1500, 1000, 2000, 1)
        self.cw.setSuffix(" m/s")
        form.addRow("Sound speed", self.cw)
        self.rho = _spin(994, 500, 2000, 1)
        self.rho.setSuffix(" kg/m³")
        form.addRow("Density", self.rho)
        self.opt = QComboBox()
        self.opt.addItems(["lsqnonlin", "fexminimize", "globalsearch"])
        form.addRow("Optimiser", self.opt)
        self.pol = QCheckBox("polish")
        form.addRow("", self.pol)
        self.pop = QSpinBox(minimum=10, maximum=100000)
        form.addRow("FEX population", self.pop)
        self.mfe = QSpinBox(minimum=100, maximum=100000000)
        form.addRow("Max. evaluations", self.mfe)
        self.seed = QSpinBox(minimum=0, maximum=2**31 - 1)
        form.addRow("Random seed", self.seed)
        form.addRow(_header("Desired profiles"))
        self.tol = _spin(0.01, 0, 10, 3)
        self.tol.setSuffix(" mm")
        self.tol.setToolTip("Tolerance: a desired depth this close to a measured one uses the measurement")
        form.addRow("Same-depth tolerance", self.tol)
        form.addRow(_header("k-Wave drive check"))
        self.grid = QLineEdit()
        self.grid.setPlaceholderText("blank = study config; auto; or nx ny nz")
        self.grid.setToolTip("GridDims of the free-water grid: 'auto' (sized to the aperture and focus), "
                             "[nx nz] / [nx ny nz], or blank for the study config's grid")
        form.addRow("Grid", self.grid)
        self.apply = QCheckBox("rescale the stored drives to the simulated ISPPA")
        self.apply.setToolTip("ApplyCorrection; off: report only")
        form.addRow("", self.apply)
        box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
                               | QDialogButtonBox.StandardButton.Reset)
        box.accepted.connect(self.accept)
        box.rejected.connect(self.reject)
        box.button(QDialogButtonBox.StandardButton.Reset).clicked.connect(self._reset)
        outer.addWidget(box)
        self._reset_all = False
        self._fill(dict(cal.ADVANCED_DEFAULTS, **current), dict(cal.EXTRA_DEFAULTS, **extra))

    def _fill(self, cur: dict, ext: dict) -> None:
        self.pc.setCurrentText(str(cur["PhaseConstraint"]))
        self.sb.setText(_pair_text(cur["PhaseStepBounds"]))
        self.lam.setValue(float(cur["Lambda"]))
        self.lama.setValue(float(cur["LambdaAmp"]))
        self.ab.setText(_pair_text(cur["AmpBounds"]))
        self.pg.setText(_pair_text(cur["PlausibleGain"]))
        self.af.setCurrentText(str(cur["AmplitudeFit"]))
        self.sv.setChecked(bool(cur["SharedVelocity"]))
        self.fw.setValue(float(cur["FocalWeight"]))
        self.win.setText(_pair_text(cur["FitWindowMm"]))
        self.nfp.setValue(float(cur["NearFieldPenalty"]))
        self.nfm.setValue(float(cur["NearFieldMargin"]))
        self.nfz.setValue(float(cur["NearFieldMinMm"]))
        self.cw.setValue(float(cur["CWater"]))
        self.rho.setValue(float(cur["Rho"]))
        self.opt.setCurrentText(str(cur["Optimizer"]))
        self.pol.setChecked(bool(cur["Polish"]))
        self.pop.setValue(int(cur["PopSize"]))
        self.mfe.setValue(int(cur["MaxFunEvals"]))
        self.seed.setValue(int(cur["Seed"]))
        self.tol.setValue(float(ext["Tolerance"]))
        g = ext["GridDims"]
        self.grid.setText(" ".join(f"{float(x):g}" for x in g) if isinstance(g, (list, tuple)) else str(g or ""))
        self.apply.setChecked(bool(ext["ApplyCorrection"]))

    def _load_preset(self, name: str) -> None:
        if name not in cal.PRESETS:
            return
        p = cal.PRESETS[name]
        self._fill(dict(cal.ADVANCED_DEFAULTS, **{k: v for k, v in p.items() if k in cal.ADVANCED_DEFAULTS}),
                   dict(cal.EXTRA_DEFAULTS, **self.extra()))

    def preset(self) -> dict | None:
        """The preset chosen in this dialog (its main-panel values apply too), or None."""
        return None if self._reset_all else cal.PRESETS.get(self.pre.currentText())

    def _reset(self) -> None:
        self._reset_all = True
        self.accept()

    def grid_dims(self):
        t = self.grid.text().strip()
        if not t or t == "auto":
            return t
        try:
            v = [int(float(x)) for x in t.replace(",", " ").replace("[", " ").replace("]", " ").split()]
        except ValueError:
            return None
        return v if len(v) in (2, 3) else None

    def values(self) -> dict:
        if self._reset_all:
            return {}
        v = {"PhaseConstraint": self.pc.currentText(), "PhaseStepBounds": _pair(self.sb.text()),
             "AmplitudeFit": self.af.currentText(), "SharedVelocity": self.sv.isChecked(),
             "Optimizer": self.opt.currentText(), "PopSize": self.pop.value(),
             "MaxFunEvals": self.mfe.value(), "Polish": self.pol.isChecked(),
             "FocalWeight": self.fw.value(), "FitWindowMm": _pair(self.win.text()),
             "Lambda": self.lam.value(), "LambdaAmp": self.lama.value(),
             "AmpBounds": _pair(self.ab.text()) or cal.ADVANCED_DEFAULTS["AmpBounds"],
             "PlausibleGain": _pair(self.pg.text()) or cal.ADVANCED_DEFAULTS["PlausibleGain"],
             "NearFieldPenalty": self.nfp.value(), "NearFieldMargin": self.nfm.value(),
             "NearFieldMinMm": self.nfz.value(), "CWater": self.cw.value(), "Rho": self.rho.value(),
             "Seed": self.seed.value()}
        return cal.non_default(v, cal.ADVANCED_DEFAULTS)

    def extra(self) -> dict:
        if self._reset_all:
            return {}
        g = self.grid_dims()
        v = {"Tolerance": self.tol.value(), "GridDims": "" if g is None else g,
             "ApplyCorrection": self.apply.isChecked()}
        return cal.non_default(v, cal.EXTRA_DEFAULTS)


class RingOptionsDialog(QDialog):
    """fit_ring_profile options and plane scans (read_offaxis_scans), as "Ring profile options" in
    calibration_gui.m; only values that differ from the defaults are kept."""

    def __init__(self, current: dict, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Ring profile options")
        cur = dict(cal.RING_DEFAULTS, **current)
        form = QFormLayout(self)
        form.addRow(_header("Profile model (used when Ring profile = fit)"))
        self.nsub = QSpinBox(minimum=1, maximum=32, value=int(cur["NSub"]))
        form.addRow("Sub-rings per ring", self.nsub)
        self.sm = _spin(cur["Smooth"], 0, 100, 3)
        self.sm.setToolTip("Smooth: penalty between neighbouring sub-rings; 0 = none")
        form.addRow("Smoothness", self.sm)
        self.skip = _spin(cur["NearFieldSkipMm"], 0, 50, 1)
        self.skip.setSuffix(" mm")
        form.addRow("Near-field skip", self.skip)
        self.ab = QLineEdit(_pair_text(cur["AmpBounds"]))
        self.ab.setToolTip("AmpBounds: bounds on sub-ring weights and channel gains")
        form.addRow("Weight bounds", self.ab)
        form.addRow(_header("Optimiser"))
        self.it = QSpinBox(minimum=1, maximum=100, value=int(cur["Iterations"]))
        form.addRow("Alternations", self.it)
        self.ns = QSpinBox(minimum=1, maximum=200, value=int(cur["NStarts"]))
        form.addRow("Starts per depth", self.ns)
        self.seed = QSpinBox(minimum=0, maximum=2**31 - 1, value=int(cur["Seed"]))
        form.addRow("Random seed", self.seed)
        self.ho = QCheckBox("leave-one-depth-out check (slower)", checked=bool(cur["HeldOut"]))
        form.addRow("", self.ho)
        form.addRow(_header("Plane / volume scans (off-axis data)"))
        f_row = QHBoxLayout()
        self.file = QLineEdit(str(cur["OffAxis"] or ""))
        self.file.setPlaceholderText("blank = axial profiles only")
        self.file.setToolTip("Sonic Concepts FUS Metrology export (.mat, folder, or .zip) or a CSV with "
                             "focus_ep_mm, x_mm, y_mm, z_ep_mm, intensity_wcm2 (read_offaxis_scans)")
        f_row.addWidget(self.file, 1)
        b = icons.PathButton("file", "Choose the off-axis scan file (.mat, folder, .zip, or CSV)")
        b.clicked.connect(self._browse)
        f_row.addWidget(b)
        form.addRow("File", f_row)
        self.kind = QComboBox()
        self.kind.addItems(["all", "transverse", "axial", "volume"])
        self.kind.setCurrentText(str(cur["Kind"]))
        self.kind.setToolTip("transverse: XY planes; axial: XZ / YZ planes; volume: 3D scans")
        form.addRow("Scans used", self.kind)
        self.wt = _spin(cur["Weight"], 0, 100, 2)
        self.wt.setToolTip("Weight of each plane relative to one axial profile")
        form.addRow("Weight per plane", self.wt)
        self.rmax = _spin(cur["RMaxMm"], 1, 200, 1)
        self.rmax.setSuffix(" mm")
        form.addRow("Radial range", self.rmax)
        self.bin = _spin(cur["BinMm"], 0.05, 10, 2)
        self.bin.setSuffix(" mm")
        form.addRow("Radial bin", self.bin)
        self.zbin = _spin(cur["ZBinMm"], 0.05, 10, 2)
        self.zbin.setSuffix(" mm")
        self.zbin.setToolTip("Axial bin of XZ / YZ planes and volumes: coarser = faster")
        form.addRow("Axial-plane / volume z bin", self.zbin)
        self.cen = QComboBox()
        self.cen.addItems(["origin", "centroid"])
        self.cen.setCurrentText(str(cur["Centre"]))
        form.addRow("Transverse / volume centre", self.cen)
        box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
                               | QDialogButtonBox.StandardButton.Reset)
        box.accepted.connect(self._accept)
        box.rejected.connect(self.reject)
        box.button(QDialogButtonBox.StandardButton.Reset).clicked.connect(self._reset)
        form.addRow(box)
        self._reset_all = False

    def _browse(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Plane scans", self.file.text() or str(Path.home()),
                                              "Plane scans (*.mat *.zip *.csv)")
        if path:
            self.file.setText(path)

    def _accept(self) -> None:
        f = self.file.text().strip()
        if f and not Path(f).exists():
            QMessageBox.warning(self, "Ring profile options", f"Not found: {f}")
            return
        self.accept()

    def _reset(self) -> None:
        self._reset_all = True
        self.accept()

    def values(self) -> dict:
        if self._reset_all:
            return {}
        v = {"NSub": self.nsub.value(), "Smooth": self.sm.value(), "NearFieldSkipMm": self.skip.value(),
             "AmpBounds": _pair(self.ab.text()) or cal.RING_DEFAULTS["AmpBounds"],
             "Iterations": self.it.value(), "NStarts": self.ns.value(), "Seed": self.seed.value(),
             "HeldOut": self.ho.isChecked(), "OffAxis": self.file.text().strip(),
             "Kind": self.kind.currentText(), "Weight": self.wt.value(), "RMaxMm": self.rmax.value(),
             "BinMm": self.bin.value(), "ZBinMm": self.zbin.value(), "Centre": self.cen.currentText()}
        return cal.non_default(v, cal.RING_DEFAULTS)


def _header(text: str) -> QLabel:
    lab = QLabel(f"<b>{text}</b>")
    return lab


def _step(number: int, title: str, sub: str = "") -> tuple[QFrame, QVBoxLayout]:
    """A numbered step group of the controls (#SubGroup look from style.py)."""
    box = QFrame()
    box.setObjectName("SubGroup")
    v = QVBoxLayout(box)
    v.setContentsMargins(10, 8, 10, 10)
    v.setSpacing(6)
    head = QLabel(f"{number}&nbsp;&nbsp;{title}" + (f"&nbsp;&nbsp;<span style='font-weight:400'>{sub}</span>"
                                                     if sub else ""))
    head.setObjectName("SubGroupTitle")
    head.setTextFormat(Qt.TextFormat.RichText)
    v.addWidget(head)
    return box, v


def _form() -> QFormLayout:
    """A form that wraps a long label above its field instead of widening the column."""
    f = QFormLayout()
    f.setContentsMargins(0, 0, 0, 0)
    f.setHorizontalSpacing(8)
    f.setVerticalSpacing(6)
    f.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
    f.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
    return f


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


def _theme_axes(axes, c: dict) -> None:
    """Axes, ticks, grid, title and legend in the active palette (light or dark)."""
    for ax in axes:
        ax.set_facecolor(c["surface"])
        ax.tick_params(colors=c["ink3"], labelsize=style.PLOT_PT["small"])
        ax.xaxis.label.set_color(c["ink2"])
        ax.yaxis.label.set_color(c["ink2"])
        ax.title.set_color(c["ink"])
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(c["line2"])
        for line in ax.get_xgridlines() + ax.get_ygridlines():
            line.set_color(c["line"])
        leg = ax.get_legend()
        if leg is not None:
            leg.get_frame().set_facecolor(c["surface"])
            leg.get_frame().set_edgecolor(c["line"])
            for txt in leg.get_texts():
                txt.set_color(c["ink2"])


class CalibrationPage(QWidget):
    def __init__(self, window, parent: QWidget | None = None):
        super().__init__(parent)
        self.window = window
        self.session_id = f"nexus-{id(self):x}"
        self.opened: tuple | None = None       # (serial, ds, library_path) open in MATLAB
        self.info: dict = {}                   # last 'open' / 'load' answer
        self.view: dict | None = None          # last 'view' / 'fit' answer
        self.advanced: dict = {}
        self.extra: dict = {}                  # desired-profile tolerance, k-Wave options (non-defaults)
        self.ring_options: dict = {}           # fit_ring_profile / plane-scan options (non-defaults)
        self.auto_depths: str | None = None    # desired-depths text this page filled in itself
        self.has_fit = False
        self._busy = False
        # MATLAB starts only on an action that needs it (Load, Predict, Evaluate, Fit, Save,
        # k-Wave); until then the page works from the library entry read in Python.
        self._device: tuple[str, str, str] | None = None   # (serial, ds, library_path) shown
        self._record_from: tuple[str, str] | None = None   # device whose meta.fit settings are applied
        self.auto_csv: str | None = None       # CSV path this page filled in itself
        self.entry: dict | None = None         # parsed library entry of the selected device

        outer = QVBoxLayout(self)
        outer.setContentsMargins(12, 8, 12, 8)
        top = QHBoxLayout()
        top.addWidget(QLabel("Transducer"))
        self.tx = QComboBox()
        self.tx.setMinimumWidth(260)
        self.tx.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
        self.tx.currentIndexChanged.connect(lambda _i: self._select())
        top.addWidget(self.tx)
        # any library transducer and driving system, without a study (calibration_gui's lists)
        self.lib_serial = QComboBox()
        self.lib_serial.setMinimumWidth(200)
        self.lib_ds = QComboBox()
        self.lib_ds.setMinimumWidth(160)
        self.lib_serial.currentIndexChanged.connect(lambda _i: self._library_serial_changed())
        self.lib_ds.currentIndexChanged.connect(lambda _i: self.refresh())
        top.addWidget(self.lib_serial)
        top.addWidget(self.lib_ds)
        self.any_lib = QCheckBox("Any library transducer")
        self.any_lib.setToolTip("Calibrate any transducer and driving system of the equipment library, not "
                                "only this study's (as calibration_gui); entries are saved to a folder you "
                                "choose, without changing the config.")
        self.any_lib.setChecked(bool(window.session.get("calibration_any_library", False)))
        self.any_lib.toggled.connect(self._toggle_any_library)
        top.addWidget(self.any_lib)
        geo = QPushButton("View geometry")
        geo.setToolTip("Side and face view of the selected transducer (Equipment sub-tab)")
        geo.clicked.connect(self._view_geometry)
        top.addWidget(geo)
        self.tx_files = FileBadge()             # the device and library entry files, on hover
        top.addWidget(self.tx_files)
        top.addStretch(1)
        outer.addLayout(top)
        # the entry note and MATLAB's state on their own row: the selectors above keep their
        # full width (the combos never shrink below their longest item)
        info = QHBoxLayout()
        self.coverage = QLabel()
        self.coverage.setObjectName("CardSub")
        self.coverage.setWordWrap(True)
        info.addWidget(self.coverage, 1)
        bridge = window.bridge()
        self.status_widget = MatlabStatus(bridge)
        info.addWidget(self.status_widget)
        bridge.state_changed.connect(self._idle_status)     # after MatlabStatus._update
        self._idle_status(getattr(bridge, "state", ""))
        outer.addLayout(info)

        # plots on top at full width (they scroll rather than squash when short), steps below
        split = QSplitter(Qt.Orientation.Vertical)
        outer.addWidget(split, 1)
        top_scroll = QScrollArea()
        top_scroll.setWidgetResizable(True)
        top_scroll.setFrameShape(QFrame.Shape.NoFrame)
        top_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        top_scroll.setWidget(self._results())
        split.addWidget(top_scroll)
        split.addWidget(self._controls())
        split.setStretchFactor(0, 1)
        split.setStretchFactor(1, 0)
        split.setChildrenCollapsible(False)
        self._split = split
        self._split_sized = False
        self._size_timer = QTimer(self, singleShot=True, interval=0)
        self._size_timer.timeout.connect(self._size_split)
        self._entry_timer = QTimer(self, singleShot=True, interval=300)
        self._entry_timer.timeout.connect(self._draw_entry)
        self.depths.textChanged.connect(lambda _t: self._entry_timer.start())
        self._update_buttons()

    # ── layout ──────────────────────────────────────────────────────────
    def showEvent(self, ev) -> None:
        super().showEvent(ev)
        if not self._split_sized:                   # once the width (row or stack) is known
            self._size_timer.start()

    def _size_split(self) -> None:
        """The steps get the height they need (at most 45 %), the plots the rest."""
        total = sum(self._split.sizes())
        if total <= 0 or not self.isVisible():
            return
        self._split_sized = True
        need = self._controls_body.sizeHint().height() + 4
        low = min(need, int(0.45 * total))
        self._split.setSizes([total - low, low])

    def _controls(self) -> QWidget:
        """Three numbered steps below the plots in two columns: Data with the log under it, then
        Fit and Save. The fit options most users keep at their defaults sit in a collapsed
        "Fit options" section, the rest under Advanced…"""
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        body = QWidget()
        cols = QHBoxLayout(body)
        cols.setContentsMargins(0, 6, 6, 0)
        cols.setSpacing(12)
        left, right = QVBoxLayout(), QVBoxLayout()
        for c in (left, right):
            c.setSpacing(8)
            cols.addLayout(c, 1)
        scroll.setWidget(body)

        # 1 Data
        box, v = _step(1, "Data", "hydrophone profiles and the depths you need")
        left.addWidget(box)
        form = _form()
        v.addLayout(form)
        csv_row = QHBoxLayout()
        self.csv = QLineEdit()
        self.csv.setPlaceholderText("hydrophone axial-profile CSV")
        self.csv.textChanged.connect(lambda _t: self._update_buttons())
        csv_row.addWidget(self.csv, 1)
        b = icons.PathButton("file", "Choose a hydrophone axial-profile CSV")
        b.clicked.connect(self._browse_csv)
        csv_row.addWidget(b)
        self.load_btn = QPushButton("Load")
        self.load_btn.clicked.connect(self._load_csv)
        csv_row.addWidget(self.load_btn)
        form.addRow("Profile CSV", csv_row)

        self.measured = QListWidget()
        self.measured.setMaximumHeight(96)
        self.measured.setFlow(QListWidget.Flow.LeftToRight)
        self.measured.setWrapping(True)
        self.measured.setResizeMode(QListWidget.ResizeMode.Adjust)
        fm = self.measured.fontMetrics()
        self.measured.setGridSize(QSize(fm.horizontalAdvance("000.00") + 36, fm.height() + 6))
        self.measured.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.measured.setToolTip("Measured depths to fit (Calibrate at: measured depths); untick to leave "
                                 "a depth out")
        form.addRow("Measured depths (mm)", self.measured)

        d_row = QHBoxLayout()
        self.depths = QLineEdit()
        self.depths.setPlaceholderText("20:5:60 or 30, 42.5, 50")
        self.depths.setToolTip("Depths from the exit plane")
        self.depths.textChanged.connect(lambda _t: self._update_buttons())
        d_row.addWidget(self.depths, 1)
        sd = QPushButton("Study")
        sd.setToolTip("The depths this study uses with this device: focal depth settings (config and "
                      "manifest) and the target depth of each placement (exit plane → focus_pos)")
        sd.clicked.connect(self.use_study_depths)
        d_row.addWidget(sd)
        self.predict_btn = QPushButton("Predict")
        self.predict_btn.setToolTip("Profiles at these depths (measured, or interpolated from the measured "
                                    "ones) and the geometric, library, and fitted models as a study "
                                    "resolves them")
        self.predict_btn.clicked.connect(self.predict)
        d_row.addWidget(self.predict_btn)
        form.addRow("Desired depths (mm)", d_row)

        # 2 Fit
        box, v = _step(2, "Fit")
        right.addWidget(box)
        form = _form()
        v.addLayout(form)
        self.calat = QComboBox()
        self.calat.addItem("Measured depths", "measured")
        self.calat.addItem("Desired depths", "desired")
        self.calat.setToolTip("Fit at the ticked measured depths, or at the desired depths (profiles "
                              "interpolated from the measured ones)")
        form.addRow("Calibrate at", self.calat)
        self.isppa = _spin(30, 0.1, 1000, 1)
        self.isppa.setSuffix(" W/cm²")
        self.isppa.setToolTip("Reference I<sub>SPPA</sub> of the stored drives; defaults to transducer.target_isppa_wcm2")
        form.addRow("Reference I<sub>SPPA</sub>", self.isppa)
        a_row = QHBoxLayout()
        self.fit_btn = QPushButton("Fit")
        self.fit_btn.setProperty("primary", True)
        self.fit_btn.clicked.connect(self.fit)
        a_row.addWidget(self.fit_btn, 1)
        self.eval_btn = QPushButton("Evaluate")
        self.eval_btn.setToolTip("Geometric, library, and fit at the measured depths (R² against the hydrophone)")
        self.eval_btn.clicked.connect(self.evaluate_measured)
        a_row.addWidget(self.eval_btn, 1)
        self.opt_btn = QToolButton(text="Fit options ▸", checkable=True)
        self.opt_btn.setToolTip("Ring profile, free parameters, forward model, loss, and starts. The defaults "
                                "suit most transducers.")
        a_row.addWidget(self.opt_btn)
        v.addLayout(a_row)

        self.options = QWidget()
        self.options.setVisible(False)
        self.opt_btn.toggled.connect(self._toggle_options)
        v.addWidget(self.options)
        form = _form()
        self.options.setLayout(form)
        self.ring = QComboBox()
        self.ring.addItem("keep the entry's", "auto")
        self.ring.addItem("fit (shared by all depths)", "fit")
        self.ring.addItem("uniform rings", "none")
        self.ring.setToolTip(
            "How each ring vibrates across its width (fit_ring_profile). Fit: every ring split into sub-rings "
            "with weights shared by all measured depths, then phases, gains, and drive per depth; the near field "
            "from 4 mm is then fitted too. Improves transducers whose near-field / focal peak ratio uniform rings "
            "cannot reproduce (e.g. CTX). Axial profiles alone do not determine the off-axis field: add plane "
            "scans under Options… where available.")
        r_row = QHBoxLayout()
        r_row.addWidget(self.ring, 1)
        self.rpo_btn = QPushButton("Options…")
        self.rpo_btn.setToolTip("Sub-rings, smoothness, optimiser, and plane scans (off-axis data) for the "
                                "ring-profile fit")
        self.rpo_btn.clicked.connect(self._ring_options)
        r_row.addWidget(self.rpo_btn)
        form.addRow("Ring profile", r_row)
        fp = QHBoxLayout()
        self.fp_phase = QCheckBox("phase", checked=True)
        self.fp_amp = QCheckBox("gain", checked=True)
        for w in (self.fp_phase, self.fp_amp):
            fp.addWidget(w)
        fp.addStretch(1)
        form.addRow("Free parameters", fp)
        self.fwd = QComboBox()
        self.fwd.addItems(["rayleigh", "oneil"])
        form.addRow("Forward model", self.fwd)
        self.skip = _spin(10, 0, 50, 1)
        self.skip.setSuffix(" mm")
        self.skip.setToolTip("Profile samples closer to the exit plane than this are left out of the loss")
        form.addRow("Near-field skip", self.skip)
        self.domain = QComboBox()
        self.domain.addItems(["intensity", "amplitude"])
        form.addRow("Loss on", self.domain)
        self.huber = _spin(0.8, 0.001, 100, 3)
        self.huber.setToolTip("Pseudo-Huber scale on the peak-normalised profile; large (e.g. 10) = least squares")
        form.addRow("Huber δ", self.huber)
        self.nstarts = QSpinBox(minimum=1, maximum=200, value=16)
        self.nstarts.setToolTip("Random starts per depth")
        form.addRow("Starts", self.nstarts)
        p_row = QHBoxLayout()
        self.align = QComboBox()
        self.align.addItems(["label", "peak"])
        self.align.setToolTip("Focus alignment when interpolating between measured depths: column "
                              "focal-depth labels (default) or intensity peaks beyond the near field")
        p_row.addWidget(self.align)
        p_row.addWidget(QLabel("extrap. ≤"))
        self.extrap = _spin(0, 0, 50, 1)
        self.extrap.setSuffix(" mm")
        self.extrap.setToolTip("Desired depths up to this far outside the measured range get the nearest "
                               "measured profile shifted; 0 = no extrapolation")
        p_row.addWidget(self.extrap)
        form.addRow("Interpolation", p_row)
        self.adv_btn = QPushButton("Advanced…")
        self.adv_btn.setToolTip("Phase constraints, absolute-intensity fit, global optimisers, focal weighting, "
                                "near-field penalty, presets")
        self.adv_btn.clicked.connect(self._advanced)
        form.addRow("", self.adv_btn)

        # 3 Save
        box, v = _step(3, "Save")
        right.addWidget(box)
        self.lib_label = QLabel()
        self.lib_label.setWordWrap(True)
        self.lib_label.setObjectName("CardSub")
        self.lib_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        v.addWidget(self.lib_label)
        self.lib_folder_btn = QPushButton("Folder…")
        self.lib_folder_btn.setToolTip("Folder of your own library for the entries of library transducers")
        self.lib_folder_btn.clicked.connect(self._choose_library_folder)
        v.addWidget(self.lib_folder_btn, 0, Qt.AlignmentFlag.AlignLeft)
        s_row = QHBoxLayout()
        self.save_btn = QPushButton("Save entry")
        self.save_btn.clicked.connect(self.save)
        s_row.addWidget(self.save_btn, 1)
        self.export_btn = QPushButton("Export")
        menu = QMenu(self.export_btn)
        menu.addAction("Profiles (CSV)…", self._export_profiles)
        menu.addAction("Fit settings (CSV)…", self._export_settings)
        self.export_btn.setMenu(menu)
        s_row.addWidget(self.export_btn, 1)
        self.report_btn = QPushButton("Report")
        self.report_btn.clicked.connect(self.open_report)
        s_row.addWidget(self.report_btn, 1)
        v.addLayout(s_row)
        # grid choice on its own row: beside the button, the pair was the widest row of the step
        k_row = QHBoxLayout()
        k_row.addWidget(QLabel("k-Wave grid"))
        self.kw_grid = QComboBox()
        self.kw_grid.addItems(["3D", "axisymmetric"])
        self.kw_grid.setToolTip("Free-water grid, unless the study config sets grid.axisymmetric")
        k_row.addWidget(self.kw_grid)
        k_row.addStretch(1)
        v.addLayout(k_row)
        self.kw_btn = QPushButton("k-Wave drive check…")
        self.kw_btn.clicked.connect(self.kwave_check)
        v.addWidget(self.kw_btn)
        self._tips = {
            self.load_btn: "Read the profiles in MATLAB (starts MATLAB on first use)",
            self.predict_btn: self.predict_btn.toolTip(),
            self.eval_btn: self.eval_btn.toolTip(),
            self.fit_btn: "Fit phases, gains, and drive per depth to the hydrophone profiles",
            self.save_btn: "Write the fitted entry to transducer.library_path (with meta.fit)",
            self.export_btn: "Export the shown profiles or the fit settings as CSV",
            self.report_btn: "Open the HTML report of the saved library entry (generate_calibration_report): "
                             "provenance, fit quality, the settings a study gets, and the amplitude lookup",
            self.kw_btn: "Per-depth fits: simulate each refined depth in free-water k-Wave with the "
                         "study config's grid and platform, correct the drive to the simulated ISPPA, "
                         "and save the entry. Runs on the Run tab (minutes per depth).",
        }

        # config sections that belong to calibration (e.g. the focal-depth reference), placed
        # here by the Transducer tab (TabPage.section_columns) rather than on its Settings sub-tab
        self.settings_column = QVBoxLayout()
        self.settings_column.setSpacing(8)
        self.log = QPlainTextEdit(readOnly=True)
        self.log.setMinimumHeight(90)
        self.log.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Ignored)   # leftover height only
        self.log.setFont(QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont))
        self.log.setPlaceholderText("Messages from the calibration steps appear here.")
        right.addLayout(self.settings_column)
        right.addStretch(1)
        left.addWidget(self.log, 1)
        self._controls_body = body
        return scroll

    def _toggle_options(self, on: bool) -> None:
        self.options.setVisible(on)
        self.opt_btn.setText("Fit options ▾" if on else "Fit options ▸")

    def _results(self) -> QWidget:
        """Two views: the library entry (CalibrationCoverage, no MATLAB needed) and, once
        MATLAB has answered, the profiles with a clickable per-depth R² strip."""
        w = QWidget()
        outer = QVBoxLayout(w)
        outer.setContentsMargins(0, 0, 0, 0)
        bar = QHBoxLayout()
        bar.setSpacing(0)
        self.view_btns: dict[str, QToolButton] = {}
        for key, text, tip in (("entry", "Library entry", "What the current library entry holds (no MATLAB)"),
                               ("profiles", "Profiles", "Models against the hydrophone profiles (MATLAB)")):
            b = QToolButton(text=text, checkable=True)
            b.setToolTip(tip)
            b.clicked.connect(lambda _c, k=key: self._show_view(k))
            bar.addWidget(b)
            self.view_btns[key] = b
        bar.addSpacing(12)
        self.depth_label = QLabel("Show depth")
        bar.addWidget(self.depth_label)
        self.show_depth = QComboBox()                   # kept in sync with the R² strip (hidden)
        self.show_depth.currentIndexChanged.connect(lambda _i: self._depth_changed())
        bar.addWidget(self.show_depth)
        self.depth_label.hide()
        self.show_depth.hide()
        bar.addStretch(1)
        outer.addLayout(bar)

        self.views = QStackedWidget()
        outer.addWidget(self.views, 1)
        page = QWidget()
        v = QVBoxLayout(page)
        v.setContentsMargins(0, 6, 0, 0)
        self.hint = QLabel()
        self.hint.setWordWrap(True)
        self.hint.setObjectName("TabIntro")
        v.addWidget(self.hint)
        self.entry_plot = cplots.CalibrationCoverage(height=230)
        v.addWidget(self.entry_plot, 1)
        self.views.addWidget(page)

        page = QWidget()
        v = QVBoxLayout(page)
        v.setContentsMargins(0, 6, 0, 0)
        self.strip_title = QLabel()
        self.strip_title.setObjectName("CardSub")
        v.addWidget(self.strip_title)
        self.strip = cplots.R2Strip()
        self.strip.selected.connect(self._strip_selected)
        v.addWidget(self.strip)
        if FigureCanvasQTAgg is not None:
            self.fig = Figure(figsize=(7, 6), layout="constrained")
            self.ax_prof, self.ax_r2 = self.fig.subplots(2, 1, height_ratios=[1.3, 1])
            self.canvas = FigureCanvasQTAgg(self.fig, export_name="calibration_profiles")
            self.canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
            self.canvas.setMinimumHeight(400)       # two axes with legends: below this they overlap
            v.addWidget(self.canvas, 3)
        else:
            self.fig = self.canvas = None
            v.addWidget(QLabel("Install matplotlib to see the profiles."), 3)
        self.table = QTableWidget()
        self.table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.table.setMinimumHeight(120)
        v.addWidget(self.table, 1)
        self.views.addWidget(page)
        self._restyle_views()
        self._show_view("entry")
        return w

    def _restyle_views(self) -> None:
        css = cplots._segmented_css(style.active())
        keys = list(self.view_btns)
        for i, k in enumerate(keys):
            r = "border-top-left-radius: 5px; border-bottom-left-radius: 5px;" if i == 0 else \
                "border-top-right-radius: 5px; border-bottom-right-radius: 5px;"
            self.view_btns[k].setStyleSheet(css.replace("border-radius: 0;", "border-radius: 0; " + r))

    def _show_view(self, key: str) -> None:
        if key == "profiles" and self.view is None:
            key = "entry"
        self.views.setCurrentIndex(0 if key == "entry" else 1)
        for k, b in self.view_btns.items():
            b.setChecked(k == key)
        self.view_btns["profiles"].setEnabled(self.view is not None)
        self.view_btns["profiles"].setToolTip(
            "Models against the hydrophone profiles" if self.view is not None else
            "Available after Predict, Evaluate, or Fit (they start MATLAB)")

    def _depth_changed(self) -> None:
        sel = self.show_depth.currentText()
        try:
            self.strip.select(None if sel in ("All depths", "") else float(sel))
        except ValueError:  # best-effort: not a depth: the selection stays
            pass
        self._draw()

    def _strip_selected(self, depth: float) -> None:
        if math.isnan(depth):
            i = self.show_depth.findText("All depths")
        else:
            i = next((j for j in range(self.show_depth.count())
                      if self.show_depth.itemText(j) not in ("All depths", "")
                      and abs(float(self.show_depth.itemText(j)) - depth) < 0.01), -1)
        if i >= 0:
            self.show_depth.setCurrentIndex(i)          # draws via _depth_changed

    def _idle_status(self, state: str) -> None:
        """While MATLAB is off, say which actions start it."""
        if state == "off":
            c = style.active()
            self.status_widget.label.setText(f"<span style='color:{c['ink3']}'>●</span> MATLAB starts when "
                                             "you load, predict, or fit")

    def _update_buttons(self) -> None:
        """Enable the step buttons for what is possible now; a disabled button says why."""
        if not hasattr(self, "_tips"):
            return
        busy = self._busy
        has_csv = bool(self.csv.text().strip()) or bool(self.info.get("csv"))
        has_depths = bool(self.depths.text().strip())
        why = {
            self.load_btn: None if has_csv else "Enter or choose a profile CSV first",
            self.predict_btn: None if has_depths else "Enter desired depths first (or click Study)",
            self.eval_btn: None if has_csv or has_depths else
            "Load a hydrophone profile CSV, or enter desired depths to predict",
            self.fit_btn: None if has_csv else "Load a hydrophone profile CSV first (step 1)",
            self.save_btn: None if self.has_fit else "Run Fit first: there is nothing to save yet",
            self.export_btn: None if self.view is not None else
            "Predict, Evaluate, or Fit first: there are no profiles to export yet",
            self.kw_btn: None if self.has_fit else "Needs a per-depth fit: run Fit first",
            self.report_btn: None if self.entry is not None or self.info.get("library_file") else
            "No library entry yet: save a fit first",
        }
        for b, reason in why.items():
            if busy:
                reason = "MATLAB is working…"
            b.setEnabled(reason is None)
            b.setToolTip(self._tips[b] if reason is None else f"{reason}.\n\n{self._tips[b]}")
        self.tx.setEnabled(not busy)

    # ── state ───────────────────────────────────────────────────────────
    @property
    def model(self):
        return self.window.model

    def _say(self, msg: str) -> None:
        self.log.appendPlainText(f"[{_dt.datetime.now():%H:%M:%S}] {msg}")

    def _set_busy(self, busy: bool, msg: str = "") -> None:
        self._busy = busy
        self._update_buttons()
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

    def library_mode(self) -> bool:
        """Any library transducer (not the study's): ti is -1 and nothing reads the config."""
        return self.any_lib.isChecked()

    def current(self) -> tuple[int, str, str] | None:
        if self.library_mode():
            serial = self.lib_serial.currentData()
            return (-1, str(serial), str(self.lib_ds.currentData() or "")) if serial else None
        data = self.tx.currentData()
        return tuple(data) if data else None

    def library_path(self) -> str:
        if self.library_mode():
            return str(self.window.session.get("calibration_library_out", "") or "")
        cur = self.current()
        lib = self.model.get(f"{TX}.library_path", cur[0] if cur else 0)
        return "" if is_empty(lib) else str(lib)

    def _toggle_any_library(self, on: bool) -> None:
        self.window.session.set("calibration_any_library", on)
        self.refresh()

    def _library_serial_changed(self) -> None:
        """The driving systems of the chosen library transducer (and none)."""
        serial = self.lib_serial.currentData()
        dev = self.model.library.get(serial) if (self.model.library and serial) else None
        keep = self.lib_ds.currentData()
        self.lib_ds.blockSignals(True)
        self.lib_ds.clear()
        self.lib_ds.addItem("(no driving system)", "")
        for ds in (dev.combos if dev is not None else {}):
            self.lib_ds.addItem(ds, ds)
        idx = self.lib_ds.findData(keep)
        self.lib_ds.setCurrentIndex(idx if idx >= 0 else (1 if self.lib_ds.count() > 1 else 0))
        self.lib_ds.blockSignals(False)
        self.refresh()

    def _view_geometry(self) -> None:
        """The selected transducer on the Equipment sub-tab (calibration_gui's View geometry)."""
        from .main_window import EQUIPMENT
        cur = self.current()
        tab = self.window.tab_pages.get("transducer")
        if tab is not None and not tab.advanced.isChecked():
            tab.advanced.setChecked(True)               # the Equipment sub-tab is an advanced one
        page = self.window.show_special(EQUIPMENT)
        page.refresh()
        page.keep_serial.setChecked(False)
        if cur is not None:
            idx = page.start.findData(cur[1])
            if idx >= 0:
                page.start.setCurrentIndex(idx)

    def _choose_library_folder(self) -> None:
        start = self.library_path() or str(Path.home())
        path = QFileDialog.getExistingDirectory(self, "Folder of your own library", start)
        if path:
            self.window.session.set("calibration_library_out", path)
            self.refresh()

    def _fill_library_lists(self) -> None:
        lib = self.model.library
        keep = self.lib_serial.currentData()
        self.lib_serial.blockSignals(True)
        self.lib_serial.clear()
        for s, d in (lib.devices.items() if lib else []):
            if d.type == "annular" and not d.excluded:
                self.lib_serial.addItem(s, s)
        fit_combo_to_items(self.lib_serial, 200)
        idx = self.lib_serial.findData(keep)
        self.lib_serial.setCurrentIndex(max(idx, 0))
        self.lib_serial.blockSignals(False)
        if self.lib_ds.count() == 0 or keep != self.lib_serial.currentData():
            self.lib_ds.blockSignals(True)
            self.lib_ds.clear()
            dev = lib.get(self.lib_serial.currentData()) if lib else None
            self.lib_ds.addItem("(no driving system)", "")
            for ds in (dev.combos if dev is not None else {}):
                self.lib_ds.addItem(ds, ds)
            self.lib_ds.setCurrentIndex(1 if self.lib_ds.count() > 1 else 0)
            self.lib_ds.blockSignals(False)

    def refresh(self) -> None:
        """Called when the sub-tab is shown: follow the config's transducers and show the
        selected one's library entry. MATLAB is not started here; once a user action has
        opened the device in MATLAB, a change of the study's depths is predicted again."""
        lib_mode = self.library_mode()
        for w in (self.lib_serial, self.lib_ds, self.lib_folder_btn):
            w.setVisible(lib_mode)
        self.tx.setVisible(not lib_mode)
        if lib_mode:
            self._fill_library_lists()
        items = self._transducers() if not lib_mode else (
            [self.current()] if self.current() is not None else [])
        keep = self.current() if not lib_mode else None
        self.tx.blockSignals(True)
        self.tx.clear()
        for ti, serial, ds in (items if not lib_mode else []):
            dev = self.model.library.get(serial) if self.model.library else None
            kind = f" ({dev.type})" if dev is not None and dev.type else ""
            self.tx.addItem(f"Transducer {ti + 1}: {serial}{' / ' + ds if ds else ''}{kind}", [ti, serial, ds])
        fit_combo_to_items(self.tx, 260)
        if keep is not None:
            idx = self.tx.findData(list(keep))
            if idx >= 0:
                self.tx.setCurrentIndex(idx)
        self.tx.blockSignals(False)
        lib = self.library_path()
        if lib_mode:
            self.lib_label.setText(f"Save to <tt>{lib}</tt> (your own library; the config is not changed)" if lib
                                   else f"<span style='color:{style.active()['amber']}'>Choose the folder of your "
                                   "own library (Folder…); Save asks for one.</span>")
        else:
            self.lib_label.setText(f"Save to <tt>{lib}</tt>" if lib else
                                   f"<span style='color:{style.active()['amber']}'>transducer.library_path is "
                                   "not set: Save asks for a folder and sets it.</span>")
        if not items:
            self.tx_files.set_files([])
            self.coverage.setText("No transducer from the equipment library in this config "
                                  "(set transducer.serial on the Settings tab).")
            self._device = None
            self.entry = None
            self.entry_plot.set_entry(None, empty_message="No transducer from the equipment library in this "
                                      "config: set transducer.serial on the Settings tab.")
            self.hint.setText("")
            self._show_view("entry")
            return
        changed = self._select()                # the config or placement changed the depths
        cur = self.current()
        if not self._busy and cur is not None and self.opened == (cur[1], cur[2], self.library_path()):
            if changed:
                self.predict()

    def _select(self) -> bool:
        """The selected device, from the library alone (no MATLAB): its entry, the fit
        settings recorded in it, its hydrophone CSV and the study's depths. True when the
        desired depths changed."""
        cur = self.current()
        if cur is None:
            return False
        key = (cur[1], cur[2], self.library_path())
        if key != self._device:
            if self._device is None or self._device[:2] != key[:2]:
                self.auto_depths = None         # another device: the old depths do not apply
                self.depths.clear()
                if self.csv.text().strip() == (self.auto_csv or ""):
                    self.csv.clear()
                self.auto_csv = None
                self.info, self.view, self.has_fit = {}, None, False
                self.measured.clear()
                target = self.model.get(f"{TX}.target_isppa_wcm2", cur[0]) if cur[0] >= 0 else None
                if isinstance(target, (int, float)) and not isinstance(target, bool) and target > 0:
                    self.isppa.setValue(float(target))
            self._device = key
            self.entry, path = cplots.load_entry(self.model.library, cur[1], cur[2], key[2])
            meta = (self.entry or {}).get("meta") or {}
            desc = (f"Library entry: {meta.get('calibration_type') or '?'}"
                    f"{', ' + str(meta['last_updated']) if meta.get('last_updated') else ''}"
                    if path else "No library entry yet")
            self.coverage.setText(desc)             # the selector names the device
            dev = self.model.library.get(cur[1]) if self.model.library else None
            self.tx_files.set_files([("Device", dev.file if dev is not None else ""),
                                     ("Library entry", path or "none yet")])
            rec = meta.get("fit")
            if isinstance(rec, dict) and self._record_from != key[:2]:
                s, _depths = cal.settings_from_record(rec, self.settings())
                self.apply_settings(s)
                self._record_from = key[:2]
                self._say(f"Fit settings taken from the entry (meta.fit: {rec.get('mode', '')}, "
                          f"{rec.get('tool', '')}).")
            csv = cplots.hydrophone_csv(self.entry, self.window.prestus_root,
                                        getattr(self.window.model, "library", None))
            if csv is not None and not self.csv.text().strip():
                self.csv.setText(str(csv))
                self.auto_csv = str(csv)
                if not self.measured.count():
                    self._fill_measured(cplots.csv_focal_depths(csv))
            if self.opened is None or self.opened != key:
                self._show_view("entry")
        changed = self._prefill_depths()
        self._draw_entry()
        self._update_buttons()
        return changed

    def _draw_entry(self) -> None:
        """The library entry against the desired depths, and what to do next."""
        cur = self.current()
        if cur is None:
            return
        depths = cal.parse_depths(self.depths.text()) or []
        name = f"{cur[1]}{' with ' + cur[2] if cur[2] else ''}"
        self.entry_plot.set_entry(self.entry, depths, name,
                                  f"No calibration entry for {name}: an annular transducer needs one (or phases and amplitude set inline) before it can run; fit one below.")
        if self.csv.text().strip():
            nxt = ("Next: <b>Evaluate</b> compares the geometric, library, and fitted models with the "
                   "hydrophone profiles; <b>Fit</b> makes a new entry. Both start MATLAB.")
        else:
            nxt = ("Load a hydrophone profile CSV (step 1) to fit a new entry, or <b>Predict</b> the "
                   "profiles at the desired depths (starts MATLAB).")
        self.hint.setText(nxt)

    def _open(self, then=None) -> None:
        cur = self.current()
        if cur is None:
            return
        key = (cur[1], cur[2], self.library_path())
        if key == self.opened and then is None:
            return
        self._set_busy(True, f"Opening {cur[1]}{' / ' + cur[2] if cur[2] else ''} in MATLAB…")
        self.window.bridge().call(
            "nexus_calibration",
            {"action": "open", "session": self.session_id, "serial": cur[1], "driving_system": cur[2],
             "library_path": key[2]},
            on_done=lambda out: self._opened(key, out, then), on_error=self._failed)

    def _ensure_open(self, then) -> bool:
        """True if the selected device is open in MATLAB; else open it (starting MATLAB) and
        call then() afterwards."""
        cur = self.current()
        if cur is None:
            return False
        if self.opened == (cur[1], cur[2], self.library_path()):
            return True
        self._open(then=then)
        return False

    def _opened(self, key, out: dict, then=None) -> None:
        self._set_busy(False)
        self.opened = key
        self.info = out
        self.has_fit = False
        self.view = None
        if not self.csv.text().strip() or self.csv.text().strip() == (self.auto_csv or ""):
            self.csv.setText(out.get("csv") or "")
            self.auto_csv = out.get("csv") or None
        self._fill_measured()
        entry = (f"Library entry: {out.get('calibration_type') or '?'}"
                 f"{', ' + out['last_updated'] if out.get('last_updated') else ''}"
                 if out.get("library_file") else "No library entry")
        self.coverage.setText(entry)                # the selector names the device
        if out.get("library_file"):
            self.tx_files.set_files([*(f for f in self.tx_files.files if f[0] != "Library entry"),
                                     ("Library entry", out["library_file"])])
        self._say(f"{out['combo']}: {entry}")
        self._update_buttons()
        if out.get("matrix"):
            self._say("Matrix transducer: geometry only. Prediction and fitting support annular arrays.")
            return
        rec = out.get("fit_record")
        fallback = as_floats(out.get("default_depths"))
        if isinstance(rec, dict):
            s, depths = cal.settings_from_record(rec, self.settings())
            fallback = depths or fallback
            if self._record_from != key[:2]:    # not yet applied from the YAML in _select
                self.apply_settings(s)
                self._record_from = key[:2]
                self._say(f"Fit settings taken from the entry (meta.fit: {rec.get('mode', '')}, "
                          f"{rec.get('tool', '')}).")
        self._prefill_depths(fallback)
        if then is not None:
            then()
        elif out.get("csv"):
            self.evaluate_measured()
        else:
            self.predict()

    def _fill_measured(self, depths: list[float] | None = None) -> None:
        """The measured depths (MATLAB's answer, else depths read from the CSV header); depths
        the user unticked stay unticked when the same depths come back."""
        off = {round(self.measured.item(i).data(Qt.ItemDataRole.UserRole), 2) for i in range(self.measured.count())
               if self.measured.item(i).checkState() != Qt.CheckState.Checked}
        self.measured.clear()
        for d in (as_floats(self.info.get("measured_depths")) if depths is None else depths):
            it = QListWidgetItem(f"{d:.2f}")
            it.setData(Qt.ItemDataRole.UserRole, d)
            it.setFlags(it.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            it.setCheckState(Qt.CheckState.Unchecked if round(d, 2) in off else Qt.CheckState.Checked)
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
            calibrate_at=self.calat.currentData(), advanced=dict(self.advanced), extra=dict(self.extra),
            ring_options=dict(self.ring_options))

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
        self.extra = dict(s.extra)
        self.ring_options = dict(s.ring_options)
        self._update_advanced_button()

    def _advanced(self) -> None:
        dlg = AdvancedDialog(self.advanced, self.extra, self)
        if dlg.exec():
            if dlg.grid_dims() is None:
                self._say("k-Wave grid: enter auto, nx nz, or nx ny nz (or leave it blank); grid unchanged.")
            p = dlg.preset()
            if p:
                self.fp_phase.setChecked("phase" in p["FreeParams"])
                self.fp_amp.setChecked("amp" in p["FreeParams"])
                self.fwd.setCurrentText(p["ForwardModel"])
                self.domain.setCurrentText(p["Domain"])
                self.huber.setValue(float(p["HuberDelta"]))
                self.skip.setValue(float(p["NearFieldSkipMm"]))
                self._say(f"Preset {dlg.pre.currentText()} loaded.")
            self.advanced = dlg.values()
            ex = dlg.extra()
            if dlg.grid_dims() is None and "GridDims" in self.extra:
                ex["GridDims"] = self.extra["GridDims"]
            self.extra = ex
            self._update_advanced_button()

    def _ring_options(self) -> None:
        dlg = RingOptionsDialog(self.ring_options, self)
        if dlg.exec():
            self.ring_options = dlg.values()
            self._update_advanced_button()
            if self.ring_options.get("OffAxis") and self.ring.currentData() != "fit":
                self._say('Plane scans are used when Ring profile is set to "fit".')

    def _update_advanced_button(self) -> None:
        n = len(self.advanced) + len(self.extra)
        self.adv_btn.setText(f"Advanced ({n})" if n else "Advanced…")
        n = len(self.ring_options)
        self.rpo_btn.setText(f"Options ({n})" if n else "Options…")

    def _study_depths(self) -> tuple[list[float], str]:
        """Desired depths for the open device from the config, manifest and placement, and
        where they came from."""
        cur = self.current()
        if cur is None or cur[0] < 0:            # a library transducer: no study depths
            return [], ""
        settings = cal.study_depths(self.model, cur[1], cur[2], self._manifest())
        targets = cal.target_depths(self.model, cur[1], cur[2], self._target_depth)
        return cal.desired_depths(settings, targets), cal.describe_depths(settings, targets)

    def _target_depth(self, ti: int) -> float | None:
        from .main_window import VIEWER
        viewer = self.window.pages.get(VIEWER)
        try:
            return viewer.target_depth_ep(ti) if viewer is not None else None
        except Exception:  # best-effort: unreadable T1 etc.: no target depth
            log_current("calibration: cannot compute the target depth")
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

        self.auto_csv = None                    # chosen by the user now: keep it

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
        if not self._ensure_open(self.evaluate_measured):
            return
        if not self.info.get("csv"):
            self.predict()
            return
        self._call("view", {"predicted": False, "settings": self.settings().bridge_settings()}, self._show)

    def fit(self) -> None:
        if not self._ensure_open(self.fit):
            return
        if not self.info.get("csv"):
            self._say("Fitting needs measured profiles: load a CSV first.")
            return
        s = self.settings()
        args = {"settings": s.bridge_settings(), "free_params": s.free_params, "ring_profile": s.ring_profile,
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
        if not out_dir and self.library_mode():
            out_dir = QFileDialog.getExistingDirectory(self, "Folder of your own library", str(Path.home()))
            if not out_dir:
                return
            self.window.session.set("calibration_library_out", out_dir)
            self.opened = (cur[1], cur[2], out_dir)     # same MATLAB state; only the folder is new
            self.refresh()
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
            self.info["library_file"] = out["file"]     # Report has an entry now
            self._update_buttons()
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

    def open_report(self) -> None:
        """Write the saved library entry's HTML report to a temporary folder and open it."""
        def written(out):
            path = Path(str(out.get("file") or ""))
            self._say(f"Calibration report: {path}")
            if path.is_file():
                QDesktopServices.openUrl(QUrl.fromLocalFile(str(path)))
        self._call("report", {}, written, "Writing the report…")

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
        lib = self.library_path()
        if self.library_mode():
            if not lib:
                self._say("Choose the folder of your own library (Folder…) first: the corrected entry is "
                          "saved there.")
                return
            sim_config = self._library_sim_config()
            if sim_config is None:
                return
            folder = Path(lib) / "calibration"
        else:
            if self.model.path is None and not self.window.save():
                return
            if not lib:
                self._say("Set transducer.library_path (or Save the analytical fit once) first: the "
                          "corrected entry is saved there.")
                return
            sim_config = str(self.model.path)
            folder = self.model.path.parent / "calibration"
        combo = self.info.get("combo", "calibration")
        mat = folder / f"{combo}_fit.mat"
        script = folder / f"kwave_drive_{combo}.m"

        def written(_out):
            folder.mkdir(parents=True, exist_ok=True)
            ko = self.settings().kwave_options()
            axisym = self.kw_grid.currentText() == "axisymmetric"
            launcher.write_kwave_drive_script(
                script, mat, self.window.prestus_root, lib, sim_config or None, axisym,
                self.skip.value(), grid_dims=ko["GridDims"], apply_correction=ko["ApplyCorrection"])
            box = QMessageBox(self)
            box.setWindowTitle("k-Wave drive check")
            box.setText(f"Script written:\n{script}\n\nIt runs one free-water simulation per refined depth "
                        "(minutes each, with the study config's grid and platform) and saves the corrected "
                        f"entry to {lib}.\n\nSubmit it as a job (it waits for its simulations on a compute "
                        "node), run it here on the Run tab, or run the script yourself?")
            job = box.addButton("Submit as a job", QMessageBox.ButtonRole.AcceptRole)
            run = box.addButton("Run on the Run tab", QMessageBox.ButtonRole.ActionRole)
            box.addButton("Only write the script", QMessageBox.ButtonRole.RejectRole)
            box.exec()
            if box.clickedButton() is job:
                self._submit_kwave_drive(script, mat, lib, axisym, ko, sim_config)
            elif box.clickedButton() is run:
                page = self.window.show_special(RUN)
                page.run_script(f"k-Wave drive check ({combo})", script)
            else:
                self._say(f"k-Wave drive check script: {script}")
        self._call("write_fit", {"path": str(mat)}, written)

    def _submit_kwave_drive(self, script: Path, mat: Path, lib: str, axisym: bool, ko: dict,
                            sim_config: str = "") -> None:
        """Submit the drive check as a scheduler job from the Run tab: here, or over SSH with a
        copy of the script that uses the host's paths."""
        from .. import launcher
        from .. import remote as rmt
        from .main_window import RUN
        page = self.window.pages[RUN]
        remote = page.remote_host() if page._remote_on() else None
        if page._remote_on() and remote is None:
            return
        if remote is not None:
            try:
                cluster = script.with_name(f"{script.stem}_cluster.m")
                launcher.write_kwave_drive_script(
                    cluster, remote.to_remote(mat), remote.remote_root(self.window.prestus_root),
                    remote.to_remote(lib), remote.to_remote(sim_config) if sim_config else None, axisym,
                    self.skip.value(),
                    grid_dims=ko["GridDims"], apply_correction=ko["ApplyCorrection"])
            except rmt.PathMapError as e:
                self._say(f"Cannot submit over SSH: {e}.")
                return
            script = cluster
        self.window.show_special(RUN)
        if page.submit_kwave_drive(script, mat, Path(sim_config) if sim_config else None, remote):
            self._say(f"k-Wave drive check submitted as a job (Run tab): {script}")

    def _library_sim_config(self) -> str | None:
        """Simulation config for the k-Wave drive check of a library transducer (grid, platform,
        hpc; calibration_gui's SimConfig): asked once and remembered. '' = PRESTUS defaults;
        None = cancelled."""
        cfg = str(self.window.session.get("calibration_library_simconfig", "") or "")
        box = QMessageBox(self)
        box.setWindowTitle("k-Wave drive check")
        box.setText("Simulation config for the check (grid, platform and HPC settings; its transducer block "
                    f"is ignored):\n\n{cfg or 'none: PRESTUS defaults'}")
        keep = box.addButton("Use this", QMessageBox.ButtonRole.AcceptRole)
        choose = box.addButton("Choose…", QMessageBox.ButtonRole.ActionRole)
        none = box.addButton("PRESTUS defaults", QMessageBox.ButtonRole.ActionRole)
        box.addButton(QMessageBox.StandardButton.Cancel)
        box.exec()
        clicked = box.clickedButton()
        if clicked is keep:
            return cfg
        if clicked is none:
            self.window.session.set("calibration_library_simconfig", "")
            return ""
        if clicked is choose:
            path, _ = QFileDialog.getOpenFileName(self, "Simulation config", cfg or str(Path.home()),
                                                  "YAML (*.yaml *.yml)")
            if path:
                self.window.session.set("calibration_library_simconfig", path)
                return path
        return None

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
        self._fill_strip()
        self._show_view("profiles")
        self._update_buttons()
        self._draw()
        self._fill_table()

    def _fill_strip(self) -> None:
        """Per-depth R² of the best model shown (fit, else library, else geometric) as
        clickable bars; a click selects the depth in show_depth."""
        v = self.view or {}
        res = v.get("results") or {}
        kind = next((k for k in ("fit", "library", "geometric") if res.get(k) and as_floats(res[k].get("r2"))), None)
        depths = as_floats(v.get("depths"))
        r2 = []
        if kind is not None:
            by = dict(zip(as_floats(res[kind].get("depths")), as_floats(res[kind].get("r2"))))
            r2 = [next((y for x, y in by.items() if abs(x - d) < 0.01), float("nan")) for d in depths]
        sel = self.show_depth.currentText()
        cur = None if sel in ("All depths", "") else float(sel)
        name = {"fit": "the fit", "library": "the library entry", "geometric": "the geometric model"}.get(kind, "")
        self.strip.set_data(depths, r2, cur, allow_all=bool(v.get("predicted")), caption=name)
        self.strip_title.setText(
            f"R² of {name} per depth (green ≥ {cplots.R2_GOOD:g}, amber ≥ {cplots.R2_FAIR:g}, red below); "
            "click a bar to show that depth" if kind else "Click a depth to show its profile")

    def _draw(self) -> None:
        if self.canvas is None or self.view is None:
            return
        v = self.view
        c = style.active()
        ax, ax2 = self.ax_prof, self.ax_r2
        ax.clear()
        ax2.clear()
        self.fig.set_facecolor(c["surface"])
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
                        ax.plot(z, col, ls, color=c["ink"], lw=1.6, label=label if first else None)
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
            ax.axvline(x, ls=":", color=c["line2"], lw=0.8)
        ax.set_xlabel("Distance from exit plane [mm]")
        ax.set_ylabel("Normalised intensity")
        ax.set_title(f"{v.get('combo', '')} — {'all depths' if show_all else sel + ' mm'}"
                     f"{' (prediction, runtime library path)' if v.get('predicted') else ''}", fontsize=style.PLOT_PT["title"])
        ax.grid(True, alpha=0.3)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=style.PLOT_PT["normal"], loc="upper right")

        predicted = bool(v.get("predicted"))
        for kind, colour, name in STYLES:
            r = res.get(kind)
            if not r:
                continue
            y = as_floats(r.get("focus" if predicted else "r2"))
            ax2.plot(as_floats(r.get("depths")), y, "-o", color=colour, ms=4, lw=1.1, label=name)
        geo = v.get("dist_geom_mm")
        if isinstance(geo, (int, float)):
            ax2.axvline(geo, ls=":", color=c["ink4"], lw=0.8)
            ax2.annotate("natural focus", (geo, 0.02), xycoords=("data", "axes fraction"), fontsize=style.PLOT_PT["small"],
                         color=c["ink3"], rotation=90)
        if predicted and depths:
            lim = [min(depths), max(depths)]
            ax2.plot(lim, lim, "--", color=c["ink"], lw=0.9, label="Target")
            fm = as_floats(v.get("measured_focus"))
            for ptype, face, label in (("measured", c["ink"], "Hydrophone"),
                                       ("interpolated", "none", "Interpolated profile"),
                                       ("extrapolated", "none", "Extrapolated profile")):
                pts = [(depths[j], fm[j]) for j in range(min(len(fm), len(depths)))
                       if j < len(types) and types[j] == ptype and math.isfinite(fm[j])]
                if pts:
                    ax2.plot(*zip(*pts), "*", ms=9, mfc=face, mec=c["ink"], ls="none", label=label)
            span = as_floats(v.get("measured_span"))
            if len(span) == 2:
                ax2.axvspan(span[0], span[1], color=c["brand"], alpha=0.08, label="Measured range")
            ax2.set_xlabel("Target depth from exit plane [mm]")
            ax2.set_ylabel("Focus, FLHM centre [mm]")
        else:
            lo = min([0.0] + [y for k, *_ in STYLES if res.get(k) for y in as_floats(res[k].get("r2"))
                              if math.isfinite(y)])
            ax2.set_ylim(max(-0.2, lo), 1)
            ax2.set_xlabel("Focal depth from exit plane [mm]")
            ax2.set_ylabel("R² (profile shape)")
        ax2.grid(True, alpha=0.3)
        if ax2.get_legend_handles_labels()[0]:
            ax2.legend(fontsize=style.PLOT_PT["normal"], loc="best")
        _theme_axes((ax, ax2), c)
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
