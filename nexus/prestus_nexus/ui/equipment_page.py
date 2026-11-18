"""Transducer tab ▸ Equipment (advanced): view any library transducer, add annular or
matrix transducers, and add driving-system combinations to the equipment library, as
calibration_gui's 'View geometry' and 'New…' do (functions/gui/calibration_gui.m,
on_view_transducer / on_new_transducer). The fields, rules and file layout are in
..equipment_editor. MATLAB (nexus_equipment, through the bridge) computes matrix element
layouts, as the simulation places them, and checks the entry (check_equipment_config)
before it is saved. Entries go into your own library folders; maintainers may save into
the shipped config/equipment (calibration_gui's AllowShipped).
"""
from __future__ import annotations

import datetime as _dt
import json
import math
import tempfile
from pathlib import Path

from PyQt6.QtCore import QPointF, QRectF, Qt, QTimer
from PyQt6.QtGui import QColor, QPainter, QPen
from PyQt6.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout, QFrame,
                             QGridLayout, QHBoxLayout, QLabel, QLayout, QLineEdit, QMessageBox,
                             QPlainTextEdit, QPushButton, QScrollArea, QSizePolicy, QSpinBox, QSplitter,
                             QStackedWidget, QVBoxLayout, QWidget)

from .. import equipment_editor as ed
from .. import yamlio
from . import style
from .calibration_page import _step
from .transducer_view import TransducerSchematic

BLANK = "(blank)"
CHANNEL_COLOURS = ["#0072b2", "#e69f00", "#009e73", "#cc79a7", "#56b4e9", "#d55e00", "#f0e442", "#999999"]


def _form_layout() -> QFormLayout:
    """Not calibration_page._form: its wrapped rows are squeezed inside a scroll area."""
    f = QFormLayout()
    f.setContentsMargins(0, 0, 0, 0)
    f.setHorizontalSpacing(8)
    f.setVerticalSpacing(6)
    f.setRowWrapPolicy(QFormLayout.RowWrapPolicy.DontWrapRows)
    f.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
    return f


def _dspin(value: float, lo: float, hi: float, decimals: int = 2, suffix: str = "") -> QDoubleSpinBox:
    s = QDoubleSpinBox()
    s.setRange(lo, hi)
    s.setDecimals(decimals)
    s.setValue(value)
    s.setKeyboardTracking(False)
    if suffix:
        s.setSuffix(suffix)
    return s


def _ispin(value: int, lo: int, hi: int, suffix: str = "") -> QSpinBox:
    s = QSpinBox()
    s.setRange(lo, hi)
    s.setValue(value)
    s.setKeyboardTracking(False)
    if suffix:
        s.setSuffix(suffix)
    return s


def _row(*widgets) -> QWidget:
    w = QWidget()
    h = QHBoxLayout(w)
    h.setContentsMargins(0, 0, 0, 0)
    for x in widgets:
        h.addWidget(QLabel(x) if isinstance(x, str) else x)
    return w


class EditorSchematic(TransducerSchematic):
    """Side view of the transducer being edited (geometry from the form). While editing, the
    scale and the bowl apex stay put on screen, so changing the curvature radius or the exit
    plane distance moves only the marks that depend on it. The view refits when a transducer
    is loaded, the widget is resized, or the drawing would no longer fit."""

    def __init__(self, model, parent=None):
        super().__init__(model, parent)
        self.geo: dict | None = None
        self._view: tuple | None = None         # (W, H, s, x0) of the held frame

    def set_geometry(self, geo: dict | None) -> None:
        self.geo = geo
        self.update()

    def reset_frame(self) -> None:
        self._view = None

    def _frame(self, z_hi, z_max, a, W, H, on_screen: bool) -> tuple[float, float]:
        fit = self._fit(z_hi, z_max, a, W, H)
        if not on_screen:                       # exported figures are always fitted
            return fit
        if self._view is not None:
            vW, vH, s, x0 = self._view
            pad = self.PAD
            fits = (x0 + z_max * s <= W - pad / 2 and 2.4 * a * s <= H - 2 * pad - 18)
            # keep the frame unless the drawing overflows it or shrinks to under half of it
            if (vW, vH) == (W, H) and fits and fit[0] <= 2 * s:
                return s, x0
        self._view = (W, H, *fit)
        return fit

    def _geometry(self) -> dict | None:
        return self.geo


class FaceView(QWidget):
    """Face view: annular rings coloured by channel, or matrix element footprints, inside the
    aperture (as calibration_gui's face axes)."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(220, 220)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.rings: list[tuple[float, float, int]] = []      # (inner r, outer r, channel)
        self.elems: list[tuple[float, float]] = []           # matrix centres (x, y) [mm]
        self.elem_size = (0.0, 0.0, "rect")                  # width, height, shape
        self.aperture = 0.0                                  # radius [mm]
        self.note = ""

    def set_annular(self, rings, aperture) -> None:
        self.rings, self.elems, self.aperture, self.note = rings, [], aperture, ""
        self.update()

    def set_matrix(self, elems, w, h, shape, aperture, note="") -> None:
        self.rings, self.elems, self.elem_size, self.aperture, self.note = [], elems, (w, h, shape), aperture, note
        self.update()

    def paintEvent(self, _e) -> None:
        C = style.active()
        qp = QPainter(self)
        qp.setRenderHint(QPainter.RenderHint.Antialiasing)
        qp.fillRect(self.rect(), QColor(C["surface"]))
        if self.aperture <= 0:
            qp.setPen(QColor(C["ink3"]))
            qp.drawText(QRectF(self.rect()), Qt.AlignmentFlag.AlignCenter, self.note or "Face view")
            return
        side = min(self.width(), self.height()) - 20
        s = side / (2.3 * self.aperture)                       # px per mm
        c = QPointF(self.width() / 2, self.height() / 2)
        for ri, ro, ch in sorted(self.rings, key=lambda r: -r[1]):
            col = QColor(CHANNEL_COLOURS[ch % len(CHANNEL_COLOURS)])
            qp.setPen(QPen(QColor(C["ink3"]), 0.5))
            qp.setBrush(col)
            qp.drawEllipse(c, ro * s, ro * s)
            qp.setBrush(QColor(C["surface"]))
            if ri > 0:
                qp.drawEllipse(c, ri * s, ri * s)
        if self.elems:
            w, h, shape = self.elem_size
            qp.setPen(Qt.PenStyle.NoPen)
            qp.setBrush(QColor(CHANNEL_COLOURS[0]))
            if shape == "rect":
                for x, y in self.elems:
                    qp.drawRect(QRectF(c.x() + (x - w / 2) * s, c.y() - (y + h / 2) * s, w * s, h * s))
            else:
                d = math.sqrt(4 * w * h / math.pi)                 # equal-area diameter (create_matrix_karray)
                for x, y in self.elems:
                    qp.drawEllipse(QPointF(c.x() + x * s, c.y() - y * s), d / 2 * s, d / 2 * s)
        pen = QPen(QColor(C["ink2"]), 1, Qt.PenStyle.DashLine)
        qp.setPen(pen)
        qp.setBrush(Qt.BrushStyle.NoBrush)
        qp.drawEllipse(c, self.aperture * s, self.aperture * s)
        if self.note:
            qp.setPen(QColor(C["ink3"]))
            qp.drawText(QRectF(4, 4, self.width() - 8, 20), Qt.AlignmentFlag.AlignLeft, self.note)


class ComboRow(QFrame):
    """One driving-system combination: DS serial (existing or new), profile CSV (annular)."""

    def __init__(self, page: "EquipmentPage", ds_serials: list[str], combo: ed.Combo | None = None):
        super().__init__()
        self.page = page
        self.keep = dict(combo.keep) if combo else {}
        g = QGridLayout(self)
        g.setContentsMargins(0, 0, 0, 4)
        g.setHorizontalSpacing(4)
        self.ds = QComboBox()
        self.ds.setEditable(True)
        self.ds.addItems(ds_serials)
        self.ds.setCurrentText(combo.ds_serial if combo else "")
        self.ds.lineEdit().setPlaceholderText("DS serial")
        self.csv = QLineEdit(combo.csv if combo else "")
        self.csv.setPlaceholderText("Axial profile CSV (optional)")
        self.browse = QPushButton("Browse…")
        remove = QPushButton("Remove")
        g.addWidget(self.ds, 0, 0)
        g.addWidget(self.csv, 0, 1)
        g.addWidget(self.browse, 0, 2)
        g.addWidget(remove, 0, 3)
        g.setColumnStretch(1, 1)
        # a driving system not in the library: defined here, written as ds-<DS>_driver.yaml
        self.new_box = QWidget()
        nb = QHBoxLayout(self.new_box)
        nb.setContentsMargins(0, 0, 0, 0)
        nb.addWidget(QLabel("New driving system:"))
        self.ds_name = QLineEdit()
        self.ds_name.setPlaceholderText("name")
        self.ds_manu = QLineEdit()
        self.ds_manu.setPlaceholderText("manufacturer")
        self.ds_ch = _ispin(max(1, page.n_elem()), 1, 4096, " channels")
        for w in (self.ds_name, self.ds_manu, self.ds_ch):
            nb.addWidget(w)
        g.addWidget(self.new_box, 1, 0, 1, 4)
        self.ds.currentTextChanged.connect(self._changed)
        self.csv.textChanged.connect(page.update_view)
        for w in (self.ds_name, self.ds_manu):
            w.textChanged.connect(page.update_view)
        self.ds_ch.valueChanged.connect(page.update_view)
        self.browse.clicked.connect(self._browse)
        remove.clicked.connect(lambda: page.remove_combo(self))
        self._changed()

    def is_new(self) -> bool:
        s = self.ds.currentText().strip()
        lib = self.page.library()
        return bool(s) and (lib is None or s not in lib.driving_systems)

    def set_matrix(self, matrix: bool) -> None:
        """Profiles are annular only (calibration_gui disables the CSV for matrix entries)."""
        self.csv.setEnabled(not matrix)
        self.browse.setEnabled(not matrix)

    def _changed(self, *_):
        self.new_box.setVisible(self.is_new())
        self.page.update_view()

    def _browse(self) -> None:
        start = self.csv.text() or str(self.page.profile_root())
        path, _ = QFileDialog.getOpenFileName(self, "Axial profile CSV", start, "CSV (*.csv)")
        if path:
            self.csv.setText(path)

    def combo(self) -> ed.Combo:
        new = None
        if self.is_new():
            new = ed.NewDrivingSystem(self.ds.currentText().strip(), self.ds_name.text(),
                                      self.ds_manu.text(), self.ds_ch.value())
        return ed.Combo(self.ds.currentText().strip(), self.csv.text().strip(), new, dict(self.keep))


class EquipmentPage(QWidget):
    def __init__(self, window):
        super().__init__()
        self.window = window
        self.model = window.model
        self.base: dict | None = None            # raw entry being edited (Edit this entry)
        self.base_dev = None
        self.rows: list[ComboRow] = []
        self._filling = False
        self.mlayout_key: str | None = None       # matrix block the last layout answer is for
        self.mlayout: dict | None = None          # MATLAB element layout of that block
        self._layout_pending: str | None = None
        self._layout_timer = QTimer(self)
        self._layout_timer.setSingleShot(True)
        self._layout_timer.setInterval(400)      # while typing, ask MATLAB once the edits pause
        self._layout_timer.timeout.connect(self._request_layout)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 6, 0, 0)
        intro = QLabel("View any library transducer, or add an annular or matrix transducer and driving-system "
                       "combinations to your equipment library (as calibration_gui's View geometry and New…). Saved "
                       "into your own library folder, an entry with a shipped serial replaces the shipped one; saved "
                       "entries appear in the device lists and on the Calibration sub-tab when the folder is one of "
                       "your libraries.")
        intro.setObjectName("CardSub")
        intro.setWordWrap(True)
        outer.addWidget(intro)
        split = QSplitter(Qt.Orientation.Horizontal)
        outer.addWidget(split, 1)

        left = QWidget()
        self.left = left
        lv = QVBoxLayout(left)
        lv.setContentsMargins(0, 0, 6, 0)
        # keep the widget at least its layout's minimum: otherwise the scroll area squeezes the
        # form rows below their height instead of scrolling
        lv.setSizeConstraint(QLayout.SizeConstraint.SetMinimumSize)

        box, v = _step(1, "Start from", "a library transducer, to view or copy")
        self.start = QComboBox()
        self.keep_serial = QCheckBox("Edit this entry: keep its serial (add driving systems or change it)")
        v.addWidget(self.start)
        v.addWidget(self.keep_serial)
        lv.addWidget(box)

        box, v = _step(2, "Transducer")
        f = _form_layout()
        self.kind = QComboBox()
        self.kind.addItems(["annular", "matrix"])
        self.serial = QLineEdit()
        self.model_serial = QLineEdit()
        self.model_serial.setPlaceholderText("optional: serial of its device model entry")
        self.name = QLineEdit()
        self.manu = QLineEdit()
        self.freq = _dspin(250, 1, 10000, 1, " kHz")
        self.rc = _dspin(63, 1, 1000, 2, " mm")
        self.dep = _dspin(52, 0, 1000, 2, " mm")
        self.dep.setToolTip("dist_geom_ep_mm: distance from the exit plane to the centre of curvature")
        self.fmin = _dspin(10, 0, 1000)
        self.fmax = _dspin(60, 0, 1000)
        f.addRow("Type", self.kind)
        f.addRow("Serial", self.serial)
        f.addRow("Model", self.model_serial)
        f.addRow("Name", self.name)
        f.addRow("Manufacturer", self.manu)
        f.addRow("Frequency", self.freq)
        f.addRow("Curvature radius", self.rc)
        f.addRow("Exit plane → geom. focus", self.dep)
        f.addRow("Focal range (from EP)", _row(self.fmin, "to", self.fmax, "mm"))
        v.addLayout(f)

        # type-specific geometry
        self.geom = QStackedWidget()
        ann = QWidget()
        af = _form_layout()
        ann.setLayout(af)
        self.nch = _ispin(4, 1, 1024)
        self.idd = QLineEdit()
        self.idd.setPlaceholderText("e.g. 0, 30.2, 42.1, 51.1")
        self.od = QLineEdit()
        self.od.setPlaceholderText("e.g. 29.6, 41.6, 50.6, 57.9")
        af.addRow("Physical channels", self.nch)
        af.addRow("Ring inner Ø (mm)", self.idd)
        af.addRow("Ring outer Ø (mm)", self.od)
        hint = QLabel("Ring Ø: one value per ring, comma- or space-separated;\n"
                      "rings group into channels in order.")
        hint.setObjectName("CardSub")
        af.addRow(hint)
        self.geom.addWidget(ann)
        mat = QWidget()
        mf = _form_layout()
        mat.setLayout(mf)
        self.m_steer = QComboBox()
        self.m_steer.addItems(["3D", "1D"])
        self.m_curved = QCheckBox("curved (spherical)")
        self.m_curved.setChecked(True)
        self.m_D = _dspin(60, 1, 1000, 2, " mm")
        self.m_shape = QComboBox()
        self.m_shape.addItems(["rect", "disc", "bowl"])
        self.m_w = _dspin(3, 0.01, 100)
        self.m_h = _dspin(3, 0.01, 100)
        self.m_grid = QComboBox()
        self.m_grid.addItems(["rect", "fibonacci", "fermat"])
        mf.addRow("Steering", _row(self.m_steer, self.m_curved))
        mf.addRow("Aperture Ø", self.m_D)
        mf.addRow("Element shape", self.m_shape)
        mf.addRow("Element width × height (mm)", _row(self.m_w, "×", self.m_h))
        mf.addRow("Grid layout", self.m_grid)
        self.grid_stack = QStackedWidget()
        g_rect = QWidget()
        gr = _form_layout()
        g_rect.setLayout(gr)
        self.r_nrow = _ispin(16, 1, 256)
        self.r_ncol = _ispin(16, 1, 256)
        self.r_sh = _dspin(0.5, 0.001, 100, 3)
        self.r_sw = _dspin(0.5, 0.001, 100, 3)
        self.r_sp = _dspin(1, 0.01, 1)
        gr.addRow("Rows × columns", _row(self.r_nrow, "×", self.r_ncol))
        gr.addRow("Spacing h, w (mm); sparsity", _row(self.r_sh, self.r_sw, self.r_sp))
        g_fib = QWidget()
        gf = _form_layout()
        g_fib.setLayout(gf)
        self.f_n = _ispin(256, 1, 4096)
        self.f_kerf = _dspin(0.5, 0, 100)
        gf.addRow("Elements; kerf (mm)", _row(self.f_n, self.f_kerf))
        g_fm = QWidget()
        gm = _form_layout()
        g_fm.setLayout(gm)
        self.fm_n = _ispin(256, 1, 4096)
        gm.addRow("Elements", self.fm_n)
        for w in (g_rect, g_fib, g_fm):
            self.grid_stack.addWidget(w)
        mf.addRow(self.grid_stack)
        self.geom.addWidget(mat)
        v.addWidget(self.geom)
        lv.addWidget(box)

        box, v = _step(3, "Driving systems", "combinations with this transducer")
        self.rows_box = QVBoxLayout()
        self.rows_box.setSpacing(2)
        v.addLayout(self.rows_box)
        add = QPushButton("Add driving system")
        add.clicked.connect(lambda: self.add_combo())
        v.addWidget(add, 0, Qt.AlignmentFlag.AlignLeft)
        lv.addWidget(box)

        box, v = _step(4, "Save")
        row = QHBoxLayout()
        self.folder = QComboBox()
        self.folder.setEditable(True)
        self.folder.lineEdit().setPlaceholderText("your own library folder")
        self.fb = QPushButton("Browse…")
        self.fb.clicked.connect(self._browse_folder)
        row.addWidget(self.folder, 1)
        row.addWidget(self.fb)
        v.addLayout(row)
        self.shipped = QCheckBox("Shipped library (maintainers): save into config/equipment")
        self.shipped.setToolTip("As calibration_gui's AllowShipped: the entry's existing file, else "
                                "model-<MODEL>/tran-<SERIAL>_transducer.yaml, else <SERIAL>.yaml; new driving "
                                "systems in ds/.")
        v.addWidget(self.shipped)
        self.save_btn = QPushButton("Save to library")
        self.save_btn.setDefault(True)
        self.save_btn.clicked.connect(self.save)
        v.addWidget(self.save_btn, 0, Qt.AlignmentFlag.AlignLeft)
        lv.addWidget(box)
        lv.addStretch(1)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setWidget(left)
        split.addWidget(scroll)

        right = QWidget()
        rv = QVBoxLayout(right)
        rv.setContentsMargins(6, 0, 0, 0)
        views = QSplitter(Qt.Orientation.Horizontal)
        self.schematic = EditorSchematic(self.model)
        self.face = FaceView()
        views.addWidget(self.schematic)
        views.addWidget(self.face)
        views.setSizes([420, 260])
        rv.addWidget(views, 3)
        self.msg = QPlainTextEdit()
        self.msg.setReadOnly(True)
        rv.addWidget(self.msg, 2)
        split.addWidget(right)
        split.setSizes([520, 680])

        for w in (self.serial, self.model_serial, self.name, self.manu, self.idd, self.od):
            w.textChanged.connect(self.update_view)
        for w in (self.freq, self.rc, self.dep, self.fmin, self.fmax, self.m_D, self.m_w, self.m_h,
                  self.r_sh, self.r_sw, self.r_sp, self.f_kerf):
            w.valueChanged.connect(self.update_view)
        for w in (self.nch, self.r_nrow, self.r_ncol, self.f_n, self.fm_n):
            w.valueChanged.connect(self.update_view)
        for w in (self.kind, self.m_steer, self.m_shape, self.m_grid):
            w.currentTextChanged.connect(self.update_view)
        self.m_curved.toggled.connect(self.update_view)
        self.start.currentTextChanged.connect(self._start_from)
        self.keep_serial.toggled.connect(lambda _: self._start_from(self.start.currentText()))
        self.folder.currentTextChanged.connect(self.update_view)
        self.shipped.toggled.connect(self.update_view)
        self.refresh()

    # ── helpers ─────────────────────────────────────────────────────────
    def library(self):
        return self.model.library

    def profile_root(self) -> Path:
        return Path(self.window.prestus_root) / "local" / "profiles"

    def shipped_dir(self) -> Path:
        return Path(self.window.prestus_root) / "config" / "equipment"

    def n_elem(self) -> int:
        if self.kind.currentText() == "matrix":
            return int(self.mlayout["n"]) if self.mlayout and self.mlayout_key == self._matrix_key() else 1
        return self.nch.value()

    def _say(self, text: str) -> None:
        self.msg.appendPlainText(f"[{_dt.datetime.now():%H:%M:%S}] {text}")

    def refresh(self) -> None:
        """Called when the sub-tab is shown: library transducers to start from, own library folders."""
        lib = self.library()
        if lib is None:
            return
        self.model._sync_library()
        keep = self.start.currentText()
        serials = [s for s, d in lib.devices.items() if d.type in ("annular", "matrix") and not d.excluded]
        self.start.blockSignals(True)
        self.start.clear()
        self.start.addItem(BLANK)
        for s in serials:
            d = lib.devices[s]
            self.start.addItem(f"{s} ({d.type})" if d.type == "matrix" else s, s)
        idx = self.start.findData(keep) if keep != BLANK else 0
        self.start.setCurrentIndex(max(idx, 0))
        self.start.blockSignals(False)
        cur = self.folder.currentText()
        self.folder.blockSignals(True)
        self.folder.clear()
        self.folder.addItems([str(p) for p in self.model.user_library_folders()])
        if cur:
            self.folder.setCurrentText(cur)
        self.folder.blockSignals(False)
        self.update_view()

    def _browse_folder(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Your own library folder", self.folder.currentText())
        if path:
            self.folder.setCurrentText(path)

    # ── combinations ────────────────────────────────────────────────────
    def add_combo(self, combo: ed.Combo | None = None) -> None:
        lib = self.library()
        row = ComboRow(self, sorted(lib.driving_systems) if lib else [], combo)
        row.set_matrix(self.kind.currentText() == "matrix")
        self.rows.append(row)
        self.rows_box.addWidget(row)
        self.update_view()

    def remove_combo(self, row: ComboRow) -> None:
        self.rows.remove(row)
        row.setParent(None)
        row.deleteLater()
        self.update_view()

    def _clear_combos(self) -> None:
        for row in list(self.rows):
            row.setParent(None)
            row.deleteLater()
        self.rows = []

    # ── form ────────────────────────────────────────────────────────────
    def _start_from(self, _text: str = "") -> None:
        lib = self.library()
        serial = self.start.currentData()
        dev = lib.get(serial) if lib and serial else None
        keep = self.keep_serial.isChecked() and dev is not None
        raw = None
        if keep and dev.file:
            try:
                raw = yamlio.load_file(dev.file)
            except Exception:
                raw = None
        form, base = ed.form_from_device(dev, raw, keep) if dev else (ed.TransducerForm(), None)
        self.base, self.base_dev = base, (dev if keep else None)
        self._fill(form)
        if dev is not None:
            self._say(f"Editing {serial}: keep its serial to add driving systems or change it." if keep else
                      f"Showing {serial}; copied its geometry: enter a new serial to add a transducer.")

    def _fill(self, form: ed.TransducerForm) -> None:
        self._filling = True
        self.kind.setCurrentText(form.kind)
        self.serial.setText(form.serial)
        self.model_serial.setText(form.model)
        self.name.setText(form.name)
        self.manu.setText(form.manufact)
        self.freq.setValue(form.freq_khz)
        self.rc.setValue(form.curv_radius_mm)
        self.dep.setValue(form.dist_geom_ep_mm)
        self.fmin.setValue(form.min_foc)
        self.fmax.setValue(form.max_foc)
        self.nch.setValue(form.n_channels)
        self.idd.setText(form.elem_id_mm)
        self.od.setText(form.elem_od_mm)
        m = form.matrix
        self.m_steer.setCurrentText(m.steering)
        self.m_curved.setChecked(m.is_curved)
        self.m_D.setValue(m.outer_diameter_mm)
        self.m_shape.setCurrentText(m.elem_shape)
        self.m_w.setValue(m.elem_width_mm)
        self.m_h.setValue(m.elem_height_mm)
        self.m_grid.setCurrentText(m.grid)
        self.r_nrow.setValue(m.n_row)
        self.r_ncol.setValue(m.n_col)
        self.r_sh.setValue(m.spacing_height_mm)
        self.r_sw.setValue(m.spacing_width_mm)
        self.r_sp.setValue(m.sparsity)
        self.f_n.setValue(m.fib_n)
        self.f_kerf.setValue(m.fib_kerf_mm)
        self.fm_n.setValue(m.fermat_n)
        self._clear_combos()
        for c in form.combos:
            self.add_combo(c)
        self._filling = False
        self.schematic.reset_frame()
        self.update_view()

    def form(self) -> ed.TransducerForm:
        m = ed.MatrixFields(
            steering=self.m_steer.currentText(), is_curved=self.m_curved.isChecked(),
            outer_diameter_mm=self.m_D.value(), elem_shape=self.m_shape.currentText(),
            elem_width_mm=self.m_w.value(), elem_height_mm=self.m_h.value(), grid=self.m_grid.currentText(),
            n_row=self.r_nrow.value(), n_col=self.r_ncol.value(), spacing_height_mm=self.r_sh.value(),
            spacing_width_mm=self.r_sw.value(), sparsity=self.r_sp.value(), fib_n=self.f_n.value(),
            fib_kerf_mm=self.f_kerf.value(), fermat_n=self.fm_n.value())
        return ed.TransducerForm(
            kind=self.kind.currentText(), serial=self.serial.text(), name=self.name.text(),
            manufact=self.manu.text(), model=self.model_serial.text(), freq_khz=self.freq.value(),
            curv_radius_mm=self.rc.value(), dist_geom_ep_mm=self.dep.value(), min_foc=self.fmin.value(),
            max_foc=self.fmax.value(), n_channels=self.nch.value(), elem_id_mm=self.idd.text(),
            elem_od_mm=self.od.text(), matrix=m, combos=[r.combo() for r in self.rows])

    def _matrix_key(self) -> str:
        return json.dumps(ed.matrix_block(self.form()), sort_keys=True)

    def _layout_n(self) -> int | None:
        if self.mlayout is not None and self.mlayout_key == self._matrix_key() and not self.mlayout.get("error"):
            return int(self.mlayout.get("n", 0))
        return None

    def update_view(self, *_) -> None:
        if self._filling:
            return
        is_m = self.kind.currentText() == "matrix"
        self.geom.setCurrentIndex(1 if is_m else 0)
        self.grid_stack.setCurrentIndex(["rect", "fibonacci", "fermat"].index(self.m_grid.currentText()))
        for stack in (self.geom, self.grid_stack):      # size to the shown page, not the largest
            for i in range(stack.count()):
                pol = QSizePolicy.Policy.Preferred if i == stack.currentIndex() else QSizePolicy.Policy.Ignored
                stack.widget(i).setSizePolicy(pol, pol)
            stack.adjustSize()
        for r in self.rows:
            r.set_matrix(is_m)
        form = self.form()
        n_layout = self._layout_n() if is_m else None
        entry, issues = ed.build_entry(form, self.profile_root(), self.base, n_layout)
        lines: list[str] = []
        if is_m:
            m = entry["transducer"]["matrix"]
            a = m["outer_diameter_mm"] / 2
            self.schematic.set_geometry({"a": a, "R": m["curv_radius_mm"] if m["is_curved"] else None,
                                         "dist_geom_ep": m["dist_geom_ep_mm"], "rings": [(0.0, a)], "fd": None,
                                         "min_foc": form.min_foc, "max_foc": form.max_foc,
                                         "name": entry.get("name", "")})
            if self.mlayout_key != self._matrix_key():
                self._layout_timer.start()
                self.face.set_matrix([], m["elem_width_mm"], m["elem_height_mm"], m["elem_shape"], a,
                                     "Element layout: asking MATLAB…")
                lines.append("Element layout: asking MATLAB (convert_to_element_pos)…")
            elif self.mlayout.get("error"):
                self.face.set_matrix([], m["elem_width_mm"], m["elem_height_mm"], m["elem_shape"], a,
                                     "No layout")
                lines.append(f"Element layout failed: {self.mlayout['error']}")
            else:
                xy = list(zip(self.mlayout.get("x_mm") or [], self.mlayout.get("y_mm") or []))
                self.face.set_matrix(xy, m["elem_width_mm"], m["elem_height_mm"], m["elem_shape"], a)
                lines += self._matrix_info(xy, m)
            issues = [i for i in issues if "not computed yet" not in i]
        else:
            a = entry["transducer"]["annular"]
            rings = [(i / 2, o / 2) for i, o in zip(a["elem_id_mm"], a["elem_od_mm"])]
            geo = None
            if rings:
                geo = {"a": max(o for _, o in rings), "R": a["curv_radius_mm"], "dist_geom_ep": a["dist_geom_ep_mm"],
                       "rings": rings, "fd": None, "min_foc": form.min_foc, "max_foc": form.max_foc,
                       "name": entry.get("name", "")}
            self.schematic.set_geometry(geo)
            n, nch = len(rings), form.n_channels
            per = n // nch if nch and n % nch == 0 else 1
            self.face.set_annular([(ri, ro, k // per if per else k) for k, (ri, ro) in enumerate(rings)],
                                  geo["a"] if geo else 0.0)
            if geo:
                lines.append(f"{n} rings in {nch} channels, aperture Ø {2 * geo['a']:.1f} mm (colour = channel).")
                bad = [k + 1 for k in range(n - 1) if rings[k][1] > rings[k + 1][0]]
                if bad:
                    lines.append(f"Rings overlap after ring(s) {bad}.")
            else:
                lines.append("Enter the ring diameters to see the transducer.")
        lines += issues
        lines.append(self._target_note(entry))
        self.msg.setPlainText("\n".join(lines))
        self.save_btn.setEnabled(not issues and (not is_m or n_layout is not None))
        self.folder.setEnabled(not self.shipped.isChecked())
        self.fb.setEnabled(not self.shipped.isChecked())

    @staticmethod
    def _matrix_info(xy, m) -> list[str]:
        out = []
        n = len(xy)
        pitch = float("nan")
        if n > 1:
            nn = []
            for i, (x, y) in enumerate(xy):
                nn.append(min(math.hypot(x - u, y - v) for j, (u, v) in enumerate(xy) if j != i))
            nn.sort()
            pitch = nn[n // 2] if n % 2 else 0.5 * (nn[n // 2 - 1] + nn[n // 2])
        gs = m["matrix_shape"]["define_here"]["grid_shape"]
        out.append(f"{n} elements ({gs['type']} grid, {m['elem_shape']} elements), median pitch {pitch:.2f} mm.")
        w, h = m["elem_width_mm"], m["elem_height_mm"]
        esize = max(w, h) if m["elem_shape"] == "rect" else math.sqrt(4 * w * h / math.pi)
        if n > 1 and esize > pitch:
            out.append(f"Elements ({esize:.2f} mm) are larger than the pitch: they overlap.")
        if gs["type"] == "rect" and gs["rect"]["sparsity_factor"] < 1:
            out.append("sparsity_factor is validated but not applied by PRESTUS; all elements are shown.")
        return out

    def _target_note(self, entry) -> str:
        name = f"tran-{entry['serial'] or '<serial>'}_transducer.yaml"
        if self.shipped.isChecked():
            return f"Saves into the shipped {self.shipped_dir()} (maintainers)."
        folder = self.folder.currentText().strip()
        return f"Saves {name} to {folder or 'your own library folder (step 4)'}."

    # ── matrix layout (MATLAB) ──────────────────────────────────────────
    def _request_layout(self) -> None:
        if self.kind.currentText() != "matrix":
            return
        key = self._matrix_key()
        if key == self._layout_pending:
            return
        self._layout_pending = key

        def done(out):
            self._layout_pending = None
            self.mlayout_key, self.mlayout = key, out if isinstance(out, dict) else {"error": "no answer"}
            self.update_view()

        def failed(err):
            self._layout_pending = None
            self.mlayout_key, self.mlayout = key, {"error": err.message}
            self.update_view()

        self.window.bridge().call("nexus_equipment", {"action": "layout", "matrix": json.loads(key)},
                                  on_done=done, on_error=failed)

    # ── saving ──────────────────────────────────────────────────────────
    def save(self) -> None:
        form = self.form()
        n_layout = self._layout_n() if form.kind == "matrix" else None
        entry, issues = ed.build_entry(form, self.profile_root(), self.base, n_layout)
        if issues:
            self._say("Not saved: " + " ".join(issues))
            return
        shipped = self.shipped_dir()
        lib = self.library()
        existing = lib.get(entry["serial"]) if lib else None
        new_ds = [r.combo().new_ds for r in self.rows if r.is_new()]
        if self.shipped.isChecked():
            if QMessageBox.question(self, "Shipped library", "This writes into the SHIPPED config/equipment. "
                                    "Use your own library folder unless you maintain the library. Continue?") \
                    != QMessageBox.StandardButton.Yes:
                return
            dev = self.base_dev if self.base_dev is not None else existing
            files = ed.target_files(entry, shipped, dev, new_ds, shipped=shipped)
            folder = shipped
        else:
            folder = Path(self.folder.currentText().strip()).expanduser()
            if not self.folder.currentText().strip():
                self._say("Not saved: choose the folder of your own library first (step 4).")
                return
            try:
                inside = folder.resolve() == shipped.resolve() or shipped.resolve() in folder.resolve().parents
            except OSError:
                inside = False
            if inside:
                self._say("Not saved: the shipped config/equipment is not changed (tick Shipped library if you "
                          "maintain it). Choose your own library folder.")
                return
            if existing is not None and self.base_dev is None and QMessageBox.question(
                    self, "Serial exists",
                    f"{entry['serial']} is already in the library ({'your own' if existing.user else 'shipped'} "
                    "entry). Your entry replaces it. Continue?") != QMessageBox.StandardButton.Yes:
                return
            dev = self.base_dev if self.base_dev is not None else (existing if existing and existing.user else None)
            files = ed.target_files(entry, folder, dev, new_ds)
        present = [p for p in files if p.exists()]
        if present and QMessageBox.question(
                self, "Overwrite", "These files exist and are replaced:\n\n" + "\n".join(map(str, present))
                + "\n\nOverwrite them?") != QMessageBox.StandardButton.Yes:
            return
        self._check_then_write(files, folder)

    def _check_then_write(self, files: dict, folder: Path) -> None:
        """check_equipment_config in MATLAB on a scratch copy (warns only, as calibration_gui)."""
        tmp = Path(tempfile.mkdtemp(prefix="prestus_equipment_"))
        ed.write({tmp / (p.name if p.name.endswith(("_transducer.yaml", "_driver.yaml"))
                         else f"tran-{d.get('serial', 'x')}_transducer.yaml"): d for p, d in files.items()})
        self._say("Checking the entry in MATLAB (check_equipment_config)…")
        self.save_btn.setEnabled(False)

        def done(out):
            self.save_btn.setEnabled(True)
            issues = [str(x) for x in (out or {}).get("issues", []) or []] if isinstance(out, dict) else []
            if issues and QMessageBox.warning(
                    self, "Equipment check", "check_equipment_config:\n\n" + "\n".join(issues),
                    QMessageBox.StandardButton.Save | QMessageBox.StandardButton.Cancel,
                    QMessageBox.StandardButton.Cancel) != QMessageBox.StandardButton.Save:
                self._say("Not saved (equipment check).")
                return
            self._write(files, folder)

        def failed(err):
            self.save_btn.setEnabled(True)
            if QMessageBox.question(self, "Equipment check",
                                    f"MATLAB could not check the entry:\n{err.message}\n\nSave without the check?") \
                    == QMessageBox.StandardButton.Yes:
                self._write(files, folder)

        self.window.bridge().call("nexus_equipment", {"action": "check", "folder": str(tmp)},
                                  on_done=done, on_error=failed)

    def _write(self, files: dict, folder: Path) -> None:
        written = ed.write(files)
        lib = self.library()
        if lib is not None:
            if self.shipped.isChecked():
                lib._shipped_entries = None              # read the shipped folder again
            lib.reload()
        for p in written:
            self._say(f"Saved {p}")
        if not self.shipped.isChecked() and \
                folder.resolve() not in [p.resolve() for p in self.model.user_library_folders()]:
            self._say("This folder is not one of this study's libraries: set transducer.library_path or "
                      "path.transducer_library to it (or PRESTUS_TRANSDUCER_LIBRARY) so studies and the "
                      "Calibration sub-tab see the entry.")
        self.refresh()
