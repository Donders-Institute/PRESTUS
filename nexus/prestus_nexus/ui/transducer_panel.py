"""Device picker for the transducer page: equipment library or custom geometry."""

from __future__ import annotations

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (QButtonGroup, QComboBox, QGridLayout, QGroupBox, QHBoxLayout,
                             QLabel, QPushButton, QRadioButton, QWidget)

from .. import placement as pl
from ..model import TX, ConfigModel, is_empty, iter_leaves
from . import style
from .readout import muted, readout, status

GEOMETRY_GROUPS = ("annular", "matrix")
# The panel sits inside a card, indented by the group box frame; its label column is that much
# narrower than the card's (section.LABEL_WIDTH = 290), so its inputs start where the card's do.
LABEL_COLUMN = 278
GRID_SPACING = 14


class DevicePanel(QGroupBox):
    """Choose a library device and driving system for one transducer."""

    changed = pyqtSignal()

    def __init__(self, model: ConfigModel, parent: QWidget | None = None):
        super().__init__("Device", parent)
        self.model = model
        self.ti = 0
        self._setting = False
        form = QGridLayout(self)
        form.setHorizontalSpacing(GRID_SPACING)
        form.setVerticalSpacing(8)
        form.setColumnMinimumWidth(0, LABEL_COLUMN)
        form.setColumnStretch(1, 1)

        def add_row(label, item) -> None:
            r = form.rowCount()
            if label is None:                       # spans both columns
                if isinstance(item, QWidget):
                    form.addWidget(item, r, 0, 1, 2)
                else:
                    form.addLayout(item, r, 0, 1, 2)
                return
            form.addWidget(QLabel(label), r, 0)
            form.addWidget(item, r, 1)

        self.rb_lib = QRadioButton("Equipment library")
        self.rb_custom = QRadioButton("Custom geometry")
        group = QButtonGroup(self)
        group.addButton(self.rb_lib)
        group.addButton(self.rb_custom)
        row = QHBoxLayout()
        row.addWidget(self.rb_lib)
        row.addWidget(self.rb_custom)
        row.addStretch(1)
        add_row(None, row)

        self.device = QComboBox()
        self.ds = QComboBox()
        add_row("Transducer", self.device)
        add_row("Driving system", self.ds)
        self.info = QLabel()
        self.info.setWordWrap(True)
        add_row(None, self.info)
        self.inline = QLabel()
        self.inline.setWordWrap(True)
        self.inline.setStyleSheet(f"color: {style.active()['amber']};")
        self.clear_btn = QPushButton("Clear inline geometry")
        self.clear_btn.setToolTip("Remove geometry values set in this study so the library values apply")
        irow = QHBoxLayout()
        irow.addWidget(self.inline, 1)
        irow.addWidget(self.clear_btn)
        add_row(None, irow)

        devices = model.library.devices if model.library else {}
        for serial, dev in devices.items():
            if not dev.excluded:
                self.device.addItem(dev.label(), serial)
        if not devices:
            self.rb_lib.setEnabled(False)
            self.rb_lib.setToolTip("No equipment library found (config/equipment)")

        self.rb_lib.toggled.connect(self._source_changed)
        self.device.currentIndexChanged.connect(self._device_changed)
        self.ds.currentIndexChanged.connect(self._ds_changed)
        self.clear_btn.clicked.connect(self._clear_inline)

    # ── model → widgets ─────────────────────────────────────────────────
    def refresh(self, ti: int) -> None:
        self.ti = ti
        m = self.model
        serial = m.get(f"{TX}.serial", ti)
        dev = m.device(ti)
        own = m.user_library_folders() if m.library is not None else []
        if m.library is not None:          # the libraries checked, most specific first
            self.rb_lib.setToolTip("Libraries checked, in order:\n" + "\n".join(
                [f"own: {p}" for p in own] + ["shipped: config/equipment, config/calibration"]))
        self._setting = True
        uses_lib = not is_empty(serial)
        self.rb_lib.setChecked(uses_lib)
        self.rb_custom.setChecked(not uses_lib)
        self.device.setEnabled(uses_lib)
        self.ds.setEnabled(uses_lib)
        if uses_lib:
            idx = self.device.findData(str(serial))
            if idx < 0:
                known = m.library.get(serial) if m.library else None
                if known is not None:              # an excluded entry the study names
                    self.device.addItem(known.label(), str(serial))
                else:
                    fallback = m.library.model_for(serial) if m.library else None
                    note = f"uses model {fallback}" if fallback else "not in library"
                    self.device.addItem(f"{serial}  ({note})", str(serial))
                idx = self.device.count() - 1
            self.device.setCurrentIndex(idx)
        self._fill_ds(dev)
        self._setting = False
        self._update_info()

    def _fill_ds(self, dev) -> None:
        current = self.model.get(f"{TX}.combo.ds_serial", self.ti)
        self.ds.clear()
        self.ds.addItem("(none: DS-agnostic entry)", "")
        if dev is not None:
            for key, ds in dev.combos.items():
                self.ds.addItem(f"{key} — {ds.name}" if ds.name else key, key)
        if not is_empty(current):
            idx = self.ds.findData(str(current))
            if idx < 0:
                self.ds.addItem(f"{current}  (not listed)", str(current))
                idx = self.ds.count() - 1
            self.ds.setCurrentIndex(idx)

    def _update_info(self) -> None:
        m, ti = self.model, self.ti
        dev = m.device(ti)
        if dev is None:
            self.info.setText("Geometry is taken from the fields below." if self.rb_custom.isChecked()
                              else status("Serial not found in config/equipment.", "bad"))
            self.inline.setText("")
            self.clear_btn.setVisible(False)
            return
        freq = f"{dev.freq_hz / 1e3:g} kHz" if dev.freq_hz else "?"
        ds = m.get(f"{TX}.combo.ds_serial", ti)
        rows = [("Device", f"<b>{dev.manufact}</b> · {dev.type or '?'} · {dev.n_elem or '?'} elements "
                           f"· {freq}")]
        serial = m.get(f"{TX}.serial", ti)
        if not is_empty(serial) and str(serial) != dev.serial:
            rows.insert(0, ("Library entry", status(f"model {dev.serial}", "warn")
                            + muted(f" ({serial} is not in the library: model geometry, "
                                    "geometric calibration)")))
        if dev.min_foc is not None and dev.max_foc is not None:
            rng = f"{dev.min_foc:g}–{dev.max_foc:g} mm" + muted(" from the exit plane")
            ep = m.get(f"{TX}.focal_distance_ep", ti)
            if isinstance(ep, (int, float)) and not isinstance(ep, bool):
                ok = dev.min_foc <= ep <= dev.max_foc
                rng += "<br>focal depth " + status(f"{ep:g} mm · {'inside' if ok else 'outside'}",
                                                   "ok" if ok else "bad")
        else:
            rng = muted("not specified")
        rows.append(("Steering range", rng))
        if dev.type == "annular":
            ann_set = all(not is_empty(m.get(f"{TX}.annular.{k}", ti))
                          for k in ("elem_phase_deg", "elem_p_amp"))
            if ann_set:
                rows.append(("Calibration", "phases and amplitude set inline" + muted(" (library not used)")))
            else:
                cal = m.library.calibration(dev.serial, None if is_empty(ds) else str(ds),
                                            m.get(f"{TX}.library_path", ti) or None)
                level = "ok" if cal.found and not cal.fallback else ("warn" if cal.found else "bad")
                rows.append(("Calibration", status(cal.summary(), level)))
                rows.append(("At the focal depth", self._depth_status_line(m, ti, cal)))
        self.info.setText(readout([("Equipment", rows)], label_width=LABEL_COLUMN + GRID_SPACING - 4))
        overrides = self._inline_geometry()
        if overrides:
            self.inline.setText(f"{len(overrides)} geometry value(s) set in this study override the "
                                f"library: {', '.join(p.split('.', 1)[1] for p in overrides[:4])}"
                                f"{' …' if len(overrides) > 4 else ''}")
        else:
            self.inline.setText("")
        self.clear_btn.setVisible(bool(overrides))

    def _depth_status_line(self, m: ConfigModel, ti: int, cal) -> str:
        """Coverage of the calibration entry at this transducer's resolved
        focal_distance_ep -- pure arithmetic against what's already in the YAML
        (Calibration.at_depth), not a re-run of the forward model."""
        txs = m.effective_transducers()
        bowl = pl.Bowl.from_transducer(txs[ti]) if ti < len(txs) else None
        ep = pl.resolve_focal_distance_ep(txs[ti], bowl) if ti < len(txs) else None
        if ep is None:
            return "set focal_distance_ep (or _bowl) to see coverage"
        cov = cal.at_depth(ep)
        if cov is None:
            return f"{ep:g} mm: no calibrated depths in this entry"
        if cov.in_range:
            lo, hi = cov.bracket
            r2 = f", R² ≈ {cov.nearest_r2:.2f} nearby" if cov.nearest_r2 is not None else ""
            return (f"{ep:g} mm: {status('within range', 'ok')} "
                    f"(between {lo:g}–{hi:g} mm, {cov.phase_model}{r2})")
        level = "warn" if cov.phase_model == "blend" else "bad"
        return (f"{ep:g} mm: "
                + status(f"{cov.extrapolation_mm:.1f} mm past the calibrated range", level) + " "
                f"(nearest {cov.nearest_mm:g} mm, {cov.phase_model})")

    def _inline_geometry(self) -> list[str]:
        t = self.model.tx[self.ti]
        paths = [f"{TX}.{p}" for p, _ in iter_leaves(t)
                 if p.split(".")[0] in GEOMETRY_GROUPS or p == "freq_hz"]
        return [p for p in paths if self.model.is_modified(p, self.ti)]

    # ── widgets → model ─────────────────────────────────────────────────
    def _source_changed(self, library: bool) -> None:
        if self._setting:
            return
        if library:
            self._device_changed(self.device.currentIndex())
        else:
            self.model.reset(f"{TX}.serial", self.ti)
            self.model.reset(f"{TX}.combo.ds_serial", self.ti)
            self.changed.emit()

    def _device_changed(self, index: int) -> None:
        if self._setting or index < 0 or not self.rb_lib.isChecked():
            return
        self.model.set(f"{TX}.serial", self.device.itemData(index), self.ti)
        self.model.reset(f"{TX}.combo.ds_serial", self.ti)
        dev = self.model.device(self.ti)
        if dev is not None and len(dev.combos) == 1:
            self.model.set(f"{TX}.combo.ds_serial", next(iter(dev.combos)), self.ti)
        self.changed.emit()

    def _ds_changed(self, index: int) -> None:
        if self._setting or index < 0:
            return
        value = self.ds.itemData(index)
        if value:
            self.model.set(f"{TX}.combo.ds_serial", value, self.ti)
        else:
            self.model.reset(f"{TX}.combo.ds_serial", self.ti)
        self.changed.emit()

    def _clear_inline(self) -> None:
        for path in self._inline_geometry():
            self.model.reset(path, self.ti)
        self.changed.emit()
