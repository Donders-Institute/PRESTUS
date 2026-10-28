"""Device picker for the transducer page: equipment library or custom geometry."""

from __future__ import annotations

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (QButtonGroup, QComboBox, QFormLayout, QGroupBox, QHBoxLayout,
                             QLabel, QPushButton, QRadioButton, QWidget)

from .. import placement as pl
from ..model import TX, ConfigModel, is_empty, iter_leaves
from . import style

GEOMETRY_GROUPS = ("annular", "matrix")


class DevicePanel(QGroupBox):
    """Choose a library device and driving system for one transducer."""

    changed = pyqtSignal()

    def __init__(self, model: ConfigModel, parent: QWidget | None = None):
        super().__init__("Device", parent)
        self.model = model
        self.ti = 0
        self._setting = False
        form = QFormLayout(self)

        self.rb_lib = QRadioButton("Equipment library")
        self.rb_custom = QRadioButton("Custom geometry")
        group = QButtonGroup(self)
        group.addButton(self.rb_lib)
        group.addButton(self.rb_custom)
        row = QHBoxLayout()
        row.addWidget(self.rb_lib)
        row.addWidget(self.rb_custom)
        row.addStretch(1)
        form.addRow(row)

        self.device = QComboBox()
        self.ds = QComboBox()
        form.addRow("Transducer", self.device)
        form.addRow("Driving system", self.ds)
        self.info = QLabel()
        self.info.setWordWrap(True)
        form.addRow("", self.info)
        self.inline = QLabel()
        self.inline.setWordWrap(True)
        self.inline.setStyleSheet(f"color: {style.active()['amber']};")
        self.clear_btn = QPushButton("Clear inline geometry")
        self.clear_btn.setToolTip("Remove geometry values set in this study so the library values apply")
        irow = QHBoxLayout()
        irow.addWidget(self.inline, 1)
        irow.addWidget(self.clear_btn)
        form.addRow(irow)

        devices = model.library.devices if model.library else {}
        for serial, dev in devices.items():
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
        self._setting = True
        uses_lib = not is_empty(serial)
        self.rb_lib.setChecked(uses_lib)
        self.rb_custom.setChecked(not uses_lib)
        self.device.setEnabled(uses_lib)
        self.ds.setEnabled(uses_lib)
        if uses_lib:
            idx = self.device.findData(str(serial))
            if idx < 0:
                self.device.addItem(f"{serial}  (not in library)", str(serial))
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
                              else f"<span style='color:{style.active()['red']}'>Serial not found in config/equipment.</span>")
            self.inline.setText("")
            self.clear_btn.setVisible(False)
            return
        freq = f"{dev.freq_hz / 1e3:g} kHz" if dev.freq_hz else "?"
        rng = (f"{dev.min_foc:g}–{dev.max_foc:g} mm" if dev.min_foc is not None and dev.max_foc is not None
               else "not specified")
        ds = m.get(f"{TX}.combo.ds_serial", ti)
        lines = [f"<b>{dev.manufact}</b> · {dev.type or '?'} · {dev.n_elem or '?'} elements · {freq}",
                 f"Focal range (exit plane): {rng}"]
        ep = m.get(f"{TX}.focal_distance_ep", ti)
        if isinstance(ep, (int, float)) and dev.min_foc is not None and dev.max_foc is not None:
            ok = dev.min_foc <= ep <= dev.max_foc
            colour = style.active()["green"] if ok else style.active()["red"]
            lines[-1] += f" — requested <span style='color:{colour}'>{ep:g} mm</span>"
        if dev.type == "annular":
            ann_set = all(not is_empty(m.get(f"{TX}.annular.{k}", ti))
                          for k in ("elem_phase_deg", "elem_p_amp"))
            if ann_set:
                lines.append("Calibration: phases and amplitude set inline (library not used)")
            else:
                cal = m.library.calibration(dev.serial, None if is_empty(ds) else str(ds),
                                            m.get(f"{TX}.library_path", ti) or None)
                mark = "✓" if cal.found and not cal.fallback else ("⚠" if cal.found else "✗")
                lines.append(f"Calibration: {mark} {cal.summary()}")
                lines.append(self._depth_status_line(m, ti, cal))
        self.info.setText("<br>".join(lines))
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
            return "&nbsp;&nbsp;at depth: set focal_distance_ep (or _bowl) to see coverage"
        status = cal.at_depth(ep)
        if status is None:
            return f"&nbsp;&nbsp;at {ep:g} mm: no calibrated depths in this entry"
        if status.in_range:
            lo, hi = status.bracket
            r2 = f", R² ≈ {status.nearest_r2:.2f} nearby" if status.nearest_r2 is not None else ""
            return (f"&nbsp;&nbsp;at {ep:g} mm: <span style='color:{style.active()['green']}'>within range</span> "
                    f"(between {lo:g}–{hi:g} mm, {status.phase_model}{r2})")
        colour = style.active()["amber"] if status.phase_model == "blend" else style.active()["red"]
        return (f"&nbsp;&nbsp;at {ep:g} mm: <span style='color:{colour}'>"
                f"{status.extrapolation_mm:.1f} mm past the calibrated range</span> "
                f"(nearest {status.nearest_mm:g} mm, {status.phase_model})")

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
