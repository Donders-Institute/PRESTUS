"""Tab 3 – Transducer."""

from __future__ import annotations
import os
from pathlib import Path

from PyQt6.QtWidgets import QWidget, QStackedWidget, QLabel
from widgets import (
    TextEdit, FloatEdit, Dropdown, Checkbox, CSVArrayEdit, XYZWidget, BrowseEdit,
    make_note,
)
from tabs.base import BaseTab
import style


def _load_transducer_serials(prestus_root: str) -> list[str]:
    items = ["(manual)"]
    eq_dir = os.path.join(prestus_root, "config", "equipment")
    if not os.path.isdir(eq_dir):
        return items
    try:
        import yaml
        for f in sorted(Path(eq_dir).glob("*.yaml")):
            with open(f) as fh:
                d = yaml.safe_load(fh)
            if isinstance(d, dict) and d.get("type") == "transducer":
                serial = d.get("serial", f.stem)
                items.append(str(serial))
    except Exception:
        pass
    return items


def _load_calibration_coverage(prestus_root: str, serial: str) -> str:
    lib_path = os.path.join(prestus_root, "config", "calibration", serial + ".yaml")
    if not os.path.isfile(lib_path):
        return "Library: not found"
    try:
        import yaml
        with open(lib_path) as f:
            lib = yaml.safe_load(f)
        if isinstance(lib, dict):
            gm = lib.get("global_model", {})
            depths = gm.get("depths_ep_mm", None)
            if depths is not None:
                return "Library: " + ", ".join(str(d) for d in depths) + " mm"
            cal = lib.get("calibration", {})
            fd = cal.get("focal_depths", {})
            if fd:
                import re
                vals = sorted(float(re.sub(r"[fp]", ".", k.lstrip("f")))
                              for k in fd.keys())
                return "Library: " + ", ".join(f"{v:.0f}" for v in vals) + " mm"
        return "Library: (no depth info)"
    except Exception:
        return "Library: (error reading)"


class TabTransducer(BaseTab):
    def __init__(self, registry: dict, prestus_root: str, parent=None):
        super().__init__(registry, parent)
        self._prestus_root = prestus_root
        r = registry

        # ── Library selection ─────────────────────────────────────────────────
        self._add_section("Transducer Library")
        serials = _load_transducer_serials(prestus_root)
        self._serial_dd = Dropdown("transducer.serial", serials, serials[0], registry=r)
        self._serial_dd.currentTextChanged.connect(self._on_serial_changed)
        self._add_row("Serial", self._serial_dd)
        self._add_row("Driving system serial",
                      TextEdit("transducer.combo.ds_serial", "", registry=r),
                      note="Optional. Specify for driving-system-specific calibration.")
        self._lib_cov_lbl = QLabel("Library: —")
        self._lib_cov_lbl.setStyleSheet(f"color:{style.TEXT_SUB}; font-size:{style.FS_SM}px;")
        self._layout.addWidget(self._lib_cov_lbl)

        # ── General ───────────────────────────────────────────────────────────
        self._add_section("General")
        self._type_dd = Dropdown("transducer.type", ["annular", "matrix"], "annular", registry=r)
        self._type_dd.currentTextChanged.connect(self._on_type_changed)
        self._add_row("* Type", self._type_dd)
        self._add_row("* Frequency",
                      FloatEdit("transducer.freq_hz", 500000, registry=r))
        self._add_row("Focal distance (exit plane)",
                      FloatEdit("transducer.focal_distance_ep", None, registry=r))
        self._add_row("Focal distance (bowl)",
                      FloatEdit("transducer.focal_distance_bowl", None, registry=r))

        # ── Position ──────────────────────────────────────────────────────────
        self._add_section("Transducer Position (T1 voxels)")
        self._add_row("Transducer position",
                      XYZWidget("transducer.trans_pos", (None, None, None), registry=r))
        self._add_row("Focus position",
                      XYZWidget("transducer.focus_pos", (None, None, None), registry=r))

        # ── Annular panel ─────────────────────────────────────────────────────
        self._stacked = QStackedWidget()
        self._layout.addWidget(self._stacked)

        ann = QWidget()
        ann.setStyleSheet("background: transparent;")
        from tabs.base import BaseTab as _B
        from PyQt6.QtWidgets import QVBoxLayout
        ann_layout = QVBoxLayout(ann)
        ann_layout.setContentsMargins(0, 0, 0, 0)
        ann_layout.setSpacing(2)

        def _sec(text):
            from widgets import make_section
            ann_layout.addSpacing(6)
            ann_layout.addWidget(make_section(text))

        def _row(label, widget, note=""):
            from PyQt6.QtWidgets import QGridLayout, QWidget as QW
            from widgets import make_label, make_note as mn
            row = QW()
            row.setStyleSheet("background: transparent;")
            gl = QGridLayout(row)
            gl.setContentsMargins(0, 1, 0, 1)
            gl.setColumnStretch(1, 1)
            gl.setColumnMinimumWidth(0, 210)
            gl.addWidget(make_label(label), 0, 0)
            gl.addWidget(widget, 0, 1)
            if note:
                gl.addWidget(mn(note), 1, 1)
            ann_layout.addWidget(row)

        _sec("Annular Array")
        _row("* Elements",
             FloatEdit("transducer.annular.elem_n", 4, registry=r))
        _row("* Inner diameters (comma-sep)",
             CSVArrayEdit("transducer.annular.elem_id_mm", "0,34,53,70", registry=r),
             note="mm")
        _row("* Outer diameters (comma-sep)",
             CSVArrayEdit("transducer.annular.elem_od_mm", "32,51,67,80", registry=r),
             note="mm")
        _row("* Radius of curvature",
             FloatEdit("transducer.annular.curv_radius_mm", 63.2, registry=r))
        _row("* Pressure amplitude",
             FloatEdit("transducer.annular.elem_amp", 1, registry=r))
        _row("Phase per element (°)",
             CSVArrayEdit("transducer.annular.elem_phase_deg", "0", registry=r))
        _row("Geometric focus-to-EP dist.",
             FloatEdit("transducer.annular.dist_geom_ep_mm", None, registry=r))
        _row("Visualization depth",
             FloatEdit("transducer.annular.depth_mm", 16, registry=r))
        ann_layout.addStretch()
        self._stacked.addWidget(ann)   # index 0 = annular

        # ── Matrix panel ──────────────────────────────────────────────────────
        mat = QWidget()
        mat.setStyleSheet("background: transparent;")
        mat_layout = QVBoxLayout(mat)
        mat_layout.setContentsMargins(0, 0, 0, 0)
        mat_layout.setSpacing(2)

        def _msec(text):
            from widgets import make_section
            mat_layout.addSpacing(6)
            mat_layout.addWidget(make_section(text))

        def _mrow(label, widget, note=""):
            from PyQt6.QtWidgets import QGridLayout, QWidget as QW
            from widgets import make_label, make_note as mn
            row = QW()
            row.setStyleSheet("background: transparent;")
            gl = QGridLayout(row)
            gl.setContentsMargins(0, 1, 0, 1)
            gl.setColumnStretch(1, 1)
            gl.setColumnMinimumWidth(0, 210)
            gl.addWidget(make_label(label), 0, 0)
            gl.addWidget(widget, 0, 1)
            if note:
                gl.addWidget(mn(note), 1, 1)
            mat_layout.addWidget(row)

        _msec("Matrix Array")
        _mrow("Element shape",
              Dropdown("transducer.matrix.elem_shape", ["rect", "disc", "bowl"], "rect", registry=r))
        _mrow("Element height",
              FloatEdit("transducer.matrix.elem_height_mm", 1, registry=r))
        _mrow("Element width",
              FloatEdit("transducer.matrix.elem_width_mm", 1, registry=r))
        _mrow("Outer diameter",
              FloatEdit("transducer.matrix.outer_diameter_mm", 70, registry=r))
        _mrow("Pressure amplitude",
              FloatEdit("transducer.matrix.elem_amp", 1, registry=r))
        _mrow("Curved surface",
              Checkbox("transducer.matrix.is_curved", "Curved aperture", False, registry=r))
        _mrow("Radius of curvature",
              FloatEdit("transducer.matrix.curv_radius_mm", None, registry=r))
        _mrow("Dist. geom. focus-to-EP",
              FloatEdit("transducer.matrix.dist_geom_ep_mm", None, registry=r))
        _mrow("Steering",
              Dropdown("transducer.matrix.steering", ["1D", "3D"], "3D", registry=r))
        _mrow("Visualization depth",
              FloatEdit("transducer.matrix.depth_mm", 16, registry=r))
        _mrow("Clover multi-aperture",
              Checkbox("transducer.matrix.is_clover_setup", "Enable Clover replication", False, registry=r))
        mat_layout.addStretch()
        self._stacked.addWidget(mat)   # index 1 = matrix

        self._finish()

    def _on_type_changed(self, val: str):
        self._stacked.setCurrentIndex(0 if val == "annular" else 1)

    def _on_serial_changed(self, serial: str):
        if serial == "(manual)":
            self._lib_cov_lbl.setText("Library: —")
            return
        self._lib_cov_lbl.setText(_load_calibration_coverage(self._prestus_root, serial))
        # Auto-populate geometry from equipment YAML
        try:
            import yaml, os
            eq_dir = os.path.join(self._prestus_root, "config", "equipment")
            for f in Path(eq_dir).glob("*.yaml"):
                with open(f) as fh:
                    d = yaml.safe_load(fh)
                if isinstance(d, dict) and str(d.get("serial", "")) == serial:
                    geom = d.get("transducer", {})
                    if "freq_hz" in geom and "transducer.freq_hz" in self.registry:
                        self.registry["transducer.freq_hz"].set_value(geom["freq_hz"])
                    if "annular" in geom:
                        ann = geom["annular"]
                        self._type_dd.set_value("annular")
                        self._on_type_changed("annular")
                        for field in ("elem_n", "curv_radius_mm", "dist_geom_ep_mm"):
                            tag = f"transducer.annular.{field}"
                            if tag in self.registry and field in ann:
                                self.registry[tag].set_value(ann[field])
                        for field in ("elem_id_mm", "elem_od_mm"):
                            tag = f"transducer.annular.{field}"
                            if tag in self.registry and field in ann:
                                self.registry[tag].set_value(ann[field])
                    elif "matrix" in geom:
                        self._type_dd.set_value("matrix")
                        self._on_type_changed("matrix")
                    break
        except Exception as exc:
            self._lib_cov_lbl.setText(f"Library: (error — {exc})")
