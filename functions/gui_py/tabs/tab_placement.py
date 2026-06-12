"""Tab 9 – Placement (with interactive MRI viewer for manual mode)."""

from __future__ import annotations
import os

from PyQt6.QtWidgets import (
    QWidget, QStackedWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QButtonGroup, QRadioButton,
)
from PyQt6.QtCore import Qt

from widgets import (
    TextEdit, FloatEdit, Dropdown, Checkbox, CSVArrayEdit,
    BrowseEdit, XYZWidget, make_note,
)
from tabs.base import BaseTab
import style


class TabPlacement(BaseTab):
    def __init__(self, registry: dict, parent=None):
        super().__init__(registry, parent)
        r = registry

        self._add_section("Transducer Placement")
        self._mode_dd = Dropdown(
            "placement.mode",
            ["manual", "localite", "heuristic", "plantus"],
            "manual", registry=r
        )
        self._mode_dd.currentTextChanged.connect(self._on_mode_changed)
        self._add_row("Placement mode", self._mode_dd,
                      note="manual: use trans_pos/focus_pos as-is.  "
                           "localite: read from XML.  "
                           "heuristic: sphere-expansion search.  "
                           "plantus: multi-objective optimisation.")

        # ── Stacked sub-panels ────────────────────────────────────────────────
        self._stacked = QStackedWidget()
        self._layout.addWidget(self._stacked)

        self._stacked.addWidget(self._build_manual_panel(r))   # 0
        self._stacked.addWidget(self._build_localite_panel(r)) # 1
        self._stacked.addWidget(self._build_heuristic_panel(r))# 2
        self._stacked.addWidget(self._build_plantus_panel(r))  # 3

        self._on_mode_changed("manual")
        self._finish()

    # ── Mode sub-panels ───────────────────────────────────────────────────────

    def _build_manual_panel(self, r):
        w = QWidget()
        w.setStyleSheet("background: transparent;")
        vl = QVBoxLayout(w)
        vl.setContentsMargins(0, 8, 0, 0)
        vl.setSpacing(4)

        note = make_note(
            "Transducer and focus positions are taken directly from the Transducer tab "
            "(transducer.trans_pos / transducer.focus_pos). Use the MRI viewer below to "
            "click-to-place them."
        )
        vl.addWidget(note)

        # ── MRI Viewer ─────────────────────────────────────────────────────
        self._viewer = None
        self._place_mode_group = None
        self._overlay_note = None

        try:
            from nifti_viewer import NiftiViewer
            self._viewer = NiftiViewer()
            self._viewer.setMinimumHeight(300)
            self._viewer.cursor_moved.connect(self._on_cursor_moved)

            # Load anatomy button + mode toggle
            tb = QHBoxLayout()
            btn_load = QPushButton("Load T1w NIfTI…")
            btn_load.setStyleSheet(
                f"background:{style.ACCENT}; color:white; border-radius:4px; padding:5px 10px;"
            )
            btn_load.clicked.connect(self._load_anatomy)
            tb.addWidget(btn_load)

            self._place_mode_group = QButtonGroup()
            for label, tag in (("Set transducer pos", "trans"), ("Set focus pos", "focus")):
                rb = QRadioButton(label)
                rb.setStyleSheet(f"font-size:{style.FS_SM}px;")
                if tag == "trans":
                    rb.setChecked(True)
                self._place_mode_group.addButton(rb, 0 if tag == "trans" else 1)
                tb.addWidget(rb)
            tb.addStretch()

            self._overlay_note = QLabel("Click any slice to place the selected position.")
            self._overlay_note.setStyleSheet(
                f"color:{style.TEXT_SUB}; font-size:{style.FS_SM}px;"
            )
            tb.addWidget(self._overlay_note)
            vl.addLayout(tb)
            vl.addWidget(self._viewer)
        except Exception:
            vl.addWidget(make_note(
                "NIfTI viewer unavailable — install matplotlib and nibabel to enable."
            ))
            self._viewer = None
        return w

    def _build_localite_panel(self, r):
        w = QWidget()
        w.setStyleSheet("background: transparent;")
        vl = QVBoxLayout(w)
        vl.setContentsMargins(0, 8, 0, 0)
        vl.setSpacing(2)

        def _row(label, widget, note=""):
            from PyQt6.QtWidgets import QGridLayout
            from widgets import make_label, make_note as mn
            row = QWidget()
            row.setStyleSheet("background: transparent;")
            gl = QGridLayout(row)
            gl.setContentsMargins(0, 1, 0, 1)
            gl.setColumnStretch(1, 1)
            gl.setColumnMinimumWidth(0, 220)
            gl.addWidget(make_label(label), 0, 0)
            gl.addWidget(widget, 0, 1)
            if note:
                gl.addWidget(mn(note), 1, 1)
            vl.addWidget(row)

        from widgets import make_section
        vl.addWidget(make_section("Localite Config"))
        _row("Localite XML file",
             BrowseEdit("placement.localite.file", "file", registry=r))
        _row("Session index",
             FloatEdit("placement.localite.session", 1, registry=r))
        _row("Marker type",
             Dropdown("placement.localite.markertype",
                       ["TriggerMarkers", "NavigationMarkers"],
                       "TriggerMarkers", registry=r))
        _row("Marker position index",
             FloatEdit("placement.localite.position", 1, registry=r))
        _row("Reference distance",
             FloatEdit("placement.localite.reference_distance_mm", 15, registry=r),
             note="mm — distance between IR trackers and transducer exit plane.")
        _row("Save Localite-aligned T1",
             Checkbox("placement.localite.save_localite_t1",
                      "Save T1 aligned to Localite header for QC", False, registry=r))
        vl.addStretch()
        return w

    def _build_heuristic_panel(self, r):
        w = QWidget()
        w.setStyleSheet("background: transparent;")
        vl = QVBoxLayout(w)
        vl.setContentsMargins(0, 8, 0, 0)
        vl.setSpacing(2)

        def _row(label, widget, note=""):
            from PyQt6.QtWidgets import QGridLayout
            from widgets import make_label, make_note as mn
            row = QWidget()
            row.setStyleSheet("background: transparent;")
            gl = QGridLayout(row)
            gl.setContentsMargins(0, 1, 0, 1)
            gl.setColumnStretch(1, 1)
            gl.setColumnMinimumWidth(0, 220)
            gl.addWidget(make_label(label), 0, 0)
            gl.addWidget(widget, 0, 1)
            if note:
                gl.addWidget(mn(note), 1, 1)
            vl.addWidget(row)

        from widgets import make_section
        vl.addWidget(make_section("Heuristic (Sphere-Expansion) Config"))
        _row("MNI target (mm)",
             XYZWidget("placement.heuristic.mni_target_mm", (None, None, None), registry=r))
        _row("Target name",
             TextEdit("placement.heuristic.target_name", "", registry=r))
        _row("Dist. close",
             FloatEdit("placement.heuristic.dist_close", None, registry=r), note="mm")
        _row("Ear exclusion radius",
             FloatEdit("placement.heuristic.ear_radius", 35, registry=r), note="mm")
        _row("Save Localite-aligned T1",
             Checkbox("placement.heuristic.save_localite_t1",
                      "Save T1 aligned to heuristic result for QC", False, registry=r))
        vl.addStretch()
        return w

    def _build_plantus_panel(self, r):
        from PyQt6.QtWidgets import QScrollArea
        outer = QScrollArea()
        outer.setWidgetResizable(True)
        outer.setStyleSheet("border: none;")
        w = QWidget()
        w.setStyleSheet("background: transparent;")
        vl = QVBoxLayout(w)
        vl.setContentsMargins(0, 8, 0, 8)
        vl.setSpacing(2)
        outer.setWidget(w)

        def _row(label, widget, note=""):
            from PyQt6.QtWidgets import QGridLayout
            from widgets import make_label, make_note as mn
            row = QWidget()
            row.setStyleSheet("background: transparent;")
            gl = QGridLayout(row)
            gl.setContentsMargins(0, 1, 0, 1)
            gl.setColumnStretch(1, 1)
            gl.setColumnMinimumWidth(0, 220)
            gl.addWidget(make_label(label), 0, 0)
            gl.addWidget(widget, 0, 1)
            if note:
                gl.addWidget(mn(note), 1, 1)
            vl.addWidget(row)

        from widgets import make_section, make_note as mn
        vl.addWidget(make_section("PlanTUS Config"))
        _row("* PlanTUS script path",
             BrowseEdit("placement.plantus.script_path", "dir", registry=r, mandatory=True))
        _row("* SimNIBS env path",
             BrowseEdit("placement.plantus.env_path", "dir", registry=r, mandatory=True))
        _row("* MNI target (mm)",
             XYZWidget("placement.plantus.mni_target_mm", (None, None, None), registry=r))
        _row("* Target name",
             TextEdit("placement.plantus.target_name", "", registry=r))
        _row("* Focal distance list",
             CSVArrayEdit("placement.plantus.focal_distance_list", "", registry=r),
             note="mm — comma-separated calibration focal distances from exit plane.")
        _row("* FLHM list",
             CSVArrayEdit("placement.plantus.flhm_list", "", registry=r),
             note="mm — comma-separated full-width half-max focal lengths.")
        _row("Max steering angle",
             FloatEdit("placement.plantus.max_angle_deg", 10, registry=r), note="deg")
        _row("Additional offset",
             FloatEdit("placement.plantus.additional_offset_mm", 0, registry=r), note="mm")
        _row("Mask radius",
             FloatEdit("placement.plantus.mask_radius_mm", 2, registry=r), note="mm")
        _row("Connectome WB path",
             BrowseEdit("placement.plantus.connectome_wb_path", "dir", registry=r),
             note="Optional — required for interactive GUI visualisation.")

        vl.addWidget(make_section("Optimisation Weights (must sum to 1.0)"))
        weights = [
            ("placement.plantus.weights.skin_target_distances",     "Skin-target distance",      0.2),
            ("placement.plantus.weights.skin_target_angles",        "Skin-target angle",         0.2),
            ("placement.plantus.weights.skin_target_intersections", "Skin-target intersection",  0.2),
            ("placement.plantus.weights.skin_skull_angles",         "Skin-skull angle",          0.2),
            ("placement.plantus.weights.skull_thickness",           "Skull thickness",           0.2),
        ]
        for tag, label, default in weights:
            _row(label, FloatEdit(tag, default, registry=r))

        vl.addStretch()
        return outer

    # ── Callbacks ─────────────────────────────────────────────────────────────

    def _on_mode_changed(self, mode: str):
        idx = {"manual": 0, "localite": 1, "heuristic": 2, "plantus": 3}.get(mode, 0)
        self._stacked.setCurrentIndex(idx)

    def _load_anatomy(self):
        from PyQt6.QtWidgets import QFileDialog
        path, _ = QFileDialog.getOpenFileName(
            self, "Load T1w NIfTI",
            filter="NIfTI (*.nii *.nii.gz);;All files (*)"
        )
        if path and self._viewer:
            self._viewer.load_anatomy(path)

    def _on_cursor_moved(self, x: int, y: int, z: int):
        """Update trans_pos or focus_pos fields depending on radio button state."""
        if not self._place_mode_group:
            return
        mode_id = self._place_mode_group.checkedId()
        if mode_id == 0:
            base = "transducer.trans_pos"
        else:
            base = "transducer.focus_pos"
        for i, val in enumerate((x, y, z)):
            tag = f"{base}_{i}"
            if tag in self.registry:
                self.registry[tag].set_value(val)
