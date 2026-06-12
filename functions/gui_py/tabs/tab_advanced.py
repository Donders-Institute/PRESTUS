"""Tab 8 – Advanced."""

from __future__ import annotations
from PyQt6.QtWidgets import QPushButton, QLabel, QHBoxLayout, QWidget, QFileDialog
from widgets import TextEdit, FloatEdit, Dropdown, Checkbox, CSVArrayEdit, BrowseEdit
from tabs.base import BaseTab
import style


class TabAdvanced(BaseTab):
    def __init__(self, registry: dict, parent=None):
        super().__init__(registry, parent)
        r = registry

        # ── Segmentation ──────────────────────────────────────────────────────
        self._add_section("Segmentation")
        self._add_row("Force qform reorientation",
                      Checkbox("segmentation.use_qform",
                               "Fix qform/sform mismatch before charm", False, registry=r))
        self._add_row("Debug segmentation",
                      Checkbox("segmentation.debug",
                               "Pass --debug flag to charm", False, registry=r))

        # ── Head Model ────────────────────────────────────────────────────────
        self._add_section("Head Model Processing")
        self._add_row("Head pad",
                      FloatEdit("headmodel.head_pad_mm", 0, registry=r), note="mm")
        self._add_row("CSF expansion",
                      FloatEdit("headmodel.csf_expansion", 40, registry=r), note="voxels")
        self._add_row("Smoothing method",
                      Dropdown("headmodel.smooth_method", ["gaussian", "box", "off"],
                                "gaussian", registry=r))
        self._add_row("Smoothing FWHM",
                      FloatEdit("headmodel.smooth_fwhm_mm", 1, registry=r), note="mm")
        self._add_row("Skull threshold",
                      FloatEdit("headmodel.smooth_threshold_skull", 0.5, registry=r),
                      note="Higher = thinner skull mask after smoothing.")
        self._add_row("Other tissue threshold",
                      FloatEdit("headmodel.smooth_threshold_other", 0.5, registry=r))
        self._add_row("Smooth acoustic properties",
                      Checkbox("headmodel.smooth_properties",
                               "Apply smoothing kernel to acoustic property maps", False, registry=r))
        self._add_row("Skull fill method",
                      Dropdown("headmodel.skull_fill_method",
                                ["rubberwrap", "imclose"], "rubberwrap", registry=r))
        self._add_row("Rubber-wrap radius",
                      FloatEdit("headmodel.skull_wrap_radius", 10, registry=r), note="voxels")
        self._add_row("Visualize rubber-wrap",
                      Checkbox("headmodel.skull_wrap_visualize",
                               "Show rubber-wrap result (disable on HPC)", False, registry=r))

        # ── Pseudo-CT ─────────────────────────────────────────────────────────
        self._add_section("Pseudo-CT Skull Mapping")
        self._pct_enabled = Checkbox("pct.enabled",
                                     "Use CT/pseudo-CT to inform skull properties",
                                     False, registry=r)
        self._pct_enabled.stateChanged.connect(self._toggle_pct)
        self._add_row("Enable pCT", self._pct_enabled)

        self._skull_map = Dropdown("pct.skull_mapping",
                                    ["kosciessa", "miscouridou", "carpino",
                                     "wiesinger", "treeby"],
                                    "kosciessa", registry=r)
        self._pct_debug = Checkbox("pct.debug",
                                   "Keep intermediate NIfTI files in pseudoCT/ subfolder",
                                   True, registry=r)
        self._density = Dropdown("pct.mapping_density",
                                  ["k-plan", "k-wave", "marsac", "aubry", "none"],
                                  "k-plan", registry=r)
        self._soundspeed = Dropdown("pct.mapping_soundspeed",
                                     ["k-plan", "marsac", "aubry", "none"],
                                     "k-plan", registry=r)
        self._attenuation = Dropdown("pct.mapping_attenuation",
                                      ["k-plan", "mueller", "aubry", "none"],
                                      "k-plan", registry=r)
        self._add_row("UTE→HU skull mapping", self._skull_map)
        self._add_row("Debug (keep intermediates)", self._pct_debug)
        self._add_row("Density mapping", self._density)
        self._add_row("Sound speed mapping", self._soundspeed)
        self._add_row("Attenuation mapping", self._attenuation)
        self._toggle_pct(0)  # start disabled

        # ── Analysis ──────────────────────────────────────────────────────────
        self._add_section("Analysis")
        self._add_row("Focus area radius",
                      FloatEdit("analysis.focus_area_radius", 5, registry=r),
                      note="mm — radius around focus for ISPPA averaging in output metrics.")

        # ── Neuronavigation Ingestion ─────────────────────────────────────────
        self._add_section("Neuronavigation Ingestion")
        self._add_row("Subject ID",
                      TextEdit("neuronav.sub_id", "", registry=r))
        self._add_row("Session ID",
                      TextEdit("neuronav.ses_id", "", registry=r))
        self._add_note("Marker type (TriggerMarkers / GUMMarkers) is auto-detected from the folder.")
        self._add_row("Raw Localite path",
                      BrowseEdit("neuronav.raw_path", "dir", registry=r))
        self._add_row("Output path",
                      BrowseEdit("neuronav.out_path", "dir", registry=r))

        self._ingest_status = QLabel("")
        self._ingest_status.setStyleSheet(
            f"color:{style.TEXT_SUB}; font-size:{style.FS_SM}px;"
        )
        btn_ingest = QPushButton("▶  Run Ingestion")
        btn_ingest.setProperty("class", "success")
        btn_ingest.setStyleSheet(
            f"background:{style.ACCENT}; color:white; border-radius:4px; padding:6px 12px;"
        )
        btn_ingest.clicked.connect(self._run_ingest)
        self._layout.addWidget(btn_ingest)
        self._layout.addWidget(self._ingest_status)

        # ── Simulation Layers ─────────────────────────────────────────────────
        self._add_section("Simulation Layers")
        self._add_note("Comma-separated SimNIBS label indices assigned to each tissue compartment.")
        layer_defaults = {
            "water": "0,3,6,9,10", "brain": "1,2", "skin": "5",
            "skull": "4", "skull_cortical": "7", "skull_trabecular": "8",
        }
        for tissue, default in layer_defaults.items():
            self._add_row(tissue,
                          CSVArrayEdit(f"layers.{tissue}", default, registry=r))

        self._finish()

    def _toggle_pct(self, state):
        enabled = bool(state)
        for w in (self._skull_map, self._pct_debug, self._density,
                  self._soundspeed, self._attenuation):
            w.setEnabled(enabled)

    def _run_ingest(self):
        r = self.registry
        sub_id = r.get("neuronav.sub_id", None)
        ses_id = r.get("neuronav.ses_id", None)
        raw_path = r.get("neuronav.raw_path", None)
        out_path = r.get("neuronav.out_path", None)

        def _val(w):
            return w.value() if w and hasattr(w, "value") else ""

        sub_id   = _val(sub_id)
        ses_id   = _val(ses_id)
        raw_path = _val(raw_path)
        out_path = _val(out_path)

        if not all([sub_id, ses_id, raw_path, out_path]):
            self._ingest_status.setStyleSheet(f"color:{style.MANDATORY};")
            self._ingest_status.setText("✗ Subject ID, Session ID, Raw path and Output path are all required.")
            return

        self._ingest_status.setStyleSheet(f"color:{style.ACCENT};")
        self._ingest_status.setText("⟳ Running ingestion…")

        # Emit signal to main app to run worker (avoids circular imports)
        import os
        prestus_root = os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__)
        )))
        from matlab_bridge import IngestWorker
        self._ingest_worker = IngestWorker(
            prestus_root, sub_id, ses_id, raw_path, out_path
        )
        self._ingest_worker.finished.connect(self._on_ingest_done)
        self._ingest_worker.error.connect(lambda msg: self._ingest_status.setText(f"✗ {msg}"))
        self._ingest_worker.start()

    def _on_ingest_done(self, success: bool, xml_path: str):
        if success:
            self._ingest_status.setStyleSheet(f"color:{style.SUCCESS};")
            self._ingest_status.setText(f"✓ Done. XML: {xml_path}")
        else:
            self._ingest_status.setStyleSheet(f"color:{style.MANDATORY};")
            self._ingest_status.setText("✗ Ingestion failed.")
