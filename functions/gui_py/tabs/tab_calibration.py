"""Tab 11 – Calibration."""

from __future__ import annotations
from PyQt6.QtWidgets import QPushButton, QLabel, QPlainTextEdit
from widgets import TextEdit, FloatEdit, Dropdown, Checkbox, CSVArrayEdit, BrowseEdit
from tabs.base import BaseTab
import style


class TabCalibration(BaseTab):
    def __init__(self, registry: dict, parent=None):
        super().__init__(registry, parent)
        r = registry

        # ── Input Data ────────────────────────────────────────────────────────
        self._add_section("Input Data")
        self._add_row("Equipment name",
                      TextEdit("calibration.equipment_name", "", registry=r),
                      note="Identifier matched against characterisation CSV files (e.g. CTX500).")
        self._add_row("Axial profile folder",
                      BrowseEdit("calibration.path_input_axial", "dir", registry=r),
                      note="Directory containing per-equipment axial hydrophone CSVs.")
        self._add_row("Phase table folder",
                      BrowseEdit("calibration.path_input_phase", "dir", registry=r),
                      note="Directory with manufacturer phase tables (Sonic Concepts / Imasonic).")

        # ── Calibration Targets ───────────────────────────────────────────────
        self._add_section("Calibration Targets")
        self._add_row("Focal depths (ep)",
                      CSVArrayEdit("calibration.focal_depths_wrt_exit_plane", "", registry=r),
                      note="mm — comma-separated target focal distances from the exit plane.")
        self._add_row("Desired intensities",
                      CSVArrayEdit("calibration.desired_intensities", "", registry=r),
                      note="W/cm² — comma-separated target Isppa values.")
        self._add_row("Add focal-distance offset",
                      Checkbox("calibration.add_FDO",
                               "Add bowl-offset to focal distance", False, registry=r),
                      note="Shifts focal target by (curvature radius − exit-plane distance).")

        # ── Optimisation ──────────────────────────────────────────────────────
        self._add_section("Optimisation Settings")
        self._add_row("Calibration mode",
                      Dropdown("calibration.Mode",
                                ["single_ref", "multi_depth"], "single_ref", registry=r),
                      note='"single_ref" optimises at one reference depth and extrapolates; '
                           '"multi_depth" fits all depths jointly.')
        self._add_row("Forward model",
                      Dropdown("calibration.ForwardModel",
                                ["oneil", "rayleigh"], "oneil", registry=r),
                      note='"oneil" is faster; "rayleigh" is more accurate for ring gaps.')
        self._add_row("Reference depth",
                      FloatEdit("calibration.RefDepth", None, registry=r),
                      note="mm — leave blank for median of available depths.")
        self._add_row("Initial velocity",
                      FloatEdit("calibration.initial_velocity", 0.05, registry=r),
                      note="m/s — starting guess for element normal velocity.")
        self._add_row("Upper velocity",
                      FloatEdit("calibration.opt_upper_velocity", 0.2, registry=r),
                      note="m/s — upper bound during optimisation.")
        self._add_row("Phase precession",
                      Dropdown("calibration.opt_phase_precession",
                                ["none", "linear", "monotonic"], "none", registry=r),
                      note='"none" = unconstrained per-element phases; '
                           '"linear"/"monotonic" add phase ramp constraints.')
        self._add_row("Regularisation λ",
                      FloatEdit("calibration.opt_regularization_lambda", 0, registry=r),
                      note="L2 penalty weight on hardware correction. 0 = none. Typical: 1e-4–1e-2.")
        self._add_row("Opt. limits",
                      CSVArrayEdit("calibration.opt_limits", "", registry=r),
                      note="mm — two-element [min max] range for profile error computation.")
        self._add_row("Profile weights",
                      FloatEdit("calibration.opt_weights", 0, registry=r),
                      note="0 = uniform; ≥1 = Gaussian centred on FLHM.")
        self._add_row("Opt. method",
                      Dropdown("calibration.opt_method",
                                ["FEXminimize", "GlobalSearch"], "FEXminimize", registry=r),
                      note='"FEXminimize" (no toolbox required); "GlobalSearch" requires Global Optimisation Toolbox.')
        self._add_row("Random seed",
                      FloatEdit("calibration.opt_seed", None, registry=r),
                      note="Integer seed for reproducibility. Leave blank for non-deterministic run.")

        # ── Simulation Control ────────────────────────────────────────────────
        self._add_section("Simulation Control")
        self._add_row("Run free-water sim",
                      Checkbox("calibration.run_free_water_sim",
                               "Simulate correction in free water", True, registry=r))
        self._add_row("Force kWaveArray",
                      Checkbox("calibration.force_kwavearray",
                               "Use kWaveArray transducer backend", False, registry=r))
        self._add_row("Axisymmetric 2D",
                      Checkbox("calibration.axisymmetric2D",
                               "Use 2-D axisymmetric grid (annular only)", False, registry=r))
        self._add_row("Amplitude validation",
                      Dropdown("calibration.opt_amp_validation",
                                ["always", "initial", "final", "none"], "final", registry=r),
                      note="When to run a validation simulation after amplitude calibration.")

        # ── Per-Element Phase Corrections ─────────────────────────────────────
        self._add_section("Per-Element Phase Corrections")
        self._add_row("Correction offsets",
                      CSVArrayEdit("calibration.elem_phase_correction_deg", "", registry=r),
                      note="deg — comma-separated hardware phase offsets per element.")
        self._add_row("Save recovered offsets",
                      Checkbox("calibration.save_elem_correction",
                               "Write recovered corrections back to transducer YAML",
                               False, registry=r))

        # ── Output ────────────────────────────────────────────────────────────
        self._add_section("Output")
        self._add_row("Output folder",
                      BrowseEdit("calibration.path_output", "dir", registry=r))
        self._add_row("Profiles output folder",
                      BrowseEdit("calibration.path_output_profiles", "dir", registry=r))
        self._add_row("Calibrated CSV filename",
                      TextEdit("calibration.filename_calibrated_CSV", "", registry=r),
                      note="Output filename for the calibrated intensity profile CSV.")
        self._add_row("Save in calibration folder",
                      Checkbox("calibration.save_in_calibration_folder",
                               "Route all outputs to the calibration output folder",
                               False, registry=r))
        self._add_row("BabelBrain export path",
                      BrowseEdit("calibration.ExportBabelBrain", "file", registry=r),
                      note="Optional .h5 output path for BabelBrain element weights.")

        # ── Run ───────────────────────────────────────────────────────────────
        self._add_section("Run Calibration")
        self._status_lbl = QLabel("Ready.")
        self._status_lbl.setStyleSheet(f"color:{style.TEXT_SUB}; font-size:{style.FS_SM}px;")

        btn = QPushButton("▶   Run Calibration")
        btn.setStyleSheet(
            "background:#2d8f57; color:white; font-weight:bold; "
            f"font-size:{style.FS_LG}px; border-radius:4px; padding:10px 16px;"
        )
        btn.clicked.connect(self._run_calibration)
        self._layout.addWidget(btn)
        self._layout.addWidget(self._status_lbl)

        # Log area
        self._log = QPlainTextEdit()
        self._log.setReadOnly(True)
        self._log.setMaximumHeight(160)
        self._log.setStyleSheet(
            "background:#191f18; color:#ccf0c0; "
            "font-family:'Courier New',monospace; font-size:11px; border:none;"
        )
        self._layout.addWidget(self._log)
        self._finish()

    def _run_calibration(self):
        import os, tempfile
        self._status_lbl.setStyleSheet(f"color:{style.ACCENT};")
        self._status_lbl.setText("⟳ Running calibration…")
        self._log.clear()

        prestus_root = os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__)
        )))
        # Write a minimal config for calibration
        try:
            import yaml
            from param_io import collect_params
            params = collect_params(self.registry)
            with tempfile.NamedTemporaryFile(suffix=".yaml", delete=False, mode="w") as f:
                yaml.dump(params, f, default_flow_style=False)
                tmp = f.name
        except Exception as exc:
            self._status_lbl.setText(f"✗ Config error: {exc}")
            return

        from matlab_bridge import CalibrationWorker
        self._worker = CalibrationWorker(tmp, prestus_root)
        self._worker.log_line.connect(lambda line: self._log.appendPlainText(line))
        self._worker.finished.connect(self._on_cal_done)
        self._worker.error.connect(lambda msg: self._status_lbl.setText(f"✗ {msg}"))
        self._worker.start()

    def _on_cal_done(self, success: bool):
        if success:
            self._status_lbl.setStyleSheet(f"color:{style.SUCCESS};")
            self._status_lbl.setText("✓ Calibration complete.")
        else:
            self._status_lbl.setStyleSheet(f"color:{style.MANDATORY};")
            self._status_lbl.setText("✗ Calibration failed.")
