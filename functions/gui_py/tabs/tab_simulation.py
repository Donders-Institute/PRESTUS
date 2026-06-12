"""Tab 2 – Simulation."""

from __future__ import annotations
from PyQt6.QtWidgets import (
    QWidget, QHBoxLayout, QPushButton, QFileDialog, QLabel,
)
from widgets import Dropdown, Checkbox, CSVArrayEdit, TextEdit
from tabs.base import BaseTab
import style


class TabSimulation(BaseTab):
    def __init__(self, registry: dict, parent=None):
        super().__init__(registry, parent)
        r = registry

        # ── Simulation Type ───────────────────────────────────────────────────
        self._add_section("Simulation Type")
        self._add_row("* Medium",
                      Dropdown("simulation.medium",
                                ["layered", "water", "phantom"], "layered", registry=r))
        self._add_row("* Code backend",
                      Dropdown("simulation.code_type",
                                ["matlab_gpu", "matlab_cpu", "cpp_gpu", "cpp_cpu"],
                                "matlab_gpu", registry=r))
        self._add_row("Precision",
                      Dropdown("simulation.precision",
                                ["single", "double"], "single", registry=r))
        self._add_row("Platform",
                      Dropdown("platform",
                                ["auto", "matlab", "slurm", "qsub"], "auto", registry=r))
        self._add_row("Acoustic backend",
                      Dropdown("simulation.acoustic_backend",
                                ["kwave", "babelbrain", "jwave"], "kwave", registry=r))
        self._add_row("Thermal backend",
                      Dropdown("simulation.thermal_backend",
                                ["kwave", "babelbrain"], "kwave", registry=r))
        self._add_row("Interactive mode",
                      Checkbox("simulation.interactive",
                               "Ask for user input and plot evolving figures (desktop only)",
                               False, registry=r))
        self._add_row("Debug mode",
                      Checkbox("simulation.debug",
                               "Enable additional intermediate diagnostic outputs",
                               False, registry=r))
        self._add_row("Uncertainty mode",
                      Checkbox("simulation.uncertainty",
                               "Run default / liberal / conservative variants + combined report",
                               False, registry=r))

        # ── Pipeline Modules ─────────────────────────────────────────────────
        self._add_section("Pipeline Modules")
        self._add_note(
            "Toggle individual pipeline stages. "
            "Disable stages only if you know what you are doing."
        )
        modules = [
            ("modules.run_grid_setup",         "Grid setup & head preprocessing",       True),
            ("modules.run_medium_setup",        "Medium acoustic property mapping",       True),
            ("modules.run_source_setup",        "Acoustic source setup",                  True),
            ("modules.run_acoustic_sims",       "Acoustic simulation (k-Wave)",           True),
            ("modules.run_acoustic_analysis",   "Acoustic analysis (ISPPA, MI, etc.)",   True),
            ("modules.run_heating_sims",        "Thermal simulation",                     False),
            ("modules.run_thermal_analysis",    "Thermal analysis (CEM43, maxT)",         True),
            ("modules.run_nifti_creation",      "NIfTI export",                           True),
            ("modules.run_posthoc_water_sims",  "Post-hoc free-water simulation",         True),
            ("modules.generate_report",         "Generate HTML report",                   True),
            ("modules.segmentation_only",       "Segmentation only — stop after SimNIBS", False),
        ]
        for tag, label, default in modules:
            self._add_row_wide(Checkbox(tag, label, default, registry=r))

        # ── Multi-ISPPA Sweep ─────────────────────────────────────────────────
        self._add_section("Multi-ISPPA Intensity Sweep")
        self._add_row("Target Isppa (W/cm²)",
                      CSVArrayEdit("transducer.target_isppa_wcm2", "", registry=r),
                      note="One value for a single run; comma-separated list triggers intensity "
                           "sweep (separate thermal sim per value, shared acoustic sim).")

        # ── Sequential Simulations ───────────────────────────────────────────
        self._add_section("Sequential Simulations")
        self._add_row_wide(
            Checkbox("options.sequential.enabled",
                     "Chain follow-up sonications (thermal state inherited between runs)",
                     False, registry=r)
        )
        self._add_note(
            "Enter one YAML config path per line. Each follow-up inherits the "
            "temperature and CEM43 maps from the preceding run."
        )

        from PyQt6.QtWidgets import QPlainTextEdit, QVBoxLayout
        self._seq_ta = QPlainTextEdit()
        self._seq_ta.setPlaceholderText("path/to/followup1.yaml")
        self._seq_ta.setMaximumHeight(120)
        registry["sequential_configs_textarea"] = self._seq_ta

        btn_add = QPushButton("+ Add YAML")
        btn_add.setProperty("class", "secondary")
        btn_add.setStyleSheet(
            f"background:#dde3f0; color:{style.TEXT}; border-radius:4px; padding:4px 8px;"
        )
        btn_add.clicked.connect(self._browse_seq)

        w = QWidget()
        vl = QHBoxLayout(w)
        vl.setContentsMargins(0, 0, 0, 0)
        vl.addWidget(self._seq_ta, stretch=1)
        vl.addWidget(btn_add)
        self._layout.addWidget(QLabel("Follow-up configs:"))
        self._layout.addWidget(w)

        self._finish()

    def _browse_seq(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select follow-up config YAML",
            filter="YAML (*.yaml *.yml);;All files (*)"
        )
        if path:
            cur = self._seq_ta.toPlainText().strip()
            self._seq_ta.setPlainText((cur + "\n" + path).strip())
