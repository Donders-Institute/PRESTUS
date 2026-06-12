"""Main application window for the PRESTUS Python GUI."""

from __future__ import annotations
import os
import sys
import tempfile
from pathlib import Path

from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QTabWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QFileDialog, QMessageBox, QSizePolicy,
    QSplitter,
)
from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QFont

import style
import param_io
from matlab_bridge import MatlabWorker

from tabs.tab_io            import TabIO
from tabs.tab_simulation    import TabSimulation
from tabs.tab_transducer    import TabTransducer
from tabs.tab_grid          import TabGrid
from tabs.tab_medium        import TabMedium
from tabs.tab_thermal       import TabThermal
from tabs.tab_hpc           import TabHPC
from tabs.tab_advanced      import TabAdvanced
from tabs.tab_placement     import TabPlacement
from tabs.tab_multitransducer import TabMultiTransducer
from tabs.tab_calibration   import TabCalibration
from tabs.tab_tx_library    import TabTxLibrary
from tabs.tab_run           import TabRun
from tabs.tab_results       import TabResults


def _get_prestus_root() -> str:
    """Return the PRESTUS repo root (two levels up from functions/gui_py/)."""
    return str(Path(__file__).resolve().parent.parent.parent)


class PrestusApp(QMainWindow):
    def __init__(self, init_params: dict | None = None):
        super().__init__()
        self._prestus_root = _get_prestus_root()
        self._registry: dict = {}
        self._last_params: dict | None = None
        self._config_path: str | None = None
        self._worker: MatlabWorker | None = None

        self.setWindowTitle("PRESTUS")
        self.resize(900, 720)
        self.setStyleSheet(style.QSS_BASE)

        central = QWidget()
        self.setCentralWidget(central)
        root_vl = QVBoxLayout(central)
        root_vl.setContentsMargins(0, 0, 0, 0)
        root_vl.setSpacing(0)

        # ── Header ────────────────────────────────────────────────────────────
        root_vl.addWidget(self._build_header())

        # ── Tab group ─────────────────────────────────────────────────────────
        self._tabs = QTabWidget()
        self._tabs.setTabPosition(QTabWidget.TabPosition.West)
        self._tabs.setDocumentMode(True)
        root_vl.addWidget(self._tabs, stretch=1)

        r = self._registry
        pr = self._prestus_root

        self._tab_io           = TabIO(r)
        self._tab_simulation   = TabSimulation(r)
        self._tab_transducer   = TabTransducer(r, pr)
        self._tab_grid         = TabGrid(r)
        self._tab_medium       = TabMedium(r)
        self._tab_thermal      = TabThermal(r)
        self._tab_hpc          = TabHPC(r)
        self._tab_advanced     = TabAdvanced(r)
        self._tab_placement    = TabPlacement(r)
        self._tab_multitrans   = TabMultiTransducer(r)
        self._tab_calibration  = TabCalibration(r)
        self._tab_tx_library   = TabTxLibrary(r, pr)
        self._tab_run          = TabRun(self._cb_run)
        self._tab_results      = TabResults()

        tabs = [
            ("  I/O & Paths   ",      self._tab_io),
            ("  Simulation    ",      self._tab_simulation),
            ("  Transducer    ",      self._tab_transducer),
            ("  Grid          ",      self._tab_grid),
            ("  Medium        ",      self._tab_medium),
            ("  Thermal       ",      self._tab_thermal),
            ("  HPC           ",      self._tab_hpc),
            ("  Advanced      ",      self._tab_advanced),
            ("  Placement     ",      self._tab_placement),
            ("  Multi-Transducer  ",  self._tab_multitrans),
            ("  Calibration   ",      self._tab_calibration),
            ("  Tx Library    ",      self._tab_tx_library),
            ("  ▶  Run        ",      self._tab_run),
            ("  Results       ",      self._tab_results),
        ]
        for title, widget in tabs:
            self._tabs.addTab(widget, title)

        # Wire thermal toggle
        heating_cb = r.get("modules.run_heating_sims")
        if heating_cb:
            heating_cb.stateChanged.connect(
                lambda state: self._tab_thermal.set_timing_enabled(bool(state))
            )

        # Load defaults on startup
        QTimer.singleShot(0, self._load_defaults)

        if init_params:
            QTimer.singleShot(100, lambda: param_io.apply_params(
                init_params, r, getattr(self._tab_medium, "table", None)
            ))

    # ── Header ────────────────────────────────────────────────────────────────

    def _build_header(self) -> QWidget:
        hdr = QWidget()
        hdr.setFixedHeight(46)
        hdr.setStyleSheet(f"background:{style.BG_DARK};")
        hl = QHBoxLayout(hdr)
        hl.setContentsMargins(14, 0, 14, 0)
        hl.setSpacing(12)

        title = QLabel("PRESTUS")
        title.setStyleSheet(
            f"color:white; font-size:{style.FS_XL}px; font-weight:bold;"
        )
        hl.addWidget(title)

        subtitle = QLabel("PREprocessing & Simulations for Transcranial Ultrasound Stimulation")
        subtitle.setStyleSheet(f"color:#a7bfc4; font-size:{style.FS_SM}px;")
        hl.addWidget(subtitle)

        hl.addStretch()

        self._config_lbl = QLabel("Config: defaults")
        self._config_lbl.setStyleSheet(f"color:#7a9097; font-size:{style.FS_SM}px;")
        hl.addWidget(self._config_lbl)

        for text, cb, bg in (
            ("⬆ Load",     self._cb_load_yaml,  "#0b7b99"),
            ("⬇ Save",     self._cb_save_yaml,  "#0b7b99"),
            ("↺ Defaults", self._load_defaults, "#075266"),
        ):
            btn = QPushButton(text)
            btn.setStyleSheet(
                f"background:{bg}; color:white; border-radius:4px; padding:5px 12px;"
                f"font-size:{style.FS_SM}px;"
            )
            btn.clicked.connect(cb)
            hl.addWidget(btn)
        return hdr

    # ── Load / save YAML ──────────────────────────────────────────────────────

    def _load_defaults(self):
        default_yaml = os.path.join(self._prestus_root, "config", "config_default.yaml")
        if os.path.isfile(default_yaml):
            self._load_yaml(default_yaml, label="defaults")

    def _cb_load_yaml(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load YAML config",
            filter="YAML (*.yaml *.yml);;All files (*)"
        )
        if path:
            self._load_yaml(path)

    def _load_yaml(self, path: str, label: str | None = None):
        try:
            import yaml
            with open(path) as f:
                params = yaml.safe_load(f) or {}
            self._last_params = params
            param_io.apply_params(
                params, self._registry,
                medium_table=getattr(self._tab_medium, "table", None)
            )
            display = label or os.path.basename(path)
            self._config_lbl.setText(f"Config: {display}")
            self._config_path = path
            # Sync thermal toggle
            self._sync_thermal_toggle()
        except Exception as exc:
            QMessageBox.critical(self, "Load Error", str(exc))

    def _cb_save_yaml(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save parameters as YAML",
            "config_study.yaml",
            "YAML (*.yaml *.yml);;All files (*)"
        )
        if not path:
            return
        try:
            import yaml
            params = param_io.collect_params(
                self._registry, self._last_params,
                medium_table=getattr(self._tab_medium, "table", None)
            )
            with open(path, "w") as f:
                yaml.dump(params, f, default_flow_style=False, allow_unicode=True)
            QMessageBox.information(self, "Saved", f"Saved to:\n{path}")
        except Exception as exc:
            QMessageBox.critical(self, "Save Error", str(exc))

    # ── Run simulation ────────────────────────────────────────────────────────

    def _cb_run(self):
        try:
            params = param_io.collect_params(
                self._registry, self._last_params,
                medium_table=getattr(self._tab_medium, "table", None)
            )
        except Exception as exc:
            QMessageBox.critical(self, "Parameter Error", str(exc))
            return

        # Write temp YAML
        try:
            import yaml
            tmp_f = tempfile.NamedTemporaryFile(
                suffix=".yaml", delete=False, mode="w"
            )
            yaml.dump(params, tmp_f, default_flow_style=False)
            tmp_f.close()
            tmp_path = tmp_f.name
        except Exception as exc:
            QMessageBox.critical(self, "YAML Error", f"Could not write config: {exc}")
            return

        # Store output dir for results tab
        out_dir = param_io.get_nested(params, "io.dir_output") or ""
        self._tab_results.set_output_dir(str(out_dir))

        self._tab_run.set_running(True)
        self._tab_run.set_status("⟳  Running…", style.ACCENT)
        self._tab_run.clear_log()
        self._tab_run.append_log("▶ Starting simulation…")

        self._worker = MatlabWorker(tmp_path, self._prestus_root)
        self._worker.log_line.connect(self._tab_run.append_log)
        self._worker.log_line.connect(self._update_stage)
        self._worker.finished.connect(self._on_sim_finished)
        self._worker.error.connect(
            lambda msg: self._tab_run.set_status(f"✗  {msg}", style.MANDATORY)
        )
        self._worker.start()

    def _on_sim_finished(self, success: bool):
        self._tab_run.set_running(False)
        if success:
            self._tab_run.set_status("✓  Done", style.SUCCESS)
            self._tab_run.append_log("✓ Simulation completed.")
            self._tab_results.refresh()
        else:
            self._tab_run.set_status("✗  Error — see log", style.MANDATORY)

    def _update_stage(self, line: str):
        """Parse pipeline section banners from log and update status label."""
        if not hasattr(self, "_last_was_banner"):
            self._last_was_banner = False
        if "========" in line:
            self._last_was_banner = True
        elif self._last_was_banner and line.strip():
            self._tab_run.set_status(f"⟳  {line.strip()}", style.ACCENT)
            self._last_was_banner = False
        else:
            self._last_was_banner = False

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _sync_thermal_toggle(self):
        heating_cb = self._registry.get("modules.run_heating_sims")
        if heating_cb:
            self._tab_thermal.set_timing_enabled(heating_cb.isChecked())
