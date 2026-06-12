"""Tab 14 – Results.

Shows:
  • Acoustic sub-tab: 3-plane NIfTI viewer (anatomy + intensity/pressure overlay)
  • Thermal sub-tab:  3-plane NIfTI viewer (anatomy + heating/CEM43 overlay)
  • Report sub-tab:   embedded HTML report (via QWebEngineView if available, else open browser)
  • Data Table sub-tab: CSV output table
"""

from __future__ import annotations
import os
import glob
from pathlib import Path

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QTabWidget,
    QPushButton, QLabel, QComboBox, QTableWidget, QTableWidgetItem,
    QHeaderView, QSizePolicy,
)
from PyQt6.QtCore import Qt
import style


def _find_nifti(output_dir: str, pattern: str) -> str | None:
    """Return the most recent NIfTI matching a glob pattern."""
    matches = sorted(glob.glob(os.path.join(output_dir, "nii", pattern)))
    if matches:
        return matches[-1]
    matches = sorted(glob.glob(os.path.join(output_dir, pattern)))
    return matches[-1] if matches else None


class _NiftiPanel(QWidget):
    """A viewer panel: dropdown to choose overlay file + 3-plane NiftiViewer."""

    def __init__(self, overlay_patterns: dict[str, str], parent=None):
        super().__init__(parent)
        self._output_dir = None
        self._patterns   = overlay_patterns  # label → glob pattern

        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(4)

        # Toolbar
        tb = QHBoxLayout()
        tb.addWidget(QLabel("Overlay:"))
        self._overlay_cb = QComboBox()
        self._overlay_cb.addItems(list(overlay_patterns.keys()))
        self._overlay_cb.setFixedWidth(220)
        self._overlay_cb.currentTextChanged.connect(self._on_overlay_changed)
        tb.addWidget(self._overlay_cb)
        tb.addStretch()
        layout.addLayout(tb)

        try:
            from nifti_viewer import NiftiViewer
            self._viewer = NiftiViewer()
            layout.addWidget(self._viewer, stretch=1)
        except Exception:
            self._viewer = None
            layout.addWidget(QLabel(
                "NIfTI viewer unavailable — install matplotlib and nibabel."
            ))

    def load(self, output_dir: str):
        self._output_dir = output_dir
        if not self._viewer:
            return
        # Load anatomy
        anat = _find_nifti(output_dir, "*T1w.nii*")
        if not anat:
            anat = _find_nifti(output_dir, "*.nii*")
        if anat:
            self._viewer.load_anatomy(anat)
        self._on_overlay_changed(self._overlay_cb.currentText())

    def _on_overlay_changed(self, label: str):
        if not self._output_dir or not self._viewer:
            return
        pattern = self._patterns.get(label, "")
        path = _find_nifti(self._output_dir, pattern) if pattern else None
        if path:
            self._viewer.load_overlay(path)
        else:
            self._viewer.clear_overlay()


class TabResults(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._output_dir = None

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 8, 12, 12)
        layout.setSpacing(6)

        # ── Toolbar ───────────────────────────────────────────────────────────
        tb = QHBoxLayout()
        btn_refresh = QPushButton("↻  Refresh")
        btn_refresh.setStyleSheet(
            f"background:#dde3f0; color:{style.TEXT}; border-radius:4px; padding:5px 10px;"
        )
        btn_refresh.clicked.connect(self.refresh)
        tb.addWidget(btn_refresh)

        btn_html = QPushButton("🌐  Open HTML Report")
        btn_html.setStyleSheet(
            f"background:#dde3f0; color:{style.TEXT}; border-radius:4px; padding:5px 10px;"
        )
        btn_html.clicked.connect(self._open_html)
        tb.addWidget(btn_html)

        btn_folder = QPushButton("📁  Open Output Folder")
        btn_folder.setStyleSheet(
            f"background:#dde3f0; color:{style.TEXT}; border-radius:4px; padding:5px 10px;"
        )
        btn_folder.clicked.connect(self._open_folder)
        tb.addWidget(btn_folder)

        tb.addStretch()
        self._status_lbl = QLabel("No output loaded.")
        self._status_lbl.setStyleSheet(f"color:{style.TEXT_SUB}; font-size:{style.FS_SM}px;")
        tb.addWidget(self._status_lbl)
        layout.addLayout(tb)

        # ── Sub-tabs ──────────────────────────────────────────────────────────
        self._subtabs = QTabWidget()
        self._subtabs.setTabPosition(QTabWidget.TabPosition.North)
        layout.addWidget(self._subtabs, stretch=1)

        # Acoustic
        self._acoustic_panel = _NiftiPanel({
            "Intensity (ISPPA)":   "*intensity*.nii*",
            "Pressure (peak)":     "*pressure*.nii*",
            "Mechanical Index":    "*MI*.nii*",
            "Phase (magnitude)":   "*p_complex_magnitude*.nii*",
        })
        self._subtabs.addTab(self._acoustic_panel, "Acoustic")

        # Thermal
        self._thermal_panel = _NiftiPanel({
            "Max temperature":     "*heating_end*.nii*",
            "CEM43":               "*CEM43_end*.nii*",
        })
        self._subtabs.addTab(self._thermal_panel, "Thermal")

        # Report
        self._report_tab = self._build_report_tab()
        self._subtabs.addTab(self._report_tab, "Report")

        # CSV
        self._csv_tab = self._build_csv_tab()
        self._subtabs.addTab(self._csv_tab, "Data Table")

    # ── Report tab ────────────────────────────────────────────────────────────

    def _build_report_tab(self) -> QWidget:
        w = QWidget()
        vl = QVBoxLayout(w)
        vl.setContentsMargins(0, 0, 0, 0)
        # Try QWebEngineView
        try:
            from PyQt6.QtWebEngineWidgets import QWebEngineView
            self._web = QWebEngineView()
            vl.addWidget(self._web)
        except ImportError:
            self._web = None
            note = QLabel(
                "PyQt6-WebEngine not installed — click 'Open HTML Report' to view in browser."
            )
            note.setAlignment(Qt.AlignmentFlag.AlignCenter)
            note.setStyleSheet(f"color:{style.TEXT_SUB}; font-size:{style.FS_SM}px;")
            vl.addWidget(note)
        return w

    # ── CSV tab ───────────────────────────────────────────────────────────────

    def _build_csv_tab(self) -> QWidget:
        w = QWidget()
        vl = QVBoxLayout(w)
        vl.setContentsMargins(4, 4, 4, 4)
        self._csv_table = QTableWidget()
        self._csv_table.setStyleSheet(
            f"QTableWidget {{ background:{style.BG_EDIT}; gridline-color:{style.BORDER}; }}"
            f"QHeaderView::section {{ background:{style.BG_FIG}; border:1px solid {style.BORDER}; "
            f"padding:4px; font-size:{style.FS_SM}px; }}"
        )
        self._csv_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        vl.addWidget(self._csv_table)
        return w

    # ── Public API ────────────────────────────────────────────────────────────

    def set_output_dir(self, path: str):
        self._output_dir = path

    def refresh(self):
        if not self._output_dir or not os.path.isdir(self._output_dir):
            self._status_lbl.setText("Output directory not set or does not exist.")
            return
        self._acoustic_panel.load(self._output_dir)
        self._thermal_panel.load(self._output_dir)
        self._load_report()
        self._load_csv()
        self._status_lbl.setText(f"Loaded: {self._output_dir}")

    def _load_report(self):
        if not self._output_dir:
            return
        reports = sorted(glob.glob(os.path.join(self._output_dir, "*report*.html")))
        if not reports:
            return
        if self._web is not None:
            from PyQt6.QtCore import QUrl
            self._web.setUrl(QUrl.fromLocalFile(reports[-1]))

    def _load_csv(self):
        if not self._output_dir:
            return
        csvs = sorted(glob.glob(os.path.join(self._output_dir, "*.csv")))
        if not csvs:
            return
        try:
            import csv
            with open(csvs[-1], newline="") as f:
                reader = csv.reader(f)
                rows = list(reader)
            if not rows:
                return
            headers = rows[0]
            data    = rows[1:]
            self._csv_table.setColumnCount(len(headers))
            self._csv_table.setRowCount(len(data))
            self._csv_table.setHorizontalHeaderLabels(headers)
            for ri, row in enumerate(data):
                for ci, val in enumerate(row):
                    item = QTableWidgetItem(val)
                    item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                    self._csv_table.setItem(ri, ci, item)
        except Exception:
            pass

    def _open_html(self):
        if not self._output_dir:
            return
        reports = sorted(glob.glob(os.path.join(self._output_dir, "*report*.html")))
        if reports:
            import subprocess, sys
            if sys.platform == "darwin":
                subprocess.run(["open", reports[-1]])
            elif sys.platform.startswith("win"):
                os.startfile(reports[-1])
            else:
                subprocess.run(["xdg-open", reports[-1]])

    def _open_folder(self):
        if not self._output_dir or not os.path.isdir(self._output_dir):
            return
        import subprocess, sys
        if sys.platform == "darwin":
            subprocess.run(["open", self._output_dir])
        elif sys.platform.startswith("win"):
            os.startfile(self._output_dir)
        else:
            subprocess.run(["xdg-open", self._output_dir])
