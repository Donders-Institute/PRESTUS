"""Results tab: acoustic/thermal NIfTI viewer, the self-contained HTML report, and the
output CSV, for the subject/study config currently open. Ported from the deprecated
functions/gui_py/tabs/tab_results.py onto this package's OrthoViewer/placement.Volume
stack instead of the old standalone nibabel viewer, and onto outputs.dir_output()
instead of a path set imperatively by the run tab."""

from __future__ import annotations

import csv
import os
import subprocess
import sys
from pathlib import Path

from PyQt6.QtCore import QUrl
from PyQt6.QtWidgets import (QComboBox, QHBoxLayout, QLabel, QPushButton, QTableWidget,
                             QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget)

from .. import placement as pl
from ..model import is_empty
from ..outputs import dir_output
from .ortho_view import OrthoViewer

ACOUSTIC_OVERLAYS = {
    "Intensity (ISPPA)": "*intensity*.nii*",
    "Pressure (peak)": "*pressure*.nii*",
    "Mechanical Index": "*MI*.nii*",
    "Phase (magnitude)": "*p_complex_magnitude*.nii*",
}
THERMAL_OVERLAYS = {
    "Max temperature": "*heating_end*.nii*",
    "CEM43": "*CEM43_end*.nii*",
}


def _find(base: Path, pattern: str) -> Path | None:
    matches = sorted((base / "nii").glob(pattern)) if (base / "nii").is_dir() else []
    matches += sorted(base.glob(pattern))
    return matches[-1] if matches else None


class _OverlayPanel(QWidget):
    """A viewer with a dropdown that picks which output NIfTI to overlay."""

    def __init__(self, patterns: dict[str, str], parent: QWidget | None = None):
        super().__init__(parent)
        self._base: Path | None = None
        self._patterns = patterns

        v = QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 0)
        tb = QHBoxLayout()
        tb.addWidget(QLabel("Overlay"))
        self.combo = QComboBox()
        self.combo.addItems(list(patterns.keys()))
        self.combo.currentTextChanged.connect(self._load_overlay)
        tb.addWidget(self.combo)
        tb.addStretch(1)
        v.addLayout(tb)
        self.viewer = OrthoViewer()
        v.addWidget(self.viewer, 1)

    def load(self, base: Path, anat: Path | None) -> None:
        self._base = base
        self.viewer.set_anatomy(anat)
        self._load_overlay(self.combo.currentText())

    def _load_overlay(self, label: str) -> None:
        if self._base is None:
            return
        pattern = self._patterns.get(label, "")
        self.viewer.set_overlay(_find(self._base, pattern) if pattern else None)


class ResultsPage(QWidget):
    """Acoustic / Thermal / Report / Data Table sub-tabs for the current study's output."""

    def __init__(self, window, parent: QWidget | None = None):
        super().__init__(parent)
        self.window = window

        v = QVBoxLayout(self)
        title = QLabel("Results")
        title.setObjectName("CardTitle")
        v.addWidget(title)

        tb = QHBoxLayout()
        refresh = QPushButton("Refresh")
        refresh.clicked.connect(self.refresh)
        tb.addWidget(refresh)
        open_html = QPushButton("Open HTML report")
        open_html.clicked.connect(self._open_html)
        tb.addWidget(open_html)
        open_folder = QPushButton("Open output folder")
        open_folder.clicked.connect(self._open_folder)
        tb.addWidget(open_folder)
        tb.addStretch(1)
        self.status = QLabel("No output loaded.")
        self.status.setStyleSheet("color: gray;")
        tb.addWidget(self.status)
        v.addLayout(tb)

        self.subtabs = QTabWidget()
        v.addWidget(self.subtabs, 1)
        self.acoustic = _OverlayPanel(ACOUSTIC_OVERLAYS)
        self.subtabs.addTab(self.acoustic, "Acoustic")
        self.thermal = _OverlayPanel(THERMAL_OVERLAYS)
        self.subtabs.addTab(self.thermal, "Thermal")
        self.report_tab, self.web = self._build_report_tab()
        self.subtabs.addTab(self.report_tab, "Report")
        self.table = QTableWidget()
        self.subtabs.addTab(self.table, "Data Table")

    def _build_report_tab(self) -> tuple[QWidget, object | None]:
        w = QWidget()
        vl = QVBoxLayout(w)
        vl.setContentsMargins(0, 0, 0, 0)
        try:
            from PyQt6.QtWebEngineWidgets import QWebEngineView
            web = QWebEngineView()
            vl.addWidget(web)
            return w, web
        except ImportError:
            note = QLabel("PyQt6-WebEngine not installed — use 'Open HTML report' to view it "
                          "in a browser instead.")
            note.setWordWrap(True)
            vl.addWidget(note)
            return w, None

    # ── model-derived paths ─────────────────────────────────────────────
    def _base(self) -> Path | None:
        base = dir_output(self.window.model)
        return base if base is not None and base.is_dir() else None

    def _anat(self) -> Path | None:
        m = self.window.model
        anat, pattern, sid = m.get("path.anat"), m.get("path.t1_pattern"), m.get("subject_id")
        if is_empty(anat) or not isinstance(sid, (int, float)) or isinstance(sid, bool):
            return None
        return pl.find_t1(anat, pattern or "", int(sid))

    def refresh(self) -> None:
        base = self._base()
        if base is None:
            self.status.setText("Output folder not found — check path.sim and subject_id, "
                                "and run the pipeline first.")
            return
        anat = self._anat()
        self.acoustic.load(base, anat)
        self.thermal.load(base, anat)
        self._load_report(base)
        self._load_table(base)
        self.status.setText(f"Loaded: {base}")

    def _reports(self, base: Path) -> list[Path]:
        return sorted(base.glob("*report*.html"))

    def _load_report(self, base: Path) -> None:
        reports = self._reports(base)
        if reports and self.web is not None:
            self.web.setUrl(QUrl.fromLocalFile(str(reports[-1])))

    def _load_table(self, base: Path) -> None:
        csvs = sorted(base.glob("*.csv"))
        if not csvs:
            self.table.setRowCount(0)
            self.table.setColumnCount(0)
            return
        with open(csvs[-1], newline="", encoding="utf-8") as fh:
            rows = list(csv.reader(fh))
        if not rows:
            return
        header, data = rows[0], rows[1:]
        self.table.setColumnCount(len(header))
        self.table.setRowCount(len(data))
        self.table.setHorizontalHeaderLabels(header)
        for ri, row in enumerate(data):
            for ci, value in enumerate(row):
                item = QTableWidgetItem(value)
                self.table.setItem(ri, ci, item)

    def _open_html(self) -> None:
        base = self._base()
        reports = self._reports(base) if base is not None else []
        if reports:
            _open_in_os(reports[-1])

    def _open_folder(self) -> None:
        base = self._base()
        if base is not None:
            _open_in_os(base)


def _open_in_os(path: Path) -> None:
    if sys.platform == "darwin":
        subprocess.run(["open", str(path)])
    elif sys.platform.startswith("win"):
        os.startfile(str(path))       # noqa: S606 (Windows-only API)
    else:
        subprocess.run(["xdg-open", str(path)])
