"""Tab 12 – Transducer Library browser."""

from __future__ import annotations
import os
from pathlib import Path

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter,
    QListWidget, QListWidgetItem, QPushButton, QLabel,
    QTextEdit,
)
from PyQt6.QtCore import Qt
from tabs.base import BaseTab
import style


class TabTxLibrary(QWidget):
    """Two-column browser: calibration YAML list on left, raw content on right."""

    def __init__(self, registry: dict, prestus_root: str, parent=None):
        super().__init__(parent)
        self.registry     = registry
        self._prestus_root = prestus_root
        self._lib_path     = os.path.join(prestus_root, "config", "calibration")

        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 12, 16, 16)
        layout.setSpacing(8)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        layout.addWidget(splitter)

        # ── Left panel ────────────────────────────────────────────────────────
        left = QWidget()
        left_vl = QVBoxLayout(left)
        left_vl.setContentsMargins(0, 0, 0, 0)
        left_vl.setSpacing(6)

        hdr = QLabel("Library entries")
        hdr.setStyleSheet(f"font-weight:bold; color:{style.TEXT};")
        left_vl.addWidget(hdr)

        self._list = QListWidget()
        self._list.currentTextChanged.connect(self._show_entry)
        left_vl.addWidget(self._list, stretch=1)

        btn_refresh = QPushButton("↻  Refresh")
        btn_refresh.setStyleSheet(
            f"background:#dde3f0; color:{style.TEXT}; border-radius:4px; padding:5px 10px;"
        )
        btn_refresh.clicked.connect(self._refresh)
        left_vl.addWidget(btn_refresh)
        splitter.addWidget(left)

        # ── Right panel ───────────────────────────────────────────────────────
        right = QWidget()
        right_vl = QVBoxLayout(right)
        right_vl.setContentsMargins(0, 0, 0, 0)
        right_vl.setSpacing(6)

        self._entry_name = QLabel("Select an entry to inspect")
        self._entry_name.setStyleSheet(f"font-weight:bold; color:{style.TEXT};")
        right_vl.addWidget(self._entry_name)

        self._detail = QTextEdit()
        self._detail.setReadOnly(True)
        self._detail.setStyleSheet(
            f"font-family:'Courier New',monospace; font-size:12px; "
            f"background:{style.BG_EDIT}; border:1px solid {style.BORDER};"
        )
        right_vl.addWidget(self._detail, stretch=1)

        self._btn_load = QPushButton("⬆  Load into active config")
        self._btn_load.setStyleSheet(
            f"background:{style.ACCENT}; color:white; border-radius:4px; padding:5px 10px;"
        )
        self._btn_load.setEnabled(False)
        self._btn_load.clicked.connect(self._load_entry)
        right_vl.addWidget(self._btn_load)

        self._status = QLabel("")
        self._status.setStyleSheet(f"color:{style.TEXT_SUB}; font-size:{style.FS_SM}px;")
        right_vl.addWidget(self._status)

        splitter.addWidget(right)
        splitter.setSizes([220, 580])

        self._refresh()

    def _refresh(self):
        self._list.clear()
        if not os.path.isdir(self._lib_path):
            self._list.addItem("(library folder not found)")
            self._status.setText("Library folder missing: config/calibration/")
            return
        files = sorted(Path(self._lib_path).glob("*.yaml"))
        if not files:
            self._list.addItem("(no entries yet)")
            self._status.setText("Run update_transducer_library() to populate.")
        else:
            for f in files:
                self._list.addItem(f.stem)
            self._status.setText(f"{len(files)} entr{'y' if len(files)==1 else 'ies'} found.")

    def _show_entry(self, name: str):
        if not name or name.startswith("("):
            self._entry_name.setText("Select an entry to inspect")
            self._detail.clear()
            self._btn_load.setEnabled(False)
            return
        yaml_path = os.path.join(self._lib_path, name + ".yaml")
        if not os.path.isfile(yaml_path):
            self._detail.setPlainText("File not found.")
            self._btn_load.setEnabled(False)
            return
        with open(yaml_path) as f:
            self._detail.setPlainText(f.read())
        self._entry_name.setText(name)
        self._btn_load.setEnabled(True)

    def _load_entry(self):
        name = self._list.currentItem()
        if name is None or name.text().startswith("("):
            return
        combo = name.text()
        yaml_path = os.path.join(self._lib_path, combo + ".yaml")
        if not os.path.isfile(yaml_path):
            self._status.setText(f"File not found: {yaml_path}")
            return
        try:
            import yaml
            from param_io import apply_params
            with open(yaml_path) as f:
                lib = yaml.safe_load(f)
            if lib is None:
                return
            apply_params(lib, self.registry)
            self._status.setText(f"✓ Loaded '{combo}' into active config.")
        except Exception as exc:
            self._status.setText(f"✗ Error: {exc}")
