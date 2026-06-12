"""Tab 13 – Run."""

from __future__ import annotations
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QPlainTextEdit,
)
from PyQt6.QtCore import Qt
import style


class TabRun(QWidget):
    """Run panel: big green button, status label, terminal-style log."""

    def __init__(self, run_callback, parent=None):
        super().__init__(parent)
        self._run_callback = run_callback

        layout = QVBoxLayout(self)
        layout.setContentsMargins(20, 16, 20, 16)
        layout.setSpacing(8)

        # ── Run button ────────────────────────────────────────────────────────
        self.btn_run = QPushButton("▶   Run Simulation")
        self.btn_run.setStyleSheet(
            f"background:{style.SUCCESS}; color:white; font-weight:bold; "
            f"font-size:{style.FS_LG}px; border-radius:6px; padding:14px 20px;"
        )
        self.btn_run.clicked.connect(self._on_run)
        layout.addWidget(self.btn_run)

        # ── Status label ──────────────────────────────────────────────────────
        self.status_lbl = QLabel("Ready.")
        self.status_lbl.setStyleSheet(f"color:{style.TEXT_SUB}; font-size:{style.FS}px;")
        layout.addWidget(self.status_lbl)

        # ── Log area ──────────────────────────────────────────────────────────
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setPlaceholderText("Simulation log will appear here.")
        self.log.setStyleSheet(
            "background:#131a12; color:#ccf0c0; "
            "font-family:'Courier New',monospace; font-size:12px; border:none; "
            "border-radius:4px;"
        )
        layout.addWidget(self.log, stretch=1)

    def _on_run(self):
        self._run_callback()

    def set_status(self, text: str, color: str = style.TEXT_SUB):
        self.status_lbl.setText(text)
        self.status_lbl.setStyleSheet(f"color:{color}; font-size:{style.FS}px;")

    def append_log(self, text: str):
        self.log.appendPlainText(text)
        sb = self.log.verticalScrollBar()
        sb.setValue(sb.maximum())

    def clear_log(self):
        self.log.clear()

    def set_running(self, running: bool):
        self.btn_run.setEnabled(not running)
        if running:
            self.btn_run.setText("⟳   Running…")
        else:
            self.btn_run.setText("▶   Run Simulation")
