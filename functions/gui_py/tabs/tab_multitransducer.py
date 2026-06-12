"""Tab 10 – Multi-Transducer."""

from __future__ import annotations
from PyQt6.QtWidgets import (
    QWidget, QHBoxLayout, QPlainTextEdit, QPushButton, QLabel, QFileDialog,
)
from widgets import Dropdown
from tabs.base import BaseTab
import style


class TabMultiTransducer(BaseTab):
    def __init__(self, registry: dict, parent=None):
        super().__init__(registry, parent)
        r = registry

        self._add_section("Multi-Transducer Mode")
        self._add_note(
            "Primary transducer is defined in the Transducer tab. "
            "Secondary transducers are loaded from YAML config files listed below."
        )
        self._add_row(
            "Coupling mode",
            Dropdown("simulation.transducer_coupling",
                      ["coherent", "async"], "coherent", registry=r),
            note="coherent: all transducers in one k-Wave run (interference modelled).  "
                 "async: each simulated independently, intensities summed incoherently."
        )

        self._add_section("Secondary Transducer Configs")
        self._add_note(
            "Each YAML fully defines one additional transducer (type, frequency, position, "
            "focus, intensity). They are appended to parameters.transducer(2..N)."
        )

        self._ta = QPlainTextEdit()
        self._ta.setPlaceholderText("path/to/transducer2.yaml\npath/to/transducer3.yaml")
        self._ta.setMaximumHeight(160)
        registry["multitrans_configs_textarea"] = self._ta

        btn_add = QPushButton("+ Add YAML")
        btn_add.setStyleSheet(
            f"background:#dde3f0; color:{style.TEXT}; border-radius:4px; padding:4px 8px;"
        )
        btn_add.clicked.connect(self._browse)

        row = QWidget()
        hl = QHBoxLayout(row)
        hl.setContentsMargins(0, 0, 0, 0)
        hl.addWidget(self._ta, stretch=1)
        hl.addWidget(btn_add)

        self._layout.addWidget(QLabel("Config paths:"))
        self._layout.addWidget(row)
        self._add_note("One YAML path per line. Leave empty for single-transducer mode.")
        self._finish()

    def _browse(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select secondary transducer config YAML",
            filter="YAML (*.yaml *.yml);;All files (*)"
        )
        if path:
            cur = self._ta.toPlainText().strip()
            self._ta.setPlainText((cur + "\n" + path).strip())
