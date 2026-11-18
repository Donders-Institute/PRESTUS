"""Expected load at the bottom of the tabs whose settings change it (workload.predict):
one line with the simulation time per stage, total and GPU memory; the full prediction
as tooltip; a link to Run ▸ Plan. Recomputed shortly after the last edit, because the
prediction checks for existing inputs on disk."""

from __future__ import annotations

import html

from PyQt6.QtCore import Qt, QTimer, pyqtSignal
from PyQt6.QtWidgets import QFrame, QHBoxLayout, QLabel

from ..model import ConfigModel
from ..workload import predict

DELAY_MS = 300


class LoadFooter(QFrame):
    """'Expected run load: …' strip under a tab's settings: the study's simulation run, never
    the transducer calibration fit."""

    details = pyqtSignal()          # 'Details' clicked: show Run ▸ Plan

    def __init__(self, model: ConfigModel, parent=None):
        super().__init__(parent)
        self.model = model
        self.setObjectName("LoadFooter")
        row = QHBoxLayout(self)
        row.setContentsMargins(12, 6, 12, 6)
        title = QLabel("<b>Expected run load</b>")
        row.addWidget(title)
        self.text = QLabel()
        self.text.setObjectName("CardSub")
        self.text.setTextFormat(Qt.TextFormat.PlainText)
        self.text.setWordWrap(True)
        row.addWidget(self.text, 1)
        link = QLabel("<a href='plan'>Details</a>")
        link.setToolTip("Run ▸ Plan: the full prediction and what the run does")
        link.linkActivated.connect(lambda _href: self.details.emit())
        row.addWidget(link)
        self._timer = QTimer(self, singleShot=True, interval=DELAY_MS)
        self._timer.timeout.connect(self._update)

    def refresh(self) -> None:
        self._timer.start()

    def showEvent(self, event) -> None:       # noqa: N802
        super().showEvent(event)
        self._timer.start()

    def _update(self) -> None:
        if not self.isVisible():
            return
        try:
            pred = predict(self.model)
        except Exception:                    # a half-edited config must not break the tab
            self.text.setText("Not available for the current settings.")
            self.setToolTip("")
            return
        self.text.setText(pred.summary())
        tip = (["The study's simulation run (not the transducer calibration):"] + pred.lines()
               + ["Simulation time only, summed over the run's simulations; "
                  "segmentation, pseudo-CT, and queue wait not included."])
        self.setToolTip("<br>".join(html.escape(s) for s in tip))
