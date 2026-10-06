"""State of the shared MATLAB session (bridge.MatlabBridge), for pages that use it."""

from __future__ import annotations

from PyQt6.QtWidgets import QHBoxLayout, QLabel, QMessageBox, QToolButton, QWidget

from . import style

TEXT = {
    "off": "MATLAB: starts on first use",
    "starting": "MATLAB: starting… (first call 15–60 s)",
    "ready": "MATLAB: ready",
    "busy": "MATLAB: working…",
    "failed": "MATLAB: could not start",
}
COLOUR = {"off": "ink3", "starting": "amber", "ready": "green", "busy": "accent", "failed": "red"}


class MatlabStatus(QWidget):
    def __init__(self, bridge, parent: QWidget | None = None):
        super().__init__(parent)
        self.bridge = bridge
        h = QHBoxLayout(self)
        h.setContentsMargins(0, 0, 0, 0)
        self.label = QLabel()
        self.label.setObjectName("CardSub")
        h.addWidget(self.label)
        self.details = QToolButton(text="?")
        self.details.setToolTip("Why MATLAB could not start")
        self.details.clicked.connect(self._show_error)
        h.addWidget(self.details)
        self.stop_btn = QToolButton(text="Stop")
        self.stop_btn.setToolTip("End the MATLAB session now (it also ends after 30 min without use)")
        self.stop_btn.clicked.connect(bridge.stop)
        h.addWidget(self.stop_btn)
        bridge.state_changed.connect(self._update)
        self._update(bridge.state)

    def _update(self, state: str) -> None:
        colour = style.active()[COLOUR.get(state, "ink3")]
        self.label.setText(f"<span style='color:{colour}'>●</span> {TEXT.get(state, state)}")
        self.details.setVisible(state == "failed")
        self.stop_btn.setVisible(state in ("ready", "busy", "starting"))

    def _show_error(self) -> None:
        QMessageBox.warning(self, "MATLAB", self.bridge.error or "MATLAB could not start.")
