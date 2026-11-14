"""The segmentation step: its status and 'Run segmentation' (a segmentation-only run, as a
cluster job where the Run tab submits there). Shown on the Medium tab's Segmentation card, on
the Placement tab when the placement method needs a segmentation that is not there yet, and
(compact) wherever a page waits for one."""

from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QFrame, QHBoxLayout, QLabel, QMessageBox, QPushButton, QWidget

from .. import checks, launcher
from .. import segmentation as sg
from . import style
from .partial_run import RunFollower, run_partial

TASK = "Segmentation"


def run_segmentation(window, rerun: bool = False) -> str | None:
    """Start a segmentation-only run ('job', 'local' or None). A re-run over an existing m2m
    folder asks first: charm deletes it, and what was derived from it goes stale."""
    st = sg.status(window.model)
    if st.state in ("found", "incomplete"):
        answer = QMessageBox.warning(
            window, "Segment again",
            f"SimNIBS charm deletes {st.m2m} and segments the T1 (and T2) again, which takes an hour "
            "or more. Placements, head models and MNI conversions made from the current segmentation "
            "no longer match it.\n\nSegment again?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Cancel)
        if answer != QMessageBox.StandardButton.Yes:
            return None
        rerun = True
    return run_partial(window, TASK, "segmentation", "the SimNIBS segmentation",
                       launcher.segmentation_overrides(rerun))


class SegmentationPanel(QFrame):
    """Status line and Run / Segment again. ``when_needed``: shown only while something in the
    config needs a segmentation that is not complete (the Placement tab's banner)."""

    def __init__(self, window, when_needed: bool = False, parent: QWidget | None = None):
        super().__init__(parent)
        self.window = window
        self.when_needed = when_needed
        self.follower = RunFollower(window, TASK, self._run_done)
        self._note = ""
        h = QHBoxLayout(self)
        if when_needed:                         # a header row of its own: drawn as a card
            self.setObjectName("Card")
            h.setContentsMargins(14, 8, 14, 8)
        else:
            h.setContentsMargins(0, 0, 0, 0)
        self.label = QLabel()
        self.label.setWordWrap(True)
        self.label.setTextFormat(Qt.TextFormat.RichText)
        h.addWidget(self.label, 1)
        self.button = QPushButton("Run segmentation")
        self.button.clicked.connect(self.run)
        h.addWidget(self.button)
        self.refresh()

    def run(self) -> None:
        how = run_segmentation(self.window)
        if how is not None:
            self._note = ("Submitted as a job (Run ▸ Jobs); this updates when charm has written the "
                          "segmentation." if how == "job" else "Running in MATLAB here (Run tab).")
            self.follower.start(how)
            self.refresh()

    def _run_done(self, ok) -> None:
        checks.clear_status_cache()             # placement modes re-check the m2m folder
        if ok is not None:
            self._note = "" if ok else "Segmentation failed: see the Run tab's log."
        self.refresh()

    def refresh(self) -> None:
        m = self.window.model
        st = sg.status(m)
        need = sg.needed_by(m)
        if self.when_needed:
            self.setVisible(bool(need) and st.state != "found")
        c = style.active()
        if st.state == "found":
            text = f"<span style='color:{c['green']}'>✓</span> Segmentation found: {st.m2m}"
            self.button.setText("Segment again…")
            self.button.setToolTip("Run SimNIBS charm again (deletes the current segmentation first)")
        else:
            colour = c["red"] if need else c["amber"]
            uses = f" Needed by {', '.join(need)}." if need else ""
            text = f"<span style='color:{colour}'>✗</span> {st.message}{uses}"
            if not st.can_run:
                text += " Set the SimNIBS folder (startup.simnibs_bin_path) to create it."
            self.button.setText("Run segmentation")
            self.button.setToolTip("Run only SimNIBS charm for this subject (an hour or more): as a cluster "
                                   "job where the Run tab submits there, else in MATLAB here. Nothing else "
                                   "runs (modules.segmentation_only).")
        self.button.setEnabled(st.state != "unset" and st.can_run and not self.follower.waiting)
        if self._note:
            text += f"<br><span style='color:{c['ink3']}'>{self._note}</span>"
        self.label.setText(text)
