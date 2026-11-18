"""Placement ▸ TUS entry: run the TUS_entry placement (Atkinson et al.) and read its result.

TUS_entry (an R tool) chooses the scalp entry point from which a transducer disk sits most
perpendicular on the scalp, within a target-distance window. position_transducer_tus_entry.m
writes its result, tus_entry_report.csv (entry and target voxel, depth, angle score), to the
placement run folder. This page summarises the settings, runs only the placement (here, or
as a cluster job where the Run tab submits there) and shows the report; the Viewer shows the
transducer it placed.
"""

from __future__ import annotations

import csv
import io
from pathlib import Path

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QHBoxLayout, QLabel, QPushButton, QVBoxLayout, QWidget

from .. import checks
from ..model import is_empty
from ..outputs import dir_placement_run
from .empty_state import EmptyStack
from .medium_selector import setup_note
from .plantus_page import run_placement

TASK = "TUS_entry placement"
REPORT = "tus_entry_report.csv"
SETTINGS = (
    ("placement.tus_entry.script_path", "TUS_entry.R"),
    ("placement.tus_entry.env_path", "R environment"),
    ("placement.tus_entry.target_name", "Target name"),
    ("placement.tus_entry.target_mask", "Target mask"),
    ("placement.tus_entry.minimal_distance_cm", "Min. distance [cm]"),
    ("placement.tus_entry.maximal_distance_cm", "Max. distance [cm]"),
    ("placement.tus_entry.transducer_size_cm", "Transducer diameter [cm]"),
)


def report_path(model) -> Path | None:
    run = dir_placement_run(model, "tus_entry")
    return None if run is None else run / REPORT


def read_report(path: str | Path) -> dict | None:
    """The first row of tus_entry_report.csv as {column: float or text}; None if unreadable."""
    try:
        text = Path(path).read_text(encoding="utf-8")
    except OSError:
        return None
    rows = list(csv.DictReader(io.StringIO(text)))
    if not rows:
        return None
    out = {}
    for k, v in rows[0].items():
        try:
            out[k.strip()] = float(v)
        except (TypeError, ValueError):
            out[k.strip()] = v
    return out


def report_rows(rep: dict) -> list[tuple[str, str]]:
    """(label, value) rows of a TUS_entry report, voxel indices 1-based as TUS_entry writes."""
    def vox(prefix):
        try:
            return "[" + " ".join(f"{rep[f'{prefix}_{a}']:.0f}" for a in "xyz") + "]"
        except (KeyError, TypeError, ValueError):
            return "n/a"

    def num(key, fmt):
        v = rep.get(key)
        return fmt.format(v) if isinstance(v, float) else "n/a"

    return [("Scalp entry (voxel)", vox("scalp")),
            ("Target centre (voxel)", vox("target")),
            ("Entry–target distance", num("distance_mm", "{:.1f} mm")),
            ("Angle score", num("angle_score", "{:.3f}") + " (SD of disk-edge-to-scalp distances; lower is flatter)")]


class TusEntryPage(QWidget):
    def __init__(self, window, parent: QWidget | None = None):
        super().__init__(parent)
        self.window = window
        self._waiting = False
        self._watching = False
        self._stale = False

        v = QVBoxLayout(self)
        intro = QLabel("TUS_entry chooses the scalp entry point from which the transducer sits flattest on "
                       "the scalp, within the distance window to the target. Run it here, then check the "
                       "transducer it placed in the Viewer.")
        intro.setObjectName("CardSub")
        intro.setWordWrap(True)
        v.addWidget(intro)
        self.setup = QLabel()                   # what the method still needs (Placement method row)
        self.setup.setTextFormat(Qt.TextFormat.RichText)
        self.setup.setWordWrap(True)
        self.setup.hide()
        v.addWidget(self.setup)
        bar = QHBoxLayout()
        self.run_btn = QPushButton("Run TUS_entry")
        self.run_btn.setToolTip("Run only the TUS_entry placement (no simulation): in MATLAB here, or as a "
                                "cluster job where the Run tab submits there")
        self.run_btn.clicked.connect(self.run)
        bar.addWidget(self.run_btn)
        viewer = QPushButton("Open in Viewer")
        viewer.clicked.connect(self._viewer)
        bar.addWidget(viewer)
        bar.addStretch(1)
        v.addLayout(bar)

        self._content = QWidget()
        c = QVBoxLayout(self._content)
        c.setContentsMargins(0, 0, 0, 0)
        self.settings = QLabel()
        self.settings.setTextFormat(Qt.TextFormat.RichText)
        self.settings.setWordWrap(True)
        c.addWidget(self.settings)
        self.result = QLabel()
        self.result.setTextFormat(Qt.TextFormat.RichText)
        self.result.setWordWrap(True)
        c.addWidget(self.result)
        c.addStretch(1)
        self.stack = EmptyStack(self._content)
        v.addWidget(self.stack, 1)

    @property
    def model(self):
        return self.window.model

    def refresh(self) -> None:
        if self.window.isVisible() and not self.isVisible():
            self._stale = True
            return
        is_mode = self.model.get("placement.mode") == "tus_entry"
        note = setup_note(checks.placement_mode_status(self.model).get("tus_entry")) if is_mode else ""
        self.setup.setText(note)
        self.setup.setVisible(bool(note))
        self.run_btn.setEnabled(is_mode and not self._waiting)
        if not is_mode:
            self.stack.show_empty("TUS_entry is not the placement method",
                                  "Choose 'tus_entry' as the placement method (top of this tab) to run it.")
            return
        self.stack.show_content()
        rows = "".join(
            f"<tr><td>{label}</td><td style='padding-left:8px'>"
            f"{'—' if is_empty(self.model.get(key)) else self.model.get(key)}</td></tr>"
            for key, label in SETTINGS)
        self.settings.setText(f"<b>Settings</b> (Placement ▸ Settings)<table>{rows}</table>")
        path = report_path(self.model)
        rep = read_report(path) if path is not None else None
        if rep is None:
            where = f"<br>Looked for {path}" if path is not None else \
                "<br>The run folder needs path.sim (or path.placement) and the subject."
            self.result.setText(f"<b>Result</b><br>No TUS_entry result yet: click <i>Run TUS_entry</i>.{where}")
            return
        rows = "".join(f"<tr><td>{k}</td><td style='padding-left:8px'><b>{val}</b></td></tr>"
                       for k, val in report_rows(rep))
        self.result.setText(f"<b>Result</b> ({path.parent})<table>{rows}</table>")

    def showEvent(self, event) -> None:        # noqa: N802
        super().showEvent(event)
        if self._stale:
            self._stale = False
            self.refresh()

    def run(self) -> None:
        from .main_window import RUN
        how = run_placement(self.window, "tus_entry", TASK)
        run = self.window.pages[RUN]
        if how == "local" and not self._waiting:
            self._waiting = True
            run.finished.connect(self._done)
        elif how == "job" and not self._watching:
            self._watching = True
            run.jobs_changed.connect(self._jobs_changed)

    def _jobs_changed(self) -> None:
        if self.isVisible():
            self.refresh()

    def _done(self, task: str, ok: bool) -> None:
        from .main_window import RUN
        if task != TASK:
            return
        self.window.pages[RUN].finished.disconnect(self._done)
        self._waiting = False
        self.window.statusBar().showMessage(
            "TUS_entry finished: the Viewer shows the placement." if ok else
            "TUS_entry failed: see the Run tab's log.", 10000)
        self.refresh()

    def _viewer(self) -> None:
        from .main_window import VIEWER
        self.window.show_special(VIEWER)
