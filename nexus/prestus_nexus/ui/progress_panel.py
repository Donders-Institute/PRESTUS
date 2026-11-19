"""Progress of a run (progress.RunProgress), shown as

- ProgressPanel: a compact card with the stage, the percentage, a bar, and one line of
  detail (solver, time left, temperature), plus a "Live view" button once frames exist.
  Used by the Run page (local runs) and the jobs view (the selected cluster job).
- LiveViewPage: the Run section's "Live view" sub-tab, shown once a run has frames: the
  run's progress card and, filling the rest, the newest frame of each field that PRESTUS
  saved with simulation.live_view (LiveViewPanel: Pressure, Temperature tabs). The running
  stage's tab is selected and marked live; a finished stage keeps its last frame, marked
  final.
"""

from __future__ import annotations

import os
import time
from pathlib import Path
from typing import Callable

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QColor, QImage, QPalette, QPixmap
from PyQt6.QtWidgets import (QFrame, QHBoxLayout, QLabel, QProgressBar, QPushButton,
                             QSizePolicy, QStackedWidget, QVBoxLayout, QWidget)

from ..progress import RunProgress
from .tabbar import EmphasisTabBar

FIELD_NAMES = {"acoustic": "Pressure", "thermal": "Temperature"}

Frames = dict  # stage -> (local path, mtime)


def local_frames(p: RunProgress | None, to_local: Callable[[str], str | None] | None) -> tuple[Frames, dict]:
    """({stage: (local path, mtime)} of the frames visible here, {stage: path as written}
    of the frames that are not (no folder mapping, or not written yet))."""
    seen, missing = {}, {}
    for stage, path in (p.live.items() if p is not None else ()):
        local = to_local(path) if to_local else path
        try:
            seen[stage] = (local, os.path.getmtime(local))
        except (OSError, TypeError):
            missing[stage] = path
    return seen, missing


def blend_margin(img: QImage, colour: QColor, tol: int = 6) -> QImage:
    """The frame with its margin in colour: the pixels of the corner colour that connect to
    the image edge (the figure background of MATLAB or matplotlib) are replaced, so the
    frame sits on the GUI's background. Pixels of that colour inside the plot (enclosed by
    the axes) stay. A frame whose corners differ (no uniform margin) is returned as is."""
    import numpy as np
    img = img.convertToFormat(QImage.Format.Format_RGB32)
    w, h = img.width(), img.height()
    if w < 3 or h < 3:
        return img
    ptr = img.constBits()
    ptr.setsize(h * img.bytesPerLine())
    a = np.frombuffer(ptr, np.uint8).reshape(h, img.bytesPerLine() // 4, 4)[:, :w, :].copy()
    rgb = a[..., :3].astype(np.int16)
    bg = rgb[0, 0]
    corners = rgb[[0, 0, -1, -1], [0, -1, 0, -1]]
    if np.abs(corners - bg).max() > tol:
        return img
    match = np.abs(rgb - bg).max(axis=2) <= tol
    from scipy import ndimage                      # a Nexus dependency (pyproject.toml)
    labels, _n = ndimage.label(match)              # 4-connected regions of the margin colour
    edge = np.unique(np.concatenate([labels[0, :], labels[-1, :], labels[:, 0], labels[:, -1]]))
    filled = np.isin(labels, edge[edge > 0])
    a[filled, 0], a[filled, 1], a[filled, 2] = colour.blue(), colour.green(), colour.red()   # BGRA
    out = QImage(a.tobytes(), w, h, 4 * w, QImage.Format.Format_RGB32)
    return out.copy()


def live_stage(p: RunProgress | None) -> str:
    """The field the run is computing now ('acoustic', 'thermal', or '')."""
    if p is None or p.finished or p.failed or p.solver is None:
        return ""
    return p.solver.stage


class ProgressPanel(QFrame):
    live_requested = pyqtSignal()

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.setObjectName("Card")
        v = QVBoxLayout(self)
        v.setContentsMargins(12, 8, 12, 10)
        v.setSpacing(6)
        top = QHBoxLayout()
        self.title = QLabel()
        self.title.setObjectName("CardTitle")
        self.percent = QLabel()
        self.percent.setObjectName("CardTitle")
        self.b_live = QPushButton("Live view")
        self.b_live.setToolTip("Show the evolving pressure and temperature fields (Live view tab)")
        self.b_live.clicked.connect(self.live_requested)
        top.addWidget(self.title)
        top.addStretch(1)
        top.addWidget(self.percent)
        top.addWidget(self.b_live)
        v.addLayout(top)
        self.bar = QProgressBar()
        self.bar.setRange(0, 1000)
        self.bar.setTextVisible(False)
        self.bar.setFixedHeight(8)
        v.addWidget(self.bar)
        self.detail = QLabel()
        self.detail.setObjectName("CardSub")
        self.detail.setWordWrap(True)
        self.detail.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        v.addWidget(self.detail)
        self.whole = QLabel()                            # whole run, from the predicted stage times
        self.whole.setObjectName("CardSub")
        self.whole.setToolTip("The whole run, from the predicted time of each stage (Plan, 'Expected "
                              "run load') and the progress of the current one.")
        v.addWidget(self.whole)
        self.warning = QLabel()                          # e.g. 'will not finish before the time limit'
        self.warning.setWordWrap(True)
        v.addWidget(self.warning)
        self.frames: Frames = {}
        self.missing: dict = {}
        self.clear()

    def clear(self) -> None:
        self.frames, self.missing = {}, {}
        self.b_live.hide()
        self.warning.hide()
        self.setVisible(False)

    def set_warning(self, text: str) -> None:
        from . import style
        self.warning.setText(text)
        self.warning.setStyleSheet(f"color: {style.active()['red']};")
        self.warning.setVisible(bool(text))

    def show_progress(self, p: RunProgress | None, now: float,
                      to_local: Callable[[str], str | None] | None = None,
                      expected: dict | None = None) -> None:
        """Show p at time now; to_local maps a frame path, as the run wrote it, to this
        computer (None: not visible here); expected: predicted seconds per stage label
        (workload.expected_seconds) for the whole-run line."""
        self.overall = None
        if p is None or not (p.has_news() or p.live):
            self.clear()
            return
        self.setVisible(True)
        est = p.estimate(now)
        if est.fraction is None and not (p.finished or p.failed):
            self.bar.setRange(0, 0)                      # busy: no fraction known
        else:
            self.bar.setRange(0, 1000)
            self.bar.setValue(int(1000 * (est.fraction or 0.0)))
        self.title.setText(est.title + (f"  ·  {est.outer}" if est.outer else ""))
        self.percent.setText(est.percent)
        self.detail.setText(est.detail)
        self.detail.setVisible(bool(est.detail))
        self.overall = p.overall(expected, now)
        if self.overall is not None and not (p.finished or p.failed):
            from ..progress import format_duration
            frac, left = self.overall
            self.whole.setText(f"Whole run (predicted): ~{100 * frac:.0f} % · about {format_duration(left)} left")
        self.whole.setVisible(self.overall is not None and not (p.finished or p.failed))
        tip = ("Estimated from k-Wave's time estimate: k-Wave's MATLAB solvers print no "
               "progress while they run." if est.estimated else "")
        done = [s for s in p.done_stages if s != "prestus_pipeline"]
        if done:
            tip += ("\n" if tip else "") + "Done: " + ", ".join(done)
        for w in (self.title, self.percent, self.bar, self.detail):
            w.setToolTip(tip)
        self.frames, self.missing = local_frames(p, to_local)
        self.b_live.setVisible(bool(self.frames or self.missing))


class LiveViewPanel(QWidget):
    """The newest live-view frames of one run: field tabs, the image, a caption."""

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        v = QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 0)
        self.heading = QLabel()
        self.heading.setObjectName("CardTitle")
        v.addWidget(self.heading)
        self.tabs = EmphasisTabBar()                     # the Nexus's tab bar (dividers)
        self.tabs.setDrawBase(False)
        self.tabs.currentChanged.connect(lambda _i: self._show_current())
        self.tabs.tabBarClicked.connect(lambda _i: setattr(self, "_user_picked", True))
        v.addWidget(self.tabs)
        self.image = QLabel()
        self.image.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.image.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Ignored)
        self.stack = QStackedWidget()                   # figure frames (PNG) | data frames (drawn)
        self.stack.addWidget(self.image)
        self.field = None                               # FieldCanvas, made on the first data frame
        v.addWidget(self.stack, 1)
        self.caption = QLabel()
        self.caption.setObjectName("CardSub")
        self.caption.setWordWrap(True)
        self.caption.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        v.addWidget(self.caption)
        self._stages: list[str] = []
        self._frames: Frames = {}
        self._missing: dict = {}
        self._live = ""
        self._pix: dict[str, tuple[float, QPixmap]] = {}   # stage -> (mtime, pixmap)
        self._user_picked = False
        self._curve: list = []
        self._data: dict = {}                            # stage -> ((mtime, curve length), meta, plane)
        self._drawn = ""                                 # stage on the canvas

    # ── content ─────────────────────────────────────────────────────────
    def update_frames(self, frames: Frames, missing: dict, live: str, heading: str = "",
                      curve: list | None = None) -> None:
        """Show the frames; live is the stage being computed now ('' when none). The live
        stage's tab is selected unless the user picked a tab."""
        self._frames, self._missing, self._live = frames, missing, live
        self._curve = list(curve or [])
        self.heading.setText(heading)
        stages = [s for s in ("acoustic", "thermal") if s in frames or s in missing]
        stages += [s for s in list(frames) + list(missing) if s not in stages]
        if stages != self._stages:
            self.tabs.blockSignals(True)
            while self.tabs.count():
                self.tabs.removeTab(0)
            for s in stages:
                self.tabs.addTab(FIELD_NAMES.get(s, s.capitalize()))
            self._stages = stages
            self.tabs.blockSignals(False)
            self._user_picked = False
        for i, s in enumerate(stages):
            self.tabs.setTabText(i, f"{FIELD_NAMES.get(s, s.capitalize())} · {'live' if s == live else 'final'}")
        if live in stages and not self._user_picked:
            self.tabs.setCurrentIndex(stages.index(live))
        self._show_current()

    def current_stage(self) -> str:
        i = self.tabs.currentIndex()
        return self._stages[i] if 0 <= i < len(self._stages) else ""

    def _show_current(self) -> None:
        stage = self.current_stage()
        if not stage:
            self.image.clear()
            self.caption.setText("No frames yet.")
            return
        if stage not in self._frames:
            self.image.clear()
            self.caption.setText(f"The frame {self._missing.get(stage, '')} is not visible from this "
                                 "computer yet (for a cluster job: set the folder mapping under Submit over SSH).")
            return
        path, mtime = self._frames[stage]
        age = max(0, int(time.time() - mtime))
        state = ("updated every 5 s while the solver runs" if stage == self._live
                 else "last frame of this stage")
        if path.endswith(".json"):
            self._show_data(stage, path, mtime)
            meta = self._data.get(stage, (None, None, None))[1] or {}
            self.caption.setText(f"{meta.get('label') or FIELD_NAMES.get(stage, stage)} · {state} · "
                                 f"saved {age} s ago")
            return
        self.stack.setCurrentWidget(self.image)
        cached = self._pix.get(stage)
        if cached is None or cached[0] != mtime:
            img = QImage(path)
            if not img.isNull():                         # being replaced: keep the last one
                bg = self.palette().color(QPalette.ColorRole.Window)
                if bg.lightness() > 128:                 # light theme: margin in the GUI colour;
                    img = blend_margin(img, bg)          # dark theme: keep it as light 'paper'
                self._pix[stage] = (mtime, QPixmap.fromImage(img))
        self._scale()
        self.caption.setText(f"{FIELD_NAMES.get(stage, stage)} field · {state} · "
                             f"saved {age} s ago · {Path(path).name}")

    def _show_data(self, stage: str, path: str, mtime: float) -> None:
        """Draw a data frame (live_field) when it or the temperature curve changed."""
        from .live_field import FieldCanvas, load_frame
        if self.field is None:
            self.field = FieldCanvas()
            self.stack.addWidget(self.field)
        self.stack.setCurrentWidget(self.field)
        key = (mtime, len(self._curve) if stage == "thermal" else 0)
        cached = self._data.get(stage)
        if cached is not None and cached[0] == key and self._drawn == stage:
            return
        frame = load_frame(path) if cached is None or cached[0][0] != mtime else (cached[1], cached[2])
        if frame is None:
            if cached is None:
                return                                  # being written: the next tick
            frame = (cached[1], cached[2])
        meta, plane = frame
        self._data[stage] = (key, meta, plane)
        self.field.draw_frame(meta, plane, self._curve if stage == "thermal" else None)
        self._drawn = stage

    def _scale(self) -> None:
        cached = self._pix.get(self.current_stage())
        if cached is None:
            self.image.clear()
            return
        size = self.image.size()
        self.image.setPixmap(cached[1].scaled(max(1, size.width() - 16), max(1, size.height() - 16),
                                              Qt.AspectRatioMode.KeepAspectRatio,
                                              Qt.TransformationMode.SmoothTransformation))

    def resizeEvent(self, event) -> None:                # noqa: N802 (Qt)
        super().resizeEvent(event)
        self._scale()


class LiveViewPage(QWidget):
    """The Run section's Live view sub-tab: one run's progress card above its fields."""

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        v = QVBoxLayout(self)
        self.card = ProgressPanel()
        v.addWidget(self.card)
        self.view = LiveViewPanel()
        v.addWidget(self.view, 1)
        self.source = None                # 'local', or the job id shown

    def show_run(self, source, p: RunProgress | None, now: float, heading: str,
                 to_local: Callable[[str], str | None] | None = None) -> bool:
        """Show this run; returns whether it has frames."""
        self.source = source
        self.card.show_progress(p, now, to_local)
        self.card.b_live.hide()                       # this is the live view
        self.view.update_frames(self.card.frames, self.card.missing, live_stage(p), heading,
                                p.temperature if p is not None else None)
        return bool(self.card.frames or self.card.missing)
