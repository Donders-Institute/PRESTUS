"""Sonication protocol diagram (pulse / pulse train / train repetition) and TUS calculator exchange."""

from __future__ import annotations

import math
from pathlib import Path

from PyQt6.QtCore import QPointF, QRectF, Qt, QUrl, pyqtSignal
from PyQt6.QtGui import QColor, QDesktopServices, QFont, QFontInfo, QFontMetricsF, QPainter, QPainterPath, QPen
from PyQt6.QtWidgets import (QCheckBox, QDialog, QDialogButtonBox, QFileDialog, QHBoxLayout, QLabel,
                             QMessageBox, QPushButton, QSizePolicy, QToolButton, QVBoxLayout, QWidget)

from .. import timing as tm
from . import icons, style
from .export import PaintedExportable, attach_painted


def _fmt_s(v: float) -> str:
    if v <= 0 or not math.isfinite(v):
        return "–"
    if v < 1e-3:
        return f"{v * 1e6:.3g} µs"
    if v < 1:
        return f"{v * 1e3:.3g} ms"
    if v < 120:
        return f"{v:.3g} s"
    return f"{v / 60:.3g} min"


def _s(n: int) -> str:
    return "" if n == 1 else "s"


class _Handle:
    """A grabbable part of the diagram: an edge (drag) or a dimension label (click to type)."""

    def __init__(self, name: str, rect: QRectF, lane: int, x0: float, sx: float, edge: bool):
        self.name, self.rect, self.lane, self.x0, self.sx, self.edge = name, rect, lane, x0, sx, edge


class TimingDiagram(PaintedExportable, QWidget):
    """Three lanes: one pulse, one pulse train, the train repetitions.

    Edges can be dragged (pulse end, ramp, next pulse, train end, next train, repetition end)
    and the dimension labels clicked to type a value; `edited` carries {config key: seconds}."""

    FONT_DELTA = 0.0        # labels at the application font, as the other plots' text
    TITLES = ("Pulse", "Pulse train", "Repetitions", "Firing schedule")
    MAX_ROWS = 8            # transducers drawn in the async lane; more are listed as a count
    ROW_H = 16              # px per transducer row
    HEADROOM = 1.12         # lane span beyond its outermost value, so that edge can be dragged longer
    GRAB = 5                # px either side of an edge
    HEAD = 92               # px, side of the head sketch beside the async lane
    EDIT_TIP = "Drag an edge in the diagram, or click a value to type it." 
    HEAD_TIP = ("Seen from above, front at the top, R on the right: where each transducer sits "
                "and points (its focus), in its colour and number. Schematic: centred on the "
                "targets and scaled to the outermost transducer. Positions: as configured "
                "(voxels taken as RAS+-stored, as for SimNIBS T1 images).")

    edited = pyqtSignal(dict)
    export_name = "timing"

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self._c, self._w = style.active(), 0.0      # palette and width being painted (paint)
        self.p = tm.Protocol()
        self.sched: tm.AsyncSchedule | None = None
        self.heads: list = []                      # mt.schematic_top_view, beside the async lane
        self._head_rect: QRectF | None = None
        self._handles: list[_Handle] = []
        self._drag: _Handle | None = None
        self._drag_moved = False
        self._spans: dict[int, float] = {}        # lane -> span frozen while dragging
        self._editor = None
        self._layout_metrics()
        self.setMinimumHeight(self.LANE_H + 20)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.setMouseTracking(True)
        self.export_button = attach_painted(self)

    def _label_font(self) -> QFont:
        f = QFont(self.font())
        size = QFontInfo(f).pointSizeF()        # the font may be set in pixels (style sheet)
        f.setPointSizeF(max(9.0, size + self.FONT_DELTA))
        return f

    def _layout_metrics(self) -> None:
        """Line height, label column and lane height from the label font: every text box and
        lane grows with it, so a larger font never overlaps the next row."""
        f = self._label_font()
        self.lh = math.ceil(QFontMetricsF(f).height()) + 2
        f.setBold(True)
        self.label_w = max(QFontMetricsF(f).horizontalAdvance(t) for t in self.TITLES) + 20
        self.wave_h = 30                                   # pulse height above the baseline
        self.row1 = 8                                      # first arrow row below the baseline
        self.row2 = self.row1 + self.lh + 6                # second arrow row
        self.LANE_H = int(self.lh + 8 + self.wave_h + self.row2 + self.lh + 10)

    def set_protocol(self, p: tm.Protocol, sched: tm.AsyncSchedule | None = None,
                     heads: list | None = None) -> None:
        """``heads``: multitransducer.schematic_top_view of the transducers, drawn beside the
        async lane when any transducer has a position."""
        if self._drag is not None:              # a refresh mid-drag must not undo the drag
            return
        self.p = p
        self.sched = sched if (sched is not None and sched.n > 1) else None
        self.heads = list(heads or [])
        self._layout_metrics()
        lanes = 1 + int(p.has_train) + int(p.has_repetition)
        self.setFixedHeight(lanes * self.LANE_H + self._async_h() + 20)
        self.update()

    def _async_h(self) -> int:
        """Height of the async lane: one row per drawn transducer, and the caption."""
        if self.sched is None or not self.p.has_pulse:
            return 0
        rows = int(self.lh + 8 + min(self.sched.n, self.MAX_ROWS) * self.ROW_H + 6 + self.lh + 10)
        return max(rows, self.HEAD + 12) if self._has_heads() else rows

    def _has_heads(self) -> bool:
        return any(h is not None for h in self.heads)

    def _span(self, lane: int) -> float:
        """Time shown across a lane: its outermost value plus headroom (frozen during a drag)."""
        if lane in self._spans:
            return self._spans[lane]
        p = self.p
        if lane == 1:
            return (p.pri if p.has_train else p.pd) * self.HEADROOM
        if lane == 2:
            return (p.ptri if p.has_repetition else p.ptd) * self.HEADROOM
        return p.ptrd * self.HEADROOM

    # ── drawing helpers ─────────────────────────────────────────────────
    def _envelope(self, t: float) -> float:
        """Amplitude of one pulse at time t in [0, pd]."""
        p = self.p
        if t < 0 or t > p.pd:
            return 0.0
        r = min(p.ramp, p.pd / 2)
        if r <= 0:
            return 1.0
        x = min(t, p.pd - t) / r
        if x >= 1:
            return 1.0
        return 0.5 * (1 - math.cos(math.pi * x)) if p.ramp_shape == "tukey" else x

    def _arrow(self, qp: QPainter, x0: float, x1: float, y: float, text: str) -> QRectF:
        qp.setPen(QPen(QColor(self._c["ink3"]), 1))
        qp.drawLine(QPointF(x0, y), QPointF(x1, y))
        for x, d in ((x0, 1), (x1, -1)):
            qp.drawLine(QPointF(x, y), QPointF(x + 5 * d, y - 3))
            qp.drawLine(QPointF(x, y), QPointF(x + 5 * d, y + 3))
        qp.setPen(QColor(self._c["ink2"]))
        fm = QFontMetricsF(qp.font())
        w = max(fm.horizontalAdvance(text) + 8, x1 - x0)
        box = QRectF((x0 + x1) / 2 - w / 2, y + 2, w, self.lh)
        qp.drawText(box, Qt.AlignmentFlag.AlignHCenter, text)
        tw = fm.horizontalAdvance(text) + 8
        return QRectF((x0 + x1) / 2 - tw / 2, y + 2, tw, self.lh)      # the text itself: click target

    def _lane_frame(self, qp: QPainter, top: float, title: str) -> tuple[float, float, float, float]:
        qp.setPen(QColor(self._c["ink2"]))
        f = qp.font()
        f.setBold(True)
        qp.setFont(f)
        qp.drawText(QRectF(8, top + 4, self.label_w, self.lh), Qt.AlignmentFlag.AlignLeft, title)
        f.setBold(False)
        qp.setFont(f)
        x0, x1 = 8.0 + self.label_w, self._w - 16.0
        y_top = top + self.lh + 8
        y_base = y_top + self.wave_h
        qp.setPen(QPen(QColor(self._c["line2"]), 1))
        qp.drawLine(QPointF(x0, y_base), QPointF(x1, y_base))
        return x0, x1, y_base, y_top

    def _next_marker(self, qp: QPainter, x: float, yt: float, yb: float) -> None:
        """Dashed line where the next pulse (or train) starts: the repetition interval's edge."""
        pen = QPen(QColor(self._c["ink3"]), 1, Qt.PenStyle.DashLine)
        qp.setPen(pen)
        qp.drawLine(QPointF(x, yt), QPointF(x, yb))

    def _edge(self, name: str, x: float, yt: float, yb: float, lane: int, x0: float, sx: float) -> None:
        g = self.GRAB
        self._handles.append(_Handle(name, QRectF(x - g, yt - 4, 2 * g, yb - yt + 8), lane, x0, sx, True))

    def _label(self, name: str, rect: QRectF, lane: int, x0: float, sx: float) -> None:
        self._handles.append(_Handle(name, rect, lane, x0, sx, False))

    def paint(self, qp: QPainter, rect: QRectF, c: dict) -> None:
        """The diagram in palette ``c`` (style.active() on screen, export.py for files)."""
        self._c, self._w = c, rect.width()
        qp.fillRect(rect, QColor(c["surface"]))
        qp.setFont(self._label_font())
        self._handles = []
        p = self.p
        if self.export_button is not None:
            self.export_button.set_available(p.has_pulse)
        if not p.has_pulse:
            qp.setPen(QColor(c["ink3"]))
            qp.drawText(rect, Qt.AlignmentFlag.AlignCenter,
                        "Set the pulse duration (and train settings) to see the protocol.")
            return
        top = 6.0
        acc = QColor(self._c["accent"])
        soft = QColor(self._c["accent_soft"])

        # Lane 1: one pulse inside one PRI
        x0, x1, yb, yt = self._lane_frame(qp, top, "Pulse")
        sx = (x1 - x0) / self._span(1)
        path = QPainterPath(QPointF(x0, yb))
        n = 200
        for i in range(n + 1):
            t = p.pd * i / n
            path.lineTo(QPointF(x0 + t * sx, yb - (yb - yt) * self._envelope(t)))
        path.lineTo(QPointF(x0 + p.pd * sx, yb))
        qp.fillPath(path, soft)
        qp.setPen(QPen(acc, 1.6))
        qp.drawPath(path)
        self._edge("pd", x0 + p.pd * sx, yt, yb, 1, x0, sx)
        g = self.GRAB        # ramp: the corner where the rise ends (at the pulse start when there is none)
        xr = x0 + min(p.ramp, p.pd / 2) * sx
        self._handles.insert(0, _Handle("ramp", QRectF(xr - g, yt - g, 2 * g, 2 * g), 1, x0, sx, True))
        self._label("pd", self._arrow(qp, x0, x0 + p.pd * sx, yb + self.row1, f"PD {_fmt_s(p.pd)}"), 1, x0, sx)
        if p.ramp > 0:
            qp.setPen(QColor(self._c["ink3"]))
            text = f"{p.ramp_shape} ramp {_fmt_s(p.ramp)} (Ipa = {p.ramp_factor:.3f} × Itp)"
            box = QRectF(x0 + p.pd * sx + 8, yt, x1 - x0, self.lh)
            qp.drawText(box, Qt.AlignmentFlag.AlignLeft, text)
            box.setWidth(QFontMetricsF(qp.font()).horizontalAdvance(text))
            self._label("ramp", box, 1, x0, sx)
        if p.has_train:
            xp = x0 + p.pri * sx
            self._next_marker(qp, xp, yt, yb)
            self._edge("pri", xp, yt, yb, 1, x0, sx)
            self._label("pri", self._arrow(qp, x0, xp, yb + self.row2,
                                           f"PRI {_fmt_s(p.pri)}  ·  PRF {p.prf_hz:.4g} Hz  ·  "
                                           f"duty cycle {100 * p.duty_cycle_pulse_train:.3g} %"), 1, x0, sx)
        top += self.LANE_H

        # Lane 2: pulses in one train
        if p.has_train:
            x0, x1, yb, yt = self._lane_frame(qp, top, "Pulse train")
            sx = (x1 - x0) / self._span(2)
            npulse = p.pulses_per_train
            if npulse * 3 <= (p.ptd * sx):
                qp.setPen(Qt.PenStyle.NoPen)
                qp.setBrush(acc)
                for k in range(npulse):
                    xs = x0 + k * p.pri * sx
                    qp.drawRect(QRectF(xs, yt, max(1.0, p.pd * sx), yb - yt))
            else:   # too many pulses to draw: a block at the duty-cycle height
                qp.fillRect(QRectF(x0, yt, p.ptd * sx, yb - yt), soft)
                qp.setPen(QPen(acc, 1))
                qp.drawRect(QRectF(x0, yt, p.ptd * sx, yb - yt))
                qp.setPen(QColor(self._c["ink2"]))
                qp.drawText(QRectF(x0, yt, p.ptd * sx, yb - yt), Qt.AlignmentFlag.AlignCenter,
                            f"{npulse} pulse{_s(npulse)}")
            self._edge("ptd", x0 + p.ptd * sx, yt, yb, 2, x0, sx)
            self._label("ptd", self._arrow(qp, x0, x0 + p.ptd * sx, yb + self.row1,
                                           f"PTD {_fmt_s(p.ptd)}  ·  {npulse} pulse{_s(npulse)}"), 2, x0, sx)
            if p.has_repetition:
                xn = x0 + p.ptri * sx
                self._next_marker(qp, xn, yt, yb)
                self._edge("ptri", xn, yt, yb, 2, x0, sx)
                self._label("ptri", self._arrow(qp, x0, xn, yb + self.row2,
                                                f"PTRI {_fmt_s(p.ptri)}  ·  train duty cycle "
                                                f"{100 * p.duty_cycle_repetition:.3g} %"), 2, x0, sx)
            top += self.LANE_H

        # Lane 3: repetitions
        if p.has_repetition:
            x0, x1, yb, yt = self._lane_frame(qp, top, "Repetitions")
            sx = (x1 - x0) / self._span(3)
            # Each train as the pulse-train block above (soft fill, accent edge): the same
            # quantity at a coarser scale, not a warning.
            for k in range(p.trains):
                r = QRectF(x0 + k * p.ptri * sx, yt, max(1.0, p.ptd * sx), yb - yt)
                qp.fillRect(r, QColor(self._c["accent_soft"]))
                qp.setPen(QPen(QColor(self._c["accent"]), 1))
                qp.drawRect(r)
            xe = x0 + p.ptrd * sx
            self._next_marker(qp, xe, yt, yb)
            self._edge("ptrd", xe, yt, yb, 3, x0, sx)
            self._label("ptrd", self._arrow(qp, x0, xe, yb + self.row1,
                                            f"PTRD {_fmt_s(p.ptrd)}  ·  {p.trains} train{_s(p.trains)}"), 3, x0, sx)
            top += self.LANE_H

        # Lane 4 (async coupling): which transducer fires when, over one cycle
        if self.sched is not None:
            self._paint_async(qp, top)
        if self._drag is not None:          # the dragged edge, highlighted
            h = next((k for k in self._handles if k.edge and k.name == self._drag.name), None)
            if h is not None:
                qp.setPen(QPen(acc, 2))
                xc = h.rect.center().x()
                qp.drawLine(QPointF(xc, h.rect.top()), QPointF(xc, h.rect.bottom()))

    def _paint_async(self, qp: QPainter, top: float) -> None:
        """One row of on-blocks per transducer (up to MAX_ROWS) over one cycle of the schedule:
        two PRIs ('pulses'), one PTRI ('trains') or until every transducer started ('custom')."""
        s, c = self.sched, self._c
        qp.setPen(QColor(c["ink2"]))
        f = qp.font()
        f.setBold(True)
        qp.setFont(f)
        qp.drawText(QRectF(8, top + 4, self.label_w, self.lh), Qt.AlignmentFlag.AlignLeft, "Firing schedule")
        f.setBold(False)
        qp.setFont(f)
        x0, x1 = 8.0 + self.label_w, self._w - 16.0
        self._head_rect = None
        if self._has_heads():
            self._head_rect = QRectF(x1 - self.HEAD, top + 4, self.HEAD, self.HEAD)
            self._paint_heads(qp, self._head_rect)
            x1 -= self.HEAD + 16
        t_win = s.cycle()
        if not t_win > 0:
            return
        sx = (x1 - x0) / t_win
        y = top + self.lh + 8
        rows = min(s.n, self.MAX_ROWS)
        levels = set()
        for i in range(rows):
            yr = y + i * self.ROW_H
            qp.setPen(QColor(c["ink3"]))
            qp.drawText(QRectF(x0 - 34, yr, 30, self.ROW_H),
                        Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter, f"T{i + 1}")
            qp.setPen(QPen(QColor(c["line2"]), 1))
            qp.drawLine(QPointF(x0, yr + self.ROW_H - 2), QPointF(x1, yr + self.ROW_H - 2))
            blocks, level = s.blocks(i, t_win)
            levels.add(level)
            col = QColor(style.series(i))
            for a, b in blocks:             # pulses, or whole trains: solid in the transducer's colour
                r = QRectF(x0 + a * sx, yr + 2, max(1.5, (b - a) * sx), self.ROW_H - 5)
                qp.fillRect(r, col)
        what = {"pulses": "pulses take turns within each PRI",
                "trains": "pulse trains take turns within each PTRI",
                "custom": "custom offsets"}[s.mode]

        def off(o: float) -> str:
            return "0" if o <= 0 else f"+{_fmt_s(o)}"
        starts = " · ".join(f"T{i + 1} {off(src.offset)}" for i, src in enumerate(s.sources[:rows]))
        more = f" · {s.n - rows} more not shown" if s.n > rows else ""
        shown = "pulse trains shown as blocks" if levels == {"trains"} else ""
        text = "  ·  ".join(t for t in (what, starts + more, shown) if t)
        caption = QFontMetricsF(qp.font()).elidedText(f"{_fmt_s(t_win)}:  {text}",
                                                      Qt.TextElideMode.ElideRight, x1 - x0)
        self._arrow(qp, x0, x1, y + rows * self.ROW_H + 6, caption)

    def _paint_heads(self, qp: QPainter, box: QRectF) -> None:
        from .head_sketch import paint_heads
        paint_heads(qp, box, self.heads, self._c)

    # ── editing ─────────────────────────────────────────────────────────
    def _hit(self, pos) -> _Handle | None:
        for h in self._handles:             # ramp first (inserted at 0), then edges and labels in draw order
            if h.rect.contains(QPointF(pos)):
                return h
        return None

    def _tip(self, h: _Handle) -> str:
        what = tm.FIELD_NAMES[h.name]
        how = "Drag to change the" if h.edge else "Click to type the"
        return f"{how} {what} ({tm.FIELD_KEYS[h.name]})" + (" · click to type a value" if h.edge else "")

    def mouseMoveEvent(self, event) -> None:  # noqa: N802 (Qt API)
        if self._drag is not None:
            h = self._drag
            if abs(event.position().x() - self._press_x) >= 2:
                self._drag_moved = True
            if not self._drag_moved:
                return
            value = tm.dragged(self._drag_base, h.name, (event.position().x() - h.x0) / h.sx)
            p = tm.Protocol(**{**self._drag_base.__dict__, h.name: value})
            self.p = p
            self.setToolTip(f"{tm.FIELD_NAMES[h.name]}: {_fmt_s(value)}" if value > 0 else "")
            self.update()
            return
        h = self._hit(event.position())
        self.setCursor(Qt.CursorShape.SizeHorCursor if h is not None and h.edge
                       else Qt.CursorShape.PointingHandCursor if h is not None else Qt.CursorShape.ArrowCursor)
        on_head = self._head_rect is not None and self._head_rect.contains(event.position())
        self.setToolTip(self._tip(h) if h is not None else self.HEAD_TIP if on_head else
                        self.EDIT_TIP if self.p.has_pulse else "")

    def mousePressEvent(self, event) -> None:  # noqa: N802 (Qt API)
        h = self._hit(event.position()) if event.button() == Qt.MouseButton.LeftButton else None
        if h is None:
            return super().mousePressEvent(event)
        self._drag, self._drag_moved = h, False
        self._press_x = event.position().x()
        self._drag_base = tm.Protocol(**self.p.__dict__)
        self._spans = {lane: self._span(lane) for lane in (1, 2, 3)}

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802 (Qt API)
        h, moved = self._drag, self._drag_moved
        self._drag, self._spans = None, {}
        if h is None:
            return super().mouseReleaseEvent(event)
        if moved:
            value = getattr(self.p, h.name)
            if value != getattr(self._drag_base, h.name):
                self.edited.emit({tm.FIELD_KEYS[h.name]: value})
            self.update()
        else:
            self.edit_value(h.name, h.rect)

    def edit_value(self, name: str, near: QRectF | None = None) -> None:
        """A small text field over the plot for `name` (e.g. '20 ms'; the PRI also takes 'Hz')."""
        from PyQt6.QtWidgets import QLineEdit
        if self._editor is not None:
            self._editor.deleteLater()
        ed = QLineEdit(self)
        ed.setText(_fmt_s(getattr(self.p, name)).replace("–", ""))
        ed.setToolTip(f"{tm.FIELD_NAMES[name]}: a duration with unit (s, ms, µs, min)"
                      + ("; or a frequency (Hz)" if name == "pri" else "") + ". Enter applies, Esc cancels.")
        w = max(110, int(QFontMetricsF(ed.font()).horizontalAdvance("0000.00 ms") + 24))
        r = near or QRectF(self.label_w, 4, w, self.lh)
        x = int(min(max(r.center().x() - w / 2, 0), self.width() - w))
        ed.setGeometry(x, int(r.top()), w, ed.sizeHint().height())
        ed.selectAll()
        done = {"v": False}

        def commit():
            if done["v"]:
                return
            done["v"] = True
            value = tm.parse_duration(ed.text(), hz=(name == "pri"))
            ed.deleteLater()
            self._editor = None
            if value is not None and value != getattr(self.p, name) and (value > 0 or name == "ramp"):
                self.edited.emit({tm.FIELD_KEYS[name]: value})

        def cancel():
            done["v"] = True
            ed.deleteLater()
            self._editor = None

        ed.returnPressed.connect(commit)
        ed.editingFinished.connect(commit)
        from PyQt6.QtGui import QKeySequence, QShortcut
        QShortcut(QKeySequence(Qt.Key.Key_Escape), ed, activated=cancel)
        self._editor = ed
        ed.show()
        ed.setFocus()


class TimingPanel(QWidget):
    """Diagram, derived numbers and the TUS calculator buttons, above the timing fields."""

    changed = pyqtSignal()
    fold_toggled = pyqtSignal(bool)        # show (True) or hide the fields the diagram edits
    FOLDED = tuple(tm.FIELD_KEYS.values())

    def __init__(self, model, parent: QWidget | None = None):
        super().__init__(parent)
        self.model = model
        v = QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 8)
        self.diagram = TimingDiagram()
        from .export import subject_stem
        self.diagram.export_button.filename_hint = subject_stem(model, "timing")
        self.diagram.edited.connect(self._apply)
        v.addWidget(self.diagram)
        hint_row = QHBoxLayout()
        hint_row.addStretch(1)                  # how to edit: the diagram's hover text
        self.fold_button = QToolButton()
        self.fold_button.setCheckable(True)
        self.fold_button.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.fold_button.setAutoRaise(True)
        self.fold_button.toggled.connect(self._fold_toggled)
        self._fold_toggled(False)
        hint_row.addWidget(self.fold_button)
        v.addLayout(hint_row)
        self.summary = QLabel()
        self.summary.setObjectName("CardSub")
        self.summary.setWordWrap(True)
        self.summary.setTextFormat(Qt.TextFormat.RichText)
        v.addWidget(self.summary)
        row = QHBoxLayout()
        b_exp = QPushButton("Export for TUS calculator…")
        b_exp.setToolTip("Write this protocol (and frequency, target Isppa, aperture, focal depth) as a "
                         "settings file for the TUS calculator's “Upload settings”.")
        b_exp.clicked.connect(lambda: self.export_settings())
        b_imp = QPushButton("Import from TUS calculator…")
        b_imp.setToolTip("Read a settings file downloaded from the TUS calculator (JSON or YAML).")
        b_imp.clicked.connect(lambda: self.import_settings())
        b_open = QPushButton("Open TUS calculator ↗")
        b_open.setToolTip(tm.TUS_CALCULATOR_URL)
        b_open.clicked.connect(lambda: QDesktopServices.openUrl(QUrl(tm.TUS_CALCULATOR_URL)))
        for b in (b_exp, b_imp, b_open):
            row.addWidget(b)
        row.addStretch(1)
        v.addLayout(row)
        self.refresh()

    def _fold_toggled(self, shown: bool) -> None:
        self.fold_button.setText("Hide the values as fields" if shown else "Show the values as fields")
        self.fold_button.setIcon(icons.icon("collapse" if shown else "expand"))
        self.fold_toggled.emit(shown)

    def show_step(self, k: int) -> None:
        """Draw (and edit) step k of the sonication sequence (sequence card selection)."""
        self.step = max(0, int(k))
        self.refresh()

    def _seq_step(self):
        """(index, Step) shown for a sonication sequence, else None."""
        from .. import sequence as sq
        steps = sq.steps(self.model)
        if not steps:
            return None
        k = min(getattr(self, "step", 0), len(steps) - 1)
        return k, steps[k]

    def _apply(self, values: dict) -> None:
        """Values edited in the diagram ({config key: seconds}) into the config. In a sonication
        sequence the timing goes to the step shown, and the later steps move with it."""
        from .. import sequence as sq
        shown = self._seq_step()
        step_values = {k.split(".", 1)[1]: v for k, v in values.items()
                       if k.split(".", 1)[1] in tm.TIMING_OVERRIDES}
        if shown is not None and step_values:
            k, st = shown
            if "ptrd" in step_values:            # a step lasts whole pulse trains
                ptri = step_values.get("ptri", st.protocol.ptri)
                step_values["ptrd"] = sq.whole_trains(step_values["ptrd"], ptri)
            sq.set_timing(self.model, k, step_values)
            values = {k2: v for k2, v in values.items() if k2.split(".", 1)[1] not in step_values}
        for key, value in values.items():        # the pulse ramp: one for every step
            self.model.set(key, value)
        self.refresh()
        self.changed.emit()

    def protocol(self) -> tm.Protocol:
        shown = self._seq_step()
        base = tm.Protocol.from_model(self.model.get)
        if shown is None:
            return base
        p = shown[1].protocol
        return tm.Protocol(pd=p.pd, pri=p.pri, ptd=p.ptd, ptri=p.ptri, ptrd=p.ptrd,
                           ramp=base.ramp, ramp_shape=base.ramp_shape)

    def refresh(self) -> None:
        p = self.protocol()
        sched = tm.AsyncSchedule.from_model(self.model)
        heads = None
        if sched is not None:
            from .. import multitransducer as mt
            heads = mt.schematic_top_view(mt.model_transducers(self.model))
        self.diagram.set_protocol(p, sched, heads)
        parts = []
        shown = self._seq_step()
        if shown is not None:
            k, st = shown
            names = ", ".join(f"T{i + 1}" for i in st.transducers)
            parts.append(f"<b>Showing step {k + 1} ({names})</b>")   # the sequence table selects it
        if p.has_pulse and sched is None:
            parts.append(f"On-time {_fmt_s(p.sonication_time_s)} of {_fmt_s(p.total_duration_s)} "
                         f"(overall duty cycle {100 * p.overall_duty_cycle:.3g} %)")
        from .. import sequence as sq
        chained = sq.has_together_steps(self.model)      # one coherent run per step: no clash
        problems = p.problems() + (sched.problems(chained) if sched is not None and p.has_pulse else [])
        parts += [f"<span style='color:{style.active()['red']}'>{m}</span>" for m in problems]
        if sched is not None and p.has_pulse:
            parts += [f"<span style='color:{style.active()['amber']}'>{m}</span>" for m in sched.notes()]
        self.summary.setText(" · ".join(parts[:1]) + ("<br>" + "<br>".join(parts[1:]) if len(parts) > 1 else ""))

    # ── TUS calculator ──────────────────────────────────────────────────
    def _first(self, path: str):
        v = self.model.get(path)
        if isinstance(v, list):
            v = v[0] if v else None
        return v if isinstance(v, (int, float)) and not isinstance(v, bool) else None

    def _aperture_mm(self):
        od = self.model.get("transducer.annular.elem_od_mm")
        if isinstance(od, list) and od:
            return max(od)
        return self._first("transducer.matrix.outer_diameter_mm")

    def export_settings(self, path: str | Path | None = None) -> Path | None:
        settings = tm.to_tus_calculator(self.protocol(), self._first("transducer.freq_hz"),
                                        self._first("transducer.target_isppa_wcm2"),
                                        self._aperture_mm(), self._first("transducer.focal_distance_ep"))
        if path is None:
            start = str((self.model.path.parent if self.model.path else Path.home()) / "tus_calculator_settings.json")
            path, _ = QFileDialog.getSaveFileName(self, "TUS calculator settings", start, "JSON (*.json)")
            if not path:
                return None
        out = tm.write_settings(path, settings)
        w = self.window()
        if hasattr(w, "statusBar"):
            w.statusBar().showMessage(f"Wrote {out} — load it with “Upload settings” in the TUS calculator", 8000)
        return out

    def import_settings(self, path: str | Path | None = None, apply_optional: bool | None = None) -> bool:
        if path is None:
            start = str(self.model.path.parent if self.model.path else Path.home())
            path, _ = QFileDialog.getOpenFileName(self, "TUS calculator settings", start,
                                                  "Settings (*.json *.yml *.yaml)")
            if not path:
                return False
        res = tm.from_tus_calculator(tm.read_settings(path))
        if not res.timing:
            QMessageBox.warning(self, "TUS calculator", "\n".join(res.notes) or "Nothing to import.")
            return False
        if apply_optional is None:
            dlg = _ImportDialog(res, self)
            if dlg.exec() != QDialog.DialogCode.Accepted:
                return False
            apply_optional = dlg.optional.isChecked()
        for key, value in res.timing.items():
            self.model.set(key, value)
        if apply_optional:
            for key, value in res.optional.items():
                self.model.set(key, value)
        self.refresh()
        self.changed.emit()
        return True


class _ImportDialog(QDialog):
    def __init__(self, res: tm.Imported, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Import from TUS calculator")
        v = QVBoxLayout(self)
        lines = [f"<tt>{k}</tt> = {v:g}" if isinstance(v, float) else f"<tt>{k}</tt> = {v}"
                 for k, v in res.timing.items()]
        text = "<b>Protocol values to set:</b><br>" + "<br>".join(lines)
        if res.notes:
            text += "<br><br>" + "<br>".join(res.notes)
        lab = QLabel(text)
        lab.setWordWrap(True)
        v.addWidget(lab)
        self.optional = QCheckBox("Also set " + ", ".join(
            f"{k} = {val:g}" for k, val in res.optional.items()))
        self.optional.setVisible(bool(res.optional))
        v.addWidget(self.optional)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        v.addWidget(bb)
