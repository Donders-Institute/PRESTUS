"""A small head seen from above with each transducer (numbered dot in its series colour) and
its beam to the focus (multitransducer.schematic_top_view). Drawn beside the async lane of
the timing diagram and in the transducers overview on Placement ▸ Settings."""

from __future__ import annotations

from PyQt6.QtCore import QPointF, QRectF, QSize, Qt
from PyQt6.QtGui import QColor, QFont, QPainter, QPainterPath, QPen
from PyQt6.QtWidgets import QSizePolicy, QWidget

from ..multitransducer import HEAD_RX, HEAD_RY
from . import style

TIP = ("Seen from above, front at the top, R on the right: where each transducer sits and "
       "points (its focus), in its colour and number. Schematic: centred on the targets and "
       "scaled to the outermost transducer.")


def paint_heads(qp: QPainter, box: QRectF, heads: list, c: dict) -> None:
    """The sketch inside ``box`` in palette ``c``; ``heads`` from schematic_top_view."""
    qp.save()
    qp.setRenderHint(QPainter.RenderHint.Antialiasing)
    r = (min(box.height(), box.width() / HEAD_RX * HEAD_RY) / 2 - 8) / HEAD_RY   # px per unit
    cx, cy = box.center().x(), box.center().y() + 3
    pt = lambda v: QPointF(cx + v[0] * r, cy - v[1] * r)     # noqa: E731 (A up, R right)
    qp.setPen(QPen(QColor(c["line2"]), 1.5))
    qp.setBrush(QColor(c["surface2"]))
    top = cy - HEAD_RY * r
    nose = QPainterPath(QPointF(cx - 5, top + 1))
    nose.lineTo(QPointF(cx, top - 6))
    nose.lineTo(QPointF(cx + 5, top + 1))
    qp.drawPath(nose)
    qp.drawEllipse(QPointF(cx, cy), HEAD_RX * r, HEAD_RY * r)
    f = QFont(qp.font())
    f.setPointSizeF(max(6.5, f.pointSizeF() - 3))
    qp.setFont(f)
    qp.setPen(QColor(c["ink3"]))
    qp.drawText(QRectF(cx + HEAD_RX * r - 2, cy - 6, 12, 12), Qt.AlignmentFlag.AlignCenter, "R")
    qp.drawText(QRectF(cx - HEAD_RX * r - 10, cy - 6, 12, 12), Qt.AlignmentFlag.AlignCenter, "L")
    dot = 6.5 if r < 60 else 8.0
    for i, h in enumerate(heads):
        if h is None:
            continue
        t, fo = h
        col = QColor(style.series(i))
        if fo is not None:
            qp.setPen(QPen(col, 1.5))
            qp.drawLine(pt(t), pt(fo))
            qp.setBrush(col)
            qp.drawEllipse(pt(fo), 1.8, 1.8)
        qp.setPen(QPen(QColor(c["surface"]), 1))
        qp.setBrush(col)
        qp.drawEllipse(pt(t), dot, dot)
        qp.setPen(QColor("#ffffff"))
        qp.drawText(QRectF(pt(t).x() - 8, pt(t).y() - 8, 16, 16), Qt.AlignmentFlag.AlignCenter, str(i + 1))
    qp.restore()


def has_heads(heads: list) -> bool:
    return any(h is not None for h in heads)


class HeadSketch(QWidget):
    """The sketch as a widget of fixed size (hidden while no transducer has a position)."""

    def __init__(self, side: int = 150, parent: QWidget | None = None):
        super().__init__(parent)
        self.heads: list = []
        self.side = side
        self.setFixedSize(QSize(side, side))
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self.setToolTip(TIP)

    def set_heads(self, heads: list, how: str = "") -> None:
        self.heads = list(heads or [])
        self.setVisible(has_heads(self.heads))
        self.setToolTip(TIP + (f" Positions: {how}." if how else ""))
        self.update()

    def paintEvent(self, _ev) -> None:  # noqa: N802 (Qt API)
        if not has_heads(self.heads):
            return
        qp = QPainter(self)
        paint_heads(qp, QRectF(0, 0, self.width(), self.height()), self.heads, style.active())
        qp.end()
