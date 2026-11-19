"""Transducer cross-section: bowl and elements, exit plane, geometric focus, chosen focal depth
and the device's focal range."""

from __future__ import annotations

import math

from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import QColor, QPainter, QPainterPath, QPen
from PyQt6.QtWidgets import QSizePolicy, QWidget

from .export import PaintedExportable, attach_painted


def _num(v):
    if isinstance(v, list):
        v = v[0] if v else None
    return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else None


def _nums(v) -> list[float]:
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        return [float(v)]
    if isinstance(v, list):
        return [float(x) for x in v if isinstance(x, (int, float)) and not isinstance(x, bool)]
    return []


class TransducerSchematic(PaintedExportable, QWidget):
    export_name = "transducer"

    def __init__(self, model, parent: QWidget | None = None):
        super().__init__(parent)
        self.model = model
        self.ti = 0
        self.setMinimumSize(320, 220)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        from .export import subject_stem
        self.export_button = attach_painted(self, filename_hint=subject_stem(model, "transducer"))

    def set_transducer(self, ti: int) -> None:
        self.ti = ti
        self.update()

    def _geometry(self) -> dict | None:
        g = self.model.get
        ti = self.ti
        kind = g("transducer.type", ti)
        if kind == "matrix":
            a = (_num(g("transducer.matrix.outer_diameter_mm", ti)) or 0) / 2
            curved = bool(g("transducer.matrix.is_curved", ti))
            R = _num(g("transducer.matrix.curv_radius_mm", ti)) if curved else None
            dg = _num(g("transducer.matrix.dist_geom_ep_mm", ti))
            rings = [(0.0, a)] if a else []
        else:
            od = _nums(g("transducer.annular.elem_od_mm", ti))
            idd = _nums(g("transducer.annular.elem_id_mm", ti))
            rings = [(i / 2, o / 2) for i, o in zip(idd, od)] if len(idd) == len(od) else [(0.0, o / 2) for o in od[-1:]]
            a = max((o for _, o in rings), default=0.0)
            R = _num(g("transducer.annular.curv_radius_mm", ti))
            dg = _num(g("transducer.annular.dist_geom_ep_mm", ti))
        if not a:
            return None
        dev = self.model.device(ti)
        return {"a": a, "R": R, "dist_geom_ep": dg, "rings": rings,
                "fd": _num(g("transducer.focal_distance_ep", ti)),
                "min_foc": getattr(dev, "min_foc", None), "max_foc": getattr(dev, "max_foc", None),
                "name": (dev.label() if dev else g("transducer.name", ti)) or ""}

    PAD = 28

    def _fit(self, z_hi: float, z_max: float, a: float, W: float, H: float) -> tuple[float, float]:
        """(mm-per-pixel factor for both axes, so the drawing is to scale; x of the bowl apex)
        that fit the drawn extent (apex at z = 0 to the farthest mark z_hi), centred."""
        pad = self.PAD
        s = min((W - 2 * pad) / z_max, (H - 2 * pad - 18) / (2.4 * a))
        return s, max(pad, (W - z_hi * s) / 2)

    def _frame(self, z_hi, z_max, a, W, H, on_screen: bool) -> tuple[float, float]:
        """Scale and apex position for this paint: fitted to the drawing each time."""
        return self._fit(z_hi, z_max, a, W, H)

    def paint(self, qp: QPainter, rect: QRectF, C: dict) -> None:
        """The schematic in palette ``C`` (style.active() on screen, export.py for files)."""
        qp.fillRect(rect, QColor(C["surface"]))
        geo = self._geometry()
        if self.export_button is not None:
            self.export_button.set_available(geo is not None)
        if geo is None:
            qp.setPen(QColor(C["ink3"]))
            qp.drawText(rect, Qt.AlignmentFlag.AlignCenter,
                        "Choose a device or enter the element\ndiameters to see the transducer.")
            return
        a, R = geo["a"], geo["R"]

        def sag(r: float) -> float:          # axial depth of the bowl surface at radius r
            if not R or R <= a:
                return 0.0
            return R - math.sqrt(R * R - r * r)

        z_rim = sag(a)
        z_ep = R - geo["dist_geom_ep"] if R and geo["dist_geom_ep"] is not None else z_rim
        marks = [z_rim + 5, (R or 0), z_ep + (geo["fd"] or 0), z_ep + (geo["max_foc"] or 0)]
        z_hi = max(marks)
        z_max = z_hi * 1.1 + 5
        W, H = rect.width(), rect.height()
        s, x0 = self._frame(z_hi, z_max, a, W, H, qp.device() is self)
        X = lambda z: x0 + z * s  # noqa: E731
        yc = H / 2 + 6
        Y = lambda r: yc - r * s  # noqa: E731

        # axis
        qp.setPen(QPen(QColor(C["line2"]), 1, Qt.PenStyle.DashLine))
        qp.drawLine(QPointF(X(0) - 8, yc), QPointF(W - 8, yc))

        # focal range of the device (its label is drawn below the focus label, further down)
        box_top = box_bottom = yc
        steer = None
        if geo["min_foc"] is not None and geo["max_foc"] is not None:
            r = QRectF(QPointF(X(z_ep + geo["min_foc"]), Y(a * 0.35)), QPointF(X(z_ep + geo["max_foc"]), Y(-a * 0.35)))
            qp.fillRect(r, QColor(C["accent_soft"]))
            box_top, box_bottom = r.top(), r.bottom()
            steer = (r, f"steering range {geo['min_foc']:g}–{geo['max_foc']:g} mm")

        # bowl and elements (both sides)
        # alternate rings: the accent and a lighter tint of it, in either theme
        colours = [QColor(C["accent"]), QColor(C["accent"]).lighter(135)]
        for k, (ri, ro) in enumerate(geo["rings"]):
            pen = QPen(colours[k % 2], 5)
            pen.setCapStyle(Qt.PenCapStyle.FlatCap)
            qp.setPen(pen)
            for sign in (1, -1):
                path = QPainterPath()
                steps = 24
                for i in range(steps + 1):
                    r = ri + (ro - ri) * i / steps
                    pt = QPointF(X(sag(r)), Y(sign * r))
                    path.moveTo(pt) if i == 0 else path.lineTo(pt)
                qp.drawPath(path)

        # exit plane
        qp.setPen(QPen(QColor(C["ink3"]), 1, Qt.PenStyle.DotLine))
        qp.drawLine(QPointF(X(z_ep), Y(a * 1.1)), QPointF(X(z_ep), Y(-a * 1.1)))
        qp.drawText(QPointF(X(z_ep) + 3, Y(a * 1.1) + 10), "exit plane")

        # geometric focus and chosen focus
        if R:
            qp.setPen(QPen(QColor(C["ink3"]), 1))
            qp.drawEllipse(QPointF(X(R), yc), 3, 3)
            # above the steering box, clear of the focus marker and its label below
            qp.drawText(QPointF(X(R) - 20, min(yc - 14, box_top - 4)), "geometric")
        if geo["fd"]:
            zf = z_ep + geo["fd"]
            qp.setPen(QPen(QColor(C["thermal"]), 2))
            qp.drawLine(QPointF(X(zf) - 6, yc - 6), QPointF(X(zf) + 6, yc + 6))
            qp.drawLine(QPointF(X(zf) - 6, yc + 6), QPointF(X(zf) + 6, yc - 6))
            qp.setPen(QColor(C["thermal"]))
            focus_base = max(yc + 22, box_bottom + 14)     # below the steering box
            qp.drawText(QPointF(X(zf) - 30, focus_base), f"focus {geo['fd']:g} mm")
            label_top = focus_base + 6
        else:
            label_top = box_bottom + 2
        if steer is not None:
            r, label = steer
            qp.setPen(QColor(C["ink3"]))
            tw = qp.fontMetrics().horizontalAdvance(label) + 4
            qp.drawText(QRectF(min(r.left(), W - 8 - tw), label_top, tw, 14), Qt.AlignmentFlag.AlignLeft, label)

        qp.setPen(QColor(C["ink2"]))
        # On screen the device selector above names the device; an exported figure names it
        # itself (eliding the name, not the dimensions, when the header is too wide).
        dims = f"aperture {2 * a:.3g} mm" + (f"  ·  ROC {R:g} mm" if R else "")
        if qp.device() is not self:
            dims = "  ·  " + dims
            fm = qp.fontMetrics()
            dims = fm.elidedText(geo["name"], Qt.TextElideMode.ElideRight,
                                 max(0, int(W - 16 - fm.horizontalAdvance(dims)))) + dims
        qp.drawText(QRectF(8, 4, W - 16, 16), Qt.AlignmentFlag.AlignLeft, dims)
