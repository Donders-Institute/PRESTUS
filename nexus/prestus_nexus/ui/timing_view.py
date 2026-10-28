"""Sonication protocol diagram (pulse / pulse train / train repetition) and TUS calculator exchange."""

from __future__ import annotations

import math
from pathlib import Path

from PyQt6.QtCore import QPointF, QRectF, Qt, QUrl, pyqtSignal
from PyQt6.QtGui import QColor, QDesktopServices, QFont, QPainter, QPainterPath, QPen
from PyQt6.QtWidgets import (QCheckBox, QDialog, QDialogButtonBox, QFileDialog, QHBoxLayout, QLabel,
                             QMessageBox, QPushButton, QSizePolicy, QVBoxLayout, QWidget)

from .. import timing as tm
from . import style


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


class TimingDiagram(QWidget):
    """Three lanes: one pulse, one pulse train, the train repetitions."""

    LANE_H = 74

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.p = tm.Protocol()
        self.setMinimumHeight(self.LANE_H + 20)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

    def set_protocol(self, p: tm.Protocol) -> None:
        self.p = p
        lanes = 1 + int(p.has_train) + int(p.has_repetition)
        self.setFixedHeight(lanes * self.LANE_H + 20)
        self.update()

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

    def _arrow(self, qp: QPainter, x0: float, x1: float, y: float, text: str) -> None:
        qp.setPen(QPen(QColor(style.active()["ink3"]), 1))
        qp.drawLine(QPointF(x0, y), QPointF(x1, y))
        for x, d in ((x0, 1), (x1, -1)):
            qp.drawLine(QPointF(x, y), QPointF(x + 5 * d, y - 3))
            qp.drawLine(QPointF(x, y), QPointF(x + 5 * d, y + 3))
        qp.setPen(QColor(style.active()["ink2"]))
        w = max(60.0, x1 - x0)
        qp.drawText(QRectF((x0 + x1) / 2 - w / 2, y + 2, w, 16), Qt.AlignmentFlag.AlignHCenter, text)

    def _lane_frame(self, qp: QPainter, top: float, title: str) -> tuple[float, float, float, float]:
        qp.setPen(QColor(style.active()["ink2"]))
        f = qp.font()
        f.setBold(True)
        qp.setFont(f)
        qp.drawText(QRectF(8, top + 4, 120, 18), Qt.AlignmentFlag.AlignLeft, title)
        f.setBold(False)
        qp.setFont(f)
        x0, x1 = 128.0, self.width() - 16.0
        y_base, y_top = top + 40, top + 12
        qp.setPen(QPen(QColor(style.active()["line2"]), 1))
        qp.drawLine(QPointF(x0, y_base), QPointF(x1, y_base))
        return x0, x1, y_base, y_top

    def paintEvent(self, _event) -> None:  # noqa: N802 (Qt API)
        qp = QPainter(self)
        qp.setRenderHint(QPainter.RenderHint.Antialiasing)
        qp.fillRect(self.rect(), QColor(style.active()["surface"]))
        f = QFont(qp.font())
        f.setPointSizeF(max(8.0, f.pointSizeF() - 1))
        qp.setFont(f)
        p = self.p
        if not p.has_pulse:
            qp.setPen(QColor(style.active()["ink3"]))
            qp.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter,
                        "Set the pulse duration (and train settings) to see the protocol.")
            return
        top = 6.0
        acc = QColor(style.active()["accent"])
        soft = QColor(style.active()["accent_soft"])

        # Lane 1: one pulse inside one PRI
        x0, x1, yb, yt = self._lane_frame(qp, top, "Pulse")
        span = p.pri if p.has_train else p.pd * 1.25
        sx = (x1 - x0) / span
        path = QPainterPath(QPointF(x0, yb))
        n = 200
        for i in range(n + 1):
            t = p.pd * i / n
            path.lineTo(QPointF(x0 + t * sx, yb - (yb - yt) * self._envelope(t)))
        path.lineTo(QPointF(x0 + p.pd * sx, yb))
        qp.fillPath(path, soft)
        qp.setPen(QPen(acc, 1.6))
        qp.drawPath(path)
        self._arrow(qp, x0, x0 + p.pd * sx, yb + 8, f"PD {_fmt_s(p.pd)}")
        if p.ramp > 0:
            qp.setPen(QColor(style.active()["ink3"]))
            qp.drawText(QRectF(x0 + p.pd * sx + 8, yt + 2, 320, 14), Qt.AlignmentFlag.AlignLeft,
                        f"{p.ramp_shape} ramp {_fmt_s(p.ramp)} (Ipa = {p.ramp_factor:.3f} × Itp)")
        if p.has_train:
            self._arrow(qp, x0, x1, yb + 24, f"PRI {_fmt_s(p.pri)}  ·  PRF {p.prf_hz:.4g} Hz  ·  "
                        f"duty cycle {100 * p.duty_cycle_pulse_train:.3g} %")
        top += self.LANE_H

        # Lane 2: pulses in one train
        if p.has_train:
            x0, x1, yb, yt = self._lane_frame(qp, top, "Pulse train")
            span = p.ptri if p.has_repetition else p.ptd * 1.1
            sx = (x1 - x0) / span
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
                qp.setPen(QColor(style.active()["ink2"]))
                qp.drawText(QRectF(x0, yt, p.ptd * sx, yb - yt), Qt.AlignmentFlag.AlignCenter,
                            f"{npulse} pulses")
            self._arrow(qp, x0, x0 + p.ptd * sx, yb + 8, f"PTD {_fmt_s(p.ptd)}  ·  {npulse} pulses")
            if p.has_repetition:
                self._arrow(qp, x0, x1, yb + 24, f"PTRI {_fmt_s(p.ptri)}  ·  train duty cycle "
                            f"{100 * p.duty_cycle_repetition:.3g} %")
            top += self.LANE_H

        # Lane 3: repetitions
        if p.has_repetition:
            x0, x1, yb, yt = self._lane_frame(qp, top, "Repetitions")
            sx = (x1 - x0) / p.ptrd
            # Each train as the pulse-train block above (soft fill, accent edge): the same
            # quantity at a coarser scale, not a warning.
            for k in range(p.trains):
                r = QRectF(x0 + k * p.ptri * sx, yt, max(1.0, p.ptd * sx), yb - yt)
                qp.fillRect(r, QColor(style.active()["accent_soft"]))
                qp.setPen(QPen(QColor(style.active()["accent"]), 1))
                qp.drawRect(r)
            self._arrow(qp, x0, x1, yb + 8, f"PTRD {_fmt_s(p.ptrd)}  ·  {p.trains} trains")
        qp.end()


class TimingPanel(QWidget):
    """Diagram, derived numbers and the TUS calculator buttons, above the timing fields."""

    changed = pyqtSignal()

    def __init__(self, model, parent: QWidget | None = None):
        super().__init__(parent)
        self.model = model
        v = QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 8)
        self.diagram = TimingDiagram()
        v.addWidget(self.diagram)
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

    def protocol(self) -> tm.Protocol:
        return tm.Protocol.from_model(self.model.get)

    def refresh(self) -> None:
        p = self.protocol()
        self.diagram.set_protocol(p)
        parts = []
        if p.has_pulse:
            parts.append(f"On-time {_fmt_s(p.sonication_time_s)} of {_fmt_s(p.total_duration_s)} "
                         f"(overall duty cycle {100 * p.overall_duty_cycle:.3g} %)")
        parts += [f"<span style='color:{style.active()['red']}'>{m}</span>" for m in p.problems()]
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
