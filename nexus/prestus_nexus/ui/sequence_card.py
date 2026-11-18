"""Sonication › Timing › Sonication sequence: stances fired one after another, with a pause
between them (sequence.py). One row per step; duration and pause are typed in the table, the
pulse timing in a small dialog, positions in the placement viewer. The strip shows the whole
protocol to scale. The card says how the sequence is simulated (one joint run or a chain of
runs) and why; the choice is a Nexus setting stored per study in the session."""

from __future__ import annotations

from typing import Callable

from PyQt6.QtCore import QPointF, QRectF, Qt, pyqtSignal
from PyQt6.QtGui import QColor, QPainter, QPen
from PyQt6.QtWidgets import (QAbstractItemView, QComboBox, QDialog, QDialogButtonBox, QFormLayout,
                             QFrame, QHBoxLayout, QHeaderView, QLabel, QLineEdit, QPushButton,
                             QTableWidget, QTableWidgetItem, QToolButton, QVBoxLayout, QWidget)

from .. import sequence as sq
from .. import timing as tm
from ..model import ConfigModel, is_empty
from . import style

COLUMNS = ("Step", "Stances", "Starts at", "Duration", "Pulse timing", "Pause after")
COL_DURATION, COL_PAUSE = 3, 5
INTRO = ("Sonicate from several stances one after another. Each step has its own position, "
         "duration, and pulse timing; the pause between steps is modelled as off time, and every "
         "step starts from the temperature and thermal dose the previous ones left.")
EMPTY_ONE = "One stance. To add a second stance that sonicates after this one:"
EMPTY_MANY = ("These transducers fire together or interleaved ({how}). To fire them one after "
              "another instead:")


def _timing_text(p: tm.Protocol) -> str:
    return (f"PD {sq.fmt_s(p.pd)} · PRI {sq.fmt_s(p.pri)} · PTD {sq.fmt_s(p.ptd)} · "
            f"PTRI {sq.fmt_s(p.ptri)}")


class SequenceStrip(QWidget):
    """The protocol to scale: one block per step, gaps for the pauses, then the steady-state."""

    H = 58

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.plan: sq.Plan | None = None
        self.selected = -1
        self.setMinimumHeight(self.H)
        self.setMaximumHeight(self.H)

    def set_plan(self, pl: sq.Plan | None, selected: int = -1) -> None:
        self.plan, self.selected = pl, selected
        self.update()

    def paintEvent(self, _ev) -> None:  # noqa: N802 (Qt API)
        pl = self.plan
        if pl is None or not pl.steps:
            return
        c = style.active()
        qp = QPainter(self)
        qp.setRenderHint(QPainter.RenderHint.Antialiasing)
        x0, x1 = 6.0, self.width() - 6.0
        total = max(pl.total_s, 1e-9)
        sx = (x1 - x0) / total
        top, h = 4.0, 26.0
        qp.setPen(QPen(QColor(c["line2"]), 1))
        qp.drawLine(QPointF(x0, top + h), QPointF(x1, top + h))
        for k, s in enumerate(pl.steps):
            r = QRectF(x0 + s.start_s * sx, top, max(2.0, s.duration_s * sx), h)
            # each transducer's colour, as in the async lanes and the head sketch; a step that
            # fires several transducers gets one band per transducer
            n = len(s.transducers)
            for j, ti in enumerate(s.transducers):
                band = QRectF(r.x(), r.y() + j * r.height() / n, r.width(), r.height() / n)
                qp.fillRect(band, QColor(style.series(ti)))
            if k == self.selected:
                qp.setPen(QPen(QColor(c["ink"]), 2))
                qp.drawRect(r)
            label = f"{k + 1}"
            if r.width() > 60:
                label += " · " + ", ".join(f"T{i + 1}" for i in s.transducers)
            qp.setPen(QColor("#ffffff"))
            qp.drawText(r, Qt.AlignmentFlag.AlignCenter, label)
        if pl.post_s > 0:
            end = max(s.end_s for s in pl.steps)
            r = QRectF(x0 + end * sx, top + h - 6, pl.post_s * sx, 6)
            soft = QColor(c["accent"])
            soft.setAlpha(60)
            qp.fillRect(r, soft)
        qp.setPen(QColor(c["ink3"]))
        y = top + h + 4
        qp.drawText(QRectF(x0, y, 80, 18), Qt.AlignmentFlag.AlignLeft, "0")
        what = f"{sq.fmt_s(pl.total_s)}" + (f" incl. {sq.fmt_s(pl.post_s)} steady-state" if pl.post_s > 0 else "")
        qp.drawText(QRectF(x1 - 260, y, 260, 18), Qt.AlignmentFlag.AlignRight, what)
        qp.end()


class TimingDialog(QDialog):
    """Pulse timing of one step (PD, PRI, PTD, PTRI), typed with units."""

    FIELDS = (("pd", "Pulse duration"), ("pri", "Pulse repetition interval (or PRF in Hz)"),
              ("ptd", "Pulse train duration"), ("ptri", "Pulse train repetition interval"))

    def __init__(self, step: int, p: tm.Protocol, parent: QWidget | None = None):
        super().__init__(parent)
        self.setWindowTitle(f"Pulse timing of step {step + 1}")
        form = QFormLayout(self)
        self.edits: dict[str, QLineEdit] = {}
        for f, label in self.FIELDS:
            e = QLineEdit(sq.fmt_s(getattr(p, f)))
            e.setToolTip("A duration with a unit: 20 ms, 0.5 s, 2 min")
            self.edits[f] = e
            form.addRow(label, e)
        self.error = QLabel()
        self.error.setStyleSheet(f"color: {style.active()['red']}")
        form.addRow(self.error)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        bb.accepted.connect(self._ok)
        bb.rejected.connect(self.reject)
        form.addRow(bb)
        self.values: dict[str, float] = {}

    def _ok(self) -> None:
        out = {}
        for f, e in self.edits.items():
            v = tm.parse_duration(e.text(), hz=(f == "pri"))
            if v is None or v <= 0:
                self.error.setText(f"Cannot read '{e.text()}'.")
                return
            out[f] = v
        self.values = out
        self.accept()


class SequenceCard(QFrame):
    """The card. ``session_getter()`` returns (session, study path) for the backend choice."""

    changed = pyqtSignal()
    edit_in_viewer = pyqtSignal(int)       # transducer index to place
    step_selected = pyqtSignal(int)        # step shown in the timing diagram

    def __init__(self, model: ConfigModel, session_getter: Callable[[], tuple] | None = None,
                 parent: QWidget | None = None):
        super().__init__(parent)
        self.model = model
        self.session_getter = session_getter or (lambda: (None, None))
        self._filling = False
        self.setObjectName("Card")
        v = QVBoxLayout(self)
        v.setContentsMargins(14, 12, 14, 12)
        v.setSpacing(8)
        head = QHBoxLayout()
        t = QLabel("Sonication sequence")
        t.setObjectName("CardTitle")
        t.setToolTip(INTRO)
        head.addWidget(t, 1)
        self.toggle = QToolButton()               # collapse / expand, as the other cards
        self.toggle.setAutoRaise(True)
        self.toggle.clicked.connect(self._toggle_clicked)
        head.addWidget(self.toggle)
        v.addLayout(head)
        self.content = QWidget()                  # everything below the title
        cv = QVBoxLayout(self.content)
        cv.setContentsMargins(0, 0, 0, 0)
        cv.setSpacing(8)
        v.addWidget(self.content)
        self._user_open: bool | None = None       # the user's choice, once they toggled

        # no sequence yet
        self.empty = QWidget()
        ev = QHBoxLayout(self.empty)
        ev.setContentsMargins(0, 0, 0, 0)
        self.empty_text = QLabel()
        self.empty_text.setWordWrap(True)
        ev.addWidget(self.empty_text, 1)
        self.start_btn = QPushButton()
        self.start_btn.clicked.connect(self._start)
        ev.addWidget(self.start_btn)
        cv.addWidget(self.empty)

        # the sequence
        self.body = QWidget()
        bv = QVBoxLayout(self.body)
        bv.setContentsMargins(0, 0, 0, 0)
        bv.setSpacing(8)
        self.strip = SequenceStrip()
        bv.addWidget(self.strip)
        self.table = QTableWidget(0, len(COLUMNS))
        self.table.setHorizontalHeaderLabels(COLUMNS)
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.DoubleClicked |
                                   QAbstractItemView.EditTrigger.EditKeyPressed |
                                   QAbstractItemView.EditTrigger.SelectedClicked)
        hh = self.table.horizontalHeader()
        hh.setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        hh.setStretchLastSection(True)
        self.table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.table.itemChanged.connect(self._item_changed)
        self.table.currentCellChanged.connect(lambda *_: self._buttons())
        self.table.currentCellChanged.connect(
            lambda row, *_: (row >= 0 and not self._filling) and self.step_selected.emit(row))
        self.table.cellDoubleClicked.connect(self._double_clicked)
        bv.addWidget(self.table)

        row = QHBoxLayout()
        self.add_btn = QPushButton("Add step")
        self.add_btn.setToolTip("A new stance that sonicates after the last step, with its timing")
        self.add_btn.clicked.connect(self._add)
        self.place_btn = QPushButton("Place in Viewer")
        self.place_btn.setToolTip("Show the selected step's stance in the placement viewer")
        self.place_btn.clicked.connect(self._place)
        self.timing_btn = QPushButton("Pulse timing…")
        self.timing_btn.setToolTip("PD, PRI, PTD, and PTRI of the selected step")
        self.timing_btn.clicked.connect(self._edit_timing)
        self.join_btn = QPushButton("Fire with previous step")
        self.join_btn.setToolTip("The selected step's transducers fire together with the previous "
                                 "step's: their fields are summed coherently, so they interfere")
        self.join_btn.clicked.connect(lambda: self._op(sq.join_previous))
        self.split_btn = QPushButton("Separate")
        self.split_btn.setToolTip("Give each transducer of the selected step its own step")
        self.split_btn.clicked.connect(lambda: self._op(sq.separate))
        self.remove_btn = QPushButton("Remove step")
        self.remove_btn.setToolTip("Remove the selected step and its transducers (not step 1's)")
        self.remove_btn.clicked.connect(lambda: self._op(sq.remove_step))
        for b in (self.add_btn, self.place_btn, self.timing_btn, self.join_btn, self.split_btn,
                  self.remove_btn):
            row.addWidget(b)
        row.addStretch(1)
        bv.addLayout(row)

        self.place_note = QLabel()
        self.place_note.setWordWrap(True)
        bv.addWidget(self.place_note)

        brow = QHBoxLayout()
        brow.addWidget(QLabel("Simulate as"))
        self.backend = QComboBox()
        for key in sq.BACKENDS:
            self.backend.addItem(sq.BACKEND_LABELS[key], key)
        self.backend.setToolTip("Automatic: one joint run, or a chain of runs when a step fires "
                                "several transducers together on per-transducer grids. Stored with "
                                "the Nexus for this study.")
        self.backend.currentIndexChanged.connect(self._backend_changed)
        brow.addWidget(self.backend)
        brow.addStretch(1)
        bv.addLayout(brow)
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.status.setTextFormat(Qt.TextFormat.RichText)
        bv.addWidget(self.status)
        cv.addWidget(self.body)
        self.refresh()

    # ── state ───────────────────────────────────────────────────────────
    def preference(self) -> str:
        session, study = self.session_getter()
        return sq.load_backend(session, study)

    def plan(self) -> sq.Plan | None:
        return sq.plan(self.model, self.preference())

    def _row(self) -> int:
        return self.table.currentRow()

    # ── collapse ────────────────────────────────────────────────────────
    def is_open(self) -> bool:
        """Closed until the user opens it (the toggle, or adding a later sonication)."""
        return bool(self._user_open)

    def _toggle_clicked(self) -> None:
        self._user_open = not self.is_open()
        self._show_open()

    def _show_open(self) -> None:
        shown = self.is_open()
        self.content.setVisible(shown)
        from . import icons
        self.toggle.setIcon(icons.icon("collapse" if shown else "expand"))
        self.toggle.setToolTip("Collapse this card" if shown else "Expand this card")

    # ── refresh ─────────────────────────────────────────────────────────
    def refresh(self) -> None:
        self._show_open()
        pl = self.plan()
        if pl is None:
            self.body.hide()
            self.empty.show()
            n = self.model.n_transducers
            if n <= 1:
                self.empty_text.setText(EMPTY_ONE)
                self.start_btn.setText("Add a later sonication")
            else:
                how = str(self.model.get("simulation.transducer_coupling") or "coherent")
                if how == "async":
                    how += f", '{self.model.get('timing.async_schedule') or 'pulses'}' schedule"
                self.empty_text.setText(EMPTY_MANY.format(how=how + " coupling"))
                self.start_btn.setText("Fire one after another")
            return
        self.empty.hide()
        self.body.show()
        row = self._row()
        self._filling = True
        self.table.setRowCount(pl.n)
        first = pl.steps[0].protocol
        for k, s in enumerate(pl.steps):
            gap = pl.pauses[k] if k < len(pl.pauses) else None
            cells = [str(k + 1), ", ".join(f"T{i + 1}" for i in s.transducers),
                     sq.fmt_s(s.start_s), sq.fmt_s(s.duration_s),
                     "as step 1" if k and all(abs(getattr(s.protocol, f) - getattr(first, f)) < 1e-12
                                              for f in ("pd", "pri", "ptd", "ptri"))
                     else _timing_text(s.protocol),
                     "–" if gap is None else sq.fmt_s(gap)]
            for col, text in enumerate(cells):
                it = QTableWidgetItem(text)
                editable = col == COL_DURATION or (col == COL_PAUSE and gap is not None)
                flags = Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled
                if editable:
                    flags |= Qt.ItemFlag.ItemIsEditable
                    it.setToolTip("Double-click to type a duration: 30 s, 2 min, 500 ms")
                it.setFlags(flags)
                self.table.setItem(k, col, it)
        self.table.setFixedHeight(self.table.horizontalHeader().height() +
                                  sum(self.table.rowHeight(r) for r in range(pl.n)) + 4)
        self._filling = False
        if 0 <= row < pl.n:
            self.table.selectRow(row)
        self.strip.set_plan(pl, self._row())
        self._set_backend(pl.preference)
        self._place_note()
        self._status(pl)
        self._buttons()

    def _set_backend(self, pref: str) -> None:
        session, study = self.session_getter()
        self.backend.setEnabled(session is not None and study is not None)
        if not self.backend.isEnabled():
            self.backend.setToolTip("Save the config first: the choice is stored per study")
        self.backend.blockSignals(True)
        self.backend.setCurrentIndex(max(0, self.backend.findData(pref)))
        self.backend.blockSignals(False)

    def _status(self, pl: sq.Plan) -> None:
        c = style.active()
        lines = []
        if pl.backend:
            what = {"joint": "One joint run", "chain": f"A chain of {pl.n} runs"}[pl.backend]
            lines.append(f"<b>{what}.</b> {pl.reason}.")
        colour = {"error": c["red"], "warning": c["amber"], "info": c["ink3"]}
        lines += [f"<span style='color:{colour.get(i.level, c['ink3'])}'>{i.message}.</span>"
                  for i in pl.issues]
        self.status.setText("<br>".join(lines))

    def _place_note(self) -> None:
        """Stances at transducer 1's position: a new step inherits it, and a step that
        returns to stance 1 keeps it (its acoustic field is then simulated once,
        acoustic_store_key.m)."""
        same = [ti + 1 for ti in range(1, self.model.n_transducers)
                if not self.model.is_overridden("transducer.trans_pos", ti)
                and not is_empty(self.model.get("transducer.trans_pos", 0))]
        if not same:
            self.place_note.hide()
            return
        names = ", ".join(f"T{i}" for i in same)
        self.place_note.setText(f"{names} {'has' if len(same) == 1 else 'have'} transducer 1's "
                                "position: select the step and place it in the Viewer, or keep it to "
                                "return to stance 1 (its acoustic field is then simulated once).")
        self.place_note.setStyleSheet(f"color: {style.active()['amber']}")
        self.place_note.show()

    def _buttons(self) -> None:
        pl = sq.plan(self.model)
        k = self._row()
        stps = pl.steps if pl else []
        has = 0 <= k < len(stps)
        self.place_btn.setEnabled(has)
        self.timing_btn.setEnabled(has)
        self.join_btn.setEnabled(has and k > 0)
        self.split_btn.setEnabled(has and stps[k].together)
        self.remove_btn.setEnabled(has and 0 not in stps[k].transducers)
        self.strip.set_plan(self.strip.plan, k)

    # ── actions ─────────────────────────────────────────────────────────
    def _done(self) -> None:
        self.refresh()
        self.changed.emit()

    def _start(self) -> None:
        self._user_open = True
        if self.model.n_transducers <= 1:
            ti = sq.add_step(self.model)
            self._done()
            self.table.selectRow(len(sq.steps(self.model)) - 1)
            self.edit_in_viewer.emit(ti)
            return
        sq.make_sequence(self.model)
        self._done()

    def _add(self) -> None:
        ti = sq.add_step(self.model)
        self._done()
        self.table.selectRow(self.table.rowCount() - 1)
        self._buttons()
        self.edit_in_viewer.emit(ti)

    def _place(self) -> None:
        stps = sq.steps(self.model)
        k = self._row()
        if 0 <= k < len(stps):
            self.edit_in_viewer.emit(stps[k].transducers[0])

    def _op(self, fn) -> None:
        k = self._row()
        if k < 0:
            return
        fn(self.model, k)
        self._done()

    def _edit_timing(self) -> None:
        stps = sq.steps(self.model)
        k = self._row()
        if not 0 <= k < len(stps):
            return
        dlg = TimingDialog(k, stps[k].protocol, self)
        if dlg.exec() == QDialog.DialogCode.Accepted:
            values = dict(dlg.values)
            values["ptrd"] = sq.whole_trains(stps[k].duration_s, values["ptri"])
            sq.set_timing(self.model, k, values)
            self._done()

    def _double_clicked(self, row: int, col: int) -> None:
        if col == 4:
            self.table.selectRow(row)
            self._edit_timing()

    def _item_changed(self, item: QTableWidgetItem) -> None:
        if self._filling:
            return
        k, col = item.row(), item.column()
        v = tm.parse_duration(item.text())
        stps = sq.steps(self.model)
        if v is None or v < 0 or not 0 <= k < len(stps):
            self.refresh()
            return
        if col == COL_DURATION and v > 0:
            ptri = stps[k].protocol.ptri
            sq.set_timing(self.model, k, {"ptrd": sq.whole_trains(v, ptri)})
        elif col == COL_PAUSE:
            sq.set_pause(self.model, k, v)
        self._done()

    def _backend_changed(self, _i: int) -> None:
        session, study = self.session_getter()
        if session is not None and study is not None:
            sq.save_backend(session, study, self.backend.currentData())
            if hasattr(session, "save"):
                session.save()
        self.refresh()
        self.changed.emit()
