"""Results ▸ Heating: a run's thermal time courses per tissue, from the HEAT file
(cache/<tag>_<medium>_HEAT<affix>.mat, see results.read_heat). Pick the metric and the
tissues to show; the choice of tissues is kept across metrics and runs."""

from __future__ import annotations

from pathlib import Path

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (QCheckBox, QComboBox, QHBoxLayout, QLabel, QListWidget,
                             QListWidgetItem, QPushButton, QVBoxLayout, QWidget)

from .. import results as rs
from .empty_state import EmptyStack
from . import style

try:
    from .mpl_theme import ThemedCanvas as FigureCanvasQTAgg   # drawn in the Nexus theme
    from matplotlib.figure import Figure
    HAVE_MPL = True
except ImportError:          # pragma: no cover
    HAVE_MPL = False

PALETTE = style.SERIES          # colour-blind-safe categorical set, one colour per tissue


class HeatPanel(QWidget):
    """Metric selector, tissue checklist and a time-course plot."""

    def __init__(self, parent: QWidget | None = None, filename_hint=None):
        super().__init__(parent)
        self.data: rs.HeatData | None = None
        self._key: tuple | None = None            # (path, mtime) of the loaded file
        self._hidden: set[str] = set()            # tissues unticked by the user (lower case)
        self._colours: dict[str, str] = {}        # tissue (lower case) -> colour, stable

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        content = QWidget()
        self.stack = EmptyStack(content)
        outer.addWidget(self.stack)
        v = QVBoxLayout(content)
        v.setContentsMargins(0, 0, 0, 0)
        tb = QHBoxLayout()
        tb.addWidget(QLabel("Metric"))
        self.metric = QComboBox()
        self.metric.setMinimumWidth(240)
        self.metric.currentTextChanged.connect(self._metric_changed)
        tb.addWidget(self.metric)
        self.log = QCheckBox("Log scale")
        self.log.setToolTip("Logarithmic y axis (useful for CEM43)")
        self.log.toggled.connect(self._plot)
        tb.addWidget(self.log)
        self.shade = QCheckBox("Shade sonication")
        self.shade.setChecked(True)
        self.shade.setToolTip("Grey bands where the protocol sonicates")
        self.shade.toggled.connect(self._plot)
        tb.addWidget(self.shade)
        tb.addStretch(1)
        self.note = QLabel()
        self.note.setObjectName("CardSub")
        self.note.setWordWrap(True)
        tb.addWidget(self.note, 2)
        v.addLayout(tb)

        body = QHBoxLayout()
        side = QVBoxLayout()
        side.addWidget(QLabel("Tissues"))
        self.tissues = QListWidget()
        self.tissues.setMaximumWidth(200)
        self.tissues.itemChanged.connect(self._tissue_toggled)
        side.addWidget(self.tissues, 1)
        row = QHBoxLayout()
        for text, state in (("All", True), ("None", False)):
            b = QPushButton(text)
            b.clicked.connect(lambda _c, s=state: self._set_all(s))
            row.addWidget(b)
        side.addLayout(row)
        body.addLayout(side)

        plot = QVBoxLayout()
        if HAVE_MPL:
            self.fig = Figure(figsize=(8, 4), tight_layout=True)
            # export through the canvas's own button (export.py); no matplotlib toolbar
            self.canvas = FigureCanvasQTAgg(self.fig, export_name="heating", filename_hint=filename_hint)
            self.ax = self.fig.add_subplot(1, 1, 1)
            plot.addWidget(self.canvas, 1)
        else:
            plot.addWidget(QLabel("Install matplotlib to plot the time courses."), 1)
        body.addLayout(plot, 1)
        v.addLayout(body, 1)
        self.stack.show_empty("No run selected")

    # ── loading ─────────────────────────────────────────────────────────
    def load(self, run: rs.Run | None) -> None:
        path = run.heat if run is not None else None
        try:
            key = (str(path), path.stat().st_mtime_ns) if path is not None else None
        except OSError:
            key = (str(path), None)
        if key == self._key and self.data is not None:
            return
        self._key = key
        self.data = None
        if path is None:
            self.note.setText("")
            self.stack.show_empty("No run selected")
        else:
            try:
                self.data = rs.read_heat(Path(path))
                self.note.setText(path.name)
                self.stack.show_content()
            except rs.HeatReadError as exc:
                self.note.setText(str(exc))
                if not Path(path).is_file():
                    self.stack.show_empty("No heating time courses for this run",
                                          "They come from the thermal stage "
                                          "(modules.run_heating_sims), which writes "
                                          f"{Path(path).name} to the cache folder.")
                else:
                    self.stack.show_empty("Cannot read this run's heating time courses", str(exc))
        current = self.metric.currentText()
        labels = list(self.data.metrics) if self.data is not None else []
        self.metric.blockSignals(True)
        self.metric.clear()
        self.metric.addItems(labels)
        if current in labels:
            self.metric.setCurrentText(current)
        self.metric.blockSignals(False)
        self._metric_changed()

    def _series(self) -> rs.HeatSeries | None:
        if self.data is None:
            return None
        return self.data.metrics.get(self.metric.currentText())

    def _metric_changed(self, *_args) -> None:
        series = self._series()
        self.tissues.blockSignals(True)
        self.tissues.clear()
        for name in (series.tissues if series is not None else {}):
            item = QListWidgetItem(name.replace("_", " "))
            item.setData(Qt.ItemDataRole.UserRole, name)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Unchecked if name.lower() in self._hidden
                               else Qt.CheckState.Checked)
            self.tissues.addItem(item)
        self.tissues.blockSignals(False)
        self._plot()

    def _tissue_toggled(self, item: QListWidgetItem) -> None:
        name = item.data(Qt.ItemDataRole.UserRole).lower()
        if item.checkState() == Qt.CheckState.Checked:
            self._hidden.discard(name)
        else:
            self._hidden.add(name)
        self._plot()

    def _set_all(self, state: bool) -> None:
        self.tissues.blockSignals(True)
        for i in range(self.tissues.count()):
            item = self.tissues.item(i)
            item.setCheckState(Qt.CheckState.Checked if state else Qt.CheckState.Unchecked)
            name = item.data(Qt.ItemDataRole.UserRole).lower()
            (self._hidden.discard if state else self._hidden.add)(name)
        self.tissues.blockSignals(False)
        self._plot()

    def shown_tissues(self) -> list[str]:
        return [self.tissues.item(i).data(Qt.ItemDataRole.UserRole)
                for i in range(self.tissues.count())
                if self.tissues.item(i).checkState() == Qt.CheckState.Checked]

    def _colour(self, name: str) -> str:
        key = name.lower()
        if key not in self._colours:
            self._colours[key] = PALETTE[len(self._colours) % len(PALETTE)]
        return self._colours[key]

    # ── plotting ────────────────────────────────────────────────────────
    def _plot(self, *_args) -> None:
        if not HAVE_MPL:
            return
        ax = self.ax
        ax.clear()
        series = self._series()
        if series is None:
            self.canvas.draw_idle()
            return
        shown = self.shown_tissues()
        for name in shown:
            y = series.tissues[name]
            ax.plot(series.time, y, lw=1.8, color=self._colour(name),
                    label=f"{name.replace('_', ' ')} (max {_fmt(y)} {series.unit})")
        if self.shade.isChecked() and self.data is not None:
            for t0, t1 in self.data.on:
                ax.axvspan(t0, t1, color="0.5", alpha=0.12, lw=0)
        self.shade.setEnabled(bool(self.data is not None and self.data.on))
        self.shade.setToolTip("Grey bands where the protocol sonicates" if self.shade.isEnabled()
                              else "This run's HEAT file has no protocol steps "
                                   "(written by PRESTUS since 2026-10-05)")
        ax.set_yscale("log" if self.log.isChecked() and shown else "linear")
        ax.set_xlabel("Time [s]", fontsize=style.PLOT_PT["normal"])
        ax.set_ylabel(f"{self.metric.currentText()} [{series.unit}]", fontsize=style.PLOT_PT["normal"])
        ax.tick_params(labelsize=style.PLOT_PT["small"])
        if len(series.time):
            ax.set_xlim(series.time[0], series.time[-1])
        ax.grid(True, alpha=0.3)
        if shown:
            ax.legend(frameon=False, fontsize=style.PLOT_PT["small"])
        self.canvas.draw_idle()


def _fmt(y) -> str:
    import numpy as np
    finite = y[np.isfinite(y)]
    return f"{finite.max():.3g}" if finite.size else "n/a"
