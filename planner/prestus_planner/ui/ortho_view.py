"""Reusable read-only 3-plane orthogonal viewer: an anatomy volume with an optional
colour overlay, both loaded as ``placement.Volume`` and sampled by RAS coordinate
(like PlacementPage._overlay), so the overlay need not share the anatomy's grid or
orientation. Used by the Results tab (acoustic/thermal maps) and the pseudo-CT
viewer — kept separate from PlacementPage, which also draws bowl/transducer geometry
and writes clicks back into the study config."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (QComboBox, QHBoxLayout, QLabel, QSlider, QVBoxLayout, QWidget)

from .. import placement as pl

try:
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
    from matplotlib.figure import Figure
    HAVE_MPL = True
except ImportError:          # pragma: no cover
    HAVE_MPL = False

VIEWS = [("Axial", 2, (0, 1)), ("Coronal", 1, (0, 2)), ("Sagittal", 0, (1, 2))]
COLORMAPS = ["hot", "plasma", "viridis", "jet", "coolwarm"]


class OrthoViewer(QWidget):
    """Axial/coronal/sagittal slices through an anatomy volume, with an optional
    colour overlay (own colormap + opacity). No click-to-place; view only."""

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.anat: pl.Volume | None = None
        self.overlay: pl.Volume | None = None
        self.cursor = np.zeros(3)

        v = QVBoxLayout(self)
        v.setContentsMargins(4, 4, 4, 4)
        v.setSpacing(4)

        tb = QHBoxLayout()
        tb.addWidget(QLabel("Overlay colormap"))
        self.cmap = QComboBox()
        self.cmap.addItems(COLORMAPS)
        self.cmap.setFixedWidth(110)
        self.cmap.currentTextChanged.connect(self._draw)
        tb.addWidget(self.cmap)
        tb.addWidget(QLabel("Opacity"))
        self.alpha = QSlider(Qt.Orientation.Horizontal, minimum=0, maximum=100, value=50)
        self.alpha.setMaximumWidth(140)
        self.alpha.valueChanged.connect(self._draw)
        tb.addWidget(self.alpha)
        tb.addStretch(1)
        self.status = QLabel()
        self.status.setStyleSheet("color: gray;")
        tb.addWidget(self.status)
        v.addLayout(tb)

        if HAVE_MPL:
            self.fig = Figure(figsize=(9, 3.2), tight_layout=True)
            self.canvas = FigureCanvasQTAgg(self.fig)
            self.axes = [self.fig.add_subplot(1, 3, i + 1) for i in range(3)]
            self.canvas.mpl_connect("scroll_event", self._scrolled)
            v.addWidget(self.canvas, 1)
        else:
            v.addWidget(QLabel("Install matplotlib and nibabel to use this viewer."), 1)

    # ── loading ─────────────────────────────────────────────────────────
    def set_anatomy(self, path: str | Path | None) -> None:
        self.anat = None
        if path is not None and Path(path).is_file():
            try:
                self.anat = pl.Volume(path)
                self.cursor = self.anat.to_ras(np.array(self.anat.shape) / 2)[0]
            except Exception as exc:
                self.status.setText(f"Cannot read anatomy: {exc}")
                self._draw()
                return
        self._draw()

    def set_overlay(self, path: str | Path | None) -> None:
        self.overlay = None
        if path is not None and Path(path).is_file():
            try:
                self.overlay = pl.Volume(path)
            except Exception as exc:
                self.status.setText(f"Cannot read overlay: {exc}")
                self._draw()
                return
        self._draw()

    # ── drawing ─────────────────────────────────────────────────────────
    def _overlay_slice(self, plane: int, k: int, shape2d: tuple[int, int]) -> np.ndarray | None:
        if self.overlay is None or self.anat is None:
            return None
        ii, jj = np.meshgrid(np.arange(shape2d[0]), np.arange(shape2d[1]), indexing="ij")
        vox = np.zeros((ii.size, 3))
        others = [a for a in range(3) if a != plane]
        vox[:, others[0]] = ii.ravel()
        vox[:, others[1]] = jj.ravel()
        vox[:, plane] = k
        values = self.overlay.sample(self.anat.canon_to_ras(vox)).astype(float)
        return values.reshape(shape2d)

    def _draw(self) -> None:
        if not HAVE_MPL:
            return
        for ax in self.axes:
            ax.clear()
            ax.set_xticks([])
            ax.set_yticks([])
        if self.anat is None:
            self.canvas.draw_idle()
            return
        data = self.anat.canon
        vs = np.sqrt((self.anat.canon_affine[:3, :3] ** 2).sum(axis=0))
        lo, hi = np.percentile(data[::4, ::4, ::4], [1, 99.5])
        c = self.anat.ras_to_canon(self.cursor)[0]
        over_lo = over_hi = None
        if self.overlay is not None:
            nz = self.overlay.data[self.overlay.data != 0]
            if nz.size:
                over_lo, over_hi = np.percentile(nz, [1, 99])
        for ax, (title, plane, (h, w)) in zip(self.axes, VIEWS):
            k = int(np.clip(round(c[plane]), 0, data.shape[plane] - 1))
            sl = np.take(data, k, axis=plane)
            ax.imshow(sl.T, cmap="gray", origin="lower", vmin=lo, vmax=hi,
                      aspect=vs[w] / vs[h], interpolation="nearest")
            over = self._overlay_slice(plane, k, sl.shape) if self.alpha.value() > 0 else None
            if over is not None and over_hi is not None and over_hi > over_lo:
                masked = np.where(over != 0, over, np.nan)
                ax.imshow(masked.T, origin="lower", cmap=self.cmap.currentText(),
                          vmin=over_lo, vmax=over_hi, alpha=self.alpha.value() / 100,
                          aspect=vs[w] / vs[h], interpolation="nearest")
            ax.set_title(title, fontsize=9)
            ax.set_xlim(-0.5, sl.shape[0] - 0.5)
            ax.set_ylim(-0.5, sl.shape[1] - 0.5)
        self.canvas.draw_idle()

    def _scrolled(self, event) -> None:
        if self.anat is None or event.inaxes not in self.axes:
            return
        _, plane, _ = VIEWS[self.axes.index(event.inaxes)]
        c = self.anat.ras_to_canon(self.cursor)[0]
        c[plane] += 1 if event.button == "up" else -1
        self.cursor = self.anat.canon_to_ras(c[None])[0]
        self._draw()
