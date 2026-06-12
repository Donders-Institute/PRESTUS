"""Interactive 3-plane NIfTI viewer embedded in a PyQt6 widget.

Shows axial / coronal / sagittal slices of an anatomical volume with
an optional overlay (e.g. intensity map).  Clicking any slice moves the
crosshair and emits a ``cursor_moved`` signal with the voxel coordinates.
"""

from __future__ import annotations
import numpy as np

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QSlider,
    QComboBox, QPushButton, QSizePolicy,
)
from PyQt6.QtCore import Qt, pyqtSignal

try:
    import matplotlib
    matplotlib.use("QtAgg")
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.figure import Figure
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

try:
    import nibabel as nib
    HAS_NIB = True
except ImportError:
    HAS_NIB = False


COLORMAPS_OVERLAY = ["hot", "plasma", "viridis", "jet", "coolwarm"]


class NiftiViewer(QWidget):
    """Three-plane slice viewer with optional overlay.

    Signals
    -------
    cursor_moved(x, y, z) : voxel indices (int, int, int)
    """

    cursor_moved = pyqtSignal(int, int, int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._anat_data  = None   # (X, Y, Z) float32 array, normalised 0-1
        self._over_data  = None   # same shape, values = metric (or None)
        self._anat_affine = None
        self._cx = 0
        self._cy = 0
        self._cz = 0
        self._alpha = 0.5

        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(4)

        # ── toolbar ──────────────────────────────────────────────────────────
        tb = QHBoxLayout()
        tb.setSpacing(8)

        self._cmap_cb = QComboBox()
        self._cmap_cb.addItems(COLORMAPS_OVERLAY)
        self._cmap_cb.setFixedWidth(110)
        self._cmap_cb.currentTextChanged.connect(self._redraw)
        tb.addWidget(QLabel("Overlay colormap:"))
        tb.addWidget(self._cmap_cb)

        self._alpha_sl = QSlider(Qt.Orientation.Horizontal)
        self._alpha_sl.setRange(0, 100)
        self._alpha_sl.setValue(50)
        self._alpha_sl.setFixedWidth(120)
        self._alpha_sl.valueChanged.connect(self._on_alpha_changed)
        tb.addWidget(QLabel("Opacity:"))
        tb.addWidget(self._alpha_sl)

        self._coord_lbl = QLabel("Cursor: —")
        self._coord_lbl.setStyleSheet("font-size:11px; color:#666;")
        tb.addStretch()
        tb.addWidget(self._coord_lbl)
        layout.addLayout(tb)

        # ── matplotlib canvas ─────────────────────────────────────────────────
        if HAS_MPL:
            self._fig = Figure(facecolor="#1a1a1a", tight_layout=True)
            self._canvas = FigureCanvas(self._fig)
            self._canvas.setSizePolicy(
                QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
            )
            layout.addWidget(self._canvas, stretch=1)
            self._axes = [
                self._fig.add_subplot(1, 3, i + 1) for i in range(3)
            ]
            for ax in self._axes:
                ax.set_facecolor("#1a1a1a")
                ax.set_xticks([])
                ax.set_yticks([])
            self._canvas.mpl_connect("button_press_event", self._on_click)
            self._im_anat  = [None, None, None]
            self._im_over  = [None, None, None]
            self._hlines   = [None, None, None]
            self._vlines   = [None, None, None]
            self._init_axes()
        else:
            no_mpl = QLabel("matplotlib not available — install it to use the viewer.")
            no_mpl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            layout.addWidget(no_mpl)

    # ── Public API ────────────────────────────────────────────────────────────

    def load_anatomy(self, path: str):
        if not HAS_NIB:
            return
        img = nib.load(path)
        data = np.asanyarray(img.dataobj).astype(np.float32)
        # Normalise
        vmin, vmax = np.percentile(data[data > 0], [1, 99]) if data.any() else (0, 1)
        data = np.clip((data - vmin) / max(vmax - vmin, 1e-9), 0, 1)
        self._anat_data   = data
        self._anat_affine = img.affine
        self._cx = data.shape[0] // 2
        self._cy = data.shape[1] // 2
        self._cz = data.shape[2] // 2
        self._redraw()

    def load_overlay(self, path: str):
        if not HAS_NIB:
            return
        img = nib.load(path)
        data = np.asanyarray(img.dataobj).astype(np.float32)
        self._over_data = data
        self._redraw()

    def clear_overlay(self):
        self._over_data = None
        self._redraw()

    def set_cursor(self, x: int, y: int, z: int):
        if self._anat_data is None:
            return
        sx, sy, sz = self._anat_data.shape
        self._cx = int(np.clip(x, 0, sx - 1))
        self._cy = int(np.clip(y, 0, sy - 1))
        self._cz = int(np.clip(z, 0, sz - 1))
        self._update_crosshairs()
        self._update_coord_label()

    def cursor(self):
        return self._cx, self._cy, self._cz

    # ── Internal ─────────────────────────────────────────────────────────────

    def _init_axes(self):
        titles = ["Sagittal (YZ)", "Coronal (XZ)", "Axial (XY)"]
        for i, ax in enumerate(self._axes):
            ax.set_title(titles[i], color="#aaa", fontsize=10)

    def _redraw(self):
        if not HAS_MPL or self._anat_data is None:
            return
        cmap = self._cmap_cb.currentText() if hasattr(self, "_cmap_cb") else "hot"
        d = self._anat_data
        o = self._over_data

        slices_anat = [
            d[self._cx, :, :],   # sagittal
            d[:, self._cy, :],   # coronal
            d[:, :, self._cz],   # axial
        ]
        slices_over = None
        if o is not None:
            slices_over = [
                o[self._cx, :, :],
                o[:, self._cy, :],
                o[:, :, self._cz],
            ]

        for i, ax in enumerate(self._axes):
            ax.cla()
            ax.set_facecolor("#1a1a1a")
            ax.set_xticks([])
            ax.set_yticks([])
            sa = np.rot90(slices_anat[i])
            ax.imshow(sa, cmap="gray", vmin=0, vmax=1, aspect="auto",
                      interpolation="nearest")
            if slices_over is not None:
                so = np.rot90(slices_over[i])
                vmax = np.nanpercentile(o[o > 0], 99) if o.any() else 1
                ax.imshow(so, cmap=cmap, alpha=self._alpha,
                          vmin=0, vmax=vmax, aspect="auto",
                          interpolation="nearest")

        self._add_crosshairs()
        titles = ["Sagittal (YZ)", "Coronal (XZ)", "Axial (XY)"]
        for i, ax in enumerate(self._axes):
            ax.set_title(titles[i], color="#aaa", fontsize=10)
        self._canvas.draw_idle()

    def _add_crosshairs(self):
        if not self._axes:
            return
        # Sagittal: Y horizontal, Z vertical  → crosshair at (cy, cz)
        # Coronal:  X horizontal, Z vertical  → crosshair at (cx, cz)
        # Axial:    X horizontal, Y vertical  → crosshair at (cx, cy)
        shapes = [
            self._anat_data.shape[1:],   # sagittal: YZ
            (self._anat_data.shape[0], self._anat_data.shape[2]),  # coronal: XZ
            self._anat_data.shape[:2],   # axial: XY
        ]
        h_pos = [self._cy, self._cx, self._cx]
        v_pos = [self._cz, self._cz, self._cy]
        for i, ax in enumerate(self._axes):
            h, w = shapes[i]
            y_cross = h - 1 - v_pos[i]   # because rot90 flips axis
            x_cross = h_pos[i]
            ax.axhline(y_cross, color="#0b7b99", lw=0.8, alpha=0.8)
            ax.axvline(x_cross, color="#0b7b99", lw=0.8, alpha=0.8)

    def _update_crosshairs(self):
        self._redraw()

    def _update_coord_label(self):
        self._coord_lbl.setText(
            f"Cursor: x={self._cx}, y={self._cy}, z={self._cz}"
        )

    def _on_alpha_changed(self, val: int):
        self._alpha = val / 100.0
        self._redraw()

    def _on_click(self, event):
        if event.inaxes is None or self._anat_data is None:
            return
        ax_idx = self._axes.index(event.inaxes) if event.inaxes in self._axes else -1
        if ax_idx < 0:
            return
        shapes = [
            self._anat_data.shape[1:],
            (self._anat_data.shape[0], self._anat_data.shape[2]),
            self._anat_data.shape[:2],
        ]
        h, w = shapes[ax_idx]
        px = int(round(event.xdata))
        py = int(round(event.ydata))
        # Undo the rot90 that was applied when displaying
        actual_v = h - 1 - py   # vertical click → original axis

        if ax_idx == 0:    # sagittal: horiz=Y, vert=Z
            self._cy = np.clip(px, 0, self._anat_data.shape[1] - 1)
            self._cz = np.clip(actual_v, 0, self._anat_data.shape[2] - 1)
        elif ax_idx == 1:  # coronal: horiz=X, vert=Z
            self._cx = np.clip(px, 0, self._anat_data.shape[0] - 1)
            self._cz = np.clip(actual_v, 0, self._anat_data.shape[2] - 1)
        else:              # axial: horiz=X, vert=Y
            self._cx = np.clip(px, 0, self._anat_data.shape[0] - 1)
            self._cy = np.clip(actual_v, 0, self._anat_data.shape[1] - 1)

        self._update_crosshairs()
        self._update_coord_label()
        self.cursor_moved.emit(self._cx, self._cy, self._cz)
