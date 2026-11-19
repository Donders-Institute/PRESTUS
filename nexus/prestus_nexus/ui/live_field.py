"""Live-view data frames (simulation.live_view), drawn in the Nexus theme.

A data frame is a 2D plane of the field written by live_view.m ('data'), run_jwave_sim.py,
or bb_progress.py: <name>.f32 (little-endian float32) and <name>.json describing it
(quantity, units, label, spacing_mm, shape, order 'F' or 'C', min, max, markers as 0-based
plane indices, data: the .f32 file name). The plane's first axis is drawn across, its second
(the beam axis in PRESTUS grids) downwards, in mm. Beside a temperature field, the peak
temperature over the run (progress.RunProgress.temperature: the hottest voxel of the grid,
not split by tissue) is drawn as a curve.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from matplotlib.figure import Figure
from PyQt6.QtGui import QPalette
from PyQt6.QtWidgets import QVBoxLayout, QWidget

from . import style
from .mpl_theme import ThemedCanvas

MARKERS = {"trans": ("^", "transducer"), "focus": ("x", "focus")}


def load_frame(json_path: str) -> tuple[dict, np.ndarray] | None:
    """(description, plane) of a data frame; None while it is incomplete or unreadable."""
    try:
        meta = json.loads(Path(json_path).read_text(encoding="utf-8"))
        raw = Path(json_path).with_name(str(meta["data"]))
        shape = tuple(int(n) for n in meta["shape"])
        data = np.fromfile(raw, dtype="<f4")
        if data.size != int(np.prod(shape)):
            return None                         # replaced while reading: the next one will do
        return meta, data.reshape(shape, order="F" if meta.get("order") == "F" else "C")
    except (OSError, ValueError, KeyError, TypeError):
        return None


class FieldCanvas(QWidget):
    """One data frame (and, for temperature, the peak temperature over time)."""

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        v = QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 0)
        self.figure = Figure(figsize=(6, 6), layout="constrained")
        self.canvas = ThemedCanvas(self.figure, export_name="live_view")
        v.addWidget(self.canvas)

    def draw_frame(self, meta: dict, plane: np.ndarray, curve: list | None = None) -> None:
        fig = self.figure
        fig.clear()
        # the page's own background (Qt's window colour, as around the tab), not a theme role:
        # mpl_theme leaves it alone (role False)
        fig.set_facecolor(self.palette().color(QPalette.ColorRole.Window).name())
        fig.patch._nexus_bg = False
        has_curve = bool(curve) and meta.get("quantity") == "temperature"
        if has_curve:                                   # the field keeps the full height; the
            gs = fig.add_gridspec(1, 2, width_ratios=[3, 2], wspace=0.15)   # curve sits beside it
            ax, axc = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
            axc.set_facecolor(style.active()["surface"])
            axc.set_box_aspect(0.75)
            axc.set_anchor("C")
        else:
            ax, axc = fig.add_subplot(111), None
        dx, dz = (list(meta.get("spacing_mm") or [1.0, 1.0]) + [1.0, 1.0])[:2]
        nx, nz = plane.shape[0], plane.shape[1] if plane.ndim > 1 else 1
        img = np.atleast_2d(plane).T                     # 2nd axis downwards
        vmin, vmax = float(np.nanmin(img)), float(np.nanmax(img))
        if meta.get("quantity") == "temperature":
            cmap = "inferno"
        elif vmin < 0 < vmax:                             # signed pressure: symmetric scale
            cmap, m = "RdBu_r", max(abs(vmin), abs(vmax))
            vmin, vmax = -m, m
        else:
            cmap = "viridis"
        im = ax.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax, origin="upper", aspect="equal",
                       extent=(0, nx * dx, nz * dz, 0), interpolation="nearest")
        ax.set_anchor("C")
        ax.set_xlabel("mm")
        ax.set_ylabel("mm (depth)")
        units = meta.get("units") or ""
        label = f"{str(meta.get('quantity') or '').capitalize()}" + (f" [{units}]" if units else "")
        cax = ax.inset_axes([1.04, 0.0, 0.035, 1.0])    # attached to the field, whatever its size
        fig.colorbar(im, cax=cax, label=label)
        for key, (marker, name) in MARKERS.items():
            pos = (meta.get("markers") or {}).get(key)
            if isinstance(pos, (list, tuple)) and len(pos) >= 2:
                ax.plot((pos[0] + 0.5) * dx, (pos[1] + 0.5) * dz, marker, color="#ffffff",
                        markeredgecolor="#000000", markersize=9, label=name, linestyle="none")
        if any((meta.get("markers") or {}).get(k) for k in MARKERS):
            ax.legend(loc="lower right", fontsize=8, framealpha=0.7)
        if axc is not None:
            xs = [c[0] for c in curve]
            ys = [c[1] for c in curve]
            axc.plot(xs, ys, color="#D55E00", linewidth=1.5)
            axc.set_xlabel("simulated time [s]" if curve[-1][2] else "progress [%]")
            kind = curve[-1][3]
            axc.set_ylabel("peak [°C]" if kind == "Tmax" else "max now [°C]")
        self.canvas.draw_idle()
