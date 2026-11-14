"""Multi-transducer coupling card (Sonication tab): an optional analytic prediction of the
coherent and independent fields (coupling.py), from the placement viewer's positions."""

from __future__ import annotations

import numpy as np
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import QCheckBox, QLabel, QVBoxLayout, QWidget

from .. import coupling as cp
from . import style
from .empty_state import EmptyStack
from .export import subject_stem

try:
    from .mpl_theme import ThemedCanvas as FigureCanvasQTAgg   # drawn in the Nexus theme
    from matplotlib.figure import Figure
    HAVE_MPL = True
except ImportError:          # pragma: no cover
    HAVE_MPL = False


def _num(v):
    if isinstance(v, (list, tuple)):
        v = v[0] if v else None
    return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else None


class CouplingPrediction(QWidget):
    """The prediction card; with no beams to draw, an empty state whose one action asks for
    the placement viewer (``open_placement_requested``, connected by the main window)."""

    open_placement_requested = pyqtSignal()

    def __init__(self, model, viewer_getter, parent=None):
        super().__init__(parent)
        self.model = model
        self.viewer_getter = viewer_getter          # -> PlacementPage or None (positions in RAS mm)
        self._shown_for = None
        v = QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 8)
        self.toggle = QCheckBox("Show an analytic prediction of coherent vs. async fields")
        self.toggle.setToolTip("Free water, Gaussian focal spots from each transducer's aperture and "
                               "focal distance, at its target free-water ISPPA. A guide to the choice "
                               "below, not a simulation: no skull, no side lobes.")
        self.toggle.toggled.connect(lambda _on: self.refresh())
        v.addWidget(self.toggle)
        self.text = QLabel()
        self.text.setObjectName("CardSub")
        self.text.setWordWrap(True)
        self.text.setTextFormat(Qt.TextFormat.RichText)
        v.addWidget(self.text)
        self.canvas = None
        self.stack = None
        if HAVE_MPL:
            self.canvas = FigureCanvasQTAgg(Figure(figsize=(8, 3.2), layout="constrained"),
                                            export_name="coupling",
                                            filename_hint=subject_stem(model, "coupling"))
            self.canvas.setMinimumHeight(260)
            self.stack = EmptyStack(self.canvas)
            v.addWidget(self.stack)
        self.refresh()

    # ── inputs ──────────────────────────────────────────────────────────
    def beams(self) -> tuple[list[cp.Beam], str]:
        """Beams from the viewer's positions and the effective transducers; why not, if empty."""
        txs = self.model.effective_transducers()
        isppas = [_num(t.get("target_isppa_wcm2")) for t in txs]
        viewer = self.viewer_getter() if self.viewer_getter else None
        if viewer is None:
            return [], "positions need the placement viewer"
        out = []
        for i, t in enumerate(txs):
            trans, focus = viewer._positions_for(i)
            bowl = viewer._bowl_for(i)
            f = _num(t.get("freq_hz"))
            if trans is None or focus is None or bowl is None or f is None or not isppas[i]:
                return [], (f"transducer {i + 1}: position, geometry, frequency or target ISPPA missing "
                            "(Placement and Transducer tabs)")
            if np.linalg.norm(np.asarray(focus) - np.asarray(trans)) <= 0:
                return [], f"transducer {i + 1}: transducer and focus positions coincide"
            out.append(cp.Beam(np.asarray(trans, float), np.asarray(focus, float),
                               float(bowl.aperture), f, float(isppas[i])))
        return out, ""

    # ── drawing ─────────────────────────────────────────────────────────
    def refresh(self) -> None:
        n = self.model.n_transducers
        self.setVisible(n > 1)
        on = self.toggle.isChecked()
        if self.stack is not None:
            self.stack.setVisible(on)
        if n < 2 or not on:
            self.text.setText("")
            return
        txs = self.model.effective_transducers()
        isppas = [_num(t.get("target_isppa_wcm2")) for t in txs]
        lines = []
        if all(isppas):
            t = cp.at_target(isppas)
            lines.append(
                f"Where the foci coincide, from target ISPPAs of {', '.join(f'{x:g}' for x in isppas)} W/cm²: "
                f"<b>async</b> {t['independent']:.3g} W/cm²; <b>coherent</b> "
                f"{t['coherent_min']:.3g}–{t['coherent_max']:.3g} W/cm², from opposite to equal phase. "
                "PRESTUS's coherent run starts every source at t = 0, so the phase at the target follows "
                "from the path lengths and is known only after the simulation.")
        beams, why = self.beams()
        if not beams:
            self.text.setText("<br>".join(lines))
            if self.stack is not None:
                self.canvas.figure.clear()
                self._shown_for = None
                self.stack.show_empty("No map of the fields yet", f"The map needs {why}.",
                                      "Open placement viewer", self.open_placement_requested.emit,
                                      "Positions come from Placement ▸ Viewer")
            return
        if self.stack is not None:
            self.stack.show_content()
        theta = cp.beam_angle_deg(beams)
        lam = beams[0].wavelength_mm()
        fringe = cp.fringe_spacing_mm(theta, lam)
        spread = max(float(np.linalg.norm(b.focus - beams[0].focus)) for b in beams)
        lines.append(f"Beams {theta:.0f}° apart, foci up to {spread:.1f} mm apart; λ = {lam:.2f} mm in water"
                     + (f"; coherent fringes about {fringe:.1f} mm apart." if fringe else "."))
        self.text.setText("<br>".join(lines))
        sig = (tuple((tuple(b.apex), tuple(b.focus), b.aperture_mm, b.freq_hz, b.isppa_wcm2) for b in beams),
               style.active()["surface"])
        if self.canvas is None or sig == self._shown_for:
            return
        self._shown_for = sig
        extent = max(15.0, 1.5 * max(b.widths_mm()[1] for b in beams) / 2 + spread)
        maps = cp.field_maps(beams, extent_mm=extent)
        fig = self.canvas.figure
        fig.clear()
        vmax = float(max(maps["independent"].max(), maps["coherent"].max())) or 1.0
        axes = fig.subplots(1, 2, sharey=True)
        for ax, key, title in ((axes[0], "independent", "Async: Σ|p|²"),
                               (axes[1], "coherent", "Coherent, in phase at the foci: |Σp|²")):
            im = ax.imshow(maps[key], origin="lower", cmap="magma", vmin=0, vmax=vmax,
                           extent=(-extent, extent, -extent, extent))
            # text and axes colours come from the theme (mpl_theme), so exports re-theme too
            ax.set_title(f"{title}\npeak {maps[key].max():.3g} W/cm²", fontsize=style.PLOT_PT["title"])
            ax.set_xlabel("along beam 1 [mm]", fontsize=style.PLOT_PT["normal"])
            ax.tick_params(labelsize=style.PLOT_PT["small"])
        axes[0].set_ylabel("across, in the beams' plane [mm]", fontsize=style.PLOT_PT["normal"])
        cb = fig.colorbar(im, ax=list(axes), shrink=0.9)
        cb.ax.tick_params(labelsize=style.PLOT_PT["small"])
        cb.set_label("Intensity [W/cm²]", fontsize=style.PLOT_PT["normal"])
        self.canvas.draw_idle()
