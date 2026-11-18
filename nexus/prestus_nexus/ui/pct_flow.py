"""Skull CT sub-tab card: the steps from the input image to the skull's acoustic properties, as
the run will take them, with the mappings drawn:

  measured CT:  CT (HU) → resampled to the segmentation grid → HU → density, sound speed, attenuation
  pseudo-CT:    UTE → bias correction and soft-tissue normalisation → UTE → pseudo-HU (line)
                → pseudo-CT (HU) → HU → density, sound speed, attenuation

The UTE → HU line comes from pct_create_pseudoCT.sh (pct_lines.py). The HU → property curves are
evaluated by MATLAB with the configured mappings (nexus_pct 'curves' → pct_mapping_curves.m), so
they are the run's own code; they are requested when MATLAB is already running, or on the button.
"""

from __future__ import annotations

from pathlib import Path

from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QPixmap
from PyQt6.QtWidgets import QFrame, QHBoxLayout, QLabel, QPushButton, QVBoxLayout, QWidget

from .. import run_plan as rp
from ..model import is_empty
from ..outputs import dir_pct
from ..pct_lines import read_script, ute_line
from . import icons, style

try:
    from matplotlib.figure import Figure

    from .mpl_theme import ThemedCanvas
    HAVE_MPL = True
except ImportError:          # pragma: no cover
    HAVE_MPL = False

STATE_ICON = {rp.RUN: ("ok", "green"), rp.REUSE: ("reset", "info"), rp.SKIP: ("minus", "ink4"),
              rp.STOP: ("warning", "red")}
PROPS = [("density", "Density", "kg/m³"), ("sound_speed", "Sound speed", "m/s"),
         ("alpha_coeff", "Attenuation", "dB/(MHz^y cm)")]
_CURVES: dict = {}           # request signature -> curves (MATLAB results of this session)


class _Step(QWidget):
    """One step: state icon, an icon of what the step does, title and a line of detail."""

    def __init__(self, parent=None):
        super().__init__(parent)
        h = QHBoxLayout(self)
        h.setContentsMargins(0, 2, 0, 2)
        self.icon = QLabel()
        self.icon.setFixedWidth(20)
        self.icon.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.what = QLabel()
        self.what.setFixedWidth(26)
        self.what.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.glyph = ""
        self.text = QLabel()
        self.text.setWordWrap(True)
        self.text.setTextFormat(Qt.TextFormat.RichText)
        h.addWidget(self.icon)
        h.addWidget(self.what)
        h.addWidget(self.text, 1)

    def set(self, state: str, title: str, detail: str, glyph: str = "") -> None:
        name, colour = STATE_ICON[state]
        c = style.active()
        self.icon.setPixmap(icons.icon(name, c[colour]).pixmap(16, 16))
        self.glyph = glyph
        self.what.setPixmap(icons.icon(glyph, c["ink3"] if state == rp.SKIP else c["ink2"]).pixmap(22, 22)
                            if glyph else QPixmap())
        self.text.setText(f"<b>{title}</b><br><span style='color:{c['ink3']}'>{detail}</span>")


class PctFlow(QFrame):
    """'From image to skull properties'."""

    def __init__(self, window, parent=None):
        super().__init__(parent)
        self.window = window
        self.setObjectName("Card")
        v = QVBoxLayout(self)
        t = QLabel("From image to skull properties")
        t.setObjectName("CardTitle")
        v.addWidget(t)
        self.sub = QLabel()
        self.sub.setObjectName("CardSub")
        self.sub.setWordWrap(True)
        v.addWidget(self.sub)
        self.steps_box = QVBoxLayout()
        v.addLayout(self.steps_box)
        self.steps = [_Step() for _ in range(5)]
        self.arrows = [QLabel("↓") for _ in range(4)]
        self.ute_canvas = self.prop_canvas = None
        if HAVE_MPL:
            self.ute_fig = Figure(figsize=(4.6, 1.9), tight_layout=True)
            self.ute_canvas = ThemedCanvas(self.ute_fig, export_name="pct_ute_to_hu")
            self.ute_canvas.setMinimumHeight(170)
            self.prop_fig = Figure(figsize=(4.6, 4.6), tight_layout=True)
            self.prop_canvas = ThemedCanvas(self.prop_fig, export_name="pct_hu_to_properties")
            self.prop_canvas.setMinimumHeight(380)
        self.curves_note = QLabel()
        self.curves_note.setObjectName("CardSub")
        self.curves_note.setWordWrap(True)
        self.curves_btn = QPushButton("Draw the mapping curves")
        self.curves_btn.setToolTip("Evaluates the configured HU → density, sound speed, and attenuation "
                                   "mappings in MATLAB (pct_mapping_curves); starts MATLAB if it is off")
        self.curves_btn.clicked.connect(lambda: self.request_curves(force=True))
        for i, s in enumerate(self.steps):
            self.steps_box.addWidget(s)
            if i == 2 and self.ute_canvas is not None:
                self.steps_box.addWidget(self.ute_canvas)
            if i == 4:
                if self.prop_canvas is not None:
                    self.steps_box.addWidget(self.prop_canvas)
                row = QHBoxLayout()
                row.addWidget(self.curves_note, 1)
                row.addWidget(self.curves_btn)
                self.steps_box.addLayout(row)
            if i < 4:
                a = self.arrows[i]
                a.setObjectName("CardSub")
                a.setContentsMargins(4, 0, 0, 0)
                self.steps_box.addWidget(a)
        self._timer = QTimer(self, singleShot=True, interval=700)
        self._timer.timeout.connect(self.request_curves)
        self._pending: tuple | None = None
        self._shown_sig = None

    @property
    def model(self):
        return self.window.model

    # ── steps ───────────────────────────────────────────────────────────
    def refresh(self) -> None:
        m = self.model
        g = m.get
        enabled = g("pct.enabled") in (1, True)
        measured = not is_empty(g("path.ct_pattern"))
        ct = rp._skull_ct(m, True)
        self.sub.setText("CT-based skull mapping is off: the skull keeps the uniform properties of the "
                         "Acoustic properties table. The steps below show what turning it on would do."
                         if not enabled else
                         ("From a measured CT." if measured else "From a pseudo-CT generated from the UTE image."))
        state_ct = ct.state if enabled else rp.SKIP
        if measured:
            ct_file = ct.detail
            self.steps[0].set(state_ct, "Measured CT", f"{_esc(g('path.ct_pattern'))} in path.anat, in Hounsfield units"
                              + ("" if state_ct != rp.STOP else f" — {_esc(ct_file)}"), "ct")
            self.steps[1].set(state_ct if state_ct != rp.STOP else rp.SKIP, "Resample to the segmentation grid",
                              "through both images' world coordinates; the CT must be co-registered to the T1w", "grid")
            self.steps[2].setVisible(False)
            self.steps[3].set(state_ct if state_ct != rp.STOP else rp.SKIP, "Skull voxels",
                              "the skull layer of the SimNIBS segmentation, as for a pseudo-CT", "tissue_skull")
            for a in self.arrows[2:3]:
                a.setVisible(False)
            if self.ute_canvas is not None:
                self.ute_canvas.setVisible(False)
        else:
            gen = rp.RUN if state_ct == rp.RUN else (rp.REUSE if state_ct == rp.REUSE else rp.SKIP)
            reuse = " (already done: pseudoCT.nii.gz exists)" if gen == rp.REUSE else ""
            t2 = g("path.t2_pattern")
            self.steps[0].set(gen, "UTE / PETRA image",
                              (f"{_esc(t2)} in path.anat (path.t2_pattern)" if not is_empty(t2)
                               else "set path.t2_pattern to the UTE image") + reuse, "mri")
            self.steps[1].set(gen, "Bias correction and normalisation",
                              "N4 bias-field correction, then division by the soft-tissue peak, so soft tissue ≈ 1",
                              "histogram")
            self.steps[2].setVisible(True)
            self.arrows[2].setVisible(True)
            line = ute_line(m, self.window.prestus_root, dir_pct(m))
            eq = (f"HU = {line.slope:g} × UTE + {line.intercept:g}" if line.slope is not None
                  else "line not known yet")
            self.steps[2].set(gen, f"UTE → pseudo-HU: {_esc(line.mapping)}",
                              f"{eq} in the skull ({_esc(line.source)}); soft tissue and air get fixed values",
                              "line_map")
            self.steps[3].set(gen, "Pseudo-CT", "pseudoCT.nii.gz in Hounsfield units, in the pseudo-CT folder", "ct")
            self._draw_ute(line)
        mapped = enabled and state_ct != rp.STOP
        cal = g("pct.density_calibration")
        density_cal = (f" (calibration: {_esc(cal) if not is_empty(cal) else 'k-Plan default'})"
                       if g("pct.mapping_density") == "k-plan" else "")
        # one property per line; the k-Plan calibration stays on the density line it belongs to
        names = "<br>".join(f"{label}: {_esc(str(g(k) or 'none'))}{extra}" for k, label, extra in
                            (("pct.mapping_density", "Density", density_cal),
                             ("pct.mapping_soundspeed", "Sound speed", ""),
                             ("pct.mapping_attenuation", "Attenuation", "")))
        self.steps[4].set(rp.RUN if mapped else rp.SKIP, "HU → skull properties", names, "table")
        self._schedule_curves()

    def _draw_ute(self, line) -> None:
        if self.ute_canvas is None:
            return
        self.ute_canvas.setVisible(True)
        fig = self.ute_fig
        fig.clear()
        ax = fig.add_subplot(1, 1, 1)
        fs = style.PLOT_PT
        try:
            consts = read_script(self.window.prestus_root)
        except OSError:
            consts = {"soft": None, "air": None}
        if line.slope is not None:
            import numpy as np
            u = np.linspace(0, 1.2, 50)
            ax.plot(u, line.slope * u + line.intercept, lw=1.6, label="skull")
        else:
            ax.text(0.5, 0.5, "fitted per subject when the pseudo-CT is made", ha="center", va="center",
                    transform=ax.transAxes, fontsize=fs["small"])
        if consts.get("soft") is not None:
            ax.axhline(consts["soft"], ls="--", lw=1, color="0.5", label=f"soft tissue {consts['soft']:g} HU")
        ax.set_xlabel("normalised UTE (soft tissue ≈ 1)", fontsize=fs["small"])
        ax.set_ylabel("HU", fontsize=fs["small"])
        ax.tick_params(labelsize=fs["small"])
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=fs["small"], frameon=False)
        self.ute_canvas.draw()

    # ── HU → property curves (MATLAB) ───────────────────────────────────
    def _signature(self) -> tuple:
        g = self.model.get
        keys = ("pct.mapping_density", "pct.mapping_soundspeed", "pct.mapping_attenuation",
                "pct.density_calibration", "pct.skull_mapping")
        path = self.model.path
        mtime = path.stat().st_mtime_ns if path is not None and Path(path).is_file() else None
        return (str(path or ""), mtime, tuple(repr(g(k)) for k in keys), repr(g("transducer.freq_hz")))

    def _args(self) -> dict:
        g = self.model.get
        pct = {k: g(f"pct.{k}") for k in ("mapping_density", "mapping_soundspeed", "mapping_attenuation",
                                          "density_calibration", "skull_mapping")}
        pct = {k: ("" if v is None else v) for k, v in pct.items()}
        args = {"action": "curves", "config": str(self.model.path or ""), "pct": pct}
        freq = g("transducer.freq_hz")
        if isinstance(freq, (int, float)) and not isinstance(freq, bool):
            args["freq_hz"] = freq
        return args

    def _schedule_curves(self) -> None:
        sig = self._signature()
        if sig in _CURVES:
            self._draw_curves(_CURVES[sig])
            return
        self._draw_curves(None)
        if getattr(self.window.bridge(), "state", "") in ("ready", "busy"):
            self._timer.start()            # MATLAB is up: follow the edits

    def request_curves(self, force: bool = False) -> None:
        sig = self._signature()
        if sig in _CURVES or self._pending == sig:
            return
        if not force and getattr(self.window.bridge(), "state", "") not in ("ready", "busy"):
            return
        if self.window.matlab_path() is None:
            self.curves_note.setText("MATLAB not found; set its path on the Run tab.")
            return
        self._pending = sig
        self.curves_note.setText("⟳ Evaluating the mappings in MATLAB…")
        self.window.bridge().call("nexus_pct", self._args(),
                                  on_done=lambda out, s=sig: self._got_curves(s, out),
                                  on_error=lambda err, s=sig: self._curves_failed(s, err))

    def _got_curves(self, sig, out: dict) -> None:
        self._pending = None
        _CURVES[sig] = out
        if sig == self._signature():
            self._draw_curves(out)

    def _curves_failed(self, sig, err) -> None:
        self._pending = None
        c = style.active()
        self.curves_note.setText(f"<span style='color:{c['red']}'>Could not evaluate the mappings: "
                                 f"{_esc(str(err))}</span>")

    def _draw_curves(self, curves: dict | None) -> None:
        if self.prop_canvas is None:
            self.curves_note.setText("Install matplotlib to see the mapping curves.")
            return
        self.curves_btn.setVisible(curves is None)
        if curves is None:
            self.prop_canvas.setVisible(False)
            if self._pending is None:
                self.curves_note.setText("The curves are evaluated by the run's own MATLAB code.")
            return
        self.prop_canvas.setVisible(True)
        mp = curves.get("mapping", {}) or {}
        cal = mp.get("density_calibration")
        self.curves_note.setText("Evaluated in MATLAB (pct_mapping_curves) with these settings"
                                 + (f"; density table: {_esc(str(cal))}" if cal else "")
                                 + ". Dashed: the uniform skull value used without CT mapping.")
        fig = self.prop_fig
        fig.clear()
        fs = style.PLOT_PT
        hu = curves.get("hu", [])
        ref = (curves.get("reference", {}) or {}).get("skull", {}) or {}
        for i, (key, label, unit) in enumerate(PROPS):
            ax = fig.add_subplot(3, 1, i + 1)
            y = curves.get(key, [])
            ax.plot(hu, y, lw=1.6)
            if isinstance(ref.get(key), (int, float)):
                ax.axhline(ref[key], ls="--", lw=1, color="0.5")
            name = {"density": mp.get("density"), "sound_speed": mp.get("sound_speed"),
                    "alpha_coeff": mp.get("attenuation")}[key]
            ax.set_title(f"{label} ({name}) [{unit}]", fontsize=fs["small"], loc="left")
            ax.tick_params(labelsize=fs["small"])
            if i == 2:
                ax.set_xlabel("HU", fontsize=fs["small"])
            else:
                ax.tick_params(labelbottom=False)
        self.prop_canvas.draw()
        self._shown_sig = self._signature()


def _esc(text) -> str:
    import html
    return html.escape(str(text))
