"""Plots of a transducer calibration library entry, readable without MATLAB.

CalibrationCoverage draws what an entry (config/calibration/<serial>_<ds>.yaml) holds:
fit quality, the drive per depth and the per-ring phase corrections, against the depths this
study uses, so a user sees at a glance whether the study's focal depths fall inside the
calibrated range. R2Strip is the clickable per-depth R² bar strip of the Calibration sub-tab.

entry_arrays() and summary_text() are pure Python (no Qt) so they can be unit-tested.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

from PyQt6.QtCore import QEvent, QRectF, QSize, Qt, pyqtSignal
from PyQt6.QtGui import QColor, QFont, QPainter, QPen
from PyQt6.QtWidgets import (QButtonGroup, QHBoxLayout, QLabel, QSizePolicy, QStackedLayout,
                             QToolButton, QVBoxLayout, QWidget)

from . import style
from .export import PaintedExportable, attach_painted

try:
    from .mpl_theme import ThemedCanvas as FigureCanvasQTAgg   # drawn in the Nexus theme
    from matplotlib.figure import Figure
except ImportError:                                    # pragma: no cover
    FigureCanvasQTAgg = None

# R² thresholds for the green / amber / red colouring of fit quality. 0.9 is where a fitted
# axial profile is visually indistinguishable from the hydrophone at the focus; below 0.8 the
# focal position or width is typically off by more than a millimetre (shipped CTX global fits).
R2_GOOD = 0.9
R2_FAIR = 0.8
# Phase steps between neighbouring calibrated depths larger than this (after unwrapping) are
# marked: interpolation between those depths may then be poorly constrained.
PHASE_JUMP_DEG = 30.0

VIEWS = ("quality", "drive", "phases")
# Categorical ring colours (light, dark), fixed order; status colours stay reserved for R².
# Checked with the dataviz palette validator (CVD and lightness band, both themes).
RING_COLOURS = [("#1f6fb2", "#4b95d1"), ("#c75b12", "#c4773a"), ("#7a4fb0", "#9f7ad6"),
                ("#00927c", "#2ea893")]


def r2_colour_key(r2: float | None) -> str:
    """style.active() key for an R² value (thresholds R2_GOOD / R2_FAIR)."""
    if r2 is None or not math.isfinite(r2):
        return "ink4"
    return "green" if r2 >= R2_GOOD else "amber" if r2 >= R2_FAIR else "red"


# ── data (no Qt) ─────────────────────────────────────────────────────────

def _floats(v: Any) -> list[float]:
    if v is None:
        return []
    if not isinstance(v, (list, tuple)):
        v = [v]
    out = []
    for x in v:
        try:
            out.append(float(x))
        except (TypeError, ValueError):
            out.append(float("nan"))
    return out


def _num(v: Any) -> float | None:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def _rows(v: Any) -> list[list[float]]:
    """Depth × element table from YAML (list of lists; a flat list is one depth)."""
    if not isinstance(v, list) or not v:
        return []
    if not isinstance(v[0], list):
        return [_floats(v)]
    return [_floats(r) for r in v]


def _unwrap_deg(col: list[float]) -> list[float]:
    """numpy.unwrap in degrees, the first value wrapped to (-180, 180] so offsets near 360°
    read as small negative offsets."""
    import numpy as np
    a = np.asarray(col, dtype=float)
    ok = np.isfinite(a)
    if not ok.any():
        return list(a)
    out = a.copy()
    out[ok] = np.rad2deg(np.unwrap(np.deg2rad(a[ok])))
    first = out[ok][0]
    out[ok] -= 360.0 * np.round(first / 360.0)
    return [float(x) for x in out]


def entry_arrays(entry: dict | None) -> dict:
    """The plottable arrays of a parsed calibration entry; missing blocks give empty lists /
    None, so legacy, geometric and DS-agnostic entries degrade gracefully.

    Keys: fit_depths, fit_r2 (global_model per-depth R²), geometric_r2_median,
    in_sample_r2_median, in_sample_r2_min, flags, cal_depths (refinement depths, else the
    legacy global_model.depths_ep_mm, else fit depths), cal_range (lo, hi) or None, blend_mm,
    drive_depths, drive_pa, drive_isppa, drive_source ('refinement' | 'global' | ''),
    global_drive_depths, global_drive_pa, amp_source, phase_depths, phase_deg (unwrapped,
    one list per channel; identical channel columns -- sub-rings of one ring -- merged),
    channel_labels, phase_jumps [(depth index, channel index, step °)], calibration_type, updated.
    """
    entry = entry if isinstance(entry, dict) else {}
    meta = entry.get("meta") or {}
    g = entry.get("global_model") or {}
    r = entry.get("refinement") or {}
    q = meta.get("quality") if isinstance(meta.get("quality"), dict) else {}
    out: dict[str, Any] = {
        "calibration_type": str(meta.get("calibration_type") or ""),
        "updated": str(meta.get("last_updated") or ""),
        "fit_depths": _floats(g.get("fit_depths_ep_mm")),
        "fit_r2": _floats(g.get("fit_r2")),
        "geometric_r2_median": _num(q.get("geometric_r2_median")),
        "in_sample_r2_median": _num(q.get("in_sample_r2_median")),
        "in_sample_r2_min": _num(q.get("in_sample_r2_min")),
        "flags": [str(f) for f in (q.get("flags") or [])] if isinstance(q.get("flags"), list) else [],
        "blend_mm": _num(r.get("blend_mm")) or 5.0,
    }
    n = min(len(out["fit_depths"]), len(out["fit_r2"]))
    out["fit_depths"], out["fit_r2"] = out["fit_depths"][:n], out["fit_r2"][:n]

    cal = _floats(r.get("depths_ep_mm")) or _floats(g.get("depths_ep_mm")) or out["fit_depths"]
    out["cal_depths"] = cal
    fin = [d for d in cal if math.isfinite(d)]
    out["cal_range"] = (min(fin), max(fin)) if fin else None

    abd = g.get("amplitude_by_depth") or {}
    out["global_drive_depths"] = _floats(abd.get("depths_ep_mm"))
    out["global_drive_pa"] = _floats(abd.get("elem_p_amp"))
    rd, rp = _floats(r.get("depths_ep_mm")), _floats(r.get("elem_p_amp"))
    if rd and len(rp) == len(rd) and any(math.isfinite(x) for x in rp):
        out.update(drive_depths=rd, drive_pa=rp, drive_isppa=_num(r.get("isppa_wcm2")),
                   drive_source="refinement")
    elif out["global_drive_depths"] and len(out["global_drive_pa"]) == len(out["global_drive_depths"]):
        out.update(drive_depths=out["global_drive_depths"], drive_pa=out["global_drive_pa"],
                   drive_isppa=_num(abd.get("isppa_wcm2")), drive_source="global")
    else:
        out.update(drive_depths=[], drive_pa=[], drive_isppa=None, drive_source="")
    src = r.get("amp_source")
    out["amp_source"] = [str(s) for s in src] if isinstance(src, list) else \
        ([str(src)] * len(rd) if isinstance(src, str) else [])

    rows = _rows(r.get("elem_phase_correction_deg"))
    pdepths = rd if rows and len(rows) == len(rd) else []
    out["phase_depths"], out["phase_deg"], out["channel_labels"], out["phase_jumps"] = pdepths, [], [], []
    if pdepths:
        n_el = min(len(row) for row in rows)
        cols = [[row[e] for row in rows] for e in range(n_el)]
        merged, groups = [], []
        for e, col in enumerate(cols):          # sub-rings of a ring profile repeat their ring's value
            if merged and all((a == b) or (math.isnan(a) and math.isnan(b)) for a, b in zip(col, merged[-1])):
                groups[-1].append(e)
                continue
            merged.append(col)
            groups.append([e])
        out["phase_deg"] = [_unwrap_deg(c) for c in merged]
        out["channel_labels"] = [f"ring {i + 1}" for i in range(len(merged))]
        for c, col in enumerate(out["phase_deg"]):
            for j in range(1, len(col)):
                step = col[j] - col[j - 1]
                if math.isfinite(step) and abs(step) > PHASE_JUMP_DEG:
                    out["phase_jumps"].append((j, c, step))
    return out


def depth_status(arrays: dict, depth: float) -> str:
    """'within' the calibrated range, 'blend' (outside, within blend_mm: fades to the global
    model) or 'outside'; 'none' without a calibrated range. Mirrors equipment.Calibration.at_depth."""
    rng = arrays.get("cal_range")
    if not rng:
        return "none"
    lo, hi = rng
    if lo - 1e-6 <= depth <= hi + 1e-6:
        return "within"
    extra = lo - depth if depth < lo else depth - hi
    return "blend" if extra <= (arrays.get("blend_mm") or 5.0) else "outside"


STATUS_TEXT = {"within": "within range", "blend": "just outside (blended to the global model)",
               "outside": "outside the calibrated range", "none": "no calibrated range"}
STATUS_COLOUR = {"within": "green", "blend": "amber", "outside": "red", "none": "ink3"}


def summary_text(arrays: dict, study_depths: list[float] | None = None) -> str:
    """One-line caption, e.g. '12 depths, 39.6–78.8 mm · median R² 0.97 · study depth 50 mm:
    within range'."""
    parts = []
    rng = arrays.get("cal_range")
    n = len(arrays.get("cal_depths") or [])
    if rng:
        parts.append(f"{n} depth{'s' if n != 1 else ''}, {rng[0]:.1f}–{rng[1]:.1f} mm")
    med = arrays.get("in_sample_r2_median")
    if med is None:
        r2 = sorted(x for x in arrays.get("fit_r2") or [] if math.isfinite(x))
        med = r2[len(r2) // 2] if r2 else None
    if med is not None:
        parts.append(f"median R² {med:.2f}")
    for d in study_depths or []:
        parts.append(f"study depth {d:g} mm: {STATUS_TEXT[depth_status(arrays, d)]}")
    if arrays.get("flags"):
        parts.append("flags: " + ", ".join(arrays["flags"]))
    return " · ".join(parts) if parts else "no calibrated depths in this entry"


def load_entry(library, serial: str, ds: str = "", library_path: str = "") -> tuple[dict | None, Path | None]:
    """The parsed entry PRESTUS would use for serial / driving system (library_path first,
    then config/calibration; DS-agnostic fallback), or (None, None)."""
    if library is None or not serial:
        return None, None
    from .. import yamlio
    cal = library.calibration(serial, ds or None, library_path or None)
    if cal.path is None:
        return None, None
    try:
        data = yamlio.load_file_cached(cal.path)    # the project volume can be slow: read once
        return (data if isinstance(data, dict) else None), cal.path
    except Exception:
        return None, cal.path




def hydrophone_csv(entry: dict | None, prestus_root: Path | None) -> Path | None:
    """meta.char_data resolved as build_transducer_library.m does (CHAR_DATA_ROOT =
    <PRESTUS>/local/profiles), if that file exists."""
    rel = ((entry or {}).get("meta") or {}).get("char_data")
    if not rel or prestus_root is None:
        return None
    p = Path(str(rel))
    if not p.is_absolute():
        p = Path(prestus_root) / "local" / "profiles" / p
    return p if p.is_file() else None


def csv_focal_depths(path: str | Path) -> list[float]:
    """Focal depths from the header of a hydrophone axial CSV (row 2: 'Distance WRT Exit
    plane (mm)', then one focal depth per column), for listing them before MATLAB has read
    the file; [] if the layout differs (MATLAB's reader then decides)."""
    import csv
    try:
        with open(path, newline="", encoding="utf-8-sig") as f:
            rows = [row for _, row in zip(range(2), csv.reader(f))]
    except (OSError, UnicodeDecodeError, csv.Error):
        return []
    if len(rows) < 2 or "focal" not in ",".join(rows[0]).lower():
        return []
    out = []
    for x in rows[1][1:]:
        try:
            out.append(float(x))
        except ValueError:
            continue
    return out


# ── Qt widgets ───────────────────────────────────────────────────────────

def _segmented_css(c: dict) -> str:
    return (f"QToolButton {{ border: 1px solid {c['line2']}; padding: 2px 10px; border-radius: 0; "
            f"background: {c['surface']}; color: {c['ink2']}; }}"
            f"QToolButton:checked {{ background: {c['brand']}; color: white; border-color: {c['brand']}; }}"
            f"QToolButton:hover:!checked {{ background: {c['surface3']}; }}")


class CalibrationCoverage(QWidget):
    """Coverage of one library entry: a "Quality · Drive · Phases" toggle over one compact plot,
    a caption (summary_text) and an empty state when there is no entry.

    Use: w = CalibrationCoverage(); w.set_entry(entry_dict_or_None, study_depths=[47.5],
    name="CTX_250_001 / SC_105_010"). set_device(library, serial, ds, library_path,
    study_depths) loads the entry itself. Redraws on theme (palette) changes."""

    def __init__(self, parent: QWidget | None = None, height: int = 230):
        super().__init__(parent)
        self.arrays: dict | None = None
        self.study_depths: list[float] = []
        self.name = ""
        self.empty_message = ""
        self.view = "quality"
        v = QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(4)
        bar = QHBoxLayout()
        bar.setSpacing(0)
        self.title = QLabel()
        self.title.setObjectName("CardSub")
        bar.addWidget(self.title, 1)
        self.buttons = QButtonGroup(self)
        self.buttons.setExclusive(True)
        self._seg: list[QToolButton] = []
        for key, text, tip in (("quality", "Quality", "Fit R² per calibrated depth"),
                               ("drive", "Drive", "Drive amplitude (elem_p_amp) per depth"),
                               ("phases", "Phases", "Per-ring phase correction per depth (unwrapped)")):
            b = QToolButton(text=text, checkable=True)
            b.setToolTip(tip)
            b.setProperty("view", key)
            b.clicked.connect(lambda _c, k=key: self.set_view(k))
            self.buttons.addButton(b)
            bar.addWidget(b)
            self._seg.append(b)
        self._seg[0].setChecked(True)
        v.addLayout(bar)

        self.stack = QStackedLayout()
        holder = QWidget()
        holder.setLayout(self.stack)
        holder.setMinimumHeight(height)
        self._holder, self._height = holder, height
        holder.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.empty = QLabel()
        self.empty.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.empty.setWordWrap(True)
        self.empty.setObjectName("CardSub")
        self.stack.addWidget(self.empty)
        if FigureCanvasQTAgg is not None:
            self.fig = Figure(figsize=(4, 2.4), layout="constrained")
            self.ax = self.fig.add_subplot(111)
            self.canvas = FigureCanvasQTAgg(self.fig, export_name="calibration_coverage",
                                            filename_hint=lambda: f"calibration_{self.view}")
            self.canvas.setMinimumHeight(140)
            self.stack.addWidget(self.canvas)
        else:                                               # pragma: no cover
            self.fig = self.ax = self.canvas = None
        v.addWidget(holder, 1)
        self.caption = QLabel()
        self.caption.setObjectName("CardSub")
        self.caption.setWordWrap(True)
        self.caption.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        v.addWidget(self.caption)
        self._restyle()
        self.set_entry(None)

    # public API
    def set_entry(self, entry: dict | None, study_depths: list[float] | None = None, name: str = "",
                  empty_message: str = "") -> None:
        """Show entry (parsed YAML) against study_depths (mm from the exit plane). With
        entry None, show empty_message (default: no calibration entry; the geometric model
        is used)."""
        self.arrays = entry_arrays(entry) if isinstance(entry, dict) else None
        self.study_depths = [float(d) for d in (study_depths or []) if d is not None and math.isfinite(float(d))]
        self.name = name
        self.empty_message = empty_message
        self.redraw()

    def set_device(self, library, serial: str, ds: str = "", library_path: str = "",
                   study_depths: list[float] | None = None) -> Path | None:
        """Load and show the entry PRESTUS would use for serial / ds; returns its path."""
        entry, path = load_entry(library, serial, ds, library_path)
        name = f"{serial}{' with ' + ds if ds else ''}"
        msg = "" if entry is not None else \
            (f"No calibration entry for {serial}{' with ' + ds if ds else ''}; the geometric model is used."
             if serial else "No transducer selected.")
        self.set_entry(entry, study_depths, name, msg)
        return path

    def set_view(self, view: str) -> None:
        if view not in VIEWS:
            return
        self.view = view
        for b in self._seg:
            b.setChecked(b.property("view") == view)
        self.redraw()

    def summary_text(self) -> str:
        return summary_text(self.arrays, self.study_depths) if self.arrays else self.empty.text()

    def sizeHint(self) -> QSize:
        return QSize(420, 300)

    # drawing
    def changeEvent(self, ev) -> None:
        if ev.type() in (QEvent.Type.PaletteChange, QEvent.Type.StyleChange):
            self._restyle()
            self.redraw()
        super().changeEvent(ev)

    def _restyle(self) -> None:
        css = _segmented_css(style.active())
        for i, b in enumerate(self._seg):
            r = "border-top-left-radius: 5px; border-bottom-left-radius: 5px;" if i == 0 else \
                "border-top-right-radius: 5px; border-bottom-right-radius: 5px;" if i == len(self._seg) - 1 else ""
            b.setStyleSheet(css.replace("border-radius: 0;", "border-radius: 0; " + r))

    def _has_data(self) -> bool:
        a = self.arrays or {}
        return bool(a.get("fit_depths") or a.get("drive_depths") or a.get("phase_depths") or a.get("cal_range"))

    def redraw(self) -> None:
        a = self.arrays
        self.title.setText(f"<b>{self.name}</b>" if self.name else "")
        for b in self._seg:
            b.setVisible(a is not None and self._has_data())
        if a is None or not self._has_data() or self.canvas is None:
            msg = self.empty_message or "No calibration entry; the geometric model is used."
            if a is not None and not self._has_data():
                msg = (f"{self.name or 'This entry'} is a {a.get('calibration_type') or 'geometric'} entry "
                       "without calibrated depths; the geometric model is used.")
            self.empty.setText(msg)
            self.stack.setCurrentWidget(self.empty)
            self._holder.setMinimumHeight(48)               # a one-line message, not empty axes
            self._holder.setMaximumHeight(96)
            self.caption.setText("")
            return
        self.stack.setCurrentWidget(self.canvas)
        self._holder.setMinimumHeight(self._height)
        self._holder.setMaximumHeight(16777215)
        c = style.active()
        self._caption(c)
        fig, ax = self.fig, self.ax
        fig.set_facecolor(c["surface"])
        ax.clear()
        ax.set_facecolor(c["surface"])
        getattr(self, f"_draw_{self.view}")(ax, a, c)
        self._decorate(ax, a, c)
        self.canvas.draw_idle()

    def _caption(self, c: dict) -> None:
        a = self.arrays or {}
        base = summary_text({**a, "flags": a.get("flags")}, [])
        parts = [base] if base else []
        for d in self.study_depths:
            st = depth_status(a, d)
            parts.append(f"study depth {d:g} mm: <span style='color:{c[STATUS_COLOUR[st]]}'>"
                         f"{STATUS_TEXT[st]}</span>")
        self.caption.setText(" · ".join(parts))

    def _decorate(self, ax, a: dict, c: dict) -> None:
        """Calibrated range, study depths, axes styling (shared by all views)."""
        rng = a.get("cal_range")
        if rng:
            ax.axvspan(rng[0], rng[1], color=c["brand"], alpha=0.08, lw=0, zorder=0)
        xs = [x for x in (a.get("cal_depths") or []) if math.isfinite(x)] + self.study_depths
        if xs:
            lo, hi = min(xs), max(xs)
            pad = max(2.0, 0.06 * (hi - lo))
            ax.set_xlim(lo - pad, hi + pad)
        for d in self.study_depths:
            col = c[STATUS_COLOUR[depth_status(a, d)]]
            ax.axvline(d, color=col, lw=1.4, ls="-", zorder=3)
            ax.annotate(f"{d:g}", (d, 1.0), xycoords=("data", "axes fraction"), xytext=(2, -2),
                        textcoords="offset points", ha="left", va="top", fontsize=style.PLOT_PT["small"], color=col,
                        fontweight="bold")
        ax.set_xlabel("Distance from exit plane [mm]", fontsize=style.PLOT_PT["normal"], color=c["ink2"])
        ax.tick_params(labelsize=style.PLOT_PT["small"], colors=c["ink3"], length=3)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(c["line2"])
        ax.grid(True, color=c["line"], lw=0.6)
        ax.set_axisbelow(True)
        leg = ax.get_legend_handles_labels()
        if leg[0]:
            lg = ax.legend(fontsize=style.PLOT_PT["small"], loc="best", frameon=True, framealpha=0.85, borderpad=0.3,
                           handlelength=1.6, labelcolor=c["ink2"])
            lg.get_frame().set_facecolor(c["surface"])
            lg.get_frame().set_edgecolor(c["line"])

    def _draw_quality(self, ax, a: dict, c: dict) -> None:
        d, r2 = a["fit_depths"], a["fit_r2"]
        if d:
            ax.plot(d, r2, "-", color=c["ink3"], lw=1.2, zorder=2)
            ax.scatter(d, r2, s=26, c=[c[r2_colour_key(x)] for x in r2], edgecolors=c["surface"],
                       linewidths=1.0, zorder=4, label="Global model R²")
        geo = a.get("geometric_r2_median")
        if geo is not None:
            ax.axhline(geo, color=c["ink3"], ls="--", lw=1.0, label=f"Geometric, median {geo:.2f}")
        ins = a.get("in_sample_r2_median")
        if ins is not None:
            ax.axhline(ins, color=c["accent"], ls=":", lw=1.3, label=f"Per-depth fit, median {ins:.2f}")
        vals = [x for x in r2 + [geo, ins] if x is not None and math.isfinite(x)]
        ax.set_ylim(max(-0.2, min(vals + [R2_FAIR]) - 0.05) if vals else 0, 1.02)
        ax.set_ylabel("R² (profile shape)", fontsize=style.PLOT_PT["normal"], color=c["ink2"])
        if not d:
            ax.text(0.5, 0.5, "No per-depth R² stored in this entry", transform=ax.transAxes,
                    ha="center", va="center", fontsize=style.PLOT_PT["normal"], color=c["ink3"])

    def _draw_drive(self, ax, a: dict, c: dict) -> None:
        d, p = a["drive_depths"], a["drive_pa"]
        if not d:
            ax.text(0.5, 0.5, "No drive per depth in this entry", transform=ax.transAxes,
                    ha="center", va="center", fontsize=style.PLOT_PT["normal"], color=c["ink3"])
            return
        kpa = [x / 1000 for x in p]
        isppa = a.get("drive_isppa")
        label = ("Per-depth fit" if a["drive_source"] == "refinement" else "Global model")
        ax.plot(d, kpa, "-o", color=c["accent"], ms=4, lw=1.6, mec=c["surface"], label=label, zorder=3)
        kw = [j for j, s in enumerate(a.get("amp_source") or []) if s == "kwave" and j < len(d)]
        if kw:
            ax.plot([d[j] for j in kw], [kpa[j] for j in kw], "s", ms=6, mfc="none", mec=c["green"],
                    label="k-Wave calibrated", zorder=4)
        if a["drive_source"] == "refinement" and a["global_drive_depths"]:
            ax.plot(a["global_drive_depths"], [x / 1000 for x in a["global_drive_pa"]], "--",
                    color=c["ink3"], lw=1.0, label="Global model", zorder=2)
        ax.set_ylabel(f"Drive{f' for {isppa:g} W/cm²' if isppa else ''} [kPa]", fontsize=style.PLOT_PT["normal"], color=c["ink2"])

    def _draw_phases(self, ax, a: dict, c: dict) -> None:
        d, cols = a["phase_depths"], a["phase_deg"]
        if not d:
            ax.text(0.5, 0.5, "No per-depth phase corrections in this entry", transform=ax.transAxes,
                    ha="center", va="center", fontsize=style.PLOT_PT["normal"], color=c["ink3"])
            return
        n = len(cols)
        if n <= len(RING_COLOURS):                  # identity: fixed categorical order
            colours = [RING_COLOURS[i][1 if style.is_dark() else 0] for i in range(n)]
        else:                                       # many rings: ordered inner → outer
            from matplotlib.colors import LinearSegmentedColormap
            cmap = LinearSegmentedColormap.from_list("rings", [c["brand"], c["accent"], c["ink3"]])
            colours = [cmap(i / max(1, n - 1)) for i in range(n)]
        for i, col in enumerate(cols):
            ax.plot(d, col, "-", color=colours[i], lw=1.4, marker="o", ms=2.5,
                    label=a["channel_labels"][i] if n <= 4 else None)
        if a["phase_jumps"]:
            xs = [d[j] for j, _ch, _s in a["phase_jumps"]]
            ys = [cols[ch][j] for j, ch, _s in a["phase_jumps"]]
            ax.scatter(xs, ys, s=40, marker="o", facecolors="none", edgecolors=c["red"], linewidths=1.4,
                       zorder=5, label=f"step > {PHASE_JUMP_DEG:g}°")
        ax.axhline(0, color=c["line2"], lw=0.8)
        ax.set_ylabel(f"Phase offset{'' if n <= 4 else f', {n} rings'} [°]", fontsize=style.PLOT_PT["normal"], color=c["ink2"])


class R2Strip(PaintedExportable, QWidget):
    """Clickable per-depth R² bars (green ≥ R2_GOOD, amber ≥ R2_FAIR, red below, grey: no R²).
    Emits selected(depth) on a click on a bar and selected(nan) on "All"."""

    selected = pyqtSignal(float)
    export_name = "calibration_r2"
    BUTTON_W = 30                                   # right margin kept clear for the export button

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.depths: list[float] = []
        self.r2: list[float] = []
        self.current: float | None = None          # None: all depths
        self.allow_all = True
        self.caption = ""
        self.setMinimumHeight(64)
        self.setMouseTracking(True)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.export_button = attach_painted(self)

    def sizeHint(self) -> QSize:
        return QSize(600, 70)

    def set_data(self, depths: list[float], r2: list[float], current: float | None = None,
                 allow_all: bool = True, caption: str = "") -> None:
        self.depths = [float(x) for x in depths]
        self.r2 = [float(x) for x in r2] + [float("nan")] * max(0, len(depths) - len(r2))
        self.current = current
        self.allow_all = allow_all
        self.caption = caption
        self.setToolTip(f"R² of {caption} against the measured profile per depth; click a bar to show "
                        f"that depth. Green ≥ {R2_GOOD:g}, amber ≥ {R2_FAIR:g}, red below." if caption else "")
        self.update()

    def select(self, depth: float | None) -> None:
        self.current = depth
        self.update()

    def _geometry(self, w: float | None = None, h: float | None = None) -> tuple[QRectF | None, list[QRectF]]:
        w = (self.width() if w is None else w) - self.BUTTON_W
        h = self.height() if h is None else h
        x0 = 4.0
        all_rect = None
        if self.allow_all:
            all_rect = QRectF(x0, 14, 34, h - 30)
            x0 += 40
        n = max(1, len(self.depths))
        slot = (w - x0 - 4) / n
        bw = max(4.0, min(36.0, slot - 3))
        return all_rect, [QRectF(x0 + i * slot + (slot - bw) / 2, 14, bw, h - 30) for i in range(len(self.depths))]

    def paint(self, p: QPainter, rect_all: QRectF, c: dict) -> None:
        """The strip in palette ``c`` (style.active() on screen, export.py for files)."""
        p.fillRect(rect_all, QColor(c["surface"]))
        f = QFont(self.font())
        f.setPointSizeF(max(7.0, f.pointSizeF() * 0.8))
        p.setFont(f)
        all_rect, rects = self._geometry(rect_all.width(), rect_all.height())
        if all_rect is not None:
            on = self.current is None
            p.setPen(QPen(QColor(c["brand"] if on else c["line2"]), 1.5 if on else 1))
            p.setBrush(QColor(c["brand_soft"] if "brand_soft" in c else c["surface2"]) if on else QColor(c["surface"]))
            p.drawRoundedRect(all_rect, 4, 4)
            p.setPen(QColor(c["ink"] if on else c["ink3"]))
            p.drawText(all_rect, Qt.AlignmentFlag.AlignCenter, "All")
        show_values = bool(rects) and rects[0].width() >= 24
        for d, r2, rect in zip(self.depths, self.r2, rects):
            ok = math.isfinite(r2)
            frac = min(1.0, max(0.03, r2)) if ok else 0.03
            bar = QRectF(rect.left(), rect.bottom() - rect.height() * frac, rect.width(), rect.height() * frac)
            sel = self.current is not None and abs(self.current - d) < 0.01
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(QColor(c["surface2"]))
            p.drawRoundedRect(rect, 3, 3)
            p.setBrush(QColor(c[r2_colour_key(r2 if ok else None)]))
            p.drawRoundedRect(bar, 3, 3)
            if sel:
                p.setBrush(Qt.BrushStyle.NoBrush)
                p.setPen(QPen(QColor(c["ink"]), 2))
                p.drawRoundedRect(rect.adjusted(-2, -2, 2, 2), 4, 4)
            p.setPen(QColor(c["ink"] if sel else c["ink3"]))
            p.drawText(QRectF(rect.center().x() - 30, rect.bottom() + 1, 60, 14),
                       Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop, f"{d:.4g}")
            if show_values or sel:
                p.setPen(QColor(c["ink2"]))
                p.drawText(QRectF(rect.center().x() - 30, 0, 60, 13),
                           Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignBottom,
                           f"{r2:.2f}" if ok else "—")
        if not rects:
            p.setPen(QColor(c["ink3"]))
            p.drawText(rect_all, Qt.AlignmentFlag.AlignCenter, "No depths yet")
        self.export_button.set_available(bool(rects))

    def mousePressEvent(self, ev) -> None:
        pos = ev.position()
        all_rect, rects = self._geometry()
        if all_rect is not None and all_rect.adjusted(-2, -14, 2, 16).contains(pos):
            self.select(None)
            self.selected.emit(float("nan"))
            return
        for d, rect in zip(self.depths, rects):
            if rect.adjusted(-2, -14, 2, 16).contains(pos):
                self.select(d)
                self.selected.emit(d)
                return
