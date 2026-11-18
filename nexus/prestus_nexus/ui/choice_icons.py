"""Schematic icons for drop-down choices (overlay 'choice_icons') and the larger figures
shown under a selected choice (overlay 'choice_figures').

Glyphs are drawn on the 24 x 24 grid of icons.py; figures on a 96 x 56 grid. Bodies use
three colour roles, written as classes and turned into SVG attributes by `qt()`:

    a   stroke in the accent colour (what the choice changes)
    f   filled in the ink colour        af  filled in the accent colour
    faint  ink stroke at reduced opacity     soft  accent fill at reduced opacity

Only primitives QSvgRenderer handles are used (path, circle, rect, ellipse, dash arrays,
opacity); figure captions are Qt labels, not SVG text, so they use the platform font.
"""

from __future__ import annotations

import math
import re

# ── Colour roles → SVG attributes ─────────────────────────────────────

_ROLE = {
    "a": 'stroke="accentColor"',
    "f": 'fill="currentColor" stroke="none"',
    "af": 'fill="accentColor" stroke="none"',
    "faint": 'stroke-opacity="0.45"',
    "soft": 'fill-opacity="0.35"',
}


def qt(body: str) -> str:
    """Replace class="…" with the attributes of its roles (placeholders currentColor and
    accentColor are substituted by icons.svg)."""
    return re.sub(r'class="([^"]*)"', lambda m: " ".join(_ROLE[r] for r in m.group(1).split()), body)


# ── Geometry helpers ──────────────────────────────────────────────────

def _f(x: float) -> str:
    return f"{x:.2f}".rstrip("0").rstrip(".")


def _poly(pts) -> str:
    return "M" + " L".join(f"{_f(x)} {_f(y)}" for x, y in pts)


def _wave(x0, x1, y, amp, cycles, phase=0.0, n=64):
    return [(x0 + (x1 - x0) * i / n, y - amp * math.sin(2 * math.pi * cycles * i / n + phase))
            for i in range(n + 1)]


def _quarter_arc(cx, cy, r, side):
    """Wavefront arc from straight up to horizontal, towards the other source."""
    if side == "L":
        return f"M{_f(cx)} {_f(cy - r)} A{_f(r)} {_f(r)} 0 0 1 {_f(cx + r)} {_f(cy)}"
    return f"M{_f(cx)} {_f(cy - r)} A{_f(r)} {_f(r)} 0 0 0 {_f(cx - r)} {_f(cy)}"


def _crossings(c1, c2, r1, r2):
    (x1, y1), (x2, y2) = c1, c2
    d = math.hypot(x2 - x1, y2 - y1)
    if d > r1 + r2 or d < abs(r1 - r2):
        return []
    a = (r1 ** 2 - r2 ** 2 + d ** 2) / (2 * d)
    h = math.sqrt(max(r1 ** 2 - a ** 2, 0))
    xm, ym = x1 + a * (x2 - x1) / d, y1 + a * (y2 - y1) / d
    return [(xm + h * (y2 - y1) / d, ym - h * (x2 - x1) / d),
            (xm - h * (y2 - y1) / d, ym + h * (x2 - x1) / d)]


def _arrow(x0, y0, x1, y1, cls="", head=2.6):
    ang = math.atan2(y1 - y0, x1 - x0)
    h1 = (x1 - head * math.cos(ang - 0.45), y1 - head * math.sin(ang - 0.45))
    h2 = (x1 - head * math.cos(ang + 0.45), y1 - head * math.sin(ang + 0.45))
    c = f' class="{cls}"' if cls else ""
    return (f'<path{c} d="M{_f(x0)} {_f(y0)} L{_f(x1)} {_f(y1)} '
            f'M{_f(h1[0])} {_f(h1[1])} L{_f(x1)} {_f(y1)} L{_f(h2[0])} {_f(h2[1])}"/>')


# ── Transducer coupling ───────────────────────────────────────────────
# coherent: pressure fields superposed in one run (interference);
# async: one run per transducer, intensities summed (config_default.yaml).

def _traces_glyph(coherent: bool) -> str:
    """Drive of transducer 1 (top) and 2 (bottom): one shared window vs. staggered bursts."""
    if coherent:
        return ('<path class="a" d="M6 3v18M18 3v18" stroke-dasharray="1.6 2.2" stroke-width="1.3"/>'
                '<path d="M2 8H6 q1.5 -3.4 3 0 t3 0 t3 0 t3 0 H22"/>'
                '<path d="M2 16H6 q1.5 -3.4 3 0 t3 0 t3 0 t3 0 H22"/>')
    return ('<path d="M2 8H4 q1.5 -3.4 3 0 t3 0 H22"/>'
            '<path d="M2 16H14 q1.5 -3.4 3 0 t3 0 H22"/>')


def _traces_figure(coherent: bool) -> str:
    b = ""
    if coherent:
        phi = math.pi / 2       # coherent = fixed phase relation, not equal phase
        b += '<path class="a" d="M34 3V36M66 3V36" stroke-dasharray="2 2.5" stroke-width="1"/>'
        for y, ph in ((12, 0.0), (28, phi)):
            b += f'<path d="{_poly([(8, y), (34, y)] + _wave(34, 66, y, 4, 4, ph) + [(66, y), (92, y)])}"/>'
        s = [(x, 46 - 4 * (math.sin(2 * math.pi * 4 * (x - 34) / 32) + math.sin(2 * math.pi * 4 * (x - 34) / 32 + phi)))
             for x in [34 + 32 * i / 96 for i in range(97)]]
        b += f'<path class="a" d="{_poly([(8, 46), (34, 46)] + s + [(66, 46), (92, 46)])}"/>'
    else:
        for y, x0 in ((12, 14), (28, 54)):
            b += f'<path d="{_poly([(8, y), (x0, y)] + _wave(x0, x0 + 28, y, 4, 3.5) + [(x0 + 28, y), (92, y)])}"/>'
        b += '<path class="a" d="M8 46H14V40H42V46H54V40H82V46H92"/>'
    return b


def _wavefront_figure(coherent: bool) -> str:
    s1, s2, radii = (12, 52), (84, 52), (12, 24, 36, 48)
    dash = "" if coherent else ' stroke-dasharray="3 3"'
    b = "".join(f'<path d="{_quarter_arc(*s1, r, "L")}"/><path d="{_quarter_arc(*s2, r, "R")}"{dash}/>'
                for r in radii)
    b += '<circle class="f" cx="12" cy="52" r="2.4"/><circle class="f" cx="84" cy="52" r="2.4"/>'
    if coherent:     # antinodes: where crests of both sources meet
        for r1 in radii:
            for r2 in radii:
                for x, y in _crossings(s1, s2, r1, r2):
                    if 3 < y < 51 and 3 < x < 93:
                        b += f'<circle class="af" cx="{_f(x)}" cy="{_f(y)}" r="2"/>'
    return b


# ── Grid orientation (grid.mode, grid.async_layout) ──────────
# preproc_head.m: focal_axis rotates the mean beam axis onto grid z, transducer_axis
# transducer(1)'s beam, ras_plus does not rotate. The grid box is drawn upright (z up);
# the accent marks what lies on z. Transducers sit above their focus.

def _bowl_beam(fx, fy, angle_deg, length=30.0, half=6.0, cls=""):
    """A bowl (thick arc, concave towards the focus) `length` from the focus at
    (fx, fy), tilted `angle_deg` from straight up, and its beam line to the focus."""
    a = math.radians(angle_deg)
    d = (-math.sin(a), math.cos(a))               # unit vector apex -> focus (SVG y down)
    p = (d[1], -d[0])
    ax, ay = fx - length * d[0], fy - length * d[1]
    e1 = (ax + 1.2 * d[0] + half * p[0], ay + 1.2 * d[1] + half * p[1])
    e2 = (ax + 1.2 * d[0] - half * p[0], ay + 1.2 * d[1] - half * p[1])
    ctl = (ax - 1.2 * d[0], ay - 1.2 * d[1])
    c = f' class="{cls}"' if cls else ""
    return (f'<path{c} d="M{_f(e1[0])} {_f(e1[1])} Q{_f(ctl[0])} {_f(ctl[1])} {_f(e2[0])} {_f(e2[1])}" '
            f'stroke-width="2.2"/>'
            f'<path{c} d="M{_f(ax)} {_f(ay)} L{_f(fx)} {_f(fy)}"/>')


def _rot_rect(cx, cy, w, h, angle_deg, cls=""):
    """Rectangle w x h centred at (cx, cy), its long side tilted angle_deg from vertical."""
    a = math.radians(angle_deg)
    u, v = (math.cos(a), math.sin(a)), (-math.sin(a), math.cos(a))
    pts = [(cx + sx * w / 2 * u[0] + sy * h / 2 * v[0], cy + sx * w / 2 * u[1] + sy * h / 2 * v[1])
           for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
    c = f' class="{cls}"' if cls else ""
    return f'<path{c} d="{_poly(pts)} Z"/>'


_GRID_BOX = '<rect class="faint" x="26" y="3" width="44" height="50" rx="1.5"/>'
_Z_AXIS = '<path class="a" d="M48 2V54" stroke-dasharray="2 2.2"/>'


def _gridmode_figure(mode: str) -> str:
    focus = '<circle class="af" cx="48" cy="40" r="1.8"/>'
    if mode == "focal_axis":
        return (_GRID_BOX + _Z_AXIS + _bowl_beam(48, 40, -28) + _bowl_beam(48, 40, 28) + focus)
    if mode == "transducer_axis":
        return (_GRID_BOX + _Z_AXIS + _bowl_beam(48, 40, 0, cls="a") + _bowl_beam(48, 40, 45, length=26) + focus)
    # ras_plus: no rotation; R (right) and S (up) axes of the T1 in the corner
    return (_GRID_BOX + _bowl_beam(48, 40, -34, length=27) + _bowl_beam(48, 40, 16) + focus
            + _arrow(8, 50, 8, 36, "a") + _arrow(8, 50, 22, 50, "a")
            + '<circle class="af" cx="8" cy="50" r="1.2"/>')


def _asyncgrid_figure(per_transducer: bool) -> str:
    f1, f2 = (36, 40), (60, 40)
    if not per_transducer:          # one grid, oriented on the mean axis, holding both beams
        return ('<rect class="faint" x="14" y="3" width="68" height="50" rx="1.5"/>'
                + _bowl_beam(*f1, -24, length=28) + _bowl_beam(*f2, 24, length=28)
                + f'<circle class="af" cx="{f1[0]}" cy="{f1[1]}" r="1.8"/>'
                  f'<circle class="af" cx="{f2[0]}" cy="{f2[1]}" r="1.8"/>')
    # one grid per transducer, each aligned with its own beam; grid 1 (accent) is the one
    # the fields are resampled onto (combine_async_intensity.m)
    f1, f2 = (30, 40), (66, 40)
    b = ""
    for (fx, fy), ang, box in ((f1, -24, "a"), (f2, 24, "faint")):
        a = math.radians(ang)
        cx, cy = fx - 13 * -math.sin(a), fy - 13 * math.cos(a)
        b += _rot_rect(cx, cy, 20, 46, ang, box) + _bowl_beam(fx, fy, ang, length=28)
        b += f'<circle class="af" cx="{fx}" cy="{fy}" r="1.8"/>'
    return b


# ── Glyphs (24 x 24) ──────────────────────────────────────────────────

_SPIRAL = "".join(  # golden-angle spiral, 21 points, as convert_to_element_pos.m 'fibonacci'
    f'<circle class="f" cx="{_f(12 + 8.5 * math.sqrt((i - .5) / 21) * math.cos(math.pi * (3 - math.sqrt(5)) * (i - 1)))}" '
    f'cy="{_f(12 + 8.5 * math.sqrt((i - .5) / 21) * math.sin(math.pi * (3 - math.sqrt(5)) * (i - 1)))}" r="0.95"/>'
    for i in range(1, 22))
_GRID4 = "".join(f'<circle class="f" cx="{x}" cy="{y}" r="1.1"/>'
                 for x in (6.75, 10.25, 13.75, 17.25) for y in (6.75, 10.25, 13.75, 17.25))
_BOWL = '<path d="M5 6Q12 1.5 19 6" stroke-width="2.2"/>'     # transducer: top, concave side to the focus below
_PROFILE = '<path d="M3 19C6.5 19 7.5 5 10.5 5C14 5 16 17 21 19"/>'
_SCALP = '<path d="M2 19Q12 11 22 19"/>'
_BOX = '<rect class="faint" x="3" y="3" width="18" height="18" rx="1"/>'
_FOCUS_ARROWS = _arrow(12, 14.5, 12, 9, "a") + _arrow(12, 14.5, 12, 20, "a")
_SMALL_BOWL_BEAM = '<path d="M9 7.5Q12 5 15 7.5" stroke-width="2.2"/><path d="M12 8.5V17"/>'
_WATER_LINE = "q2.25-2.4 4.5 0t4.5 0t4.5 0t4.5 0"

GLYPHS: dict[str, str] = {
    # simulation.transducer_coupling
    "coupling_coherent": _traces_glyph(True),
    "coupling_async": _traces_glyph(False),
    # simulation.medium: curved layers = segmented head; flat layers in a box = test object
    "medium_water": f'<path d="M3 7{_WATER_LINE}M3 12{_WATER_LINE}M3 17{_WATER_LINE}"/>',
    "medium_layered": '<path d="M3 9Q12 2 21 9"/><path class="a" d="M3 14Q12 7 21 14"/><path d="M3 19Q12 12 21 19"/>',
    "medium_phantom": '<rect x="3.5" y="3.5" width="17" height="17" rx="1.5"/><path class="a" d="M3.5 10h17M3.5 14h17"/>',
    # transducer.type (aperture, face view)
    "tx_annular": '<circle cx="12" cy="12" r="9"/><circle cx="12" cy="12" r="6"/><circle cx="12" cy="12" r="3"/>',
    "tx_matrix": '<circle cx="12" cy="12" r="9"/>' + "".join(
        f'<circle class="f" cx="{x}" cy="{y}" r="1.2"/>' for x in (8, 12, 16) for y in (8, 12, 16)),
    # transducer.matrix.steering: focus moves along the axis vs. in all directions
    "steer_1d": _BOWL + '<path class="faint" d="M12 7V21" stroke-dasharray="1.4 2"/>' + _FOCUS_ARROWS
                + '<circle class="af" cx="12" cy="14.5" r="1.8"/>',
    "steer_3d": _BOWL + _FOCUS_ARROWS + _arrow(12, 14.5, 6.5, 14.5, "a") + _arrow(12, 14.5, 17.5, 14.5, "a")
                + '<circle class="af" cx="12" cy="14.5" r="1.8"/>',
    # transducer.matrix.elem_shape (one element, face view)
    "elem_rect": '<rect x="6.5" y="6.5" width="11" height="11" rx="0.8"/>',
    "elem_disc": '<circle cx="12" cy="12" r="6"/>',
    "elem_bowl": '<circle cx="12" cy="12" r="6"/><path class="a" d="M8 10.5Q12 15.5 16 10.5"/>',
    # …grid_shape.type: fibonacci is laid out flat then curved; fermat is k-Wave's sampling
    # on the bowl surface, so it is drawn side-on
    "layout_rect": '<circle class="faint" cx="12" cy="12" r="9.5"/>' + _GRID4,
    "layout_fibonacci": '<circle class="faint" cx="12" cy="12" r="9.5"/>' + _SPIRAL,
    "layout_fermat": '<path class="faint" d="M3 16Q12 2 21 16"/>' + "".join(       # dots on the curve itself
        f'<circle class="f" cx="{_f((1 - u) ** 2 * 3 + 2 * u * (1 - u) * 12 + u * u * 21)}" '
        f'cy="{_f((1 - u) ** 2 * 16 + 2 * u * (1 - u) * 2 + u * u * 16)}" r="1.35"/>'
        for u in [0.06 + 0.88 * k / 8 for k in range(9)]),
    # grid.async_layout
    "asyncgrid_shared": '<rect x="2.5" y="4.5" width="19" height="15" rx="1.2"/><path d="M6 8L12 15M18 8L12 15"/>'
                        '<circle class="af" cx="12" cy="15" r="1.4"/>',
    "asyncgrid_per_transducer": '<rect x="2" y="4.5" width="9" height="15" rx="1.2"/>'
                                '<rect x="13" y="4.5" width="9" height="15" rx="1.2"/><path d="M4.5 8L8 15M19.5 8L16 15"/>',
    # grid.mode: the grid box is upright (z vertical); accent = what is aligned with z
    "gridmode_focal_axis": _BOX + '<path d="M6 5L12 16M18 5L12 16"/><path class="a" d="M12 3v18" stroke-dasharray="1.6 2"/>',
    "gridmode_transducer_axis": _BOX + '<path class="a" d="M12 4.5V16"/><path d="M19 7L12 16"/>',
    "gridmode_ras_plus": _BOX + '<path d="M6 5L13 16M19 7L13 16"/><path class="a" d="M5.5 18.5h4M5.5 18.5v-4"/>',
    # calibration.focal_reference: axial profile; peak, or centre of the length above peak/2, peak/4
    "focalref_peak": _PROFILE + '<path class="a" d="M10.5 5V19" stroke-dasharray="1.4 2"/><circle class="af" cx="10.5" cy="5" r="1.7"/>',
    "focalref_flhm_3db": _PROFILE + '<path class="a" d="M7.3 12H15.1M11.2 10.5v3"/>',
    "focalref_flhm_6db": _PROFILE + '<path class="a" d="M6.4 15.5H17.4M11.9 14v3"/>',
    # headmodel.smooth_method (kernel shape; 'off' = the untouched mask edge)
    "smooth_gaussian": '<path d="M3 18C8 18 9 6 12 6S16 18 21 18"/>',
    "smooth_box": '<path d="M3 18H8V7H16V18H21"/>',
    "smooth_off": '<path d="M3 18H12V6H21"/>',
    # headmodel.skull_fill_method
    "fill_rubberwrap": '<path d="M4 19L6 12L8.5 15L10.5 8.5L13 13L15 9.5L17.5 14L20 19"/>'
                       '<path class="a" d="M3 19C4 10 8 6.5 12 6.5S20 10 21 19"/>',
    "fill_imclose": '<path d="M4 18C6 11 8.5 10.5 10 12M14 12C15.5 10.5 18 11 20 18"/><path class="a" d="M10 12Q12 13.4 14 12"/>'
                    '<circle class="a" cx="12" cy="7" r="2.4" stroke-dasharray="1.2 1.4"/>',
    "fill_off": '<path d="M4 18C6 11 8.5 10.5 10 12M14 12C15.5 10.5 18 11 20 18"/>',
    # timing.pulse_ramp_shape (pulse envelope)
    "ramp_linear": '<path d="M3 18L7.5 6H16.5L21 18"/>',
    "ramp_tukey": '<path d="M3 18C5.6 18 5.6 6 8.2 6H15.8C18.4 6 18.4 18 21 18"/>',
    # grid.axisymmetric_output ('auto' stays text)
    "axisym_2d": '<path class="a" d="M12 3v18" stroke-dasharray="1.6 2"/><rect x="12" y="4" width="8" height="16"/>'
                 '<rect class="faint" x="4" y="4" width="8" height="16" stroke-dasharray="1.6 2"/>',
    "axisym_3d": '<path class="a" d="M12 3v18" stroke-dasharray="1.6 2"/><ellipse cx="12" cy="6" rx="7" ry="2.4"/>'
                 '<path d="M5 6v12M19 6v12M5 18A7 2.4 0 0 0 19 18"/>',
    # grid.water_baseline_grid
    "baseline_beam": '<rect class="faint" x="3" y="3" width="18" height="18" rx="1" stroke-dasharray="1.6 2"/>'
                     '<rect class="a" x="8" y="4.5" width="8" height="16" rx="0.8"/>' + _SMALL_BOWL_BEAM,
    "baseline_full": '<rect class="a" x="3" y="3" width="18" height="18" rx="1"/>' + _SMALL_BOWL_BEAM,
    # time_reversal.pose_workflow.aim: corrected aim points past the target by the predicted shift
    # (conventions for all glyphs: transducer = thick arc at the top, concave side towards its
    # focus; target = ring; focus = filled dot; accent = what the choice changes)
    "aim_geometric": _BOWL + '<circle cx="12" cy="16" r="3.2"/><path d="M12 5V16"/>',
    "aim_corrected": _BOWL + '<circle cx="12" cy="16" r="3.2"/><path d="M12 5L15.5 19"/><path class="a" d="M14 17.5L17 20.5M17 17.5L14 20.5"/>',
    "aim_both": _BOWL + '<circle cx="12" cy="16" r="3.2"/><path d="M12 5V16"/><path class="a" d="M12 5L15.5 19"/>',
    # time_reversal.placement_search.jitter_mode: ± grid axes (in-plane view) vs. 8 moves on the scalp
    "jitter_grid": '<path class="faint" d="M12 4v16M4 12h16"/><circle class="f" cx="12" cy="12" r="1.6"/>' + "".join(
        f'<circle class="af" cx="{x}" cy="{y}" r="1.4"/>' for x, y in ((12, 5), (12, 19), (5, 12), (19, 12))),
    "jitter_scalp": _SCALP + '<circle class="f" cx="12" cy="15" r="1.6"/>' + "".join(
        f'<circle class="af" cx="{_f(12 + 7 * math.cos(a))}" cy="{_f(15 + 2.4 * math.sin(a))}" r="1.2"/>'
        for a in [i * math.pi / 4 for i in range(8)]),
    # time_reversal.pose_workflow.rank_by
    "rank_target_intensity": '<circle cx="12" cy="12" r="7.5"/><circle class="af" cx="12" cy="12" r="3.5"/>',
    "rank_focus_error": '<circle cx="9" cy="14" r="3.2"/><circle class="af" cx="16.5" cy="7.5" r="2"/>'
                        '<path class="a" d="M11.6 11.8L14.9 8.9" stroke-dasharray="1.4 1.6"/>',
    # placement.mode
    "place_manual": '<circle cx="12" cy="12" r="6.5"/><path d="M12 3v5M12 16v5M3 12h5M16 12h5"/><circle class="af" cx="12" cy="12" r="1.5"/>',
    "place_localite": '<rect x="3" y="3.5" width="18" height="5" rx="2"/><circle class="f" cx="7" cy="6" r="1.1"/>'
                      '<circle class="f" cx="17" cy="6" r="1.1"/><path class="a" d="M7 9L12 18M17 9L12 18" stroke-dasharray="1.4 2"/>'
                      '<circle class="af" cx="12" cy="18.5" r="1.7"/>',
    "place_heuristic": '<circle class="af" cx="12" cy="12" r="1.6"/><circle cx="12" cy="12" r="5"/>'
                       '<circle cx="12" cy="12" r="9" stroke-dasharray="1.6 2"/>',
    "place_mni": '<path d="M4 13C4 7 8 4 12 4S20 7 20 13C20 17 17 20 12 20S4 17 4 13Z"/>'
                 '<path class="faint" d="M12 4V20M4.5 10H19.5M4.5 15.5H19.5M8 5V19.5M16 5V19.5"/>',
    "place_plantus": '<path d="M4 3V20H21"/><path d="M6 6C8 13 13 17 19 18"/><circle class="f" cx="7.2" cy="9.8" r="1.2"/>'
                     '<circle class="af" cx="11" cy="14.6" r="1.6"/><circle class="f" cx="16" cy="17.2" r="1.2"/>',
    "place_tus_entry": _SCALP + '<path class="a" d="M8.5 15.4L15.5 15.4" stroke-width="2.4"/>'
                       '<path d="M8 11.5Q12 7.5 16 11.5" stroke-width="2.2"/>',
}

# ── Figures (96 x 56): body and caption ───────────────────────────────

FIGURES: dict[str, tuple[str, str]] = {
    "coupling_wavefronts_coherent": (
        _wavefront_figure(True),
        "Both transducers fire in one run. Their wavefronts superpose; dots mark where crests reinforce."),
    "coupling_wavefronts_async": (
        _wavefront_figure(False),
        "Each transducer is simulated on its own (dashed: the other run). No interference."),
    "coupling_traces_coherent": (
        _traces_figure(True),
        "Drives 1 and 2 share one window with a fixed phase relation; pressures add (blue: p₁ + p₂)."),
    "coupling_traces_async": (
        _traces_figure(False),
        "Drives fire in separate windows; time-averaged intensities add (blue: I₁ + I₂)."),
    "gridmode_focal_axis_fig": (
        _gridmode_figure("focal_axis"),
        "Grid z (blue, dashed) follows the mean beam axis: here the bisector, each beam tilted "
        "by half the angle between them."),
    "gridmode_transducer_axis_fig": (
        _gridmode_figure("transducer_axis"),
        "Grid z follows transducer 1's beam (blue); the other beams cross the grid at their "
        "full angle."),
    "gridmode_ras_plus_fig": (
        _gridmode_figure("ras_plus"),
        "No rotation: the grid keeps the T1's R and S axes (blue), so every beam is oblique. "
        "Cannot be combined with 2-D axisymmetric simulations (grid.axisymmetric), which need "
        "the beam along the grid z axis."),
    "asyncgrid_shared_fig": (
        _asyncgrid_figure(False),
        "One grid for all transducers, oriented by grid.mode and cropped to hold every beam; "
        "off-axis beams are voxelised at an angle."),
    "asyncgrid_per_transducer_fig": (
        _asyncgrid_figure(True),
        "One grid per transducer, aligned with its beam as in a solo run; each field is "
        "resampled onto transducer 1's grid (blue) to combine them."),
}


# ── Tissues (medium_properties table headers) ──────────────────────────
# A dome cross-section of the head layers; the tissue's own layer is accented, the rest faint.
HEAD_LAYERS = {
    "skin": 'd="M3 7.5Q12 -0.5 21 7.5"',
    "cortical_out": 'd="M3 11.5Q12 3.5 21 11.5"',
    "trabecular": 'd="M3 14.5Q12 6.5 21 14.5" stroke-dasharray="1.4 1.8"',
    "cortical_in": 'd="M3 17.5Q12 9.5 21 17.5"',
    "brain": 'd="M4 21.5q2 -2.5 4 0t4 0t4 0t4 0"',
}
TISSUE_LAYERS = {
    "skin": ("skin",),
    "skull": ("cortical_out", "trabecular", "cortical_in"),
    "skull_cortical": ("cortical_out", "cortical_in"),
    "skull_trabecular": ("trabecular",),
    "brain": ("brain",),
}


def head_layers(accent: tuple[str, ...] = (), faint: bool = True) -> str:
    """The head-layer dome; layers in `accent` drawn in the accent colour."""
    return "".join(f'<path class="{"a" if k in accent else ("faint" if faint else "")}" {d}/>'
                   for k, d in HEAD_LAYERS.items()).replace(' class=""', "")


def tissue_glyph(tissue: str) -> str | None:
    """Glyph name for a medium_properties tissue column, or None."""
    name = f"tissue_{tissue}"
    return name if name in GLYPHS else None


GLYPHS.update({f"tissue_{t}": head_layers(layers) for t, layers in TISSUE_LAYERS.items()})
GLYPHS["tissue_water"] = GLYPHS["medium_water"]


def glyph_body(name: str) -> str:
    return qt(GLYPHS[name])


def figure_body(name: str) -> str:
    return qt(FIGURES[name][0])


def caption(name: str) -> str:
    return FIGURES[name][1]
