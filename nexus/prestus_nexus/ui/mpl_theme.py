"""Matplotlib figures in the Nexus theme (light or dark): ThemedCanvas applies the active
palette (style.active()) to its figure just before every draw.

Only matplotlib's own defaults are recoloured: a white figure or axes background becomes the
card surface, black text, spines and ticks become the theme's ink. Colours a plot sets on
purpose (data, markers, warnings, image backgrounds) are kept, except colours taken from a
Nexus palette (style.C or style.C_DARK, e.g. style.active()["ink3"]): those follow the theme by
their role. Each artist's role is decided once, from its colour when first seen, so switching
themes (or exporting in the light theme, export.py) recolours it again.

ThemedCanvas also carries the figure export button (export.ExportButton), unless created with
``export=False``.
"""

from __future__ import annotations

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.colors import to_hex, to_rgba

from . import style
from .errors import log_current

_BG, _INK = "_nexus_bg", "_nexus_ink"
# containers drawn as cards (style.py: a surface background); a canvas elsewhere is on the page
CARD_NAMES = {"Card", "SubGroup", "LoadFooter"}

# hex colour -> palette role, for colours a plot took from either Nexus palette
_ROLES: dict[str, str] = {}
for _pal in (style.C_DARK, style.C):
    for _k, _v in _pal.items():
        if _k != "on_brand":                      # same white as 'surface' in the light theme
            _ROLES[_v.lower()] = _k


def _is_default_bg(color) -> bool:
    try:
        r, g, b, a = to_rgba(color)
    except (TypeError, ValueError):
        return False
    return a > 0 and min(r, g, b) > 0.94


def _is_default_ink(color) -> bool:
    try:
        r, g, b, a = to_rgba(color)
    except (TypeError, ValueError):
        return False
    return a > 0 and max(r, g, b) < 0.25


def _palette_role(color) -> str | None:
    try:
        r, g, b, a = to_rgba(color)
    except (TypeError, ValueError):
        return None
    if a == 0:
        return None
    return _ROLES.get(to_hex((r, g, b)))


def _recolour(artist, flag: str, test, get, set_, colour, c: dict | None = None) -> None:
    """Set artist's colour to ``colour`` if it is a matplotlib default (``test``), or to the
    same role in palette ``c`` if it came from a Nexus palette; leave chosen colours alone."""
    role = getattr(artist, flag, None)
    if role is None:
        current = get()
        role = _palette_role(current) or bool(test(current))
        try:
            setattr(artist, flag, role)
        except AttributeError:          # pragma: no cover
            return
    if role is True:
        set_(colour)
    elif role and c is not None and role in c:
        set_(c[role])


def _single(colours):
    """The one colour of a collection drawn in a single colour, else None."""
    try:
        return colours[0] if len(colours) == 1 else None
    except TypeError:
        return None


def _palette_only(artist, flag: str, get, set_, c: dict) -> None:
    """Recolour by role only a colour taken from a Nexus palette (data colours are kept)."""
    _recolour(artist, flag, lambda _x: False, get, set_, None, c)


def _text(t, c) -> None:
    _recolour(t, _INK, _is_default_ink, t.get_color, t.set_color, c["ink2"], c)


def _legend(leg, c) -> None:
    if leg is None:
        return
    frame = leg.get_frame()
    _recolour(frame, _BG, _is_default_bg, frame.get_facecolor, frame.set_facecolor, c["surface"], c)
    frame.set_edgecolor(c["line"] if "line" in c else c["ink4"])
    for t in leg.get_texts():
        _text(t, c)
    if leg.get_title() is not None:
        _text(leg.get_title(), c)


def _axis_ticks(ax, axis, c) -> None:
    """Tick marks and labels, also those matplotlib creates while drawing (tick_params): the
    theme's ink unless the plot set a tick colour itself; a colour it took from a Nexus palette
    follows the theme by its role."""
    roles = getattr(axis, "_nexus_ticks", None)
    if roles is None:
        kw = getattr(axis, "_major_tick_kw", {})
        tick = _palette_role(kw["color"]) if "color" in kw else "ink3"
        label = _palette_role(kw["labelcolor"]) if "labelcolor" in kw else ("ink2" if "color" not in kw else tick)
        roles = axis._nexus_ticks = (tick, label)
    tick, label = roles
    kw = {}
    if tick in c:
        kw["colors"] = c[tick]
    if label in c:
        kw["labelcolor"] = c[label]
    if kw:
        ax.tick_params(axis=axis.axis_name, which="both", **kw)


def apply(fig, c: dict | None = None) -> None:
    """Recolour fig's default background, text, spines, ticks and legends to the palette c
    (default: the active Nexus theme)."""
    c = c or style.active()
    _recolour(fig.patch, _BG, _is_default_bg, fig.patch.get_facecolor, fig.patch.set_facecolor, c["surface"], c)
    for t in fig.texts:
        _text(t, c)
    if getattr(fig, "_suptitle", None) is not None:
        _text(fig._suptitle, c)
    for leg in fig.legends:
        _legend(leg, c)
    for ax in fig.axes:
        _recolour(ax.patch, _BG, _is_default_bg, ax.patch.get_facecolor, ax.patch.set_facecolor, c["surface"], c)
        for sp in ax.spines.values():
            _recolour(sp, _INK, _is_default_ink, sp.get_edgecolor, sp.set_edgecolor, c["ink3"], c)
        axes = [ax.xaxis, ax.yaxis] + ([ax.zaxis] if hasattr(ax, "zaxis") else [])
        for axis in axes:
            _axis_ticks(ax, axis, c)
            _text(axis.label, c)
            if hasattr(axis, "set_pane_color"):         # 3D: the panes behind the data
                pane = axis.pane
                _recolour(pane, _BG, _is_default_bg, pane.get_facecolor, axis.set_pane_color, c["surface2"], c)
        for t in (ax.title, ax._left_title, ax._right_title, *ax.texts):
            _text(t, c)
        for ln in ax.lines:                     # lines and markers drawn in a palette colour
            _palette_only(ln, "_nexus_c", ln.get_color, ln.set_color, c)
            _palette_only(ln, "_nexus_mfc", ln.get_markerfacecolor, ln.set_markerfacecolor, c)
            _palette_only(ln, "_nexus_mec", ln.get_markeredgecolor, ln.set_markeredgecolor, c)
        for coll in ax.collections:             # vlines / hlines / scatter in one palette colour
            _palette_only(coll, "_nexus_fc", lambda k=coll: _single(k.get_facecolor()), coll.set_facecolor, c)
            _palette_only(coll, "_nexus_ec", lambda k=coll: _single(k.get_edgecolor()), coll.set_edgecolor, c)
        _legend(ax.get_legend(), c)


class ThemedCanvas(FigureCanvasQTAgg):
    """FigureCanvasQTAgg drawn in the Nexus theme (see apply), with the export button
    (export.ExportButton) in its top-right corner unless ``export=False``.

    ``export_name`` names the figure (the default filename); ``filename_hint`` (callable ->
    filename stem, e.g. export.subject_stem(model, "beam")) overrides it. The button hides
    while the figure has no axes."""

    def __init__(self, figure=None, *, export: bool = True, export_name: str = "figure",
                 filename_hint=None):
        super().__init__(figure)
        self.export_button = None
        if export:
            from .export import ExportButton, MplExporter
            self.export_button = ExportButton(self, MplExporter(self), export_name, filename_hint)

    def _in_card(self) -> bool:
        """Whether a card (white panel) holds this canvas: then the figure stays card white."""
        w = self.parentWidget()
        while w is not None:
            if w.objectName() in CARD_NAMES:
                return True
            w = w.parentWidget()
        return False

    def wheelEvent(self, event) -> None:
        """Scroll the page under a plot: matplotlib accepts every wheel event, so a canvas inside
        a scroll area stops it scrolling. Only canvases with a scroll_event handler (slice and
        zoom views) keep the wheel (matplotlib's own handlers do not count)."""
        if getattr(self, "_uses_wheel", False):
            super().wheelEvent(event)
        else:
            event.ignore()

    def mpl_connect(self, s, func):
        if s == "scroll_event":
            self._uses_wheel = True
        return super().mpl_connect(s, func)

    def draw(self) -> None:
        try:
            apply(self.figure)
            # a figure on the page itself (not inside a card) takes the page's colour instead of
            # the card white; exports keep the card colour (export.MplExporter: white)
            # (role True: matplotlib's default white; 'surface': the card colour set on purpose)
            if getattr(self.figure.patch, _BG, None) in (True, "surface") and not self._in_card():
                self.figure.patch.set_facecolor(style.active()["bg"])
        except Exception:  # best-effort: theming must never stop a plot
            log_current("plot theme")
        super().draw()
        if self.export_button is not None:
            self.export_button.set_available(any(ax.get_visible() for ax in self.figure.axes))
