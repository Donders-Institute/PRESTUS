"""Line icons for the Nexus chrome, drawn here as small SVGs on a 24 × 24 grid and tinted
with the active palette, so they follow the light and dark theme.

    icon("advanced")                 # ink colour of the active palette
    icon("save", style.active()["on_brand"])
    icon("coupling_async")           # a choice glyph (choice_icons.py): ink + accent
    figure("coupling_traces_async")  # a 96 x 56 choice figure as a pixmap
"""

from __future__ import annotations

from PyQt6.QtCore import QByteArray, QEvent, QRectF, Qt, QUrl
from PyQt6.QtGui import QIcon, QPainter, QPixmap, QTextDocument
from PyQt6.QtWidgets import QToolButton

from . import choice_icons, style

# Stroke paths (no fill) unless a shape sets fill="currentColor" itself.
_SHAPES = {
    # three sliders: the advanced settings
    "advanced": '<path d="M4 6h16M4 12h16M4 18h16"/>'
                '<circle cx="9" cy="6" r="2.2" fill="currentColor"/>'
                '<circle cx="15" cy="12" r="2.2" fill="currentColor"/>'
                '<circle cx="7" cy="18" r="2.2" fill="currentColor"/>',
    "help": '<circle cx="12" cy="12" r="9"/>'
            '<path d="M9.2 9.3a2.9 2.9 0 0 1 5.6 1c0 1.9-2.8 2.3-2.8 4.2"/>'
            '<circle cx="12" cy="17.6" r="0.9" fill="currentColor" stroke="none"/>',
    "reset": '<path d="M4 12a8 8 0 1 0 2.4-5.7L4 8.6"/><path d="M4 3.6v5h5"/>',
    "lock": '<rect x="5" y="11" width="14" height="10" rx="2"/><path d="M8 11V7.5a4 4 0 0 1 8 0V11"/>',
    "unlock": '<rect x="5" y="11" width="14" height="10" rx="2"/><path d="M8 11V7.5a4 4 0 0 1 7.6-1.7"/>',
    "folder": '<path d="M3 7a2 2 0 0 1 2-2h4l2 2h8a2 2 0 0 1 2 2v8a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2z"/>',
    # a sheet with a folded corner: a file
    "file": '<path d="M6.5 3h7.5l4.5 4.5V19a2 2 0 0 1-2 2h-10a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2z"/><path d="M14 3v4.5h4.5"/>',
    "save": '<path d="M5 3h11l3 3v13a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2z"/>'
            '<path d="M8 3v5h7V3M7 21v-7h10v7"/>',
    "collapse": '<path d="M6 15l6-6 6 6"/>',
    "expand": '<path d="M6 9l6 6 6-6"/>',
    # an arrow into a tray: export a figure
    "download": '<path d="M12 4v11"/><path d="M7.5 10.5L12 15l4.5-4.5"/><path d="M4.5 15.5v3a1.5 1.5 0 0 0 1.5 1.5h12a1.5 1.5 0 0 0 1.5-1.5v-3"/>',
    "plus": '<path d="M12 5v14M5 12h14"/>',
    "undo": '<path d="M9 14L4 9l5-5"/><path d="M4 9h10.5a5.5 5.5 0 0 1 0 11H11"/>',
    "redo": '<path d="M15 14l5-5-5-5"/><path d="M20 9H9.5a5.5 5.5 0 0 0 0 11H13"/>',
    "minus": '<path d="M5 12h14"/>',
    "ok": '<circle cx="12" cy="12" r="9"/><path d="M8 12.2l2.8 2.8L16.2 9.5"/>',
    "warning": '<path d="M12 3.5l9.5 16.5h-19z"/><path d="M12 10v4.5"/>'
               '<circle cx="12" cy="17.2" r="0.9" fill="currentColor" stroke="none"/>',
    # section cards (overlay 'sections: <name>: icon')
    # the head's layers (skin, cortical / trabecular / cortical skull, brain): the SimNIBS
    # segmentation; the same dome as the tissue glyphs (choice_icons.HEAD_LAYERS)
    "segmentation": choice_icons.head_layers(faint=False),
    # two error bars: uncertainty bounds (simulation.uncertainty)
    "range": '<path d="M7 4.5v15M5 4.5h4M5 19.5h4"/><circle cx="7" cy="11" r="1.5" fill="currentColor" stroke="none"/>'
             '<path d="M17 7.5v10M15 7.5h4M15 17.5h4"/><circle cx="17" cy="13" r="1.5" fill="currentColor" stroke="none"/>',
    # a wand with sparks: refining the segmentation
    "refine": '<path d="M4 20L14.5 9.5"/><path d="M13 8l3 3"/>'
              '<path d="M18 3v3.5M16.25 4.75h3.5M20.5 10v2.5M19.25 11.25h2.5M9 3.5v2M8 4.5h2"/>',
    # a jagged edge becoming a smooth curve: smoothing
    "smooth": '<path d="M3 19l2-3 2 3 2-3 2 3"/><path d="M3 11c3 0 3.5-6 7-6s4 6 7 6 2.5-2 4-2"/>',
    # a grid with one finer quadrant: grid refinement
    "grid": '<rect x="4" y="4" width="16" height="16" rx="1.5"/><path d="M4 12h16M12 4v16"/>'
            '<path d="M4 8h8M8 4v8" stroke-width="1.1"/>',
    # stacked sheets: tissue layers
    "layers": '<path d="M12 4l8 4-8 4-8-4z"/><path d="M4 12l8 4 8-4"/><path d="M4 16l8 4 8-4"/>',
    # a table: tissue properties
    "table": '<rect x="4" y="5" width="16" height="14" rx="1.5"/><path d="M4 10h16M4 14.5h16M10 5v14"/>',
    # a CT gantry: pseudo-CT
    "ct": '<circle cx="12" cy="11" r="7.5"/><circle cx="12" cy="11" r="3"/><path d="M5 20.5h14"/>',
    # section cards and sub-tabs
    # a person: the subject
    "subject": '<circle cx="12" cy="8" r="3.6"/><path d="M5 20c.6-3.8 3.4-6 7-6s6.4 2.2 7 6"/>',
    # a wrench: the environment (external tools)
    "wrench": '<path d="M14.5 4.2a4.6 4.6 0 0 0-4.9 6.2L4 16v4h4l5.6-5.6a4.6 4.6 0 0 0 6.2-4.9l-2.9 2.9-3.2-.8-.8-3.2z"/>',
    # arrows in and out of a tray: input / output
    "io": '<path d="M4.5 15.5v3a1.5 1.5 0 0 0 1.5 1.5h12a1.5 1.5 0 0 0 1.5-1.5v-3"/><path d="M8.5 4v9M5.5 10l3 3 3-3"/><path d="M15.5 13V4M12.5 7l3-3 3 3"/>',
    # a pulse train: sonication timing
    "pulse": '<path d="M3 17h3V7h4v10h4V7h4v10h3"/>',
    # a chip: the compute device
    "chip": '<rect x="6.5" y="6.5" width="11" height="11" rx="1.5"/><rect x="9.5" y="9.5" width="5" height="5"/><path d="M9.5 3.5v3M14.5 3.5v3M9.5 17.5v3M14.5 17.5v3M3.5 9.5h3M3.5 14.5h3M17.5 9.5h3M17.5 14.5h3"/>',
    # two people: a study of several runs
    "people": '<circle cx="9" cy="8.5" r="3"/><path d="M3.5 19c.5-3.2 2.7-5 5.5-5s5 1.8 5.5 5"/><circle cx="16.5" cy="9.5" r="2.4"/><path d="M15.5 14.2c2.6-.2 4.5 1.4 5 4.3"/>',
    # boxes joined by arrows: pipeline stages
    "stages": '<rect x="3" y="9" width="5" height="6" rx="1"/><rect x="9.5" y="9" width="5" height="6" rx="1"/><rect x="16" y="9" width="5" height="6" rx="1"/><path d="M8 12h1.5M14.5 12H16"/>',
    # a server rack: the cluster
    "server": '<rect x="4" y="4" width="16" height="6.5" rx="1.5"/><rect x="4" y="13.5" width="16" height="6.5" rx="1.5"/><circle cx="7.5" cy="7.25" r=".9" fill="currentColor" stroke="none"/><circle cx="7.5" cy="16.75" r=".9" fill="currentColor" stroke="none"/><path d="M11 7.25h6M11 16.75h6"/>',
    # a play triangle: launch a run
    "play": '<circle cx="12" cy="12" r="9"/><path d="M10 8.5v7l6-3.5z"/>',
    # a cube: the 3D view
    "cube": '<path d="M12 3l8 4.5v9L12 21l-8-4.5v-9z"/><path d="M4 7.5l8 4.5 8-4.5M12 12v9"/>',
    # a box inside a dashed margin: padding around the head
    "pad": '<rect x="3.5" y="3.5" width="17" height="17" rx="2" stroke-dasharray="2.2 2"/><rect x="8" y="8" width="8" height="8" rx="1"/>',
    # a list with markers: the job list
    "jobs": '<circle cx="5.5" cy="7" r="1" fill="currentColor" stroke="none"/><circle cx="5.5" cy="12" r="1" fill="currentColor" stroke="none"/><circle cx="5.5" cy="17" r="1" fill="currentColor" stroke="none"/><path d="M9 7h10M9 12h10M9 17h10"/>',
    # pipeline steps (Run ▸ Plan, Skull CT flow)
    # an MR slice: a head outline in a scan frame
    "mri": '<rect x="3.5" y="3.5" width="17" height="17" rx="2"/><path d="M12 6.5c-3 0-5 2.4-5 5.4 0 2.2 1 3.6 2.4 4.6v1.5h5.2v-1.5c1.4-1 2.4-2.4 2.4-4.6 0-3-2-5.4-5-5.4z"/>',
    # a histogram with its peak marked: soft-tissue normalisation
    "histogram": '<path d="M4 20h16"/><path d="M6.5 20v-3M9.5 20v-7M12.5 20v-12M15.5 20v-6M18.5 20v-2"/><path d="M12.5 4.5v1" stroke-dasharray="1 1.5"/>',
    # axes with a straight line: a value-to-value mapping
    "line_map": '<path d="M4.5 4v15.5H20"/><path d="M7.5 16.5L18.5 6.5"/><circle cx="7.5" cy="16.5" r="1.2" fill="currentColor" stroke="none"/><circle cx="18.5" cy="6.5" r="1.2" fill="currentColor" stroke="none"/>',
    # wavefronts leaving a source: the acoustic simulation
    "wave": '<path d="M4 12h2"/><path d="M9 7.5a6.5 6.5 0 0 1 0 9"/><path d="M13 5a10 10 0 0 1 0 14"/><path d="M17 3a13 13 0 0 1 0 18"/>',
    # a target: targeting accuracy
    "target": '<circle cx="12" cy="12" r="8"/><circle cx="12" cy="12" r="4"/><circle cx="12" cy="12" r="1" fill="currentColor" stroke="none"/>',
    # a thermometer: the thermal simulation
    "thermometer": '<path d="M10 14.5V5a2 2 0 0 1 4 0v9.5a4 4 0 1 1-4 0z"/><path d="M12 9v7.5"/>',
    # bars on axes: an analysis
    "chart": '<path d="M4.5 4v15.5H20"/><path d="M8.5 16v-4M12.5 16v-7M16.5 16v-5"/>',
    # a sheet with text lines: a report
    "report": '<path d="M6.5 3h11a2 2 0 0 1 2 2v14a2 2 0 0 1-2 2h-11a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2z"/><path d="M8 8h8M8 11.5h8M8 15h5"/>',
}

_CACHE: dict[tuple[str, str, int], QIcon] = {}


def _wrap(body: str, colour: str, accent: str, box: str, width: float) -> str:
    body = body.replace("currentColor", colour).replace("accentColor", accent)
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{box}" fill="none" '
            f'stroke="{colour}" stroke-width="{width}" stroke-linecap="round" '
            f'stroke-linejoin="round">{body}</svg>')


def svg(name: str, colour: str, accent: str | None = None) -> str:
    """SVG text of a chrome icon or a choice glyph; `accent` (default: `colour`) tints the
    part of a choice glyph that marks what the choice changes."""
    body = _SHAPES[name] if name in _SHAPES else choice_icons.glyph_body(name)
    return _wrap(body, colour, accent or colour, "0 0 24 24", 1.8)


def icon(name: str, colour: str | None = None, size: int = 16, accent: str | None = None) -> QIcon:
    """The named icon in `colour` (default: the active palette's secondary ink) and, for
    choice glyphs, `accent` (default: the palette accent), sharp on high-resolution
    screens; a disabled state is drawn in the faint ink."""
    c = style.active()
    colour = colour or c["ink2"]
    accent = accent or c["accent"]
    key = (name, colour, accent, size)
    if key not in _CACHE:
        ic = QIcon()
        ic.addPixmap(_render(svg(name, colour, accent), size, size), QIcon.Mode.Normal)
        ic.addPixmap(_render(svg(name, c["ink4"]), size, size), QIcon.Mode.Disabled)
        _CACHE[key] = ic
    return _CACHE[key]


def for_choice(node, value) -> QIcon | None:
    """The glyph of `value` of a choice key (overlay 'choice_icons'), or None."""
    name = (node.choice_icons or {}).get(value) if node is not None and isinstance(value, str) else None
    return icon(name) if name else None


def figure(name: str, width: int = 160) -> QPixmap:
    """A choice figure (96 x 56 grid, choice_icons.FIGURES) at `width` px, in the active
    palette's secondary ink and accent."""
    c = style.active()
    text = _wrap(choice_icons.figure_body(name), c["ink2"], c["accent"], "0 0 96 56", 1.4)
    return _render(text, width, round(width * 56 / 96))


def _render(text: str, width: int, height: int) -> QPixmap:
    from PyQt6.QtSvg import QSvgRenderer
    scale = 3                                     # drawn at 3x, so it stays sharp on retina
    pix = QPixmap(width * scale, height * scale)
    pix.fill(Qt.GlobalColor.transparent)
    painter = QPainter(pix)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    QSvgRenderer(QByteArray(text.encode())).render(painter, QRectF(0, 0, width * scale, height * scale))
    painter.end()
    pix.setDevicePixelRatio(scale)
    return pix


def names() -> list[str]:
    """Chrome icons and choice glyphs."""
    return list(_SHAPES) + list(choice_icons.GLYPHS)


def figure_names() -> list[str]:
    return list(choice_icons.FIGURES)


def html_img(doc: QTextDocument, name: str, colour: str | None = None, size: int = 14) -> str:
    """An <img> tag showing icon `name` inside a QTextDocument (QTextBrowser, QLabel rich
    text): the icon is drawn from its SVG and registered as an image resource of `doc`
    (no encoded data in the HTML)."""
    colour = colour or style.active()["ink2"]
    url = f"nexus-icon:{name}/{colour.lstrip('#')}/{size}"
    if doc.resource(QTextDocument.ResourceType.ImageResource.value, QUrl(url)) is None:
        pix = _render(svg(name, colour), size, size)
        doc.addResource(QTextDocument.ResourceType.ImageResource.value, QUrl(url), pix.toImage())
    return f'<img src="{url}" width="{size}" height="{size}" style="vertical-align: middle">'


class PathButton(QToolButton):
    """A browse button with a folder ('dir') or file icon in the active palette, redrawn when
    the theme changes."""

    def __init__(self, kind: str = "dir", tooltip: str | None = None, parent=None):
        super().__init__(parent)
        self._name = "file" if kind == "file" else "folder"
        self.setIcon(icon(self._name))
        self.setToolTip(tooltip or ("Choose a file" if self._name == "file" else "Choose a folder"))

    def changeEvent(self, event) -> None:
        if event.type() in (QEvent.Type.StyleChange, QEvent.Type.PaletteChange):
            self.setIcon(icon(self._name))
        super().changeEvent(event)

