"""Hover help for a config field: the config_default.yaml comment reflowed into
paragraphs, its choices as a list (with their glyphs), and the key, default and
tuning goals in a footer.

Help text is written as YAML comments, wrapped at ~80 columns, with choices as
lines like  'coherent'  — text  or inline as  'a' text; 'b' text  or  'a' | 'b'.
"""

from __future__ import annotations

import html
import re
import urllib.parse

from .. import docs, labels
from ..schema import Node

WIDTH = 420                 # px; tooltips wider than a short sentence are wrapped to this
FAINT = "#9aa7b4"          # footer ink: the tooltip background is dark in both themes (style.py)
TIP_BG = "#10161d"         # the QToolTip background of style.py, for the clickable help popup
LINK = "#5aa8e6"           # link ink on that background (the glyphs' accent)
_CHOICE_LINE = re.compile(r"^'([^']+)'\s*(\(default\))?\s*(?:—|–|-{1,2}|:)?\s*(.*)$")
_DEFAULT = re.compile(r"\s*\(default\)\s*")


def _clean(text: str) -> tuple[str, bool]:
    """Item text without a '(default)' marker or wrapping parentheses; whether it was marked."""
    marked = bool(_DEFAULT.search(text))
    text = _DEFAULT.sub(" ", text).strip()
    if text.startswith("(") and text.endswith(")") and text.count("(") == 1:
        text = text[1:-1].strip()
    return text, marked


def _quoted(choice) -> str:
    return re.escape(f"'{choice}'")


def _split_inline(par: str, choices: list) -> tuple[str, list[tuple[str, str]]] | None:
    """'intro 'a' text; 'b' text' → (intro, [(a, text), (b, text)]) when at least two
    choices are listed inline, separated by ';' or '|'."""
    names = [str(c) for c in choices or []]
    hits = sorted((m.start(), n) for n in names for m in re.finditer(_quoted(n), par))
    if len(hits) < 2:
        return None
    cuts = [p for p, _ in hits]
    for a, b in zip(cuts, cuts[1:]):
        if not re.search(r"[;|]\s*$", par[a:b].rstrip()):
            return None
    items = []
    for (p, n), end in zip(hits, cuts[1:] + [len(par)]):
        text = par[p + len(n) + 2:end].strip().rstrip(";|").strip()
        items.append((n, text.lstrip(":—– ").strip()))
    return par[:cuts[0]].strip(), items


def parse(text: str, choices: list | None = None) -> list:
    """Blocks of the help text: ('p', paragraph) and ('items', [(choice, text, is_default)])."""
    blocks: list = []
    par: list[str] = []
    items: list[list] = []

    def flush_par():
        if par:
            blocks.append(("p", " ".join(par)))
            par.clear()

    def flush_items():
        if items:
            blocks.append(("items", [tuple(i) for i in items]))
            items.clear()

    for raw in (text or "").splitlines():
        line = raw.strip()
        if not line:
            flush_par(); flush_items()
            continue
        m = _CHOICE_LINE.match(line)
        if m and (not choices or m.group(1) in [str(c) for c in choices]):
            flush_par()
            items.append([m.group(1), m.group(3), bool(m.group(2))])
        elif items and not par:
            items[-1][1] = f"{items[-1][1]} {line}".strip()     # continuation of the item
        else:
            par.append(line)
    flush_par(); flush_items()
    blocks = [(k, [(c, *_cleaned(t, mk)) for c, t, mk in v]) if k == "items" else (k, v) for k, v in blocks]

    out = []
    for kind, val in blocks:            # inline lists inside paragraphs
        split = _split_inline(val, choices) if kind == "p" else None
        if split:
            intro, its = split
            if intro:
                out.append(("p", intro))
            out.append(("items", [(c, *_clean(t)) for c, t in its]))
        else:
            out.append((kind, val))
    return out


def _cleaned(text: str, marked: bool) -> tuple[str, bool]:
    text, m = _clean(text)
    return text, marked or m


def _esc(s: str) -> str:
    s = html.escape(s, quote=False)
    s = s.replace("&gt;=", "≥").replace("&lt;=", "≤").replace(" -- ", " – ")
    return labels.rich(s)               # ITRUSST notation: I<sub>SPPA</sub>, I<sub>SPTA</sub>


def _glyph_img(name: str) -> str:
    """The glyph as SVG markup in a data URL (text, no encoded binary), drawn at 48 px
    and shown at 16 px so it stays sharp on high-resolution screens."""
    from . import icons
    svg = icons.svg(name, "#e8edf2", "#5aa8e6").replace("<svg ", '<svg width="48" height="48" ', 1)
    return f'<img src="data:image/svg+xml;utf8,{urllib.parse.quote(svg)}" width="16" height="16">'


def html_for(node: Node, title: str, value_text: str | None, icons_ok: bool = True, note: str = "",
             links: bool = False, doc: str = "") -> str:
    """Rich-text tooltip for `node` (title = its display label, value_text = its default as
    shown in the editor, None for groups; note = a context line after the help). links: doc_*.md
    mentions become links and `doc` (a doc reference) is added as a 'Docs' row, for the
    clickable popup (show_help); a hover tooltip cannot be clicked."""
    def text(s: str) -> str:
        return docs.linkify(_esc(s), LINK) if links else _esc(s)

    blocks = parse(node.help, node.choices)
    body = [f"<div style='font-weight:600; font-size:13px'>{_esc(title)}</div>"]
    default = node.default if not node.is_group else None
    if blocks and blocks[0][0] == "p" and blocks[0][1].rstrip(" :").lower() == title.lower():
        blocks = blocks[1:]                 # an intro that only repeats the title
    for kind, val in blocks:
        if kind == "p":
            body.append(f"<p style='margin-top:4px'>{text(val)}</p>")
            continue
        rows = []
        for choice, txt, _ in val:     # '(default)' in the help is dropped: the value below is authoritative
            icon = node.choice_icons.get(choice) if icons_ok else None
            img = _glyph_img(icon) if icon else ""
            tag = (" <span style='color:%s'>(default)</span>" % FAINT
                   if default is not None and str(choice) == str(default) else "")
            rows.append(f"<tr><td style='padding:2px 6px 2px 0' valign='top'>{img}</td>"
                        f"<td style='padding:2px 10px 2px 0' valign='top'><b>{_esc(choice)}</b></td>"
                        f"<td style='padding:2px 0' valign='top'>{text(txt)}{tag}</td></tr>")
        body.append("<table cellspacing='0' style='margin-top:4px'>" + "".join(rows) + "</table>")

    if note:
        body.append(f"<p style='margin-top:4px'><i>{_esc(note)}</i></p>")
    meta = [("Key", f"<code>{_esc(node.path)}</code>")]
    if value_text is not None:
        meta.append(("Default", _esc(value_text or "(empty)")))
    if node.tune:
        meta.append(("Tuning", "affects " + _esc(", ".join(node.tune))))
    if node.extra:
        meta.append(("Source", "read by PRESTUS code; not in config_default.yaml"))
    if links and doc:
        meta.append(("Docs", f"<a href='{_esc(doc)}' style='color:{LINK}'>{_esc(docs.parse_ref(doc)[0])} ↗</a>"))
    foot = "".join(f"<tr><td style='padding-right:10px'>{k}</td><td>{v}</td></tr>" for k, v in meta)
    body.append(f"<table cellspacing='0' style='margin-top:8px; color:{FAINT}'>{foot}</table>")

    text_len = len(node.help or "") + len(title)
    width = f" width='{WIDTH}'" if text_len > 70 else ""
    return f"<table{width} cellspacing='0'><tr><td>{''.join(body)}</td></tr></table>"


def show_help(html_text: str, anchor=None) -> None:
    """Show help HTML in a popup at the cursor whose doc links can be clicked (a QToolTip
    closes before a click reaches it). Closes on a click outside, Escape, or a followed link."""
    from PyQt6.QtCore import QPoint, Qt
    from PyQt6.QtGui import QCursor
    from PyQt6.QtWidgets import QFrame, QLabel, QVBoxLayout

    pop = QFrame(anchor, Qt.WindowType.Popup)
    pop.setObjectName("HelpPopup")
    pop.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
    pop.setStyleSheet(f"#HelpPopup {{ background: {TIP_BG}; border-radius: 4px; }}"
                      "#HelpPopup QLabel { color: #ffffff; background: transparent; }")
    lay = QVBoxLayout(pop)
    lay.setContentsMargins(8, 6, 8, 6)
    label = QLabel(html_text)
    label.setTextFormat(Qt.TextFormat.RichText)
    label.setTextInteractionFlags(Qt.TextInteractionFlag.TextBrowserInteraction)
    label.setOpenExternalLinks(False)

    def follow(link: str) -> None:
        pop.close()
        docs.open_doc(link)
    label.linkActivated.connect(follow)
    lay.addWidget(label)
    pop.adjustSize()
    pos = QCursor.pos() + QPoint(10, 12)
    screen = (anchor.screen() if anchor is not None else None)
    if screen is not None:
        geo = screen.availableGeometry()
        pos.setX(max(geo.left(), min(pos.x(), geo.right() - pop.width())))
        pos.setY(max(geo.top(), min(pos.y(), geo.bottom() - pop.height())))
    pop.move(pos)
    pop.show()


class MarkerTips:
    """Hover text for points drawn on a matplotlib canvas (the transducer markers of the
    slice views). The host registers each marker's data position while drawing; a mouse move
    compares the pointer with those few points in display pixels, so hovering stays cheap."""

    RADIUS_PX = 8           # pick radius in logical pixels

    def __init__(self, canvas):
        self.canvas = canvas
        self.points: dict = {}          # axes -> [(x, y, text)] in data coordinates
        self._shown = ""
        canvas.mpl_connect("motion_notify_event", self._moved)
        canvas.mpl_connect("figure_leave_event", lambda _e: self._show(""))

    def clear(self, ax=None) -> None:
        if ax is None:
            self.points.clear()
        else:
            self.points.pop(ax, None)

    def add(self, ax, x, y, text: str) -> None:
        x, y = float(x), float(y)
        if text and x == x and y == y:          # NaN-free
            self.points.setdefault(ax, []).append((x, y, text))

    def at(self, ax, px: float, py: float) -> str:
        """Text of the marker nearest to display point (px, py) within the pick radius."""
        pts = self.points.get(ax)
        if not pts:
            return ""
        import numpy as np
        disp = ax.transData.transform(np.array([(x, y) for x, y, _ in pts]))
        d = np.hypot(disp[:, 0] - px, disp[:, 1] - py)
        i = int(np.argmin(d))
        ratio = getattr(self.canvas, "device_pixel_ratio", 1) or 1
        return pts[i][2] if d[i] <= self.RADIUS_PX * ratio else ""

    def _moved(self, event) -> None:
        self._show(self.at(event.inaxes, event.x, event.y) if event.inaxes is not None else "")

    def _show(self, text: str) -> None:
        if text == self._shown:
            return
        self._shown = text
        from PyQt6.QtGui import QCursor
        from PyQt6.QtWidgets import QToolTip
        if text:
            QToolTip.showText(QCursor.pos(), text, self.canvas)
        else:
            QToolTip.hideText()
