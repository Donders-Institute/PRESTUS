"""Skull CT sub-tab card: the equations of the selected UTE→pHU (pseudo-CT only), density, sound
speed and attenuation mappings, read from the PRESTUS docs so that they stay in step with them:

- pct.skull_mapping       doc/doc_pseudoCT_generation.md, table under 'Choosing the UTE-to-HU Mapping'
- pct.mapping_density     doc/doc_pseudoCT_property-mappings.md, 'Density equations'
- pct.mapping_soundspeed  same page, 'Sound speed equations'
- pct.mapping_attenuation same page, 'Attenuation equations'

Display equations ($$…$$) and inline maths ($…$) are rendered with matplotlib mathtext; without
matplotlib they are shown as their LaTeX source.
"""

from __future__ import annotations

import io
import re
from pathlib import Path

from PyQt6.QtCore import QUrl
from PyQt6.QtGui import QImage, QTextDocument
from PyQt6.QtWidgets import QFrame, QHBoxLayout, QLabel, QTextBrowser, QToolButton, QVBoxLayout

from ..model import is_empty
from . import icons, style

GENERATION = "doc_pseudoCT_generation.md"
MAPPINGS = "doc_pseudoCT_property-mappings.md"
PROPERTIES = [("pct.mapping_density", "Density", "Mapping Skull Density"),
              ("pct.mapping_soundspeed", "Sound speed", "Mapping Skull Sound Speed"),
              ("pct.mapping_attenuation", "Attenuation", "Mapping Skull Attenuation")]


def _section(text: str, heading: str) -> str:
    m = re.search(rf"^## {re.escape(heading)}\s*$(.*?)(?=^## |\Z)", text, re.M | re.S)
    return m.group(1) if m else ""


def table_rows(text: str) -> dict[str, str]:
    """{value: description} from '| `value` ... | description |' rows."""
    out = {}
    for line in text.splitlines():
        m = re.match(r"\|\s*`([^`]+)`[^|]*\|\s*(.*?)\s*\|\s*$", line)
        if m:
            out[m.group(1)] = m.group(2)
    return out


def equation_blocks(text: str) -> dict[str, str]:
    """{value: text} from the '**`value`.**' paragraphs of an equations note."""
    lines = [ln[4:] if ln.startswith("    ") else ln for ln in text.splitlines()]
    body = "\n".join(lines)
    out = {}
    parts = re.split(r"^\*\*`([^`]+)`\.\*\*", body, flags=re.M)
    for name, chunk in zip(parts[1::2], parts[2::2]):
        chunk = re.split(r"^(?:\?\?\?|## )", chunk, maxsplit=1, flags=re.M)[0]
        out[name] = chunk.strip()
    return out


def load(doc_dir: Path) -> dict:
    """{'skull_mapping': {value: desc}, key: {value: (desc, equations)}} from the docs."""
    out: dict = {"skull_mapping": {}}
    try:
        gen = (doc_dir / GENERATION).read_text(encoding="utf-8")
        out["skull_mapping"] = table_rows(_section(gen, "Choosing the UTE-to-HU Mapping"))
    except OSError:
        pass
    try:
        maps = (doc_dir / MAPPINGS).read_text(encoding="utf-8")
    except OSError:
        maps = ""
    for key, _label, heading in PROPERTIES:
        sec = _section(maps, heading)
        rows, eqs = table_rows(sec), equation_blocks(sec)
        out[key] = {v: (rows.get(v, ""), eqs.get(v, "")) for v in {*rows, *eqs}}
    return out


def _mathtext(tex: str) -> str:
    """LaTeX of the docs → matplotlib mathtext."""
    tex = re.sub(r"\\text\{([^}]*)\}", lambda m: r"\mathrm{" + m.group(1).replace("_", r"\_") + "}", tex)
    return tex.replace(r"\,", r"\ ")


def render_math(tex: str, colour: str, size: int) -> QImage | None:
    """A LaTeX line rendered with matplotlib mathtext; None when it cannot be."""
    return _render(tex, colour, size)


def _screen_scale() -> tuple[float, float]:
    """(logical px per point, device pixel ratio) of the primary screen; (1, 1) without one."""
    from PyQt6.QtGui import QGuiApplication
    screen = QGuiApplication.primaryScreen() if QGuiApplication.instance() else None
    if screen is None:
        return 1.0, 1.0
    return screen.logicalDotsPerInch() / 72.0, max(1.0, screen.devicePixelRatio())


def _render(tex: str, colour: str, size: int) -> QImage | None:
    """The formula as vector graphics (matplotlib mathtext → SVG, glyphs as paths), drawn by
    Qt's SVG renderer at the screen's device pixel ratio: sharp on high-DPI screens, and
    ``size`` points as large as text of that point size. The image carries its device pixel
    ratio, so labels and text documents show it at its logical size."""
    try:
        from matplotlib import rc_context
        from matplotlib.backends.backend_svg import FigureCanvasSVG
        from matplotlib.figure import Figure
        from PyQt6.QtCore import QByteArray, QRectF
        from PyQt6.QtGui import QPainter
        from PyQt6.QtSvg import QSvgRenderer
    except ImportError:          # pragma: no cover
        return None
    try:
        with rc_context({"svg.fonttype": "path"}):
            fig = Figure(figsize=(0.01, 0.01), dpi=72)
            FigureCanvasSVG(fig)
            fig.text(0, 0, f"${_mathtext(tex)}$", fontsize=size, color=colour)
            buf = io.BytesIO()
            fig.savefig(buf, format="svg", transparent=True, bbox_inches="tight", pad_inches=0.02)
        renderer = QSvgRenderer(QByteArray(buf.getvalue()))
        box = renderer.viewBoxF()                      # points (dpi 72)
        if not renderer.isValid() or box.isEmpty():
            return None
        per_pt, dpr = _screen_scale()
        w, h = box.width() * per_pt, box.height() * per_pt
        img = QImage(max(1, round(w * dpr)), max(1, round(h * dpr)),
                     QImage.Format.Format_ARGB32_Premultiplied)
        img.fill(0)
        painter = QPainter(img)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        renderer.render(painter, QRectF(0, 0, img.width(), img.height()))
        painter.end()
        img.setDevicePixelRatio(dpr)
        return img
    except Exception:            # unsupported LaTeX: show the source instead
        return None


def _inline_md(s: str) -> str:
    s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    s = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", s)                 # links: their text
    s = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", s)
    s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
    return s


class PctFormulas(QFrame):
    def __init__(self, model, prestus_root: str | Path, parent=None):
        super().__init__(parent)
        self.setObjectName("Card")
        self.model = model
        self.docs = load(Path(prestus_root) / "doc")
        self._shown_for: tuple | None = None
        self._n = 0
        v = QVBoxLayout(self)
        v.setContentsMargins(14, 10, 14, 12)
        head = QHBoxLayout()
        t = QLabel("Mapping equations")
        t.setObjectName("CardTitle")
        head.addWidget(t, 1)
        # collapsed by default: the flow card above draws the mappings
        self.toggle = QToolButton()
        self.toggle.setCheckable(True)
        self.toggle.setToolTip("Show the equations of the selected mappings")
        self.toggle.toggled.connect(self._toggled)
        head.addWidget(self.toggle)
        v.addLayout(head)
        self.sub = QLabel("The equations of the selected mappings, from the PRESTUS documentation "
                          "(Skull Property Mappings, From UTE to Pseudo-HU).")
        self.sub.setObjectName("CardSub")
        self.sub.setWordWrap(True)
        v.addWidget(self.sub)
        self.view = QTextBrowser()
        self.view.setOpenLinks(False)
        self.view.setMinimumHeight(360)
        v.addWidget(self.view)
        self._toggled(False)
        self.refresh()

    def _toggled(self, on: bool) -> None:
        self.toggle.setIcon(icons.icon("collapse" if on else "expand"))
        self.sub.setVisible(on)
        self.view.setVisible(on)

    def _math(self, tex: str, display: bool) -> str:
        c = style.active()
        img = _render(tex, c["ink"], 12 if display else 10)
        if img is None:
            return f"<code>{_inline_md(tex)}</code>"
        self._n += 1
        url = QUrl(f"math://{self._n}")
        self.view.document().addResource(QTextDocument.ResourceType.ImageResource, url, img)
        dpr = img.devicePixelRatio()                       # logical size, so it is not drawn at 2×
        tag = (f"<img src='{url.toString()}' width='{round(img.width() / dpr)}' "
               f"height='{round(img.height() / dpr)}' style='vertical-align:middle'>")
        return f"<p align='center'>{tag}</p>" if display else tag

    def _body(self, text: str) -> str:
        """Markdown paragraphs with $$…$$ / $…$ maths → HTML."""
        html = []
        text = re.sub(r"```.*?\n(.*?)```", lambda m: "<pre>" + _inline_md(m.group(1)) + "</pre>", text, flags=re.S)
        for chunk in re.split(r"(\$\$.*?\$\$)", text, flags=re.S):
            if chunk.startswith("$$"):
                html.append(self._math(chunk.strip("$ \n"), True))
                continue
            for para in re.split(r"\n\s*\n", chunk):
                para = para.strip()
                if not para:
                    continue
                if para.startswith("<pre>"):
                    html.append(para)
                    continue
                pieces = re.split(r"(\$[^$]+\$)", para)
                html.append("<p>" + "".join(self._math(p.strip("$"), False) if p.startswith("$")
                                            else _inline_md(p) for p in pieces) + "</p>")
        return "".join(html)

    def _generation_html(self) -> str:
        """The UTE → pseudo-HU step of a pseudo-CT (path.ct_pattern empty)."""
        sm = str(self.model.get("pct.skull_mapping") or "")
        desc = self.docs["skull_mapping"].get(sm, "").replace(" (see below)", "")
        if sm == "custom":
            desc = (f"pHU = {self.model.get('pct.custom_slope')} UTE + "
                    f"{self.model.get('pct.custom_intercept')}")
        return (f"<h3>Pseudo-CT, UTE → pseudo-HU: <code>{sm}</code></h3><p>{_inline_md(desc) or '(not documented)'}</p>"
                "<p>Soft tissue is set to 42 HU and air to −1000 HU for every mapping.</p>")

    def refresh(self) -> None:
        keys = ["pct.enabled", "path.ct_pattern", "pct.skull_mapping", *(k for k, _, _ in PROPERTIES),
                "pct.custom_slope", "pct.custom_intercept"]
        sig = tuple(repr(self.model.get(k)) for k in keys) + (style.active()["ink"],)
        if sig == self._shown_for:
            return
        self._shown_for = sig
        self.view.clear()
        self._n = 0
        if self.model.get("pct.enabled") not in (1, True):
            self.view.setHtml("<p>CT-based skull mapping is off: the skull keeps the uniform skull "
                              "properties of the Acoustic properties table. Turn on <b>Enabled</b> "
                              "to map skull properties from a measured CT or a pseudo-CT.</p>")
            return
        parts = []
        if not is_empty(self.model.get("path.ct_pattern")):
            parts.append("<h3>CT: measured</h3><p>The HU of the measured CT (path.ct_pattern) are "
                         "mapped directly; no UTE → pseudo-HU step.</p>")
        else:
            parts.append(self._generation_html())
        for key, label, _ in PROPERTIES:
            value = str(self.model.get(key) or "none")
            parts.append(f"<h3>{label}: <code>{value}</code></h3>")
            if value == "none":
                parts.append("<p>Not mapped: the skull keeps its uniform value.</p>")
                continue
            desc, eq = self.docs.get(key, {}).get(value, ("", ""))
            if desc:
                parts.append(f"<p>{_inline_md(desc)}</p>")
            parts.append(self._body(eq) if eq else "<p>(no equation in the documentation)</p>")
        self.view.setHtml("".join(parts))
