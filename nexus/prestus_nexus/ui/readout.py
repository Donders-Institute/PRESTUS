"""Readouts: short titled sections of label/value rows, as rich text for a QLabel.

One layout for the side panels that report measured or resolved values (placement viewer,
transducer panel), so they scan the same way: muted labels in a left column, values on the
right, status words coloured from the active palette."""

from __future__ import annotations

from html import escape

from . import style

Row = tuple[str, str]                 # (label, value as rich text)
Section = tuple[str, list[Row]]       # (title, rows); a section without rows is skipped


GLYPH = {"ok": "✓", "warn": "▲", "bad": "✗"}   # status shape, so colour is never the only cue


def status(text: str, level: str) -> str:
    """text coloured for level 'ok' | 'warn' | 'bad' | 'muted'; the first three also lead
    with a glyph, readable with red-green colour deficits and in greyscale."""
    c = style.active()
    colour = {"ok": c["green"], "warn": c["amber"], "bad": c["red"], "muted": c["ink2"]}[level]
    glyph = GLYPH.get(level)
    return f"<span style='color:{colour}'>{glyph + ' ' if glyph else ''}{text}</span>"


def muted(text: str) -> str:
    return status(text, "muted")


def readout(sections: list[Section], footer: str = "", label_width: int | None = None) -> str:
    """Sections in one two-column table (so the columns line up across sections), each
    under a small capitalised title row; footer as a muted last line. label_width [px] fixes
    the label column, e.g. to line the values up with the inputs of a form."""
    c = style.active()
    out = ["<table width='100%' cellspacing='0' cellpadding='2'>"]
    first = True
    for title, rows in sections:
        if not rows:
            continue
        pad = "2px" if first else "12px"
        first = False
        out.append(f"<tr><td colspan='2' style='padding-top:{pad}; color:{c['ink2']}; "
                   f"font-weight:600; font-size:small'>{escape(title.upper())}</td></tr>")
        for label, value in rows:
            width = f" width='{label_width}'" if label_width else ""
            out.append(f"<tr><td nowrap{width} style='color:{c['ink2']}; padding-right:10px' valign='top'>{escape(label)}</td>"
                       f"<td valign='top'>{value}</td></tr>")
    out.append("</table>")
    if footer:
        out.append(f"<p style='margin-top:10px; color:{c['ink2']}; font-size:small'>{footer}</p>")
    return "".join(out)
