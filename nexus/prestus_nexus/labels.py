"""Readable field labels for schema nodes.

An explicit ``label`` from overlay.yaml wins; otherwise the key is turned
into words (``elem_od_mm`` -> "Element outer diameter"), with units and a
leading "is"/"do" dropped (the editor shows the unit). Abbreviations are
spelled out or kept in their usual capitals. The full help text is shown on
hover.
"""

from __future__ import annotations

import re

from .schema import Node

_WORDS = {
    "elem": "element", "n": "number of", "id": "inner diameter", "od": "outer diameter",
    "pos": "position", "freq": "frequency", "ep": "exit plane", "dist": "distance",
    "geom": "geometric", "curv": "curvature", "p": "pressure", "amp": "amplitude",
    "isppa": "ISPPA", "ispta": "ISPTA", "mni": "MNI", "ras": "RAS", "t1": "T1", "t2": "T2",
    "pct": "pseudo-CT", "fov": "FOV", "pml": "PML", "ppw": "PPW", "cfl": "CFL", "gpu": "GPU",
    "cpu": "CPU", "hpc": "HPC", "roi": "ROI", "csf": "CSF", "fwhm": "FWHM", "flhm": "FLHM",
    "tx": "transducer", "ds": "driving system", "roc": "radius of curvature", "wb": "Workbench",
    "cem43": "CEM43", "gmres": "GMRES", "maxiter": "max. iterations", "q": "Q", "xy": "xy", "kwavearray": "kWaveArray", "simnibs": "SimNIBS",
    "ct": "CT", "hu": "HU", "ute": "UTE", "io": "input/output", "col": "column", "rot": "rotation",
    "rect": "rectangular", "posthoc": "post-hoc", "timelimit": "time limit", "memorylimit": "memory limit",
    "fd": "finite-difference", "tol": "tolerance", "coeff": "coefficient", "long": "longitudinal",
    "az": "azimuthal", "sep": "separation", "env": "environment", "exclu": "exclusion",
    "soundspeed": "sound speed", "preproc": "preprocessing", "min": "min.", "max": "max.",
    "iso": "(ISO)", "rayleigh": "Rayleigh", "poisson": "Poisson", "babelbrain": "BabelBrain", "jwave": "jwave",
}
_UNITS = {"mm", "cm", "hz", "s", "us", "deg", "rad", "wcm2", "pa", "ms", "mpa", "c", "gb"}
_LEAD = {"is", "do"}          # 'is_curved' -> "Curved", 'do_refocusing' -> "Refocusing"


def plural(n: int, word: str, many: str | None = None) -> str:
    """'1 error', '3 errors' (`many` for irregular plurals: plural(n, 'is', 'are'))."""
    return f"{n} {word if n == 1 else (many or word + 's')}"


def humanize(key: str) -> str:
    words = [w for w in re.split(r"_+", key) if w]
    if len(words) > 1 and words[-1].lower() in _UNITS:
        words = words[:-1]
    if len(words) > 1 and words[0].lower() in _LEAD:
        words = words[1:]
    out = []
    for w in words:
        lw = w.lower()
        out.append(_WORDS.get(lw, w if any(c.isupper() for c in w[1:]) else lw))
    text = " ".join(out).replace("number of number", "number")
    text = re.sub(r"^element number of$", "number of elements", text)
    return text[:1].upper() + text[1:] if text else key


def label(node: Node) -> str:
    """Plain label (combo items, file names, Markdown, window titles): 'ISPPA', 'ISPTA'."""
    return node.label or humanize(node.name)


def choice(node: Node | None, value) -> str:
    """Drop-down text of a choice value: its overlay 'choice_labels' entry, else the value."""
    labels = node.choice_labels if node is not None else {}
    return str(labels.get(value, value))


# ITRUSST notation: the intensity acronyms as I with a subscript where text is rich. A
# pipeline name ('multi-ISPPA') or a file name keeps the plain form.
_INTENSITY = re.compile(r"(?<![\w-])I(SPPA|SPTA)\b")


def rich(text: str) -> str:
    """Qt rich text with I<sub>SPPA</sub> / I<sub>SPTA</sub>. `text` must already be HTML
    (escaped where needed); plain text goes through rich_label/html.escape first."""
    return _INTENSITY.sub(r"I<sub>\1</sub>", text)


def rich_label(node: Node) -> str:
    """The label as escaped rich text with the intensity subscripts (form labels, tooltips)."""
    import html
    return rich(html.escape(label(node), quote=False))


def mathtext(text: str) -> str:
    """Matplotlib label with the intensity subscripts as mathtext: I$_{\\mathrm{SPPA}}$."""
    return _INTENSITY.sub(r"I$_{\\mathrm{\1}}$", text)

