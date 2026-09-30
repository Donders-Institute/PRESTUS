"""Readable field labels for schema nodes.

An explicit ``label`` from overlay.yaml wins; otherwise the key is turned
into words (``elem_od_mm`` -> "Element outer diameter"), with units dropped
because the editor shows them. The full help text is shown on hover.
"""

from __future__ import annotations

import re

from .schema import Node

_WORDS = {
    "elem": "element", "n": "number of", "id": "inner diameter", "od": "outer diameter",
    "pos": "position", "freq": "frequency", "ep": "exit plane", "dist": "distance",
    "geom": "geometric", "curv": "curvature", "p": "pressure", "amp": "amplitude",
    "isppa": "Isppa", "ispta": "Ispta", "mni": "MNI", "ras": "RAS", "t1": "T1", "t2": "T2",
    "pct": "pseudo-CT", "fov": "FOV", "pml": "PML", "ppw": "PPW", "cfl": "CFL", "gpu": "GPU",
    "cpu": "CPU", "hpc": "HPC", "roi": "ROI", "csf": "CSF", "fwhm": "FWHM", "flhm": "FLHM",
    "tx": "transducer", "ds": "driving system", "roc": "radius of curvature", "wb": "Workbench",
    "cem43": "CEM43", "gmres": "GMRES", "maxiter": "max. iterations", "q": "Q", "xy": "xy", "kwavearray": "kWaveArray", "simnibs": "SimNIBS",
}
_UNITS = {"mm", "hz", "s", "deg", "rad", "wcm2", "pa", "ms", "mpa", "c", "gb"}


def humanize(key: str) -> str:
    words = [w for w in re.split(r"_+", key) if w]
    if len(words) > 1 and words[-1].lower() in _UNITS:
        words = words[:-1]
    out = []
    for w in words:
        lw = w.lower()
        out.append(_WORDS.get(lw, w if any(c.isupper() for c in w[1:]) else lw))
    text = " ".join(out).replace("number of number", "number")
    return text[:1].upper() + text[1:] if text else key


def label(node: Node) -> str:
    return node.label or humanize(node.name)

