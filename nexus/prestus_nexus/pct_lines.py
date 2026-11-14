"""The UTE -> pseudo-HU step of the pseudo-CT, read from functions/pct/pct_create_pseudoCT.sh
(the code that applies it), so the Nexus never restates the coefficients:

- the linear map of each pct.skull_mapping, from the fslmaths chain of its 'case' branch
  (e.g. miscouridou: fslmaths "$ute_norm" -mul -2085 -add 2329 ...) as HU = slope * UTE + intercept
  on the soft-tissue-normalised UTE;
- 'custom': pct.custom_slope / pct.custom_intercept;
- 'kosciessa': fitted per subject; the fit is read from <pct folder>/debug_pct/pCT_skull_mapping.txt
  (slope on line 1, intercept on line 2) once a pseudo-CT has been generated;
- the fixed soft-tissue and air values (pCT_soft_val, pCT_air_val).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

SCRIPT = Path("functions") / "pct" / "pct_create_pseudoCT.sh"


@dataclass
class UteLine:
    mapping: str
    slope: float | None            # HU per normalised UTE unit; None: not known yet (kosciessa before a fit)
    intercept: float | None
    source: str                    # where the numbers come from, for the caption


def _fslmaths_linear(ops: str) -> tuple[float, float] | None:
    """slope, intercept of an fslmaths chain on x made of -mul/-add/-sub/-div steps."""
    a, b = 1.0, 0.0
    for op, val in re.findall(r"-(mul|add|sub|div)\s+(-?[0-9.]+)", ops):
        v = float(val)
        if op == "mul":
            a, b = a * v, b * v
        elif op == "add":
            b += v
        elif op == "sub":
            b -= v
        else:
            a, b = a / v, b / v
    return a, b


def read_script(prestus_root: str | Path) -> dict:
    """{'lines': {mapping: (slope, intercept)}, 'soft': HU, 'air': HU} from the script."""
    text = (Path(prestus_root) / SCRIPT).read_text(encoding="utf-8", errors="replace")
    lines: dict[str, tuple[float, float]] = {}
    # a branch:  name)\n  fslmaths "$ute_norm" -mul ... -mul "$skull_mask_PVC_dil" "$pCT_skull" ;;
    for name, ops in re.findall(r'^\s*([a-z]+)\)\s*\n\s*fslmaths "\$ute_norm"((?:\s+-(?:mul|add|sub|div)\s+-?[0-9.]+)+)',
                                text, flags=re.M):
        lines[name] = _fslmaths_linear(ops)
    soft = re.search(r"pCT_soft_val=(-?[0-9.]+)", text)
    air = re.search(r"pCT_air_val=(-?[0-9.]+)", text)
    return {"lines": lines, "soft": float(soft.group(1)) if soft else None,
            "air": float(air.group(1)) if air else None}


def ute_line(model, prestus_root: str | Path, pct_dir: Path | None) -> UteLine:
    """The UTE -> HU line the configured pct.skull_mapping applies."""
    mapping = str(model.get("pct.skull_mapping") or "kosciessa")
    if mapping == "custom":
        s, i = model.get("pct.custom_slope"), model.get("pct.custom_intercept")
        ok = all(isinstance(x, (int, float)) and not isinstance(x, bool) for x in (s, i))
        return UteLine(mapping, float(s) if ok else None, float(i) if ok else None,
                       "pct.custom_slope / pct.custom_intercept")
    if mapping == "kosciessa":
        fit = pct_dir / "debug_pct" / "pCT_skull_mapping.txt" if pct_dir is not None else None
        if fit is not None and fit.is_file():
            try:
                vals = [float(ln.split()[0]) for ln in fit.read_text().splitlines() if ln.strip()][:2]
                return UteLine(mapping, vals[0], vals[1], f"this subject's fit ({fit.name})")
            except (ValueError, IndexError):
                pass
        return UteLine(mapping, None, None, "fitted per subject from the UTE histograms when the pseudo-CT is made")
    try:
        line = read_script(prestus_root)["lines"].get(mapping)
    except OSError:
        line = None
    if line is None:
        return UteLine(mapping, None, None, f"not found in {SCRIPT.name}")
    return UteLine(mapping, line[0], line[1], SCRIPT.name)
