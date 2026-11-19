#!/usr/bin/env python3
"""Generate doc/doc_functions.md from the MATLAB sources in functions/.

One table per subsystem folder (functions/<folder>/, recursively), one row per
function file. Each row's description is the first help line of its file
('% NAME  Description'), so a description is edited in the source, never in
doc/doc_functions.md. Only the intro block (everything before the first
'#### ' heading) is edited by hand.

Usage:
    python tests/docs/gen_function_index.py           # rewrite doc/doc_functions.md
    python tests/docs/gen_function_index.py --check   # exit 1 if it is out of date
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FUNCTIONS = ROOT / "functions"
DOC = ROOT / "doc" / "doc_functions.md"

# Section titles that differ from the upper-cased folder name
TITLES = {"pct": "PSEUDO-CT", "babelbrain": "BABELBRAIN", "jwave": "JWAVE"}
# Files at the top of functions/ are listed under this section
TOP_LEVEL_SECTION = "core"
HEADER = ("| **Function Name**                       | **Description**"
          "                                                                                 |\n"
          "|-----------------------------------------|-------------------------------------"
          "------------------------------------------------------------|")


def help_line(path: Path) -> str:
    """First help-text line of a MATLAB file, without the function name."""
    name = path.stem
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines()[:40]:
        line = raw.strip()
        if not line.startswith("%"):
            continue
        body = line.lstrip("%").strip()
        if not body:
            continue
        # '% NAME  Description' / '% NAME - Description' / '% Description'
        m = re.match(rf"^{re.escape(name)}\b[\s:\-–—]*(.*)$", body, flags=re.IGNORECASE)
        desc = m.group(1).strip() if m else body
        return desc.replace("|", "\\|") or "(no description)"
    return "(no description)"


def is_function_file(path: Path) -> bool:
    """MATLAB function files (not scripts, classdefs are listed as well)."""
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("%"):
            continue
        return line.startswith("function") or line.startswith("classdef")
    return False


def collect() -> dict[str, list[Path]]:
    groups: dict[str, list[Path]] = {}
    for path in sorted(FUNCTIONS.rglob("*.m")):
        rel = path.relative_to(FUNCTIONS)
        if any(part.startswith(".") or part in ("private", "__pycache__") for part in rel.parts[:-1]):
            continue
        if not is_function_file(path):
            continue
        group = rel.parts[0] if len(rel.parts) > 1 else TOP_LEVEL_SECTION
        groups.setdefault(group, []).append(path)
    return groups


def render(intro: str) -> str:
    parts = [intro.rstrip("\n"), ""]
    for group, paths in sorted(collect().items()):
        parts.append(f"#### {TITLES.get(group, group.upper())}")
        parts.append("")
        parts.append(HEADER)
        for path in sorted(paths, key=lambda p: p.stem.lower()):
            name = path.stem
            desc = help_line(path)
            parts.append(f"| `{name}`".ljust(42) + f"| {desc} |")
        parts.append("")
    return "\n".join(parts).rstrip("\n") + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="fail if doc/doc_functions.md is out of date")
    args = ap.parse_args()

    text = DOC.read_text(encoding="utf-8")
    first = text.find("\n#### ")
    intro = text if first < 0 else text[: first + 1]
    new = render(intro)
    if args.check:
        if new != text:
            print("doc/doc_functions.md is out of date; run: python tests/docs/gen_function_index.py")
            return 1
        print("doc/doc_functions.md is up to date")
        return 0
    DOC.write_text(new, encoding="utf-8")
    print(f"wrote {DOC.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
