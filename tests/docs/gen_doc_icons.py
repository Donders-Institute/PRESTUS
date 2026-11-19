#!/usr/bin/env python3
"""Write the option-table icons of the docs from the Nexus drawings.

The docs show the same pictures as the PRESTUS Nexus drop-downs: every glyph of
nexus/prestus_nexus/ui/choice_icons.py becomes doc/assets/icons/icon_<name>_glyph.svg
(24 x 24, text SVG), coloured with the docs palette (ink #6e7781, accent #0e7490; the
dark scheme inverts them through stylesheets/extra.css, class "choice-icon").
choice_icons.py is plain Python, so no Qt is needed.

Usage:
    python tests/docs/gen_doc_icons.py           # (re)write doc/assets/icons/
    python tests/docs/gen_doc_icons.py --check   # exit 1 if an icon is missing or out of date
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "nexus" / "prestus_nexus" / "ui" / "choice_icons.py"
OUT = ROOT / "doc" / "assets" / "icons"
INK = "#6e7781"      # secondary text (docs figure palette)
ACCENT = "#0e7490"   # PRESTUS/Nexus teal: what the choice changes


def load_glyphs():
    spec = importlib.util.spec_from_file_location("choice_icons", SOURCE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def render(mod) -> dict[str, str]:
    files = {}
    for name in mod.GLYPHS:
        body = mod.glyph_body(name).replace("currentColor", INK).replace("accentColor", ACCENT)
        files[f"icon_{name.replace('_', '-')}_glyph.svg"] = (
            '<svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" '
            f'fill="none" stroke="{INK}" stroke-width="1.8" stroke-linecap="round" '
            f'stroke-linejoin="round">{body}</svg>\n')
    return files


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="only check; exit 1 if out of date")
    args = ap.parse_args()
    files = render(load_glyphs())
    stale = [n for n, s in files.items()
             if not (OUT / n).is_file() or (OUT / n).read_text() != s]
    if args.check:
        for n in stale:
            print(f"out of date: doc/assets/icons/{n}")
        return 1 if stale else 0
    OUT.mkdir(parents=True, exist_ok=True)
    for n in stale:
        (OUT / n).write_text(files[n])
    print(f"{len(files)} icons in {OUT.relative_to(ROOT)} ({len(stale)} written)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
