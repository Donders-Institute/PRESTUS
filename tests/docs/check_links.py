#!/usr/bin/env python3
"""Check relative links in the documentation.

Scans doc/*.md for Markdown links [text](target) with a relative target and
checks that the file exists, lies inside doc/ (MkDocs serves only those), and,
for 'file.md#anchor' / '#anchor', that the
anchor matches a heading of that file (MkDocs' default slug rule: lower case,
punctuation removed, spaces to hyphens) or an explicit <a id/name="..."> tag.
External links (http, https, mailto) are not checked.

Usage:
    python tests/docs/check_links.py
"""
from __future__ import annotations

import re
import sys
import unicodedata
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DOC = ROOT / "doc"

LINK = re.compile(r"(?<!!)\[[^\]]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")
HEADING = re.compile(r"^(#{1,6})\s+(.*?)\s*#*\s*$")
EXPLICIT = re.compile(r"""<a\s+(?:id|name)=["']([^"']+)["']""")
ATTR_LIST = re.compile(r"\{\s*([^}]*)\}\s*$")         # attr_list: { #id .class key=value }
FENCE = re.compile(r"^\s*(```|~~~)")


def slugify(text: str) -> str:
    """Python-Markdown toc default slugify (as used by MkDocs)."""
    text = re.sub(r"<[^>]+>", "", text)                       # inline HTML
    text = re.sub(r"`([^`]*)`", r"\1", text)                  # code spans
    text = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", text)      # links
    text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode("ascii")
    text = re.sub(r"[^\w\s-]", "", text).strip().lower()
    return re.sub(r"[-\s]+", "-", text)


def anchors(path: Path) -> set[str]:
    found: set[str] = set()
    seen: dict[str, int] = {}
    in_fence = False
    for line in path.read_text(encoding="utf-8").splitlines():
        if FENCE.match(line):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        found.update(EXPLICIT.findall(line))
        m = HEADING.match(line)
        if m:
            title = m.group(2)
            attr = ATTR_LIST.search(title)
            if attr:
                found.update(t[1:] for t in attr.group(1).split() if t.startswith("#"))
                title = ATTR_LIST.sub("", title)
            slug = slugify(title)
            if slug in seen:                       # duplicates get _1, _2, ...
                seen[slug] += 1
                found.add(f"{slug}_{seen[slug]}")
            else:
                seen[slug] = 0
                found.add(slug)
    return found


def links(path: Path):
    in_fence = False
    for n, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if FENCE.match(line):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        for target in LINK.findall(line):
            yield n, target


def main() -> int:
    cache: dict[Path, set[str]] = {}
    problems: list[str] = []
    files = sorted(DOC.glob("*.md"))
    for md in files:
        for n, target in links(md):
            if re.match(r"^[a-z][a-z0-9+.\-]*:", target, flags=re.IGNORECASE):
                continue                                   # http:, mailto:, ...
            file_part, _, anchor = target.partition("#")
            dest = md if not file_part else (md.parent / file_part).resolve()
            where = f"{md.relative_to(ROOT)}:{n}"
            if not dest.exists():
                problems.append(f"{where}: missing target {target}")
                continue
            if DOC.resolve() not in dest.resolve().parents and dest.resolve() != DOC.resolve():
                problems.append(f"{where}: {target} is outside doc/ (MkDocs cannot serve it; "
                                "use a code span or a GitHub URL)")
                continue
            if anchor and dest.suffix == ".md":
                if dest not in cache:
                    cache[dest] = anchors(dest)
                if anchor not in cache[dest]:
                    problems.append(f"{where}: no heading for #{anchor} in {dest.name}")
    if problems:
        print(f"{len(problems)} broken link(s):")
        print("\n".join(f"  {p}" for p in problems))
        return 1
    print(f"all relative links in {len(files)} files resolve")
    return 0


if __name__ == "__main__":
    sys.exit(main())
