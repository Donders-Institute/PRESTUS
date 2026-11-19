#!/usr/bin/env python3
"""Check that every image and PRESTUS_bin file the documentation uses exists.

Two kinds of reference are checked:

- local images, ![alt](assets/...) in doc/*.md: the file must exist inside doc/;
- files hosted in the PRESTUS_bin repository (any URL into
  github.com/jkosciessa/PRESTUS_bin or raw.githubusercontent.com/jkosciessa/PRESTUS_bin)
  in doc/*.md, doc/stylesheets/*.css, and mkdocs.yaml (logo, favicon).

Hosted files are looked up in the committed tree (HEAD) of a local PRESTUS_bin clone with
--bin DIR, so a file not yet committed counts as missing; without it, each URL is requested
over HTTP (as in CI).

Usage:
    python tests/docs/check_images.py
    python tests/docs/check_images.py --bin ../PRESTUS_bin
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DOC = ROOT / "doc"

IMAGE = re.compile(r"!\[[^\]]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")
BIN_URL = re.compile(
    r"https://(?:github\.com/jkosciessa/PRESTUS_bin/raw|raw\.githubusercontent\.com/jkosciessa/PRESTUS_bin)"
    r"/(?:refs/heads/)?main/([^)\s\"'<>]+)")
RAW = "https://raw.githubusercontent.com/jkosciessa/PRESTUS_bin/main/"


def sources() -> list[Path]:
    return sorted(DOC.glob("*.md")) + sorted((DOC / "stylesheets").glob("*.css")) + [ROOT / "mkdocs.yaml"]


def exists_remote(path: str) -> bool:
    req = urllib.request.Request(RAW + path, method="HEAD", headers={"User-Agent": "prestus-doc-check"})
    try:
        with urllib.request.urlopen(req, timeout=20) as r:
            return r.status == 200
    except urllib.error.HTTPError:
        return False


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--bin", type=Path, help="local PRESTUS_bin clone to check hosted files against")
    args = ap.parse_args()

    problems: list[str] = []
    hosted: dict[str, list[str]] = {}
    for src in sources():
        for n, line in enumerate(src.read_text(encoding="utf-8").splitlines(), 1):
            where = f"{src.relative_to(ROOT)}:{n}"
            for path in BIN_URL.findall(line):
                hosted.setdefault(path, []).append(where)
            if src.suffix != ".md":
                continue
            for target in IMAGE.findall(line):
                if re.match(r"^[a-z][a-z0-9+.\-]*:", target, flags=re.IGNORECASE):
                    continue                               # external: hosted ones checked above
                dest = (src.parent / target.split("#")[0]).resolve()
                if not dest.is_file():
                    problems.append(f"{where}: missing image {target}")

    committed: set[str] = set()
    if args.bin:
        committed = set(subprocess.run(["git", "-C", str(args.bin), "ls-tree", "-r", "--name-only", "HEAD"],
                                       capture_output=True, text=True, check=True).stdout.splitlines())
    for path, wheres in sorted(hosted.items()):
        ok = path in committed if args.bin else exists_remote(path)
        if not ok:
            problems.append(f"{wheres[0]}: PRESTUS_bin has no {path}"
                            + (f" (also {len(wheres) - 1} more use(s))" if len(wheres) > 1 else ""))

    if problems:
        print(f"{len(problems)} missing image(s) or file(s):")
        print("\n".join(f"  {p}" for p in problems))
        return 1
    print(f"all images and files resolve ({len(hosted)} hosted in PRESTUS_bin"
          f"{', checked against ' + str(args.bin) if args.bin else ''})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
