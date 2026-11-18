#!/usr/bin/env python3
"""Place one built docs version on the gh-pages tree (versioned site, mike layout).

The docs workflow (.github/workflows/mkdocs.yml) builds one version with
build_docs.sh and calls this with a checkout of the gh-pages branch. The tree
then holds

    main/                 the main branch (default)
    development/          the development branch
    v0.6.1/ ...           release tags, frozen as tagged
    versions.json         the list the theme's version selector reads (mike format)
    index.html, <page>/   redirects into the default version, so links to the
                          unversioned site (https://.../PRESTUS/doc_hpc/) keep working
    .nojekyll

Usage:
    python .github/scripts/deploy_docs_version.py --site <built site> --pages <gh-pages tree> \
        --version main --title "main (v0.6.1)" [--default main]
"""

import argparse
import json
import re
import shutil
from pathlib import Path

ROOT_MANIFEST = ".root-redirects"  # root redirect stubs written by the last default deploy

STUB = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Redirecting</title>
<link rel="canonical" href="{target}">
<meta http-equiv="refresh" content="0; url={target}">
<script>location.replace("{target}" + location.search + location.hash)</script>
</head>
<body><a href="{target}">Continue to the documentation</a></body>
</html>
"""


def load_versions(pages: Path) -> list:
    f = pages / "versions.json"
    if f.is_file():
        try:
            data = json.loads(f.read_text())
            if isinstance(data, list):
                return data
        except json.JSONDecodeError:
            pass
    return []


def upsert(versions: list, version: str, title: str, aliases: list, pages: Path) -> list:
    out = []
    for v in versions:
        if v.get("version") == version:
            # aliases dropped since the last deploy (e.g. mike's latest/dev) leave no folder behind
            for a in set(v.get("aliases", [])) - set(aliases):
                if a and (pages / a).is_dir() and a != version:
                    shutil.rmtree(pages / a)
            continue
        out.append(v)
    out.append({"version": version, "title": title, "aliases": aliases})
    return out


def version_key(name: str):
    """Release tags newest first: v0.10.0 before v0.6.1 before v0.6.0."""
    nums = [int(n) for n in re.findall(r"\d+", name)]
    return [-n for n in nums]


def order_versions(versions: list, first: list) -> list:
    """The names in first (main, development) lead; release tags follow, newest first."""
    rank = {name: i for i, name in enumerate(first)}
    return sorted(versions, key=lambda v: (rank.get(v["version"], len(rank)), version_key(v["version"]), v["version"]))


def write_root_redirects(site: Path, pages: Path, default: str, protected: set) -> int:
    """One stub per page of the default version, at its unversioned path."""
    manifest = pages / ROOT_MANIFEST
    if manifest.is_file():
        for rel in manifest.read_text().split("\n"):
            stub = pages / rel
            if rel and stub.is_file():
                stub.unlink()
                d = stub.parent
                while d != pages and d.is_dir() and not any(d.iterdir()):
                    d.rmdir()
                    d = d.parent
    written = []
    for page in sorted(site.rglob("index.html")):
        rel = page.relative_to(site)
        if rel.parts[0] in protected:
            continue
        depth = len(rel.parts) - 1
        target = "../" * depth + default + "/" + "".join(p + "/" for p in rel.parts[:-1])
        dest = pages / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(STUB.format(target=target))
        written.append(rel.as_posix())
    manifest.write_text("\n".join(written) + "\n")
    return len(written)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--site", required=True, type=Path, help="built site (zensical build)")
    ap.add_argument("--pages", required=True, type=Path, help="checkout of the gh-pages branch")
    ap.add_argument("--version", required=True, help="folder name, e.g. main or development")
    ap.add_argument("--title", default=None, help="label in the version selector")
    ap.add_argument("--alias", action="append", default=[], help="extra name in versions.json")
    ap.add_argument("--default", default="main", help="version the site root redirects to")
    ap.add_argument("--order", nargs="*", default=["main", "development"], help="selector order")
    a = ap.parse_args()

    site, pages, version = a.site.resolve(), a.pages.resolve(), a.version
    if not (site / "index.html").is_file():
        raise SystemExit(f"no built site in {site}")
    if not version or "/" in version or version.startswith("."):
        raise SystemExit(f"invalid version name {version!r}")
    pages.mkdir(parents=True, exist_ok=True)

    dest = pages / version
    if dest.exists():
        shutil.rmtree(dest)
    shutil.copytree(site, dest)

    versions = upsert(load_versions(pages), version, a.title or version, a.alias, pages)
    versions = order_versions(versions, a.order)
    (pages / "versions.json").write_text(json.dumps(versions, indent=2) + "\n")
    (pages / ".nojekyll").touch()

    msg = f"deployed {version} ({sum(1 for _ in dest.rglob('*.html'))} pages)"
    if version == a.default:
        protected = {v["version"] for v in versions} | {x for v in versions for x in v.get("aliases", [])}
        n = write_root_redirects(site, pages, a.default, protected)
        msg += f"; {n} root redirects into {a.default}/"
    elif not (pages / "index.html").is_file():
        # first deploy is not the default: the root still needs an entry point
        (pages / "index.html").write_text(STUB.format(target=a.default + "/"))
    print(msg)
    print("versions:", ", ".join(v["version"] for v in versions))


if __name__ == "__main__":
    main()
