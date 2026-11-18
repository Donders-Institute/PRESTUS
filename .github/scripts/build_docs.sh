#!/usr/bin/env bash
# Build the docs site of one PRESTUS checkout for the versioned site.
#
#   build_docs.sh <checkout> <out_dir> [banner]
#
# A checkout whose requirements-docs.txt lists zensical is built with Zensical,
# any other (main, release tags up to v0.6.x) with MkDocs + Material. Both read
# versions.json for the version selector (extra.version, added here when the
# config has none). [banner] is the text shown on every page of a version other
# than the default (main): it goes into extra.version_banner, and checkouts
# without theme overrides get a minimal override that shows it. All changes are
# made to the given checkout only (a temporary worktree in the workflow).
# The tools (zensical, mkdocs, mkdocs-material, ...) must be on PATH.
set -euo pipefail

src=$1
out=$2
banner=${3:-}
here=$(cd "$(dirname "$0")" && pwd)

cd "$src"
if [ -f mkdocs.yaml ]; then cfg=mkdocs.yaml
elif [ -f mkdocs.yml ]; then cfg=mkdocs.yml
else echo "no mkdocs.yaml in $src" >&2; exit 1; fi

mkdir -p _version_overrides
cp "$here/version_banner.html" _version_overrides/main.html
python3 "$here/docs_version_config.py" "$cfg" "$banner"

if grep -qi zensical requirements-docs.txt 2>/dev/null; then
  zensical build --clean -f "$cfg"       # writes site/ inside the checkout
  rm -rf "$out"
  cp -r site "$out"
else
  mkdocs build --clean -f "$cfg" -d "$out"
fi
echo "built $(find "$out" -name '*.html' | wc -l | tr -d ' ') pages from $src into $out"
