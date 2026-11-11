#!/usr/bin/env python3
"""Check that every key of config/config_default.yaml is documented in doc/doc_parameters.md.

A key counts as documented when its dotted path (e.g. `io.overwrite_files`), or
a trailing part of it in backticks (e.g. `overwrite_files`, `subset.random_seed`),
appears in doc/doc_parameters.md. Keys of list entries (transducer elements) are checked
through the first entry.

The Default column of the parameter tables is also compared with the value in
config_default.yaml, for rows under a heading that names the section in
backticks (e.g. "### I/O management (`io`)" or "#### `placement.mtnu`") and
whose Default cell is a single code literal (`0`, `'auto'`, `[1, 2]`, `[] → 1`
reads as `[]`). Tables inside collapsible blocks (indented by four spaces) are
checked as well; there a bold-only line naming the section in backticks
(e.g. "**`placement.manual`**") sets the section like a heading. Exit status 1 lists undocumented keys and differing defaults.

Usage:
    python tests/docs/check_param_docs.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / "config" / "config_default.yaml"
DOC = ROOT / "doc" / "doc_parameters.md"


def leaf_keys(node, prefix: str = "") -> list[str]:
    keys: list[str] = []
    if isinstance(node, dict):
        for key, value in node.items():
            path = f"{prefix}.{key}" if prefix else str(key)
            if isinstance(value, dict) and value:
                keys.extend(leaf_keys(value, path))
            elif isinstance(value, list) and value and isinstance(value[0], dict):
                keys.extend(leaf_keys(value[0], path))
            else:
                keys.append(path)
    return keys


_MISSING = object()


def lookup(config, path: str):
    node = config
    for part in path.split("."):
        if isinstance(node, list) and node:
            node = node[0]
        if not isinstance(node, dict) or part not in node:
            return _MISSING
        node = node[part]
    return node


def normalise(value) -> str:
    """Canonical text of a YAML value."""
    if isinstance(value, bool):
        return str(int(value))
    if value is None:
        return "[]"
    if isinstance(value, (int, float)):
        return repr(float(value))
    if isinstance(value, list):
        return "[" + ",".join(normalise(v) for v in value) + "]"
    text = str(value).strip()
    if text.lower() in ("true", "false"):
        return str(int(text.lower() == "true"))
    try:
        return repr(float(text))
    except ValueError:
        return text


def documented_default(cell: str):
    """Canonical text of a Default cell, or None when it is not a plain code literal."""
    cell = cell.strip()
    m = re.fullmatch(r"`([^`]*)`(\s*→.*)?", cell)
    if not m:
        return None
    text = m.group(1).strip().strip("'\"")
    if text.startswith("[") and text.endswith("]"):
        items = [x.strip().strip("'\"") for x in re.split(r"[,;]", text[1:-1]) if x.strip()]
        return "[" + ",".join(normalise(x) for x in items) + "]"
    if text in ("''", ""):
        return ""
    if text in ("~", "null"):
        return "[]"
    return normalise(text)


def default_mismatches(config, doc: str) -> list[str]:
    out: list[str] = []
    prefix = None
    default_col = None
    for line in doc.splitlines():
        # Folded content (`??? note` blocks) is indented. A bold-only line
        # (e.g. "**`placement.manual`**") acts as a heading.
        line = line.lstrip()
        heading = re.match(r"^#{2,6}\s", line) or re.fullmatch(r"\*\*[^*]+\*\*:?", line.rstrip())
        if heading:
            m = (re.search(r"\(`([\w.]+)`", line) or re.match(r"^#{2,6}\s+`([\w.]+)`", line)
                 or re.match(r"^\*\*`([\w.]+)`", line))
            prefix = m.group(1) if m else None
            default_col = None
            continue
        if not prefix or not line.startswith("|"):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if any("Default" in c for c in cells):
            default_col = next(i for i, c in enumerate(cells) if "Default" in c)
            continue
        if default_col is None or len(cells) <= default_col:
            continue
        name = re.match(r"`([\w.]+)`", cells[0])
        if not name:
            continue
        key = name.group(1)
        path = key if key.startswith(prefix.split(".")[0] + ".") else f"{prefix}.{key}"
        value = lookup(config, path)
        if value is _MISSING or isinstance(value, dict):
            continue
        doc_value = documented_default(cells[default_col])
        if doc_value is None:
            continue
        yaml_value = "" if value == "" else ("[]" if value is None else normalise(value))
        if doc_value != yaml_value:
            out.append(f"{path}: documented {cells[default_col]}, config_default.yaml has {value!r}")
    return out


def main() -> int:
    config = yaml.safe_load(DEFAULT.read_text(encoding="utf-8"))
    doc = DOC.read_text(encoding="utf-8")
    status = 0
    ticked = set(re.findall(r"`([\w.]+)`", doc))
    missing = [k for k in leaf_keys(config)
               if k not in doc and not any(".".join(k.split(".")[i:]) in ticked
                                           for i in range(len(k.split("."))))]
    if missing:
        print(f"{len(missing)} key(s) of config_default.yaml are not documented in doc/doc_parameters.md:")
        for k in missing:
            print(f"  {k}")
        status = 1
    else:
        print(f"all {len(leaf_keys(config))} keys of config_default.yaml are documented")
    wrong = default_mismatches(config, doc)
    if wrong:
        print(f"{len(wrong)} documented default(s) differ from config_default.yaml:")
        for w in wrong:
            print(f"  {w}")
        status = 1
    else:
        print("documented defaults match config_default.yaml")
    return status


if __name__ == "__main__":
    sys.exit(main())
