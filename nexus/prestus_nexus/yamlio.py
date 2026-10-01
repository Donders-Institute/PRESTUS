"""YAML reading and writing compatible with PRESTUS' MATLAB loader.

PRESTUS parses configs with ``yaml.loadFile`` (SnakeYAML 1.30, YAML 1.1).
Two differences from PyYAML matter:

* SnakeYAML accepts floats such as ``2.17e-3`` (no sign in the exponent,
  used in config_default.yaml); PyYAML reads those as strings. ``load``
  uses a loader with SnakeYAML's float rule.
* YAML 1.1 turns unquoted ``yes``/``no``/``on``/``off`` into booleans and
  ``04:00:00`` into a sexagesimal number. ``dump`` therefore single-quotes
  every string and writes floats so that both parsers read them back as
  floats.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

import yaml


class _Loader(yaml.SafeLoader):
    """SafeLoader with SnakeYAML's (lenient) float resolver."""


# Remove PyYAML's float resolver and install SnakeYAML's, verbatim from
# org.yaml.snakeyaml.resolver.Resolver.FLOAT (SnakeYAML 1.30).
_Loader.yaml_implicit_resolvers = {
    ch: [(tag, rx) for tag, rx in resolvers if tag != "tag:yaml.org,2002:float"]
    for ch, resolvers in yaml.SafeLoader.yaml_implicit_resolvers.items()
}
_Loader.add_implicit_resolver(
    "tag:yaml.org,2002:float",
    re.compile(
        r"^([-+]?(?:[0-9][0-9_]*)\.[0-9_]*(?:[eE][-+]?[0-9]+)?"
        r"|[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)"
        r"|[-+]?\.[0-9_]+(?:[eE][-+]?[0-9]+)?"
        r"|[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\.[0-9_]*"
        r"|[-+]?\.(?:inf|Inf|INF)"
        r"|\.(?:nan|NaN|NAN))$"
    ),
    list("-+0123456789."),
)


def load(text: str) -> Any:
    """Parse YAML text the way PRESTUS' MATLAB loader would."""
    return yaml.load(text, Loader=_Loader)


def load_file(path: str | Path) -> Any:
    with open(path, encoding="utf-8") as fh:
        return load(fh.read())


# ── Writing ──────────────────────────────────────────────────────────────

def format_scalar(value: Any) -> str:
    """Format a scalar so SnakeYAML and PyYAML both read it back unchanged."""
    if value is None:
        return "~"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if math.isnan(value):
            return ".nan"
        if math.isinf(value):
            return ".inf" if value > 0 else "-.inf"
        text = repr(value)                      # shortest round-trip repr
        mantissa, _, exponent = text.partition("e")
        if "." not in mantissa:
            mantissa += ".0"
        if exponent and exponent[0] not in "+-":
            exponent = "+" + exponent
        return mantissa + ("e" + exponent if exponent else "")
    if isinstance(value, str):
        return "'" + value.replace("'", "''") + "'"
    raise TypeError(f"Unsupported YAML scalar type: {type(value).__name__}")


def _is_scalar(value: Any) -> bool:
    return value is None or isinstance(value, (bool, int, float, str))


def _flow(value: Any) -> str:
    """Flow-style text for a scalar or a (nested) list of scalars."""
    if isinstance(value, list):
        return "[" + ", ".join(_flow(v) for v in value) + "]"
    return format_scalar(value)


def _is_flow_list(value: Any) -> bool:
    return isinstance(value, list) and all(
        _is_scalar(v) or _is_flow_list(v) for v in value)


def _emit(value: Any, indent: int, lines: list[str]) -> None:
    pad = " " * indent
    for key, val in value.items():
        if isinstance(val, dict):
            if val:
                lines.append(f"{pad}{key}:")
                _emit(val, indent + 2, lines)
            else:
                lines.append(f"{pad}{key}: {{}}")
        elif _is_flow_list(val):
            lines.append(f"{pad}{key}: {_flow(val)}")
        elif isinstance(val, list):
            # List of mappings (e.g. multiple transducers): block style
            lines.append(f"{pad}{key}:")
            for item in val:
                if not isinstance(item, dict):
                    raise TypeError(f"Mixed list under '{key}' is not supported")
                sub: list[str] = []
                _emit(item, indent + 4, sub)
                if not sub:
                    lines.append(f"{pad}  - {{}}")
                    continue
                lines.append(f"{pad}  - {sub[0].lstrip()}")
                lines.extend(sub[1:])
        else:
            lines.append(f"{pad}{key}: {format_scalar(val)}")


def dump(data: dict, header: str = "") -> str:
    """Serialise a nested dict as MATLAB-compatible block YAML."""
    lines: list[str] = []
    if header:
        lines.extend("# " + h if h else "#" for h in header.splitlines())
        lines.append("")
    _emit(data, 0, lines)
    return "\n".join(lines) + "\n"


def dump_file(path: str | Path, data: dict, header: str = "") -> None:
    Path(path).write_text(dump(data, header), encoding="utf-8")
