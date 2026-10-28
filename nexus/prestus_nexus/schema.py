"""Config schema derived from config_default.yaml plus overlay.yaml.

The key tree, key order, default values and help texts all come from
config_default.yaml (help = the inline comment plus its indented
continuation lines). overlay.yaml adds widget kinds, choices, units,
visibility rules and the keys PRESTUS reads that the default omits.
"""

from __future__ import annotations

import copy
import fnmatch
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator

from . import yamlio
from .merge import is_empty, mergestruct
from .paths import package_dir

OVERLAY_PATH = Path(__file__).with_name("overlay.yaml")

KINDS = {"str", "int", "float", "number", "numbers", "vec3", "flag", "bool",
         "tristate", "dir", "file", "path", "strlist", "yaml"}

# Top-level scalar keys are shown in a pseudo-section of this name
GENERAL = "general"


@dataclass
class Node:
    path: str                       # dotted path, e.g. 'transducer.annular.elem_n'
    default: Any = None
    is_group: bool = False
    kind: str = "yaml"
    choices: list | None = None
    unit: str = ""
    help: str = ""
    visible_if: dict | None = None
    derived: bool = False
    extra: bool = False             # read by PRESTUS but absent from config_default
    minimum: float | None = None
    maximum: float | None = None
    label: str = ""                 # display label ('' = derived from the key, see labels.py)
    advanced: bool = False          # hidden unless 'Show advanced options' (or when set)
    tune: list[str] = field(default_factory=list)   # goals this key trades off (overlay 'tuning')
    group: str | None = None        # show in another section's form (top-level sections only)
    group_first: bool = False       # …at the top of that form (default: after its own fields)
    children: list["Node"] = field(default_factory=list)

    @property
    def name(self) -> str:
        return self.path.rsplit(".", 1)[-1]

    @property
    def section(self) -> str:
        return self.path.split(".", 1)[0] if self.is_group or "." in self.path else GENERAL


# ── Help text from comments ──────────────────────────────────────────────

_KEY_LINE = re.compile(r"^(?P<indent> *)(?P<key>[A-Za-z_][\w]*):(?P<rest>.*)$")


def _split_comment(rest: str) -> str:
    """Return the comment part of the text after 'key:' (quote-aware)."""
    quote = None
    for i, ch in enumerate(rest):
        if quote:
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
        elif ch == "#" and (i == 0 or rest[i - 1] in " \t"):
            return rest[i + 1:].strip()
    return ""


_BANNER = re.compile(r"^#\s*-{3,}")
_COMMENTED_KEY = re.compile(r"^#\s*[A-Za-z_]\w*:(\s|$)")


def extract_help(text: str) -> dict[str, str]:
    """Map dotted key paths to their comment text in a YAML config file.

    Help is the inline comment plus following comment-only lines whose '#'
    sits right of the key's own indentation (continuation lines). Keys
    without an inline comment use the comment block directly above them at
    the same indentation, unless that block holds commented-out keys.
    Section banners ('#-----') are ignored.
    """
    helps: dict[str, str] = {}
    stack: list[tuple[int, str]] = []       # (indent, key)
    current: tuple[str, int] | None = None  # (path, key indent) for continuations
    block: list[str] = []                   # comment lines above the next key
    block_col = -1
    block_ok = True
    for line in text.splitlines():
        stripped = line.lstrip(" ")
        if stripped.startswith("#"):
            col = len(line) - len(stripped)
            if current and col > current[1]:
                extra = stripped[1:].strip()
                if extra:
                    helps[current[0]] = (helps.get(current[0], "") + "\n" + extra).strip()
                continue
            current = None
            if _BANNER.match(stripped):
                block, block_ok = [], True
                continue
            if block and col != block_col:
                block, block_ok = [], True
            block_col = col
            block_ok = block_ok and not _COMMENTED_KEY.match(stripped)
            block.append(stripped[1:].strip())
            continue
        if not stripped:
            current = None
            block, block_ok = [], True
            continue
        m = _KEY_LINE.match(line)
        if not m:
            current = None
            block, block_ok = [], True
            continue
        indent = len(m["indent"])
        while stack and stack[-1][0] >= indent:
            stack.pop()
        stack.append((indent, m["key"]))
        path = ".".join(k for _, k in stack)
        comment = _split_comment(m["rest"])
        if comment:
            helps[path] = comment
        elif block and block_ok and block_col == indent:
            helps[path] = "\n".join(b for b in block if b).strip()
        current = (path, indent)
        block, block_ok = [], True
    return helps


# ── Kind inference ───────────────────────────────────────────────────────

def infer_kind(value: Any) -> str:
    if isinstance(value, bool):
        return "bool"
    if isinstance(value, (int, float)):
        return "number"     # MATLAB has no int/float distinction; strict ints come from the overlay
    if isinstance(value, str):
        return "str"
    if isinstance(value, list) and value and all(
            isinstance(v, (int, float)) and not isinstance(v, bool) for v in value):
        return "numbers"
    return "yaml"   # [], None, nested lists: free YAML until the overlay says otherwise


# ── Schema ───────────────────────────────────────────────────────────────

TOOLBOX_DEFAULT = package_dir().parents[1] / "config" / "config_default.yaml"


def _same_file(a: Path, b: Path) -> bool:
    try:
        return a.resolve() == b.resolve() or a.samefile(b)
    except OSError:
        return False


def _leaf_paths(data: Any, prefix: str = "") -> list[str]:
    out: list[str] = []
    if isinstance(data, dict):
        for k, v in data.items():
            path = f"{prefix}.{k}" if prefix else str(k)
            if isinstance(v, dict) and v:
                out.extend(_leaf_paths(v, path))
            else:
                out.append(path)
    return out


def _get(data: Any, path: str) -> Any:
    for part in path.split("."):
        if not isinstance(data, dict) or part not in data:
            return _ABSENT
        data = data[part]
    return data


_ABSENT = object()


def _missing_keys(toolbox: dict, project: dict) -> list[str]:
    """Leaf keys of the toolbox default absent from a project copy."""
    return [p for p in _leaf_paths(toolbox) if _get(project, p) is _ABSENT]


def _differing_values(project: dict, toolbox: dict) -> dict[str, tuple[Any, Any]]:
    """Leaf values a project copy sets differently: {key: (project, toolbox)}."""
    out: dict[str, tuple[Any, Any]] = {}
    for p in _leaf_paths(project):
        tv = _get(toolbox, p)
        pv = _get(project, p)
        if tv is not _ABSENT and pv != tv and not (is_empty(pv) and is_empty(tv)):
            out[p] = (pv, tv)
    return out


def _ordered_merge(toolbox: dict, project: dict) -> dict:
    """mergestruct(toolbox, project), keeping the toolbox key order."""
    return mergestruct(toolbox, project)


class Schema:
    """Key tree for one config_default.yaml."""

    def __init__(self, default_path: str | Path, overlay_path: str | Path = OVERLAY_PATH,
                 toolbox_path: str | Path | None = None):
        """default_path: the config_default.yaml a study is merged onto (a project
        copy or the toolbox file). As in load_parameters, a project copy is merged
        on top of the toolbox default, so keys added to the toolbox since the copy
        was made keep their defaults; stale_keys lists them and pinned the values
        the copy sets differently."""
        self.default_path = Path(default_path)
        self.toolbox_path = Path(toolbox_path) if toolbox_path else TOOLBOX_DEFAULT
        text = self.default_path.read_text(encoding="utf-8")
        own: dict = yamlio.load(text) or {}
        helps = extract_help(text)
        self.project_path: Path | None = None
        self.stale_keys: list[str] = []
        self.pinned: dict[str, tuple[Any, Any]] = {}
        if self.toolbox_path.is_file() and not _same_file(self.toolbox_path, self.default_path):
            tb_text = self.toolbox_path.read_text(encoding="utf-8")
            toolbox = yamlio.load(tb_text) or {}
            self.project_path = self.default_path
            self.stale_keys = _missing_keys(toolbox, own)
            self.pinned = _differing_values(own, toolbox)
            self.defaults: dict = _ordered_merge(toolbox, own)
            helps = {**helps, **extract_help(tb_text)}
        else:
            self.defaults = own
        overlay = yamlio.load_file(overlay_path) or {}
        self.sections: dict[str, dict] = overlay.get("sections", {}) or {}
        self.tabs: list[dict] = overlay.get("tabs", []) or []
        self.field_rules: dict[str, dict] = overlay.get("fields", {}) or {}
        self.display_rules: dict[str, dict] = overlay.get("display", {}) or {}
        self.extra_rules: dict[str, dict] = overlay.get("extra", {}) or {}
        self.deprecated: dict[str, dict] = overlay.get("deprecated", {}) or {}
        self.tune_goals: dict[str, dict] = overlay.get("tune_goals", {}) or {}
        self.tuning_rules: dict[str, list] = overlay.get("tuning", {}) or {}
        self.workflow_keys = self._calibration_keys()

        self.root = Node(path="", is_group=True)
        self._index: dict[str, Node] = {}
        self._build(self.root, self.defaults, helps)
        self._add_extras()
        for node in self._index.values():
            self._apply_rules(node)

    # construction ---------------------------------------------------------
    def _calibration_keys(self) -> set[str]:
        """Dotted keys of config_calibration.yaml (calibration workflow settings).

        The file has no top-level 'calibration:' key; the calibration workflow
        loads it into parameters.calibration (examples/calibration_standalone.m).
        """
        for folder in (self.default_path.parent, package_dir().parents[1] / "config"):
            cal = folder / "config_calibration.yaml"
            if cal.is_file():
                data = yamlio.load_file(cal) or {}
                keys: set[str] = set()

                def walk(d: dict, prefix: str) -> None:
                    for k, v in d.items():
                        path = f"{prefix}.{k}" if prefix else k
                        keys.add(path)
                        if isinstance(v, dict):
                            walk(v, path)
                walk(data, "calibration")
                return keys
        return set()

    def _build(self, parent: Node, data: dict, helps: dict[str, str]) -> None:
        for key, value in data.items():
            path = f"{parent.path}.{key}" if parent.path else key
            node = Node(path=path, default=copy.deepcopy(value),
                        is_group=isinstance(value, dict), help=helps.get(path, ""))
            if not node.is_group:
                node.kind = infer_kind(value)
            parent.children.append(node)
            self._index[path] = node
            if node.is_group:
                self._build(node, value, helps)

    def _add_extras(self) -> None:
        for path, rule in self.extra_rules.items():
            if path in self._index:
                continue    # the drift test reports extras that exist in the default
            parent = self.root
            parts = path.split(".")
            for i, part in enumerate(parts[:-1]):
                sub = ".".join(parts[:i + 1])
                if sub not in self._index:
                    group = Node(path=sub, is_group=True, default={}, extra=True)
                    parent.children.append(group)
                    self._index[sub] = group
                parent = self._index[sub]
            node = Node(path=path, default=copy.deepcopy(rule.get("default")),
                        extra=True, kind=rule.get("kind", "yaml"))
            parent.children.append(node)
            self._index[path] = node

    def _apply_rules(self, node: Node) -> None:
        rules: list[dict] = [r for pat, r in self.field_rules.items()
                             if fnmatch.fnmatchcase(node.path, pat)]
        if node.extra and node.path in self.extra_rules:
            rules.append(self.extra_rules[node.path])
        rules += [r for pat, r in self.display_rules.items() if fnmatch.fnmatchcase(node.path, pat)]
        for rule in rules:
            if "kind" in rule and not node.is_group:
                node.kind = rule["kind"]
            if "choices" in rule:
                node.choices = list(rule["choices"])
            if "unit" in rule:
                node.unit = rule["unit"]
            if "help" in rule:
                node.help = rule["help"]
            if "visible_if" in rule:
                node.visible_if = rule["visible_if"]
            if rule.get("derived"):
                node.derived = True
            if "min" in rule:
                node.minimum = rule["min"]
            if "max" in rule:
                node.maximum = rule["max"]
            if "label" in rule:
                node.label = rule["label"]
            if "advanced" in rule:
                node.advanced = bool(rule["advanced"])
            if "group" in rule:
                node.group = rule["group"]
            if "group_first" in rule:
                node.group_first = bool(rule["group_first"])
        for pat, goals in self.tuning_rules.items():
            if fnmatch.fnmatchcase(node.path, pat):
                node.tune = list(goals or [])
        if node.choices and node.kind == "yaml":
            node.kind = "str"

    # queries --------------------------------------------------------------
    def get(self, path: str) -> Node | None:
        return self._index.get(path)

    def __contains__(self, path: str) -> bool:
        return path in self._index

    def leaves(self) -> Iterator[Node]:
        for node in self._index.values():
            if not node.is_group:
                yield node

    def nodes(self) -> Iterator[Node]:
        yield from self._index.values()

    def section_order(self) -> list[str]:
        """Section keys: 'general' first, then top-level groups in default order."""
        order = [GENERAL]
        order += [n.name for n in self.root.children if n.is_group]
        return order

    def section_nodes(self, section: str) -> list[Node]:
        if section == GENERAL:
            return [n for n in self.root.children if not n.is_group]
        node = self._index.get(section)
        return node.children if node else []

    def display_nodes(self, section: str) -> list[Node]:
        """Nodes shown in a section's form: its own (minus those moved elsewhere)
        followed by nodes other sections moved here with ``group:``."""
        own = [n for n in self.section_nodes(section) if not n.group or n.group == section]
        moved = [n for n in self._index.values() if n.group == section and n.section != section]
        return [n for n in moved if n.group_first] + own + [n for n in moved if not n.group_first]

    def tab_of(self, section: str) -> str | None:
        for tab in self.tabs:
            if section in tab.get("sections", []):
                return tab["key"]
        return None

    def section_label(self, section: str) -> str:
        return self.sections.get(section, {}).get("label", section.replace("_", " ").capitalize())

    def section_is_advanced(self, section: str) -> bool:
        return bool(self.sections.get(section, {}).get("advanced", False))

    def section_is_locked(self, section: str) -> bool:
        return bool(self.sections.get(section, {}).get("locked", False))

    def tune_label(self, goal: str) -> str:
        return self.tune_goals.get(goal, {}).get("label", goal.capitalize())

    def ancestors(self, path: str) -> list[Node]:
        """Nodes from the top-level group down to (and including) path."""
        parts = path.split(".")
        return [self._index[".".join(parts[:i + 1])] for i in range(len(parts))
                if ".".join(parts[:i + 1]) in self._index]

    def suggest(self, path: str) -> str | None:
        """Closest known key for a probably misspelt path (same parent)."""
        import difflib
        parent, _, name = path.rpartition(".")
        siblings = [n.name for n in (self._index[parent].children if parent in self._index
                                     else self.root.children if not parent else [])]
        match = difflib.get_close_matches(name, siblings, n=1, cutoff=0.6)
        if not match:
            return None
        return f"{parent}.{match[0]}" if parent else match[0]
