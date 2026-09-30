"""Study configuration: defaults + sparse overrides, merged like PRESTUS.

The effective configuration is computed exactly as load_parameters does:
``mergestruct(defaults, overrides)`` followed by ``restore_struct_arrays``.
The latter merges transducer 2..N onto transducer 1's *merged* values, so
later transducers inherit transducer 1's overrides, not just the defaults.

Only overrides are written to the study YAML.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterable

from . import yamlio
from .schema import Node, Schema

if TYPE_CHECKING:
    from .equipment import Device, EquipmentLibrary

TX = "transducer"


# ── Merge semantics (functions/helper/mergestruct.m) ─────────────────────

def is_empty(value: Any) -> bool:
    """MATLAB isempty() for values read from YAML ('' and [] and ~ are empty)."""
    return value is None or value == [] or value == "" or value == {}


def mergestruct(a: Any, b: Any) -> Any:
    """Python port of PRESTUS mergestruct: b's non-struct values replace a's."""
    if not isinstance(b, dict):
        return copy.deepcopy(a)
    res = copy.deepcopy(a) if isinstance(a, dict) else {}
    for key, new in b.items():
        old = res.get(key) if isinstance(a, dict) else None
        if is_empty(old) or not isinstance(new, dict):
            res[key] = copy.deepcopy(new)
        else:
            res[key] = mergestruct(old, new)
    return res


def _merge_unset(base: dict, override: dict) -> dict:
    """mergestruct_shallow in resolve_transducer_from_serial.m: override wins,
    except that empty override values do not replace non-empty base values."""
    out = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _merge_unset(out[key], value)
        elif is_empty(value) and not is_empty(out.get(key)):
            continue
        else:
            out[key] = copy.deepcopy(value)
    return out


def values_equal(a: Any, b: Any) -> bool:
    """Equality as MATLAB sees values after loading.

    Numbers compare by value (1 == 1.0). Booleans only equal booleans and
    '' / [] / ~ only equal their own kind, because yaml.loadFile maps them to
    different MATLAB types (logical vs double, string vs [] vs yaml.Null).
    """
    if isinstance(a, bool) or isinstance(b, bool):
        return type(a) is type(b) and a == b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return float(a) == float(b)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(values_equal(x, y) for x, y in zip(a, b))
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(values_equal(a[k], b[k]) for k in a)
    return type(a) is type(b) and a == b


# ── Nested dict helpers ──────────────────────────────────────────────────

_MISSING = object()     # get_in default: raise KeyError
_ABSENT = object()      # get_in default: key not present


def get_in(data: Any, parts: list[str], default: Any = _MISSING) -> Any:
    for part in parts:
        if not isinstance(data, dict) or part not in data:
            if default is _MISSING:
                raise KeyError(".".join(parts))
            return default
        data = data[part]
    return data


def set_in(data: dict, parts: list[str], value: Any) -> None:
    for part in parts[:-1]:
        if not isinstance(data.get(part), dict):
            data[part] = {}
        data = data[part]
    data[parts[-1]] = value


def del_in(data: dict, parts: list[str]) -> None:
    """Delete a nested key and prune dicts left empty."""
    trail = []
    node = data
    for part in parts[:-1]:
        if not isinstance(node.get(part), dict):
            return
        trail.append((node, part))
        node = node[part]
    node.pop(parts[-1], None)
    for parent, key in reversed(trail):
        if parent[key] == {}:
            del parent[key]
        else:
            break


def iter_leaves(data: dict, prefix: str = "") -> Iterable[tuple[str, Any]]:
    for key, value in data.items():
        path = f"{prefix}.{key}" if prefix else key
        if isinstance(value, dict) and value:
            yield from iter_leaves(value, path)
        else:
            yield path, value


@dataclass
class Issue:
    level: str          # 'error' | 'warning' | 'info'
    path: str           # dotted key ('' for general)
    message: str
    transducer: int | None = None   # 0-based transducer index, if any

    def where(self) -> str:
        if self.transducer is not None and self.path.startswith(TX + "."):
            return f"transducer[{self.transducer + 1}].{self.path[len(TX) + 1:]}"
        return self.path


class ConfigModel:
    """Defaults + sparse overrides for one study config."""

    def __init__(self, schema: Schema, library: "EquipmentLibrary | None" = None):
        self.schema = schema
        self.library = library            # equipment library for transducer.serial
        self.overrides: dict = {}         # everything except the transducer block
        self.tx: list[dict] = [{}]        # sparse override per transducer
        self.tx_as_list = False           # write transducer as a YAML list even if N == 1
        self.path: Path | None = None     # study YAML file

    # ── loading / saving ────────────────────────────────────────────────
    @classmethod
    def from_file(cls, schema: Schema, path: str | Path,
                  library: "EquipmentLibrary | None" = None) -> "ConfigModel":
        model = cls(schema, library)
        data = yamlio.load_file(path) or {}
        if not isinstance(data, dict):
            raise ValueError(f"{path}: top level must be a mapping")
        model.load_dict(data)
        model.path = Path(path)
        return model

    def load_dict(self, data: dict) -> None:
        data = copy.deepcopy(data)
        tx = data.pop(TX, None)
        self.overrides = data
        if isinstance(tx, list):
            self.tx = [t if isinstance(t, dict) else {} for t in tx] or [{}]
            self.tx_as_list = True
        elif isinstance(tx, dict):
            self.tx = [tx]
            self.tx_as_list = False
        else:
            self.tx = [{}]
            self.tx_as_list = False

    def override_dict(self) -> dict:
        """Sparse overrides in config_default order (unknown keys last)."""
        combined = copy.deepcopy(self.overrides)
        tx = self._aligned_transducers()
        if any(tx) or len(tx) > 1:
            combined[TX] = tx if (len(tx) > 1 or self.tx_as_list) else tx[0]
        return self._ordered(combined, self.schema.root)

    def _aligned_transducers(self) -> list[dict]:
        """Transducer overrides with the same top-level keys in every element.

        yaml.loadFile (ConvertToArray) concatenates a list of mappings into a
        struct array, which fails unless all elements have identical fields.
        Keys transducer 1 sets are therefore repeated on later elements: plain
        values with the value they inherit anyway, groups as {} (mergestruct
        then keeps the inherited sub-fields). The meaning is unchanged.
        """
        tx = copy.deepcopy(self.tx)
        if len(tx) < 2:
            return tx
        first = mergestruct(self.schema.defaults.get(TX, {}), self.tx[0])
        for t in tx[1:]:
            for key, value in tx[0].items():
                if key not in t:
                    t[key] = {} if isinstance(value, dict) else copy.deepcopy(first[key])
        return tx

    def _ordered(self, data: dict, node: Node) -> dict:
        out: dict = {}
        known = {c.name: c for c in node.children}
        for child in node.children:
            if child.name not in data:
                continue
            value = data[child.name]
            if child.is_group and isinstance(value, dict):
                out[child.name] = self._ordered(value, child)
            elif child.is_group and isinstance(value, list):
                out[child.name] = [self._ordered(v, child) if isinstance(v, dict) else v
                                   for v in value]
            else:
                out[child.name] = value
        for key, value in data.items():
            if key not in known:
                out[key] = value
        return out

    def to_yaml(self, header: str = "") -> str:
        return yamlio.dump(self.override_dict(), header)

    def save(self, path: str | Path | None = None, header: str = "") -> Path:
        target = Path(path) if path else self.path
        if target is None:
            raise ValueError("No file name given")
        yamlio.dump_file(target, self.override_dict(), header)
        self.path = target
        return target

    # ── effective values ────────────────────────────────────────────────
    def effective(self) -> dict:
        """The merged configuration, as load_parameters would produce it.

        Transducers with a serial also get the equipment-library geometry
        (resolve_transducer_from_serial); calibrated phases and amplitudes
        are computed by MATLAB and not modelled here.
        """
        merged = mergestruct(self.schema.defaults, self.overrides)
        base_tx = merged.get(TX, {})
        first = mergestruct(base_tx, self.tx[0])
        txs = [first] + [mergestruct(first, t) for t in self.tx[1:]]
        txs = [self._resolve_serial(t) for t in txs]
        merged[TX] = txs if len(txs) > 1 else txs[0]
        return merged

    def effective_transducers(self) -> list[dict]:
        tx = self.effective()[TX]
        return tx if isinstance(tx, list) else [tx]

    # ── equipment library layer ─────────────────────────────────────────
    def device(self, ti: int = 0) -> "Device | None":
        """Library device named by transducer ti's serial, if any."""
        if self.library is None:
            return None
        return self.library.get(self._tx_value("transducer.serial", ti))

    def _resolve_serial(self, tx: dict) -> dict:
        dev = self.library.get(tx.get("serial")) if self.library else None
        if dev is None:
            return tx
        out = _merge_unset(dev.geometry, tx)
        if is_empty(out.get("type")) and dev.type:
            out["type"] = dev.type
        return out

    def _library_value(self, path: str, ti: int) -> Any:
        dev = self.device(ti)
        if dev is None:
            return _ABSENT
        parts = path.split(".")[1:]
        if parts == ["type"]:
            return dev.type or _ABSENT
        value = get_in(dev.geometry, parts, _ABSENT)
        return _ABSENT if value is _ABSENT or is_empty(value) else copy.deepcopy(value)

    # ── values ──────────────────────────────────────────────────────────
    def _tx_value(self, path: str, ti: int) -> Any:
        """Transducer value from overrides and inheritance only (no library)."""
        parts = path.split(".")
        if get_in(self.tx[ti], parts[1:], _ABSENT) is not _ABSENT:
            return copy.deepcopy(get_in(self.tx[ti], parts[1:]))
        return self._config_baseline(path, ti)

    def _config_baseline(self, path: str, ti: int) -> Any:
        parts = path.split(".")
        if parts[0] == TX:
            if ti == 0:
                base = self.schema.defaults.get(TX, {})
                return get_in(base, parts[1:], self._node_default(path))
            first = mergestruct(self.schema.defaults.get(TX, {}), self.tx[0])
            return get_in(first, parts[1:], self._node_default(path))
        return get_in(self.schema.defaults, parts, self._node_default(path))

    def baseline(self, path: str, ti: int = 0) -> Any:
        """Value a key has when not overridden: default, transducer 1, or library."""
        value = self._config_baseline(path, ti)
        if path.startswith(TX + ".") and is_empty(value):
            lib = self._library_value(path, ti)
            if lib is not _ABSENT:
                return lib
        return value

    def _node_default(self, path: str) -> Any:
        node = self.schema.get(path)
        return copy.deepcopy(node.default) if node else None

    def get(self, path: str, ti: int = 0) -> Any:
        parts = path.split(".")
        if parts[0] == TX:
            if self.is_overridden(path, ti):
                value = copy.deepcopy(get_in(self.tx[ti], parts[1:]))
                if is_empty(value):             # empty study values count as unset
                    lib = self._library_value(path, ti)
                    return value if lib is _ABSENT else lib
                return value
            return self.baseline(path, ti)
        if self.is_overridden(path):
            return copy.deepcopy(get_in(self.overrides, parts))
        return self.baseline(path, ti)

    def from_library(self, path: str, ti: int = 0) -> bool:
        """The shown value comes from the equipment library."""
        return path.startswith(TX + ".") and not self.is_modified(path, ti) and \
            self._library_value(path, ti) is not _ABSENT and \
            is_empty(self._config_baseline(path, ti))

    def is_overridden(self, path: str, ti: int = 0) -> bool:
        parts = path.split(".")
        if parts[0] == TX:
            return get_in(self.tx[ti], parts[1:], _ABSENT) is not _ABSENT
        return get_in(self.overrides, parts, _ABSENT) is not _ABSENT

    def is_modified(self, path: str, ti: int = 0) -> bool:
        """Overridden with a value that differs from its baseline."""
        return self.is_overridden(path, ti) and not values_equal(
            self.get(path, ti), self.baseline(path, ti))

    def set(self, path: str, value: Any, ti: int = 0) -> None:
        """Set a value; setting it back to the baseline removes the override."""
        parts = path.split(".")
        target, key_parts = (self.tx[ti], parts[1:]) if parts[0] == TX else (self.overrides, parts)
        if values_equal(value, self.baseline(path, ti)):
            del_in(target, key_parts)
        else:
            set_in(target, key_parts, copy.deepcopy(value))

    def reset(self, path: str, ti: int = 0) -> None:
        parts = path.split(".")
        if parts[0] == TX:
            del_in(self.tx[ti], parts[1:])
        else:
            del_in(self.overrides, parts)

    def prune_defaults(self) -> int:
        """Remove overrides whose value equals their baseline; return how many."""
        removed = 0
        for path, _ in list(iter_leaves(self.overrides)):
            if path in self.schema and not self.is_modified(path):
                self.reset(path)
                removed += 1
        for ti in range(len(self.tx)):
            for sub, _ in list(iter_leaves(self.tx[ti])):
                path = f"{TX}.{sub}"
                if path in self.schema and not self.is_modified(path, ti):
                    self.reset(path, ti)
                    removed += 1
        return removed

    def modified_count(self) -> int:
        n = sum(1 for _ in iter_leaves(self.overrides))
        n += sum(1 for t in self.tx for _ in iter_leaves(t))
        return n

    # ── transducers ─────────────────────────────────────────────────────
    @property
    def n_transducers(self) -> int:
        return len(self.tx)

    def add_transducer(self) -> int:
        """Append a transducer. It starts from transducer 1's values (PRESTUS semantics)."""
        self.tx.append({})
        return len(self.tx) - 1

    def remove_transducer(self, ti: int) -> None:
        if len(self.tx) == 1:
            self.tx = [{}]
            return
        del self.tx[ti]

    # ── visibility ──────────────────────────────────────────────────────
    def condition_holds(self, cond: dict | None, ti: int = 0) -> bool:
        if not cond:
            return True
        value = self.get(cond["key"], ti)
        if "in" in cond:
            return any(values_equal(value, c) for c in cond["in"])
        if "not_in" in cond:
            return not any(values_equal(value, c) for c in cond["not_in"])
        return True

    def is_visible(self, path: str, ti: int = 0) -> bool:
        return all(self.condition_holds(n.visible_if, ti) for n in self.schema.ancestors(path))

    # ── unknown keys ────────────────────────────────────────────────────
    def unknown_keys(self) -> list[Issue]:
        issues: list[Issue] = []

        def walk(data: dict, prefix: str, ti: int | None) -> None:
            for key, value in data.items():
                path = f"{prefix}.{key}" if prefix else key
                node = self.schema.get(path)
                if path in self.schema.deprecated:
                    dep = self.schema.deprecated[path]
                    note = "no longer read by PRESTUS" if dep.get("ignored") else "PRESTUS still converts it"
                    issues.append(Issue("warning", path, f"Deprecated; use {dep['use']} ({note})", ti))
                elif node is None and path in self.schema.workflow_keys:
                    issues.append(Issue("info", path,
                                        "Calibration-workflow setting (config_calibration.yaml); "
                                        "not used by the simulation pipeline", ti))
                    continue
                elif node is None:
                    hint = self.schema.suggest(path)
                    msg = "Unknown key; PRESTUS will ignore it"
                    msg += f" (did you mean '{hint.rsplit('.', 1)[-1]}'?)" if hint else ""
                    issues.append(Issue("warning", path, msg, ti))
                elif node.is_group and isinstance(value, dict):
                    walk(value, path, ti)

        walk(self.overrides, "", None)
        for ti, t in enumerate(self.tx):
            walk(t, TX, ti)
        return issues
