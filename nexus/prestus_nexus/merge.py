"""PRESTUS merge semantics (functions/helper/mergestruct.m), shared by schema and model."""

from __future__ import annotations

import copy
from typing import Any


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
