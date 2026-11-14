"""Refreshing a project's copy of config_default.yaml (prestus_config_upgrade.m).

A project copy pins the defaults a study is merged onto. When the toolbox
gains keys, the copy lacks them (Schema.stale_keys); load_parameters then
fills them from the toolbox default and warns. Replacing the copy with the
toolbox default ends the warning but also drops the values the copy pinned
differently (Schema.pinned); carry_pinned_values moves those into the study
config first, so the effective configuration stays the same.
"""

from __future__ import annotations

import copy
import time
from pathlib import Path
from typing import TYPE_CHECKING

from .fsio import copy_atomic
from .model import TX, get_in, set_in

if TYPE_CHECKING:
    from .model import ConfigModel

_ABSENT = object()


def carry_pinned_values(model: "ConfigModel", pinned: dict) -> list[str]:
    """Write the project copy's pinned values into the study's overrides.

    Keys the study already sets are left alone. Transducer defaults go to
    transducer 1, which later transducers inherit from. Returns the keys
    written. The values are written directly (not via ConfigModel.set, which
    would drop them as equal to the current, still pinned, baseline).
    """
    with model.edit("Carry pinned defaults"):        # one undo step
        return _carry(model, pinned)


def _carry(model: "ConfigModel", pinned: dict) -> list[str]:
    written: list[str] = []
    for key, (project_value, _toolbox_value) in pinned.items():
        parts = key.split(".")
        if parts[0] == TX:
            if not model.tx:
                continue
            target, key_parts = model.tx[0], parts[1:]
        else:
            target, key_parts = model.overrides, parts
        if not key_parts or get_in(target, key_parts, _ABSENT) is not _ABSENT:
            continue
        set_in(target, key_parts, copy.deepcopy(project_value))
        written.append(key)
    return written


def replace_project_default(project: Path, toolbox: Path) -> Path:
    """Replace the project copy by the toolbox default; return the backup's path."""
    project, toolbox = Path(project), Path(toolbox)
    backup = project.with_name(f"{project.name}.bak-{time.strftime('%Y%m%d-%H%M%S')}")
    copy_atomic(project, backup)
    copy_atomic(toolbox, project)            # the project copy stays intact if the copy fails
    return backup
