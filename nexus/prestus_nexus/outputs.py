"""Resolves the on-disk locations the Results tab and pseudo-CT viewer read from,
the same way io_preview.py previews them (path.sim / path.pct + sub-XXX) and
path_log_setup.m lays them out (nii/, cache/ under dir_output). Read-only: this
never creates directories, only says where the pipeline would have written them.
"""

from __future__ import annotations

from pathlib import Path

from .model import ConfigModel, is_empty
from .paths import m2m_dir


def _subject_id(model: ConfigModel) -> int | None:
    sid = model.get("subject_id")
    if isinstance(sid, (int, float)) and not isinstance(sid, bool) and sid == int(sid) and 0 <= sid <= 999:
        return int(sid)
    return None


def dir_output(model: ConfigModel) -> Path | None:
    """<path.sim>/sub-XXX, where the CSV/report/nii/cache/debug for this subject land."""
    sim, sid = model.get("path.sim"), _subject_id(model)
    if is_empty(sim) or sid is None:
        return None
    return Path(sim) / f"sub-{sid:03d}"


def dir_pct(model: ConfigModel) -> Path | None:
    """Where pseudoCT.nii.gz lives: path.pct/sub-XXX if set, else the SimNIBS m2m folder."""
    sid = model.get("subject_id")
    pct, seg_root = model.get("path.pct"), model.get("path.seg") or model.get("path.anat")
    sid = _subject_id(model)
    if sid is None:
        return None
    if is_empty(pct):
        return None if is_empty(seg_root) else m2m_dir(seg_root, sid)
    return Path(pct) / f"sub-{sid:03d}"


def placement_record(model: ConfigModel) -> dict | None:
    """Positions the placement stage resolved, from <dir_output>/sub-XXX_placement<affix>.json
    (write_placement_positions.m): {'mode', 'created', 'transducers': [{'trans_pos', 'focus_pos'}]}
    with 1-based T1 voxel positions; None when the file is missing or unreadable."""
    import json

    from .study import sanitize_affix
    base, sid = dir_output(model), _subject_id(model)
    if base is None or sid is None:
        return None
    affix = sanitize_affix(model.get("io.output_affix") or "")
    path = base / f"sub-{sid:03d}_placement{affix}.json"
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if not isinstance(data, dict) or not isinstance(data.get("transducers"), list):
        return None
    data["path"] = str(path)
    return data