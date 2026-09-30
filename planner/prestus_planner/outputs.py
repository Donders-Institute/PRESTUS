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
