"""Resolves the on-disk locations the Results tab and pseudo-CT viewer read from,
the same way io_preview.py previews them (path.sim / path.pct + sub-XXX) and
path_log_setup.m lays them out (nii/, cache/ under dir_output). Read-only: this
never creates directories, only says where the pipeline would have written them.
"""

from __future__ import annotations

from pathlib import Path

from .model import ConfigModel, is_empty
from .paths import m2m_dir


_TEXT_CACHE: dict[str, tuple[tuple[int, int], str]] = {}


def _read_text(path: Path) -> str:
    """The file's text, reused while its size and modification time are unchanged: a stat
    on the project volume is quick, but each open costs a round trip (about 0.2 s), and
    the placement views read these files on every redraw. Raises OSError like read_text."""
    st = path.stat()
    stamp = (st.st_size, st.st_mtime_ns)
    hit = _TEXT_CACHE.get(str(path))
    if hit is not None and hit[0] == stamp:
        return hit[1]
    text = path.read_bytes().decode("utf-8")
    if len(_TEXT_CACHE) > 64:
        _TEXT_CACHE.clear()
    _TEXT_CACHE[str(path)] = (stamp, text)
    return text


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


def dir_placement_run(model: ConfigModel, mode: str) -> Path | None:
    """Folder of one placement run (placement_run_dir.m): placement.<mode>.work_dir if set,
    else <path.placement or path.sim/placement>/sub-XXX/<mode>/<target_name>."""
    work_dir = model.get(f"placement.{mode}.work_dir")
    if not is_empty(work_dir):
        return Path(work_dir)
    sid = _subject_id(model)
    root = model.get("path.placement")
    if is_empty(root):
        sim = model.get("path.sim")
        root = None if is_empty(sim) else Path(sim) / "placement"
    if root is None or sid is None:
        return None
    target = model.get(f"placement.{mode}.target_name")
    return Path(root) / f"sub-{sid:03d}" / mode / ("target" if is_empty(target) else str(target))


def dir_pct(model: ConfigModel) -> Path | None:
    """Where the pseudo-CT is written (get_pct_dir.m): path.pct/sub-XXX if set, else
    <parent of path.seg>/pct/sub-XXX, next to (never inside) the SimNIBS folder."""
    pct, seg_root = model.get("path.pct"), model.get("path.seg")
    sid = _subject_id(model)
    if sid is None:
        return None
    if not is_empty(pct):
        return Path(pct) / f"sub-{sid:03d}"
    if is_empty(seg_root):
        return None
    return Path(str(seg_root).rstrip("/\\")).parent / "pct" / f"sub-{sid:03d}"


def pseudoct_file(model: ConfigModel) -> tuple[Path | None, bool]:
    """(existing pseudoCT.nii.gz, found in the legacy SimNIBS m2m location?), as
    find_pseudoct_dir.m; (None, False) when there is none yet."""
    pct_dir = dir_pct(model)
    if pct_dir is not None and (pct_dir / "pseudoCT.nii.gz").is_file():
        return pct_dir / "pseudoCT.nii.gz", False
    seg_root, sid = model.get("path.seg"), _subject_id(model)
    if not is_empty(seg_root) and sid is not None:
        legacy = m2m_dir(seg_root, sid) / "pseudoCT.nii.gz"
        if legacy.is_file():
            return legacy, True
    return None, False


def placement_record(model: ConfigModel) -> dict | None:
    """Positions the placement stage resolved, from <dir_output>/sub-XXX_placement<affix>.json
    (write_placement_positions.m): {'mode', 'created', 'transducers': [{'trans_pos', 'focus_pos'}]}
    with 1-based T1 voxel positions; None when the file is missing or unreadable."""
    import json
def measured_ct_file(model: ConfigModel) -> Path | None:
    """The measured CT path.ct_pattern names (as pct_ct_file.m), or None when it is empty or
    matches nothing."""
    from .paths import resolve_pattern
    pattern, anat, sid = model.get("path.ct_pattern"), model.get("path.anat"), _subject_id(model)
    if is_empty(pattern) or is_empty(anat) or sid is None:
        return None
    hits = resolve_pattern(anat, pattern, sid)
    return hits[0] if hits else None


_EXISTING: dict = {}


def existing_outputs(model: ConfigModel) -> dict[str, bool] | None:
    """Which optional output categories (predict_outputs.save_items) already have files in this
    subject's output folder, for this medium and output affix: {label: exists}. None when
    path.sim or subject_id is unset. Cached until a folder it reads changes (mtime)."""
    from .predict_outputs import save_items
    out = dir_output(model)
    if out is None:
        return None
    items = save_items(model)
    folders = (out, out / "cache", out / "nii", out / "nii" / "properties", out / "debug")
    stamp = tuple(p.stat().st_mtime_ns if p.is_dir() else None for p in folders)
    key = (str(out), tuple(i.pattern for i in items), stamp)
    if key not in _EXISTING:
        _EXISTING.clear()
        _EXISTING[key] = {i.label: bool(i.pattern) and out.is_dir() and any(out.glob(i.pattern))
                          for i in items}
    return dict(_EXISTING[key])


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