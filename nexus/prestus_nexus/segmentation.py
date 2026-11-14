"""The SimNIBS segmentation as a step of its own: where it is, whether it is complete, and
which settings use it. Many stages read it (the layered head model; heuristic, MNI, PlanTUS
and TUS_entry placement; the Viewer's scalp and beam checks; the 3D head), and SimNIBS charm
takes an hour or more, so the Nexus offers to run it on its own (modules.segmentation_only,
launcher.segmentation_overrides) before anything else.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .model import ConfigModel, is_empty
from .paths import m2m_dir

# what charm writes that PRESTUS reads: tissue labels and the MNI warps (checks._segmentation_status)
REQUIRED = ("final_tissues.nii.gz", "toMNI/MNI2Conform_nonl.nii.gz")
MODES = ("heuristic", "mni", "plantus", "tus_entry")      # placement modes that read it


@dataclass(frozen=True)
class SegStatus:
    state: str              # 'found' | 'incomplete' | 'missing' | 'unset'
    m2m: Path | None
    message: str
    can_run: bool           # startup.simnibs_bin_path is set (where charm runs may be the cluster)


def _subject(model: ConfigModel) -> int | None:
    sid = model.get("subject_id")
    return int(sid) if isinstance(sid, (int, float)) and not isinstance(sid, bool) else None


def m2m_path(model: ConfigModel) -> Path | None:
    """<path.seg or path.anat>/m2m_sub-XXX, without touching the disk."""
    sid = _subject(model)
    root = model.get("path.seg")
    root = model.get("path.anat") if is_empty(root) else root
    if sid is None or is_empty(root):
        return None
    return m2m_dir(str(root), sid)


def needed_by(model: ConfigModel) -> list[str]:
    """What in this config uses the segmentation (empty: nothing needs it)."""
    out = []
    if str(model.get("simulation.medium") or "layered") == "layered":
        out.append("the layered head model")
    mode = str(model.get("placement.mode") or "manual")
    if mode in MODES:
        out.append(f"{mode} placement")
    return out


def status(model: ConfigModel) -> SegStatus:
    simnibs = model.get("startup.simnibs_bin_path")
    can_run = not is_empty(simnibs)
    m2m = m2m_path(model)
    if m2m is None:
        return SegStatus("unset", None, "Set subject_id and path.seg (or path.anat) to find the segmentation.",
                         can_run)
    if not m2m.is_dir():
        return SegStatus("missing", m2m, f"No segmentation yet ({m2m.name}).", can_run)
    missing = [r for r in REQUIRED if not (m2m / r).is_file()]
    if missing:
        return SegStatus("incomplete", m2m, f"{m2m.name} lacks {', '.join(missing)} (charm unfinished "
                         "or failed?).", can_run)
    return SegStatus("found", m2m, f"Found: {m2m}", can_run)
