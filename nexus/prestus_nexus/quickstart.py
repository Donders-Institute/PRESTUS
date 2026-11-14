"""QuickStart tab logic: everything a first run needs, derived from one T1 image.

PRESTUS finds its inputs through path.anat + path.t1_pattern (sprintf'd with subject_id),
looks for the segmentation in <path.seg or path.anat>/m2m_sub-XXX, and needs path.sim for
its outputs. A user who picks a T1 file should not have to know any of that: from_t1()
proposes all of these keys, and the QuickStart tab writes them into the model.
"""

from __future__ import annotations

import re
import shutil
from dataclasses import dataclass, field
from pathlib import Path

from . import simnibs
from .paths import m2m_dir, matlab_sprintf, resolve_pattern

NIFTI_FILTER = "NIfTI images (*.nii *.nii.gz)"
_SUB = re.compile(r"sub-(\d+)")
_T1_TAGS = (("T1w", "T2w"), ("T1", "T2"), ("t1", "t2"))


@dataclass
class T1Inputs:
    anat: str
    t1_pattern: str
    subject_id: int | None          # None: the file name has no sub-<n>; keep the config's ID
    t2_pattern: str = ""            # '' = no T2 found next to the T1
    seg: str = ""                   # '' = the segmentation is (or will be) in path.anat
    m2m_found: bool = False
    notes: list[str] = field(default_factory=list)

    def as_settings(self) -> dict:
        """Config keys and values to set (subject_id only when the file name gives one)."""
        out = {"path.anat": self.anat, "path.t1_pattern": self.t1_pattern,
               "path.t2_pattern": self.t2_pattern, "path.seg": self.seg}
        if self.subject_id is not None:
            out["subject_id"] = self.subject_id
        return out


def _pattern_for(name: str, sid: int | None) -> str:
    """A t1_pattern that gives back `name` for subject `sid`: 'sub-012_T1w.nii.gz' becomes
    'sub-%1$03d_T1w.nii.gz', so a study manifest can reuse it for other subjects. Without a
    subject number, the file name itself (with % escaped for sprintf)."""
    literal = name.replace("%", "%%")
    m = _SUB.search(name)
    if m is None or sid is None:
        return literal
    digits = m.group(1)
    spec = f"%1$0{len(digits)}d" if digits.startswith("0") else "%1$d"
    head = name[:m.start(1)].replace("%", "%%")
    tail = name[m.end(1):].replace("%", "%%")
    pattern = head + spec + tail
    return pattern if matlab_sprintf(pattern, sid) == name else literal


def _t2_next_to(t1: Path) -> Path | None:
    for a, b in _T1_TAGS:
        if a in t1.name:
            cand = t1.with_name(t1.name.replace(a, b))
            if cand != t1 and cand.is_file():
                return cand
    return None


def _find_m2m(anat: Path, sid: int, siblings: bool) -> Path | None:
    """The folder that holds m2m_sub-XXX: the anatomy folder itself or one of its subfolders,
    and with `siblings` also a folder next to it (one level only, so this stays quick on network
    shares). Siblings are searched only when the file name gave the subject number: otherwise
    the ID is a guess, and a neighbour's m2m_sub-XXX would belong to someone else."""
    if m2m_dir(str(anat), sid).is_dir():
        return anat
    name = m2m_dir("", sid).name
    for base in ((anat, anat.parent) if siblings else (anat,)):
        try:
            subs = [p for p in base.iterdir() if p.is_dir()]
        except OSError:
            continue
        for d in sorted(subs):
            if (d / name).is_dir():
                return d
    return None


def from_t1(t1: str | Path, current_subject: int | None = None) -> T1Inputs:
    """Propose path.anat, path.t1_pattern, subject_id, path.t2_pattern and path.seg for `t1`."""
    t1 = Path(t1).expanduser().resolve()
    anat = t1.parent
    m = _SUB.search(t1.name)
    sid = int(m.group(1)) if m else None
    pattern = _pattern_for(t1.name, sid)
    notes: list[str] = []
    check_sid = sid if sid is not None else (current_subject or 0)
    if [p.resolve() for p in resolve_pattern(str(anat), pattern, check_sid)] != [t1]:
        pattern = t1.name.replace("%", "%%")          # the generic pattern matches other files too
    t2 = _t2_next_to(t1)
    t2_pattern = _pattern_for(t2.name, sid) if t2 is not None else ""
    if t2 is not None:
        notes.append(f"T2 found: {t2.name}")
    seg_root = _find_m2m(anat, check_sid, siblings=sid is not None)
    seg = "" if seg_root in (None, anat) else str(seg_root)
    return T1Inputs(str(anat), pattern, sid, t2_pattern, seg, seg_root is not None, notes)


def default_output(anat: str, config_path: Path | None) -> Path:
    """Where results go when the user has not said: next to the study config if it is saved,
    else next to the anatomy folder (never inside it, which may be a read-only dataset)."""
    if config_path is not None:
        return Path(config_path).resolve().parent / "prestus_sim"
    return Path(anat).parent / "prestus_sim"


def is_simnibs_bin(folder: str | Path) -> bool:
    return simnibs.is_simnibs_bin(folder)


def find_simnibs_bin() -> Path | None:
    """A SimNIBS 4 bin folder (one that holds charm) on this computer, or None: on PATH first,
    then the newest install found (simnibs.find_installs)."""
    installs = [i for i in simnibs.find_installs() if i.supported is not False]
    exe = shutil.which("charm")
    if exe:
        on_path = Path(exe).resolve().parent
        for inst in installs:
            if inst.bin.resolve() == on_path:
                return inst.bin
    return installs[0].bin if installs else None
