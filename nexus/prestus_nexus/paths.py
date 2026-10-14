"""Locating PRESTUS, its config_default.yaml, and subject files."""

from __future__ import annotations

import glob
import os
import re
from pathlib import Path

DEFAULT_NAME = "config_default.yaml"


def is_prestus_root(path: str | Path) -> bool:
    p = Path(path)
    return (p / "config" / DEFAULT_NAME).is_file() and (p / "functions").is_dir()


def package_dir() -> Path:
    """This package's folder in the PRESTUS checkout (nexus/prestus_nexus), also when the
    Nexus runs from a local copy of it (localcopy.py), which sets $PRESTUS_NEXUS_SOURCE."""
    source = os.environ.get("PRESTUS_NEXUS_SOURCE")
    return Path(source) if source else Path(__file__).resolve().parent


def find_prestus_root(hint: str | Path | None = None) -> Path | None:
    """PRESTUS root: explicit hint, $PRESTUS_ROOT, or the checkout this package lives in."""
    candidates = [hint, os.environ.get("PRESTUS_ROOT"), package_dir().parents[1]]
    for cand in candidates:
        if cand and is_prestus_root(cand):
            return Path(cand).resolve()
    return None


def toolbox_default(prestus_root: str | Path) -> Path:
    return Path(prestus_root) / "config" / DEFAULT_NAME


def base_default_for(study_file: str | Path | None, prestus_root: str | Path) -> Path:
    """config_default.yaml that load_parameters(file, dir) will merge onto.

    The launcher calls load_parameters(<file>, <dir>), which prefers a
    config_default.yaml in <dir> over the toolbox copy.
    """
    if study_file:
        local = Path(study_file).resolve().parent / DEFAULT_NAME
        if local.is_file():
            return local
    return toolbox_default(prestus_root)


# ── MATLAB sprintf subset used by path.t1_pattern / t2_pattern ───────────

_FMT = re.compile(r"%(?:(\d+)\$)?([-+ 0#]*)(\d*)(?:\.(\d+))?([diuscfg%])")


def matlab_sprintf(fmt: str, *args) -> str:
    """Enough of MATLAB sprintf for subject patterns ('%03d', '%1$03d', '%s')."""
    counter = iter(range(len(args)))

    def repl(m: re.Match) -> str:
        pos, flags, width, prec, conv = m.groups()
        if conv == "%":
            return "%"
        idx = int(pos) - 1 if pos else next(counter, len(args))
        if idx >= len(args):
            return ""
        spec = "%" + flags + width + (f".{prec}" if prec else "") + ("d" if conv in "iu" else conv)
        return spec % args[idx]

    return _FMT.sub(repl, fmt)


def resolve_pattern(directory: str, pattern: str, subject_id: int) -> list[Path]:
    """Files matching fullfile(directory, sprintf(pattern, subject_id)), like dir()."""
    if not directory or not pattern:
        return []
    expanded = os.path.join(directory, matlab_sprintf(pattern, subject_id))
    return sorted(Path(p) for p in glob.glob(expanded))


def m2m_dir(seg_path: str, subject_id: int) -> Path:
    return Path(seg_path) / f"m2m_sub-{subject_id:03d}"
