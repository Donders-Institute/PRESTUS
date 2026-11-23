"""Study manifest rows from Localite ingestion outputs (neuronav_gui / neuronav_ingest_batch).

Reads what the MATLAB side writes under path.localite_post:

  neuronav_ingest_summary.csv            one row per ingested target (batch window)
  sub-XXX/ses-XX/InstrumentMarker_sub-XXX_ses-XX.json   one entry per written target
                                         (single-session "Write outputs" and batch)

and builds localite-mode runs that PRESTUS resolves itself at run time
(placement.mode = localite, placement.localite.session / .target), either one
run per target (one transducer), or one run per session with the session's
targets assigned to transducers 1..N (bilateral / multi-transducer setups).
"""

from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .study import Manifest
from .ui.errors import log_current

_SUB = re.compile(r"^sub-(\d+)$")


def _f(row: dict, key: str) -> float | None:
    """A summary-CSV numeric cell: MATLAB writetable writes NaN as the literal 'NaN'."""
    v = (row.get(key) or "").strip()
    if not v or v.upper() == "NAN":
        return None
    try:
        return float(v)
    except ValueError:
        return None


def _vec3(row: dict, prefix: str) -> np.ndarray | None:
    vals = [_f(row, f"{prefix}_{a}") for a in ("x", "y", "z")]
    return None if any(v is None for v in vals) else np.array(vals)


@dataclass(frozen=True)
class Entry:
    sub: str        # 'sub-001'
    ses: str        # 'ses-02'
    target: str
    source: str     # 'summary' | 'json'
    # QC and coordinates, from neuronav_ingest_summary.csv only (summary_vars() in
    # neuronav_ingest_batch.m); None for 'json'-sourced entries or missing columns.
    pos_dev_mm: float | None = None
    rot_dev_rad: float | None = None
    axis_dist_mm: float | None = None
    depth_mm: float | None = None
    n_pulses: float | None = None
    trans_ras: np.ndarray | None = None     # RAS mm, bowl apex
    focus_ras: np.ndarray | None = None     # RAS mm

    @property
    def subject_id(self) -> str:
        m = _SUB.match(self.sub)
        return str(int(m.group(1))) if m else self.sub


def scan(localite_post: str | Path) -> list[Entry]:
    """All ingested (subject, session, target) under localite_post."""
    root = Path(localite_post)
    found: dict[tuple[str, str, str], Entry] = {}
    summary = root / "neuronav_ingest_summary.csv"
    if summary.is_file():
        with summary.open(encoding="utf-8", newline="") as f:
            for r in csv.DictReader(f):
                t = (r.get("target") or "").strip()
                if t:
                    e = Entry(r["sub"].strip(), r["ses"].strip(), t, "summary",
                              pos_dev_mm=_f(r, "pos_dev_mm"), rot_dev_rad=_f(r, "rot_dev_rad"),
                              axis_dist_mm=_f(r, "axis_dist_mm"), depth_mm=_f(r, "depth_mm"),
                              n_pulses=_f(r, "n_pulses"), trans_ras=_vec3(r, "trans_ras"),
                              focus_ras=_vec3(r, "focus_ras"))
                    found[(e.sub, e.ses, e.target)] = e
    for js in sorted(root.glob("sub-*/ses-*/InstrumentMarker_*.json")):
        sub, ses = js.parent.parent.name, js.parent.name
        if any(k[0] == sub and k[1] == ses for k in found):
            continue          # the summary has this session (with the real target names)
        try:
            data = json.loads(js.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):     # this session is not offered
            log_current(f"neuronav import: cannot read {js}")
            continue
        if isinstance(data, dict):
            for t in data:
                found[(sub, ses, t)] = Entry(sub, ses, t, "json")
    return sorted(found.values(), key=lambda e: (e.sub, e.ses, e.target))


def coordinates_for(localite_post: str | Path, session: str, target: str,
                    subject_id: int | None = None) -> Entry | None:
    """The scan() entry for one (session, target) [, subject], for the placement viewer.

    session/target match placement.localite.session/.target as PRESTUS resolves them
    (ses_label / the raw target string); session may be given with or without 'ses-'."""
    ses = session if str(session).startswith("ses-") else f"ses-{int(session):02d}"
    sub = f"sub-{subject_id:03d}" if subject_id is not None else None
    for e in scan(localite_post):
        if e.ses == ses and e.target == target and (sub is None or e.sub == sub):
            return e
    return None


def order_targets(targets: list[str], rule: str | list[str]) -> list[str]:
    """'lr': left (L…) targets before right (R…); 'alphabetical'; or an explicit order."""
    targets = sorted(set(targets))
    if isinstance(rule, list):
        pos = {t: i for i, t in enumerate(rule)}
        return sorted(targets, key=lambda t: (pos.get(t, len(pos)), t))
    if rule == "lr":
        side = lambda t: 0 if re.match(r"^(L|left)", t, re.I) else 1 if re.match(r"^(R|right)", t, re.I) else 2  # noqa: E731
        return sorted(targets, key=lambda t: (side(t), t))
    return targets


def build_manifest(entries: list[Entry], grouping: str = "target", order: str | list[str] = "lr",
                   n_transducers: int | None = None) -> tuple[Manifest, list[str]]:
    """Manifest rows (localite placement) and notes about skipped sessions."""
    notes: list[str] = []
    m = Manifest(columns=["subject_id", "condition", "placement.mode", "placement.localite.session",
                          "placement.localite.target"])
    if grouping == "target":
        count: dict[tuple[str, str], int] = {}
        for e in entries:
            count[(e.sub, e.target)] = count.get((e.sub, e.target), 0) + 1
        for e in entries:
            cond = e.target if count[(e.sub, e.target)] == 1 else f"{e.target}-{e.ses}"
            m.add_row({"subject_id": e.subject_id, "condition": cond, "placement.mode": "localite",
                       "placement.localite.session": e.ses, "placement.localite.target": e.target})
        return m, notes

    sessions: dict[tuple[str, str], list[str]] = {}
    for e in entries:
        sessions.setdefault((e.sub, e.ses), []).append(e.target)
    for (sub, ses), targets in sessions.items():
        ordered = order_targets(targets, order)
        if isinstance(order, list):
            missing = [t for t in order if t not in ordered]
            if missing:
                notes.append(f"{sub} {ses}: no {', '.join(missing)}; skipped.")
                continue
            ordered = [t for t in ordered if t in order]
        if n_transducers and len(ordered) != n_transducers:
            notes.append(f"{sub} {ses}: {len(ordered)} target(s) ({', '.join(ordered)}) but "
                         f"{n_transducers} transducer(s) configured; skipped.")
            continue
        subject_id = Entry(sub, ses, "", "").subject_id
        m.add_row({"subject_id": subject_id, "condition": ses, "placement.mode": "localite",
                   "placement.localite.session": ses,
                   "placement.localite.target": "[" + ", ".join(ordered) + "]" if len(ordered) > 1 else ordered[0]})
    return m, notes
