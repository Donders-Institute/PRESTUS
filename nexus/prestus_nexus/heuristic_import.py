"""Reads what transducer_positioning.m / tp_select_heuristic_position.m already wrote
to disk for heuristic placement (doc_placement_heuristic.md), so the placement viewer
can show candidates and the selected position without recomputing anything.

tpars_sub-XXX_<target>.csv  - candidate positions (transducer_positioning.m:61-96),
                              voxel space, one row per candidate on the skull surface.
sub-XXX_<target>.txt        - the selected position (tp_select_heuristic_position.m:170),
                              a key/value TSV (no header) with voxel, Localite mm and
                              RAS mm coordinates for the transducer and the target.
Both live in the heuristic's placement run folder (placement_run_dir.m),
i.e. outputs.dir_placement_run(model, "heuristic").
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Candidate:
    idx: int
    trans_vox: np.ndarray        # [x, y, z], voxel
    targ_vox: np.ndarray
    dist_to_target: float
    prop_intersect: float        # fraction of the beam path intersecting excluded tissue; low is good
    mean_dist_skin: float
    var_dist_skin: float
    mean_dist_skull: float
    var_dist_skull: float


@dataclass(frozen=True)
class Selected:
    trans_vox: np.ndarray | None
    targ_vox: np.ndarray | None
    trans_ras: np.ndarray | None    # RAS mm, bowl apex
    focus_ras: np.ndarray | None    # RAS mm


def _candidates_path(run_dir: Path, subject_id: int, target_name: str) -> Path:
    return Path(run_dir) / f"tpars_sub-{subject_id:03d}_{target_name}.csv"


def _selected_path(run_dir: Path, subject_id: int, target_name: str) -> Path:
    return Path(run_dir) / f"sub-{subject_id:03d}_{target_name}.txt"


def candidates(run_dir: str | Path, subject_id: int, target_name: str) -> list[Candidate]:
    path = _candidates_path(Path(run_dir), subject_id, target_name)
    if not path.is_file():
        return []
    out = []
    with path.open(encoding="utf-8", newline="") as f:
        for i, r in enumerate(csv.DictReader(f)):
            try:
                out.append(Candidate(
                    idx=int(float(r.get("idx", i))),
                    trans_vox=np.array([float(r["trans_x"]), float(r["trans_y"]), float(r["trans_z"])]),
                    targ_vox=np.array([float(r["targ_x"]), float(r["targ_y"]), float(r["targ_z"])]),
                    dist_to_target=float(r["dist_to_target"]),
                    prop_intersect=float(r["prop_intersect"]),
                    mean_dist_skin=float(r["mean_dist_skin"]),
                    var_dist_skin=float(r["var_dist_skin"]),
                    mean_dist_skull=float(r["mean_dist_skull"]),
                    var_dist_skull=float(r["var_dist_skull"]),
                ))
            except (KeyError, ValueError):
                continue
    return out


def selected(run_dir: str | Path, subject_id: int, target_name: str) -> Selected | None:
    path = _selected_path(Path(run_dir), subject_id, target_name)
    if not path.is_file():
        return None
    kv: dict[str, float] = {}
    with path.open(encoding="utf-8", newline="") as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) < 2:
                continue
            try:
                kv[row[0].strip()] = float(row[1])
            except ValueError:
                continue       # section-header rows (e.g. "RAS_Transducer_mm") and blanks

    def vec(keys: tuple[str, str, str]) -> np.ndarray | None:
        vals = [kv.get(k) for k in keys]
        return None if any(v is None for v in vals) else np.array(vals)

    trans_ras = vec(("ras_trans_x_mm", "ras_trans_y_mm", "ras_trans_z_mm"))
    focus_ras = vec(("ras_targ_x_mm", "ras_targ_y_mm", "ras_targ_z_mm"))
    trans_vox = vec(("trans_x", "trans_y", "trans_z"))
    targ_vox = vec(("targ_x", "targ_y", "targ_z"))
    if trans_ras is None and trans_vox is None:
        return None
    return Selected(trans_vox=trans_vox, targ_vox=targ_vox, trans_ras=trans_ras, focus_ras=focus_ras)
