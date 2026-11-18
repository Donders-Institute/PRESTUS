"""Which tissues of the acoustic property table the subject's segmentation gives voxels to.

A layered medium assigns properties through ``layers`` (SimNIBS labels per tissue).
check_layers.m drops a layer none of whose labels occur in the segmentation and, with
pct.enabled, first merges skull_cortical and skull_trabecular into the single skull layer;
water takes every label no remaining layer claims. This module applies the same rules to the
whole-head segmentation (<path.seg>/m2m_sub-XXX/final_tissues.nii.gz, segmentation_file.m).

PRESTUS runs check_layers on the segmentation cropped to the simulation grid, so a tissue
present in the head can still be absent from the simulated domain: 'modelled' here is an
upper bound. The label volumes are counted once per file (cached by path, mtime and size);
nothing is written next to the segmentation.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

SKULL_SPLIT = ("skull_cortical", "skull_trabecular")
SEG_FILE = "final_tissues.nii.gz"

_VOLUMES: dict[tuple, dict[int, float]] = {}


@dataclass(frozen=True)
class TissueUse:
    state: str                  # 'modelled' | 'absent' | 'merged' | 'no_layer' | 'remainder'
    volume_ml: float | None     # voxels of its labels in the segmentation [mL]
    labels: tuple[int, ...]     # the SimNIBS labels it takes (after the pCT merge)
    note: str                   # why, in one sentence

    @property
    def modelled(self) -> bool:
        return self.state in ("modelled", "remainder")


def seg_file(m2m: str | Path | None) -> Path | None:
    return Path(m2m) / SEG_FILE if m2m is not None else None


def _labels(v) -> list[int]:
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        return [int(v)]
    if isinstance(v, (list, tuple)):
        return [int(x) for x in v if isinstance(x, (int, float)) and not isinstance(x, bool)]
    return []


def count_labels(data: np.ndarray, voxel_ml: float) -> dict[int, float]:
    """{label: volume [mL]} of the labels present in a label volume."""
    flat = np.asarray(data).ravel()
    if flat.dtype.kind in "iu" and (flat.size == 0 or flat.min() >= 0):
        counts = np.bincount(flat)
        return {int(l): float(n) * voxel_ml for l, n in enumerate(counts) if n}
    labels, counts = np.unique(np.rint(flat).astype(np.int64), return_counts=True)
    return {int(l): float(n) * voxel_ml for l, n in zip(labels, counts)}


def label_volumes(path: str | Path) -> dict[int, float]:
    """{label: volume [mL]} of a segmentation NIfTI, cached while the file is unchanged.
    Reads the whole volume: call it off the GUI thread."""
    p = Path(path)
    st = p.stat()
    key = (str(p), st.st_mtime, st.st_size)
    if key not in _VOLUMES:
        from .placement import Volume
        vol = Volume(p)
        voxel_ml = abs(float(np.linalg.det(vol.affine[:3, :3]))) / 1000.0
        if len(_VOLUMES) > 8:
            _VOLUMES.clear()
        _VOLUMES[key] = count_labels(vol.data, voxel_ml)
    return _VOLUMES[key]


def cached_volumes(path: str | Path) -> dict[int, float] | None:
    """The cached count for this file as it is now, else None (no read)."""
    p = Path(path)
    try:
        st = p.stat()
    except OSError:
        return None
    return _VOLUMES.get((str(p), st.st_mtime, st.st_size))


def tissue_use(layers: dict, medium_names: list[str], volumes: dict[int, float],
               pct: bool) -> dict[str, TissueUse]:
    """For each medium_properties tissue: does the segmentation give it voxels under these
    layers (check_layers.m, on the whole head)? Label 0 (background) is not counted in
    water's volume, which is the head's remaining labels only."""
    lay = {str(k): _labels(v) for k, v in (layers or {}).items()}
    merged_from: list[str] = []
    if pct:
        split = [n for n in SKULL_SPLIT if n in lay]
        if split:
            merged_from = split
            lay["skull"] = sorted(set(lay.get("skull", [])).union(*(lay.pop(n) for n in split)))
    present = {l for l, v in volumes.items() if v > 0}
    out: dict[str, TissueUse] = {}
    claimed: set[int] = set()
    for name in medium_names:
        if name == "water":
            continue
        if name in merged_from:
            out[name] = TissueUse("merged", None, (), "With pct.enabled, cortical and trabecular "
                                  "skull merge into the single skull layer, whose properties the "
                                  "(pseudo-)CT maps.")
            continue
        if name not in lay:
            out[name] = TissueUse("no_layer", None, (), "No entry in layers assigns labels to it.")
            continue
        labs = lay[name]
        hit = [l for l in labs if l in present]
        if hit:
            claimed.update(labs)
            vol = sum(volumes[l] for l in hit)
            note = f"Labels {', '.join(map(str, hit))} of the segmentation."
            if name == "skull" and merged_from:
                note += " Cortical and trabecular skull merged into it (pct.enabled)."
            out[name] = TissueUse("modelled", vol, tuple(labs), note)
        else:
            out[name] = TissueUse("absent", 0.0, tuple(labs),
                                  f"None of its labels ({', '.join(map(str, labs)) or 'none'}) "
                                  "occur in the segmentation; PRESTUS drops the layer.")
    if "water" in medium_names:
        rest = sorted(l for l in present - claimed if l != 0)
        out["water"] = TissueUse("remainder", sum(volumes[l] for l in rest), tuple(rest),
                                 "Everything no other layer claims: the coupling medium outside "
                                 "the head" + (f", and labels {', '.join(map(str, rest))} inside it."
                                               if rest else "."))
    return out


def volume_text(ml: float | None) -> str:
    if ml is None:
        return ""
    return f"{ml / 1000:.2f} L" if ml >= 1000 else f"{ml:.0f} mL" if ml >= 1 else f"{ml:.1f} mL"
