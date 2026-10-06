"""Localite page logic without Qt: where the raw exports are, which sessions exist, which
config values the ingestion depends on, and the target_map written for the ticked series."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .model import ConfigModel, is_empty

# Config values neuronav_detect_series / neuronav_name_series / neuronav_ingest_markers read;
# a detection is stale once one of them changes.
SETTINGS_KEYS = (
    "path.localite_raw", "path.localite", "path.localite_post",
    "placement.localite.loc_folder", "placement.localite.session_folder_pattern",
    "placement.localite.expected_segment_length", "placement.localite.voxel_size",
    "io.overwrite_files", "neuronav.coil_map", "neuronav.suggest_names",
    "neuronav.axis_tolerance_mm", "neuronav.target_depth_mm", "neuronav.target_labels",
)


def raw_root(model: ConfigModel) -> str:
    """path.localite_raw, else the deprecated path.localite (neuronav_resolve_paths.m)."""
    for key in ("path.localite_raw", "path.localite"):
        v = model.get(key)
        if not is_empty(v):
            return str(v)
    return ""


def post_root(model: ConfigModel) -> str:
    v = model.get("path.localite_post")
    return "" if is_empty(v) else str(v)


def list_sessions(raw: str | Path) -> dict[str, list[str]]:
    """sub-* folders and their ses-* folders under the raw root."""
    root = Path(raw) if raw else None
    if root is None or not root.is_dir():
        return {}
    out = {}
    for sub in sorted(p for p in root.glob("sub-*") if p.is_dir()):
        out[sub.name] = sorted(p.name for p in sub.glob("ses-*") if p.is_dir())
    return out


def settings_snapshot(model: ConfigModel) -> tuple:
    return tuple(repr(model.get(k)) for k in SETTINGS_KEYS)


def settings_summary(model: ConfigModel) -> list[tuple[str, str]]:
    """(label, value) pairs shown above the page; edited on the Settings tab."""
    seg = model.get("placement.localite.expected_segment_length")
    coil = model.get("neuronav.coil_map")
    labels = model.get("neuronav.target_labels")
    naming = "planned target nearest the beam axis" if model.get("neuronav.suggest_names") else "series_<k>"
    return [
        ("Raw exports", raw_root(model) or "— (path.localite_raw)"),
        ("Ingested markers", post_root(model) or "— (path.localite_post)"),
        ("Localite folder", str(model.get("placement.localite.loc_folder") or "(directly under ses-XX)")),
        ("Pulses per series", "accept all" if is_empty(seg) else f"{seg} (±20)"),
        ("Coils", f"{len(coil)} in neuronav.coil_map" if isinstance(coil, list) and coil else "single coil"),
        ("Names", naming + (f", {len(labels)} label rule(s)" if isinstance(labels, list) and labels else "")),
        ("Overwrite", str(model.get("io.overwrite_files"))),
    ]


def target_map(series: list[dict], ticked: list[bool], names: list[str], colours: list[str]) -> list[dict]:
    """neuronav_ingest_markers target_map for the ticked, non-empty series (checked in
    MATLAB by neuronav_check_target_map)."""
    out = []
    for s, use, name, colour in zip(series, ticked, names, colours):
        if use and not s.get("empty"):
            out.append({"series_index": int(s["index"]), "transducer_id": int(s.get("transducer_id") or 0),
                        "target_name": name.strip(), "color": colour.strip()})
    return out


def use_in_study(model: ConfigModel, ses: str, targets: list[str]) -> list[tuple[str, Any]]:
    """Point the study's single run at this session and target(s) (placement.localite.session /
    .target; several targets = one per transducer). Returns the values set."""
    changes: list[tuple[str, Any]] = [("placement.localite.session", ses)]
    if targets:
        changes.append(("placement.localite.target", targets[0] if len(targets) == 1 else list(targets)))
    for key, value in changes:
        model.set(key, value)
    return changes
