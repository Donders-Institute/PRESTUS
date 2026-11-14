"""Follow-up runs: a new study config that starts from the results of an earlier run.

A follow-up is a copy of the parent config's overrides plus the changes of one follow-up
type, saved as its own study YAML. The child records where it came from in a comment header
line (``# nexus-followup: {...}``, model.LINEAGE_TAG); MATLAB ignores comments, so no config
key is added. Qt-free.

Types and the keys each writes (verified against the MATLAB code):

- ``new_intensity``: ``transducer.target_isppa_wcm2`` and ``path.acoustic_input`` = the
  parent's T1w pressure map ``nii/<tag>_<medium>_T1w<affix>_pressure.nii.gz``
  (write_acoustic_sidecar.m:33-34). The run skips the acoustic simulation
  (acoustic_stage.m:84-90) and the free-water baseline (acoustic_water_baseline.m:33-36) and
  rescales the map with the sidecar's ``isppa.level_wcm2`` (load_acoustic_input.m:119-131).
- ``new_timing`` (thermal only): as new_intensity (target optional) plus ``timing.*`` and
  ``modules.run_heating_sims: 1``. Nothing else is simulated acoustically; the post-hoc water
  run is switched off (``modules.run_posthoc_water_sims: 0``) because posthoc_water_parameters.m
  keeps path.acoustic_input on a water-medium run, which has no T1 transform
  (load_acoustic_input.m:53-57).
- ``new_placement_same_head``: ``headmodel.reuse_across_placements: true`` plus placement
  keys (head_model_load_or_build.m: one head model per subject, rotated per placement).
- ``pose_from_time_reversal``: a row of ``<tag>_<medium><affix>_timereversal_poses_T<i>.csv``
  (time_reversal_stage.m:369,387) as a manual placement in RAS mm, as
  time_reversal_pose_parameters.m does it; or ``time_reversal.apply_phases: true``.
- ``sequential``: the child continues from the parent's thermal end state; it is added to the
  parent's launcher options (``options.sequential_configs``, sequential_pipeline.m), not run
  on its own.

Every type sets a new ``io.output_affix``.
"""

from __future__ import annotations

import copy
import csv
import datetime as _dt
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .model import TX, ConfigModel, Issue, is_empty
from .study import sanitize_affix

# ── types ──────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class FollowupType:
    key: str
    label: str
    blurb: str                  # one line, shown under the type chooser
    uses_pressure: bool         # starts from the parent's pressure map (path.acoustic_input)


TYPES: dict[str, FollowupType] = {t.key: t for t in (
    FollowupType("new_intensity", "New intensity",
                 "Same beam at another target ISPPA: reuses this run's pressure map, no new "
                 "acoustic simulation.", True),
    FollowupType("new_timing", "New timing (thermal only)",
                 "Other pulse timing on the same beam: reuses the pressure map and runs only "
                 "the thermal simulation.", True),
    FollowupType("new_placement_same_head", "New placement, same head",
                 "Another transducer position on this head: the head model is built once and "
                 "reused.", False),
    FollowupType("pose_from_time_reversal", "Pose from the time-reversal search",
                 "Simulate a pose the time-reversal pose search suggested, or apply its "
                 "refocusing phases.", False),
    FollowupType("sequential", "Sequential (heat carry-over)",
                 "A sonication that follows this one and starts from its end temperature and "
                 "thermal dose.", False),
)}

POSE_AIMS = ("corrected", "geometric", "apply_phases")

# time_reversal_pose_parameters.m:52-53: the RAS columns a pose needs
POSE_RAS_COLUMNS = ("apex_ras_x", "apex_ras_y", "apex_ras_z", "aim_ras_x", "aim_ras_y", "aim_ras_z",
                    "target_ras_x", "target_ras_y", "target_ras_z")


# ── lineage ────────────────────────────────────────────────────────────────────

def make_lineage(parent_model: ConfigModel, parent_run, ftype: str) -> dict:
    """The lineage record of a follow-up (model.lineage)."""
    parent = parent_model.path.resolve() if parent_model.path is not None else None
    return {"parent_config": str(parent) if parent else "",
            "parent_affix": getattr(parent_run, "affix", "") or "",
            "type": ftype,
            "created": _dt.datetime.now().replace(microsecond=0).isoformat()}


# ── the parent run's files ─────────────────────────────────────────────────────

def pressure_nifti(run) -> Path:
    """T1w pressure map of a run (write_acoustic_sidecar.m:33-34, nifti_acoustic.m)."""
    return run.base / "nii" / f"{run.tag}_{run.medium}_T1w{run.affix}_pressure.nii.gz"


def sidecar_of(data_file: Path) -> Path:
    """JSON sidecar of a data file: same name with .json (find_sidecar.m:22-23)."""
    stem = re.sub(r"(\.nii\.gz|\.nii|\.mat|\.h5)$", "", str(data_file), flags=re.IGNORECASE)
    return Path(stem + ".json")


def thermal_end_maps(run) -> list[Path]:
    """The T1w end-state maps a sequential follow-up starts from (sequential_pipeline.m:60-70)."""
    return [run.base / "nii" / f"{run.tag}_{run.medium}_T1w{run.affix}_{name}.nii.gz"
            for name in ("heating_end", "CEM43_end", "CEM43_iso_end")]


def pose_tables(run) -> list[Path]:
    """The run's time-reversal pose tables, one per transducer (time_reversal_stage.m:369,387;
    io.dir_tabular is the output folder, path_log_setup.m:128)."""
    stem = f"{run.tag}_{run.medium}{run.affix}_timereversal_poses_T"
    try:
        names = sorted(p.name for p in run.base.iterdir()
                       if p.name.startswith(stem) and p.name.endswith(".csv"))
    except OSError:
        return []
    out = []
    for name in names:
        if name[len(stem):-len(".csv")].isdigit():
            out.append(run.base / name)
    return sorted(out, key=lambda p: int(p.name[len(stem):-len(".csv")]))


def _cell(text: str) -> Any:
    """A CSV cell as MATLAB's writetable writes it: number, logical (true/false/1/0) or text."""
    t = text.strip()
    if t.lower() in ("true", "false"):
        return t.lower() == "true"
    if t.lower() == "nan" or t == "":
        return math.nan if t else None
    try:
        return int(t) if re.fullmatch(r"[-+]?\d+", t) else float(t)
    except ValueError:
        return t


def read_pose_table(path: Path) -> list[dict]:
    """Rows of a pose table as dicts, in file order (sorted by p_corrected_rel, best first:
    reciprocity_score_poses.m:124)."""
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    if not rows:
        return []
    head = [h.strip() for h in rows[0]]
    return [{k: _cell(v) for k, v in zip(head, r)} for r in rows[1:] if any(c.strip() for c in r)]


def read_sidecar(path: Path) -> dict | None:
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


# ── sidecar signature vs the child's settings (mirrors check_sidecar_signature.m) ──

# Transducer fields PRESTUS derives at run time (placement, library, validators), so the
# config value cannot be compared with the recorded one.
_DERIVED_TX = {"trans_pos", "focus_pos", "focal_distance_ep", "focal_distance_bowl",
               "focal_distance_offset", "align_to_focus", "elem_phase_deg", "elem_phase_rad",
               "elem_p_amp", "target_isppa_wcm2", "tx_calibration"}


def _norm_leaf(v: Any) -> Any:
    if v is None or v == "" or v == [] or v == {}:
        return None
    if isinstance(v, bool):
        return float(v)
    if isinstance(v, (int, float)):
        return float(v)
    if isinstance(v, list):
        flat: list = []

        def walk(x):
            if isinstance(x, list):
                for y in x:
                    walk(y)
            else:
                flat.append(_norm_leaf(x))
        walk(v)
        return flat[0] if len(flat) == 1 else tuple(flat)
    return v


def _same(a: Any, b: Any) -> bool:
    a, b = _norm_leaf(a), _norm_leaf(b)
    if isinstance(a, float) and isinstance(b, float):
        if math.isnan(a) and math.isnan(b):
            return True
        return math.isclose(a, b, rel_tol=1e-6, abs_tol=1e-9)
    if isinstance(a, tuple) and isinstance(b, tuple):
        return len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b))
    return a == b


def _flatten(d: Any, prefix: str = "") -> dict[str, Any]:
    out: dict[str, Any] = {}
    if isinstance(d, dict):
        for k, v in d.items():
            out.update(_flatten(v, f"{prefix}.{k}" if prefix else k))
    elif isinstance(d, list) and d and all(isinstance(x, dict) for x in d):
        for i, x in enumerate(d):
            out.update(_flatten(x, f"{prefix}({i + 1})"))
    else:
        out[prefix] = d
    return out


def signature_diffs(stored: dict, effective: dict) -> list[str]:
    """Dotted keys whose recorded value (the sidecar's ``signature``, cache_signature.m kind
    'acoustic') differs from the child's configured one. Values PRESTUS derives at run time
    (positions, library drive, 'auto' sizes) are skipped, so this is a preview of
    check_sidecar_signature.m, which stays the authority."""
    if not isinstance(stored, dict):
        return []
    eff = copy.deepcopy(effective)
    mine: dict[str, Any] = {"medium": (eff.get("simulation") or {}).get("medium"),
                            "subject_id": eff.get("subject_id")}
    for sec in ("grid", "layers", "headmodel", "pct"):
        mine[sec] = eff.get(sec)
    if isinstance(mine["grid"], dict):          # cache_signature.m:91-94: thermal_* not compared
        mine["grid"] = {k: v for k, v in mine["grid"].items() if not k.startswith("thermal_")}
    mp = eff.get("medium_properties")
    if isinstance(mp, dict):                    # cache_signature.m:89 (drop_absorption_coeff, [])
        mine["medium_properties"] = {t: ({k: v for k, v in p.items() if k != "absorption_coeff"}
                                         if isinstance(p, dict) else p) for t, p in mp.items()}
    tx = eff.get(TX)
    mine[TX] = tx if isinstance(tx, list) else [tx]
    stored = dict(stored)
    if isinstance(stored.get(TX), dict):
        stored[TX] = [stored[TX]]
    have, want = _flatten(stored), _flatten(mine)
    diffs = []
    for key, value in have.items():
        if key not in want or key.startswith("signature_version"):
            continue
        last = key.rsplit(".", 1)[-1]
        if key.startswith(TX) and last in _DERIVED_TX:
            continue
        mine_v = want[key]
        if isinstance(mine_v, str) and isinstance(_norm_leaf(value), (float, tuple)):
            continue                            # 'auto' in the config, a number at run time
        if not _same(value, mine_v):
            diffs.append(key)
    return sorted(diffs)


# ── affixes and file names ─────────────────────────────────────────────────────

def _raw_affix(affix: str) -> str:
    """'_desc-pilot' -> 'pilot' (the value io.output_affix holds before sanitising)."""
    return affix[len("_desc-"):] if affix.startswith("_desc-") else affix.lstrip("_")


def _num_label(x: float) -> str:
    return f"{x:g}".replace(".", "p").replace("-", "m")


def default_affix(parent_affix: str, ftype: str, params: dict | None = None) -> str:
    """io.output_affix for the child: the parent's label plus what changed."""
    params = params or {}
    if ftype == "new_intensity":
        v = params.get("target_isppa_wcm2")
        part = f"isppa{_num_label(float(v))}" if _num(v) else "isppa"
    elif ftype == "new_timing":
        part = "timing"
    elif ftype == "new_placement_same_head":
        part = "placement"
    elif ftype == "pose_from_time_reversal":
        aim = params.get("aim", "corrected")
        if aim == "apply_phases":
            part = "trphases"
        else:
            pose = (params.get("pose") or {}).get("pose")
            part = (f"pose{int(pose):03d}" if _num(pose) else "pose") + ("-geo" if aim == "geometric" else "")
    elif ftype == "sequential":
        part = f"seq{int(params.get('seq_index', 2))}"
    else:
        part = ftype
    base = _raw_affix(parent_affix)
    return f"{base}_{part}" if base else part


def followup_filename(parent_config: Path, ftype: str, affix: str) -> str:
    """<parent_stem>_<type>_<affix>.yaml; lowercase, hyphens within fields."""
    def field(s: str) -> str:
        return re.sub(r"[^a-z0-9-]+", "-", s.lower().replace("_", "-")).strip("-")
    parts = [Path(parent_config).stem.lower(), field(ftype), field(_raw_affix(sanitize_affix(affix)))]
    return "_".join(p for p in parts if p) + ".yaml"


def _num(v: Any) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)


# ── create ─────────────────────────────────────────────────────────────────────

def copy_model(parent: ConfigModel) -> ConfigModel:
    child = ConfigModel(parent.schema, parent.library)
    child.load_dict(parent.override_dict())
    child.tx_as_list = parent.tx_as_list
    return child


def _effective_tx(model: ConfigModel, ti: int = 0) -> dict:
    txs = model.effective_transducers()
    return txs[ti] if ti < len(txs) else {}


def _apply_pose(child: ConfigModel, pose: dict, aim: str) -> list[str]:
    """time_reversal_pose_parameters.m: apex -> trans_pos, aim point (or target) -> focus_pos,
    manual placement in RAS mm, focal depth of the pose, time reversal off. Returns notes."""
    notes = []
    apex = [float(pose[f"apex_ras_{a}"]) for a in "xyz"]
    src = "aim" if aim == "corrected" else "target"
    focus = [float(pose[f"{src}_ras_{a}"]) for a in "xyz"]
    child.set("transducer.trans_pos", apex)
    child.set("transducer.focus_pos", focus)
    # focal depth: |focus - apex| minus the bowl-to-exit-plane offset (lines 79-84)
    ann = _effective_tx(child).get("annular") or {}
    offset = 0.0
    R, d = ann.get("curv_radius_mm"), ann.get("dist_geom_ep_mm")
    if _num(R) and _num(d):                     # annular_bowl_depth_mm.m
        offset = float(R) - float(d)
    child.set("transducer.focal_distance_ep", round(math.dist(apex, focus) - offset, 4))
    child.set("transducer.focal_distance_bowl", [])
    # inline ring phases: the pose run steers geometrically to the aim point (lines 96-101);
    # cleared phases let the focal depth set the steering (fill_annular_drive.m)
    for key in ("transducer.annular.elem_phase_deg", "transducer.annular.elem_phase_rad"):
        if child.is_overridden(key) and not is_empty(child.get(key)):
            child.set(key, [])
            notes.append(f"{key} cleared: the pose's focal depth sets the ring steering")
    child.set("placement.mode", "manual")
    child.set("placement.manual.coordinate_space", "ras")
    child.set("modules.run_transducer_placement", 1)
    for key in ("time_reversal.enabled", "time_reversal.only",
                "time_reversal.placement_search.enabled", "time_reversal.pose_workflow.enabled"):
        child.set(key, False)
    return notes


def create_followup(parent_model: ConfigModel, parent_run, ftype: str,
                    params: dict | None = None) -> ConfigModel:
    """A new ConfigModel: the parent's overrides plus the changes of ``ftype``, with the
    lineage header set. ``params``: ``affix`` (raw io.output_affix; default default_affix),
    ``target_isppa_wcm2``, ``timing`` ({'pd': ..}), ``pose`` (a pose table row) and ``aim``
    ('corrected' | 'geometric' | 'apply_phases'), ``seq_index``, and ``changes``
    ({dotted key: value}, applied last)."""
    if ftype not in TYPES:
        raise ValueError(f"Unknown follow-up type '{ftype}'")
    params = dict(params or {})
    child = copy_model(parent_model)
    child.path = None
    child.lineage = make_lineage(parent_model, parent_run, ftype)

    if TYPES[ftype].uses_pressure:
        child.set("path.acoustic_input", str(pressure_nifti(parent_run)))
        child.set("modules.run_posthoc_water_sims", 0)
        if _num(params.get("target_isppa_wcm2")):
            child.set("transducer.target_isppa_wcm2", float(params["target_isppa_wcm2"]))
    if ftype == "new_timing":
        for key, value in (params.get("timing") or {}).items():
            child.set(key if key.startswith("timing.") else f"timing.{key}", value)
        child.set("modules.run_heating_sims", 1)
    elif ftype == "new_placement_same_head":
        child.set("headmodel.reuse_across_placements", True)
    elif ftype == "pose_from_time_reversal":
        aim = params.get("aim", "corrected")
        if aim == "apply_phases":
            child.set("time_reversal.enabled", True)
            child.set("time_reversal.apply_phases", True)
        elif params.get("pose"):
            _apply_pose(child, params["pose"], aim)
    elif ftype == "sequential":
        # sequential_pipeline.m:80-88 needs the T1w end maps: heating and NIfTI stages on
        child.set("modules.run_heating_sims", 1)
        child.set("modules.run_nifti_creation", 1)

    for key, value in (params.get("changes") or {}).items():
        child.set(key, value)
    affix = params.get("affix")
    if not affix:
        affix = default_affix(getattr(parent_run, "affix", ""), ftype, params)
    child.set("io.output_affix", affix)
    return child


# ── validation ─────────────────────────────────────────────────────────────────

def _issue(level: str, path: str, message: str) -> Issue:
    return Issue(level, path, message)


def pressure_issues(model: ConfigModel, nifti: Path, need_level: bool) -> list[Issue]:
    """What PRESTUS checks when a run starts from a pressure map (load_acoustic_input.m,
    load_parameters.m:343-353). ``need_level``: the map must be rescalable (a new target)."""
    issues: list[Issue] = []
    key = "path.acoustic_input"
    g = model.get
    if not nifti.is_file():
        issues.append(_issue("error", key, f"Pressure map not found: {nifti} (PRESTUS stops with "
                                           "prestus:acousticInput:notFound)"))
    if g("simulation.medium") not in ("layered", "phantom"):
        issues.append(_issue("error", "simulation.medium",
                             "A pressure map needs the T1 → grid transform: medium layered or phantom "
                             "(load_acoustic_input.m:53-57)"))
    if g("grid.axisymmetric") in (1, True):
        issues.append(_issue("error", "grid.axisymmetric",
                             "A pressure map needs a 3D grid, not axisymmetric (load_acoustic_input.m:59-62)"))
    if g("thermal_backend") == "babelbrain":
        issues.append(_issue("error", "thermal_backend",
                             "thermal_backend 'babelbrain' cannot start from a pressure map "
                             "(load_parameters.m:343-353)"))
    txs = model.effective_transducers()
    if len(txs) > 1 and any(_num(t.get("target_isppa_wcm2")) or
                            (isinstance(t.get("target_isppa_wcm2"), list) and
                             any(_num(x) for x in t["target_isppa_wcm2"])) for t in txs):
        issues.append(_issue("error", "transducer.target_isppa_wcm2",
                             "Several transducers: a combined pressure map cannot be rescaled to "
                             "per-transducer targets; set target_isppa_wcm2: [] on every transducer "
                             "(load_acoustic_input.m:112-118)"))
    if g("modules.run_source_setup") != 1 or g("modules.run_grid_setup") != 1:
        issues.append(_issue("error", "modules.run_source_setup",
                             "A pressure map needs the grid and source stages (modules.run_grid_setup "
                             "and run_source_setup: 1)"))
    if g("modules.run_posthoc_water_sims") == 1:
        issues.append(_issue("warning", "modules.run_posthoc_water_sims",
                             "The post-hoc water run keeps path.acoustic_input (posthoc_water_parameters.m) "
                             "on a water medium, which cannot read a T1w map; set run_posthoc_water_sims: 0"))
    target = _effective_tx(model).get("target_isppa_wcm2")
    if isinstance(target, list) and len([x for x in target if _num(x)]) > 1:
        issues.append(_issue("warning", "transducer.target_isppa_wcm2",
                             "A sweep of targets with a pressure map: the map is rescaled only to a "
                             "scalar target (load_acoustic_input.m:119-122)"))
    if not nifti.is_file():
        return issues
    sc_path = sidecar_of(nifti)
    sc = read_sidecar(sc_path) if sc_path.is_file() else None
    if sc is None:
        msg = (f"No PRESTUS sidecar ({sc_path.name}): the map is used as it is, its settings are "
               "not checked")
        issues.append(_issue("error" if need_level and _num(target) else "warning", key,
                             msg + (" and it cannot be rescaled to the new target" if _num(target) else "")))
        return issues
    if sc.get("kind") != "acoustic_pressure":
        issues.append(_issue("warning", key, f"{sc_path.name} is not an acoustic_pressure sidecar; "
                                             "its settings are not checked"))
        return issues
    isppa = sc.get("isppa") if isinstance(sc.get("isppa"), dict) else {}
    level = isppa.get("level_wcm2")
    if need_level and _num(target) and not (_num(level) and level > 0):
        why = (" (a combined multi-transducer field)" if isppa.get("combination_mode") else "")
        issues.append(_issue("error", "transducer.target_isppa_wcm2",
                             f"The map's free-water ISPPA is not recorded{why}, so it cannot be "
                             "rescaled to a new target (load_acoustic_input.m:124-136)"))
    diffs = signature_diffs(sc.get("signature") or {}, model.effective())
    if diffs:
        shown = ", ".join(diffs[:6]) + (f" and {len(diffs) - 6} more" if len(diffs) > 6 else "")
        check = g("io.cache_signature_check")
        stop = " PRESTUS stops (io.cache_signature_check: error)." if check == "error" else ""
        issues.append(_issue("warning", key,
                             f"Settings differ from those the map was simulated with: {shown}.{stop}"))
    return issues


def validate_followup(parent_model: ConfigModel, parent_run, ftype: str, child: ConfigModel,
                      chain: list[str] | None = None) -> list[Issue]:
    """Problems with a follow-up before it is created. Errors block creation.
    ``chain``: the parent's current sequential chain (config paths), for 'sequential'."""
    issues: list[Issue] = []
    if parent_model.path is None or not Path(parent_model.path).is_file():
        issues.append(_issue("error", "", "Save the parent config first: the follow-up refers to it"))
    if parent_run is None:
        issues.append(_issue("error", "", "Choose the run to follow up"))
        return issues
    affix = sanitize_affix(child.get("io.output_affix") or "")
    if affix == parent_run.affix:
        issues.append(_issue("error", "io.output_affix",
                             "The follow-up needs its own output affix, or it overwrites the parent's outputs"))
    elif not affix:
        issues.append(_issue("error", "io.output_affix", "Set an output affix for the follow-up"))
    else:
        csv_name = f"{parent_run.tag}_{parent_run.medium}{affix}.csv"
        if (parent_run.base / csv_name).is_file():
            issues.append(_issue("warning", "io.output_affix",
                                 f"Outputs with this affix exist ({csv_name}); the follow-up replaces them"))

    if TYPES[ftype].uses_pressure:
        issues += pressure_issues(child, pressure_nifti(parent_run), need_level=ftype == "new_intensity")
        if ftype == "new_intensity" and not _num(_effective_tx(child).get("target_isppa_wcm2")):
            issues.append(_issue("error", "transducer.target_isppa_wcm2", "Set the new target ISPPA"))
    if ftype == "new_timing" and child.get("modules.run_heating_sims") != 1:
        issues.append(_issue("error", "modules.run_heating_sims", "A timing follow-up runs the thermal stage"))

    if ftype == "new_placement_same_head":
        if child.get("grid.mode") != "ras_plus":
            issues.append(_issue("info", "grid.mode",
                                 "With grid.mode other than 'ras_plus', the shared head model is smoothed "
                                 "before rotation; thin skull layers differ slightly (preproc_head.m:194-200)"))
        res = child.get("grid.resolution_mm")
        if _num(res):
            cache = parent_run.base / "cache" / f"{parent_run.tag}_{parent_run.medium}_headmodel_{res:g}mm.mat"
            if not cache.is_file():
                issues.append(_issue("info", "headmodel.reuse_across_placements",
                                     "No shared head model yet: this run builds it, later placements "
                                     "reuse it (head_model_load_or_build.m:46)"))

    if ftype == "pose_from_time_reversal":
        tx = _effective_tx(child)
        if child.n_transducers != 1:
            issues.append(_issue("error", TX, "Pose placement supports one transducer "
                                              "(time_reversal_pose_parameters.m:57-60)"))
        ann = tx.get("annular") or {}
        R = ann.get("curv_radius_mm") if isinstance(ann, dict) else None
        if tx.get("type") != "annular" or not _num(R):
            issues.append(_issue("error", "transducer.type",
                                 "Pose placement and the pose search need a curved annular transducer "
                                 "(time_reversal_pose_parameters.m:62-66)"))
        if child.get("time_reversal.apply_phases") is True:
            n_rings = len(ann.get("elem_od_mm") or []) if isinstance(ann, dict) else 0
            if n_rings == 1:
                issues.append(_issue("warning", "time_reversal.apply_phases",
                                     "A single element has no phases to apply (doc_time_reversal.md)"))
        elif is_empty(child.get("transducer.trans_pos")) or child.get("placement.manual.coordinate_space") != "ras":
            issues.append(_issue("error", "transducer.trans_pos", "Choose a pose from the pose table"))

    if ftype == "sequential":
        g = parent_model.get
        if g("modules.run_heating_sims") != 1:
            issues.append(_issue("error", "modules.run_heating_sims",
                                 "The parent run has no thermal stage: there is no end state to carry over "
                                 "(sequential_pipeline.m:80-86 stops with missingThermalEnd)"))
        pt = _effective_tx(parent_model).get("target_isppa_wcm2")
        multi = isinstance(pt, list) and len([x for x in pt if _num(x)]) > 1
        async_ = parent_model.n_transducers > 1 and g("simulation.transducer_coupling") == "async"
        if g("simulation.uncertainty") is True and (multi or async_):
            issues.append(_issue("error", "simulation.uncertainty",
                                 "Uncertainty runs with a multi-ISPPA or async pipeline cannot chain "
                                 "sequential runs (uncertainty_pipeline.m:275-281)"))
        elif multi or async_:
            issues.append(_issue("warning", "", "The parent runs a multi-ISPPA or async pipeline; "
                                                "only a plain run hands its end state on (prestus_pipeline.m:563)"))
        if chain:
            issues.append(_issue("info", "", f"Added as step {len(chain) + 2} after the "
                                             f"{len(chain)} follow-up(s) already chained to this config"))
    issues += adopted_issues(child)
    order = {"error": 0, "warning": 1, "info": 2}
    return sorted(issues, key=lambda i: order.get(i.level, 3))


ADOPTED_KEYS = ("io.adopted_heatmap", "io.adopted_cem43", "io.adopted_cem43_iso")


def adopted_issues(model: ConfigModel) -> list[Issue]:
    """io.adopted_* set but missing: thermal_simulation.m:163,292,306 only adopt a file that
    exists and otherwise start from the baseline without a warning."""
    issues = []
    for key in ADOPTED_KEYS:
        value = model.get(key)
        if isinstance(value, str) and value:
            p = resolve(model, value)
            if not p.is_file():
                issues.append(_issue("error", key, f"File not found: {p}. PRESTUS would silently start "
                                                   "from the baseline temperature (thermal_simulation.m:163)"))
    return issues


def resolve(model: ConfigModel, value: str) -> Path:
    """A path from the config; relative ones against the config's folder (where the Nexus
    launcher starts MATLAB, launcher.run_command -sd)."""
    p = Path(value).expanduser()
    if not p.is_absolute() and model.path is not None:
        p = Path(model.path).parent / p
    return p


def lineage_issues(model: ConfigModel) -> list[Issue]:
    lin = getattr(model, "lineage", None)
    if not isinstance(lin, dict):
        return []
    parent = lin.get("parent_config")
    if parent and not Path(parent).is_file():
        return [_issue("warning", "", f"This follow-up's parent config no longer exists: {parent}")]
    return []


# ── saving ─────────────────────────────────────────────────────────────────────

def save_followup(child: ConfigModel, path: Path, header: str = "") -> Path:
    """Write the child config atomically (fsio), with the lineage header."""
    from .fsio import write_text_atomic
    path = Path(path)
    write_text_atomic(path, child.to_yaml(header))
    child.path = path
    return path
