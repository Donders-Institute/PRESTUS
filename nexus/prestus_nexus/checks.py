"""Validation of a study config before it is handed to MATLAB.

Checks mirror what PRESTUS itself enforces (with the MATLAB source noted),
so that problems show up while editing instead of at submission.
MATLAB's load_parameters stays the final authority.
"""

from __future__ import annotations

import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .model import TX, ConfigModel, Issue, is_empty
from .paths import m2m_dir, resolve_pattern
from .schema import Node


def _is_num(v: Any) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def check_value(node: Node, value: Any) -> str | None:
    """Problem with a single value for its schema kind, or None."""
    kind = node.kind
    if kind in ("yaml",) or node.derived:
        return None
    if kind == "int":
        if is_empty(value):
            return None
        if not _is_num(value) or float(value) != int(value):
            return "must be an integer"
        if node.minimum is not None and value < node.minimum:
            return f"must be ≥ {node.minimum}"
        if node.maximum is not None and value > node.maximum:
            return f"must be ≤ {node.maximum}"
        return None
    if kind in ("number", "float"):
        return None if is_empty(value) or _is_num(value) else "must be a number"
    if kind == "numbers":
        if is_empty(value) or _is_num(value):
            return None
        if isinstance(value, list) and all(_is_num(v) or isinstance(v, list) for v in value):
            return None
        return "must be a number or a list of numbers"
    if kind == "vec3":
        if is_empty(value):
            return None
        if isinstance(value, list) and len(value) == 3 and all(_is_num(v) for v in value):
            return None
        return "must be three numbers [x, y, z]"
    if kind == "flag":
        return None if value in (0, 1) else "must be 0 or 1"
    if kind == "bool":
        return None if isinstance(value, bool) or is_empty(value) else "must be true or false"
    if kind == "tristate":
        return None if value is None or isinstance(value, bool) else "must be true, false or ~ (auto)"
    if kind == "strlist":
        if is_empty(value) or isinstance(value, str):
            return None
        if isinstance(value, list) and all(isinstance(v, str) for v in value):
            return None
        return "must be a list of strings"
    if kind in ("str", "dir", "file", "path"):
        if is_empty(value):
            return None
        if not isinstance(value, str):
            return "must be text"
        if node.choices and value not in node.choices:
            return f"'{value}' is not one of: {', '.join(map(str, node.choices))}"
    return None


def value_issues(model: ConfigModel) -> list[Issue]:
    """Type/choice problems in every overridden value."""
    issues: list[Issue] = []

    def walk(data: dict, prefix: str, ti: int | None) -> None:
        for key, value in data.items():
            path = f"{prefix}.{key}" if prefix else key
            node = model.schema.get(path)
            if node is None:
                continue
            if node.is_group:
                if isinstance(value, dict):
                    walk(value, path, ti)
                continue
            problem = check_value(node, value)
            if problem:
                level = "warning" if node.choices and node.path == "io.overwrite_files" else "error"
                issues.append(Issue(level, path, problem, ti))
            if node.derived:
                issues.append(Issue("info", path, "Derived by PRESTUS; the value set here is overwritten", ti))

    walk(model.overrides, "", None)
    for ti, t in enumerate(model.tx):
        walk(t, TX, ti)
    return issues


def subject_issues(model: ConfigModel) -> list[Issue]:
    """Subject, path and environment checks (Setup tab inputs)."""
    issues: list[Issue] = []
    g = model.get
    sid = g("subject_id")
    layered = g("simulation.medium") == "layered"

    if not _is_num(sid) or float(sid) != int(sid) or sid < 0:
        issues.append(Issue("error", "subject_id", "Subject ID must be a non-negative integer"))
        sid = None
    elif (sid := int(sid)) > 999:
        issues.append(Issue("warning", "subject_id", "IDs above 999 do not fit the sub-%03d naming"))

    sim = g("path.sim")
    if is_empty(sim):
        issues.append(Issue("error", "path.sim", "Output folder is required (get_output_dir.m)"))
    elif not os.path.isabs(sim):
        issues.append(Issue("warning", "path.sim",
                            "Relative path: PRESTUS resolves it against MATLAB's working folder"))
    elif not Path(sim).is_dir():
        parent = Path(sim).parent
        level = "info" if parent.is_dir() else "warning"
        issues.append(Issue(level, "path.sim", "Folder does not exist yet; PRESTUS will create it"))

    anat = g("path.anat")
    if layered:
        if is_empty(anat):
            issues.append(Issue("error", "path.anat", "Anatomy folder is required for layered simulations"))
        elif not Path(anat).is_dir():
            issues.append(Issue("error", "path.anat", "Folder not found"))
        elif sid is not None:
            hits = resolve_pattern(anat, g("path.t1_pattern") or "", sid)
            if not hits:
                issues.append(Issue("error", "path.t1_pattern", "No T1 image matches this pattern"))
            elif len(hits) > 1:
                issues.append(Issue("warning", "path.t1_pattern",
                                    f"{len(hits)} files match; PRESTUS uses {hits[0].name}"))
            else:
                issues.append(Issue("info", "path.t1_pattern", f"T1: {hits[0]}"))
            t2 = g("path.t2_pattern")
            if not is_empty(t2) and not resolve_pattern(anat, t2, sid):
                issues.append(Issue("warning", "path.t2_pattern", "No T2 image matches this pattern"))
            ct = g("path.ct_pattern")
            if g("pct.enabled") in (1, True) and not is_empty(ct):
                ct_hits = resolve_pattern(anat, ct, sid)
                if not ct_hits:
                    issues.append(Issue("error", "path.ct_pattern", "No CT image matches this pattern"))
                else:
                    issues.append(Issue("info", "path.ct_pattern",
                                        f"Measured CT: {ct_hits[0]} (must be co-registered to the T1)"))

    seg = g("path.seg")
    seg_root = anat if is_empty(seg) else seg      # load_parameters: seg defaults to anat
    mode = g("placement.mode")
    needs_m2m = layered or mode in SEGMENTATION_MODES or \
        g("placement.manual.coordinate_space") == "mni"
    if needs_m2m and sid is not None and not is_empty(seg_root):
        m2m = m2m_dir(seg_root, sid)
        simnibs = g("startup.simnibs_bin_path")
        if m2m.is_dir():
            issues.append(Issue("info", "path.seg", f"Segmentation found: {m2m}"))
            mesh = m2m / f"sub-{sid:03d}.msh"
            if mode == "plantus" and not mesh.is_file():
                issues.append(Issue("error", "path.seg",
                                    f"PlanTUS needs the SimNIBS mesh {mesh.name}, not found in {m2m.name}"))
        elif is_empty(simnibs):
            issues.append(Issue("error", "startup.simnibs_bin_path",
                                f"{m2m.name} not found and no SimNIBS path set to run charm"))
        else:
            issues.append(Issue("info", "path.seg", f"{m2m.name} not found; charm will run"))
        if not is_empty(simnibs) and not Path(simnibs).is_dir():
            issues.append(Issue("error", "startup.simnibs_bin_path", "Folder not found"))
    return issues


# ── Placement modes: availability for the configured subject ─────────────

PLACEMENT_MODES = ("manual", "localite", "heuristic", "mni", "plantus", "tus_entry")
# Modes that read the SimNIBS m2m folder (segmentation and/or MNI transforms)
SEGMENTATION_MODES = ("mni", "heuristic", "plantus", "tus_entry")
_RANK = {"ready": 0, "pending": 1, "unavailable": 2}


@dataclass(frozen=True)
class ModeStatus:
    """Whether a placement mode can run for the configured subject.

    state  'ready'        all inputs found
           'pending'      selectable: something runs first (SimNIBS charm) or a setting
                          of the mode itself is still to be filled in (see reason)
           'unavailable'  an input is missing that PRESTUS cannot produce (no T1; no
                          segmentation and no SimNIBS to make one; no recorded poses)
    Settings that are only shown once the mode is selected (tool paths, targets)
    never make a mode unavailable, or it could not be selected to fill them in.
    """
    state: str
    reason: str = ""


def _worst(*parts: ModeStatus) -> ModeStatus:
    state = max((p.state for p in parts), key=_RANK.__getitem__)
    reasons = [p.reason for p in parts if p.state != "ready" and p.reason]
    if state == "ready":
        reasons = [p.reason for p in parts if p.reason]
    return ModeStatus(state, "; ".join(reasons))


def _subject_id(model: ConfigModel) -> int | None:
    sid = model.get("subject_id")
    return int(sid) if _is_num(sid) and float(sid) == int(sid) and sid >= 0 else None


def _t1_status(model: ConfigModel) -> ModeStatus:
    g = model.get
    sid, anat = _subject_id(model), g("path.anat")
    if sid is None or is_empty(anat):
        return ModeStatus("pending", "set subject_id and path.anat to find the T1")
    if not resolve_pattern(anat, g("path.t1_pattern") or "", sid):
        return ModeStatus("unavailable", "no T1 matches path.anat / path.t1_pattern")
    return ModeStatus("ready")


def _segmentation_status(model: ConfigModel, mesh: bool = False) -> ModeStatus:
    """The SimNIBS m2m folder: tissue labels, MNI warps and (mesh=True) the head mesh."""
    g = model.get
    sid = _subject_id(model)
    seg_root = g("path.seg") if not is_empty(g("path.seg")) else g("path.anat")
    if sid is None or is_empty(seg_root):
        return ModeStatus("pending", "set subject_id and path.seg (or path.anat) to find the segmentation")
    m2m = m2m_dir(seg_root, sid)
    need = [m2m / "final_tissues.nii.gz", m2m / "toMNI" / "MNI2Conform_nonl.nii.gz"]
    if mesh:
        need.append(m2m / f"sub-{sid:03d}.msh")
    missing = [p.relative_to(m2m).as_posix() for p in need if not p.is_file()]
    if not missing:
        return ModeStatus("ready")
    simnibs = g("startup.simnibs_bin_path")
    what = f"{m2m.name}: {', '.join(missing)} not found"
    if not is_empty(simnibs) and Path(simnibs).is_dir():
        return ModeStatus("pending", f"{what}; SimNIBS charm runs first")
    return ModeStatus("unavailable", f"{what}, and no startup.simnibs_bin_path to run charm")


def _localite_status(model: ConfigModel) -> ModeStatus:
    from .localite import post_root, raw_root   # avoids an import cycle at module load
    g = model.get
    f = g("placement.localite.file")
    if not is_empty(f):
        return ModeStatus("ready") if Path(str(f)).is_file() else \
            ModeStatus("pending", "placement.localite.file not found")
    sid = _subject_id(model)
    if sid is None:
        return ModeStatus("pending", "set subject_id to find the recorded poses")
    sub = f"sub-{sid:03d}"
    if (post := post_root(model)) and (Path(post) / sub).is_dir():
        return ModeStatus("ready", f"ingested poses in {Path(post) / sub}")
    if (raw := raw_root(model)) and (Path(raw) / sub).is_dir():
        return ModeStatus("ready", f"raw exports in {Path(raw) / sub}; ingested on first use")
    return ModeStatus("unavailable", f"no recorded poses for {sub} (placement.localite.file, "
                                     "path.localite_post or path.localite_raw)")


def _tool_status(model: ConfigModel, mode: str) -> ModeStatus:
    """External tools: only 'pending' (to be set or fixed in the mode's own settings)."""
    g = model.get
    todo = []
    if mode == "plantus":
        if is_empty(g("placement.plantus.script_path")) or \
                (is_empty(g("placement.plantus.env_path")) and is_empty(g("startup.simnibs_bin_path"))):
            todo.append("set the PlanTUS paths (placement.plantus.script_path; startup.simnibs_bin_path "
                        "or placement.plantus.env_path)")
        if is_empty(g("placement.plantus.focal_distance_list")):
            todo.append("set placement.plantus.focal_distance_list / flhm_list")
    elif mode == "tus_entry":
        if is_empty(g("placement.tus_entry.script_path")) or is_empty(g("placement.tus_entry.env_path")):
            todo.append("set the TUS_entry paths (placement.tus_entry.script_path, env_path)")
    return ModeStatus("pending", "; ".join(todo)) if todo else ModeStatus("ready")


# Settings placement_mode_status reads. Its answer is reused for _STATUS_TTL_S while
# these are unchanged: it looks for files (T1, m2m, poses), ~1 s on a network drive, and
# the Issues panel and the Placement card both ask after every edit.
_STATUS_KEYS = ("subject_id", "path.anat", "path.t1_pattern", "path.seg", "startup.simnibs_bin_path",
                "placement.manual.coordinate_space", "placement.localite.file", "path.localite_post",
                "path.localite_raw", "path.localite", "placement.plantus.script_path",
                "placement.plantus.env_path", "placement.plantus.focal_distance_list",
                "placement.tus_entry.script_path", "placement.tus_entry.env_path")
_STATUS_TTL_S = 5.0
_status_cache: dict[tuple, tuple[float, dict]] = {}


def clear_status_cache() -> None:
    _status_cache.clear()


def placement_mode_status(model: ConfigModel) -> dict[str, ModeStatus]:
    """Availability of every placement mode for the configured subject
    (preproc_transducer_placement.m and the position_transducer_* functions)."""
    key = tuple(repr(model.get(k)) for k in _STATUS_KEYS)
    hit = _status_cache.get(key)
    if hit is not None and time.monotonic() - hit[0] < _STATUS_TTL_S:
        return hit[1]
    status = _placement_mode_status(model)
    _status_cache.clear()
    _status_cache[key] = (time.monotonic(), status)
    return status


def _placement_mode_status(model: ConfigModel) -> dict[str, ModeStatus]:
    t1 = _t1_status(model)
    seg = _segmentation_status(model)
    manual = t1
    if model.get("placement.manual.coordinate_space") == "mni":
        manual = _worst(t1, seg)
    return {
        "manual":    manual,
        "localite":  _worst(t1, _localite_status(model)),
        "heuristic": _worst(t1, seg),
        "mni":       _worst(t1, seg),
        "plantus":   _worst(t1, _segmentation_status(model, mesh=True), _tool_status(model, "plantus")),
        "tus_entry": _worst(t1, seg, _tool_status(model, "tus_entry")),
    }


def placement_issues(model: ConfigModel) -> list[Issue]:
    """The selected placement mode: availability and its tool settings
    (position_transducer_plantus.m, position_transducer_tus_entry.m)."""
    issues: list[Issue] = []
    g = model.get
    mode = g("placement.mode") or "manual"
    st = placement_mode_status(model).get(mode)
    if st is not None and st.state == "unavailable":
        issues.append(Issue("error", "placement.mode", f"'{mode}' cannot run: {st.reason}"))

    if mode == "plantus":
        script, env = g("placement.plantus.script_path"), g("placement.plantus.env_path")
        if is_empty(script):
            issues.append(Issue("error", "placement.plantus.script_path", "Required for PlanTUS"))
        elif not (Path(script) / "code" / "PlanTUS.py").is_file():
            issues.append(Issue("error", "placement.plantus.script_path", "code/PlanTUS.py not found here"))
        if is_empty(env):
            issues.append(Issue("error", "placement.plantus.env_path",
                                "Required: the SimNIBS (or PRESTUS) Python environment"))
        elif not Path(env).is_dir():
            issues.append(Issue("error", "placement.plantus.env_path", "Folder not found"))
        wb = g("placement.plantus.connectome_wb_path")
        if not is_empty(wb) and not Path(wb).is_dir():
            issues.append(Issue("warning", "placement.plantus.connectome_wb_path", "Folder not found"))
        fd, flhm = g("placement.plantus.focal_distance_list"), g("placement.plantus.flhm_list")
        if is_empty(fd):
            issues.append(Issue("error", "placement.plantus.focal_distance_list",
                                "Required: focal depth(s) of the transducer [mm]"))
        elif not is_empty(flhm):
            n_fd = len(fd) if isinstance(fd, list) else 1
            n_fl = len(flhm) if isinstance(flhm, list) else 1
            if n_fd != n_fl:
                issues.append(Issue("error", "placement.plantus.flhm_list",
                                    f"{n_fl} value(s) for {n_fd} focal distance(s)"))
        lo, hi = g("placement.plantus.min_distance_mm"), g("placement.plantus.max_distance_mm")
        if _is_num(lo) and _is_num(hi) and lo >= hi:
            issues.append(Issue("error", "placement.plantus.max_distance_mm",
                                "Must be larger than min_distance_mm"))
        if is_empty(g("placement.plantus.mni_target_mm")) and \
                all(is_empty(t.get("focus_pos")) for t in model.effective_transducers()):
            issues.append(Issue("error", "placement.plantus.mni_target_mm",
                                "Set the target: mni_target_mm or transducer.focus_pos"))

    if mode == "tus_entry":
        script, env = g("placement.tus_entry.script_path"), g("placement.tus_entry.env_path")
        if is_empty(script):
            issues.append(Issue("error", "placement.tus_entry.script_path", "Required for TUS_entry"))
        elif not Path(script).is_file():
            issues.append(Issue("error", "placement.tus_entry.script_path", "File not found"))
        if is_empty(env):
            issues.append(Issue("error", "placement.tus_entry.env_path",
                                "Required: conda environment with Rscript"))
        elif not Path(env).is_dir():
            issues.append(Issue("error", "placement.tus_entry.env_path", "Folder not found"))
        elif not any((Path(env) / p).is_file() for p in ("bin/Rscript", "Rscript", "bin/Rscript.exe")):
            issues.append(Issue("warning", "placement.tus_entry.env_path", "No Rscript found in this environment"))
        for key in ("scalp_mask", "target_mask", "exclu_mask"):
            v = g(f"placement.tus_entry.{key}")
            if not is_empty(v) and not Path(str(v)).is_file():
                issues.append(Issue("error", f"placement.tus_entry.{key}", "File not found"))
        if is_empty(g("placement.tus_entry.target_mask")) and is_empty(g("placement.tus_entry.mni_target_mm")) \
                and all(is_empty(t.get("focus_pos")) for t in model.effective_transducers()):
            issues.append(Issue("error", "placement.tus_entry.mni_target_mm",
                                "Set the target: target_mask, mni_target_mm or transducer.focus_pos"))
        lo, hi = g("placement.tus_entry.minimal_distance_cm"), g("placement.tus_entry.maximal_distance_cm")
        if _is_num(lo) and _is_num(hi) and lo >= hi:
            issues.append(Issue("error", "placement.tus_entry.maximal_distance_cm",
                                "Must be larger than minimal_distance_cm"))
    return issues


def rule_issues(model: ConfigModel) -> list[Issue]:
    """Cross-field rules PRESTUS enforces at load or run time."""
    issues: list[Issue] = []
    g = model.get

    # load_parameters.m: backend combinations
    ab, tb = g("acoustic_backend"), g("thermal_backend")
    if tb == "babelbrain" and ab != "babelbrain":
        issues.append(Issue("error", "thermal_backend",
                            "thermal_backend 'babelbrain' requires acoustic_backend 'babelbrain'"))

    # babelbrain_check_geometry.m: BabelBrain simulates one beam, along the grid z axis
    if "babelbrain" in (ab, tb):
        key = "acoustic_backend" if ab == "babelbrain" else "thermal_backend"
        if len(model.effective_transducers()) > 1:
            issues.append(Issue("error", key,
                                "The BabelBrain backend supports one transducer per run; "
                                "use kwave or jwave for multi-transducer runs"))
        if g("grid.mode") == "ras_plus":
            issues.append(Issue("error", "grid.mode",
                                "grid.mode 'ras_plus' leaves the beam off the grid z axis, which the "
                                "BabelBrain backend cannot simulate; use 'focal_axis'"))

    # coherent_drive_calibration.m / multi_isppa_pipeline.m: coherent multi-transducer
    # runs set each drive from one target per transducer; a sweep needs one run per point
    txs_all = model.effective_transducers()
    if len(txs_all) > 1 and g("simulation.transducer_coupling") != "async":
        for ti, t in enumerate(txs_all):
            v = t.get("target_isppa_wcm2")
            if isinstance(v, (list, tuple)) and len([x for x in v if _is_num(x)]) > 1:
                issues.append(Issue("error", "transducer.target_isppa_wcm2",
                                    "Coherent multi-transducer runs take one target per transducer "
                                    "(each sweep point needs its own simulation); use async coupling "
                                    "for a sweep", ti))
                break

    # prestus_pipeline.m: mm positions need the placement stage
    cs = g("placement.manual.coordinate_space")
    if g("placement.mode") == "manual" and cs in ("ras", "mni") and \
            g("modules.run_transducer_placement") != 1:
        issues.append(Issue("error", "modules.run_transducer_placement",
                            f"coordinate_space '{cs}' needs the placement stage to convert positions"))

    # load_transducer_parameters.m: one frequency for all transducers
    txs = model.effective_transducers()
    freqs = {t.get("freq_hz") for t in txs if _is_num(t.get("freq_hz"))}
    if len(freqs) > 1:
        issues.append(Issue("error", "transducer.freq_hz", "All transducers must use the same frequency"))

    # load_parameters.m (restore_struct_arrays): identical fields in every element
    first_keys = set(txs[0].keys())
    for ti, t in enumerate(model.tx[1:], start=1):
        extra = sorted(set(t) - first_keys)
        if extra:
            issues.append(Issue("error", f"{TX}.{extra[0]}",
                                f"Set {', '.join(extra)} on transducer 1 as well; "
                                "PRESTUS requires the same fields for every transducer", ti))

    for ti, t in enumerate(txs):
        if is_empty(t.get("serial")) and is_empty(t.get("type")):
            issues.append(Issue("warning", f"{TX}.type",
                                "Set a transducer type (or an equipment serial)", ti))
        if is_empty(t.get("serial")) and is_empty(t.get("freq_hz")):
            issues.append(Issue("warning", f"{TX}.freq_hz", "Frequency not set", ti))
        if g("placement.mode") == "manual":
            for key in ("trans_pos", "focus_pos"):
                if is_empty(t.get(key)):
                    issues.append(Issue("warning", f"{TX}.{key}",
                                        "Needed for manual placement", ti))

    if g("placement.mode") == "plantus":
        names = ["skin_target_distances", "skin_target_angles", "skin_target_intersections",
                 "skin_skull_angles", "skull_thickness"]
        weights = [g(f"placement.plantus.weight_{n}") for n in names]
        if all(_is_num(w) for w in weights) and abs(sum(weights) - 1) > 1e-6:
            issues.append(Issue("warning", "placement.plantus.weight_skull_thickness",
                                f"PlanTUS weights sum to {sum(weights):.3g}, expected 1"))
    return issues


def library_issues(model: ConfigModel) -> list[Issue]:
    """Equipment/calibration library rules (resolve_transducer_from_serial.m)."""
    issues: list[Issue] = []
    if model.library is None:
        return issues
    for ti, t in enumerate(model.effective_transducers()):
        serial = t.get("serial")
        if is_empty(serial):
            continue
        dev = model.library.get(serial)
        if dev is None:
            issues.append(Issue("error", f"{TX}.serial",
                                f"'{serial}' is not in config/equipment", ti))
            continue
        ds = (t.get("combo") or {}).get("ds_serial") if isinstance(t.get("combo"), dict) else None
        if not is_empty(ds) and ds not in dev.combos:
            issues.append(Issue("warning", f"{TX}.combo.ds_serial",
                                f"'{ds}' is not a listed driving system for {serial} "
                                f"({', '.join(dev.combos) or 'none listed'})", ti))
        ep = t.get("focal_distance_ep")
        if _is_num(ep):
            if _is_num(dev.min_foc) and ep < dev.min_foc:
                issues.append(Issue("warning", f"{TX}.focal_distance_ep",
                                    f"{ep:g} mm is below the device's min_foc {dev.min_foc:g} mm", ti))
            if _is_num(dev.max_foc) and ep > dev.max_foc:
                issues.append(Issue("warning", f"{TX}.focal_distance_ep",
                                    f"{ep:g} mm exceeds the device's max_foc {dev.max_foc:g} mm", ti))
        if t.get("type") != "annular":
            continue
        ann = t.get("annular") or {}
        if not is_empty(ann.get("elem_phase_deg")) and not is_empty(ann.get("elem_p_amp")):
            continue            # fully specified inline: no library lookup
        isppa = t.get("target_isppa_wcm2")
        if not _is_num(ep) or is_empty(isppa):
            issues.append(Issue("error", f"{TX}.focal_distance_ep",
                                "Phases/amplitude come from the calibration library, which needs "
                                "focal_distance_ep and target_isppa_wcm2", ti))
            continue
        cal = model.library.calibration(serial, None if is_empty(ds) else ds,
                                        t.get("library_path") or None)
        if not cal.found:
            issues.append(Issue("error", f"{TX}.combo.ds_serial",
                                f"No calibration entry ({cal.key}.yaml or {serial}.yaml); "
                                "set phases and amplitude inline or build the library", ti))
        elif cal.fallback:
            issues.append(Issue("warning", f"{TX}.combo.ds_serial",
                                f"No entry for {cal.key}; PRESTUS falls back to {cal.path.name}", ti))
        elif cal.depths_mm and not min(cal.depths_mm) <= ep <= max(cal.depths_mm):
            issues.append(Issue("warning", f"{TX}.focal_distance_ep",
                                f"{ep:g} mm is outside the calibrated depths "
                                f"{min(cal.depths_mm):g}–{max(cal.depths_mm):g} mm", ti))
    return issues


def isppa_issues(model: ConfigModel) -> list[Issue]:
    """Whether transducer.target_isppa_wcm2 can be met (needs_water_baseline.m,
    isppa_scale_factor.m): the free-water baseline must run unless the library
    calibration supplies a matching ISPPA."""
    issues: list[Issue] = []
    mode = model.get("modules.run_water_baseline")
    off = mode in (0, False, "0")
    for ti, t in enumerate(model.effective_transducers()):
        target = t.get("target_isppa_wcm2")
        if not _is_num(target):
            continue                       # empty (drive level) or a multi-ISPPA vector
        ann = t.get("annular") or {}
        inline_drive = isinstance(ann, dict) and not is_empty(ann.get("elem_p_amp"))
        library_drive = not is_empty(t.get("serial")) and not inline_drive
        if library_drive:
            continue                       # calibration ISPPA from the library
        if off:
            level = "warning" if not is_empty(t.get("serial")) else "error"
            issues.append(Issue(level, "modules.run_water_baseline",
                                f"target_isppa_wcm2 = {target:g} W/cm² needs a free-water baseline "
                                "(this transducer has no library calibration), but "
                                "run_water_baseline is 0: PRESTUS stops with prestus:isppa:noBaseline. "
                                "Use 'auto', or set target_isppa_wcm2: [] to report the drive level", ti))
        elif mode != 1:
            issues.append(Issue("info", f"{TX}.target_isppa_wcm2",
                                f"An extra free-water simulation runs to scale the field to "
                                f"{target:g} W/cm² (no library calibration; run_water_baseline: auto). "
                                "Set target_isppa_wcm2: [] to report the elem_p_amp drive level", ti))
    return issues


def defaults_issues(model: ConfigModel) -> list[Issue]:
    """A project copy of config_default.yaml that predates toolbox keys (load_parameters.m)."""
    schema = model.schema
    if schema.project_path is None or not schema.stale_keys:
        return []
    shown = ", ".join(schema.stale_keys[:4]) + (", …" if len(schema.stale_keys) > 4 else "")
    return [Issue("warning", "",
                  f"The project's config_default.yaml lacks {len(schema.stale_keys)} key(s) added to "
                  f"PRESTUS since it was copied ({shown}); they use the toolbox defaults. "
                  "Tools ▸ Update project defaults… refreshes the copy")]


def uncertainty_issues(model: ConfigModel) -> list[Issue]:
    """uncertainty_pipeline.m errors when a variant's medium file is missing (configNotFound)."""
    if model.get("simulation.uncertainty") is not True:
        return []
    from .uncertainty import VARIANT_SHORT, property_sets
    sets = property_sets(model, model.schema.toolbox_path.parent.parent)
    return [Issue("error", f"simulation.uncertainty_{variant}", f"{VARIANT_SHORT[variant]} variant: {why}")
            for variant, why in sets.errors.items()]


def all_issues(model: ConfigModel) -> list[Issue]:
    order = {"error": 0, "warning": 1, "info": 2}
    issues = (value_issues(model) + model.unknown_keys()
              + subject_issues(model) + rule_issues(model) + placement_issues(model)
              + library_issues(model) + isppa_issues(model) + defaults_issues(model)
              + uncertainty_issues(model))
    return sorted(issues, key=lambda i: order.get(i.level, 3))
