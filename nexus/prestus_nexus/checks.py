"""Validation of a study config before it is handed to MATLAB.

Checks mirror what PRESTUS itself enforces (with the MATLAB source noted),
so that problems show up while editing instead of at submission.
MATLAB's load_parameters stays the final authority.
"""

from __future__ import annotations

import math

import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from . import simnibs as sn
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


class _KeyRecorder:
    """The model, recording which keys get() reads (subject_issues' cache key)."""

    def __init__(self, model: ConfigModel):
        self._model = model
        self.read: dict[tuple, str] = {}

    def get(self, path: str, ti: int = 0):
        value = self._model.get(path, ti)
        self.read[(path, ti)] = repr(value)
        return value

    def __getattr__(self, name):
        return getattr(self._model, name)


# subject_issues looks for folders and files (~1 s on a network drive) and runs after every
# edit. Its answer is reused for _STATUS_TTL_S while every setting it read is unchanged, so
# typing in an unrelated field does no file I/O.
_subject_cache: list[tuple[float, dict[tuple, str], tuple, list[Issue]]] = []


def _simnibs_context() -> tuple:
    """What the SimNIBS checks read besides the config: the SSH host and its scanned installs."""
    return sn.remote_host(), sn.remote_generation()


def subject_issues(model: ConfigModel) -> list[Issue]:
    """Subject, path and environment checks (Setup tab inputs)."""
    now = time.monotonic()
    ctx = _simnibs_context()
    for stamp, read, c, issues in _subject_cache:
        if now - stamp < _STATUS_TTL_S and c == ctx and \
                all(repr(model.get(k, ti)) == v for (k, ti), v in read.items()):
            return list(issues)
    rec = _KeyRecorder(model)
    issues = _subject_issues(rec)
    _subject_cache[:] = [(time.monotonic(), rec.read, ctx, issues)]
    return list(issues)


def _subject_issues(model: ConfigModel) -> list[Issue]:
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
        issues += _simnibs_issues(model, m2m)
    return issues


def _simnibs_issues(model: ConfigModel, m2m: Path) -> list[Issue]:
    """The SimNIBS environment (SimNIBS 4: charm), the charm settings, and an existing
    segmentation that another version or other settings made (segmentation_check_existing.m)."""
    g = model.get
    out: list[Issue] = []
    bin_path, settings = g("startup.simnibs_bin_path"), g("segmentation.settings_file")
    inst = None
    host = sn.remote_host()             # runs submitted over SSH: the install may be on the cluster
    if not is_empty(bin_path):
        if Path(bin_path).is_dir():
            inst = sn.install_info(bin_path)
            if inst is None:
                out.append(Issue("error", "startup.simnibs_bin_path", "No charm in this folder (SimNIBS bin folder)"))
        else:
            inst = sn.lookup(bin_path)
            if inst is None and host:
                out.append(Issue("warning", "startup.simnibs_bin_path",
                                 f"Not on this computer, nor among the SimNIBS installs found on {host}: "
                                 "scan the cluster on Medium › Segmentation"))
            elif inst is None:
                out.append(Issue("error", "startup.simnibs_bin_path", "Folder not found"))
        if inst is not None and inst.supported is False:
            out.append(Issue("error", "startup.simnibs_bin_path",
                             f"SimNIBS {inst.version}: PRESTUS segments with charm from SimNIBS 4 or later"))
        for key, value in sn.linked_mismatches(g, inst):
            out.append(Issue("warning", key,
                             f"Belongs to another SimNIBS install than startup.simnibs_bin_path ({inst.bin}); "
                             "leave it empty to follow the selected environment"))
        if g("segmentation.deface") in (1, True) and inst is not None and inst.pydeface is False:
            out.append(Issue("warning", "segmentation.deface",
                             f"pydeface not found in {inst.bin}: no defaced copy is made "
                             "(install it in the SimNIBS environment)"))
    if not is_empty(settings):
        sp = Path(settings)
        if not sp.is_file():
            out.append(Issue("warning" if host else "error", "segmentation.settings_file",
                             "Not found on this computer" + (f"; charm on {host} needs it there" if host else "")))
        elif sp.resolve().parent == m2m.resolve():
            out.append(Issue("error", "segmentation.settings_file",
                             "charm deletes the m2m folder when it re-runs: keep the settings file outside it"))
    state = sn.segmentation_state(m2m)
    rerun = g("io.overwrite_simnibs") in (1, True)
    if state.exists and rerun:
        out.append(Issue("warning", "io.overwrite_simnibs",
                         f"charm deletes {m2m.name} (incl. toMNI/T1_to_MNI_post-hoc) and segments again"))
    elif state.exists:
        for d in sn.differences(state, inst, None if is_empty(settings) else str(settings)):
            out.append(Issue("warning", "path.seg",
                             f"{m2m.name} {d}; it is reused as is (Re-run SimNIBS segmentation to replace it)"))
    return out


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
    _subject_cache.clear()


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
            # empty: the Python of the SimNIBS install in startup.simnibs_bin_path
            inst = sn.lookup(g("startup.simnibs_bin_path"))
            if is_empty(g("startup.simnibs_bin_path")):
                issues.append(Issue("error", "placement.plantus.env_path",
                                    "Required without startup.simnibs_bin_path: the SimNIBS Python environment"))
            elif inst is not None and not inst.python:
                issues.append(Issue("warning", "placement.plantus.env_path",
                                    f"No Python found for the SimNIBS install {inst.bin}: set the environment here"))
        elif not Path(env).is_dir() and sn.remote_host() is None:
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
    model._sync_library()                  # equipment entries in transducer.library_path
    for ti, t in enumerate(model.effective_transducers()):
        serial = t.get("serial")
        if is_empty(serial):
            continue
        dev = model.library.get(serial)
        if dev is None:
            fallback = model.library.model_for(serial)
            if fallback is None:
                issues.append(Issue("error", f"{TX}.serial",
                                    f"'{serial}' is not in config/equipment or the study's library_path", ti))
                continue
            issues.append(Issue("warning", f"{TX}.serial",
                                f"'{serial}' is not in config/equipment; PRESTUS uses its model "
                                f"'{fallback}' (model geometry, geometric calibration)", ti))
            serial, dev = fallback, model.library.get(fallback)
        if dev.deprecated or dev.excluded:
            alt = f"; use '{dev.superseded_by}' for new studies" if dev.superseded_by else ""
            issues.append(Issue("warning", f"{TX}.serial",
                                f"'{serial}' is {dev.status} (kept for earlier studies){alt}", ti))
        ds =(t.get("combo") or {}).get("ds_serial") if isinstance(t.get("combo"), dict) else None
        is_model = dev.model == dev.serial and not dev.combos   # no driving systems to list
        if not is_empty(ds) and ds not in dev.combos and not is_model:
            issues.append(Issue("warning", f"{TX}.combo.ds_serial",
                                f"'{ds}' is not a listed driving system for {serial} "
                                f"({', '.join(dev.combos) or 'none listed'})", ti))
        ep = t.get("focal_distance_ep")
        if _is_num(ep):
            if _is_num(dev.min_foc) and ep < dev.min_foc:
                issues.append(Issue("warning", f"{TX}.focal_distance_ep",
                                    f"{ep:g} mm is below the hardware steering range (min_foc {dev.min_foc:g} mm)", ti))
            if _is_num(dev.max_foc) and ep > dev.max_foc:
                issues.append(Issue("warning", f"{TX}.focal_distance_ep",
                                    f"{ep:g} mm is beyond the hardware steering range (max_foc {dev.max_foc:g} mm)", ti))
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


STAIRCASE_MM = 0.75     # coarser grids: staircased skull and tissue boundaries in head media


def grid_issues(model: ConfigModel) -> list[Issue]:
    """Grid resolution against the wavelength and the head geometry.

    PPW as source_sensor_setup.m checks it: the slowest sound speed of the modelled tissues
    over the fundamental frequency, at least grid.min_ppw (default 6) points per wavelength;
    below that PRESTUS warns, and errors for jwave (its absorption diverges). For head and
    phantom media, a grid coarser than STAIRCASE_MM staircases the tissue transitions."""
    issues: list[Issue] = []
    g = model.get
    dx = g("grid.resolution_mm")
    if not _is_num(dx) or dx <= 0:
        return issues
    freqs = [t.get("freq_hz") for t in model.effective_transducers()]
    freqs = [float(f[0] if isinstance(f, (list, tuple)) and f else f) for f in freqs
             if _is_num(f) or (isinstance(f, (list, tuple)) and f and _is_num(f[0]))]
    medium = g("simulation.medium") or "layered"
    props = g("medium_properties") or {}
    tissues = ["water"] if medium == "water" else list(props)
    speeds = [(props.get(t, {}) or {}).get("sound_speed") for t in tissues if isinstance(props.get(t), dict)]
    speeds = [(float(c), t) for c, t in zip(speeds, tissues) if _is_num(c) and c > 0]
    min_ppw = g("grid.min_ppw")
    min_ppw = float(min_ppw) if _is_num(min_ppw) and min_ppw > 0 else 6.0
    if freqs and speeds:
        f = max(freqs)
        c, tissue = min(speeds)
        ppw = c / (f * dx / 1e3)
        if ppw < min_ppw:
            need = c / (f * min_ppw) * 1e3
            jwave = str(g("acoustic_backend") or "").startswith("jwave")
            issues.append(Issue(
                "error" if jwave else "warning", "grid.resolution_mm",
                f"{ppw:.1f} points per wavelength at {f / 1e3:.0f} kHz in {tissue} ({c:.0f} m/s), "
                f"below the minimum of {min_ppw:g}: set {need:.2f} mm or finer"
                + (" (jwave diverges below the minimum)" if jwave else "")))
    if medium in ("layered", "phantom") and dx > STAIRCASE_MM:
        issues.append(Issue(
            "warning", "grid.resolution_mm",
            f"Above {STAIRCASE_MM} mm, staircasing at the head tissue transitions (skull, scalp) "
            "may affect the results; use 0.75 mm or finer for results you report"))
    return issues


def uncertainty_issues(model: ConfigModel) -> list[Issue]:
    """uncertainty_pipeline.m errors when a variant's medium file is missing (configNotFound)."""
    if model.get("simulation.uncertainty") is not True:
        return []
    from .uncertainty import VARIANT_SHORT, property_sets
    sets = property_sets(model, model.schema.toolbox_path.parent.parent)
    return [Issue("error", f"simulation.uncertainty_{variant}", f"{VARIANT_SHORT[variant]} variant: {why}")
            for variant, why in sets.errors.items()]


def drive_issues(model: ConfigModel) -> list[Issue]:
    """A custom annular transducer (no library serial) is driven by its focal depth setting
    and target ISPPA (fill_annular_drive.m): without phases it needs focal_distance_ep (or
    _bowl), without an amplitude a target ISPPA. Flat rings keep needing phases."""
    issues: list[Issue] = []
    for ti, t in enumerate(model.effective_transducers()):
        if t.get("type") != "annular" or not is_empty(t.get("serial")):
            continue
        ann = t.get("annular") or {}
        if not isinstance(ann, dict):
            continue
        no_phase = is_empty(ann.get("elem_phase_deg")) and is_empty(ann.get("elem_phase_rad"))
        depth = not is_empty(t.get("focal_distance_ep")) or not is_empty(t.get("focal_distance_bowl"))
        flat = _is_num(ann.get("curv_radius_mm")) and not math.isfinite(ann["curv_radius_mm"])
        if no_phase and (flat or not depth):
            issues.append(Issue("error", f"{TX}.focal_distance_ep",
                                "No steering: set the focal depth setting"
                                + (" (flat rings need element phases)" if flat else "")
                                + ", or element phases from calibration (advanced)", ti))
        if is_empty(ann.get("elem_p_amp")) and is_empty(t.get("target_isppa_wcm2")):
            issues.append(Issue("error", f"{TX}.target_isppa_wcm2",
                                "No drive level: set the target free-water ISPPA, or element "
                                "pressure amplitudes from calibration (advanced)", ti))
    return issues


def _on(v: Any) -> bool:
    """A 0/1 flag or bool config value that is set."""
    return v is True or (_is_num(v) and v != 0)


def multitransducer_issues(model: ConfigModel) -> list[Issue]:
    """Multi-transducer and grid-orientation combinations PRESTUS refuses or reduces."""
    issues: list[Issue] = []
    g = model.get
    n = len(model.effective_transducers())
    axisym = _on(g("grid.axisymmetric"))
    is_async = n > 1 and g("simulation.transducer_coupling") == "async"     # is_async_mode.m
    medium = g("simulation.medium") or "layered"

    # grid_transducer_location.m:50-57: non-layered media keep transducer(1) only
    if n > 1 and not is_async and "layered" not in str(medium):
        issues.append(Issue("warning", "simulation.medium",
                            f"'{medium}' simulations model only transducer 1; the other "
                            f"{n - 1} transducer(s) are dropped. Use a layered medium for several "
                            "transducers"))
    # source_create.m:241-244: axisymmetric kWaveArray sources take one transducer
    if n > 1 and axisym and not is_async and g("grid.use_kWaveArray") not in (0, False):
        issues.append(Issue("error", "grid.axisymmetric",
                            "Axisymmetric simulations support one transducer; switch "
                            "grid.axisymmetric off for several transducers"))
    # combine_async_intensity.m:134-139: async fields must be 3-D
    if is_async and axisym:
        issues.append(Issue("error", "grid.axisymmetric",
                            "Async coupling combines 3-D fields; axisymmetric (2-D) results "
                            "cannot be combined. Switch grid.axisymmetric off"))
    # grid_axisymmetry.m:84-90: axisymmetry needs the beam on the grid z axis
    if g("grid.mode") == "ras_plus" and axisym:
        issues.append(Issue("error", "grid.mode",
                            "grid.mode 'ras_plus' is incompatible with grid.axisymmetric: the "
                            "beam must lie on the grid z axis ('focal_axis')"))
    return issues


def heating_issues(model: ConfigModel) -> list[Issue]:
    """Timing and thermal grid checks load_parameters.m runs when run_heating_sims is on
    (it calls thermal_parameters.m), and the thermal FOV against the acoustic FOV."""
    issues: list[Issue] = []
    g = model.get
    if _on(g("modules.run_heating_sims")):
        # thermal_parameters.m:76 / :81: ptd/pri and ptrd/ptri must be integers; with pri
        # or ptri 0 (the defaults) both are 0/0 = NaN and round_if_integer.m stops the run
        for key, msg in (("timing.pri", "Pulses per PT must be integer"),
                         ("timing.ptri", "PTRI reps in PTRD must be integer")):
            v = g(key)
            if not _is_num(v) or v <= 0:
                issues.append(Issue("error", key,
                                    f"Heating is on but {key} is not set (> 0): PRESTUS stops "
                                    f"with '{msg}'. Set the timing protocol (Timing tab)"))
        # thermal_parameters.m:64/:76/:81: pd, ptd or ptrd 0 leaves no sonication time
        zero = [k for k in ("timing.pd", "timing.ptd", "timing.ptrd") if _is_num(g(k)) and g(k) == 0]
        if zero:
            issues.append(Issue("warning", zero[0],
                                f"Heating is on but {', '.join(zero)} = 0: the thermal model has "
                                "no sonication time, so no heating"))
        # thermal_grid_setup.m:75-78: one thermal FOV extent per grid dimension
        fov = g("grid.thermal_fov_mm")
        if not is_empty(fov) and not _on(g("grid.axisymmetric")):
            if not (isinstance(fov, list) and len(fov) == 3 and all(_is_num(x) for x in fov)):
                issues.append(Issue("error", "grid.thermal_fov_mm",
                                    "Needs three extents [x, y, z] in mm (one per grid dimension)"))
    # resolve_acoustic_fov_params.m: the acoustic FOV is on only with both diameter and
    # length set; thermal_stage.m:89-91: thermal_use_acoustic_fov needs it
    diam, length = g("grid.acoustic_fov_diameter_mm"), g("grid.acoustic_fov_length_mm")
    fov_on = not is_empty(diam) and not is_empty(length)
    if _on(g("grid.thermal_use_acoustic_fov")) and not fov_on:
        issues.append(Issue("info", "grid.thermal_use_acoustic_fov",
                            "Has no effect without an acoustic FOV (set both "
                            "grid.acoustic_fov_diameter_mm and acoustic_fov_length_mm)"))
    # thermal_grid_setup.m:79-80: thermal_fov_mm is clamped to the grid the thermal stage
    # starts from, which is the acoustic FOV box when thermal_use_acoustic_fov is on
    fov = g("grid.thermal_fov_mm")
    if fov_on and _on(g("grid.thermal_use_acoustic_fov")) and len(model.effective_transducers()) == 1 \
            and isinstance(fov, list) and len(fov) == 3 and all(_is_num(x) for x in fov) \
            and _is_num(diam) and _is_num(length) and g("grid.mode") != "ras_plus":
        if max(fov[0], fov[1]) > diam or fov[2] > length:
            issues.append(Issue("info", "grid.thermal_fov_mm",
                                f"Larger than the acoustic FOV ({diam:g} mm across, {length:g} mm "
                                "long) that the thermal model runs in here; PRESTUS clamps it to "
                                "that box"))
    return issues


def _pichardo_shear_fit(f0: float) -> tuple[float, float]:
    """medium_elastic_properties.m pichardo_shear_fit: cS = a*rho + b, linear in frequency
    through the 270 and 836 kHz fits (extrapolated outside)."""
    t = (f0 - 270e3) / (836e3 - 270e3)
    return 0.458 + t * (0.386 - 0.458), 545.52 + t * (815.51 - 545.52)


def elastic_issues(model: ConfigModel) -> list[Issue]:
    """medium_elastic_properties.m:155-165: warning when the slowest shear wave has fewer
    than elastic.min_ppw_shear (default 6) points per wavelength. Predicted from the layered
    skull densities; pseudo-CT densities are only known at run time."""
    g = model.get
    if g("acoustic_backend") != "kwave_elastic" or (g("elastic.shear_mapping") or "pichardo") != "pichardo":
        return []
    if (g("simulation.medium") or "layered") != "layered" or _on(g("pct.enabled")):
        return []
    dx = g("grid.resolution_mm")
    txs = model.effective_transducers()
    f = txs[0].get("freq_hz") if txs else None
    f = f[0] if isinstance(f, list) and f else f
    if not _is_num(dx) or dx <= 0 or not _is_num(f) or f <= 0:
        return []
    rho_min = g("elastic.shear_density_min")
    rho_min = float(rho_min) if _is_num(rho_min) else 1200.0
    nu = g("elastic.poisson_min")
    nu = float(nu) if _is_num(nu) else 0.1
    a, b = _pichardo_shear_fit(float(f))
    props = g("medium_properties") or {}
    speeds = []
    for tissue in ("skull", "skull_cortical", "skull_trabecular"):    # tissuemask_binary skull
        p = props.get(tissue)
        if not isinstance(p, dict) or not _is_num(p.get("density")) or p["density"] <= rho_min:
            continue
        cs = a * p["density"] + b
        if _is_num(p.get("sound_speed")):                               # Poisson-ratio cap
            cs = min(cs, p["sound_speed"] * math.sqrt((1 - 2 * nu) / (2 * (1 - nu))))
        speeds.append((cs, tissue))
    if not speeds:
        return []
    cs, tissue = min(speeds)
    need = g("elastic.min_ppw_shear")
    need = float(need) if _is_num(need) else 6.0
    ppw = cs / (f * dx * 1e-3)
    if ppw >= need:
        return []
    return [Issue("warning", "grid.resolution_mm",
                  f"Shear waves in {tissue} (about {cs:.0f} m/s at {f / 1e3:.0f} kHz) get {ppw:.1f} "
                  f"points per wavelength, below elastic.min_ppw_shear = {need:g}: set "
                  f"{cs / (f * need) * 1e3:.2f} mm or finer")]


def all_issues(model: ConfigModel) -> list[Issue]:
    order = {"error": 0, "warning": 1, "info": 2}
    issues = (value_issues(model) + model.unknown_keys()
              + subject_issues(model) + rule_issues(model) + placement_issues(model)
              + library_issues(model) + isppa_issues(model) + defaults_issues(model)
              + uncertainty_issues(model))
    return sorted(issues, key=lambda i: order.get(i.level, 3))
