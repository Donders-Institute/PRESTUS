"""Mirrors functions/helper/predict_outputs.m: which pipeline mode
prestus_pipeline_start.m would dispatch a study to, and the primary files
that mode writes. Pure and read-only — reads the current (possibly
unfinished) study config and never touches MATLAB or disk.

Mode detection mirrors functions/helper/is_async_mode.m and
is_multi_isppa_mode.m exactly (same field paths, same conditions) so this
can't silently drift from what MATLAB would actually do; test_predict_outputs.py
is the drift test, the same role test_schema.py plays for the config schema.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .model import ConfigModel, is_empty
from .study import sanitize_affix
from .uncertainty import PAIRS, VARIANT_SHORT
from .uncertainty import VARIANTS as UNCERTAINTY_VARIANTS

ACOUSTIC_ONLY = PAIRS["trans"]         # uncertainty variants without heating

_MODE_FUNCTION = {
    "async": "async_transducer_pipeline",
    "multi_isppa": "multi_isppa_pipeline",
    "uncertainty": "uncertainty_pipeline",
}


def _on(value, default: bool) -> bool:
    """A 0/1 config flag; unset (None or '') means the default."""
    if value is None or value == "":
        return default
    return bool(value)


def _is_num(v) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _finite_values(v) -> list[float]:
    """Mirrors MATLAB's target_isppa_wcm2(:)' plus an isfinite() filter."""
    raw = v if isinstance(v, list) else [v]
    return [float(x) for x in raw if _is_num(x) and x == x]  # x == x excludes NaN


@dataclass
class Artifact:
    name: str
    note: str


@dataclass
class Plan:
    mode: str                          # 'single' | 'async' | 'multi_isppa' | 'uncertainty'
    nested: str                        # uncertainty only: '' | 'async' | 'multi_isppa' run per variant
    platform_used: str                 # 'matlab' | 'hpc' | 'auto' (only changes the jobs wording)
    tag: str                           # 'sub-003', or 'sub-XXX' when subject_id is unset/invalid
    artifacts: list[Artifact] = field(default_factory=list)
    common_dirs: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    jobs: list[str] = field(default_factory=list)      # deployment structure: stages/jobs this run submits


def is_async_mode(model: ConfigModel) -> bool:
    """Mirrors functions/helper/is_async_mode.m."""
    return model.get("simulation.transducer_coupling") == "async" and len(model.tx) > 1


def is_multi_isppa_mode(model: ConfigModel) -> bool:
    """Mirrors functions/helper/is_multi_isppa_mode.m."""
    for ti in range(len(model.tx)):
        if len(_finite_values(model.get("transducer.target_isppa_wcm2", ti))) > 1:
            return True
    return False


def _trans_suffix(ti: int, n_trans: int) -> str:
    """Mirrors functions/helper/trans_suffix.m."""
    return f"_trans-{ti:02d}" if n_trans > 1 else ""


def _async_sweep_points(model: ConfigModel) -> int:
    """Mirrors async_transducer_pipeline.m's own n_points derivation: the
    longest target_isppa_wcm2 vector across transducers (1 when every
    transducer carries a scalar or unset target)."""
    n_points = 1
    for ti in range(len(model.tx)):
        n_points = max(n_points, len(_finite_values(model.get("transducer.target_isppa_wcm2", ti))) or 1)
    return n_points


_STANCE_FREE = ("timing_offset_s", "timing", "name", "target_isppa_wcm2", "tx_calibration")


def shared_fields(model: ConfigModel) -> list[int]:
    """For each transducer (0-based), the first transducer with the same acoustic field,
    as acoustic_store_key.m: with explicit positions (placement.mode 'manual', or no
    placement stage), entries that differ only in timing, start, name, and target share
    one acoustic simulation. Other placement modes resolve each entry on its own."""
    tx = model.effective_transducers()
    mode = model.get("placement.mode")
    placement_on = model.get("modules.run_transducer_placement") not in (0, False)
    if placement_on and mode not in (None, "", "manual"):
        return list(range(len(tx)))
    stance = [{k: v for k, v in t.items() if k not in _STANCE_FREE} for t in tx]
    return [next(j for j in range(i + 1) if stance[j] == stance[i]) for i in range(len(tx))]


def _async_affix_thermal(base_affix: str, j: int, n_points: int) -> str:
    if n_points == 1:
        return f"{base_affix}_desc-asyncJoint"
    return f"{base_affix}_desc-asyncJointS{j:02d}"


def predict(model: ConfigModel) -> Plan:
    """Mirrors functions/helper/predict_outputs.m."""
    sid_raw = model.get("subject_id")
    sid = int(sid_raw) if _is_num(sid_raw) and sid_raw == int(sid_raw) and 0 <= sid_raw <= 999 else None
    tag = f"sub-{sid:03d}" if sid is not None else "sub-XXX"

    medium = model.get("simulation.medium") or "layered"
    layered = medium == "layered"
    base_affix = sanitize_affix(model.get("io.output_affix") or "")   # as load_parameters does

    async_on = is_async_mode(model)
    multi_isppa_on = is_multi_isppa_mode(model)
    uncertainty_on = bool(model.get("simulation.uncertainty"))

    # platform='auto' is resolved inside prestus_pipeline_start.m; it only changes
    # how the jobs list describes the run.
    platform = model.get("platform")
    if platform == "matlab":
        platform_used = "matlab"
    elif platform in ("slurm", "qsub"):
        platform_used = "hpc"
    else:
        platform_used = "auto"

    # Same order as prestus_pipeline_start.m on every platform. Uncertainty
    # combined with async or a multi-ISPPA sweep nests that pipeline inside
    # each variant (uncertainty_pipeline.m).
    if uncertainty_on:
        mode = "uncertainty"
    elif async_on:
        mode = "async"
    elif multi_isppa_on:
        mode = "multi_isppa"
    else:
        mode = "single"
    nested = ""
    if uncertainty_on and async_on:
        nested = "async"
    elif uncertainty_on and multi_isppa_on:
        nested = "multi_isppa"

    warnings: list[str] = []
    if mode != "single":
        warnings.append(
            f"Assumes the default options of {_MODE_FUNCTION[mode]}; other options can change this list.")

    common_dirs = ["nii/", "img/", "log/", "cache/"]
    if bool(model.get("simulation.debug")):
        common_dirs.append("debug/")

    artifacts: list[Artifact] = []
    jobs: list[str] = []

    if mode == "single":
        artifacts.append(Artifact(f"{tag}_{medium}{base_affix}.csv", "metrics"))
        artifacts.append(Artifact(f"{tag}_{medium}_report{base_affix}.html", "HTML report"))
        if layered and bool(model.get("modules.run_targeting_analysis")):
            artifacts.append(Artifact(f"{tag}_{medium}{base_affix}_targeting.csv", "targeting accuracy"))
        if layered and _on(model.get("report.segmentation"), True):
            artifacts.append(Artifact(f"{tag}_{medium}_segmentation_report{base_affix}.html",
                                      "segmentation report (head-model postprocessing)"))
        if bool(model.get("simulation.debug")) and _on(model.get("report.debug"), True):
            artifacts.append(Artifact(f"{tag}_{medium}_debug_report{base_affix}.html", "debug report"))
        # preproc_transducer_placement.m: an empty placement.report counts as on
        report = model.get("placement.report")
        if _on(model.get("modules.run_transducer_placement"), True) and \
                (report in (None, "", []) or bool(report)):
            mode_name = str(model.get("placement.mode") or "manual").lower()
            artifacts.append(Artifact(f"{tag}_placement_{mode_name}{base_affix}.html", "placement report"))
        # time_reversal_stage.m (k-Wave source setup only; skipped with BabelBrain)
        if _on(model.get("time_reversal.enabled"), False) and \
                str(model.get("acoustic_backend") or "").lower() != "babelbrain":
            artifacts.append(Artifact(f"{tag}_{medium}{base_affix}_timereversal_summary.csv",
                                      "time reversal: refocusing and aberration per transducer"))
            artifacts.append(Artifact(f"{tag}_{medium}{base_affix}_timereversal_elements.csv",
                                      "time reversal: transfer and phases per element"))
        jobs.append("1 job — the full pipeline, start to finish (prestus_pipeline.m)")

    elif mode == "async":
        n_trans = len(model.tx)
        n_points = _async_sweep_points(model)
        rep = shared_fields(model)
        fields = [ti for ti in range(n_trans) if rep[ti] == ti]
        run_water = bool(model.get("modules.run_posthoc_water_sims"))
        if run_water:
            for ti in (i + 1 for i in fields):
                artifacts.append(Artifact(f"{tag}_water{base_affix}{_trans_suffix(ti, n_trans)}.csv",
                                           f"free-water validation — transducer {ti}"))
        for j in range(1, n_points + 1):
            affix = _async_affix_thermal(base_affix, j, n_points)
            label = "joint thermal + report"
            if n_points > 1:
                label = f"{label} (sweep point {j}/{n_points})"
            artifacts.append(Artifact(f"{tag}_{medium}{affix}.csv", label))
            artifacts.append(Artifact(f"{tag}_{medium}_report{affix}.html", label))
        grid_mode = model.get("grid.async_layout") or "shared"
        if grid_mode == "per_transducer":
            jobs.append(f"{n_trans} preprocessing jobs — grid.async_layout='per_transducer', one grid each")
        else:
            jobs.append("1 shared preprocessing stage — grid.async_layout='shared', one grid for all transducers")
        if len(fields) < n_trans:
            jobs.append(f"{len(fields)} acoustic jobs — one per distinct field ({n_trans} transducers; "
                        "transducers with the same stance and settings share one)")
        else:
            jobs.append(f"{n_trans} acoustic jobs — one per transducer")
        jobs.append("fields kept in this output folder by earlier runs with the same settings are "
                    "reused, not simulated again")
        if run_water:
            jobs.append(f"{len(fields)} free-water validation jobs — one per distinct field")
        jobs.append(f"{n_points} combine job(s)" + (" — one per sweep point" if n_points > 1 else ""))
        jobs.append(f"{n_points} joint thermal + report job(s)" + (" — one per sweep point" if n_points > 1 else ""))

    elif mode == "multi_isppa":
        targets = _finite_values(model.get("transducer.target_isppa_wcm2", 0))
        for target in targets:
            affix = f"{base_affix}_desc-isppa{round(target):03d}"
            label = f"{target:.0f} W/cm² — thermal + report"
            artifacts.append(Artifact(f"{tag}_{medium}{affix}.csv", label))
            artifacts.append(Artifact(f"{tag}_{medium}_report{affix}.html", label))
        artifacts.append(Artifact(f"{tag}_{medium}_desc-multiIsppa_report{base_affix}.html",
                                   "summary comparing all targets"))
        jobs.append("1 acoustic + free-water baseline job (stage 1)")
        jobs.append(f"{len(targets)} thermal jobs (stages 2–{len(targets) + 1}, one per target)")
        jobs.append("1 summary report job (final stage)")

    elif mode == "uncertainty":
        # Output affix suffix of each sweep point an inner pipeline writes.
        if nested == "multi_isppa":
            targets = _finite_values(model.get("transducer.target_isppa_wcm2", 0))
            point_suffixes = [f"_desc-isppa{round(t):03d}" for t in targets]
        elif nested == "async":
            n_points = _async_sweep_points(model)
            point_suffixes = [_async_affix_thermal("", j, n_points) for j in range(1, n_points + 1)]
        else:
            point_suffixes = [""]
        n_pts = len(point_suffixes)
        for variant in ("default", *UNCERTAINTY_VARIANTS):
            # output affix as in uncertainty_variants.m: _desc-<name without underscore>
            variant_affix = base_affix if variant == "default" else f"{base_affix}_desc-{variant.replace('_', '')}"
            vlabel = f"{VARIANT_SHORT[variant]} variant"
            if variant in ACOUSTIC_ONLY:
                vlabel += " (acoustic only)"
            for j, suffix in enumerate(point_suffixes, start=1):
                affix = variant_affix + suffix
                label = vlabel + (f", sweep point {j}/{n_pts}" if n_pts > 1 else "")
                artifacts.append(Artifact(f"{tag}_{medium}{affix}.csv", f"{label} metrics"))
                artifacts.append(Artifact(f"{tag}_{medium}_report{affix}.html", f"{label} report"))
            if nested == "multi_isppa":
                artifacts.append(Artifact(f"{tag}_{medium}_desc-multiIsppa_report{variant_affix}.html",
                                           f"{vlabel}: summary comparing all targets"))
        for j, suffix in enumerate(point_suffixes, start=1):
            label = "combined uncertainty report (default, Temp-H/L, Trans-H/L)"
            if n_pts > 1:
                label = f"{label}, sweep point {j}/{n_pts}"
            artifacts.append(Artifact(f"{tag}_{medium}_desc-uncertainty_report{base_affix}{suffix}.html", label))
        if nested == "multi_isppa":
            jobs.append("1 preprocessing job (stage 1)")
            jobs.append("5 multi-ISPPA runs, one per variant (default, Temp-H, Temp-L, Trans-H, Trans-L): "
                        f"default, Temp-H, and Temp-L each 1 acoustic job, {n_pts} thermal jobs and 1 summary "
                        "job; Trans-H and Trans-L run acoustic only, without thermal jobs")
            jobs.append(f"{n_pts} uncertainty report job(s), one per target (stage 7)")
        elif nested == "async":
            jobs.append("5 async multi-transducer runs, one per variant (default, Temp-H, Temp-L, Trans-H, "
                        "Trans-L), each with its own preprocessing, acoustic, combine, and joint jobs; "
                        "Trans-H and Trans-L run acoustic only")
            jobs.append(f"{n_pts} uncertainty report job(s), one per sweep point (stage 7)")
        else:
            jobs.append("1 preprocessing job (stage 1)")
            jobs.append("5 variant simulation jobs: default, Temp-H, Temp-L, Trans-H, Trans-L (stages 2-6); "
                        "Trans-H and Trans-L run acoustic only")
            jobs.append("1 uncertainty report job (stage 7)")

    if mode != "single":
        if platform_used == "hpc":
            jobs.append("On SLURM/qsub these submit as separate jobs with scheduler (afterok) dependencies.")
        elif platform_used == "matlab":
            jobs.append("In this MATLAB session these run sequentially, one after another.")
        else:
            jobs.append("platform='auto': separate SLURM/qsub jobs with dependencies if a scheduler is "
                        "found at submission time, otherwise sequential stages in one MATLAB session.")

    return Plan(mode=mode, nested=nested, platform_used=platform_used, tag=tag,
                artifacts=artifacts, common_dirs=common_dirs, warnings=warnings, jobs=jobs)


@dataclass
class SaveItem:
    label: str                 # 'source matrices'
    kind: str                  # 'file' | 'folder'
    where: str                 # path relative to the output folder, e.g. cache/sub-003_layered_results.mat
    saved: bool                # what the flags decide
    reason: str                # why: the flag that decides, or the default
    applies: bool = True       # False when the stage that writes it does not run
    note: str = ""             # what it holds
    pattern: str = ""          # glob, relative to the output folder, of this run's files (outputs.existing_outputs)


def save_items(model: ConfigModel) -> list[SaveItem]:
    """The optional outputs (io.save_*, simulation.debug) with where they go and why they are
    or are not saved, as should_save_output.m and the writers decide. The metrics CSV, the
    report, log/, img/ and the T1w NIfTI maps are always written (see predict)."""
    g = model.get
    sid_raw = g("subject_id")
    tag = f"sub-{int(sid_raw):03d}" if _is_num(sid_raw) and sid_raw == int(sid_raw) and 0 <= sid_raw <= 999 \
        else "sub-XXX"
    medium = g("simulation.medium") or "layered"
    layered = medium == "layered"
    affix = sanitize_affix(g("io.output_affix") or "")
    stem = f"{tag}_{medium}"

    def on(key, default=True) -> bool:
        v = g(key)
        return default if v is None or (isinstance(v, list) and not v) else bool(v)

    heating = on("modules.run_heating_sims", False)
    nifti = on("modules.run_nifti_creation")
    babelbrain = str(g("acoustic_backend") or "kwave").lower() == "babelbrain"
    timing_pd = not is_empty(g("timing.pd"))
    timing_pri = not is_empty(g("timing.pri"))

    items = []
    for key, label, where, note, applies, why_not in [
        ("save_source_matrices", "source matrices", f"cache/{stem}_kwave_source{affix}.mat",
         "transducer source array, shared across variants", not babelbrain, "BabelBrain sets up its own source"),
        ("save_acoustic_matrices", "acoustic matrices", f"cache/{stem}_results{affix}.mat",
         "sensor data, grid, medium, and source; can be very large", True, ""),
        ("save_thermal_matrices", "thermal matrices", f"cache/{stem}_heating_res{affix}.mat",
         "temperature time series and field maps", heating, "no thermal simulation"),
        ("save_grid_cache", "grid cache", f"cache/{stem}_grid_cache{affix}.mat",
         "preprocessed grid, reused by a re-run", True, ""),
        ("save_medium_cache", "medium cache", f"cache/{stem}_medium_cache{affix}.mat",
         "mapped medium, reused by a re-run", True, ""),
    ]:
        value = g(f"io.{key}")
        saved = _effective_save(g, key)
        if not is_empty(value):
            reason = f"io.{key}: {_flag(value)}"
        elif not is_empty(g("io.save_matrices")):
            reason = f"io.{key} empty: follows io.save_matrices: {_flag(g('io.save_matrices'))}"
        else:
            reason = f"io.{key} and io.save_matrices empty: default {'on' if saved else 'off'}"
        if not applies:
            reason = why_not
        items.append(SaveItem(label, "file", where, saved, reason, applies, note, where))

    items.append(SaveItem("complex pressure field", "file", f"cache/{stem}_p_complex{affix}.mat",
                          on("io.save_p_complex", False), f"io.save_p_complex: {_flag(g('io.save_p_complex'))}",
                          True, "magnitude and phase of the steady-state field (large)",
                          f"cache/{stem}_p_complex{affix}*.mat"))
    v = g("io.save_heatingvideo")
    video = f"{stem}_heating_animation{affix}.avi"
    items.append(SaveItem("heating video", "file", video,
                          v is None or v == 1, f"io.save_heatingvideo: {_flag(v)}" if heating else "no thermal simulation",
                          heating, "incremental heating, uncompressed AVI", video))
    for key, label, kind, where, pattern, note, applies, why_not in [
        ("save_property_maps", "property maps", "folder", "nii/properties/", "nii/properties/*.nii*",
         "one NIfTI per acoustic and thermal property", layered and nifti, "layered medium with NIfTI maps only"),
        ("save_MNI", "MNI-space maps", "file", f"nii/{stem}_MNI{affix}_*.nii.gz", f"nii/{stem}_MNI{affix}_*.nii*",
         "the NIfTI maps also in MNI space (needs the SimNIBS m2m folder)",
         layered and nifti, "layered medium with NIfTI maps only"),
        ("save_Ipa", "pulse-average intensity (Ipa) map", "file", f"nii/{stem}_T1w{affix}_Ipa.nii.gz",
         f"nii/{stem}_T1w{affix}_Ipa.nii*", "", nifti and timing_pd, "needs timing.pd" if nifti else "no NIfTI maps"),
        ("save_Ita", "time-average intensity (Ita) map", "file", f"nii/{stem}_T1w{affix}_Ita.nii.gz",
         f"nii/{stem}_T1w{affix}_Ita.nii*", "", nifti and timing_pd and timing_pri,
         "needs timing.pd and timing.pri" if nifti else "no NIfTI maps"),
    ]:
        value = g(f"io.{key}")
        saved = _effective_save(g, key)
        reason = f"io.{key}: {_flag(value)}" if not is_empty(value) else f"io.{key} empty: default {'on' if saved else 'off'}"
        items.append(SaveItem(label, kind, where, saved, reason if applies else why_not, applies, note, pattern))
    debug = bool(g("simulation.debug"))
    items.append(SaveItem("debug diagnostics", "folder", "debug/", debug, f"simulation.debug: {_flag(g('simulation.debug'))}",
                          True, "preproc/, medium/, source/, and babelbrain/ figures and data, indexed by the debug report", "debug/*"))
    return items


def _flag(v) -> str:
    if is_empty(v):
        return "empty"
    return "1" if bool(v) else "0"


def save_summary(model: ConfigModel) -> list[str]:
    """Plain-text read of save_items, for the command line and the tests."""
    g = model.get
    global_flag = g("io.save_matrices")
    if is_empty(global_flag):
        global_text = "unset — each matrix flag below uses its own default"
    elif bool(global_flag):
        global_text = "ON — saves every matrix below whose flag is left empty"
    else:
        global_text = "OFF — saves no matrix below whose flag is left empty"
    lines = ["Global default (io.save_matrices): " + global_text]
    items = save_items(model)
    width = max(len(i.label) for i in items) + 2
    for i in items:
        if i.label == "debug diagnostics":
            continue
        mark = "✓" if i.saved and i.applies else "✗"
        lines.append(f"  {mark} {i.label:<{width}}{i.where}  ({i.reason})")
    debug_on = bool(g("simulation.debug"))
    lines.append(("✓" if debug_on else "✗") +
                 " Debug diagnostics (simulation.debug): " +
                 ("ON — writes debug/{preproc,medium,source,babelbrain}/ and a debug report" if debug_on else "OFF"))
    return lines


# Built-in defaults of should_save_output.m. Matrix flags inherit
# io.save_matrices when empty; the other io.save_* flags do not.
_MATRIX_SAVE_DEFAULTS = {"save_source_matrices": True, "save_acoustic_matrices": False,
                         "save_thermal_matrices": True, "save_grid_cache": False,
                         "save_medium_cache": False}
_OTHER_SAVE_DEFAULTS = {"save_MNI": True, "save_property_maps": False, "save_kwave_medium": False}


def _effective_save(g, key: str) -> bool:
    """Mirror should_save_output.m: per-step flag, else io.save_matrices (matrix
    flags only), else the built-in default (any flag not listed there: on)."""
    value = g(f"io.{key}")
    if not is_empty(value):
        return bool(value)
    if key in _MATRIX_SAVE_DEFAULTS:
        global_flag = g("io.save_matrices")
        if not is_empty(global_flag):
            return bool(global_flag)
        return _MATRIX_SAVE_DEFAULTS[key]
    return _OTHER_SAVE_DEFAULTS.get(key, True)
