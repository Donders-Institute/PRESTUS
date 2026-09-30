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

from .model import ConfigModel

_MODE_FUNCTION = {
    "async": "async_transducer_pipeline",
    "multi_isppa": "multi_isppa_pipeline",
    "uncertainty": "uncertainty_pipeline",
}


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
    platform_used: str                 # 'matlab' | 'hpc'
    ambiguous: bool
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
    base_affix = model.get("io.output_affix") or ""

    async_on = is_async_mode(model)
    multi_isppa_on = is_multi_isppa_mode(model)
    uncertainty_on = bool(model.get("simulation.uncertainty"))

    # platform='auto' (config_default.yaml's own default!) is resolved to a concrete
    # platform by hpc_detect_system() *inside* prestus_pipeline_start.m — by the time
    # its dispatch switch runs, parameters.platform is never still 'auto'. Predicting
    # from a config alone, before that call, we genuinely don't know which way 'auto'
    # will go: keep it a distinct third state rather than guessing 'hpc'.
    platform = model.get("platform")
    if platform == "matlab":
        platform_used = "matlab"
    elif platform in ("slurm", "qsub"):
        platform_used = "hpc"
    else:
        platform_used = "auto"

    active = {"async": async_on, "multi_isppa": multi_isppa_on, "uncertainty": uncertainty_on}
    precedence_matlab = ["async", "multi_isppa", "uncertainty"]
    precedence_hpc = ["uncertainty", "async", "multi_isppa"]
    mode_matlab = next((m for m in precedence_matlab if active[m]), "single")
    mode_hpc = next((m for m in precedence_hpc if active[m]), "single")

    # A vector target_isppa_wcm2 under async coupling is NOT ambiguous: async mode
    # absorbs it as its own internal sweep (see async_transducer_pipeline.m's
    # header), identically on every platform. The one platform-DEPENDENT case is
    # uncertainty mode competing with async/multi-ISPPA — HPC checks uncertainty
    # first, MATLAB checks it last (see prestus_pipeline_start.m) — so only that
    # combination is flagged, and only then does platform_used=='auto' matter.
    ambiguous = uncertainty_on and (async_on or multi_isppa_on)
    if platform_used == "matlab":
        mode = mode_matlab
    elif platform_used == "hpc":
        mode = mode_hpc
    else:
        mode = mode_matlab   # representative pick; see the 'auto' warning below when ambiguous

    warnings: list[str] = []
    if ambiguous:
        active_names = [n for n in ("async", "multi_isppa", "uncertainty") if active[n]]
        base = ("More than one pipeline trigger is active at once (" + " + ".join(active_names) +
                "). prestus_pipeline_start.m resolves this differently by platform: async > "
                "multi-ISPPA > uncertainty on 'matlab', but uncertainty > async > multi-ISPPA on "
                "'slurm'/'qsub'.")
        if platform_used == "auto":
            warnings.append(
                f"{base} platform is 'auto': PRESTUS decides at submission time (hpc_detect_system), "
                f"so this could run as '{mode_matlab}' mode (no scheduler found) or as '{mode_hpc}' "
                f"mode (SLURM/qsub found). This preview shows '{mode}' mode; set platform explicitly "
                "to pin down which one.")
        else:
            warnings.append(
                f"{base} For platform_used='{platform_used}' this prediction assumes '{mode}' mode "
                "runs; the others are ignored here, not combined.")
    if mode != "single":
        warnings.append(
            f"Assumes {_MODE_FUNCTION[mode]} runs with its default options (see its own help text); "
            "a non-default call to prestus_pipeline_start(parameters, options) can add to or remove "
            "from this list.")

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
        jobs.append("1 job — the full pipeline, start to finish (prestus_pipeline.m)")

    elif mode == "async":
        n_trans = len(model.tx)
        n_points = _async_sweep_points(model)
        run_water = bool(model.get("modules.run_posthoc_water_sims"))
        if run_water:
            for ti in range(1, n_trans + 1):
                artifacts.append(Artifact(f"{tag}_water{base_affix}{_trans_suffix(ti, n_trans)}.csv",
                                           f"free-water validation — transducer {ti}"))
        for j in range(1, n_points + 1):
            affix = _async_affix_thermal(base_affix, j, n_points)
            label = "joint thermal + report"
            if n_points > 1:
                label = f"{label} (sweep point {j}/{n_points})"
            artifacts.append(Artifact(f"{tag}_{medium}{affix}.csv", label))
            artifacts.append(Artifact(f"{tag}_{medium}_report{affix}.html", label))
        grid_mode = model.get("simulation.async_grid_mode") or "shared"
        if grid_mode == "per_transducer":
            jobs.append(f"{n_trans} preprocessing jobs — async_grid_mode='per_transducer', one grid each")
        else:
            jobs.append("1 shared preprocessing stage — async_grid_mode='shared', one grid for all transducers")
        jobs.append(f"{n_trans} acoustic jobs — one per transducer")
        if run_water:
            jobs.append(f"{n_trans} free-water validation jobs — one per transducer")
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
        for variant in ("default", "liberal", "conservative"):
            affix = base_affix if variant == "default" else f"{base_affix}_desc-{variant}"
            artifacts.append(Artifact(f"{tag}_{medium}{affix}.csv", f"{variant} variant metrics"))
            artifacts.append(Artifact(f"{tag}_{medium}_report{affix}.html", f"{variant} variant report"))
        artifacts.append(Artifact(f"{tag}_{medium}_uncertainty_report.html",
                                   "combined uncertainty report (default/liberal/conservative)"))
        jobs.append("1 preprocessing job (stage 1)")
        jobs.append("3 variant simulation jobs: default, liberal, conservative (stages 2–4)")
        jobs.append("1 uncertainty report job (stage 5)")

    if mode != "single":
        if platform_used == "hpc":
            jobs.append("On SLURM/qsub these submit as separate jobs with scheduler (afterok) dependencies.")
        elif platform_used == "matlab":
            jobs.append("In this MATLAB session these run sequentially, one after another.")
        else:
            jobs.append("platform='auto': separate SLURM/qsub jobs with dependencies if a scheduler is "
                        "found at submission time, otherwise sequential stages in one MATLAB session.")

    return Plan(mode=mode, platform_used=platform_used, ambiguous=ambiguous, tag=tag,
                artifacts=artifacts, common_dirs=common_dirs, warnings=warnings, jobs=jobs)


def save_summary(model: ConfigModel) -> list[str]:
    """Human-readable read of io.save_* / simulation.debug: what gets kept besides
    the always-written CSV/report. Mirrors config_default.yaml's io section
    comments — a direct readout, not a prediction, so it needs no mode branching."""
    g = model.get
    global_default = bool(g("io.save_matrices"))
    lines = [
        "Global default (io.save_matrices): " +
        ("ON — saves everything below unless a flag overrides it" if global_default
         else "OFF — nothing below is saved unless a flag overrides it"),
    ]
    flags = [
        ("save_source_matrices", "source matrices", "kwave_source.mat — shared across variants"),
        ("save_acoustic_matrices", "acoustic matrices", "sensor_data + kgrid + medium + source (can be large)"),
        ("save_thermal_matrices", "thermal matrices", "heating_res.mat — temperature timeseries + field maps"),
        ("save_heatingvideo", "heating video", "incremental heating video"),
        ("save_property_maps", "property maps", "nii/properties/*.nii.gz (layered medium only)"),
        ("save_MNI", "MNI-space NIfTIs", "in addition to native T1w space"),
        ("save_Ipa", "Ipa intensity map", "pulse-average intensity NIfTI"),
        ("save_Ita", "Ita intensity map", "temporal-average intensity NIfTI"),
        ("save_p_complex", "complex pressure field", "p_complex alongside p_max_all (large)"),
    ]
    width = max(len(label) for _, label, _ in flags) + 2
    for key, label, note in flags:
        on = bool(g(f"io.{key}"))
        mark = "✓" if on else "✗"
        lines.append(f"  {mark} {label:<{width}}{note}")
    debug_on = bool(g("simulation.debug"))
    lines.append(("✓" if debug_on else "✗") +
                 " Debug diagnostics (simulation.debug): " +
                 ("ON — writes debug/{preproc,medium,source}/" if debug_on else "OFF"))
    return lines
