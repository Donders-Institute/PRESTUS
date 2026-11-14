"""What a run of the current study config will do, stage by stage, for the Run page.

Mirrors the stage gating of functions/prestus_pipeline.m (segmentation → placement → grid →
medium → source → water baseline → acoustic → analysis/NIfTI/targeting → thermal → thermal
analysis → post-hoc water → report). Pure and read-only, like predict_outputs.py: it reads
the (possibly unsaved) config and only checks whether inputs already exist on disk (SimNIBS
segmentation, pseudo-CT, measured CT), which decides between generating and reusing them.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .model import ConfigModel, is_empty

RUN, REUSE, SKIP, STOP = "run", "reuse", "skip", "stop"


@dataclass
class Stage:
    name: str
    state: str                  # RUN | REUSE | SKIP | STOP
    detail: str = ""
    icon: str = ""              # glyph of what the stage does (ui/icons.py), set by plan()


@dataclass
class RunPlan:
    stages: list[Stage] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


def _on(v, default: bool = True) -> bool:
    """modules.* flag as prestus_pipeline.m reads it: a missing flag counts as on where the
    pipeline uses ~isfield(...) || flag, so callers pass that default."""
    if v is None or (isinstance(v, list) and not v):
        return default
    if isinstance(v, str):
        return v.lower() not in ("0", "false", "off", "")
    return bool(v)


def _subject(model: ConfigModel) -> int | None:
    sid = model.get("subject_id")
    if isinstance(sid, (int, float)) and not isinstance(sid, bool) and sid == int(sid) and 0 <= sid <= 999:
        return int(sid)
    return None


def _num(v):
    if isinstance(v, list):
        v = v[0] if v else None
    return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else None


def _targets(model: ConfigModel) -> list[float]:
    v = model.get("transducer.target_isppa_wcm2")
    raw = v if isinstance(v, list) else [v]
    return [float(x) for x in raw if isinstance(x, (int, float)) and not isinstance(x, bool) and x == x]


def _segmentation(model: ConfigModel, sid: int | None, layered: bool) -> Stage:
    if not layered:
        return Stage("Segmentation", SKIP, "no head model for this medium")
    from .placement import find_segmentation
    seg = model.get("path.seg")
    seg_root = model.get("path.anat") if is_empty(seg) else seg      # load_parameters default
    found = find_segmentation(str(seg_root), sid) if not is_empty(seg_root) and sid is not None else None
    if found is not None and not _on(model.get("io.overwrite_simnibs"), False):
        return Stage("Segmentation", REUSE, f"SimNIBS charm output {found.parent.name}/")
    t2 = model.get("path.t2_pattern")
    return Stage("Segmentation", RUN, "SimNIBS charm from the T1" + ("" if is_empty(t2) else " and T2/UTE")
                 + " (takes hours)")


def _skull_ct(model: ConfigModel, layered: bool) -> Stage | None:
    if not layered:
        return None
    if not _on(model.get("pct.enabled"), False):
        return Stage("Skull CT", SKIP, "skull properties uniform, from the segmentation")
    if not is_empty(model.get("path.ct_pattern")):
        from .outputs import measured_ct_file
        f = measured_ct_file(model)
        return Stage("Skull CT", REUSE if f else STOP,
                     f"measured CT {f.name}, resampled to the segmentation grid" if f
                     else "path.ct_pattern matches no file: the run stops")
    from .outputs import pseudoct_file
    f, _legacy = pseudoct_file(model)
    if f is not None:
        return Stage("Skull CT", REUSE, f"pseudo-CT {f.parent.name}/{f.name}")
    mapping = model.get("pct.skull_mapping") or "kosciessa"
    return Stage("Skull CT", RUN, f"generate a pseudo-CT from the UTE image ({mapping} mapping; FSL, ANTs)")


def _placement(model: ConfigModel) -> Stage:
    if not _on(model.get("modules.run_transducer_placement")):
        return Stage("Placement", SKIP, "trans_pos / focus_pos used as configured")
    mode = model.get("placement.mode") or "manual"
    what = {"manual": "configured trans_pos / focus_pos",
            "localite": "Localite neuronavigation markers",
            "heuristic": "heuristic search on the scalp",
            "mni": "MNI coordinates mapped to the subject",
            "plantus": "PLANTUS planning",
            "tus_entry": "TUS_entry planning"}.get(str(mode), str(mode))
    return Stage("Placement", RUN, f"{mode}: {what}")


def _medium(model: ConfigModel, medium: str) -> str:
    if medium == "water":
        return "water"
    if medium == "phantom":
        return "phantom labels mapped to tissues"
    skull = "skull from the CT" if _on(model.get("pct.enabled"), False) else "uniform skull"
    return f"layered head model, {skull}"


def _water_baseline(model: ConfigModel, backend: str) -> Stage:
    if backend == "babelbrain":
        return Stage("Free-water baseline", SKIP, "BabelBrain backend")
    mode = model.get("modules.run_water_baseline")
    targets = _targets(model)
    target = targets[0] if targets else None
    calib = _num(model.get("transducer.tx_calibration.isppa_wcm2"))
    if mode in (1, True, "1"):
        return Stage("Free-water baseline", RUN, "always (run_water_baseline: 1)")
    if target is None:
        return Stage("Free-water baseline", SKIP, "no target ISPPA: results stay at the configured drive")
    differs = calib is None or abs(target - calib) > 0.01 * max(abs(calib), 1e-9)
    if mode in (0, False, "0") and calib is None:
        return Stage("Free-water baseline", STOP, "target set, no calibration ISPPA and run_water_baseline: 0")
    if differs:
        why = "no calibration ISPPA recorded" if calib is None else f"calibration is {calib:g} W/cm²"
        return Stage("Free-water baseline", RUN, f"to scale to {target:g} W/cm² ({why})")
    return Stage("Free-water baseline", SKIP, f"target {target:g} W/cm² equals the calibration")


# What each stage does, as a glyph (ui/icons.py); placement, medium and source follow the config
STAGE_ICONS = {
    "Segmentation": "segmentation", "Everything after segmentation": "warning", "Skull CT": "ct",
    "Grid and head model": "grid", "Free-water baseline": "tissue_water", "Time reversal": "redo",
    "Acoustic simulation": "wave", "Acoustic analysis": "chart", "NIfTI maps": "layers",
    "Targeting analysis": "target", "Thermal simulation": "thermometer", "Thermal analysis": "chart",
    "Post-hoc free-water run": "tissue_water", "Report": "report", "Group report": "report",
}


def _icons(stages: list[Stage], model: ConfigModel) -> None:
    g = model.get
    for s in stages:
        if s.name == "Placement":
            mode = str(g("placement.mode") or "manual")
            s.icon = f"place_{mode}"
        elif s.name == "Medium":
            s.icon = f"medium_{g('simulation.medium') or 'layered'}"
        elif s.name == "Source":
            s.icon = "tx_matrix" if g("transducer.type") == "matrix" else "tx_annular"
        else:
            s.icon = STAGE_ICONS.get(s.name, "")


def plan(model: ConfigModel) -> RunPlan:
    """The stages prestus_pipeline.m runs for this config, in order."""
    out = _plan(model)
    _icons(out.stages, model)
    return out


def _plan(model: ConfigModel) -> RunPlan:
    g = model.get
    out = RunPlan()
    st = out.stages
    sid = _subject(model)
    medium = g("simulation.medium") or "layered"
    layered = medium == "layered"
    head = medium in ("layered", "phantom")
    backend = str(g("acoustic_backend") or "kwave")
    thermal_backend = str(g("thermal_backend") or "kwave")

    st.append(_segmentation(model, sid, layered))
    if _on(g("modules.segmentation_only"), False):
        st.append(Stage("Everything after segmentation", STOP, "modules.segmentation_only is on"))
        return out
    ct = _skull_ct(model, layered)
    if ct is not None:
        st.append(ct)
    st.append(_placement(model))

    res = g("grid.resolution_mm")
    st.append(Stage("Grid and head model", RUN if _on(g("modules.run_grid_setup")) else SKIP,
                    f"{g('grid.mode') or 'focal_axis'} grid, {res:g} mm" if isinstance(res, (int, float))
                    else str(g("grid.mode") or "")))
    st.append(Stage("Medium", RUN if _on(g("modules.run_medium_setup")) else SKIP, _medium(model, medium)))

    acoustic_input = g("path.acoustic_input")
    if not is_empty(acoustic_input):
        st.append(Stage("Acoustic simulation", REUSE, f"start from the pressure map {acoustic_input}"))
    else:
        source_on = _on(g("modules.run_source_setup"))
        n_tx = len(model.tx)
        name = g("transducer.name")
        tx = f"{n_tx} transducers" if n_tx > 1 else (str(name) if not is_empty(name) else "1 transducer")
        st.append(Stage("Source", SKIP if backend == "babelbrain" or not source_on else RUN,
                        "BabelBrain sets up its own source" if backend == "babelbrain" else tx))
        st.append(_water_baseline(model, backend))
        if _on(g("time_reversal.enabled"), False):
            only = _on(g("time_reversal.only"), False)
            st.append(Stage("Time reversal", RUN, "refocusing / aberration report"
                            + (": no forward simulation (time_reversal.only)" if only else "")))
            if only:
                return out
        acoustic_on = _on(g("modules.run_acoustic_sims")) and (source_on or backend == "babelbrain")
        st.append(Stage("Acoustic simulation", RUN if acoustic_on else SKIP,
                        f"{backend}" + (f", {g('simulation.code_type')}" if backend.startswith("kwave")
                                        and not is_empty(g("simulation.code_type")) else "")))
    st.append(Stage("Acoustic analysis", RUN if _on(g("modules.run_acoustic_analysis")) else SKIP,
                    "intensity, pressure and MI metrics"))
    nifti = _on(g("modules.run_nifti_creation"))
    st.append(Stage("NIfTI maps", RUN if nifti else SKIP, "maps in T1w and MNI space"))
    if head:
        targeting = _on(g("modules.run_targeting_analysis")) and nifti
        st.append(Stage("Targeting analysis", RUN if targeting else SKIP,
                        "simulated vs planned focus" if targeting or not _on(g("modules.run_targeting_analysis"))
                        else "needs NIfTI maps"))

    heating = _on(g("modules.run_heating_sims"), False)
    st.append(Stage("Thermal simulation", RUN if heating else SKIP,
                    thermal_backend if heating else "modules.run_heating_sims is off"))
    st.append(Stage("Thermal analysis", RUN if heating and _on(g("modules.run_thermal_analysis")) else SKIP,
                    "heating and CEM43 thermal dose" if heating else "no heating simulated"))
    if head:
        water = _on(g("modules.run_posthoc_water_sims"), False)
        st.append(Stage("Post-hoc free-water run", RUN if water else SKIP,
                        "the same setup in water, as a separate run" if water else ""))
    st.append(Stage("Report", RUN if _on(g("modules.generate_report"), False) else SKIP,
                    f"HTML report ({g('modules.report_style') or 'instrument'})"))
    if _on(g("modules.generate_group_report"), False):
        st.append(Stage("Group report", RUN, "updated after the run"))
    return out
