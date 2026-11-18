"""MATLAB launch: find MATLAB, write a run script, build -batch commands."""

from __future__ import annotations

import datetime as _dt
import glob
import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path

from .fsio import write_text_atomic

OK_SENTINEL = "[PRESTUS Nexus] OK"
FAIL_SENTINEL = "[PRESTUS Nexus] FAILED"


def find_matlab(*hints: str | None) -> str | None:
    """MATLAB binary: explicit hints, $PRESTUS_MATLAB, PATH, then standard install folders."""
    cands = [h for h in hints if h] + [os.environ.get("PRESTUS_MATLAB"), shutil.which("matlab")]
    cands += sorted(glob.glob("/Applications/MATLAB_R20*.app/bin/matlab"), reverse=True)
    cands += sorted(glob.glob("/usr/local/MATLAB/R20*/bin/matlab"), reverse=True)
    cands += sorted(glob.glob(r"C:\Program Files\MATLAB\R20*\bin\matlab.exe"), reverse=True)
    for c in cands:
        if c and Path(c).is_file():
            return str(c)
    return None


def _q(text: str | Path) -> str:
    """MATLAB single-quoted char literal."""
    return "'" + str(text).replace("'", "''") + "'"


def launcher_path(study: Path) -> Path:
    stem = study.stem[len("config_"):] if study.stem.startswith("config_") else study.stem
    return study.with_name(f"run_{stem}.m")


# ── pipeline options (prestus_pipeline_start's second argument) ─────────────────────
# Nexus-side settings, not config keys: stored per study in the session (load_options).
# Stage limits of the multi-ISPPA pipeline (multi_isppa_pipeline.m:110-118) and of the
# uncertainty pipeline (uncertainty_pipeline.m:122-131); a memory limit is in GB like
# hpc.memorylimit, [] keeps hpc.memorylimit.
MULTI_ISPPA_OPTIONS = ("acoustic_timelimit", "thermal_timelimit", "summary_timelimit",
                       "acoustic_memorylimit", "thermal_memorylimit", "summary_memorylimit",
                       "acoustic_partition", "summary_partition")
UNCERTAINTY_OPTIONS = ("stage1_timelimit", "sim_timelimit", "report_timelimit",
                       "stage1_memorylimit", "sim_memorylimit", "report_memorylimit",
                       "stage1_partition", "report_partition")
STAGE_OPTIONS = MULTI_ISPPA_OPTIONS + UNCERTAINTY_OPTIONS
SESSION_KEY = "pipeline_options"


@dataclass
class PipelineOptions:
    """What the run script passes as ``options`` to prestus_pipeline_start.

    sequential_configs: follow-up configs run after this one with heat carry-over, in order
    (options.sequential_configs.config_<N>, N = 2, 3, ...; sequential_pipeline.m picks the
    lowest N first). stages: STAGE_OPTIONS name -> value. generate_summary_report: the
    multi-ISPPA summary report (multi_isppa_pipeline.m:110; None = PRESTUS' default, on)."""
    sequential_configs: list[str] = field(default_factory=list)
    stages: dict[str, object] = field(default_factory=dict)
    generate_summary_report: bool | None = None

    def is_empty(self) -> bool:
        return not self.sequential_configs and not self.stages and self.generate_summary_report is None

    def to_dict(self) -> dict:
        out: dict = {}
        if self.sequential_configs:
            out["sequential_configs"] = list(self.sequential_configs)
        if self.stages:
            out["stages"] = dict(self.stages)
        if self.generate_summary_report is not None:
            out["generate_summary_report"] = bool(self.generate_summary_report)
        return out

    @classmethod
    def from_dict(cls, d) -> "PipelineOptions":
        d = d if isinstance(d, dict) else {}
        seq = [str(x) for x in d.get("sequential_configs") or [] if x]
        stages = {k: v for k, v in (d.get("stages") or {}).items()
                  if k in STAGE_OPTIONS and v not in (None, "", [])}
        gsr = d.get("generate_summary_report")
        return cls(seq, stages, gsr if isinstance(gsr, bool) else None)


def _study_key(study: Path) -> str:
    return str(Path(study).resolve())


def load_options(session, study: Path | None) -> PipelineOptions:
    """The study's pipeline options from the session (empty when none are stored)."""
    if session is None or study is None:
        return PipelineOptions()
    store = session.get(SESSION_KEY) or {}
    return PipelineOptions.from_dict(store.get(_study_key(study)) if isinstance(store, dict) else None)


def save_options(session, study: Path, options: PipelineOptions) -> None:
    """Store (or, when empty, forget) the study's pipeline options; the caller saves the session."""
    store = dict(session.get(SESSION_KEY) or {})
    if options.is_empty():
        store.pop(_study_key(study), None)
    else:
        store[_study_key(study)] = options.to_dict()
    session.set(SESSION_KEY, store)


def add_to_chain(session, study: Path, child: Path) -> PipelineOptions:
    """Append a follow-up config to the study's sequential chain (no duplicates)."""
    opts = load_options(session, study)
    child_s = str(Path(child).resolve())
    if child_s not in opts.sequential_configs:
        opts.sequential_configs.append(child_s)
    save_options(session, study, opts)
    return opts


def options_lines(options: PipelineOptions | None, local: bool, load_expr=None) -> list[str]:
    """MATLAB lines that build ``options`` (none when there are no options).

    load_expr(path) -> MATLAB expression that loads a sequential config; the default loads
    it from its local folder (remote.submit_command maps it to the cluster instead)."""
    if options is None or options.is_empty():
        return []
    lines = ["options = struct();"]
    for key in STAGE_OPTIONS:
        if key in options.stages:
            lines.append(f"options.{key} = {matlab_value(options.stages[key])};")
    if options.generate_summary_report is not None:
        lines.append(f"options.generate_summary_report = {matlab_value(bool(options.generate_summary_report))};")
    if options.sequential_configs:
        lines += ["% sequential follow-ups with heat carry-over (doc_sequential.md)",
                  "options.sequential_configs = struct();"]
        for n, cfg in enumerate(options.sequential_configs, start=2):
            p = Path(cfg).resolve()
            field_ = f"options.sequential_configs.config_{n}"
            expr = load_expr(p) if load_expr else f"load_parameters({_q(p.name)}, {_q(p.parent)})"
            lines.append(f"{field_} = {expr};")
            if local:
                lines.append(f"{field_}.platform = 'matlab';")
    return lines


def launcher_script(study: Path, prestus_root: Path, local: bool = True,
                    options: PipelineOptions | None = None) -> str:
    """Script that loads the study config and starts the pipeline.

    local=True sets platform 'matlab' for this run only (PRESTUS' 'auto'
    requires sbatch or qsub); the study config is left unchanged. ``options``
    (PipelineOptions) are passed to prestus_pipeline_start when there are any.
    """
    study = Path(study).resolve()
    lines = [
        f"% Generated by the PRESTUS Nexus on {_dt.date.today().isoformat()}.",
        f"% Runs {study.name}. Edit the study config, not this script.",
        "",
        f"prestus_root = {_q(prestus_root)};",
        "addpath(genpath(fullfile(prestus_root, 'functions')));",
        "",
        f"parameters = load_parameters({_q(study.name)}, {_q(study.parent)});",
    ]
    if local:
        lines += ["parameters.platform = 'matlab';   % run in this MATLAB session"]
    opt = options_lines(options, local)
    if opt:
        lines += opt + ["prestus_pipeline_start(parameters, options);", ""]
    else:
        lines += ["prestus_pipeline_start(parameters);", ""]
    return "\n".join(lines)


def write_launcher(study: Path, prestus_root: Path, local: bool = True,
                   options: PipelineOptions | None = None) -> Path:
    path = launcher_path(study)
    write_text_atomic(path, launcher_script(study, prestus_root, local, options))
    return path


def submit_launcher_path(study: Path) -> Path:
    stem = study.stem[len("config_"):] if study.stem.startswith("config_") else study.stem
    return study.with_name(f"submit_{stem}.m")


def write_submit_launcher(study: Path, prestus_root: Path,
                          options: PipelineOptions | None = None) -> Path:
    """Launcher that keeps the config's platform, so PRESTUS submits to SLURM/qsub
    (prestus_pipeline_start -> hpc_submit_job) and returns after submission."""
    path = submit_launcher_path(study)
    write_text_atomic(path, launcher_script(study, prestus_root, local=False, options=options))
    return path


def _wrapped(body: str) -> str:
    """-batch statement that prints a sentinel on success and on failure."""
    return (f"try, {body}; fprintf('\\n{OK_SENTINEL}\\n'); "
            f"catch err, fprintf(2, '\\n{FAIL_SENTINEL}: %s\\n', err.message); "
            f"fprintf(2, '%s\\n', getReport(err, 'extended', 'hyperlinks', 'off')); exit(1); end")


def validate_command(matlab: str, study: Path, prestus_root: Path) -> list[str]:
    """Load the config with load_parameters only (no simulation)."""
    study = Path(study).resolve()
    body = (f"addpath(genpath(fullfile({_q(prestus_root)}, 'functions'))); "
            f"p = load_parameters({_q(study.name)}, {_q(study.parent)}); "
            "fprintf('Loaded: subject %d, %d transducer(s), medium %s\\n', "
            "p.subject_id, numel(p.transducer), string(p.simulation.medium))")
    return [matlab, "-batch", _wrapped(body)]


def run_command(matlab: str, launcher: Path) -> list[str]:
    launcher = Path(launcher).resolve()
    return [matlab, "-sd", str(launcher.parent), "-batch", _wrapped(f"run({_q(launcher)})")]


# ── PRESTUS' MATLAB GUIs (need an interactive MATLAB session, not -batch) ──

def gui_command(matlab: str, prestus_root: Path, call: str) -> list[str]:
    """Start MATLAB interactively and open a PRESTUS GUI."""
    return [matlab, "-nosplash", "-r",
            f"addpath(genpath(fullfile({_q(prestus_root)}, 'functions'))); {call}"]


def neuronav_gui_call(study: Path | None) -> str:
    """neuronav_gui pre-filled from the study config (Localite paths, coil map)."""
    return f"neuronav_gui({_q(Path(study).resolve())})" if study else "neuronav_gui()"


def calibration_gui_call(output_path: str | None) -> str:
    """calibration_gui saving entries to transducer.library_path, if set."""
    if output_path:
        return f"calibration_gui(struct('OutputPath', {_q(output_path)}))"
    return "calibration_gui()"


# ── study (multi-subject) runs: prestus_study_start / prestus_group_report_start ──

def study_script(study: Path, prestus_root: Path, local: bool = True, dry_run: bool = False) -> str:
    """Script that runs every manifest row of the study config (see prestus_study_start)."""
    study = Path(study).resolve()
    opts = []
    if dry_run:
        opts.append("'DryRun', true")
    if local:
        opts.append("'Platform', 'matlab'")
    call = f"prestus_study_start({_q(study)}" + (", " + ", ".join(opts) if opts else "") + ");"
    return "\n".join([
        f"% Generated by the PRESTUS Nexus on {_dt.date.today().isoformat()}.",
        f"% Runs all rows of the manifest named in {study.name}. Edit the config/manifest, not this script.",
        "",
        f"prestus_root = {_q(prestus_root)};",
        "addpath(genpath(fullfile(prestus_root, 'functions')));",
        "",
        call,
        "",
    ])


def study_launcher_path(study: Path) -> Path:
    stem = study.stem[len("config_"):] if study.stem.startswith("config_") else study.stem
    return study.with_name(f"run_study_{stem}.m")


def write_study_launcher(study: Path, prestus_root: Path, local: bool = True) -> Path:
    path = study_launcher_path(study)
    write_text_atomic(path, study_script(study, prestus_root, local))
    return path


def study_dry_run_command(matlab: str, study: Path, prestus_root: Path) -> list[str]:
    """Load and check every manifest row, print the plan; runs nothing."""
    study = Path(study).resolve()
    body = (f"addpath(genpath(fullfile({_q(prestus_root)}, 'functions'))); "
            f"prestus_study_start({_q(study)}, 'DryRun', true)")
    return [matlab, "-batch", _wrapped(body)]


def group_report_command(matlab: str, study: Path, prestus_root: Path) -> list[str]:
    """Group MNI maps + group report from existing outputs (no simulation)."""
    study = Path(study).resolve()
    body = (f"addpath(genpath(fullfile({_q(prestus_root)}, 'functions'))); "
            f"rp = prestus_group_report_start({_q(study)}); "
            "if isempty(rp), error('no group report was written'); end")
    return [matlab, "-batch", _wrapped(body)]


# ── Localite batch ingestion and the k-Wave drive check (long: Run tab, not the bridge) ──

def _stem(study: Path) -> str:
    return study.stem[len("config_"):] if study.stem.startswith("config_") else study.stem


PLACEMENT_ONLY = ("for f = fieldnames(parameters.modules)', if startsWith(f{1}, {'run_', 'generate_'}), "
                  "parameters.modules.(f{1}) = 0; end, end; parameters.modules.run_transducer_placement = 1;")


def matlab_value(v) -> str:
    """A Python value as a MATLAB literal (str, bool, number, flat list)."""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return repr(v)
    if isinstance(v, (list, tuple)):
        return "[" + " ".join(matlab_value(x) for x in v) + "]"
    return _q(v)


def placement_overrides(mode: str = "plantus", interactive: bool = False) -> list[tuple[str, object]]:
    """Config values a placement-only run sets on top of the study config: the placement
    always runs again, without interaction unless PlanTUS is to open Workbench."""
    out: list[tuple[str, object]] = [("simulation.interactive", 0)]
    if mode == "plantus":
        out.append(("placement.plantus.interactive", int(interactive)))
    return out + [("placement.reuse", "never")]


def plantus_overrides(interactive: bool) -> list[tuple[str, object]]:
    return placement_overrides("plantus", interactive)


def partial_script(study: Path, prestus_root: Path, what: str, overrides: list[tuple[str, object]],
                   extra: str = "", local: bool = True) -> str:
    """Script that runs part of the pipeline for the study config: ``overrides`` replace config
    values (load_parameters), ``extra`` MATLAB statements then adjust ``parameters`` (e.g.
    PLACEMENT_ONLY). local=False keeps the config's platform, so PRESTUS submits it as a job."""
    study = Path(study).resolve()
    cell = "; ...\n".join(f"    {_q(k)}, {matlab_value(v)}" for k, v in overrides)
    return "\n".join([
        f"% Generated by the PRESTUS Nexus on {_dt.date.today().isoformat()}.",
        f"% Runs only {what} of {study.name}. Edit the study config, not this script.",
        "",
        f"prestus_root = {_q(prestus_root)};",
        "addpath(genpath(fullfile(prestus_root, 'functions')));",
        "",
        f"parameters = load_parameters({_q(study.name)}, {_q(study.parent)}, {{ ...",
        cell + "});",
        *(["parameters.platform = 'matlab';   % run in this MATLAB session"] if local else []),
        *([f"% {what} only", extra] if extra else []),
        "prestus_pipeline_start(parameters);",
        "",
    ])


def write_partial_script(study: Path, prestus_root: Path, name: str, what: str,
                         overrides: list[tuple[str, object]], extra: str = "", local: bool = True) -> Path:
    path = Path(study).with_name(f"{name}_{_stem(Path(study))}.m")
    write_text_atomic(path, partial_script(study, prestus_root, what, overrides, extra, local))
    return path


def segmentation_overrides(rerun: bool = False) -> list[tuple[str, object]]:
    """A segmentation-only run (SimNIBS charm, modules.segmentation_only): the segmentation
    stage runs for a layered medium only; rerun=True segments again over an existing m2m
    folder (io.overwrite_simnibs), which charm deletes first."""
    return [("simulation.interactive", 0), ("simulation.medium", "layered"),
            ("modules.segmentation_only", 1), ("io.overwrite_simnibs", int(rerun))]


def write_segmentation_script(study: Path, prestus_root: Path, rerun: bool = False,
                              local: bool = True) -> Path:
    return write_partial_script(study, prestus_root, "segmentation", "the SimNIBS segmentation",
                                segmentation_overrides(rerun), local=local)


def placement_script(study: Path, prestus_root: Path, mode: str, interactive: bool = False,
                     local: bool = True) -> str:
    """Script that runs only the placement stage of the study config (placement_overrides,
    every other stage off). local=False keeps the config's platform, so PRESTUS submits it
    as a job (HPC)."""
    return partial_script(study, prestus_root, f"the {mode} placement",
                          placement_overrides(mode, interactive), PLACEMENT_ONLY, local)


def placement_script_path(study: Path, mode: str) -> Path:
    return Path(study).with_name(f"{mode}_{_stem(Path(study))}.m")


def write_placement_script(study: Path, prestus_root: Path, mode: str, interactive: bool = False,
                           local: bool = True) -> Path:
    return write_partial_script(study, prestus_root, mode, f"the {mode} placement",
                                placement_overrides(mode, interactive), PLACEMENT_ONLY, local)


def plantus_script(study: Path, prestus_root: Path, interactive: bool, local: bool = True) -> str:
    """Script that runs only the PlanTUS placement of the study config.

    interactive=True opens Connectome Workbench to pick the vertex and confirm it in a
    review dialog (no MATLAB desktop needed); False takes the composite-score optimum, or
    placement.plantus.vertex when the config sets one. The placement always runs again
    (placement.reuse: never) and writes the placement report and the positions the
    placement viewer shows; no simulation runs. local=False keeps the config's platform,
    so PRESTUS submits it as a job (HPC).
    """
    return placement_script(study, prestus_root, "plantus", interactive, local)


def plantus_script_path(study: Path) -> Path:
    return placement_script_path(study, "plantus")


def write_plantus_script(study: Path, prestus_root: Path, interactive: bool, local: bool = True) -> Path:
    return write_placement_script(study, prestus_root, "plantus", interactive, local)


def ingest_script(study: Path, prestus_root: Path, selection: list[tuple[str, str]]) -> str:
    """Script that ingests the selected sessions with the study config's settings
    (neuronav_ingest_batch; outputs and neuronav_ingest_summary.csv under path.localite_post)."""
    study = Path(study).resolve()
    rows = "; ".join(f"{_q(sub)}, {_q(ses)}" for sub, ses in selection)
    return "\n".join([
        f"% Generated by the PRESTUS Nexus on {_dt.date.today().isoformat()}.",
        f"% Ingests {len(selection)} Localite session(s) with the settings of {study.name}.",
        "",
        f"prestus_root = {_q(prestus_root)};",
        "addpath(genpath(fullfile(prestus_root, 'functions')));",
        "",
        f"parameters = load_parameters({_q(study.name)}, {_q(study.parent)});",
        f"selection = {{{rows}}};",
        "T = neuronav_ingest_batch(parameters, selection, 'progress', ...",
        "    @(k, n, sub, ses, status) fprintf('[%d/%d] %s %s: %s\\n', k, n, sub, ses, status));",
        "fprintf('%d target(s) written. Summary: %s\\n', nnz(~cellfun(@isempty, T.target)), ...",
        "    fullfile(parameters.path.localite_post, 'neuronav_ingest_summary.csv'));",
        "",
    ])


def ingest_script_path(study: Path) -> Path:
    return Path(study).with_name(f"ingest_{_stem(Path(study))}.m")


def write_ingest_script(study: Path, prestus_root: Path, selection: list[tuple[str, str]]) -> Path:
    path = ingest_script_path(study)
    write_text_atomic(path, ingest_script(study, prestus_root, selection))
    return path


def kwave_drive_script(fit_mat: Path, prestus_root: Path, library_path: str | Path,
                       sim_config: str | Path | None, axisymmetric: bool, near_field_skip_mm: float,
                       grid_dims: str | list | None = None, apply_correction: bool = True) -> str:
    """Script for the free-water k-Wave drive check of an unsaved per-depth fit (written by
    nexus_calibration 'write_fit'): kwave_calibrate_drive, then the corrected entry is saved
    to library_path with calibration_save_entry, as calibration_gui does."""
    fit_mat = Path(fit_mat).resolve()
    ko = [f"'Axisymmetric', {str(bool(axisymmetric)).lower()}",
          f"'NearFieldSkipMm', {float(near_field_skip_mm):g}"]
    if sim_config:
        ko.append(f"'SimConfig', {_q(sim_config)}")
    if isinstance(grid_dims, (list, tuple)) and grid_dims:
        ko.append(f"'GridDims', [{' '.join(f'{int(x)}' for x in grid_dims)}]")
    elif grid_dims:
        ko.append(f"'GridDims', {_q(grid_dims)}")
    if not apply_correction:
        ko.append("'ApplyCorrection', false")
    grid = "axisymmetric" if axisymmetric else "3D"
    return "\n".join([
        f"% Generated by the PRESTUS Nexus on {_dt.date.today().isoformat()}.",
        "% Free-water k-Wave drive check of a per-depth calibration fit (kwave_calibrate_drive):",
        "% one simulation per refined depth (minutes each; SLURM/qsub if the sim config or",
        "% the machine provides it), then the corrected entry is saved to the library folder.",
        "",
        f"prestus_root = {_q(prestus_root)};",
        "addpath(genpath(fullfile(prestus_root, 'functions')));",
        "addpath(genpath(fullfile(prestus_root, 'external')));",
        "",
        f"s = load({_q(fit_mat)});",
        "fit = s.fit;",
        f"ko = struct({', '.join(ko)});",
        "[lib, report] = kwave_calibrate_drive(fit.entry, fit.prof, ko);",
        "disp(report);",
        f"kw = struct('grid', {_q(grid)}, 'sim_config', {_q(sim_config or '')});",
        f"file = calibration_save_entry(lib, fit, {_q(library_path)}, struct('Kwave', kw));",
        "fprintf('Saved %s\\n', file);",
        "",
    ])


def write_kwave_drive_script(path: Path, *args, **kwargs) -> Path:
    path = Path(path)
    write_text_atomic(path, kwave_drive_script(*args, **kwargs))
    return path


def kwave_drive_submit_body(script: str | Path, fit_mat: str | Path, sim_config: str | Path | None,
                            prestus_root: str | Path) -> str:
    """MATLAB statements that submit a drive-check script as a scheduler job
    (kwave_drive_submit: a CPU job that runs the script, which submits one simulation per
    depth and waits for them). Paths as the machine running them sees them."""
    return (f"addpath(genpath(fullfile({_q(prestus_root)}, 'functions'))); "
            f"addpath(genpath(fullfile({_q(prestus_root)}, 'external'))); "
            f"kwave_drive_submit({_q(script)}, {_q(fit_mat)}, {_q(sim_config or '')})")


def write_kwave_drive_submit(script: Path, fit_mat: Path, sim_config: str | Path | None,
                             prestus_root: Path) -> Path:
    """submit_<script>.m next to the drive-check script, for MATLAB on this computer."""
    script = Path(script).resolve()
    path = script.with_name(f"submit_{script.stem}.m")
    body = kwave_drive_submit_body(script, Path(fit_mat).resolve(), sim_config, prestus_root)
    write_text_atomic(path, "\n".join([
        f"% Generated by the PRESTUS Nexus on {_dt.date.today().isoformat()}.",
        f"% Submits {script.name} as a scheduler job (kwave_drive_submit).",
        "", body.replace("; ", ";\n") + ";", ""]))
    return path
