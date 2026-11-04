# How PRESTUS Works

!!! note "At a glance"
    - **What:** How a run is put together: config assembly, where it executes, what each stage reads, writes and caches, and the conventions that keep results consistent.
    - **When you need it:** When you debug a run or change the pipeline. For the stage switches see [Modules](doc_modules.md); for every parameter see [Parameters](doc_parameters.md).
    - **Before you start:** Nothing; this is a reference page for contributors and advanced users.

## Overview

```
config YAML(s) ──► load_parameters ──► prestus_pipeline_start ──► prestus_pipeline ──► outputs
                    (merge, validate)    (where to run)             (stages 1-12)        CSV, NIfTI, HTML,
                                                                                         cache/, log/
```

Entry points:

| Function | Use |
|---|---|
| `prestus_pipeline_start(parameters)` | One subject: runs locally or submits an HPC job. |
| `prestus_study_start(config)` | Many subjects/conditions from a study manifest; one `prestus_pipeline_start` per row ([Group analysis](doc_group.md#study-setup)). |
| `prestus_group_report_start(config)` | Group maps and report from finished runs. |
| `prestus_nexus` | Desktop app that writes the study config and starts the above ([Nexus](doc_nexus.md)). |
| `prestus_pipeline(parameters)` | The stages themselves (called by the above; call it directly only for debugging). |

All MATLAB code is under `functions/<subsystem>/`; [Functions](doc_functions.md) lists every function by subsystem.

## Configuration

The user-facing description is [Configuring a study](doc_configuration.md). `load_parameters` builds one `parameters` struct from up to four layers, each merged onto the previous one (`mergestruct` + `restore_struct_arrays`, so struct arrays such as `transducer(2)` keep their own fields):

| Layer | Source | Notes |
|---|---|---|
| 1. Toolbox defaults | `config/config_default.yaml` | Always loaded first. |
| 2. Project defaults | `config_default.yaml` next to the study config (or in the folder passed as second argument, or the current folder) | Copied by `prestus_config_init` to pin defaults. Keys added to the toolbox since the copy was made fall back to layer 1, with the warning `prestus:staleProjectDefault`; `prestus_config_upgrade` refreshes the copy. Using the current folder's copy warns (`prestus:cwdDefault`). |
| 3. Study config(s) | the YAML file(s) or struct passed in | YAML keys PRESTUS does not know warn (`prestus:unknownConfigKey`, with a suggestion for typos). |
| 4. Overrides | trailing `{'dotted.key', value; ...}` cell | e.g. one row of a study manifest (`apply_config_overrides`). |

??? note "Details: known keys and post-merge steps"
    A key is *known* when the merged defaults contain it, or when it is listed under `extra:` or `deprecated:` in `nexus/prestus_nexus/overlay.yaml`: keys the code reads that `config_default.yaml` omits. That list is shared by MATLAB (`config_known_extra_keys`) and the Nexus, and the Nexus drift test checks that every entry is still read by `functions/`.

    After merging, `load_parameters` migrates deprecated keys, derives transducer parameters (`load_transducer_parameters`, which also resolves `transducer.serial` against the equipment library), validates the timing protocol, sanitises `io.output_affix` into a BIDS `_desc-<label>` entity, and validates backends. Missing tool paths (`hpc.ld_library_path`, `startup.simnibs_bin_path`) prompt only in an interactive desktop session with `simulation.interactive: 1`; otherwise they warn, so batch and HPC runs never wait on a dialog.

## Where a run executes

`prestus_pipeline_start` resolves `platform` (`auto` detects SLURM or qsub via `hpc_detect_system`, and errors when neither is on the PATH):

- **`matlab`**: runs in the current session.
- **`slurm` / `qsub`**: saves `parameters` to a temporary `.mat`, writes a small `.m` launcher that adds the PRESTUS paths and calls `prestus_pipeline`, and submits it (`hpc_submit_job`; SLURM via `sbatch --parsable`). With `hpc.wait_for_job` it blocks until the job ends (`hpc_wait_for_completion`); otherwise it submits a dependent report job (`hpc_job_report`: run time, memory, energy).

Three **chained pipelines** split one run into several jobs that pass results through caches:

| Mode | Triggered by | Function | Stages |
|---|---|---|---|
| Async multi-transducer | `simulation.transducer_coupling: async` and >1 transducer | `async_transducer_pipeline` | shared preprocessing → one acoustic job per transducer → incoherent combination → thermal |
| Multi-intensity | `transducer.target_isppa_wcm2` is a vector | `multi_isppa_pipeline` | one acoustic job with free-water baseline → one thermal job per intensity (field rescaled) → summary |
| Uncertainty | `simulation.uncertainty: 1` | `uncertainty_pipeline` | shared preprocessing → one job per medium-property variant → uncertainty report |

`resolve_chained_overwrite` maps `io.overwrite_files` onto these: `always` re-runs the output stages but resumes the intermediate ones from cache; `force` re-runs everything.

## Pipeline stages

`prestus_pipeline` calls one function per stage and passes the stage outputs explicitly; no stage reads another stage's variables from a global. Positions are converted between coordinate systems at fixed points only: RAS mm → T1 voxels once in the placement stage (`ras_to_grid`), T1 voxels → grid voxels in the grid stage, and back for NIfTI export. See [Coordinate systems](doc_coordinate_systems.md) before changing any transform.

??? note "Details: data passed between stages"
    | Variable | Contents | Produced by |
    |---|---|---|
    | `parameters` | Configuration plus derived values: `grid.dims`, transducer positions in grid voxels, `io.dir_*` output folders, `state.*` flags | `load_parameters`, then updated by most stages |
    | `planimg` | Planning image (T1) header and the transforms between T1 voxels and the simulation grid | grid stage |
    | `segmentation`, `medium_masks` | Tissue labels on the grid, and their mapping to `medium_properties` layers | grid stage |
    | `kwave_medium` | Sound speed, density, attenuation (and thermal properties in `medium_plus`) on the grid | medium stage |
    | `kgrid`, `source`, `sensor` | k-Wave grid, source and sensor | source stage |
    | `acoustic_provenance` | Free-water baseline ISPPA and its origin, calibration, FOV offsets | baseline + acoustic stage |
    | `sensor_data` | `p_max_all`, `p_min_all`, `p_final` (scaled to the target ISPPA) | acoustic stage |
    | `fov_state` | Full-grid copies when the acoustic run uses a cropped field of view | FOV crop |

    `parameters.state` records what is available downstream: `acoustics_available`, `heating_available`, `acoustic_dry_run` (a backend wrote its inputs but did not simulate), and `planned_trans_t1_vox` / `planned_focus_t1_vox` (planned positions in T1 voxels, kept for the targeting analysis because the grid stage replaces `trans_pos`/`focus_pos` with grid indices).

??? note "Details: stage functions, inputs and outputs"
    | # | Stage | Function(s) | Reads | Writes / caches |
    |---|---|---|---|---|
    | 1 | Segmentation | `preproc_segmentation` (SimNIBS `charm`) | T1 (T2) | `m2m_sub-XXX/` (governed by `io.overwrite_simnibs`) |
    | 2 | Placement | `preproc_transducer_placement` → `position_transducer_*` | `placement.*`, Localite/heuristic/MNI/PlanTUS inputs | `transducer(i).trans_pos/focus_pos` in T1 voxels; placement report |
    | 3 | Grid | `grid_stage` → `grid_tissue_setup`, `autopad_grid_for_transducers`, `grid_transducer_location`, `grid_axisymmetry` | T1, segmentation, pseudo-CT | `cache/…_grid_cache…mat` (`io.save_grid_cache`), `…_planimg…mat` sidecar |
    | 4 | Medium | `medium_stage` → `medium_setup`; `nifti_medium` | masks, `medium_properties`, pseudo-CT | `cache/…_medium_cache…mat` (`io.save_medium_cache`); property NIfTIs |
    | – | FOV crop | `acoustic_fov_crop` | `grid.acoustic_fov_*` | cropped medium and positions; `fov_state` |
    | 5 | Source | `coherent_drive_calibration` (coherent, >1 transducer); `source_sensor_setup` → `source_create`; `checkStability` | transducer geometry, drive | `cache/…_kwave_source…mat` |
    | – | Free-water baseline | `acoustic_water_baseline` → `needs_water_baseline`, `water_baseline` | source, water properties | `acoustic_provenance.freefield_isppa_wcm2`, `freefield_source` |
    | 6 | Acoustic | `acoustic_stage` → `run_acoustic_backend` (`acoustic_wrapper` / `jwave_wrapper` / `babelbrain_wrapper`); `apply_isppa_scaling` | medium, source, sensor | `cache/…_results…mat` (`cache_acoustic`) |
    | 7 | Acoustic analysis | `acoustic_analysis`; `nifti_acoustic`; `targeting_analysis` | `sensor_data`, masks | `…<affix>.csv`, intensity/pressure/MI NIfTIs, `…_targeting.csv` |
    | 8 | Thermal | `thermal_stage` → `thermal_simulation` (k-Wave `kWaveDiffusion`) or `babelbrain_thermal_wrapper` | acoustic field, timing protocol | `cache/…_heating_res…mat` |
    | 9 | Thermal analysis | `thermal_analysis`; `nifti_thermal` | heating results | thermal columns in the CSV, temperature/CEM43 NIfTIs |
    | 10 | Reports | `generate_simulation_report` (+ group, multi-intensity, uncertainty) | all outputs | `…_report…html` |
    | 11 | Post-hoc water | `prestus_pipeline_start` with `simulation.medium: water` | same transducer | a water run per transducer |
    | 12 | Sequential | `sequential_pipeline` | `options.sequential_configs` | follow-up run starting from this run's end temperature |

## Caches and reuse

Every expensive stage writes a cache under `cache/` named `sub-XXX_<medium>_<stage>…<affix>.mat` (the full list is under [Output directory structure](doc_outputs.md#output-directory-structure)). Whether a cache is reused is decided by `io.overwrite_files` (`confirm_overwriting`): `never` reuses, `always`/`force` recompute, anything else asks in interactive sessions.

Caches are for re-runs of the same run. A separate run starts from an earlier run's **stage outputs** instead: the SimNIBS segmentation (`path.seg`), the pseudo-CT (`pseudoCT.json`), a placement run (`placement_run.json`, `placement.reuse`) and the T1w pressure map (`path.acoustic_input`, `LOAD_ACOUSTIC_INPUT`). See [Starting from earlier stage outputs](doc_outputs.md#starting-from-earlier-stage-outputs).

??? note "Details: cache signatures and sidecars"
    Because cache names do not encode the settings, a reused cache could silently ignore a changed setting. Grid, medium and acoustic caches therefore store the settings they depend on (`cache_signature`) and compare them on reuse (`check_cache_signature`, `io.cache_signature_check`); see [Cache reuse checks](doc_outputs.md#cache-reuse-checks). The acoustic cache also checks that its grid size matches the current grid (`cache_acoustic:staleGrid`) and re-runs the acoustic simulation otherwise. A run that reuses another run's grid cache (`io.preproc_affix`, `io.reuse_upstream_caches`) keeps its own transducer settings (drive, phases, ISPPA targets) and takes only the grid-space positions from the cache (`load_grid_cache_file`). Likewise it keeps its own grid settings (e.g. a partial acoustic FOV, `grid.acoustic_fov_*`) and takes only what the grid stage computes (`dims`, `default_dims`, `resolution_mm`, `display_orientation`). Pipeline sub-jobs that reuse their own stage-0 caches (`io.reuse_upstream_caches`) skip the signature check.

    Each stage output, except the standard SimNIBS output, has a JSON sidecar (`write_sidecar`) recording the settings it depends on, compared on reuse by `check_sidecar_signature` under the same `io.cache_signature_check`. A new reusable stage output gets a sidecar written with `write_sidecar` and a signature of the settings it depends on.

## Intensity scaling

Simulations run at the drive `elem_p_amp`; results are reported at the free-water ISPPA `transducer.target_isppa_wcm2`, by scaling pressure with `sqrt(target / baseline)`. The baseline comes from a free-water simulation of the same source (`modules.run_water_baseline: auto` runs it when needed) or from a matching library calibration. A target that cannot be met stops the run rather than being reported unscaled. See [ISPPA targeting](doc_simulations-acoustic.md#isppa-targeting).

With several transducers each target is that transducer's own free-water ISPPA. Coherent runs cannot scale the combined field per transducer, so `coherent_drive_calibration` sets each transducer's `elem_p_amp` before the source is built (baseline: the transducer alone in water, or its matching calibration) and `apply_isppa_scaling` leaves the result unscaled (`acoustic_provenance.coherent_drive`). Async runs scale each field in `combine_async_intensity`. See [Multi-transducer: ISPPA targets](doc_multitransducer.md#isppa-targets).

## Acoustic and thermal backends

`acoustic_backend` (`kwave`, `jwave`, `babelbrain`) and `thermal_backend` (`kwave`, `babelbrain`) select the solver; everything before and after the solver call is shared. `run_acoustic_backend` dispatches the acoustic call, so a fresh run and a re-run after a stale cache use the same backend and the same scaling. `simulation.code_type` selects k-Wave's MATLAB or C++ code on CPU or GPU. See [Backends](doc_backend.md) and [k-Wave](doc_kwave.md).

## Logging and provenance

`path_log_setup` creates the per-subject output folders (`io.dir_*`), starts a MATLAB diary in `log/`, writes a parameters snapshot to `cache/` and a human-readable resolved-parameters YAML to `log/`. The acoustic cache carries `acoustic_provenance` (baseline ISPPA and its source, calibration, BabelBrain version). Errors close the diary and are rethrown with their identifier; opt-in anonymous usage statistics record the run's outcome ([Telemetry](doc_telemetry.md)).

## Changing PRESTUS

- **New config key:** add it with a comment to `config/config_default.yaml` and a row to [Parameters](doc_parameters.md). If the code reads a key that `config_default.yaml` does not define, list it under `extra:` in `nexus/prestus_nexus/overlay.yaml`, so that neither MATLAB nor the Nexus reports it as unknown.
- **New function:** start the file with a help line `% NAME  Short description`; `tools/docs/gen_function_index.py` adds it to [Functions](doc_functions.md).
- **New stage or behaviour:** add a test at the right level ([Testing](doc_testing.md)), and for anything that changes numbers a physics or golden-metric check.
- Before committing, run the checks that CI runs:

```bash
python tools/docs/gen_function_index.py --check
python tools/docs/check_param_docs.py
python tools/docs/check_links.py
mkdocs build --strict -d /tmp/prestus-site
(cd nexus && python -m pytest -q -m "not matlab")
matlab -batch "cd tests; run_all_tests('unit')"
```
