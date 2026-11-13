# Modules

A run segments the head, places the transducer and simulates acoustics and heating, and segmentation alone can take an hour or more. Each stage of `prestus_pipeline` is switched on or off with a `modules.*` flag in the config. This allows partial runs: re-running only the analysis after changing output thresholds, or stopping after segmentation before committing to a full simulation.

On this page:

- [Common partial-run patterns](#common-partial-run-patterns): the settings for frequent partial runs.
- [Pipeline stages](#pipeline-stages): the stage order and the flag of each stage.

!!! note "At a glance"
    - **When you need it:** To run only part of the pipeline, e.g. segmentation first, or heating on top of an earlier acoustic run.
    - **Before you start:** Nothing.

## Common Partial-Run Patterns

| Goal | Settings |
|---|---|
| Run [SimNIBS](https://simnibs.github.io/simnibs/build/html/index.html) segmentation only | `modules.segmentation_only = 1` |
| Reuse a heuristic / [PlanTUS](https://github.com/mlueckel/PlanTUS) / [TUS_entry](https://github.com/CyrilAtkinson/TUS_entry) placement | Nothing to set: `placement.reuse: auto` reuses a recorded run whose settings match; `always` requires one. See [Starting from earlier stage outputs](doc_outputs_reuse.md#starting-from-earlier-stage-outputs) |
| Thermal simulation (or analysis) on a completed acoustic run | `path.acoustic_input: <earlier run>/nii/sub-XXX_<medium>_T1w<affix>_pressure.nii.gz`; keep the other modules on (grid and medium take seconds to a minute). The map is rescaled to `transducer.target_isppa_wcm2` from its sidecar |
| Regenerate analysis outputs from the acoustic cache | `run_acoustic_sims = 0`, `run_heating_sims = 0` (the acoustic stage then loads `cache/…_results….mat`, which the earlier run only wrote with `io.save_acoustic_matrices: 1`). Without that cache, use `path.acoustic_input` as above |
| Skip free-water reference | `modules.run_posthoc_water_sims = 0` |
| Run thermal at a specific free-water ISPPA | `transducer.target_isppa_wcm2: 30` (with `modules.run_water_baseline: auto`, the free-water baseline runs when the target needs it); see [ISPPA targeting](doc_simulations-acoustic.md#isppa-targeting) |
| Report the field at the simulated drive (`elem_p_amp`) | `transducer.target_isppa_wcm2: []` |
| Run thermal at multiple ISPPAs (parallel jobs) | `transducer.target_isppa_wcm2: [10, 20, 30, 50]` triggers multi-ISPPA mode; see [Multi-intensity sweeps](doc_multiintensity.md) |
| Simulate two or more transducers asynchronously | `simulation.transducer_coupling: async` with multiple `transducer` entries; see [Multi-transducer](doc_multitransducer.md) |
| Always record the free-water baseline | `modules.run_water_baseline: 1` (`auto`, the default, runs it only when a target needs it) |
| Generate a cross-subject group report | `modules.generate_group_report = 1` (or run `prestus_group_report_start` manually); see [Group HTML report](doc_group.md#group-html-report) |
| Run many subjects with per-subject settings | `study.manifest` + `prestus_study_start(config)`; see [Study setup](doc_group.md#study-setup) |

!!! note
    Downstream modules depend on outputs from upstream ones. Within-run caches in `cache/` are reused subject to `io.overwrite_files`, and a reused grid or medium cache must match the current settings ([Cache reuse checks](doc_outputs_reuse.md#cache-reuse-checks)). The grid, medium and acoustic caches are only written when `io.save_grid_cache`, `io.save_medium_cache` and `io.save_acoustic_matrices` (or `io.save_matrices`) are on; by default they are not. If a stage is switched off and its cache does not exist, the stage does not stop the run: the grid and medium stages continue with empty results, and the acoustic and thermal stages report that no results are available and skip their analysis. To start a run from the results of an earlier one, use its stage outputs instead ([Starting from earlier stage outputs](doc_outputs_reuse.md#starting-from-earlier-stage-outputs)).

## Pipeline Stages

`prestus_pipeline_start` decides where a run executes (local [MATLAB](https://www.mathworks.com/products/matlab.html) or an HPC job) and hands uncertainty, pose-workflow, async multi-transducer and multi-ISPPA runs to their own pipelines (pose workflow: [Time Reversal](doc_time_reversal.md#pose-workflow)). `prestus_pipeline` then runs the stages below. [How PRESTUS works](doc_architecture.md) describes what each stage reads, writes and caches.

??? note "Stage diagram"
    ```
    load_parameters                 toolbox default → project default → study config (+ overrides)
          │
    prestus_pipeline_start          platform: matlab | slurm | qsub;  special modes, in this order:
          │                           simulation.uncertainty = 1           → uncertainty_pipeline
          │                           time_reversal.pose_workflow.enabled  → pose_workflow_pipeline
          │                           async coupling                       → async_transducer_pipeline
          │                           target_isppa_wcm2 is a vector        → multi_isppa_pipeline
          │
    prestus_pipeline
          │  special modes (run instead of the stages below):
          │    simulation.uncertainty = 1          → uncertainty_pipeline
          │    target_isppa_wcm2 is a vector       → multi_isppa_pipeline
          │    time_reversal.pose_workflow.enabled → pose_workflow_pipeline
          │    modules.combine_async = 1           → combine_async_intensity
          ▼
    [1 segmentation]        SimNIBS charm (layered media only)
          │                 modules.segmentation_only = 1 → stop here
          ▼
    [2 placement]           modules.run_transducer_placement
          │                 placement.mode → trans_pos / focus_pos in T1 voxels
          ▼
    [3 grid]                modules.run_grid_setup              (grid cache)
          │                 reorient, crop, smooth the head; positions → grid voxels
          ▼
    [4 medium]              modules.run_medium_setup            (medium cache)
          │                 tissue labels / pseudo-CT → acoustic & thermal properties
          │                 NIfTI export of masks and property maps
          ▼
    [  acoustic FOV crop]   grid.acoustic_fov_*  (optional sub-grid for the acoustic run)
          ▼
    [5 source]              modules.run_source_setup
          │                 kWaveArray / mask source, sensor, CFL check
          ▼
    [  free-water baseline] modules.run_water_baseline (auto | 1 | 0)
          │                 free-water ISPPA of the current drive, for target_isppa_wcm2 scaling
          ▼
    [  time reversal]       time_reversal.enabled  (refocusing phases, aberration report, pose search; skipped with BabelBrain)
          │                 time_reversal.only = 1 → stop here (no water baseline either)
          ▼
    [6 acoustic]            modules.run_acoustic_sims           (acoustic cache)
          │                 k-Wave / jwave / BabelBrain, or a prior pressure map (path.acoustic_input);
          │                 field scaled to target_isppa_wcm2; dry-run backends stop here
          ▼
    [7 acoustic analysis]   modules.run_acoustic_analysis
          │                 Isppa, MI, focal metrics → CSV; NIfTI (modules.run_nifti_creation);
          │                 targeting accuracy (modules.run_targeting_analysis)
          │                 BabelBrain input export (io.export_babelbrain)
          ▼
    [8 thermal]             modules.run_heating_sims
          │                 k-Wave diffusion or BabelBrain BHTE over the timing protocol
          ▼
    [9 thermal analysis]    modules.run_thermal_analysis
          │                 peak temperature, CEM43 → CSV; NIfTI
          ▼
    [10 reports]            modules.generate_report / generate_group_report /
          │                 multi_isppa_report / uncertainty_report
          ▼
    [11 post-hoc water]     modules.run_posthoc_water_sims  (one water run per transducer)
          ▼
    [12 sequential]         options.sequential_configs  (follow-up run on the same medium)
    ```
