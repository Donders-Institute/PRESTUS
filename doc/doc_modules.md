# Pipeline Modules

A run segments the head, places the transducer, and simulates acoustics and heating; segmentation alone can take an hour or more. Each stage of `prestus_pipeline` is switched on or off with a `modules.*` flag in the config. This allows partial runs: re-running only the analysis after changing output thresholds, or stopping after segmentation before committing to a full simulation.

On this page:

- [Pipeline Stages](#pipeline-stages): the stage order and the flag of each stage.
- [Partial Runs](#partial-runs): the settings for the common partial runs.

!!! note "At a glance"
    - **When you need it:** To run only part of the pipeline, e.g. segmentation first, or heating on top of an earlier acoustic run.
    - **Before you start:** Nothing.

## Pipeline Stages

`prestus_pipeline_start` decides where a run executes (local [MATLAB](https://www.mathworks.com/products/matlab.html) or an HPC job) and hands uncertainty, pose-workflow, async multi-transducer, and multi-[I<sub>SPPA</sub>](doc_glossary.md#isppa) runs to their own pipelines (pose workflow: [Time Reversal](doc_time_reversal.md#pose-workflow)). `prestus_pipeline` then runs the stages below. [How PRESTUS works](doc_architecture.md) describes what each stage reads, writes, and caches.

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
[segmentation]          SimNIBS charm (layered media only)
      │                 modules.segmentation_only = 1 → stop here
      ▼
[placement]             modules.run_transducer_placement
      │                 placement.mode → trans_pos / focus_pos in T1 voxels
      ▼
[grid]                  modules.run_grid_setup              (grid cache)
      │                 reorient, crop, smooth the head; positions → grid voxels
      ▼
[medium]                modules.run_medium_setup            (medium cache)
      │                 tissue labels / pseudo-CT → acoustic & thermal properties
      │                 NIfTI export of masks and property maps
      ▼
[acoustic FOV crop]     grid.acoustic_fov_*  (optional sub-grid for the acoustic run)
      ▼
[source]                modules.run_source_setup
      │                 kWaveArray / mask source, sensor, CFL check
      ▼
[free-water baseline]   modules.run_water_baseline (auto | 1 | 0)
      │                 free-water ISPPA of the current drive, for target_isppa_wcm2 scaling
      ▼
[time reversal]         time_reversal.enabled  (refocusing phases, aberration report, pose search; skipped with BabelBrain)
      │                 time_reversal.only = 1 → stop here (no water baseline either)
      ▼
[acoustic]              modules.run_acoustic_sims           (acoustic cache)
      │                 k-Wave / jwave / BabelBrain, or a prior pressure map (path.acoustic_input);
      │                 field scaled to target_isppa_wcm2; dry-run backends stop here
      ▼
[acoustic analysis]     modules.run_acoustic_analysis
      │                 Isppa, MI, focal metrics → CSV; NIfTI (modules.run_nifti_creation);
      │                 targeting accuracy (modules.run_targeting_analysis)
      │                 BabelBrain input export (io.export_babelbrain)
      ▼
[thermal]               modules.run_heating_sims
      │                 k-Wave diffusion or BabelBrain BHTE over the timing protocol
      ▼
[thermal analysis]      modules.run_thermal_analysis
      │                 peak temperature, CEM43 → CSV; NIfTI
      ▼
[reports]               modules.generate_report / generate_group_report /
      │                 multi_isppa_report / uncertainty_report
      ▼
[post-hoc water]        modules.run_posthoc_water_sims  (one water run per transducer)
      ▼
[sequential]            options.sequential_configs  (follow-up run on the same medium)
```

## Partial Runs

| To | Set |
|---|---|
| Segment first, simulate later | `modules.segmentation_only: 1`: the run stops after [SimNIBS](https://simnibs.github.io/simnibs/build/html/index.html) segmentation. |
| Reuse an earlier heuristic, [PlanTUS](https://github.com/mlueckel/PlanTUS), or [TUS_entry](https://github.com/CyrilAtkinson/TUS_entry) placement | Nothing: `placement.reuse: auto` reuses a recorded run whose settings match; `always` requires one. See [Starting from earlier stage outputs](doc_outputs_reuse.md#starting-from-earlier-stage-outputs). |
| Simulate heating on top of an earlier acoustic run | `path.acoustic_input: <earlier run>/nii/sub-XXX_<medium>_T1w<affix>_pressure.nii.gz`, with the other modules left on (grid and medium take seconds to a minute). The map is rescaled to `transducer.target_isppa_wcm2` from its sidecar. |
| Redo the acoustic analysis only | `modules.run_acoustic_sims: 0` and `modules.run_heating_sims: 0`. The acoustic stage then loads `cache/…_results….mat`, which the earlier run wrote only with `io.save_acoustic_matrices: 1`; without it, use `path.acoustic_input` as above. |
| Skip the post-hoc free-water runs | `modules.run_posthoc_water_sims: 0` |

For example, heating on top of an earlier acoustic run (`modules.run_heating_sims` is off by default):

```yaml
path:
  acoustic_input: /project/sim/sub-001/nii/sub-001_layered_T1w_pressure.nii.gz
modules:
  run_heating_sims: 1
```

!!! note
    Downstream modules depend on outputs from upstream ones. Within-run caches in `cache/` are reused subject to `io.overwrite_files`, and a reused grid or medium cache must match the current settings ([Cache reuse checks](doc_outputs_reuse.md#cache-reuse-checks)). The grid, medium, and acoustic caches are only written when `io.save_grid_cache`, `io.save_medium_cache`, and `io.save_acoustic_matrices` (or `io.save_matrices`) are on; by default they are not. If a stage is switched off and its cache does not exist, the stage does not stop the run: the grid and medium stages continue with empty results, and the acoustic and thermal stages report that no results are available and skip their analysis. To start a run from the results of an earlier one, use its stage outputs instead ([Starting from earlier stage outputs](doc_outputs_reuse.md#starting-from-earlier-stage-outputs)).
