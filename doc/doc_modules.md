# Modules

`prestus_pipeline` is composed of discrete, independently toggleable modules. Each module corresponds to a stage of the TUS simulation workflow and can be enabled or disabled via `modules.*` flags in the config. This allows partial runs — for example, re-running only analysis after changing output thresholds, or stopping after segmentation before committing to a full simulation.

---

## Pipeline stages

`prestus_pipeline_start` decides where a run executes (local MATLAB or an HPC job) and hands async multi-transducer runs to `async_transducer_pipeline`. `prestus_pipeline` then runs the stages below; each is a function call in `functions/prestus_pipeline.m`, and [How PRESTUS works](doc_architecture.md) describes what each stage reads, writes and caches.

```
load_parameters                 toolbox default → project default → study config (+ overrides)
      │
prestus_pipeline_start          platform: matlab | slurm | qsub;  async coupling → async_transducer_pipeline
      │
prestus_pipeline
      │  special modes (run instead of the stages below):
      │    simulation.uncertainty = 1          → uncertainty_pipeline
      │    target_isppa_wcm2 is a vector       → multi_isppa_pipeline
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
[6 acoustic]            modules.run_acoustic_sims           (acoustic cache)
      │                 k-Wave / jwave / BabelBrain, or io.external_acoustic_nifti;
      │                 field scaled to target_isppa_wcm2
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

## Common partial-run patterns

| Goal | Settings |
|---|---|
| Run SimNIBS segmentation only | `modules.segmentation_only = 1` |
| Re-run from acoustic simulation onwards | `modules.run_grid_setup = 0`, `run_medium_setup = 0`, `run_source_setup = 0` |
| Add thermal simulation to completed acoustic run | `modules.run_heating_sims = 1`; all others can remain `1` |
| Regenerate analysis outputs only | `run_acoustic_sims = 0`, `run_heating_sims = 0` |
| Skip free-water reference | `modules.run_posthoc_water_sims = 0` |
| Run thermal at a specific free-water ISPPA | `transducer.target_isppa_wcm2: 30` (with `modules.run_water_baseline: auto`, the free-water baseline runs when the target needs it); see [ISPPA targeting](doc_simulations-acoustic.md#isppa-targeting) |
| Report the field at the simulated drive (`elem_p_amp`) | `transducer.target_isppa_wcm2: []` |
| Run thermal at multiple ISPPAs (parallel jobs) | `transducer.target_isppa_wcm2: [10, 20, 30, 50]` — triggers multi-ISPPA mode; see [doc_multiintensity.md](doc_multiintensity.md) |
| Simulate two or more transducers asynchronously | `simulation.transducer_coupling: async` with multiple `transducer` entries; see [doc_multitransducer.md](doc_multitransducer.md) |
| Always record the free-water baseline | `modules.run_water_baseline: 1` (`auto`, the default, runs it only when a target needs it) |
| Generate a cross-subject group report | `modules.generate_group_report = 1` (or run `prestus_group_report_start` manually); see [doc_group.md](doc_group.md#group-html-report) |
| Run many subjects with per-subject settings | `study.manifest` + `prestus_study_start(config)`; see [doc_group.md](doc_group.md#study-setup) |

> **Note:** Downstream modules depend on outputs from upstream ones. If intermediate files already exist on disk (e.g. from a prior run), PRESTUS reuses them subject to `io.overwrite_files`; a reused grid or medium cache must match the current settings ([Cache reuse checks](doc_outputs.md#cache-reuse-checks)). If they do not exist, disabling an upstream module while enabling a downstream one will error.

---
