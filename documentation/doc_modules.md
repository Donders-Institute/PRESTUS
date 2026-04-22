# Modules

`prestus_pipeline` is composed of discrete, independently toggleable modules. Each module corresponds to a stage of the TUS simulation workflow and can be enabled or disabled via `modules.*` flags in the config. This allows partial runs — for example, re-running only analysis after changing output thresholds, or stopping after segmentation before committing to a full simulation.

---

## Pipeline stages

```
load_parameters
      │
      ▼
[segmentation]          SimNIBS charm — structural MRI → tissue labels
      │  modules.segmentation_only = 1 → stop here
      ▼
[grid setup]            modules.run_grid_setup
      │                 Crop & orient head, position transducer & target
      ▼
[medium setup]          modules.run_medium_setup
      │                 Map tissue labels → acoustic/thermal properties
      ▼
[source setup]          modules.run_source_setup
      │                 Construct k-Wave source from transducer geometry
      ▼
[amplitude calibration] modules.run_amplitude_calibration  (default: off)
      │                 Fast water sim → scale elem_amp to target pressure
      ▼
[acoustic simulation]   modules.run_acoustic_sims
      │                 k-Wave pressure field simulation
      ▼
[acoustic analysis]     modules.run_acoustic_analysis
      │                 Extract Isppa, Ispta, focal metrics; write CSV
      ▼
[thermal simulation]    modules.run_heating_sims
      │                 k-Wave heat diffusion; requires acoustic results
      ▼
[thermal analysis]      (run automatically when heating results available)
      │
      ▼
[NIfTI export]          modules.run_nifti_creation
      │                 Write pressure/temperature maps in subject space
      ▼
[report]                modules.generate_report
      │                 Self-contained HTML report
      ▼
[free-water reference]  modules.run_posthoc_water_sims
                        Re-run acoustics in water for Isppa normalisation
```

---

## Common partial-run patterns

| Goal | Settings |
|---|---|
| Run SimNIBS segmentation only | `modules.segmentation_only = 1` |
| Re-run from acoustic simulation onwards | `modules.run_grid_setup = 0`, `run_medium_setup = 0`, `run_source_setup = 0` |
| Add thermal simulation to completed acoustic run | `modules.run_heating_sims = 1`; all others can remain `1` |
| Regenerate analysis outputs only | `run_acoustic_sims = 0`, `run_heating_sims = 0` |
| Skip free-water reference | `modules.run_posthoc_water_sims = 0` |
| Scale to a target free-water ISPPA | `modules.run_amplitude_calibration = 1`, `calibration.target_isppa_wcm2: 30` |

> **Note:** Downstream modules depend on outputs from upstream ones. If intermediate files already exist on disk (e.g. from a prior run), PRESTUS reuses them subject to `io.overwrite_files`. If they do not exist, disabling an upstream module while enabling a downstream one will error.

---
