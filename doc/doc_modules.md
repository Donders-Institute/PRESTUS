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
[NIfTI export — medium] modules.run_nifti_creation
      │                 Write medium_masks + property maps (sound_speed,
      │                 density, alpha_coeff) in subject and MNI space
      ▼
[source setup]          modules.run_source_setup
      │                 Construct k-Wave source from transducer geometry
      ▼
[water baseline]        modules.run_water_baseline  (default: off)
      │                 Fast water sim → record free-water ISPPA provenance
      │                 alongside the acoustic cache; no elem_amp changes
      ▼
[acoustic simulation]   modules.run_acoustic_sims
      │                 k-Wave pressure field simulation; cache saved with
      │                 acoustic_provenance.freefield_isppa_wcm2
      ▼
[ISPPA scaling]         Applied automatically on cache load when
      │                 transducer.target_isppa_wcm2 is set;
      │                 scales p_max_all by sqrt(target/baseline)
      ▼
[acoustic analysis]     modules.run_acoustic_analysis
      │                 Extract Isppa, Ispta, focal metrics; write CSV
      ▼
[NIfTI export — acoustic] modules.run_nifti_creation
      │                 Write intensity, MI, pressure maps
      ▼
[thermal simulation]    modules.run_heating_sims
      │                 k-Wave heat diffusion; requires acoustic results
      ▼
[thermal analysis]      (run automatically when heating results available)
      │
      ▼
[NIfTI export — thermal] modules.run_nifti_creation
      │                 Write heating, heatrise, CEM43 maps
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
| Run thermal at a specific free-water ISPPA | `transducer.target_isppa_wcm2: 30` (water baseline runs automatically) |
| Run thermal at multiple ISPPAs (parallel jobs) | `transducer.target_isppa_wcm2: [10, 20, 30, 50]` — triggers multi-ISPPA mode; see [doc_multi_isppa.md](doc_multi_isppa.md) |
| Simulate two or more transducers asynchronously | `simulation.transducer_coupling: async` with multiple `transducer` entries; see [doc_async_transducer.md](doc_async_transducer.md) |
| Enable the water baseline measurement | `modules.run_water_baseline: 1` (required for ISPPA scaling) |

> **Note:** Downstream modules depend on outputs from upstream ones. If intermediate files already exist on disk (e.g. from a prior run), PRESTUS reuses them subject to `io.overwrite_files`. If they do not exist, disabling an upstream module while enabling a downstream one will error.

---
