## Transducer Calibration

!!! note "At a glance"
    - **What:** How the phases and drive of a simulated transducer are fitted to measurements of the real device, and how calibration relates to intensity targeting.
    - **When you need it:** Only if your device is not in the [transducer library](doc_transducer_library.md) or you calibrate at depths and intensities of your own.
    - **Before you start:** Hydrophone axial profiles of the device in free water.

To emulate transducers, PRESTUS calibrates the phase and amplitude settings of virtual transducers so that their analytically predicted focal axis profiles match measurements made with real driving system–transducer setups in free water (measured in-house via hydrophones or provided by device manufacturers). Calibrated settings are stored in a transducer library and loaded at runtime.

Transducer calibration relies on an additional config (`calibration_config`) that should be loaded as `parameters.calibration` (see the example `calibration_standalone`).

---

## Calibration methods

The transducer library stores, for each measured depth, per-ring phase offsets, gains and the drive that reproduce the hydrophone profile (per-depth refinement, the default). Between measured depths they are interpolated; beyond them, geometric steering is used. See [doc_transducer_library.md](doc_transducer_library.md) for the model, the schema and the entry points.

| Level | Entry point | Phase source | Amplitude | k-Wave needed | Typical use |
|---|---|---|---|---|---|
| **Per-depth refinement** (default) | `build_transducer_library`, `refine_transducer_library`, `calibration_gui`, `calibration_transducer` | Analytical fit per depth (`fit_transducer_model` via `calibrate_library_entry`) | Analytical, per depth | No | Library entries from hydrophone profiles; calibration at chosen depths and intensities |
| **k-Wave drive calibration** (optional) | `kwave_calibrate_drive` (`KWave` option / GUI checkbox / `kwave_drive_calibration`) | Unchanged | Calibrated to free-water k-Wave ISPPA | Yes | Simulation-consistent drive at the refined depths |
| **O'Neil phase-only fit** | `fit_preset: oneil_phase_only` / `Preset = 'oneil_phase_only'` | Per-depth fit with the former global-search settings (O'Neil, FEX-minimize, phases only) | Analytical, per depth | No | Comparison with results of the former `perform_global_search` pipeline |

Because the k-Wave drive calibration is optional, most entries stay at the "Analytical" amplitude level unless someone deliberately runs it. `load_transducer_from_library` does not assume that step happened: it resolves `tx_calibration.amp_kwave_calibrated` per requested depth and warns (`load_transducer_from_library:uncalibratedAmplitude`) whenever the driven amplitude has not been k-Wave-calibrated — see [doc_transducer_library.md](doc_transducer_library.md#optional-free-water-k-wave-drive-calibration).

### Focal depth convention (FLHM)

A "focal depth" quoted for a real transducer (e.g. a manufacturer steering table, or
`focal_depths_wrt_exit_plane` in a calibration config) is not the raw intensity-peak location —
it is the centre of the **Focal Length at Half Maximum** (FLHM), the midpoint of the contiguous
span around the peak where the axial profile stays at or above a threshold fraction of the peak.
Manufacturers differ on what "half maximum" means:

| `focal_reference` | Threshold | Convention |
|---|---|---|
| `flhm_3db` (default) | intensity ≥ peak / 2 | ITRUSST / Martin et al. (2025); Sonic Concepts CTX + TPO ("NeuroFUS") steering tables |
| `flhm_6db` | intensity ≥ peak / 4 (≡ pressure ≥ peak / 2) | Imasonic / IGT steering tables |
| `peak` | n/a — raw intensity-peak location | No FLHM centring |

Set `parameters.calibration.focal_reference` to match the convention the characterisation data
actually uses — mismatching it produces calibrated depths that look shifted by several mm/cm even
though the underlying phases/amplitudes are correct. The value is:

- Read by `compute_flhm_center` (via `focal_reference_to_frac`) wherever a calibration or
  post-hoc analysis needs to locate the FLHM centre of an axial profile.
- Echoed to the console during `calibration_setup` (`Focus from exit plane: 85.00 (FLHM -6dB)`).
- Recorded into the sidecar YAML's `meta.calibration.focal_reference` by `save_optimized_values`,
  so calibrated outputs stay traceable to the convention they were produced under.
- Used by `compute_focal_metrics` (and therefore `acoustic_analysis`) so post-hoc/in-pipeline
  focal-centre reporting stays consistent with how the transducer was calibrated. This is
  independent of the fixed ITRUSST −3 dB / −6 dB **volume** fields (`mask_3dB`/`mask_6dB`), which
  always use both thresholds regardless of `focal_reference`.

Defaults to `flhm_3db` when unset, preserving prior behaviour.

### Library build (per-depth refinement)

The library is built by `build_transducer_library.m` from the profiles in `local/profiles/`. Every measured depth becomes a refinement fitted by `fit_transducer_model`: phase offsets and gains on top of geometric steering (`set_real_phases`, path-length geometry from the ring dimensions), plus the drive that reproduces the measured peak. The global model is geometric steering with the measured drive per depth. `meta.quality` records in-sample, held-out interpolation and geometric R². The result is written to `config/calibration/{combo_name}.yaml` (`LIBRARY_PATH`).

At runtime, load calibrated parameters with:
```matlab
parameters = load_transducer_from_library(combo_name, focal_depth_ep, isppa, equip_param);
```

### One calibration core

All entry points run the same functions, so a calibration done with `calibration_transducer` is identical to one done with `refine_transducer_library` or the GUI for the same profiles, depths and settings:

```
calibration_setup ─► calibration_pipeline_start ─► calibration_transducer ─┐   (config_calibration.yaml;
refine_transducer_library ─────────────────────────────────────────────────┤    chosen depths and
                                                                            ▼    intensities)
                                   calibrate_library_entry(lib, tran, meas, options)
                                     ├─ interpolate_axial_profiles   profile at each requested depth
                                     ├─ calibration_fit_depths       one fit_transducer_model fit per depth
                                     │    (calibration_fit_options: 'rayleigh_phase_gain' | 'oneil_phase_only')
                                     ├─ library_refinement_from_fits phases, gains, drive → refinements
                                     ├─ kwave_calibrate_drive        optional: drive calibrated in free-water k-Wave
                                     └─ evaluate_library_entry       R² / focal error before and after
build_transducer_library ─► calibration_fit_depths (same fit loop and presets)
calibration_gui, Nexus Calibration page ─► calibration_fit ─► calibration_fit_depths
```

`calibrate_library_entry` reads and writes no files; the callers decide where the result goes.

**Unmeasured depths.** A requested depth within 0.5 mm of a measured one uses that measurement. Any other depth is interpolated from its two neighbouring measured profiles (aligned on their focus, `profile_align`) and then fitted, and is stored with `profile_type: interpolated`. Between calibrated depths, studies interpolate the stored offsets and gains (`library_model_at_depth`). On held-out measured depths both ways predict the profile about equally well when the measured depths are up to ~10 mm apart (median R² 0.975 each, vs 0.976 for a direct fit and 0.90 for geometric steering). Fitting the interpolated profile is more robust where the focal behaviour changes quickly. Neither bridges a large gap across a change in focal behaviour, so a warning is issued when the neighbouring measurements are more than `max_gap_mm` (default 10 mm) apart.

**Fit presets** (`calibration_fit_options`):

| Preset | Free parameters | Forward model | Loss | Optimiser | Pull to prior |
|---|---|---|---|---|---|
| `rayleigh_phase_gain` (default) | phase offset + gain per ring | Rayleigh | pseudo-Huber on intensity (δ 0.8) | `lsqnonlin`, 16 starts | `Lambda` 1e-3 |
| `oneil_phase_only` | phase offset per ring | O'Neil | (practically) least squares on intensity | FEX-minimize (or `GlobalSearch`) | none |

Any `fit_transducer_model` option can be set on top of a preset (`fit_options` in the config, `FitOptions` in the functions).

---

### Calibration at chosen depths and intensities (`calibration_transducer`)

**Script: `examples/calibration_standalone.m`** → `calibration_setup` → `calibration_pipeline_start` → `calibration_transducer`

For each equipment combination in `config_calibration.yaml`, `calibration_setup` locates the hydrophone profile CSV (equipment YAML `char_data_path`, relative to `path_input_axial`), drops depths outside the transducer's focal range and calibrates all requested depths and intensities in one call (locally, or as one HPC job per combination via `platform`).

```yaml
combinations:
  - IS_PCD15287_01001_IGT_32_ch_comb_10_ch
  - IS_PCD15473_01001_IGT_32_ch_comb_10_ch
focal_depths_wrt_exit_plane:   # [] = every measured depth
  - [40, 50]
  - [40]
desired_intensities:           # free-water ISPPA [W/cm²]; the first is the stored reference
  - [30, 60]
  - [30]
```

`calibration_transducer` can also be called directly:

```matlab
[lib, report, out] = calibration_transducer(csv_path, 'IS_PCD15287_01001_IGT_32_ch_comb_10_ch', ...
    [30 60], [40 50], parameters);
```

`profile_empirical` may be the profile CSV, a `read_axial_profiles` struct, or (former interface) a single profile with `axial_intensity` and `axial_distance_bowl`.

#### Two endpoints

1. **Library entry** `{library_path}/{combo}.yaml` (default `{path_output}/library`): the calibrated depths as refinements on top of the existing user entry, or else the shipped one. Studies use it by setting `transducer.library_path` to that folder.
2. **Calibration output folder** `path_output_profiles` (default `{path_output}/profiles`):
    - `PRESTUS_virtual_{combo}.csv`: phases and amplitude per intensity (rows) × depth (columns)
    - `{combo}-F{depth}mm-I{intensity}wpercm2.yaml`: one PRESTUS transducer block per depth × intensity, with `target_isppa_wcm2` and `tx_calibration.isppa_wcm2` set to that intensity (can be ingested with `update_transducer_library`)
    - `{combo}_settings_I{intensity}.csv`: per-depth phases, per-ring amplitudes, gains, model focus (`library_settings_table`)
    - `{combo}_calibration_profiles.csv`: the measured or interpolated profiles that were fitted
    - `{combo}_calibration_report.csv`, `{combo}_calibration.png`, and `{combo}_kwave_drive_calibration.csv` if the k-Wave drive calibration ran

All values in the output folder are resolved from the library entry (`library_model_at_depth`), so both endpoints agree exactly. The drive for other intensities scales as √(I / I_ref).

#### Configuration

| Key | Meaning |
|---|---|
| `phase_mode` | `single_depth` (default): fit each depth · `geometric`: no fit, geometric steering with the measured drive (the former `multi_depth` joint model has been removed) |
| `fit_preset` | `rayleigh_phase_gain` (default) · `oneil_phase_only` (see above) |
| `forward_model`, `fit_options` | Optional overrides of the preset |
| `opt_phase_precession` | `false` · `linear` · `monotonic` phase constraint (`PhaseConstraint`) |
| `opt_weights` | 0 = uniform; ≥ 1 Gaussian weights around the focus (`FocalWeight`) |
| `opt_limits` | [min max] fitted window, mm from the bowl apex (`FitWindowMm`, converted to the exit plane) |
| `skip_front_peak_mm` | Near field ignored in fit and metrics (`NearFieldSkipMm`) |
| `profile_align`, `max_extrapolation_mm`, `max_gap_mm` | Profile interpolation at unmeasured depths |
| `ring_profile` | `auto` (default: keep the entry's) · `fit`: within-ring velocity profile shared by all depths (`fit_ring_profile`) · `none`: uniform rings; `ring_profile_n_sub` (4), `ring_profile_smooth` (0.3), `ring_profile_near_field_skip_mm` (4); `ring_profile_offaxis` (plane-scan file, `''` = axial only) and `ring_profile_offaxis_kind` (`all`). See [Within-ring velocity profile](doc_transducer_library.md#optional-within-ring-velocity-profile) |
| `add_FDO` | 1 (default): profile distances are from the exit plane · 0: from the bowl apex |
| `kwave_drive_calibration` | Free-water k-Wave simulation per depth that rescales the stored drive to the simulated ISPPA (default true; former name `run_free_water_sim`); `force_kwavearray`, `axisymmetric2D` set up its grid |
| `opt_method`, `opt_seed` | Optimiser and seed of the `oneil_phase_only` preset |

Former per-element corrections in the equipment YAML (`elem_phase_correction`) are used as the prior the fits are regularised towards. `save_elem_correction` and `opt_amp_validation` are obsolete: offsets are written to the library entry, and since the drive scales exactly with √ISPPA, one k-Wave drive calibration per depth replaces the validation run per intensity. `opt_upper_velocity` is no longer used.

---

## Library YAML format

The library YAML schema (`meta`, `global_model`, `refinement`, and the legacy forms that are still read) is documented in [doc_transducer_library.md](doc_transducer_library.md#library-yaml-schema).

---

## Calibration and ISPPA targeting

Calibration and intensity targeting are separate steps:

| | Calibration (`calibration_transducer`, library) | ISPPA targeting |
|---|---|---|
| **When** | Once per transducer/depth, before any study | Every run, after the acoustic simulation |
| **What is set** | Phases **and** drive per focal depth | Neither — the simulated field is scaled |
| **Input** | Hydrophone axial profiles | `transducer.target_isppa_wcm2` (free-water ISPPA, W/cm²; default 30) |
| **Output** | Library entry (YAML) used by studies | Pressure field scaled by √(target / free-water baseline) |

The acoustic simulation runs at the drive `elem_p_amp`. A short free-water baseline of the same source (`modules.run_water_baseline: auto`, the default) gives the free-water ISPPA of that drive, and the field is scaled to the target, so changing the target never re-runs the head simulation. A target that cannot be met stops the run; it is never reported unscaled. [ISPPA targeting](doc_simulations-acoustic.md#isppa-targeting) describes the rules, the provenance fields and the error messages; [Multi-intensity](doc_multiintensity.md) runs several targets from one acoustic simulation.

### Relationship to calibration

The typical workflow is:

1. **Once per transducer (library build)**: Run `build_transducer_library.m` to fit every measured depth. Results are stored in `config/calibration/{combo_name}.yaml`.
2. **Optionally, at chosen depths**: Run `calibration_transducer` (or `refine_transducer_library`) with your own profiles or at unmeasured depths, writing a user library entry (`transducer.library_path`). Calibrated depths are used exactly there and blended into the rest of the entry within `blend_mm`.
3. **Per subject simulation**: set `transducer.serial`, `combo.ds_serial` and `focal_distance_ep` (PRESTUS loads the calibrated parameters via `load_transducer_from_library`) and `transducer.target_isppa_wcm2` to the desired free-water intensity. The water baseline and scaling handle the rest — no need to re-run the acoustic simulation when changing the target intensity.
