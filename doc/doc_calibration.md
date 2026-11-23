# Transducer Calibration

![Pipeline steps; this page belongs to 2 Transducer](assets/figures/fig_stage-strip-transducer_diagram.svg){ .stage-strip }

A simulation stands for a real transducer only if its virtual source produces the same focus and intensity as the device. Calibration fits the simulated phases and drive to hydrophone measurements of the real device in free water, so that a focal depth or intensity set in the simulation corresponds to the same setting in the lab.

On this page:

- [Focal Depth Convention (FLHM)](#focal-depth-convention-flhm): which point a quoted focal depth refers to.
- [Calibration Methods](#calibration-methods): fit presets and unmeasured depths.
- [Calibration and I<sub>SPPA</sub> Targeting](#calibration-and-isppa-targeting): how the calibrated drive relates to intensity scaling.
- [Calibrating at Chosen Depths and Intensities](#calibrating-at-chosen-depths-and-intensities): calibrating single conditions with `calibration_transducer`.

!!! note "At a glance"
    - **When you need it:** To know how a calibrated entry relates to a focal depth and a target intensity, or to calibrate a [library](doc_transducer_library.md) entry at depths and intensities of your own.
    - **Before you start:** Hydrophone axial profiles of the device in free water ([Input file formats](doc_calibration_profile-formats.md#input-file-formats)).

PRESTUS calibrates the phase and amplitude settings of a virtual transducer so that its analytically predicted focal axis profile matches free-water measurements of the real driving system–transducer setup (hydrophone profiles measured in-house or provided by the manufacturer). Calibrated settings are stored in the [transducer library](doc_transducer_library.md) and loaded at run time. The library stores, for each measured depth, per-ring phase offsets, gains, and the drive that reproduce the hydrophone profile (per-depth refinement, the default). Between measured depths they are interpolated; beyond them, the entry fades back to its global model (geometric steering) over 5 mm.

Subpages:

- [Emulating Hydrophone Profiles](doc_calibration_hydrophone-profiles.md): for anyone fitting a library entry to their own hydrophone profiles; also covers the within-ring velocity profile, the free-water k-Wave drive calibration, and the profile file formats.
- [Matrix Arrays: Steering Correction](doc_calibration_matrix-steering.md): only for users of matrix (2D) arrays.

## Focal Depth Convention (FLHM)

A "focal depth" quoted for a real transducer (e.g. a manufacturer steering table, or `focal_depths_wrt_exit_plane` in a calibration config) is not the raw intensity-peak location. It is the centre of the **Focal Length at Half Maximum** ([FLHM](doc_glossary.md#flhm)): the midpoint of the contiguous span around the peak where the axial profile stays at or above a threshold fraction of the peak. Manufacturers differ on what "half maximum" means:

| `focal_reference` | Threshold | Convention |
|---|---|---|
| ![](assets/icons/icon_focalref-flhm-3db_glyph.svg){ .choice-icon } `flhm_3db` (default) | intensity ≥ peak / 2 | [ITRUSST](https://itrusst.org/) / [Martin et al. (2024)](https://doi.org/10.1016/j.brs.2024.04.013); Sonic Concepts CTX + TPO ("NeuroFUS") steering tables |
| ![](assets/icons/icon_focalref-flhm-6db_glyph.svg){ .choice-icon } `flhm_6db` | intensity ≥ peak / 4 (≡ pressure ≥ peak / 2) | Imasonic / IGT steering tables |
| ![](assets/icons/icon_focalref-peak_glyph.svg){ .choice-icon } `peak` | n/a (raw intensity-peak location) | No FLHM centring |

Set `parameters.calibration.focal_reference` to the convention the characterisation data uses. A mismatch produces calibrated depths that look shifted by several mm, although the underlying phases and amplitudes are correct. Default: `flhm_3db`.

??? note "Where focal_reference is used"
    - Echoed to the console during `calibration_setup` (`Focus from exit plane: 85.00 (FLHM -6dB)`).
    - Recorded in the sidecar YAML's `meta.calibration.focal_reference` by `save_optimized_values`, so calibrated outputs stay traceable to their convention.
    - Used by `compute_focal_metrics` (and therefore `acoustic_analysis`), so focal-centre reporting in the pipeline matches the calibration. This is independent of the fixed ITRUSST −3 dB / −6 dB **volume** fields (`mask_3dB`/`mask_6dB`), which always use both thresholds.

??? developer "FLHM centre computation"
    `focal_reference` is read by `compute_flhm_center` (via `focal_reference_to_frac`) wherever a calibration or post-hoc analysis locates the FLHM centre of an axial profile.

## Calibration Methods

Per-depth refinement (analytical, no k-Wave) is the default; the k-Wave drive calibration optionally recalibrates its drive, and the [O'Neil](https://doi.org/10.1121/1.1906542) phase-only fit serves comparison with the `perform_global_search` pipeline.

??? note "Calibration levels compared"
    | Level | Entry point | Phase source | Amplitude | k-Wave needed | Typical use |
    |---|---|---|---|---|---|
    | **Per-depth refinement** (default) | `build_transducer_library`, `refine_transducer_library`, `calibration_gui`, `calibration_transducer` | Analytical fit per depth (`fit_transducer_model` via `calibrate_library_entry`) | Analytical, per depth | No | Library entries from hydrophone profiles; calibration at chosen depths and intensities |
    | **k-Wave drive calibration** (optional) | `kwave_calibrate_drive` (`KWave` option / GUI checkbox / `kwave_drive_calibration`) | Unchanged | Calibrated to free-water k-Wave I<sub>SPPA</sub> | Yes | Simulation-consistent drive at the refined depths |
    | **O'Neil phase-only fit** | `fit_preset: oneil_phase_only` / `Preset = 'oneil_phase_only'` | Per-depth fit with global-search settings (O'Neil, [FEX-minimize](https://www.mathworks.com/matlabcentral/fileexchange/24298-minimize), phases only) | Analytical, per depth | No | Comparison with results of the `perform_global_search` pipeline |

Most entries stay at the analytical amplitude level unless someone runs the k-Wave drive calibration. `load_transducer_from_library` resolves `tx_calibration.amp_kwave_calibrated` per requested depth and warns (`load_transducer_from_library:uncalibratedAmplitude`) whenever the driven amplitude has not been k-Wave-calibrated; see [free-water k-Wave drive calibration](doc_calibration_hydrophone-profiles.md#free-water-k-wave-drive-calibration).

??? note "Fit settings"
    | Setting | Behaviour |
    |---|---|
    | Free parameters (`FreeParams`) | `phase` and `amp` (per-ring gain). `z_shift` and `c_steer` only when calling `fit_transducer_model` directly. Rings on one channel share phase and gain (the 30 sub-rings of IS_PCD15287 fit as 10 channels). |
    | Loss | Pseudo-Huber on the peak-normalised profile beyond `NearFieldSkipMm`, every depth weighted equally. |
    | Optimiser | `lsqnonlin` from `NStarts` starts (geometry first, then random). Needs the [Optimization Toolbox](https://www.mathworks.com/products/optimization.html); loading entries does not. |
    | Near field | The first `NearFieldSkipMm` (10 mm) are not fitted. `NearFieldPenalty` discourages a model near field hotter than the measured one; every fit reports `nearfield_excess`. |
    | Plausibility | Gains range over 0.1–10× (`AmpBounds`). A gain outside 0.5–2× the mean (`PlausibleGain`), or at a bound, is flagged in `model.flags`; check the element, cable, or channel. |

    The static joint model (`LIBRARY_MODE = 'static'`, `Mode = 'static'`, `phase_mode: multi_depth`) raises an error. Entries that already hold a static offset still load.

??? developer "How the model is fitted and library model version"
    `fit_transducer_model` fits the on-axis field p(z) = H(z)·w, linear in the complex element drives w, to the measured profile at each depth.
    H comes from `axial_transfer_matrix` (Rayleigh by default).

    **Library model version.** Entries with `meta.library_model_version: 2` apply per-ring gains to `elem_p_amp`; version 1 applies a phase offset only. Version-3 entries (`library_mode: per_depth`) take phases, gains, and drive from the refinement at each measured depth.

??? note "One calibration core"
    All entry points run the same functions, so a calibration done with `calibration_transducer` is identical to one done with `refine_transducer_library` or the GUI for the same profiles, depths and settings.

??? developer "Calibration call graph"
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

??? developer "Library build and run-time loading"
    The shipped library is built by `build_transducer_library.m` from the profiles in `local/profiles/`: every measured depth becomes a refinement fitted by `fit_transducer_model` (phase offsets and gains on top of geometric steering from `set_real_phases`, plus the drive that reproduces the measured peak), written to `config/calibration/model-<MODEL>/tran-<TRAN>_ds-<DS>_calib.yaml` (`LIBRARY_PATH`). The entry schema is documented in the [library YAML schema](doc_transducer_library.md#library-yaml-schema). At run time, `load_transducer_from_library(combo_name, focal_depth_ep, isppa, equip_param)` loads the calibrated parameters.

### Fit Presets

Presets come from `calibration_fit_options`. Any `fit_transducer_model` option can be set on top of a preset (`fit_options` in the config, `FitOptions` in the functions).

??? note "Preset settings"
    | Preset | Free parameters | Forward model | Loss | Optimiser | Pull to prior |
    |---|---|---|---|---|---|
    | `rayleigh_phase_gain` (default) | phase offset + gain per ring | Rayleigh | pseudo-Huber on intensity (δ 0.8) | `lsqnonlin`, 16 starts | `Lambda` 1e-3 |
    | `oneil_phase_only` | phase offset per ring | O'Neil | (practically) least squares on intensity | FEX-minimize (or `GlobalSearch`) | none |

### Unmeasured Depths

A requested depth within 0.5 mm of a measured one uses that measurement.

??? note "Interpolated profiles, the gap warning, and interpolating fits instead"
    Any other depth is interpolated from its two neighbouring measured profiles (aligned on their focus, `profile_align`) and then fitted, and is stored with `profile_type: interpolated`. A warning is issued when the neighbouring measurements are more than `max_gap_mm` (default 10 mm) apart.

    Between calibrated depths, studies interpolate the stored offsets and gains (`library_model_at_depth`). On held-out measured depths both ways predict the profile about equally well when the measured depths are up to ~10 mm apart (median R² 0.975 each, vs 0.976 for a direct fit and 0.90 for geometric steering). Fitting the interpolated profile is more robust where the focal behaviour changes quickly. Neither bridges a large gap across a change in focal behaviour.

## Calibration and I<sub>SPPA</sub> Targeting

Calibration and intensity targeting are separate steps:

| | Calibration (`calibration_transducer`, library) | I<sub>SPPA</sub> targeting |
|---|---|---|
| **When** | Once per transducer/depth, before any study | Every run, after the acoustic simulation |
| **What is set** | Phases **and** drive per focal depth | Neither; the simulated field is scaled |
| **Input** | Hydrophone axial profiles | `transducer.target_isppa_wcm2` (free-water I<sub>SPPA</sub>, W/cm²; default 30) |
| **Output** | Library entry (YAML) used by studies | Pressure field scaled by √(target / free-water baseline) |

The acoustic simulation runs at the drive `elem_p_amp`. A short [free-water baseline](doc_glossary.md#free-water-baseline) of the same source (`modules.run_water_baseline: auto`, the default) gives the free-water [I<sub>SPPA</sub>](doc_glossary.md#isppa) of that drive, and the field is scaled to the target, so changing the target never re-runs the head simulation. A target that cannot be met stops the run; it is never reported unscaled. [I<sub>SPPA</sub> targeting](doc_simulations-acoustic.md#isppa-targeting) describes the rules, the provenance fields, and the error messages; [Multi-intensity](doc_multiintensity.md) runs several targets from one acoustic simulation.

The typical workflow:

1. **Once per transducer (library build)**: run `build_transducer_library.m` to fit every measured depth. Results are stored in `config/calibration/` ([layout](doc_transducer_library.md#using-a-library-entry)).
2. **Optionally, at chosen depths**: run `calibration_transducer` (or `refine_transducer_library`) with your own profiles or at unmeasured depths, writing a user library entry (`transducer.library_path`). Calibrated depths are used exactly there and blended into the rest of the entry within `blend_mm`.
3. **Per subject simulation**: set `transducer.serial`, `combo.ds_serial`, and `focal_distance_ep` (PRESTUS loads the calibrated parameters via `load_transducer_from_library`) and `transducer.target_isppa_wcm2` to the desired free-water intensity. The water baseline and scaling handle the rest; the acoustic simulation is not re-run when the target intensity changes.

## Calibrating at Chosen Depths and Intensities

**Script: `examples/calibration_standalone.m`** → `calibration_setup` → `calibration_pipeline_start` → `calibration_transducer`

Calibration reads its own config (`config_calibration.yaml`), loaded as `parameters.calibration` (as in the example script). For each equipment combination listed there, `calibration_setup` locates the hydrophone profile CSV (equipment YAML `char_data_path`, relative to `path_input_axial`), drops depths outside the transducer's focal range and calibrates all requested depths and intensities in one call (locally, or as one HPC job per combination via `platform`).

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

`profile_empirical` may be the profile CSV, a `read_axial_profiles` struct, or a single profile with `axial_intensity` and `axial_distance_bowl`.

### Calibration Outputs

Calibration writes a library entry `{library_path}/{combo}.yaml`, which studies use through `transducer.library_path`, and a calibration output folder `path_output_profiles`.

??? note "Files written to each endpoint"
    1. **Library entry** `{library_path}/{combo}.yaml` (default `{path_output}/library`): the calibrated depths as refinements on top of the existing user entry, or else the shipped one. Studies use it by setting `transducer.library_path` to that folder.
    2. **Calibration output folder** `path_output_profiles` (default `{path_output}/profiles`):
        - `PRESTUS_virtual_{combo}.csv`: phases and amplitude per intensity (rows) × depth (columns)
        - `{combo}-F{depth}mm-I{intensity}wpercm2.yaml`: one PRESTUS transducer block per depth × intensity, with `target_isppa_wcm2` and `tx_calibration.isppa_wcm2` set to that intensity (can be ingested with `update_transducer_library`)
        - `{combo}_settings_I{intensity}.csv`: per-depth phases, per-ring amplitudes, gains, model focus (`library_settings_table`)
        - `{combo}_calibration_profiles.csv`: the measured or interpolated profiles that were fitted
        - `{combo}_calibration_report.csv`, `{combo}_calibration.png`, and `{combo}_kwave_drive_calibration.csv` if the k-Wave drive calibration ran

All values in the output folder are resolved from the library entry (`library_model_at_depth`), so both endpoints agree exactly. The drive for other intensities scales as √(I / I_ref).

### Calibration Configuration

Keys of `config_calibration.yaml` beyond the combinations, depths, and intensities above.

??? note "Calibration config keys"
    All keys of `config_calibration.yaml` are listed in [Calibration Parameters](doc_parameters_calibration.md).

    | Key | Meaning |
    |---|---|
    | `phase_mode` | `single_depth` (default): fit each depth · `geometric`: no fit, geometric steering with the measured drive (`multi_depth` raises an error) |
    | `fit_preset` | `rayleigh_phase_gain` (default) · `oneil_phase_only` (see [Fit presets](#fit-presets)) |
    | `forward_model`, `fit_options` | Optional overrides of the preset |
    | `opt_phase_precession` | `false` · `linear` · `monotonic` phase constraint (`PhaseConstraint`) |
    | `opt_weights` | 0 = uniform; ≥ 1 Gaussian weights around the focus (`FocalWeight`) |
    | `opt_limits` | [min max] fitted window, mm from the bowl apex (`FitWindowMm`, converted to the exit plane) |
    | `skip_front_peak_mm` | Near field ignored in fit and metrics (`NearFieldSkipMm`) |
    | `profile_align`, `max_extrapolation_mm`, `max_gap_mm` | Profile interpolation at unmeasured depths |
    | `ring_profile` | `auto` (default: keep the entry's) · `fit`: within-ring velocity profile shared by all depths (`fit_ring_profile`) · `none`: uniform rings; `ring_profile_n_sub` (4), `ring_profile_smooth` (0.3), `ring_profile_near_field_skip_mm` (4); `ring_profile_offaxis` (plane-scan file, `''` = axial only) and `ring_profile_offaxis_kind` (`all`). See [Within-ring velocity profile](doc_calibration_hydrophone-profiles.md#within-ring-velocity-profile) |
    | `add_FDO` | 1 (default): profile distances are from the exit plane · 0: from the bowl apex |
    | `focal_reference` | FLHM convention of the profile headers; see [Focal depth convention](#focal-depth-convention-flhm) |
    | `kwave_drive_calibration` | Free-water k-Wave simulation per depth that rescales the stored drive to the simulated I<sub>SPPA</sub> (default true; `run_free_water_sim` is still read); `force_kwavearray`, `axisymmetric2D` set up its grid |
    | `opt_method`, `opt_seed` | Optimiser and seed of the `oneil_phase_only` preset |

    `elem_phase_correction` in the equipment YAML is used as the prior the fits are regularised towards. `save_elem_correction`, `opt_amp_validation`, and `opt_upper_velocity` are not used: offsets are written to the library entry, and since the drive scales exactly with √I<sub>SPPA</sub>, one k-Wave drive calibration per depth replaces a validation run per intensity.
