# Transducer Calibration

!!! note "At a glance"
    - **What:** How the phases and drive of a simulated transducer are fitted to hydrophone measurements of the real device, and how calibration relates to intensity targeting.
    - **When you need it:** To refine a [library](doc_transducer_library.md) entry to your own profiles, calibrate at depths and intensities of your own, or correct the steering of a matrix array.
    - **Before you start:** Hydrophone axial profiles of the device in free water ([Input file formats](#appendix-input-file-formats)).

PRESTUS calibrates the phase and amplitude settings of a virtual transducer so that its analytically predicted focal axis profile matches free-water measurements of the real driving system–transducer setup (hydrophone profiles measured in-house or provided by the manufacturer). Calibrated settings are stored in the [transducer library](doc_transducer_library.md) and loaded at run time. The library stores, for each measured depth, per-ring phase offsets, gains and the drive that reproduce the hydrophone profile (per-depth refinement, the default). Between measured depths they are interpolated; beyond them, geometric steering is used.

## Refining a library entry to your hydrophone profiles

Measure axial profiles in your water tank at the depths you will stimulate ([CSV format](#axial-profile-csv)). One call fits them and writes an entry to **your** folder, with a before/after report and, with `Plot`, a figure:

```matlab
refine_transducer_library('CTX_250_001_SC_105_010', 'my_profiles.csv', ...
    '/project/my_study/transducer_library', struct('Plot', true));
```

Then set `transducer.library_path` in the study config. `examples/demo_transducer_refinement.m` runs the whole path on demonstration profiles.

To inspect an entry, `plot_transducer_library(combo, library_path)` plots the ring phases, drive and ring gains across depths, and `generate_calibration_report(combo, library_path)` writes an HTML summary with the fit quality and the settings at the fitted depths.

**GUI.** The same steps are available interactively:

```matlab
calibration_gui(struct('OutputPath', '/project/my_study/transducer_library'))
```

**Predict** shows the fit at the depths you will stimulate before you save ([Predictions at desired depths](#predictions-at-desired-depths)). The Nexus's *Transducer ▸ Calibration* page gives the same results for the transducers of a study config and saves to the study's `transducer.library_path` ([Calibration page](doc_nexus.md#calibration-page)). A saved entry records its fit settings (mode, `FitOptions`, depths, profile CSV, reference ISPPA, interpolation and k-Wave options) in `meta.fit`; opening the entry restores them.

**Config-driven.** `calibration_transducer` (`examples/calibration_standalone.m`) runs the same calibration (`calibrate_library_entry`) and writes the library entry plus a calibration output folder; see [Calibration at chosen depths and intensities](#calibration-at-chosen-depths-and-intensities-calibration_transducer). To ingest per-run YAMLs (`{combo}-F*mm-I*wpercm2.yaml`) as refinements:

```matlab
update_transducer_library(combo_name, path_output_profiles, my_library_path)
```

??? note "Details: fit defaults and fit options"
    Per-depth fits free the phase offsets and per-ring gains (`FreeParams = {'phase','amp'}`), with an intensity-domain loss (`Domain = 'intensity'`, `HuberDelta = 0.8`). They are regularised towards geometric steering. The GUI and `build_transducer_library` use the same defaults. `z_shift` stays at 0, because a study does not apply it; correct a persistent axial misregistration in the equipment geometry or the profiles.

    **Fit options** (`FitOptions`, see `help fit_transducer_model`) in both GUIs:

    | Panel | Options |
    |---|---|
    | Main | Free parameters, forward model, near-field skip, loss domain, Huber δ, starts |
    | *Advanced…* | Phase constraint and priors (`Lambda`, `LambdaAmp`), gain bounds (`AmpBounds`, `PlausibleGain`), amplitude fit, focal weighting and fit window, near-field penalty (`NearFieldPenalty`, `NearFieldMargin`, `NearFieldMinMm`), medium (`CWater`, `Rho`), optimiser and seed, desired-profile `Tolerance`, k-Wave drive calibration `GridDims` and `ApplyCorrection` |
    | *Start from preset* | Loads a `calibration_fit_options` preset into both |

    Not offered: z shift and steering sound speed (not applied at run time), `ChannelGroups` (from the equipment), `Prior` (the entry) and `EvaluateOnly` (*Evaluate measured*).

### Predictions at desired depths

**Predict** models the profiles at the depths in **Desired depths**, for example `20:5:60` or `30, 42.5, 50` (default `min_foc:5:max_foc`). The depths need not have been measured. Each model is resolved as a study resolves it (`evaluate_library_entry` → `library_model_at_depth`):

- **Geometric**: geometric steering only.
- **Library**: the library entry with its refinements (the entry in **Save to folder** if there is one, else the shipped entry).
- **Fit**: the current unsaved fit, with k-Wave-calibrated drives if that calibration ran.

Interpolation and blending follow [How phases are computed](doc_transducer_library.md#how-phases-are-computed); drives and gains are held constant outside the refined span. `fit_z_shift_mm` is not applied.

The GUI plots modelled against measured (solid) or derived (dashed) profiles, and focus against target depth. The table lists, per desired depth: profile type, phase model (`refined`, `blend w=…`, `static` or `geometric`), drive at the reference ISPPA, predicted focus, R² and element phases. Without a profile CSV, the GUI opens in prediction view only (no R², no fitting).

**Calibrating at desired depths.** Set **Calibrate at** to *Desired depths* and press **Fit**. The GUI calibrates at every desired depth with a profile and reports skipped depths. With *Fit depths (measured)*, **Fit** uses the depths selected under **Fit depths**.

??? note "Details: profiles at desired depths, exports and headless refinement"
    **Profiles at desired depths.** With a profile CSV loaded, each desired depth gets a profile from the measurements (`interpolate_axial_profiles`):

    | Type | Source |
    |---|---|
    | Measured | The desired depth was measured |
    | Interpolated | The two neighbouring measured profiles, shifted so their foci coincide at the interpolated position, blended linearly by distance |
    | Extrapolated | Outside the measured range but within **extrapolate ≤** (default 0 mm, off): the nearest profile shifted to the desired depth |
    | None | Further outside the range |

    **Desired profiles** sets what the foci are aligned on: `label` (default, the column headers; see [focal-depth labels](#focal-depth-labels)) or `peak` (intensity peak beyond the near-field skip). The shifted near field is not physical; fits skip it.

    **Exports.**

    - **Save entry** writes `refinement.profile_type`. If any depth was fitted to a derived profile, it also writes `{combo}_calibration_profiles.csv` next to the entry and records it in `refinement.calibration_profiles_csv`.
    - **Export profiles…** writes the profiles of the last fit, or else those at the desired depths on view.
    - **Export settings…** writes per-depth settings of the fit (else the entry) at the depths on view, via `library_settings_table`: profile type, phase model, weight, drive, focus, and per ring phase, drive and gain.

    Profile CSVs use the format `read_axial_profiles` reads, with a YAML sidecar recording, for every derived column, its source depths, weights and shifts. `ExportFolder` writes both exports without a file dialog.

    **Headless.** `refine_transducer_library` accepts unmeasured depths in `Depths`: a depth within 0.5 mm of a measurement uses it; any other is interpolated (`InterpolateOptions`: `Align`, `MaxExtrapolationMm`). It writes the calibration profiles (`ExportProfiles`, on by default when any depth is derived) and optionally the settings table (`ExportSettings`, `SettingsDepths`):

    ```matlab
    refine_transducer_library('CTX_250_001_SC_105_010', 'my_profiles.csv', lib, struct( ...
        'Depths', [30 42.5 55], 'ExportSettings', true, ...
        'InterpolateOptions', struct('MaxExtrapolationMm', 3)));
    % writes {combo}.yaml, {combo}_calibration_profiles.csv (+ .yaml sidecar), {combo}_settings.csv
    T = library_settings_table(tran, lib_struct, 20:5:70, struct('IsppaWcm2', 30, 'CsvPath', 'settings.csv'));
    ```

### Optional: within-ring velocity profile

The library models each ring as a uniformly vibrating piston. For several Sonic Concepts CTX transducers no choice of phases and gains then reproduces the measured near-field to focal peak ratio. On CTX-500-026, the measured near-field peak (8–12 mm from the exit plane, 0.4–0.6 of the focal intensity at deep settings) comes out about 7× too weak. Uniform rings match Imasonic transducers well.

The ring profile splits every ring into `n_sub` sub-rings with a complex weight that is the **same at every focal depth**. Each depth keeps one phase and gain per channel. On-axis profiles alone do not determine the weights uniquely. **Without off-axis data, keep the library fit and treat the off-axis field of a ring profile as unvalidated.** Off-axis scans ([format](#plane-and-volume-scans-optional)) constrain the weights and are fitted together with the axial profiles.

??? note "Details: where to enable the ring profile"
    | Where | How to enable |
    |---|---|
    | `calibration_gui`, Nexus *Calibration* page | **Ring profile** → *fit (shared by all depths)*; **Options…** sets the `fit_ring_profile` options and a plane-scan file. The GUI shows the profile (`plot_ring_profile`) and R² per plane. |
    | `refine_transducer_library`, `calibrate_library_entry` | `RingProfile = 'fit'`; `RingProfileOptions`: `NSub`, `Smooth`, `NearFieldSkipMm`, `AmpBounds`, `Iterations`, `NStarts`, `Seed`, `HeldOut`, `OffAxis`, `OffAxisOptions`. With `Plot`, the profile figure is saved next to the entry. |
    | `config_calibration.yaml` (`calibration_transducer`) | `ring_profile: fit`, `ring_profile_n_sub`, `ring_profile_smooth`, `ring_profile_near_field_skip_mm`, `ring_profile_offaxis` (plane-scan file), `ring_profile_offaxis_kind` |
    | `build_transducer_library` | `RING_PROFILE = 'fit'` or a list of combos (default `'none'`) |

`'auto'` (the default everywhere) keeps an entry's profile if it has one; `'none'` removes it.

```matlab
ro = struct('OffAxis', 'SC CTX-500-026 (4 ch.).zip', 'OffAxisOptions', struct('Kind', 'transverse'));
[lib, fits, report] = calibrate_library_entry(lib, tran, meas, struct('RingProfile', 'fit', 'RingProfileOptions', ro));
```

??? note "Details: how the profile is fitted, what changes in the entry, validation"
    `fit_ring_profile` fits the weights jointly over all measured depths, with a smoothness penalty and the near field included from 4 mm. Each ring's area-weighted mean weight is 1.

    | | Uniform rings | Ring profile, depth held out |
    |---|---|---|
    | CTX-500-026 (4 rings), median R² (worst) | 0.957 (0.79) | 0.989 (0.93) |
    | CTX-250-009 (2 rings) | 0.894 (0.87) | 0.980 (0.95) |
    | CTX-250-001 (4 rings) | 0.978 (0.96) | 0.980 (0.96) |
    | CTX-500-026 near-field / focal ratio error | 7× | 33 % |

    Exploratory analysis: peak-normalised intensity from 4 mm, held-out depth with its channel drive refitted. "Uniform rings" here is an unregularised per-depth fit, not the shipped library fit.

    **What changes.**

    - The profile is stored in `global_model.ring_profile` (`n_sub`, `weight_abs`, `weight_phase_deg` per sub-ring, `smooth`, `near_field_skip_mm`, fit quality).
    - All per-element arrays of the entry then refer to sub-rings. Adding, removing or refitting a profile drops the entry's refinements (with a warning); only the depths calibrated then are stored.
    - The simulation gets an annular source of sub-rings, each driven with its channel's phase and drive times its weight (`ring_profile_drive`). `transducer.tx_calibration.ring_profile` records it. Sub-rings are steered with their ring (`annular.steer_r_mid_mm`).
    - Inline `elem_phase_deg` / `elem_p_amp` in a study config refer to the equipment rings and bypass the profile.
    - Dephased conditions are not supported with a profile.

    **Off-axis data.** A refit from synthetic axial profiles reproduces the profiles (R² > 0.998) but not the weights. Validation on CTX-500-026 hydrophone planes at the 40 mm setting (R² of intensity normalised by the axial peak, planes held out of the fit):

    | Model (fitted to 12 axial profiles) | Axis (40 mm) | Sagittal XZ (z < 25 mm) | XY 35 mm | XY 40 mm (focus) | −3 dB width |
    |---|---|---|---|---|---|
    | Shipped library fit (uniform rings, `rayleigh_phase_gain`) | 0.87 | 0.89 (0.26) | 0.96 | 0.97 | 2.7 mm |
    | Same, gain regularisation 10× weaker (`LambdaAmp` 1e-3) | 0.89 | 0.90 (0.33) | 0.94 | 0.96 | 2.6 mm |
    | Uniform rings, unregularised | 0.90 | 0.86 (0.31) | 0.90 | 0.90 | 2.4 mm |
    | Ring profile, axial | 0.97 | 0.88 (0.64) | 0.88 | 0.87 | 2.4 mm |
    | Ring profile, axial + XY 35 mm plane | 0.97 | 0.94 (0.65) | 0.95 (fitted) | 0.98 | 2.9 mm |
    | Hydrophone | | | | | 2.8 mm |

    The shipped fit predicts the transverse planes and focal width well but misses the axis and near field. A ring profile fitted to axial profiles alone fixes the axis and loses the planes. Adding one transverse plane gives both.

### Optional: free-water k-Wave drive calibration

`kwave_calibrate_drive` simulates every refined depth in free-water k-Wave as a study would resolve it, with ISPPA targeting off. It rescales the stored drive so the simulation hits the reference ISPPA (`elem_p_amp × √(I_ref / I_kWave)`, `amp_source: kwave`). Phases are validated, not changed.

**Why run it.** In free-water k-Wave the analytical drive overshoots the reference ISPPA, consistently per transducer and almost independently of depth: by about 3 % (IGT 32-channel, 10 rings), 6–9 % (CTX-250) and 14–16 % (CTX-500), over three depths each. A study inherits that error in every ISPPA-dependent output (intensity, MI, heating) when it skips the free-water baseline. That happens with no `transducer.target_isppa_wcm2`, or with a target equal to the entry's calibration ISPPA (within 1 %; the default 30 W/cm² is the shipped reference), because PRESTUS then trusts the entry's drive ([ISPPA targeting](doc_simulations-acoustic.md#isppa-targeting)). A different target, or `modules.run_water_baseline: 1`, removes the error at the cost of one water simulation per run. The drive calibration removes it once per entry.

**Run it.**

```matlab
refine_transducer_library('CTX_250_001_SC_105_010', 'my_profiles.csv', lib, struct( ...
    'KWaveDriveCalibration', true, 'KWaveDriveOptions', struct('SimConfig', 'my_study_config.yaml')));
```

??? note "Details: GUI, Nexus, simulation settings, grid and run time"
    - In `calibration_gui`, tick **k-Wave drive: calibrate in free water, after fit** (per-depth mode). The simulated profiles are overlaid in green, and **Save entry** stores the corrected drives.
    - The Nexus Calibration page writes the unsaved fit and a `kwave_drive_<combo>.m` script, to run on its Run tab or on the cluster.
    - `SimConfig` supplies grid, `code_type` and `platform`/HPC settings. Without it, PRESTUS runs locally with `config_default` (CPU on macOS).
    - `GridDims` (`kwave_grid_dims` in `config_calibration.yaml`) selects the free-water grid; `'auto'` sizes it to the aperture and focus (`grid.default_dims: 'auto'`).
    - Each depth is one free-water simulation: about 3–20 min per depth in 3D on 16 CPU cores, about 1 min on the **axisymmetric** grid. 3D removes the error; the axisymmetric grid leaves at most about 3 % (2D and 3D agree within 0.98–1.03 in ISPPA, with the same focus). Use 3D for final entries.

When the drive at the requested depth is not k-Wave-calibrated, loading an entry warns (`load_transducer_from_library:uncalibratedAmplitude`) and names the nearest checked depth with its R² and focal error.

??? note "Details: report fields and the uncalibrated-drive flag"
    **Report** (`kwave_drive_metrics`), which also scores the absolute level that peak-normalised R² cannot show:

    | Field | Meaning |
    |---|---|
    | `isppa_err_pct` | How far the analytical drive misses the reference ISPPA in k-Wave |
    | `i_at_target_ratio` | Simulated intensity at the measured FLHM centre relative to the reference |
    | `kwave_over_analytical` | k-Wave ISPPA over Rayleigh ISPPA at the same drive (the Rayleigh ISPPA of a shipped drive equals the reference to within 2 %) |
    | `nrmse_before`, `nrmse_after` | RMS difference from the hydrophone profile scaled to the reference ISPPA, before and after correction |
    | `sim_seconds`, `grid_points` | Cost |

    It also reports simulated ISPPA, R² and FLHM-centre error against the hydrophone profile.

    `load_transducer_from_library` sets `tx_calibration.amp_kwave_calibrated`, true only when the drive at the requested depth comes wholly from a refinement with `amp_source: kwave`. BabelBrain's `water_amp_calibration` ([doc_babelbrain.md](doc_babelbrain.md)) is independent: it corrects BabelBrain's FDTD against its own Rayleigh reference on every BabelBrain run.

## Calibration at chosen depths and intensities (`calibration_transducer`)

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

### Two endpoints

Calibration writes a library entry `{library_path}/{combo}.yaml`, which studies use through `transducer.library_path`, and a calibration output folder `path_output_profiles`.

??? note "Details: files written to each endpoint"
    1. **Library entry** `{library_path}/{combo}.yaml` (default `{path_output}/library`): the calibrated depths as refinements on top of the existing user entry, or else the shipped one. Studies use it by setting `transducer.library_path` to that folder.
    2. **Calibration output folder** `path_output_profiles` (default `{path_output}/profiles`):
        - `PRESTUS_virtual_{combo}.csv`: phases and amplitude per intensity (rows) × depth (columns)
        - `{combo}-F{depth}mm-I{intensity}wpercm2.yaml`: one PRESTUS transducer block per depth × intensity, with `target_isppa_wcm2` and `tx_calibration.isppa_wcm2` set to that intensity (can be ingested with `update_transducer_library`)
        - `{combo}_settings_I{intensity}.csv`: per-depth phases, per-ring amplitudes, gains, model focus (`library_settings_table`)
        - `{combo}_calibration_profiles.csv`: the measured or interpolated profiles that were fitted
        - `{combo}_calibration_report.csv`, `{combo}_calibration.png`, and `{combo}_kwave_drive_calibration.csv` if the k-Wave drive calibration ran

All values in the output folder are resolved from the library entry (`library_model_at_depth`), so both endpoints agree exactly. The drive for other intensities scales as √(I / I_ref).

### Configuration

Keys of `config_calibration.yaml` beyond the combinations, depths and intensities above.

??? note "Details: calibration config keys"
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
    | `ring_profile` | `auto` (default: keep the entry's) · `fit`: within-ring velocity profile shared by all depths (`fit_ring_profile`) · `none`: uniform rings; `ring_profile_n_sub` (4), `ring_profile_smooth` (0.3), `ring_profile_near_field_skip_mm` (4); `ring_profile_offaxis` (plane-scan file, `''` = axial only) and `ring_profile_offaxis_kind` (`all`). See [Within-ring velocity profile](#optional-within-ring-velocity-profile) |
    | `add_FDO` | 1 (default): profile distances are from the exit plane · 0: from the bowl apex |
    | `focal_reference` | FLHM convention of the profile headers; see [Focal depth convention](#focal-depth-convention-flhm) |
    | `kwave_drive_calibration` | Free-water k-Wave simulation per depth that rescales the stored drive to the simulated ISPPA (default true; `run_free_water_sim` is still read); `force_kwavearray`, `axisymmetric2D` set up its grid |
    | `opt_method`, `opt_seed` | Optimiser and seed of the `oneil_phase_only` preset |

    `elem_phase_correction` in the equipment YAML is used as the prior the fits are regularised towards. `save_elem_correction`, `opt_amp_validation` and `opt_upper_velocity` are not used: offsets are written to the library entry, and since the drive scales exactly with √ISPPA, one k-Wave drive calibration per depth replaces a validation run per intensity.

## Focal depth convention (FLHM)

A "focal depth" quoted for a real transducer (e.g. a manufacturer steering table, or `focal_depths_wrt_exit_plane` in a calibration config) is not the raw intensity-peak location. It is the centre of the **Focal Length at Half Maximum** (FLHM): the midpoint of the contiguous span around the peak where the axial profile stays at or above a threshold fraction of the peak. Manufacturers differ on what "half maximum" means:

| `focal_reference` | Threshold | Convention |
|---|---|---|
| `flhm_3db` (default) | intensity ≥ peak / 2 | ITRUSST / Martin et al. (2025); Sonic Concepts CTX + TPO ("NeuroFUS") steering tables |
| `flhm_6db` | intensity ≥ peak / 4 (≡ pressure ≥ peak / 2) | Imasonic / IGT steering tables |
| `peak` | n/a (raw intensity-peak location) | No FLHM centring |

Set `parameters.calibration.focal_reference` to the convention the characterisation data uses. A mismatch produces calibrated depths that look shifted by several mm, although the underlying phases and amplitudes are correct. Default: `flhm_3db`.

??? note "Details: where focal_reference is used"
    - Read by `compute_flhm_center` (via `focal_reference_to_frac`) wherever a calibration or post-hoc analysis locates the FLHM centre of an axial profile.
    - Echoed to the console during `calibration_setup` (`Focus from exit plane: 85.00 (FLHM -6dB)`).
    - Recorded in the sidecar YAML's `meta.calibration.focal_reference` by `save_optimized_values`, so calibrated outputs stay traceable to their convention.
    - Used by `compute_focal_metrics` (and therefore `acoustic_analysis`), so focal-centre reporting in the pipeline matches the calibration. This is independent of the fixed ITRUSST −3 dB / −6 dB **volume** fields (`mask_3dB`/`mask_6dB`), which always use both thresholds.

## Calibration methods

Per-depth refinement (analytical, no k-Wave) is the default; the k-Wave drive calibration optionally recalibrates its drive, and the O'Neil phase-only fit serves comparison with the `perform_global_search` pipeline.

??? note "Details: calibration levels compared"
    | Level | Entry point | Phase source | Amplitude | k-Wave needed | Typical use |
    |---|---|---|---|---|---|
    | **Per-depth refinement** (default) | `build_transducer_library`, `refine_transducer_library`, `calibration_gui`, `calibration_transducer` | Analytical fit per depth (`fit_transducer_model` via `calibrate_library_entry`) | Analytical, per depth | No | Library entries from hydrophone profiles; calibration at chosen depths and intensities |
    | **k-Wave drive calibration** (optional) | `kwave_calibrate_drive` (`KWave` option / GUI checkbox / `kwave_drive_calibration`) | Unchanged | Calibrated to free-water k-Wave ISPPA | Yes | Simulation-consistent drive at the refined depths |
    | **O'Neil phase-only fit** | `fit_preset: oneil_phase_only` / `Preset = 'oneil_phase_only'` | Per-depth fit with global-search settings (O'Neil, FEX-minimize, phases only) | Analytical, per depth | No | Comparison with results of the `perform_global_search` pipeline |

Most entries stay at the analytical amplitude level unless someone runs the k-Wave drive calibration. `load_transducer_from_library` resolves `tx_calibration.amp_kwave_calibrated` per requested depth and warns (`load_transducer_from_library:uncalibratedAmplitude`) whenever the driven amplitude has not been k-Wave-calibrated; see [free-water k-Wave drive calibration](#optional-free-water-k-wave-drive-calibration).

??? note "Details: how the model is fitted, fit settings and library model version"
    `fit_transducer_model` fits the on-axis field p(z) = H(z)·w, linear in the complex element drives w, to the measured profile at each depth.
    H comes from `axial_transfer_matrix` (Rayleigh by default).

    | Setting | Behaviour |
    |---|---|
    | Free parameters (`FreeParams`) | `phase` and `amp` (per-ring gain). `z_shift` and `c_steer` only when calling `fit_transducer_model` directly. Rings on one channel share phase and gain (the 30 sub-rings of IS_PCD15287 fit as 10 channels). |
    | Loss | Pseudo-Huber on the peak-normalised profile beyond `NearFieldSkipMm`, every depth weighted equally. |
    | Optimiser | `lsqnonlin` from `NStarts` starts (geometry first, then random). Needs the Optimization Toolbox; loading entries does not. |
    | Near field | The first `NearFieldSkipMm` (10 mm) are not fitted. `NearFieldPenalty` discourages a model near field hotter than the measured one; every fit reports `nearfield_excess`. |
    | Plausibility | Gains range over 0.1–10× (`AmpBounds`). A gain outside 0.5–2× the mean (`PlausibleGain`), or at a bound, is flagged in `model.flags`; check the element, cable or channel. |

    The static joint model (`LIBRARY_MODE = 'static'`, `Mode = 'static'`, `phase_mode: multi_depth`) raises an error. Entries that already hold a static offset still load.

    **Library model version.** Entries with `meta.library_model_version: 2` apply per-ring gains to `elem_p_amp`; version 1 applies a phase offset only. Version-3 entries (`library_mode: per_depth`) take phases, gains and drive from the refinement at each measured depth.

??? note "Details: one calibration core"
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

??? note "Details: library build and run-time loading"
    The shipped library is built by `build_transducer_library.m` from the profiles in `local/profiles/`: every measured depth becomes a refinement fitted by `fit_transducer_model` (phase offsets and gains on top of geometric steering from `set_real_phases`, plus the drive that reproduces the measured peak), written to `config/calibration/{combo_name}.yaml` (`LIBRARY_PATH`). The entry schema is documented in the [library YAML schema](doc_transducer_library.md#library-yaml-schema). At run time, `load_transducer_from_library(combo_name, focal_depth_ep, isppa, equip_param)` loads the calibrated parameters.

### Fit presets

Presets come from `calibration_fit_options`. Any `fit_transducer_model` option can be set on top of a preset (`fit_options` in the config, `FitOptions` in the functions).

??? note "Details: preset settings"
    | Preset | Free parameters | Forward model | Loss | Optimiser | Pull to prior |
    |---|---|---|---|---|---|
    | `rayleigh_phase_gain` (default) | phase offset + gain per ring | Rayleigh | pseudo-Huber on intensity (δ 0.8) | `lsqnonlin`, 16 starts | `Lambda` 1e-3 |
    | `oneil_phase_only` | phase offset per ring | O'Neil | (practically) least squares on intensity | FEX-minimize (or `GlobalSearch`) | none |

### Unmeasured depths

A requested depth within 0.5 mm of a measured one uses that measurement.

??? note "Details: interpolated profiles and the gap warning"
    Any other depth is interpolated from its two neighbouring measured profiles (aligned on their focus, `profile_align`) and then fitted, and is stored with `profile_type: interpolated`. A warning is issued when the neighbouring measurements are more than `max_gap_mm` (default 10 mm) apart.

??? note "Details: interpolating profiles vs interpolating fits"
    Between calibrated depths, studies interpolate the stored offsets and gains (`library_model_at_depth`). On held-out measured depths both ways predict the profile about equally well when the measured depths are up to ~10 mm apart (median R² 0.975 each, vs 0.976 for a direct fit and 0.90 for geometric steering). Fitting the interpolated profile is more robust where the focal behaviour changes quickly. Neither bridges a large gap across a change in focal behaviour.

## Calibration and ISPPA targeting

Calibration and intensity targeting are separate steps:

| | Calibration (`calibration_transducer`, library) | ISPPA targeting |
|---|---|---|
| **When** | Once per transducer/depth, before any study | Every run, after the acoustic simulation |
| **What is set** | Phases **and** drive per focal depth | Neither; the simulated field is scaled |
| **Input** | Hydrophone axial profiles | `transducer.target_isppa_wcm2` (free-water ISPPA, W/cm²; default 30) |
| **Output** | Library entry (YAML) used by studies | Pressure field scaled by √(target / free-water baseline) |

The acoustic simulation runs at the drive `elem_p_amp`. A short free-water baseline of the same source (`modules.run_water_baseline: auto`, the default) gives the free-water ISPPA of that drive, and the field is scaled to the target, so changing the target never re-runs the head simulation. A target that cannot be met stops the run; it is never reported unscaled. [ISPPA targeting](doc_simulations-acoustic.md#isppa-targeting) describes the rules, the provenance fields and the error messages; [Multi-intensity](doc_multiintensity.md) runs several targets from one acoustic simulation.

The typical workflow:

1. **Once per transducer (library build)**: run `build_transducer_library.m` to fit every measured depth. Results are stored in `config/calibration/{combo_name}.yaml`.
2. **Optionally, at chosen depths**: run `calibration_transducer` (or `refine_transducer_library`) with your own profiles or at unmeasured depths, writing a user library entry (`transducer.library_path`). Calibrated depths are used exactly there and blended into the rest of the entry within `blend_mm`.
3. **Per subject simulation**: set `transducer.serial`, `combo.ds_serial` and `focal_distance_ep` (PRESTUS loads the calibrated parameters via `load_transducer_from_library`) and `transducer.target_isppa_wcm2` to the desired free-water intensity. The water baseline and scaling handle the rest; the acoustic simulation is not re-run when the target intensity changes.

## Steering correction of matrix arrays

A matrix (2D) array's delivered focal peak can differ from the geometric model by a factor that changes with the steering target (element directivity, crosstalk, manufacturing spread). The steering correction measures that factor at a set of targets and scales the drive by it.

**How it is applied.** `create_matrix_karray` multiplies every element amplitude by k(target) and records it in `transducer.matrix.steering_k`. `elem_p_amp` stays the nominal drive setting, and the simulated field is what the device delivers at that setting. With [ISPPA targeting](doc_simulations-acoustic.md#isppa-targeting) the results are reported at the target ISPPA either way; the correction changes the drive that achieves them: `elem_p_amp × isppa_scale_p`, or in volts `drive_ref_v × elem_p_amp × isppa_scale_p / elem_p_amp_ref`. With `target_isppa_wcm2: []` the results are reported at the corrected drive. Switch it off with `transducer.matrix.use_steering_correction: false`.

**Building an entry.** Measure the focal peak at 30–40 targets spread through the steering range you will use, including its edges and the geometric focus, plus a few kept back to check the fit. Use one drive voltage and the geometric phases PRESTUS computes for each target. Write one CSV row per target, then build the entry with `build_steering_correction`.

??? note "Details: CSV columns and build call"
    | Column | Unit | Required | Meaning |
    |---|---|---|---|
    | `target_x_mm`, `target_y_mm`, `target_z_mm` | mm | yes | Steering target in the transducer frame: origin at the apex (curved) or array centre (flat), z along the beam |
    | `peak_pressure_mpa` | MPa | yes | Largest pressure amplitude in the plane through the target, normal to the beam axis |
    | `drive_v` | V | no | Drive voltage; rows at other voltages are scaled to the median voltage |
    | `peak_x_mm`, `peak_y_mm` | mm | no | Where the peak was found |
    | `model_pressure_mpa` | MPa | no | A modelled peak at the same drive, used instead of the Rayleigh model |

    Then, with the array's geometry in `tr` (from its equipment YAML) and `tr.matrix.elem_p_amp` the element amplitude that stands for the measurement drive:

    ```matlab
    % run from the PRESTUS folder; writes the library entry next to the shipped ones
    sc = build_steering_correction(parameters, tr, 'steering_peaks.csv', ...
        struct('out_file', fullfile('config', 'calibration', '<serial>.yaml'), 'serial', '<serial>'));
    ```

`resolve_transducer_from_serial` attaches the entry to a matrix transducer given by `serial` (`<serial>_<ds_serial>.yaml` first, then `<serial>.yaml`, or `transducer.library_path`). Without an entry the array is simulated uncorrected. The entry stores the fit, its cross-validated errors, `elem_p_amp_ref` and `drive_ref_v`. An inline `transducer.matrix.steering_correction` takes precedence.

**Limits.** Only the focal peak is corrected, not side lobes or field shape. The model assumes the beam along +z of the transducer frame and does not support Clover setups. The Rayleigh model treats each element as a baffled piston in free water; check it against k-Wave with a free-water run at a few targets.

??? note "Details: what is fitted"
    At each measured target x, `k(x) = measured / modelled focal peak`, both at the same drive. The modelled peak is the largest free-field pressure in the plane through the target, normal to the beam axis, from a Rayleigh integral over the array (`matrix_focal_peak`, `matrix_rayleigh_pressure`). The measured peak is the same quantity from a hydrophone scan of that plane. `fit_steering_correction` fits k as one constant (`scale`) or a full quadratic surface in x, y, z (`quadratic`, 10 coefficients, linear least squares), and scores both by leave-one-out cross-validation. The default `auto` keeps the quadratic surface only with at least 15 targets and a median cross-validated error at most 0.8 times that of the constant. A quadratic map is not extrapolated: targets outside the measured range use the value at its edge, with a warning.

    Dadgar-Kiani et al. (2026) fit such a correction with a Gaussian of the difference in MPa; PRESTUS fits the ratio, which applies to the drive at any level.

??? note "References"
    Dadgar-Kiani et al. (2026). *Calibrated simulations for dynamic focusing of ultrasound through the temporal window.* Cell Reports Methods 6:101568.

## Appendix: Input file formats

`refine_transducer_library` and `build_transducer_library` need one axial profile CSV per transducer–driving-system combo. Steering phases are computed from the ring geometry (`set_real_phases`), so manufacturer phase tables are not used. Plane and volume scans are optional: they constrain a [within-ring velocity profile](#optional-within-ring-velocity-profile).

### Axial profile CSV

```
,Focal depths (mm),,,,,,,
Distance WRT Exit plane (mm),40.0,50.0,60.0,70.0,80.0
0.5,0.12,0.31,0.04,0.01,0.0
1.0,0.18,0.44,0.08,0.02,0.0
...
```

| Row | Content |
|---|---|
| 1 (title) | Label strings only; skipped |
| 2 (depth) | First cell a label (ignored); then the focal depth of each column in mm from the exit plane |
| 3 onward (data) | Axial distance from the exit plane in mm; then intensity in W/cm² at each focal depth |

Rules:

- Separator: comma. Encoding: UTF-8 (a BOM is tolerated).
- Any number of focal-depth columns ≥ 1.
- Axial distance starts at or below the transducer face (0 mm is fine; negative values are ignored).
- Intensity is non-negative and always W/cm² (no units row). Trailing NaN columns are trimmed.

#### Focal-depth labels

Column headers are nominal driver setpoints. They correspond to the FLHM centre of the measured profile, at a manufacturer-dependent threshold ([focal depth convention](#focal-depth-convention-flhm)). On the Donders profiles (2026-09-30; median |header − FLHM centre|, first 10 mm skipped):

| Equipment | Convention | Median error |
|---|---|---|
| Sonic Concepts CTX + TPO | `flhm_3db` | 0.03–0.26 mm (−6 dB: 1.0–7.0 mm) |
| Imasonic + IGT 32 ch | `flhm_6db` | 0.04–0.08 mm (−3 dB: 2.1–2.6 mm) |
| Imasonic + IGT 128 ch | closest to `flhm_6db` | 0.60–0.75 mm (−3 dB: 1.2–1.7 mm) |

#### Depth coverage recommendations

Cover the full intended depth range and include the natural focal depth.

??? note "Details: depths per transducer type"
    | Transducer type | Minimum depths | Recommended |
    |---|---|---|
    | Imasonic (10-element annular) | 3 | 5–7 spaced 10 mm apart |
    | Sonic Concepts CTX (4-element) | 4 | 6–8, including natural focal depth |
    | Any | n/a | Cover full intended depth range; include the natural focal depth |

#### Profile folder and file naming

Measured profiles are site data in `local/profiles/` inside the PRESTUS root (ignored by git), one subfolder per transducer. `build_transducer_library` and `calibration_setup` read them from there.

??? note "Details: folder layout, file names and sidecar"
    ```
    local/profiles/
    ├── IS_PCD15473_01001/
    │   ├── tran-IS_PCD15473_01001_ds-IGT_32_ch_comb_10_ch_axial.csv
    │   ├── tran-IS_PCD15473_01001_ds-IGT_32_ch_comb_10_ch_axial.yaml         ← sidecar
    │   ├── tran-IS_PCD15473_01001_ds-IGT_32_ch_comb_10_ch_acq-OPM_axial.csv  ← repeat measurement
    │   ├── tran-IS_PCD15473_01001_ds-IGT_128_ch_comb_10_ch_axial.csv
    │   └── sourcedata/        ← never read: vendor spreadsheets, steering tables, deprecated files
    └── CTX_250_001/
        ├── tran-CTX_250_001_ds-SC_105_010_axial.csv
        └── tran-CTX_250_001_ds-SC_203_035_axial.csv
    ```

    File names follow `tran-{TRAN_SERIAL}_ds-{DS_SERIAL}[_acq-{LABEL}]_axial.csv`, using the serials of `config/equipment/` (underscores, not hyphens). `acq-` marks a repeat or alternative measurement of the same pair. One subfolder per transducer holds all its measurements whatever driving system recorded them.

    `char_data_path` in the equipment YAML is relative to this folder, for example `CTX_250_001/tran-CTX_250_001_ds-SC_105_010_axial.csv`. Point `calibration.path_input_axial` (for `calibration_setup`) or `CHAR_DATA_ROOT` (for `build_transducer_library`) at the folder.

    **Sidecar.** Each CSV has a YAML sidecar with the same stem, recording `quantity` and `units`, `distance_reference` (`exit_plane` | `bowl`), whether headers are `nominal` setpoints or `measured` positions (`focal_label_is`), their FLHM convention (`focal_label_reference`), acquisition conditions and provenance. PRESTUS does not read sidecars: set `calibration.add_FDO` and `calibration.focal_reference` in the calibration config to match it.

### Plane and volume scans (optional)

Off-axis hydrophone scans for the ring-profile fit, read by `read_offaxis_scans` (`OffAxis` in `RingProfileOptions`; *Ring profile ▸ Options…* in the GUIs; `ring_profile_offaxis` in `config_calibration.yaml`).

**What to measure.** One transverse plane a few millimetres before the focus, at a setting that also has an axial profile, made the ring profile identifiable on CTX-500-026 (see [Within-ring velocity profile](#optional-within-ring-velocity-profile)). Further planes or a volume add checks (e.g. of an off-centre or tilted axis) more than constraint.

Accepted files are a Sonic Concepts FUS Metrology export or a CSV with one row per point.

??? note "Details: accepted files, CSV columns and reading call"
    **Accepted files.**

    - **Sonic Concepts FUS Metrology export**: `FUS_Metrology_Export_BableBrain.mat`, its folder, or a `.zip` of the export folder with `Planes/` and `Axial/`. Each entry of the cell array `data` needs `type`, `focus` (focal setting, mm from the exit plane), `x`, `y`, `z` (mm) and `ISPPA.data` with `ISPPA.unit = 'W/cm2'`. Entries of type `Axial` are used only to rescale the scans.
    - **CSV**, one row per point:

    | Column | Meaning |
    |---|---|
    | `focus_ep_mm` | Focal setting the scan was recorded at [mm from the exit plane] |
    | `x_mm`, `y_mm` | Lateral position [mm]; the beam axis at x = y = 0 (or use `Centre = 'centroid'`) |
    | `z_ep_mm` | Distance from the exit plane [mm] |
    | `intensity_wcm2` | Intensity [W/cm²] |
    | `label` (optional) | Name of the scan; separates scans at the same setting |

    ```
    focus_ep_mm,x_mm,y_mm,z_ep_mm,intensity_wcm2,label
    40,-15.0,-15.0,35,0.012,XY 35
    40,-14.5,-15.0,35,0.014,XY 35
    ...
    ```

    ```matlab
    [oa, scans] = read_offaxis_scans('SC CTX-500-026 (4 ch.).zip', struct('Axial', meas));
    % oa: OffAxis input of fit_ring_profile; scans: maps for plotting
    ```

??? note "Details: how scans are read and reduced"
    - **Kind from the coordinates.** Points need not lie on a grid.

        | Varying coordinates | Kind | Reduced to |
        |---|---|---|
        | x, y (constant z) | Transverse plane | Azimuthally averaged radial profile |
        | x or y with z | Axial plane | (\|x\|, z) map, both sides folded together |
        | x, y and z | 3D volume | (r, z) map, averaged around one axis parallel to z in each z bin |

    - **One scan per setting and label.** Several planes at the same setting need distinct `label` values in a CSV; without them, their points read as one volume.
    - **Absolute level.** A scan is compared with the axial profile at its focal setting (within 1 mm) relative to that profile's peak, so it must be in W/cm² at the same drive. When a Sonic Concepts export includes its own axial profiles, each scan is rescaled by the ratio of the fitted axial peak to the export's at that setting.
    - **Axis.** The scan axis is x = y = 0, or with `Centre = 'centroid'` the intensity centroid of the top 50 % of a transverse plane or volume. A volume's axis is assumed parallel to z; tilt is not corrected.
    - **Reduction options** (`OffAxisOptions`): `Kind` (`all`, `transverse`, `axial`, `volume`, or a list), `Weight` (per scan, relative to one axial profile; default 1), `RMaxMm` (30), `BinMm` (radial bin, 0.5), `ZBinMm` (z bin of axial planes and volumes, 1; the fit's field matrix grows with the number of reduced points), `Centre`.

### Checklist before refining or rebuilding

- [ ] Axial profile CSV is readable by `readmatrix(path, 'NumHeaderLines', 1)`
- [ ] Row 2 contains focal depth values from column 2 onward (row 1 is a title row)
- [ ] The CSV is in `local/profiles/{TRAN_SERIAL}/`, named `tran-…_ds-…_axial.csv`, with a filled-in sidecar
- [ ] Column 1 contains axial distances in ascending order
- [ ] Intensity values are non-negative
- [ ] Equipment YAML validates: `check_equipment_config(equip_path)` passes
- [ ] Focal depths in CSV are within `[min_foc, max_foc]` for the transducer
