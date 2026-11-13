# Refining With Your Hydrophone Profiles

A shipped library entry is fitted to reference hydrophone profiles of the transducer model where such profiles exist, not to your own device. Refining fits the entry to the axial profiles you measure in your own water tank and writes it to your own folder, so a study reproduces the focus and drive of your setup. This page covers the refinement, with its optional within-ring velocity profile and free-water k-Wave drive calibration, and the profile files the fits read.

On this page:

- [Refining a Library Entry](#refining-a-library-entry): the usual path, from your measured profiles to your own library entry.
- [Input File Formats](#input-file-formats): the profile files the fits read.

!!! note "At a glance"
    - **When you need it:** You have hydrophone profiles of your own device and want a [library](doc_transducer_library.md) entry fitted to them.
    - **Before you start:** Free-water axial profiles in the [axial profile CSV](#axial-profile-csv) format; the [Optimization Toolbox](https://www.mathworks.com/products/optimization.html) for fitting. The focal depth convention and calibration methods are on [Calibration](doc_calibration.md).

## Refining a Library Entry

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

**Predict** shows the fit at the depths you will stimulate before you save ([Predictions at desired depths](#predictions-at-desired-depths)). The Nexus's *Transducer ▸ Calibration* page gives the same results for the transducers of a study config and saves to the study's `transducer.library_path` ([Calibration page](doc_nexus_calibration-localite.md#calibration-page)). A saved entry records its fit settings (mode, `FitOptions`, depths, profile CSV, reference ISPPA, interpolation and [k-Wave](https://github.com/ucl-bug/k-wave) options) in `meta.fit`; opening the entry restores them.

**Config-driven.** `calibration_transducer` (`examples/calibration_standalone.m`) runs the same calibration (`calibrate_library_entry`) and writes the library entry plus a calibration output folder; see [Calibration at chosen depths and intensities](doc_calibration.md#calibrating-at-chosen-depths-and-intensities). To ingest per-run YAMLs (`{combo}-F*mm-I*wpercm2.yaml`) as refinements:

```matlab
update_transducer_library(combo_name, path_output_profiles, my_library_path)
```

??? note "Fit defaults and fit options"
    Per-depth fits free the phase offsets and per-ring gains (`FreeParams = {'phase','amp'}`), with an intensity-domain loss (`Domain = 'intensity'`, `HuberDelta = 0.8`). They are regularised towards geometric steering. The GUI and `build_transducer_library` use the same defaults. `z_shift` stays at 0, because a study does not apply it; correct a persistent axial misregistration in the equipment geometry or the profiles.

    **Fit options** (`FitOptions`, see `help fit_transducer_model`) in both GUIs:

    | Panel | Options |
    |---|---|
    | Main | Free parameters, forward model, near-field skip, loss domain, Huber δ, starts |
    | *Advanced…* | Phase constraint and priors (`Lambda`, `LambdaAmp`), gain bounds (`AmpBounds`, `PlausibleGain`), amplitude fit, focal weighting and fit window, near-field penalty (`NearFieldPenalty`, `NearFieldMargin`, `NearFieldMinMm`), medium (`CWater`, `Rho`), optimiser and seed, desired-profile `Tolerance`, k-Wave drive calibration `GridDims` and `ApplyCorrection` |
    | *Start from preset* | Loads a `calibration_fit_options` preset into both |

    Not offered: z shift and steering sound speed (not applied at run time), `ChannelGroups` (from the equipment), `Prior` (the entry) and `EvaluateOnly` (*Evaluate measured*).

### Predictions at Desired Depths

**Predict** models the profiles at the depths in **Desired depths**, for example `20:5:60` or `30, 42.5, 50` (default `min_foc:5:max_foc`). The depths need not have been measured. Each model is resolved as a study resolves it (`evaluate_library_entry` → `library_model_at_depth`):

- **Geometric**: geometric steering only.
- **Library**: the library entry with its refinements (the entry in **Save to folder** if there is one, else the shipped entry).
- **Fit**: the current unsaved fit, with k-Wave-calibrated drives if that calibration ran.

Interpolation and blending follow [How phases are computed](doc_transducer_library.md#how-phases-are-computed); outside the refined span, the drive and gains of the nearest refinement are mixed with the global ones by the same blend weight, so beyond `blend_mm` (5 mm) the global drive and gains apply. Without a global drive, the nearest refined drive is used at any depth. `fit_z_shift_mm` is not applied.

The GUI plots modelled against measured (solid) or derived (dashed) profiles, and focus against target depth. The table lists, per desired depth: profile type, phase model (`refined`, `blend w=…`, `static` or `geometric`), drive at the reference ISPPA, predicted focus, R² and element phases. Without a profile CSV, the GUI opens in prediction view only (no R², no fitting).

**Calibrating at desired depths.** Set **Calibrate at** to *Desired depths* and press **Fit**. The GUI calibrates at every desired depth with a profile and reports skipped depths. With *Fit depths (measured)*, **Fit** uses the depths selected under **Fit depths**.

??? note "Profiles at desired depths, exports and headless refinement"
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

### Within-Ring Velocity Profile

Optional. The library models each ring as a uniformly vibrating piston. For several [Sonic Concepts](https://sonicconcepts.com/) CTX transducers no choice of phases and gains then reproduces the measured near-field to focal peak ratio. On CTX-500-026, the measured near-field peak (8–12 mm from the exit plane, 0.4–0.6 of the focal intensity at deep settings) comes out about 7× too weak. Uniform rings match [Imasonic](https://www.imasonic.com) transducers well.

The ring profile splits every ring into `n_sub` sub-rings with a complex weight that is the **same at every focal depth**. Each depth keeps one phase and gain per channel. On-axis profiles alone do not determine the weights uniquely. **Without off-axis data, keep the library fit and treat the off-axis field of a ring profile as unvalidated.** Off-axis scans ([format](#plane-and-volume-scans)) constrain the weights and are fitted together with the axial profiles.

??? note "Where to enable the ring profile"
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

??? note "How the profile is fitted, what changes in the entry, validation"
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

??? developer "Sub-ring source in the simulation"
    The simulation gets an annular source of sub-rings, each driven with its channel's phase and drive times its weight (`ring_profile_drive`). `transducer.tx_calibration.ring_profile` records it. Sub-rings are steered with their ring (`annular.steer_r_mid_mm`).

### Free-Water k-Wave Drive Calibration

Optional. `kwave_calibrate_drive` simulates every refined depth in free-water k-Wave as a study would resolve it, with ISPPA targeting off. It rescales the stored drive so the simulation hits the reference ISPPA (`elem_p_amp × √(I_ref / I_kWave)`, `amp_source: kwave`). Phases are validated, not changed.

**Why run it.** In free-water k-Wave the analytical drive overshoots the reference ISPPA, consistently per transducer and almost independently of depth: by about 3 % (IGT 32-channel, 10 rings), 6–9 % (CTX-250) and 14–16 % (CTX-500), over three depths each. A study inherits that error in every ISPPA-dependent output (intensity, MI, heating) when it skips the free-water baseline. That happens with a target equal to the entry's calibration ISPPA (within 1 %; the default 30 W/cm² is the shipped reference), because PRESTUS then trusts the entry's drive ([ISPPA targeting](doc_simulations-acoustic.md#isppa-targeting)). A different target, or `modules.run_water_baseline: 1`, removes the error at the cost of one water simulation per run. The drive calibration removes it once per entry.

**Run it.**

```matlab
refine_transducer_library('CTX_250_001_SC_105_010', 'my_profiles.csv', lib, struct( ...
    'KWaveDriveCalibration', true, 'KWaveDriveOptions', struct('SimConfig', 'my_study_config.yaml')));
```

??? note "GUI, Nexus, simulation settings, grid and run time"
    - In `calibration_gui`, tick **k-Wave drive: calibrate in free water, after fit** (per-depth mode). The simulated profiles are overlaid in green, and **Save entry** stores the corrected drives.
    - The Nexus Calibration page writes the unsaved fit and a `kwave_drive_<combo>.m` script, to run on its Run tab or on the cluster.
    - `SimConfig` supplies grid, `code_type` and `platform`/HPC settings. Without it, PRESTUS runs locally with `config_default` (CPU on macOS).
    - `GridDims` (`kwave_grid_dims` in `config_calibration.yaml`) selects the free-water grid; `'auto'` sizes it to the aperture and focus (`grid.default_dims: 'auto'`).
    - Each depth is one free-water simulation: about 3–20 min per depth in 3D on 16 CPU cores, about 1 min on the **axisymmetric** grid. 3D removes the error; the axisymmetric grid leaves at most about 3 % (2D and 3D agree within 0.98–1.03 in ISPPA, with the same focus). Use 3D for final entries.

When the drive at the requested depth is not k-Wave-calibrated, loading an entry warns (`load_transducer_from_library:uncalibratedAmplitude`) and names the nearest checked depth with its R² and focal error.

??? note "Report fields"
    **Report** (`kwave_drive_metrics`), which also scores the absolute level that peak-normalised R² cannot show:

    | Field | Meaning |
    |---|---|
    | `isppa_err_pct` | How far the analytical drive misses the reference ISPPA in k-Wave |
    | `i_at_target_ratio` | Simulated intensity at the measured FLHM centre relative to the reference |
    | `kwave_over_analytical` | k-Wave ISPPA over Rayleigh ISPPA at the same drive (the Rayleigh ISPPA of a shipped drive equals the reference to within 2 %) |
    | `nrmse_before`, `nrmse_after` | RMS difference from the hydrophone profile scaled to the reference ISPPA, before and after correction |
    | `sim_seconds`, `grid_points` | Cost |

    It also reports simulated ISPPA, R² and FLHM-centre error against the hydrophone profile.

    [BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain)'s `water_amp_calibration` ([doc_babelbrain.md](doc_babelbrain.md)) is independent: it corrects BabelBrain's FDTD against its own Rayleigh reference on every BabelBrain run.

??? developer "Uncalibrated-drive flag"
    `load_transducer_from_library` sets `tx_calibration.amp_kwave_calibrated`, true only when the drive at the requested depth comes wholly from a refinement with `amp_source: kwave`.

## Input File Formats

`refine_transducer_library` and `build_transducer_library` need one axial profile CSV per transducer–driving-system combo. Steering phases are computed from the ring geometry (`set_real_phases`), so manufacturer phase tables are not used. Plane and volume scans are optional: they constrain a [within-ring velocity profile](#within-ring-velocity-profile).

### Axial Profile CSV

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

#### Focal-Depth Labels

Column headers are nominal driver setpoints. They correspond to the FLHM centre of the measured profile, at a manufacturer-dependent threshold ([focal depth convention](doc_calibration.md#focal-depth-convention-flhm)). On the Donders profiles (2026-09-30; median |header − FLHM centre|, first 10 mm skipped):

| Equipment | Convention | Median error |
|---|---|---|
| Sonic Concepts CTX + TPO | `flhm_3db` | 0.03–0.26 mm (−6 dB: 1.0–7.0 mm) |
| Imasonic + IGT 32 ch | `flhm_6db` | 0.04–0.08 mm (−3 dB: 2.1–2.6 mm) |
| Imasonic + IGT 128 ch | closest to `flhm_6db` | 0.60–0.75 mm (−3 dB: 1.2–1.7 mm) |

#### Depth Coverage Recommendations

Cover the full intended depth range and include the natural focal depth.

??? note "Depths per transducer type"
    | Transducer type | Minimum depths | Recommended |
    |---|---|---|
    | Imasonic (10-element annular) | 3 | 5–7 spaced 10 mm apart |
    | Sonic Concepts CTX (4-element) | 4 | 6–8, including natural focal depth |
    | Any | n/a | Cover full intended depth range; include the natural focal depth |

#### Profile Folder and File Naming

Measured profiles are site data in `local/profiles/` inside the PRESTUS root (ignored by [Git](https://git-scm.com/)), one subfolder per transducer. `build_transducer_library` and `calibration_setup` read them from there.

??? note "Folder layout, file names and sidecar"
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

### Plane and Volume Scans

Optional. Off-axis hydrophone scans for the ring-profile fit, read by `read_offaxis_scans` (`OffAxis` in `RingProfileOptions`; *Ring profile ▸ Options…* in the GUIs; `ring_profile_offaxis` in `config_calibration.yaml`).

**What to measure.** One transverse plane a few millimetres before the focus, at a setting that also has an axial profile, made the ring profile identifiable on CTX-500-026 (see [Within-ring velocity profile](#within-ring-velocity-profile)). Further planes or a volume add checks (e.g. of an off-centre or tilted axis) more than constraint.

Accepted files are a Sonic Concepts FUS Metrology export or a CSV with one row per point.

??? note "Accepted files, CSV columns and reading call"
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

??? note "How scans are read and reduced"
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

### Checklist Before Refining or Rebuilding

- [ ] Axial profile CSV is readable by `readmatrix(path, 'NumHeaderLines', 1)`
- [ ] Row 2 contains focal depth values from column 2 onward (row 1 is a title row)
- [ ] The CSV is in `local/profiles/{TRAN_SERIAL}/`, named `tran-…_ds-…_axial.csv`, with a filled-in sidecar
- [ ] Column 1 contains axial distances in ascending order
- [ ] Intensity values are non-negative
- [ ] Equipment YAML validates: `check_equipment_config(equip_path)` passes
- [ ] Focal depths in CSV are within `[min_foc, max_foc]` for the transducer
