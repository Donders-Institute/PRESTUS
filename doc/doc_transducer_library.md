# Transducer Library

!!! note "At a glance"
    - **What:** The shipped library of transducers and their calibrations, and how to refine an entry to your own measurements.
    - **When you need it:** Whenever you simulate a device from the library: set `serial`, `combo.ds_serial`, `focal_distance_ep` and `target_isppa_wcm2` (Level 0 is all most users need).
    - **Before you start:** For refining (Level 2): hydrophone axial profiles of your own device.

PRESTUS ships a static transducer library and an optional path to refine it against your own hydrophone measurements.

| Folder | Content | Purpose |
|---|---|---|
| `config/equipment/` | One YAML per physical device | Ring geometry, frequency, steerable range |
| `config/calibration/` | One YAML per transducer–driving-system combo | Per-depth refinements (phases, gains, drive) on geometric steering |
| your folder (optional) | Same schema | Your own refinements; set `transducer.library_path` |

## Level 0: Use as-is (most users)

Name the transducer and driving system, the focal depth and the target intensity:

```yaml
transducer:
  serial: IS_PCD15473_01001
  combo:
    ds_serial: IGT_32_ch_comb_10_ch
  target_isppa_wcm2: 32.0
  focal_distance_ep: 60.0
  # resolves config/calibration/IS_PCD15473_01001_IGT_32_ch_comb_10_ch.yaml
```

PRESTUS resolves geometry and calibrated parameters at load time: the ring geometry from `config/equipment/`, and the phases, per-ring gains and drive at `focal_distance_ep` from the calibration entry. Fields already set in the study config override library values. [Shipped combos](#shipped-combos) lists the entries and how well each one fits its hydrophone profiles.

To use an entry you refined yourself ([Level 2](#level-2-refine-to-your-own-hydrophone-profiles-optional)), add `library_path`:

```yaml
transducer:
  serial: CTX_250_001
  combo:
    ds_serial: SC_105_010
  library_path: /project/my_study/transducer_library   # written by refine_transducer_library
  target_isppa_wcm2: 30.0
  focal_distance_ep: 35.0
```

The run records what was used: `transducer.tx_calibration.phase_model` (`geometric`, `static`, `refined` or `blend`), `amplitude_source` (`refinement`, `blend` or `global`) and `library_meta` (the library model version). Record the version with your results; results are not comparable across versions.

To bypass the library, give the full geometry, phases and drive inline: see [Transducer definition](doc_transducer.md#what-to-set).

??? note "Details: library keys and fallback chain"
    | Filename | When used |
    |---|---|
    | `{tran_serial}_{ds_serial}.yaml` | DS-specific calibration (all shipped entries) |
    | `{tran_serial}.yaml` | DS-agnostic entry (none shipped at present) |

    PRESTUS resolves in order:

    1. `transducer.library_path` (your folder), if set and it contains an entry for this transducer
    2. `config/calibration/{tran_serial}_{ds_serial}.yaml`
    3. `config/calibration/{tran_serial}.yaml`, with a warning
    4. If none exists, an error with instructions

    Without `combo.ds_serial`, PRESTUS looks only for `{tran_serial}.yaml`.

### How phases are computed

At a requested depth, the entry gives geometric steering plus a fitted per-ring correction, and the drive that reaches the reference ISPPA.

??? note "Details: phase and drive at a requested depth"
    For any requested focal depth *d* (mm from the exit plane):

    ```
    phase(d) = geometric steering phase(d; c_steer) + offset(d)
    drive_i  = elem_p_amp(d) · gain_i · exp(i · phase_i(d))
    ```

    - **Geometric steering** (`set_real_phases`) comes from the ring geometry in `config/equipment/`, at `global_model.steering_sound_speed` (default 1500 m/s).
    - **Refinement** (default): one fit per measured depth (`refinement.depths_ep_mm`). Each holds a per-ring phase offset (`elem_phase_correction_deg`), per-ring gains (`elem_amp_weight`) and the drive that gives the reference ISPPA at that depth (`elem_p_amp`). Between refined depths, offsets are interpolated per element on the unit circle; gains and drive are interpolated linearly. Outside that span, the nearest refinement fades back to the global model over `refinement.blend_mm` (default 5 mm).
    - **Global model**: geometric steering with the measured drive per depth (`global_model.amplitude_by_depth`). Some entries hold one depth-independent offset and gains instead (`global_model.elem_phase_correction_deg`, `elem_amp_weight`). These are read but cannot be fitted.

    **Amplitude.** Hydrophone profiles are recorded at a matched focal ISPPA per depth. The fit therefore uses profile *shape* only, with a free scale per depth, and derives the drive that produced the measured peak. For a refinement this is `refinement.elem_p_amp` at `refinement.isppa_wcm2`; for the global model it is `global_model.amplitude_by_depth`. Both are normalised to the reference ISPPA, interpolated in depth and scaled as √I at load time. Entries without either use the depth-independent `amplitude_scaling`. The optional [free-water k-Wave drive calibration](#optional-free-water-k-wave-drive-calibration) rescales the refinement drive to the simulated ISPPA.

### How the model is fitted

`fit_transducer_model` fits the on-axis field p(z) = H(z)·w, linear in the complex element drives w, to the measured profile at each depth.

??? note "Details: fit settings and library model version"
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

### Library YAML schema

??? note "Details: entry schema"
    ```yaml
    meta:
      tran_serial: CTX_250_001
      driving_system: SC_105_010
      calibration_type: empirical        # geometric | empirical | user_empirical
      created_at: '2026-06-12'
      prestus_version: <commit>

    global_model:                         # legacy static model (read only)
      elem_phase_correction_deg: [22.3, 56.3, 72.0, 349.9]   # one value per ring
      elem_amp_weight: [0.92, 1.03, 1.06, 0.99]              # per-ring gain, mean 1
      steering_sound_speed: 1500          # [m/s] steering law
      amplitude_scaling: {i30: 114442}    # depth-independent fallback
      amplitude_by_depth:                 # drive giving isppa_wcm2 in water per depth
        depths_ep_mm: [14, 16, ...]
        elem_p_amp: [98000, 101000, ...]
        isppa_wcm2: 30
      forward_model: rayleigh             # model the fit used
      free_params: [phase, amp, z_shift]
      fit_depths_ep_mm: [14, 16, ...]
      fit_z_shift_mm: 0.8                 # diagnostic only
      fit_r2: [0.95, 0.93, ...]           # in-sample R² per depth

    refinement:                           # per-depth (default for built and user entries)
      depths_ep_mm: [35, 40]
      elem_phase_correction_deg:
        - [25.1, 60.2, 70.3, 351.0]       # offset at 35 mm
        - [23.8, 58.9, 71.5, 350.2]       # offset at 40 mm
      elem_amp_weight:                    # per-ring gains, mean 1 (NaN: global gains)
        - [0.95, 1.02, 1.05, 0.98]
        - [0.97, 1.10, 1.03, 0.90]
      elem_p_amp: [98000, 104000]         # drive giving isppa_wcm2 in water (NaN: global drive)
      isppa_wcm2: 30
      amp_source: [analytical, kwave]     # kwave: calibrated by the free-water k-Wave drive calibration
      fit_z_shift_mm: [-1.5, -0.4]        # axial registration, diagnostic only
      profile_type: [measured, interpolated]   # profile each depth was fitted to (interpolate_axial_profiles)
      calibration_profiles_csv: /path/{combo}_calibration_profiles.csv   # written when any depth is derived
      kwave_isppa_wcm2: [.nan, 22.1]      # k-Wave ISPPA at the analytical drive (NaN: not checked)
      kwave_r2: [.nan, 0.93]              # k-Wave vs hydrophone profile shape
      kwave_focal_err_mm: [.nan, 0.4]
      blend_mm: 5
      forward_model: rayleigh
    ```

    `elem_phase_correction_deg` must have one entry per ring in the current equipment geometry. Otherwise `load_transducer_from_library` raises an error and `validate_transducer_library` reports `error`. Changing the ring geometry of a transducer invalidates its entries; rebuild them.

    Per-depth absolute phase tables (`global_model.depths_ep_mm` + `elem_phase_deg`) and parametric `calibration.focal_depths.f*` entries are converted to refinement offsets on load.

## Shipped combos

Agreement of the shipped entries (library model v3, `library_mode: per_depth`, built 2026-09-30) with the hydrophone axial profiles they were fitted to, from `meta.quality`. R² is on peak-normalised intensity beyond 10 mm from the exit plane, as a median over measured depths.

Columns: *geo*, geometric steering only; *in-sample*, the refinement at each measured depth (minimum over depths in brackets); *held-out*, each interior depth predicted from its two neighbours, as at run time between measured depths; *worst Δ*, the largest held-out loss against geometric steering at any depth.

| Combo | Manufacturer | Depths mm (n) | R² geo | R² in-sample (min) | R² held-out | worst Δ vs geo | Flags |
|---|---|---|---|---|---|---|---|
| IS_PCD15287_01001_IGT_128_ch_comb_10_ch | Imasonic | 7–96 (13) | 0.990 | 0.998 (0.930) | 0.998 | +0.000 | n/a |
| IS_PCD15287_01001_IGT_32_ch_comb_10_ch | Imasonic | 10–90 (10) | 0.986 | 0.998 (0.988) | 0.998 | -0.002 | n/a |
| IS_PCD15287_01002_IGT_128_ch_comb_10_ch | Imasonic | 7–94 (13) | 0.982 | 0.998 (0.889) | 0.998 | +0.002 | n/a |
| IS_PCD15287_01002_IGT_32_ch_comb_10_ch | Imasonic | 10–90 (10) | 0.988 | 0.999 (0.955) | 0.999 | +0.001 | n/a |
| IS_PCD15473_01001_IGT_128_ch_comb_10_ch | Imasonic | 11–93 (11) | 0.981 | 0.996 (0.980) | 0.995 | -0.001 | n/a |
| IS_PCD15473_01001_IGT_32_ch_comb_10_ch | Imasonic | 10–90 (10) | 0.983 | 0.998 (0.983) | 0.998 | -0.005 | n/a |
| IS_PCD15473_01001_OPM_IGT_32_ch_comb_10_ch | Imasonic | 10–90 (10) | 0.980 | 0.998 (0.975) | 0.997 | -0.043 | n/a |
| IS_PCD15473_01003_IGT_128_ch_comb_10_ch | Imasonic | n/a | n/a | n/a | n/a | n/a | geometric only (no profiles) |
| IS_PCD15473_01003_IGT_32_ch_comb_10_ch | Imasonic | 10–90 (10) | 0.982 | 0.997 (0.981) | 0.997 | -0.012 | n/a |
| IS_PCD15473_01003_OPM_IGT_32_ch_comb_10_ch | Imasonic | 10–90 (10) | 0.982 | 0.997 (0.984) | 0.997 | -0.003 | n/a |
| CTX_250_001_SC_105_010 | Sonic Concepts | 14–61 (12) | 0.786 | 0.984 (0.916) | 0.975 | -0.069 | n/a |
| CTX_250_001_SC_203_035 | Sonic Concepts | 14–61 (12) | 0.849 | 0.970 (0.904) | 0.938 | -0.277 | nearfield_excess |
| CTX_250_009_SC_105_010 | Sonic Concepts | 16–46 (12) | 0.838 | 0.917 (0.607) | 0.926 | -0.002 | nearfield_excess |
| CTX_250_009_SC_203_035 | Sonic Concepts | 16–38 (10) | 0.845 | 0.912 (0.243) | 0.913 | -0.004 | nearfield_excess |
| CTX_250_014_SC_105_010 | Sonic Concepts | 13–44 (12) | 0.869 | 0.906 (0.754) | 0.898 | -0.003 | n/a |
| CTX_250_014_SC_203_035 | Sonic Concepts | 13–44 (12) | 0.858 | 0.902 (0.751) | 0.892 | -0.001 | n/a |
| CTX_250_026_SC_105_010 | Sonic Concepts | 22–62 (12) | 0.852 | 0.975 (0.893) | 0.975 | +0.062 | nearfield_excess |
| CTX_250_026_SC_203_035 | Sonic Concepts | 22–61 (12) | 0.858 | 0.975 (0.911) | 0.976 | +0.056 | nearfield_excess |
| CTX_500_006_SC_105_010 | Sonic Concepts | 33–79 (10) | 0.910 | 0.921 (0.450) | 0.944 | -0.020 | n/a |
| CTX_500_006_SC_203_035 | Sonic Concepts | 36–78 (7) | 0.849 | 0.873 (0.477) | 0.924 | -0.002 | nearfield_excess |
| CTX_500_026_SC_105_010 | Sonic Concepts | 40–80 (12) | 0.909 | 0.973 (0.800) | 0.960 | +0.015 | n/a |
| CTX_500_026_SC_203_035 | Sonic Concepts | 40–79 (12) | 0.896 | 0.971 (0.832) | 0.959 | +0.013 | n/a |

`nearfield_excess`: the fitted model has a hotter near field than the measured profile at some depth.

Geometric steering alone reproduces the Imasonic profiles well (R² ≈ 0.98–0.99). For the Sonic Concepts CTX-500 / CTX-250 it is less accurate (R² ≈ 0.79–0.91), so refinement matters most there.

## Level 1: Geometry-only baseline (no hydrophone data)

For a transducer or driving system that is not in the shipped library:

```matlab
COMBOS_TO_RUN = {'MY_SERIAL_MY_DS_SERIAL'};              % a combo without char_data_path
LIBRARY_PATH  = '/project/my_study/transducer_library';  % omit to write config/calibration/
run(fullfile(prestus_root, 'functions', 'calibration', 'build_transducer_library.m'))
```

This writes an entry with an empty offset (pure geometric steering). Accuracy equals the geometric baseline in the table above. Adding the device to `config/equipment/` first is described in [Adding and maintaining library entries](doc_transducer_library_maintenance.md).

## Level 2: Refine to your own hydrophone profiles (optional)

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

**Config-driven.** `calibration_transducer` (`examples/calibration_standalone.m`) runs the same calibration (`calibrate_library_entry`) and writes the library entry plus a calibration output folder; see [Calibration at chosen depths and intensities](doc_calibration.md#calibration-at-chosen-depths-and-intensities-calibration_transducer). To ingest per-run YAMLs (`{combo}-F*mm-I*wpercm2.yaml`) as refinements:

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

### Optional: within-ring velocity profile

The library models each ring as a uniformly vibrating piston. For several Sonic Concepts CTX transducers no choice of phases and gains then reproduces the measured near-field to focal peak ratio. On CTX-500-026, the measured near-field peak (8–12 mm from the exit plane, 0.4–0.6 of the focal intensity at deep settings) comes out about 7× too weak. Uniform rings match Imasonic transducers well.

The ring profile splits every ring into `n_sub` sub-rings with a complex weight that is the **same at every focal depth**. Each depth keeps one phase and gain per channel. On-axis profiles alone do not determine the weights uniquely. **Without off-axis data, keep the library fit and treat the off-axis field of a ring profile as unvalidated.** Off-axis scans ([format](#plane-and-volume-scans-optional)) constrain the weights and are fitted together with the axial profiles.

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

### Adding transducers and rebuilding entries

Adding a new transducer or driving system (including `calibration_gui`'s **New…** form), rebuilding a shipped entry and contributing entries back are described in [Adding and maintaining library entries](doc_transducer_library_maintenance.md).

### Predictions at desired depths

**Predict** models the profiles at the depths in **Desired depths**, for example `20:5:60` or `30, 42.5, 50` (default `min_foc:5:max_foc`). The depths need not have been measured. Each model is resolved as a study resolves it (`evaluate_library_entry` → `library_model_at_depth`):

- **Geometric**: geometric steering only.
- **Library**: the library entry with its refinements (the entry in **Save to folder** if there is one, else the shipped entry).
- **Fit**: the current unsaved fit, with k-Wave-calibrated drives if that calibration ran.

Interpolation and blending follow [How phases are computed](#how-phases-are-computed); drives and gains are held constant outside the refined span. `fit_z_shift_mm` is not applied.

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

## Dephased conditions (experimental)

!!! warning "Experimental, not validated"
    No dephased hydrophone measurements have been compared with these predictions yet. Every use prints a `library_condition:experimental` warning. Do not report dephased-condition results without your own validation.

Dephased (defocused) settings drive the rings with explicit phases (e.g. alternate rings shifted by 180°) instead of steering to a focal depth. They cannot be interpolated in depth, so the library entry stores them as named **conditions**. Select one with `transducer.condition: geom_180` (with `serial`, `combo.ds_serial` and `target_isppa_wcm2`). Without `focal_distance_ep`, it defaults to the condition's `setpoint_ep_mm`.

**Adding conditions without measurements:**

```matlab
equip = load_equipment_config();
tran  = equip.trans.CTX_250_001;
lib   = yaml.loadFile('config/calibration/CTX_250_001_SC_203_035.yaml', 'ConvertToArray', true);
cond  = struct('label', 'geom_180', 'phase_shift_deg', [0 180 0 180], 'setpoint_ep_mm', 52.38);
lib   = library_add_conditions(lib, tran, cond, struct('IsppaWcm2', 30));
yaml.dumpFile('my_library/CTX_250_001_SC_203_035.yaml', lib);   % then transducer.library_path: my_library
```

**With measurements**, store them as `tran-{TRAN}_ds-{DS}_acq-dephased_axial.csv` in the usual CSV layout. The depth row holds each column's nominal setpoint. The sidecar describes the columns in order:

```yaml
columns:
  - label: geom_180
    phase_shift_deg: [0, 180, 0, 180]
  - label: geom_120
    elem_phase_deg: [0, 120, 240, 0]     # absolute driver phases
```

`library_ingest_conditions(lib_file, csv_path)` predicts every column, prints the agreement table and writes the conditions to `lib_file` (`struct('DryRun', true)` to only score). Do not pass dephased CSVs to `build_transducer_library` or `refine_transducer_library`, which read every column as a focal depth.

**Before treating a condition as validated:** compare `predicted_r2` and the focal/peak errors with those of the focused depths of the same entry, check `peak_z_ep_mm` against the measured peak, and confirm the free-water ISPPA in a k-Wave water run. Then set `validated: true` in your copy of the entry.

??? note "Details: condition schema and how phases and drive are derived"
    ```yaml
    conditions:
      geom_180:
        status: experimental
        validated: false
        phase_mode: shift                    # shift | absolute
        phase_shift_deg: [0, 180, 0, 180]    # shift mode: added to the setpoint's steering phases
        nominal_phase_deg: [...]             # driver phases before the hardware correction
        setpoint_ep_mm: 52.38                # nominal focal setpoint the condition derives from
        inherit_correction: true
        elem_phase_correction_deg: [...]     # hardware offset of the entry at the setpoint
        elem_phase_deg: [...]                # phases used in simulation
        elem_amp_weight: [...]               # per-ring gains of the setpoint
        elem_p_amp: 98000                    # drive for isppa_wcm2 at the model's free-water peak
        isppa_wcm2: 30
        peak_z_ep_mm: 18.5                   # where that peak is (can be a near-field or split focus)
        predicted_r2: .nan                   # model vs measured profile (NaN = no measurement)
        predicted_focal_err_mm: .nan
        source_csv: ''
    ```

    - `shift` mode: geometric steering phases at the setpoint plus `phase_shift_deg`.
    - `absolute` mode: the given driver phases.
    - `inherit_correction` adds the entry's hardware offset and gains at the setpoint. This assumes the correction fitted to focused profiles transfers to dephased drives; a dephased measurement has to confirm it.
    - The drive is set so that the Rayleigh free-water on-axis peak beyond 10 mm reaches `isppa_wcm2`.
    - Measured profiles only score the prediction (`predicted_r2`, focal and peak error); they are never fitted.

## Steering correction of matrix arrays

A matrix (2D) array's delivered focal peak can differ from the geometric model by a factor that changes with the steering target (element directivity, crosstalk, manufacturing spread). The steering correction measures that factor at a set of targets and scales the drive by it.

**How it is applied.** `create_matrix_karray` multiplies every element amplitude by k(target) and records it in `transducer.matrix.steering_k`. `elem_p_amp` stays the nominal drive setting, and the simulated field is what the device delivers at that setting. With [ISPPA targeting](doc_simulations-acoustic.md#isppa-targeting) the results are reported at the target ISPPA either way; the correction changes the drive that achieves them: `elem_p_amp × isppa_scale_p`, or in volts `drive_ref_v × elem_p_amp × isppa_scale_p / elem_p_amp_ref`. With `target_isppa_wcm2: []` the results are reported at the corrected drive. Switch it off with `transducer.matrix.use_steering_correction: false`.

**Building an entry.** Measure the focal peak at 30–40 targets spread through the steering range you will use, including its edges and the geometric focus, plus a few kept back to check the fit. Use one drive voltage and the geometric phases PRESTUS computes for each target. Write one CSV row per target:

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

## Library limitations

| Limitation | Effect | Mitigation |
|---|---|---|
| Refinements are exact only at measured depths | Between depths, interpolated phases fit less well (held-out interpolation in `meta.quality`) | Measure closely spaced depths around the targets |
| Axial shift fitted, not applied | A real model/hydrophone offset remains in focal placement | Correct the equipment geometry if `fit_z_shift_mm` is persistent |
| Refinement is only as good as your profiles | Offsets inherit hydrophone/tank errors | Cover the intended depth range; record conditions |
| Water-tank basis | Does not account for tissue acoustic properties | Combined with the PRESTUS forward simulation in tissue |
| Profiles recorded at matched ISPPA per depth | Depth dependence of focal gain is not calibrated | Water-baseline ISPPA rescaling |
| DS-agnostic fallback | Ignores DS hardware phase offset | Use DS-specific entries |
| Outside `[min_foc, max_foc]` | Unknown accuracy, warning emitted | Keep `focal_distance_ep` within range |

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

Column headers are nominal driver setpoints. They correspond to the FLHM centre of the measured profile, at a manufacturer-dependent threshold ([focal depth convention](doc_calibration.md#focal-depth-convention-flhm)). On the Donders profiles (2026-09-30; median |header − FLHM centre|, first 10 mm skipped):

| Equipment | Convention | Median error |
|---|---|---|
| Sonic Concepts CTX + TPO | `flhm_3db` | 0.03–0.26 mm (−6 dB: 1.0–7.0 mm) |
| Imasonic + IGT 32 ch | `flhm_6db` | 0.04–0.08 mm (−3 dB: 2.1–2.6 mm) |
| Imasonic + IGT 128 ch | closest to `flhm_6db` | 0.60–0.75 mm (−3 dB: 1.2–1.7 mm) |

#### Depth coverage recommendations

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
