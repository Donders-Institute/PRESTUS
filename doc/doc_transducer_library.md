# Transducer Library

!!! note "At a glance"
    - **What:** The shipped library of transducers and their calibrations, and how to refine an entry to your own measurements.
    - **When you need it:** Whenever you simulate a device from the library: set `serial`, `combo.ds_serial`, `focal_distance_ep` and `target_isppa_wcm2` (Level 0 below is all most users need).
    - **Before you start:** For refining: hydrophone axial profiles of your own device.

PRESTUS ships a static transducer library and an optional path to refine it against your own hydrophone measurements.

| Folder | Content | Purpose |
|---|---|---|
| `config/equipment/` | One YAML per physical device | Ring geometry, frequency, steerable range |
| `config/calibration/` | One YAML per transducer–driving-system combo | Per-depth refinements (phases, gains, drive) on geometric steering |
| your folder (optional) | Same schema | Your own refinements; set `transducer.library_path` |

When a study config contains `transducer.serial`, PRESTUS resolves geometry and calibrated parameters automatically at load time. Fields already set in the study config always override library values.

---

## How phases are computed

For any requested focal depth *d* (mm from the exit plane):

```
phase(d) = geometric steering phase(d; c_steer) + offset(d)
drive_i  = elem_p_amp(d) · gain_i · exp(i · phase_i(d))
```

- **Geometric steering** (`set_real_phases`) is computed from the ring geometry in `config/equipment/`, at the steering sound speed `global_model.steering_sound_speed` (default 1500 m/s).
- **Refinement** (default): one fit per measured depth (`refinement.depths_ep_mm`), holding a per-ring phase offset (`elem_phase_correction_deg`), per-ring gains (`elem_amp_weight`) and the drive that gives the reference ISPPA at that depth (`elem_p_amp`). Between refined depths, offsets are interpolated per element on the unit circle, and gains and drive linearly. Outside that span, the nearest refinement fades back to the global model over `refinement.blend_mm` (default 5 mm), so phases and drive stay continuous in depth.
- **Global model**: geometric steering with the measured drive per depth (`global_model.amplitude_by_depth`). Older entries may instead hold one depth-independent static offset and gains (`global_model.elem_phase_correction_deg`, `elem_amp_weight`); these are still read, but can no longer be fitted (see below).

`transducer.tx_calibration.phase_model` records what was used: `geometric`, `static`, `refined` or `blend`; `amplitude_source` records `refinement`, `blend` or `global`.

**Amplitude.** Hydrophone profiles are recorded at a matched focal ISPPA per depth. The fit therefore uses profile *shape* only, with a free scale per depth, and derives the drive that produced the measured peak at each depth. For a refinement this is the drive with its own refined phases and gains (`refinement.elem_p_amp` at `refinement.isppa_wcm2`); for the global model it is `global_model.amplitude_by_depth`. Both are normalised to the reference ISPPA, interpolated in depth and scaled as √I at load time. Entries without either use the depth-independent `amplitude_scaling`. An optional free-water k-Wave drive calibration (Level 2) rescales the refinement drive to the simulated ISPPA.

### How the model is fitted

`fit_transducer_model` uses the fact that the on-axis field is linear in the complex element drives: p(z) = H(z)·w. The transfer matrix H is computed once per transducer (`axial_transfer_matrix`, Rayleigh by default), so each evaluation is a matrix product. Details:

- **Free parameters** (`FreeParams`): `phase`, `amp` (per-ring gain), and for `fit_transducer_model` called directly `z_shift` (axial registration of the hydrophone profiles; a fit nuisance, not applied at run time) and `c_steer`; the calibration entry points fit `phase` and `amp` only. Rings driven by one channel share phase and gain; for example, the 30 sub-rings of IS_PCD15287 fit as 10 channels.
- **Loss**: pseudo-Huber on peak-normalised *pressure* amplitude beyond `NearFieldSkipMm`, with every depth weighted equally. Secondary lobes and split foci therefore count, and isolated artefacts don't dominate.
- **Optimiser**: `lsqnonlin` from `NStarts` starts (geometry first, then random), because profiles with two foci have several local optima. `lsqnonlin` needs the Optimization Toolbox; loading library entries does not.
- **Near field**: the first `NearFieldSkipMm` (10 mm) are not fitted. A soft penalty (`NearFieldPenalty`) discourages a model near field hotter than the measured one, and every fit reports `nearfield_excess`, because a hot spot near the exit plane matters in skull simulations.
- **Plausibility**: gains may range over 0.1–10× (`AmpBounds`). A gain outside 0.5–2× the mean (`PlausibleGain`), or at a bound, is flagged (`model.flags`) instead of being shipped silently: it more likely indicates a faulty element, cable or driver channel than something to calibrate away.

**Static joint model (removed 2026-10).** One depth-independent per-ring offset and gain fitted jointly over all depths (`LIBRARY_MODE = 'static'`, `Mode = 'static'`, `phase_mode: multi_depth`, GUI *Static model (joint, debug)*) was kept for debugging only. No benchmark showed it ahead of per-depth refinement, and all shipped entries are per-depth, so it was removed; those options now raise an error. Entries that already hold a static offset are still loaded.

**Library model version.** Entries with `meta.library_model_version: 2` apply per-ring gains to `elem_p_amp`, so simulated pressure and heating differ from runs with version-1 entries (phase offset only). Version-3 entries (`library_mode: per_depth`) take phases, gains and drive from the refinement at each measured depth. Results are not comparable across versions; record the version with your results (`transducer.tx_calibration.library_meta`).

### Library YAML schema

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

`elem_phase_correction_deg` must have one entry per ring in the current equipment geometry. If it doesn't, the entry was built for a different geometry: `load_transducer_from_library` raises an error and `validate_transducer_library` reports `error`.

Older entries are still read. Per-depth absolute phase tables (`global_model.depths_ep_mm` + `elem_phase_deg`) and parametric `calibration.focal_depths.f*` entries are converted to refinement offsets, relative to geometric steering at each depth.

---

## Library keys and fallback chain

| Filename | When used |
|---|---|
| `{tran_serial}_{ds_serial}.yaml` | DS-specific calibration (all shipped entries) |
| `{tran_serial}.yaml` | DS-agnostic entry (none shipped at present) |

PRESTUS resolves in order:
1. `transducer.library_path` (your folder), if set and it contains an entry for this transducer
2. `config/calibration/{tran_serial}_{ds_serial}.yaml`
3. `config/calibration/{tran_serial}.yaml`, with a warning
4. If none exists, an error with instructions

Specifying `combo.ds_serial` in the study config enables the DS-specific key. Omitting it looks only for `{tran_serial}.yaml`.

---

## Study config: minimal to fully specified

**Minimal: let the library fill in phases and amplitude:**

```yaml
transducer:
  serial: IS_PCD15473_01001
  combo:
    ds_serial: IGT_32_ch_comb_10_ch
  target_isppa_wcm2: 32.0
  focal_distance_ep: 60.0
  # resolves config/calibration/IS_PCD15473_01001_IGT_32_ch_comb_10_ch.yaml
```

**With your own refinement:**

```yaml
transducer:
  serial: CTX_250_001
  combo:
    ds_serial: SC_105_010
  library_path: /project/my_study/transducer_library   # written by refine_transducer_library
  target_isppa_wcm2: 30.0
  focal_distance_ep: 35.0
```

**Fully specified: bypass the library entirely:**

```yaml
transducer:
  serial: IS_PCD15473_01001
  type: annular
  freq_hz: 300000.0
  focal_distance_ep: 60.0
  annular:
    elem_n: 10
    elem_id_mm: [10.0, 22.1, ...]
    elem_od_mm: [21.1, 28.8, ...]
    curv_radius_mm: 100.0
    dist_geom_ep_mm: 92.7
    elem_phase_deg: [0, 45, 90, ...]
    elem_p_amp: [150000, 150000, ...]
```

---

## Shipped combos

Agreement of the shipped entries (library model v3, `library_mode: per_depth`, built 2026-09-30) with the hydrophone axial profiles they were fitted to, from each entry's `meta.quality`. R² is on peak-normalised intensity beyond 10 mm from the exit plane, as a median over measured depths. *geo* is geometric steering only; *in-sample* uses the refinement at each measured depth (minimum over depths in brackets); *held-out* predicts each interior depth from its two neighbours, as at run time between measured depths; *worst Δ* is the largest held-out loss against geometric steering at any depth.

| Combo | Manufacturer | Depths mm (n) | R² geo | R² in-sample (min) | R² held-out | worst Δ vs geo | Flags |
|---|---|---|---|---|---|---|---|
| IS_PCD15287_01001_IGT_128_ch_comb_10_ch | Imasonic | 7–96 (13) | 0.990 | 0.998 (0.930) | 0.998 | +0.000 | — |
| IS_PCD15287_01001_IGT_32_ch_comb_10_ch | Imasonic | 10–90 (10) | 0.986 | 0.998 (0.988) | 0.998 | -0.002 | — |
| IS_PCD15287_01002_IGT_128_ch_comb_10_ch | Imasonic | 7–94 (13) | 0.982 | 0.998 (0.889) | 0.998 | +0.002 | — |
| IS_PCD15287_01002_IGT_32_ch_comb_10_ch | Imasonic | 10–90 (10) | 0.988 | 0.999 (0.955) | 0.999 | +0.001 | — |
| IS_PCD15473_01001_IGT_128_ch_comb_10_ch | Imasonic | 11–93 (11) | 0.981 | 0.996 (0.980) | 0.995 | -0.001 | — |
| IS_PCD15473_01001_IGT_32_ch_comb_10_ch | Imasonic | 10–90 (10) | 0.983 | 0.998 (0.983) | 0.998 | -0.005 | — |
| IS_PCD15473_01001_OPM_IGT_32_ch_comb_10_ch | Imasonic | 10–90 (10) | 0.980 | 0.998 (0.975) | 0.997 | -0.043 | — |
| IS_PCD15473_01003_IGT_128_ch_comb_10_ch | Imasonic | — | — | — | — | — | geometric only (no profiles) |
| IS_PCD15473_01003_IGT_32_ch_comb_10_ch | Imasonic | 10–90 (10) | 0.982 | 0.997 (0.981) | 0.997 | -0.012 | — |
| IS_PCD15473_01003_OPM_IGT_32_ch_comb_10_ch | Imasonic | 10–90 (10) | 0.982 | 0.997 (0.984) | 0.997 | -0.003 | — |
| CTX_250_001_SC_105_010 | Sonic Concepts | 14–61 (12) | 0.786 | 0.984 (0.916) | 0.975 | -0.069 | — |
| CTX_250_001_SC_203_035 | Sonic Concepts | 14–61 (12) | 0.849 | 0.970 (0.904) | 0.938 | -0.277 | nearfield_excess |
| CTX_250_009_SC_105_010 | Sonic Concepts | 16–46 (12) | 0.838 | 0.917 (0.607) | 0.926 | -0.002 | nearfield_excess |
| CTX_250_009_SC_203_035 | Sonic Concepts | 16–38 (10) | 0.845 | 0.912 (0.243) | 0.913 | -0.004 | nearfield_excess |
| CTX_250_014_SC_105_010 | Sonic Concepts | 13–44 (12) | 0.869 | 0.906 (0.754) | 0.898 | -0.003 | — |
| CTX_250_014_SC_203_035 | Sonic Concepts | 13–44 (12) | 0.858 | 0.902 (0.751) | 0.892 | -0.001 | — |
| CTX_250_026_SC_105_010 | Sonic Concepts | 22–62 (12) | 0.852 | 0.975 (0.893) | 0.975 | +0.062 | nearfield_excess |
| CTX_250_026_SC_203_035 | Sonic Concepts | 22–61 (12) | 0.858 | 0.975 (0.911) | 0.976 | +0.056 | nearfield_excess |
| CTX_500_006_SC_105_010 | Sonic Concepts | 33–79 (10) | 0.910 | 0.921 (0.450) | 0.944 | -0.020 | — |
| CTX_500_006_SC_203_035 | Sonic Concepts | 36–78 (7) | 0.849 | 0.873 (0.477) | 0.924 | -0.002 | nearfield_excess |
| CTX_500_026_SC_105_010 | Sonic Concepts | 40–80 (12) | 0.909 | 0.973 (0.800) | 0.960 | +0.015 | — |
| CTX_500_026_SC_203_035 | Sonic Concepts | 40–79 (12) | 0.896 | 0.971 (0.832) | 0.959 | +0.013 | — |

`nearfield_excess`: the fitted model has a hotter near field than the measured profile at some depth (see *Near field* above).

**Accuracy notes**

- **Imasonic**: geometric steering alone already reproduces the measured profiles well (R² ≈ 0.98–0.99).
- **Sonic Concepts CTX-500 / CTX-250**: geometric steering is less accurate (R² ≈ 0.79–0.91). These are the transducers where an empirical offset, and optionally depth-specific refinement, matters most.

---

## Calibration levels

### Level 0: Use as-is (most users)

No action needed. Set `serial`, `combo.ds_serial`, `focal_distance_ep` and `target_isppa_wcm2`.

### Level 1: Geometry-only baseline (no hydrophone data)

For a transducer or driving system that isn't in the shipped library:

```matlab
COMBOS_TO_RUN = {'MY_SERIAL_MY_DS_SERIAL'};              % a combo without char_data_path
LIBRARY_PATH  = '/project/my_study/transducer_library';  % omit to write config/calibration/
run(fullfile(prestus_root, 'functions', 'calibration', 'build_transducer_library.m'))
```

This writes an entry with an empty offset (pure geometric steering). Accuracy equals the geometric baseline in the table above.

### Level 2: Refine to your own hydrophone profiles (optional)

Measure axial profiles in your own water tank at the depths you will stimulate (CSV format below). One call fits them and writes an entry to **your** folder, leaving the shipped library unchanged. It also prints a before/after report and, with `Plot`, saves a figure next to the entry:

```matlab
refine_transducer_library('CTX_250_001_SC_105_010', 'my_profiles.csv', ...
    '/project/my_study/transducer_library', struct('Plot', true));
```

Then set `transducer.library_path` in the study config. At the measured depths PRESTUS uses your refinement. Between them it interpolates, and more than `blend_mm` outside them it returns to the global model. `examples/demo_transducer_refinement.m` runs the whole path on demonstration profiles in a few seconds.

To inspect an entry, `plot_transducer_library(combo, library_path)` plots the ring phases, drive and ring gains it resolves to across depths, and `generate_calibration_report(combo, library_path)` writes an HTML summary with the fit quality and the settings at the fitted depths.

Per-depth fits free the phase offsets and per-ring gains (`FreeParams = {'phase','amp'}`), with an intensity-domain loss (`Domain = 'intensity'`, `HuberDelta = 0.8`). They are regularised towards the entry's global model, which is geometric steering for per-depth entries. The same defaults apply in the GUI and in `build_transducer_library`. The axial registration `z_shift` stays at 0: a study does not apply it, so phases fitted with it would match the hydrophone only after a shift the entry never gets. The GUI therefore does not offer it. A persistent axial misregistration points to an error in the equipment geometry or exit-plane distance, which should be corrected in the profiles or the equipment config rather than fitted.

You can also use the GUI:

```matlab
calibration_gui(struct('OutputPath', '/project/my_study/transducer_library'))
```

Before saving, **Predict** shows what the fit gives at the depths you will stimulate, including depths between and beyond the measured ones (see [Predictions at desired depths](#predictions-at-desired-depths)).

The same steps are on the PRESTUS Nexus's *Transducer ▸ Calibration* page, for the transducers of a study config: it saves to the study's `transducer.library_path` and can predict at the study's own focal depths ([Calibration page](doc_nexus.md#calibration-page)). Both use `calibration_desired_profiles`, `calibration_evaluate`, `calibration_fit` and `calibration_save_entry`, so they give the same results. A saved entry records how it was fitted in `meta.fit`: the fit mode, free parameters and `FitOptions`, whether it was calibrated at measured or desired depths (and those depths), the profile CSV, the reference ISPPA, the profile interpolation options and, if run, the k-Wave drive check settings. When `calibration_gui` or the Nexus opens an entry with `meta.fit`, they restore these settings, so a refit starts from the last one.

Fit settings are passed through `FitOptions` (see `help fit_transducer_model`).

#### Optional: within-ring velocity profile

The library models each ring as a uniformly vibrating piston. For several Sonic Concepts CTX transducers that model cannot reproduce the measured axial profiles, whatever phases and gains are chosen: in particular the ratio of the near-field peak to the focal peak. On CTX-500-026, the measured near-field peak (8–12 mm from the exit plane, 0.4–0.6 of the focal intensity at deep settings, reproduced by in-house and manufacturer measurements) comes out about 7× too weak. Imasonic transducers are matched well by uniform rings.

The ring profile describes how each ring vibrates across its width. Every ring is split into `n_sub` sub-rings with a complex weight that is the **same at every focal depth** (a property of the hardware), while each depth keeps only what the driver controls: one phase and gain per channel. `fit_ring_profile` fits the weights jointly over all measured depths, with a smoothness penalty between neighbouring sub-rings, and the near field included from 4 mm. Each ring's area-weighted mean weight is 1, so its level and phase stay with the channel drive.

| | Uniform rings | Ring profile, depth held out |
|---|---|---|
| CTX-500-026 (4 rings), median R² (worst) | 0.957 (0.79) | 0.989 (0.93) |
| CTX-250-009 (2 rings) | 0.894 (0.87) | 0.980 (0.95) |
| CTX-250-001 (4 rings) | 0.978 (0.96) | 0.980 (0.96) |
| CTX-500-026 near-field / focal ratio error | 7× | 33 % |

(Exploratory analysis, 2026-10-08; peak-normalised intensity from 4 mm, held-out depth with its channel drive refitted. "Uniform rings" here is an unregularised per-depth fit, not the shipped library fit, which regularises the ring gains and penalises a near field hotter than measured; see the off-axis comparison below. The benchmark `code/diagnostics/calib_ring_profile` repeats this for all combos, including prediction between measured depths.)

Use it from:

- **`calibration_gui`** and the Nexus *Calibration* page: **Ring profile** → *fit (shared by all depths)*; the GUI also shows the fitted profile (`plot_ring_profile`).
- **`refine_transducer_library`** / **`calibrate_library_entry`**: `RingProfile = 'fit'` (`RingProfileOptions`: `NSub`, `Smooth`, `NearFieldSkipMm`, `HeldOut`); with `Plot`, the profile figure is saved next to the entry.
- **`config_calibration.yaml`** (`calibration_transducer`): `ring_profile: fit`, `ring_profile_n_sub`, `ring_profile_smooth`, `ring_profile_near_field_skip_mm`.
- **`build_transducer_library`**: `RING_PROFILE = 'fit'` or a list of combos (default `'none'`).

`'auto'` (the default everywhere) keeps an entry's profile if it has one; `'none'` removes it.

**What changes.** The profile is stored in `global_model.ring_profile` (`n_sub`, `weight_abs`, `weight_phase_deg` per sub-ring, `smooth`, `near_field_skip_mm`, fit quality). All per-element arrays of the entry then refer to the sub-rings, so adding, removing or refitting a profile drops the entry's refinements (with a warning) and only the depths calibrated then are stored. Every function that combines an entry with the equipment geometry resolves it through `library_transducer`, so the profile is used consistently: the analytical models apply the weights in `axial_transfer_matrix`, and `load_transducer_from_library` hands the simulation an annular source of sub-rings, each driven with its channel's phase and drive times its weight (`ring_profile_drive`; `transducer.tx_calibration.ring_profile` records it). Sub-rings are steered with their ring (`annular.steer_r_mid_mm`). Inline `elem_phase_deg` / `elem_p_amp` in a study config refer to the equipment rings and bypass the profile; dephased conditions are not supported with a profile.

**Limitation and off-axis data.** On-axis profiles do not determine the surface vibration uniquely: refitting a known smooth profile from synthetic axial profiles reproduces the profiles almost perfectly (R² > 0.998) but not the weights. Off-axis measurements constrain it. `fit_ring_profile` accepts them through `OffAxis` (points r, z with intensity at a focal setting, e.g. the azimuthally averaged radial profile of a transverse hydrophone plane), modelled with `field_transfer_matrix` (the Rayleigh integral at points off the axis) and fitted together with the axial profiles. Validation on CTX-500-026 hydrophone planes at the 40 mm setting (`code/report/01_calibration/ring_profile_validation_*`, project repository; R² of intensity normalised by the axial peak, planes held out of the fit):

| Model (fitted to 12 axial profiles) | Axis (40 mm) | Sagittal XZ (z < 25 mm) | XY 35 mm | XY 40 mm (focus) | −3 dB width |
|---|---|---|---|---|---|
| Shipped library fit (uniform rings, `rayleigh_phase_gain`) | 0.87 | 0.89 (0.26) | 0.96 | 0.97 | 2.7 mm |
| Same, gain regularisation 10× weaker (`LambdaAmp` 1e-3) | 0.89 | 0.90 (0.33) | 0.94 | 0.96 | 2.6 mm |
| Uniform rings, unregularised | 0.90 | 0.86 (0.31) | 0.90 | 0.90 | 2.4 mm |
| Ring profile, axial | 0.97 | 0.88 (0.64) | 0.88 | 0.87 | 2.4 mm |
| Ring profile, axial + XY 35 mm plane | 0.97 | 0.94 (0.65) | 0.95 (fitted) | 0.98 | 2.9 mm |
| Hydrophone | | | | | 2.8 mm |

The shipped library fit predicts the transverse planes and the focal width well: its gain regularisation keeps the ring gains near equal, which costs on-axis fit but keeps the off-axis field right (unregularised fits match the axis better and the planes worse). Its shortfall is the axis and the near field. A ring profile fitted to axial profiles only improves those but loses the transverse planes; with one transverse plane added, it improves the axis and the near field and keeps the focal plane. Without off-axis data, keep the library fit and treat the off-axis field of a ring profile as unvalidated.

#### Optional: free-water k-Wave drive calibration

The analytical fit (Rayleigh) sets the phases and the drive for each refined depth. `kwave_calibrate_drive` simulates every refined depth in free-water k-Wave exactly as a study would resolve it from the entry (serial, driving system, depth, target ISPPA; in-pipeline ISPPA targeting off), and rescales the stored drive so that the simulation hits the reference ISPPA (`elem_p_amp × √(I_ref / I_kWave)`, `amp_source: kwave`). It also reports the simulated ISPPA, R² and FLHM-centre error against the hydrophone profile. Phases are not changed: k-Wave calibrates the amplitude of the analytical fit and validates its phases; it is not a phase optimiser. (Former name: `kwave_check_refinement`, kept as a deprecated wrapper.)

The report also scores the **absolute** profile (`kwave_drive_metrics`), because the peak-normalised R² cannot show a wrong level: `isppa_err_pct` is how far the analytical drive misses the reference ISPPA in k-Wave, `i_at_target_ratio` the simulated intensity at the measured FLHM centre relative to the reference, `kwave_over_analytical` the k-Wave ISPPA over the Rayleigh ISPPA at the same drive (the Rayleigh ISPPA of a shipped drive equals the reference to within 2 %), and `nrmse_before` / `nrmse_after` the RMS difference from the hydrophone profile scaled to the reference ISPPA, before and after the drive correction (what remains after is shape). `sim_seconds` and `grid_points` record the cost. `GridDims` (`kwave_grid_dims` in `config_calibration.yaml`) selects the free-water grid: `'auto'` sizes it to the aperture and focus (`grid.default_dims: 'auto'`).

```matlab
refine_transducer_library('CTX_250_001_SC_105_010', 'my_profiles.csv', lib, struct( ...
    'KWaveDriveCalibration', true, 'KWaveDriveOptions', struct('SimConfig', 'my_study_config.yaml')));
```

**Consumption-time warning.** Until the drive calibration is run, an entry's `amp_source` stays `analytical` at every depth — the drive was never validated against a real k-Wave simulation, only against the O'Neil/Rayleigh model. `load_transducer_from_library` surfaces this at the point of use: it resolves `tx_calibration.amp_kwave_calibrated` (true only when the drive at the requested depth comes wholly from a k-Wave-calibrated refinement (`amp_source: kwave`), not the always-analytical global model or a blend with it) and, if false, issues `load_transducer_from_library:uncalibratedAmplitude` naming the nearest checked depth (if any) and its R²/focal error for context. (The former names `amp_kwave_verified`, `amp_kwave_ref_depth_mm` are still written as deprecated aliases of `amp_kwave_calibrated`, `amp_kwave_nearest_depth_mm`; the warning was formerly `:unverifiedAmplitude`.) This is independent of BabelBrain's own `water_amp_calibration` (see `doc_babelbrain.md`), which corrects a *different* bias (BabelBrain's FDTD vs. its own Rayleigh reference) and runs automatically on every BabelBrain simulation regardless of whether the library entry itself was k-Wave-checked.

In the GUI, tick **k-Wave drive: calibrate in free water, after fit** (per-depth mode); the Nexus's Calibration page instead writes the unsaved fit and a `kwave_drive_<combo>.m` script, to run on its Run tab or on the cluster. The simulated profiles are overlaid in green, and **Save entry** stores the corrected drives. `SimConfig` supplies grid, `code_type` and `platform`/HPC settings. Without it, PRESTUS runs locally with `config_default` (CPU on macOS). Each depth is one full free-water simulation. On the default free-water grid (`grid.default_dims: 'auto'`, sized to the transducer) that takes about 3–20 min per depth in 3D on 16 CPU cores, and about 1 min with the **axisymmetric** grid option.

**Why run it.** The shipped drive is fitted with an analytical (Rayleigh) model. In a free-water k-Wave simulation it overshoots the reference ISPPA, consistently per transducer and almost independently of depth: by about 3 % (IGT 32-channel, 10 rings), 6–9 % (CTX-250) and 14–16 % (CTX-500) (k-Wave drive benchmark, 2026-10; three depths each). A study inherits that error in every ISPPA-dependent output (intensity, MI, heating) whenever the free-water baseline is skipped: with no `transducer.target_isppa_wcm2`, or with a target equal to the entry's calibration ISPPA (within 1 %; the default target of 30 W/cm² is the shipped reference), because PRESTUS then trusts the entry's drive ([ISPPA targeting](doc_simulations-acoustic.md#isppa-targeting)). A target that differs, or `modules.run_water_baseline: 1`, rescales to a k-Wave water simulation and removes the error, at the cost of one extra simulation on the study grid in every run; the drive calibration removes it once per entry. Calibrating in 3D removes it; calibrating on the axisymmetric grid leaves at most about 3 % (2D and 3D agree within 0.98–1.03 in ISPPA, with the same focus). Use 3D for final entries and the axisymmetric grid for quick checks or when ±3 % is acceptable. (Until 2026-10 the axisymmetric grid placed the rings of an annular array at the apex plane instead of on the bowl, which put the focus tens of mm off; fixed, see the changelog.)

`calibration_transducer` (`examples/calibration_standalone.m`, config-driven) runs the same calibration (`calibrate_library_entry`) and writes the library entry directly, plus a calibration output folder; see [doc_calibration.md](doc_calibration.md#calibration-at-chosen-depths-and-intensities-calibration_transducer). Per-run YAMLs from older versions (`{combo}-F*mm-I*wpercm2.yaml`) can still be ingested as refinements with:

```matlab
update_transducer_library(combo_name, path_output_profiles, my_library_path)
```

### Adding transducers and rebuilding entries

Adding a new transducer or driving system (including `calibration_gui`'s **New…** form), rebuilding a shipped entry and contributing entries back are described in [Adding and maintaining library entries](doc_transducer_library_maintenance.md).

### Predictions, exports and headless refinement

#### Predictions at desired depths

**Predict** models the profiles at the depths in **Desired depths**, for example `20:5:60` or `30, 42.5, 50` (default `min_foc:5:max_foc` in 5 mm steps). The depths need not have been measured. Each model is resolved the way a study resolves it (`evaluate_library_entry` → `library_model_at_depth`, the depth logic of `load_transducer_from_library`):

- **Geometric**: geometric steering only.
- **Library**: the library entry, including any refinements (the entry in **Save to folder**, if there is one, otherwise the shipped entry).
- **Fit**: the current fit, before it is saved, with k-Wave-calibrated drives if the k-Wave drive calibration ran.

Between refined depths the refinement offsets are interpolated. Outside the refined span they fade back to the global model over `blend_mm`. Drives and per-ring gains are interpolated in depth and held constant outside the span. The diagnostic `fit_z_shift_mm` is not applied, because a study does not apply it either.

**Profiles at desired depths.** With a profile CSV loaded, each desired depth also gets a profile from the measurements (`interpolate_axial_profiles`):

- **Measured**: the desired depth was measured.
- **Interpolated**: the profiles at the two neighbouring measured depths are shifted so their foci coincide at the interpolated focal position, then blended linearly by distance.
- **Extrapolated**: a desired depth outside the measured range, but within **extrapolate ≤** (default 0 mm, i.e. off). The nearest measured profile is shifted to the desired depth.
- None: further outside the range.

**Desired profiles** sets what the foci are aligned on. `label` (default) uses the column focal-depth headers, which are the nominal setpoints (≈ FLHM centre; see [focal-depth labels](#focal-depth-labels)). `peak` uses the intensity peak beyond the near-field skip. The near-field structure moves with the focus, so it is not physical; fits skip it.

The GUI plots the modelled axial profiles against the measured profile (solid) or the interpolated/extrapolated one (dashed). With **Show depth** set to one depth, the legend gives R² against that profile and the phase model; **All depths** overlays all of them. The lower panel plots the focus (FLHM centre beyond the near-field skip) against the target depth, with the measured depth range shaded. Hydrophone foci are drawn as filled stars (measured) or open stars (derived). The table lists, per desired depth: the profile type, the phase model (`refined`, `blend w=…`, `static` or `geometric`), the drive at the reference ISPPA, the predicted focus, R² and the element phases. These are the fit's values once there is a fit, otherwise the library's.

**Calibrating at desired depths.** Set **Calibrate at** to *Desired depths* and press **Fit**. The GUI then calibrates at every desired depth that has a profile (measured, interpolated or extrapolated), using the fit mode and settings above. Depths without a profile are skipped and reported. The result appears in the prediction view as a study would resolve it. In-sample fit R² is reported in the status, and the table shows R² through the runtime path. With *Fit depths (measured)*, **Fit** uses the measured depths selected under **Fit depths**, as before.

**Exports.**

- **Save entry**: writes each refinement's profile type to `refinement.profile_type`. If any depth was calibrated on a derived profile, it also writes the profiles fitted to as `{combo}_calibration_profiles.csv` next to the entry and records the path in `refinement.calibration_profiles_csv`.
- **Export profiles…**: writes the profiles of the last fit, or else those at the desired depths on view.
- **Export settings…**: writes the per-depth settings of the fit (else the library entry) at the depths on view, via `library_settings_table`. The columns are profile type, phase model, weight, drive, focus, and per ring phase, drive and gain.

Profile CSVs use the same format that `read_axial_profiles` reads, with a YAML sidecar. The sidecar records, for every derived column, the measured depths it was built from, their weights and shifts. The option `ExportFolder` writes both exports without a file dialog.

The predicted profile shape is exact for the chosen forward model: it depends only on geometry, frequency and the resolved phases and gains, while measured data only set the amplitude scale. When no profile CSV is found for the selected combination, the GUI opens in the prediction view (model only, no R², no fitting).

**Headless.** `refine_transducer_library` does the same from a script. It accepts depths that were not measured in `Depths`: a depth within 0.5 mm of a measurement uses it, and any other is interpolated (`InterpolateOptions`: `Align`, `MaxExtrapolationMm`). It writes the calibration profiles (`ExportProfiles`, on by default when any depth is derived) and optionally the settings table (`ExportSettings`, `SettingsDepths`):

```matlab
refine_transducer_library('CTX_250_001_SC_105_010', 'my_profiles.csv', lib, struct( ...
    'Depths', [30 42.5 55], 'ExportSettings', true, ...
    'InterpolateOptions', struct('MaxExtrapolationMm', 3)));
% writes {combo}.yaml, {combo}_calibration_profiles.csv (+ .yaml sidecar), {combo}_settings.csv
T = library_settings_table(tran, lib_struct, 20:5:70, struct('IsppaWcm2', 30, 'CsvPath', 'settings.csv'));
```

Changing the ring geometry of a transducer invalidates its library entries. Rebuild them.

---

## Dephased conditions (experimental)

> **Experimental — not validated.** No dephased hydrophone measurements have been compared with these predictions yet. Every use prints a `library_condition:experimental` warning. Do not report dephased-condition results without your own validation.

Dephased (defocused) settings drive the rings with explicit phases (e.g. alternate rings shifted by 180°) instead of steering to a focal depth. They cannot be interpolated in depth, so they are stored in the same library entry as separate, named **conditions**:

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

A study selects a condition with `transducer.condition: geom_180` (with `serial`, `combo.ds_serial` and `target_isppa_wcm2` as usual). `focal_distance_ep` may be omitted; it then defaults to `setpoint_ep_mm`.

**How phases and drive are derived.** In `shift` mode the phases are the geometric steering phases at the setpoint plus `phase_shift_deg`, which keeps the dephasing angles explicit. In `absolute` mode they are the given driver phases. With `inherit_correction` the entry's hardware phase offset and gains at the setpoint are added, assuming the correction fitted to focused profiles transfers to dephased drives. That assumption is what a dephased measurement has to confirm. The drive is set so that the Rayleigh free-water on-axis peak beyond 10 mm reaches `isppa_wcm2`. Measured profiles are only used to score the prediction (`predicted_r2`, focal and peak error), never to fit it.

**Adding conditions.** Without measurements:

```matlab
equip = load_equipment_config();
tran  = equip.trans.CTX_250_001;
lib   = yaml.loadFile('config/calibration/CTX_250_001_SC_203_035.yaml', 'ConvertToArray', true);
cond  = struct('label', 'geom_180', 'phase_shift_deg', [0 180 0 180], 'setpoint_ep_mm', 52.38);
lib   = library_add_conditions(lib, tran, cond, struct('IsppaWcm2', 30));
yaml.dumpFile('my_library/CTX_250_001_SC_203_035.yaml', lib);   % then transducer.library_path: my_library
```

With measurements, store them as `tran-{TRAN}_ds-{DS}_acq-dephased_axial.csv` in the usual CSV layout. The depth row holds each column's nominal setpoint. The sidecar describes the columns in order:

```yaml
columns:
  - label: geom_180
    phase_shift_deg: [0, 180, 0, 180]
  - label: geom_120
    elem_phase_deg: [0, 120, 240, 0]     # absolute driver phases
```

`library_ingest_conditions(lib_file, csv_path)` predicts every column, prints the agreement table and writes the conditions to `lib_file` (`struct('DryRun', true)` to only score). Dephased CSVs must not be passed to `build_transducer_library` or `refine_transducer_library`, which read every column as a focal depth.

**Before treating a condition as validated:** compare `predicted_r2` and the focal/peak errors with those of the focused depths of the same entry, check the location of `peak_z_ep_mm` against the measured peak, and confirm the free-water ISPPA in a k-Wave water run. Then set `validated: true` in your copy of the entry.

---

## Steering correction of matrix arrays

A matrix (2D) array is steered with geometric phases at one element amplitude. The focal peak it delivers can differ from the model by a factor that changes with the steering target: element directivity, crosstalk and manufacturing spread are not in the model. The steering correction measures that factor at a set of targets and scales the drive by it.

**What is fitted.** At each measured target x, `k(x) = measured / modelled focal peak`, both at the same drive. The modelled peak is the largest free-field pressure in the plane through the target, normal to the beam axis, from a Rayleigh integral over the array (`matrix_focal_peak`, `matrix_rayleigh_pressure`); the measured peak is the same quantity from a hydrophone scan of that plane. `fit_steering_correction` fits k as one constant (`scale`) or as a full quadratic surface in x, y, z (`quadratic`, 10 coefficients, linear least squares) and scores both by leave-one-out cross-validation. The default `auto` keeps the quadratic surface only with at least 15 targets and a median cross-validated error at most 0.8 times that of the constant. Arrays differ: on some the steering-dependent part is real, on others one constant is as good, and the cross-validation decides. A quadratic map is not extrapolated: targets outside the measured range use the value at its edge, with a warning.

**How it is applied.** `create_matrix_karray` multiplies every element amplitude by k(target) and records it in `transducer.matrix.steering_k`. `elem_p_amp` therefore stays the nominal drive setting and the simulated field is what the device delivers at that setting. With [ISPPA targeting](doc_simulations-acoustic.md#isppa-targeting) the field is reported at the target free-water ISPPA either way, so the correction does not change the results; it changes the drive that achieves them: `elem_p_amp × isppa_scale_p` (the scale recorded with the results), or in volts `drive_ref_v × elem_p_amp × isppa_scale_p / elem_p_amp_ref`. With `target_isppa_wcm2: []` the results are reported at the corrected drive. Switch it off with `transducer.matrix.use_steering_correction: false`.

**Building an entry.** Measure the focal peak at targets spread through the steering range you will use (30–40 targets, including its edges and the geometric focus, plus a few kept back to check the fit), all at one drive voltage and with the geometric phases PRESTUS computes for each target. Write one CSV row per target:

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

`resolve_transducer_from_serial` attaches the entry to a matrix transducer given by `serial` (`<serial>_<ds_serial>.yaml` first, then `<serial>.yaml`, or `transducer.library_path`); without an entry the array is simulated uncorrected. The entry stores the fit, its cross-validated errors, `elem_p_amp_ref` and `drive_ref_v`. An inline `transducer.matrix.steering_correction` takes precedence.

**Limits.** Only the focal peak is corrected, not side lobes or the shape of the field. The model assumes the beam along +z of the transducer frame and does not support Clover setups. The Rayleigh model treats each element as a baffled piston in free water; whether k-Wave reproduces it for a given array can be checked with a free-water run at a few targets. Background: Dadgar-Kiani et al. (2026), *Calibrated simulations for dynamic focusing of ultrasound through the temporal window*, Cell Reports Methods 6:101568, fit such a correction with a Gaussian of the difference in MPa; PRESTUS fits the ratio, which applies to the drive at any level.

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

---

## Appendix: Input file formats

`refine_transducer_library` (Level 2) and `build_transducer_library` (Level 3) need one axial profile CSV per transducer–driving-system combo. Steering phases are always computed from the ring geometry, so manufacturer phase tables are not needed.

---

### Axial profile CSV

#### Layout

```
,Focal depths (mm),,,,,,,
Distance WRT Exit plane (mm),40.0,50.0,60.0,70.0,80.0
0.5,0.12,0.31,0.04,0.01,0.0
1.0,0.18,0.44,0.08,0.02,0.0
...
```

**Row 1 (title row):** Label strings only; skipped.

**Row 2 (depth row):** First cell is a label string (ignored). Remaining cells give the focal depth of each column in mm from the transducer exit plane.

**Row 3 onward (data rows):** First column is axial distance from exit plane in mm. Remaining columns are measured intensity in W/cm² at each focal depth.

#### Rules

- Separator: comma
- Encoding: UTF-8 (BOM is tolerated — `readmatrix` strips it automatically)
- Focal depths: any number of columns ≥ 1; more depths improve global model accuracy
- Axial distance: must start at or below the face of the transducer (0 mm is fine; negative values are ignored)
- Intensity values: must be non-negative; `readmatrix` pads ragged CSVs with NaN — trailing NaN columns are trimmed automatically
- No units row: intensity is always interpreted as W/cm²

#### Reading in MATLAB

```matlab
raw = readmatrix(csv_path, 'NumHeaderLines', 1);
focal_depths_ep = round(raw(1, 2:end), 2);   % row 1, cols 2+
dist_ep_mm      = raw(2:end, 1);             % col 1, rows 2+
intens_data     = raw(2:end, 2:end);         % data block

% Drop NaN columns (ragged padding)
valid = isfinite(focal_depths_ep);
focal_depths_ep = focal_depths_ep(valid);
intens_data     = intens_data(:, valid);
```

#### Profile folder and file naming

Measured profiles are site data and live in `local/profiles/` inside the PRESTUS root. The folder is ignored by git. There is one subfolder per transducer, and all measurements of that transducer sit side by side, whatever driving system recorded them:

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

File names follow `tran-{TRAN_SERIAL}_ds-{DS_SERIAL}[_acq-{LABEL}]_axial.csv`, using the serials of `config/equipment/` (underscores, not hyphens). `acq-` marks a repeat or alternative measurement of the same pair.

`char_data_path` in the equipment YAML is relative to this folder, for example `CTX_250_001/tran-CTX_250_001_ds-SC_105_010_axial.csv`. Point `calibration.path_input_axial` (for `calibration_setup`) or `CHAR_DATA_ROOT` (for `build_transducer_library`) at the folder.

**Sidecar.** Each CSV has a YAML sidecar with the same stem. It records what the numbers mean: `quantity` and `units`, `distance_reference` (`exit_plane` | `bowl`), whether the column headers are `nominal` setpoints or `measured` positions (`focal_label_is`), the FLHM convention of the headers (`focal_label_reference`), the acquisition conditions and the file's provenance. PRESTUS does not read sidecars yet. Until it does, `calibration.add_FDO` and `calibration.focal_reference` must still be set in the calibration config to match the sidecar.

<a id="focal-depth-labels"></a>**Focal-depth labels.** Column headers are nominal driver setpoints. They correspond to the FLHM centre of the measured profile, but the threshold depends on the manufacturer. Checked on the Donders profiles (2026-09-30; median |header − FLHM centre|, first 10 mm skipped):

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
| Any | — | Cover full intended depth range; include the natural focal depth |

---

### Phase tables

Not used. Earlier versions read manufacturer steering tables (Imasonic `.ini`, Sonic Concepts `.csv`). PRESTUS now derives steering phases from ring geometry (`set_real_phases`), and any systematic difference from the hardware is absorbed by the fitted offset.

---

### Checklist before refining or rebuilding

- [ ] Axial profile CSV is readable by `readmatrix(path, 'NumHeaderLines', 1)`
- [ ] Row 2 contains focal depth values from column 2 onward (row 1 is a title row)
- [ ] The CSV is in `local/profiles/{TRAN_SERIAL}/`, named `tran-…_ds-…_axial.csv`, with a filled-in sidecar
- [ ] Column 1 contains axial distances in ascending order
- [ ] Intensity values are non-negative
- [ ] Equipment YAML validates: `check_equipment_config(equip_path)` passes
- [ ] Focal depths in CSV are within `[min_foc, max_foc]` for the transducer
