# Transducer Library

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
- **Global model**: geometric steering with the measured drive per depth (`global_model.amplitude_by_depth`) for per-depth entries. With `LIBRARY_MODE = 'static'` (debugging) it is instead one depth-independent static model, a per-ring phase offset and gain fitted jointly to all measured depths.

`transducer.tx_calibration.phase_model` records what was used: `geometric`, `static`, `refined` or `blend`; `amplitude_source` records `refinement`, `blend` or `global`.

**Amplitude.** Hydrophone profiles are recorded at a matched focal ISPPA per depth. The fit therefore uses profile *shape* only, with a free scale per depth, and derives the drive that produced the measured peak at each depth. For a refinement this is the drive with its own refined phases and gains (`refinement.elem_p_amp` at `refinement.isppa_wcm2`); for the global model it is `global_model.amplitude_by_depth`. Both are normalised to the reference ISPPA, interpolated in depth and scaled as √I at load time. Entries without either use the depth-independent `amplitude_scaling`. An optional free-water k-Wave drive calibration (Level 2) rescales the refinement drive to the simulated ISPPA.

### How the model is fitted

`fit_transducer_model` uses the fact that the on-axis field is linear in the complex element drives: p(z) = H(z)·w. The transfer matrix H is computed once per transducer (`axial_transfer_matrix`, Rayleigh by default), so each evaluation is a matrix product. Details:

- **Free parameters** (`FreeParams`): `phase`, `amp` (per-ring gain), `z_shift` (axial registration of the hydrophone profiles; a fit nuisance, not applied at run time) and `c_steer`. Rings driven by one channel share phase and gain; for example, the 30 sub-rings of IS_PCD15287 fit as 10 channels.
- **Loss**: pseudo-Huber on peak-normalised *pressure* amplitude beyond `NearFieldSkipMm`, with every depth weighted equally. Secondary lobes and split foci therefore count, and isolated artefacts don't dominate.
- **Optimiser**: `lsqnonlin` from `NStarts` starts (geometry first, then random), because profiles with two foci have several local optima. `lsqnonlin` needs the Optimization Toolbox; loading library entries does not.
- **Near field**: the first `NearFieldSkipMm` (10 mm) are not fitted. A soft penalty (`NearFieldPenalty`) discourages a model near field hotter than the measured one, and every fit reports `nearfield_excess`, because a hot spot near the exit plane matters in skull simulations.
- **Plausibility**: gains may range over 0.1–10× (`AmpBounds`). A gain outside 0.5–2× the mean (`PlausibleGain`), or at a bound, is flagged (`model.flags`) instead of being shipped silently: it more likely indicates a faulty element, cable or driver channel than something to calibrate away.

**Static library model (`LIBRARY_MODE = 'static'`, debugging).** Per-ring phase offset + per-ring gain (`FreeParams = {'phase','amp'}`). The axial shift and steering sound speed are fitted only as diagnostics. In the first benchmarked combo they did not improve held-out depths; confirm across all combos before changing this. A persistent axial shift points to an error in the equipment geometry or exit-plane distance and should be corrected there. Held-out results are recorded per entry in `meta.quality`.

**Library model version.** Entries with `meta.library_model_version: 2` apply per-ring gains to `elem_p_amp`, so simulated pressure and heating differ from runs with version-1 entries (phase offset only). Version-3 entries (`library_mode: per_depth`) take phases, gains and drive from the refinement at each measured depth. Results are not comparable across versions; record the version with your results (`transducer.tx_calibration.library_meta`).

### Library YAML schema

```yaml
meta:
  tran_serial: CTX_250_001
  driving_system: SC_105_010
  calibration_type: empirical        # geometric | empirical | user_empirical
  created_at: '2026-06-12'
  prestus_version: <commit>

global_model:                         # static model
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

Per-depth fits free the phase offsets and per-ring gains (`FreeParams = {'phase','amp'}`), with an intensity-domain loss (`Domain = 'intensity'`, `HuberDelta = 0.8`). They are regularised towards the entry's static model, which is geometric steering for per-depth entries. The same defaults apply in the GUI and in `build_transducer_library`. The axial registration `z_shift` stays at 0: a study does not apply it, so phases fitted with it would match the hydrophone only after a shift the entry never gets. The GUI therefore does not offer it. `build_transducer_library` still fits it once, separately, as a diagnostic (`z_shift_mm` in its quality report); a shift at its bound (±3 mm) at many depths points to an error in the equipment geometry or exit-plane distance, which should be corrected in the profiles or the equipment config rather than fitted.

You can also use the GUI:

```matlab
calibration_gui(struct('OutputPath', '/project/my_study/transducer_library'))
```

Before saving, **Predict** shows what the fit gives at the depths you will stimulate, including depths between and beyond the measured ones (see [Predictions at desired depths](#predictions-at-desired-depths)).

Fit settings are passed through `FitOptions` (see `help fit_transducer_model`). `Mode = 'static'` (GUI: **Static model (joint, debug)**) refits one depth-independent model over all selected depths instead. It is kept for debugging and comparison, not as the default.

#### Optional: free-water k-Wave drive calibration

The analytical fit (Rayleigh) sets the phases and the drive for each refined depth. `kwave_calibrate_drive` simulates every refined depth in free-water k-Wave exactly as a study would resolve it from the entry (serial, driving system, depth, target ISPPA; in-pipeline ISPPA targeting off), and rescales the stored drive so that the simulation hits the reference ISPPA (`elem_p_amp × √(I_ref / I_kWave)`, `amp_source: kwave`). It also reports the simulated ISPPA, R² and FLHM-centre error against the hydrophone profile. Phases are not changed: k-Wave calibrates the amplitude of the analytical fit and validates its phases; it is not a phase optimiser. (Former name: `kwave_check_refinement`, kept as a deprecated wrapper.)

```matlab
refine_transducer_library('CTX_250_001_SC_105_010', 'my_profiles.csv', lib, struct( ...
    'KWaveDriveCalibration', true, 'KWaveDriveOptions', struct('SimConfig', 'my_study_config.yaml')));
```

**Consumption-time warning.** Until the drive calibration is run, an entry's `amp_source` stays `analytical` at every depth — the drive was never validated against a real k-Wave simulation, only against the O'Neil/Rayleigh model. `load_transducer_from_library` surfaces this at the point of use: it resolves `tx_calibration.amp_kwave_calibrated` (true only when the drive at the requested depth comes wholly from a k-Wave-calibrated refinement (`amp_source: kwave`), not the always-analytical global model or a blend with it) and, if false, issues `load_transducer_from_library:uncalibratedAmplitude` naming the nearest checked depth (if any) and its R²/focal error for context. (The former names `amp_kwave_verified`, `amp_kwave_ref_depth_mm` are still written as deprecated aliases of `amp_kwave_calibrated`, `amp_kwave_nearest_depth_mm`; the warning was formerly `:unverifiedAmplitude`.) This is independent of BabelBrain's own `water_amp_calibration` (see `doc_babelbrain.md`), which corrects a *different* bias (BabelBrain's FDTD vs. its own Rayleigh reference) and runs automatically on every BabelBrain simulation regardless of whether the library entry itself was k-Wave-checked.

In the GUI, tick **k-Wave drive: calibrate in free water, after fit** (per-depth mode). The simulated profiles are overlaid in green, and **Save entry** stores the corrected drives. `SimConfig` supplies grid, `code_type` and `platform`/HPC settings. Without it, PRESTUS runs locally with `config_default` (CPU on macOS). Each depth is one full free-water simulation; the 3D default grid takes about an hour per depth on a CPU, so use a GPU or the cluster. The **axisymmetric** grid option is currently not reliable: on the free-water grid the simulated focus lies far short of the target (tested 2026-09-30, CTX-500-026). Use 3D until that is resolved.

`calibration_transducer` (`examples/calibration_standalone.m`, config-driven) runs the same calibration (`calibrate_library_entry`) and writes the library entry directly, plus a calibration output folder; see [doc_calibration.md](doc_calibration.md#calibration-at-chosen-depths-and-intensities-calibration_transducer). Per-run YAMLs from older versions (`{combo}-F*mm-I*wpercm2.yaml`) can still be ingested as refinements with:

```matlab
update_transducer_library(combo_name, path_output_profiles, my_library_path)
```

### Level 3: Rebuild a shipped entry

Maintainers only. Set `char_data_path` in `config/equipment/{serial}.yaml` and run `build_transducer_library`. `CHAR_DATA_ROOT` defaults to `local/profiles/` (see Appendix), and `LIBRARY_PATH` to `config/calibration/`.

By default (`LIBRARY_MODE = 'per_depth'`) every measured depth becomes a refinement (`PER_DEPTH_FIT_OPTIONS`, as in Level 2), on a geometric global model with the measured drive per depth. `meta.quality` records:
- the in-sample R² (median, min) and the geometric R²;
- the **held-out interpolation** R²: each interior depth predicted from its two neighbours, without axial shift, as at run time.

Held-out interpolation shows how much of the refinement carries over between measured depths. On CTX-250-001 + TPO-203-035 it is 0.938, against 0.970 in-sample and 0.849 for geometric steering. Space the measured depths closely where you will stimulate.

`LIBRARY_MODE = 'static'` (debugging) fits the static model with `FIT_OPTIONS` instead, behind an **acceptance gate**: fitted leave-one-depth-out, it must beat geometric steering on the held-out depths (median R² at least geometric, and no depth worse by more than `GATE_WORST` = 0.05), or a geometric entry is written. A missing profile CSV is an error in both modes; it does not fall back silently to a geometry-only entry. Then check the entry:

```matlab
validate_transducer_library(library_path, equip_path)   % must report no 'error'
```

### Level 4: Contribute back

Open a PR with:
- the rebuilt YAML (with `global_model.forward_model` set)
- `validate_transducer_library` passing with no errors
- the `meta.quality` summary (in-sample, held-out interpolation and geometric R²)
- characterisation conditions in the combo YAML `notes` field (water temperature, hydrophone model, tank geometry)

---

## Adding a new transducer or driving system

### New driving system for an existing transducer

1. Add a combo entry to `config/equipment/{tran_serial}.yaml`:
   ```yaml
   combos:
     MY_NEW_DS:
       ds:
         serial: MY_NEW_DS
         name: My Driving System
         manufact: MyCompany
       # char_data_path: add when hydrophone data is available
   ```
2. Generate a geometric baseline (Level 1).
3. Use it in study configs with `combo.ds_serial: MY_NEW_DS`.
4. Refine (Level 2) when hydrophone data is available.

### New transducer (not in equipment config)

1. Create `config/equipment/{MY_SERIAL}.yaml`, copying an existing file as a template.

   Required fields for annular arrays:
   ```yaml
   type: transducer
   serial: MY_SERIAL
   name: Human-readable name
   manufact: Imasonic    # or: Sonic Concepts
   min_foc: 30.0         # minimum steerable focal distance from exit plane [mm]
   max_foc: 90.0
   transducer:
     annular:
       elem_n: 10
       elem_id_mm: [0.0, ...]
       elem_od_mm: [21.0, ...]
       curv_radius_mm: 100.0
       dist_geom_ep_mm: 92.7
     freq_hz: 300000.0
   ```
2. Run `check_equipment_config(equip_path)` to validate the geometry and its physical plausibility.
3. Continue from Level 1.

Alternatively, click **New…** under the transducer list in `calibration_gui`. The form takes the type (annular or matrix), frequency, curvature radius, exit-plane-to-geometric-focus distance (`dist_geom_ep_mm`) and focal range, then the type-specific geometry:

- **Annular:** ring inner and outer diameters and the number of physical channels, plus optionally one driving system with its profile CSV.
- **Matrix:** steering (1D/3D), curved or flat, aperture diameter, element shape and size, and a `define_here` grid (`rect` rows × columns with spacing, `fibonacci` count and kerf, or `fermat` count). Coordinate-file layouts (`extract_from_file`) are not offered in the form.

**Copy from** fills the form from an existing transducer. A side view and a face view update as you type: rings coloured by channel, or matrix element footprints, with the exit plane, geometric focus and steering range. The message box reports the ring/channel count or the element count and median pitch, and flags overlapping rings or elements. On save the GUI:

- rejects a serial that is not a valid identifier, unequal inner/outer diameter counts, a ring whose outer diameter does not exceed its inner diameter, a ring count that is not a multiple of the physical channel count (`fit_transducer_model` groups consecutive rings into channels), a curved matrix whose curvature radius is not larger than half the aperture, and a matrix grid with no element inside the aperture;
- runs `check_equipment_config` on the entry and asks for confirmation if it reports issues;
- writes `{serial}.yaml` to the GUI's equipment folder (`EquipPath`, default `config/equipment/`) and selects the new transducer. For an annular transducer, the geometric prediction is shown right away.

The profile CSV path is stored relative to `ProfileRoot` when it lies under it, otherwise as an absolute path. For matrix entries, the top-level `n_elem` is the element count inside the aperture. Serial-based study configs (`transducer.serial`) always read `config/equipment/`, so save there if the pipeline should resolve the new serial.

**View geometry** draws the selected transducer the same way. Matrix transducers are geometry only in the calibration GUI: the axial forward model behind prediction and fitting (`axial_transfer_matrix`, `steering_path_lengths`) is defined for annular arrays.

#### Predictions at desired depths

**Predict** models the profiles at the depths in **Desired depths**, for example `20:5:60` or `30, 42.5, 50` (default `min_foc:5:max_foc` in 5 mm steps). The depths need not have been measured. Each model is resolved the way a study resolves it (`evaluate_library_entry` → `library_model_at_depth`, the depth logic of `load_transducer_from_library`):

- **Geometric**: geometric steering only.
- **Library**: the library entry, including any refinements (the entry in **Save to folder**, if there is one, otherwise the shipped entry).
- **Fit**: the current fit, before it is saved, with k-Wave-calibrated drives if the k-Wave drive calibration ran.

Between refined depths the refinement offsets are interpolated. Outside the refined span they fade back to the static model over `blend_mm`. Drives and per-ring gains are interpolated in depth and held constant outside the span. The diagnostic `fit_z_shift_mm` is not applied, because a study does not apply it either.

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

**Focal-depth labels.** Column headers are nominal driver setpoints. They correspond to the FLHM centre of the measured profile, but the threshold depends on the manufacturer. Checked on the Donders profiles (2026-09-30; median |header − FLHM centre|, first 10 mm skipped):

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
