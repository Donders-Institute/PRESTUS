# Transducer Library

!!! note "At a glance"
    - **What:** The shipped library of transducers and their calibrations, and how a study uses an entry.
    - **When you need it:** Whenever you simulate a device from the library: set `serial`, `combo.ds_serial`, `focal_distance_ep` and `target_isppa_wcm2`.
    - **Before you start:** Nothing. To fit an entry to your own hydrophone profiles, see [Calibration](doc_calibration.md); to add a device, see [Adding & Maintaining Entries](doc_transducer_library_maintenance.md).

PRESTUS ships a static transducer library. You can refine an entry against your own hydrophone measurements ([Calibration](doc_calibration.md)).

| Folder | Content | Purpose |
|---|---|---|
| `config/equipment/` | One YAML per physical device | Ring geometry, frequency, steerable range |
| `config/calibration/` | One YAML per transducer–driving-system combo | Per-depth refinements (phases, gains, drive) on geometric steering |
| your folder (optional) | Same schema | Your own refinements; set `transducer.library_path` |

## Using an entry

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

To use an entry you refined yourself ([Calibration](doc_calibration.md#refining-a-library-entry-to-your-hydrophone-profiles)), add `library_path`:

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

    **Amplitude.** Hydrophone profiles are recorded at a matched focal ISPPA per depth. The fit therefore uses profile *shape* only, with a free scale per depth, and derives the drive that produced the measured peak. For a refinement this is `refinement.elem_p_amp` at `refinement.isppa_wcm2`; for the global model it is `global_model.amplitude_by_depth`. Both are normalised to the reference ISPPA, interpolated in depth and scaled as √I at load time. Entries without either use the depth-independent `amplitude_scaling`. The optional [free-water k-Wave drive calibration](doc_calibration.md#optional-free-water-k-wave-drive-calibration) rescales the refinement drive to the simulated ISPPA.

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

## Dephased conditions (experimental)

!!! warning "Experimental, not validated"
    No dephased hydrophone measurements have been compared with these predictions yet. Every use prints a `library_condition:experimental` warning. Do not report dephased-condition results without your own validation.

Dephased (defocused) settings drive the rings with explicit phases (e.g. alternate rings shifted by 180°) instead of steering to a focal depth. They cannot be interpolated in depth, so the library entry stores them as named **conditions**. Select one with `transducer.condition: geom_180` (with `serial`, `combo.ds_serial` and `target_isppa_wcm2`). Without `focal_distance_ep`, it defaults to the condition's `setpoint_ep_mm`.

**Before treating a condition as validated:** compare `predicted_r2` and the focal/peak errors with those of the focused depths of the same entry, check `peak_z_ep_mm` against the measured peak, and confirm the free-water ISPPA in a k-Wave water run. Then set `validated: true` in your copy of the entry.

??? note "Details: adding conditions, with or without measurements"
    **Without measurements:**

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
