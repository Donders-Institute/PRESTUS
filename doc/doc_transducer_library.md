# Transducer Library

![Pipeline steps; this page belongs to 2 Transducer](assets/figures/fig_stage-strip-transducer_diagram.svg){ .stage-strip }

A simulation reproduces a real device only if it uses that device's ring geometry, phases, and drive. PRESTUS ships a static library that holds these settings for each supported transducer and driving-system combination, fitted to free-water hydrophone axial profiles where profiles exist, so that a study names the device and a focal depth instead of entering element settings. You can refine an entry against your own hydrophone measurements ([Calibration](doc_calibration.md)).

On this page:

- [Using a Library Entry](#using-a-library-entry): selecting an entry, the lookup order, how phases and drive are computed, and the entry schema.
- [Dephased Conditions](#dephased-conditions): named defocused settings, which are experimental and not validated.
- [Library Limitations](#library-limitations): where the library is less accurate and how to mitigate it.

!!! note "At a glance"
    - **When you need it:** Whenever you simulate a device from the library: set `serial`, `combo.ds_serial`, `focal_distance_ep`, and `target_isppa_wcm2`.
    - **Before you start:** Nothing. To fit an entry to your own hydrophone profiles, see [Calibration](doc_calibration.md); to add a device, see [Adding & Maintaining Entries](doc_transducer_library_maintenance.md).

| Folder | Content | Purpose |
|---|---|---|
| `config/equipment/` | One YAML per physical device | Ring geometry, frequency, steerable range |
| `config/calibration/` | One YAML per transducer–driving-system combo | Per-depth refinements (phases, gains, drive) on geometric steering |
| your folder (optional) | Same schema | Your own refinements; set `transducer.library_path` |

## Using a Library Entry

Name the transducer and driving system, the focal depth and the target intensity:

```yaml
transducer:
  serial: IS_PCD15473_01001
  combo:
    ds_serial: IGT_32_ch_comb_10_ch
  target_isppa_wcm2: 32.0
  focal_distance_ep: 60.0
  # resolves config/calibration/model-IS_PCD15473/tran-IS_PCD15473_01001_ds-IGT_32_ch_comb_10_ch_calib.yaml
```

PRESTUS resolves geometry and calibrated parameters at load time: the ring geometry from `config/equipment/`, and the phases, per-ring gains, and drive at `focal_distance_ep` from the calibration entry. Fields already set in the study config override library values. [Shipped entries](doc_transducer_library_shipped.md#shipped-entries-and-fit-quality) lists the entries and how well each one fits its hydrophone profiles.

To use an entry you refined yourself ([Calibration](doc_calibration_hydrophone-profiles.md#refining-a-library-entry)), add `library_path`:

```yaml
transducer:
  serial: CTX_250_001
  combo:
    ds_serial: SC_105_010
  library_path: /project/my_study/transducer_library   # written by refine_transducer_library
  target_isppa_wcm2: 30.0
  focal_distance_ep: 35.0
```

The run records what was used: `transducer.tx_calibration.phase_model` (`geometric`, `static`, `refined`, or `blend`), `amplitude_source` (`refinement`, `blend`, or `global`), and `library_meta` (the library model version). Record the version with your results; results are not comparable across versions.

To bypass the library, give the full geometry, phases, and drive inline: see [Transducer definition](doc_transducer.md#transducer-fields).

??? note "Library keys and layout"
    Each entry has a key: `<TRAN>_<DS>` for a transducer with a driving system, `<TRAN>` for an entry without one (the geometric baselines of the device models). A library folder holds entries in either of two layouts, and PRESTUS reads both:

    | Layout | File | Used by |
    |---|---|---|
    | Nested | `model-<MODEL>/tran-<TRAN>_ds-<DS>_calib.yaml`, `model-<MODEL>/tran-<TRAN>_calib.yaml` | the shipped library (`config/calibration/`) |
    | Flat | `<TRAN>_<DS>.yaml`, `<TRAN>.yaml` | your own libraries (`transducer.library_path`) and libraries from earlier versions |

    `library_entry_file(library_path, key)` returns the file of an entry. PRESTUS resolves in order:

    1. `transducer.library_path` (your folder), if set and it contains an entry for this transducer
    2. the shipped entry `<TRAN>_<DS>`
    3. the shipped entry `<TRAN>`, with a warning
    4. If none exists, an error with instructions

    Without `combo.ds_serial`, PRESTUS looks only for the entry `<TRAN>`.

### Your Own Transducers and Calibrations

Keep your own devices, calibrations, and hydrophone profiles in a library folder of your own, outside the PRESTUS folder, so that updating PRESTUS never touches them. Do not add them to `config/equipment` or `config/calibration`. Name the folder in one of three places; PRESTUS checks them in this order, then its shipped library:

| Setting | Scope |
|---|---|
| `transducer.library_path` | This study (per transducer) |
| `path.transducer_library` | Every study of a project: set it in the project copy of `config_default.yaml` ([project defaults](doc_configuration.md#project-config-folder)). A folder or a list |
| `PRESTUS_TRANSDUCER_LIBRARY` | This machine or user: an environment variable, folders separated by `:` (`;` on Windows). `calibration_gui` and the Nexus read it too |

The first library that has an entry wins: an equipment entry replaces entries with the same serial from later libraries and from the shipped library (with a message), and a calibration entry for the driving system is taken from the first library that has one. The folder can hold:

| File | Content |
|---|---|
| `tran-<SERIAL>_transducer.yaml` | An equipment entry: a new unit, or a replacement for a shipped entry with the same serial |
| `ds-<DS>_driver.yaml` | A driving system, optionally with its drive limit |
| `<TRAN>_<DS>.yaml` (or `tran-<TRAN>_ds-<DS>_calib.yaml`) | A calibration entry for that transducer and driving system |
| `profiles/` | Hydrophone profiles; a relative `char_data_path` is found here first, then in `local/profiles/` |

A unit of a model in the library needs only its identity and combos; it takes geometry and steering range from the model entry:

```yaml
# my_library/tran-IS_PCD15473_01001_OPM_transducer.yaml
type: transducer
serial: IS_PCD15473_01001_OPM
model: IS_PCD15473
name: Imasonic 10 ch. IS_PCD15473_01001 (OPM setup)
combos:
  IGT_32_ch_comb_10_ch:
    char_data_path: tran-IS_PCD15473_01001_ds-IGT_32_ch_comb_10_ch_acq-OPM_axial.csv   # in my_library/profiles/
```

```yaml
# study config
transducer:
  serial: IS_PCD15473_01001_OPM
  combo: {ds_serial: IGT_32_ch_comb_10_ch}
  library_path: /path/to/my_library
```

Create the calibration entry from your hydrophone profiles with the Nexus Calibration sub-tab ([Calibration page](doc_nexus_tabs.md#calibration-page)), `calibration_gui`, `refine_transducer_library`, or `calibration_transducer`. `calibration_gui` saves entries and new transducers into its output folder, by default the first folder in `PRESTUS_TRANSDUCER_LIBRARY`; it does not write into `config/equipment` (option `AllowShipped` for maintainers), and an edited shipped transducer is saved as a replacement in your library. If your driving system gives you phases per depth instead, write a short entry by hand:

```yaml
# my_library/IS_PCD15473_01001_OPM_IGT_32_ch_comb_10_ch.yaml
meta: {tran_serial: IS_PCD15473_01001_OPM, driving_system: IGT_32_ch_comb_10_ch, calibration_type: user}
global_model:
  depths_ep_mm: [30, 50, 70]
  elem_phase_deg:                 # absolute phase of each ring at each depth [deg]
  - [0, 12, 25, 37, 50, 62, 75, 87, 100, 112]
  - [...]
  - [...]
  amplitude_by_depth: {depths_ep_mm: [30, 50, 70], elem_p_amp: [88000, 103000, 132000], isppa_wcm2: 30}
```

Between the given depths, PRESTUS interpolates the phases, and outside them it extrapolates and warns. Entries fitted by PRESTUS record the geometry they were fitted for (`meta.geometry_fingerprint`). When the equipment entry in use has another geometry, for example because your folder replaces a shipped entry, the entry still loads and PRESTUS warns (`load_transducer_from_library:geometryMismatch`). The Nexus reads the equipment entries in `transducer.library_path` as well and labels them "(own library)".

### How Phases Are Computed

At a requested depth, the entry gives geometric steering plus a fitted per-ring correction, and the drive that reaches the reference [I<sub>SPPA</sub>](doc_glossary.md#isppa).

??? note "Phase and drive at a requested depth"
    For any requested focal depth *d* (mm from the exit plane):

    ```
    phase(d) = geometric steering phase(d; c_steer) + offset(d)
    drive_i  = elem_p_amp(d) · gain_i · exp(i · phase_i(d))
    ```

    - **Geometric steering** (`set_real_phases`) comes from the ring geometry in `config/equipment/`, at `global_model.steering_sound_speed` (default 1500 m/s).
    - **Refinement** (default): one fit per measured depth (`refinement.depths_ep_mm`). Each holds a per-ring phase offset (`elem_phase_correction_deg`), per-ring gains (`elem_amp_weight`), and the drive that gives the reference I<sub>SPPA</sub> at that depth (`elem_p_amp`). Between refined depths, offsets are interpolated per element on the unit circle; gains and drive are interpolated linearly. Outside that span, the nearest refinement fades back to the global model over `refinement.blend_mm` (default 5 mm).
    - **Global model**: geometric steering with the measured drive per depth (`global_model.amplitude_by_depth`). Some entries hold one depth-independent offset and gains instead (`global_model.elem_phase_correction_deg`, `elem_amp_weight`). These are read but cannot be fitted.

    **Amplitude.** Hydrophone profiles are recorded at a matched focal I<sub>SPPA</sub> per depth. The fit therefore uses profile *shape* only, with a free scale per depth, and derives the drive that produced the measured peak. For a refinement this is `refinement.elem_p_amp` at `refinement.isppa_wcm2`; for the global model it is `global_model.amplitude_by_depth`. Both are normalised to the reference I<sub>SPPA</sub>, interpolated in depth, and scaled as √I at load time. Entries without either use the depth-independent `amplitude_scaling`. The optional [free-water k-Wave drive calibration](doc_calibration_hydrophone-profiles.md#free-water-k-wave-drive-calibration) rescales the refinement drive to the simulated I<sub>SPPA</sub>.

### Library YAML Schema

??? note "Entry schema"
    ```yaml
    meta:
      tran_serial: CTX_250_001
      driving_system: SC_105_010
      calibration_type: empirical        # geometric | empirical | user_empirical
      created_at: '2026-06-12'
      prestus_version: <version>
      prestus_commit: <commit>

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

The entries shipped with PRESTUS, and how well each fits its profiles, are listed in [Shipped Library Entries](doc_transducer_library_shipped.md).

## Dephased Conditions

!!! warning "Experimental, not validated"
    No dephased hydrophone measurements have been compared with these predictions yet. Every use prints a `library_condition:experimental` warning. Do not report dephased-condition results without your own validation.

Dephased (defocused) settings drive the rings with explicit phases (e.g. alternate rings shifted by 180°) instead of steering to a focal depth. They cannot be interpolated in depth, so the library entry stores them as named **conditions**. Select one with `transducer.condition: geom_180` (with `serial`, `combo.ds_serial`, and `target_isppa_wcm2`). Without `focal_distance_ep`, it defaults to the condition's `setpoint_ep_mm`.

**Before treating a condition as validated:** compare `predicted_r2` and the focal/peak errors with those of the focused depths of the same entry, check `peak_z_ep_mm` against the measured peak, and confirm the free-water I<sub>SPPA</sub> in a [k-Wave](https://github.com/ucl-bug/k-wave) water run. Then set `validated: true` in your copy of the entry.

??? note "Adding conditions, with or without measurements"
    **Without measurements:**

    ```matlab
    equip = load_equipment_config();
    tran  = equip.trans.CTX_250_001;
    lib   = yaml.loadFile(library_entry_file('config/calibration', 'CTX_250_001_SC_203_035'), 'ConvertToArray', true);
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

    `library_ingest_conditions(lib_file, csv_path)` predicts every column, prints the agreement table, and writes the conditions to `lib_file` (`struct('DryRun', true)` to only score). Do not pass dephased CSVs to `build_transducer_library` or `refine_transducer_library`, which read every column as a focal depth.

??? note "Condition schema and how phases and drive are derived"
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

## Library Limitations

| Limitation | Effect | Mitigation |
|---|---|---|
| Refinements are exact only at measured depths | Between depths, interpolated phases fit less well (held-out interpolation in `meta.quality`) | Measure closely spaced depths around the targets |
| Axial shift fitted, not applied | A real model/hydrophone offset remains in focal placement | Correct the equipment geometry if `fit_z_shift_mm` is persistent |
| Refinement is only as good as your profiles | Offsets inherit hydrophone/tank errors | Cover the intended depth range; record conditions |
| Water-tank basis | Does not account for tissue acoustic properties | Combined with the PRESTUS forward simulation in tissue |
| Profiles recorded at matched I<sub>SPPA</sub> per depth | Depth dependence of focal gain is not calibrated | Water-baseline I<sub>SPPA</sub> rescaling |
| DS-agnostic fallback | Ignores DS hardware phase offset | Use DS-specific entries |
| Outside the measured range | Extrapolated from the fitted model, warning emitted | Measure the depths you target |
| Outside `[min_foc, max_foc]` | Beyond the hardware steering range, warning emitted | Keep `focal_distance_ep` within range |
