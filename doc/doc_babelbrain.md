# BabelBrain Backend

!!! note "At a glance"
    - **What:** BabelBrain as an alternative acoustic solver, thermal solver, or both. k-Wave stays the default.
    - **When you need it:** To cross-check a k-Wave result with an independent solver, or to use the Apple Silicon GPU.
    - **Before you start:** The `external/BabelBrain` submodule and a BabelBrain conda environment ([Installation](#installation)).

PRESTUS can run [BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain) (Samuel Pichardo, University of Calgary) in place of k-Wave. Pre-processing (segmentation, grid, medium masks, placement) and post-processing (focal metrics, CEM43, NIfTI export, HTML report) run unchanged. BabelBrain's solver, BabelViscoFDTD, solves the viscoelastic wave equation and adds a Rayleigh integral step for the water coupling layer. PRESTUS runs it in fluid-bone mode (shear off, PRESTUS tissue properties injected), so results are comparable to k-Wave and jwave. On Apple Silicon it has a native Metal backend, where k-Wave only has the slower `matlab_gpu` path.

## Quick start

Set three things: the backend, the Python environment and (optionally) the transducer preset.

```yaml
acoustic_backend: babelbrain

babelbrain:
  python_env: /opt/conda/envs/BabelBrain   # directory containing bin/python; must be set
  tx_system: R15473                        # BabelBrain preset name; auto-detected if empty
  computing_backend: auto                  # 'auto' | 'cuda' | 'opencl' | 'metal' | 'cpu'  (or 1-4)
```

`babelbrain.python_env` must be set. Empty falls back to system `python3`, which fails with `No module named 'BabelViscoFDTD'`. PRESTUS finds BabelBrain at `external/BabelBrain` (git submodule); set `babelbrain.root` to override. For `tx_system`, see [Transducer mapping](#transducer-mapping).

To run the thermal simulation in BabelBrain as well:

```yaml
acoustic_backend: babelbrain
thermal_backend:  babelbrain   # BabelBrain Pennes BHTE instead of k-Wave; needs acoustic_backend: babelbrain
```

| `acoustic_backend` | `thermal_backend` | Scenario |
|---|---|---|
| `kwave` (default) | `kwave` (default) | Standard PRESTUS, no BabelBrain |
| `babelbrain` | `kwave` | BabelBrain acoustics, PRESTUS thermal |
| `babelbrain` | `babelbrain` | Full BabelBrain pipeline, PRESTUS post-processing |

`thermal_backend: babelbrain` reads the `*DataForSim.h5` file of the BabelBrain acoustic stage, so other combinations are rejected at config load. The thermal step has no CPU backend.

??? note "Details: all `babelbrain` keys"
    ```yaml
    babelbrain:
      tx_system: R15473            # BabelBrain preset name; auto-detected if empty
      computing_backend: auto      # 'auto' | 'cuda' | 'opencl' | 'metal' | 'cpu'  (or 1-4)
      device_name: auto            # 'auto' selects device based on computing_backend; or set explicitly e.g. 'Apple'
      ppw: 12                      # points-per-wavelength floor (default 12)
      do_refocusing: false         # skull back-propagation refocusing
      z_beyond_focus_mm: 40        # axial domain truncation past the focal point [mm]
      python_env: ''               # path to conda/venv environment; empty = system python3
      root: ''                     # path to BabelBrain source; auto-detected if empty
      tx_mech_adjust_z_m: 0.0      # ADDITIONAL transducer z offset into skin [m]; a placement
                                   # correction (preset FL − |focus_pos − trans_pos|) is always added
      z_steering_m: ~              # ZSteering [m]; null = auto-derived (0 with explicit phases)
      q_correction: []             # bone Q-factor correction of the relaxation attenuation model; [] = 1.0 non-CT / 3.0 CT
      focal_calibration: true              # geometric fallback only: steer so the achieved FWHM-centre lands at focal_distance_ep;
                                           # locate-only (report) when explicit element phases exist
      focal_calibration_target: flhm_center # quantity placed at focal_distance_ep: 'flhm_center' | 'peak'
      profile_calibration: false   # annular explicit phases only: refit BabelBrain phases to the O'Neil axial profile
      water_amp_calibration: true  # scale output so the free-water peak matches Rayleigh for the same drive
      use_rayleigh_for_water: false # water-only solver: false = full-wave FDTD, true = Rayleigh integral
      ct_mapping: ''               # HU-to-property mapping for CT-based skull; empty = label-based path
      air_regions: false           # pseudo-CT path only: internal air as pressure-release reflector mask
      air_hu_window: [-1200, -400] # HU range classified as air
      air_keep_boundary: false     # keep air components touching the domain boundary
      cleanup_cache: true          # delete the per-subject BabelBrain cache after reading
      dry_run: false               # prepare inputs only; stop after the acoustic stage
      # isppa_w_per_cm2: ~         # thermal only; leave empty: derived from the scaled acoustic field (in-situ brain peak)
    ```

    Placement correction, `z_steering_m` and `focal_calibration` are explained in [Transducer placement, ZSteering, and the drive phases](doc_babelbrain_internals.md#transducer-placement-zsteering-and-the-drive-phases).

    **PPW and resolution.** PRESTUS uses `max(ppw, round(c_min / (freq × res_mm × 1e-3)))`, where `c_min` is the smaller of BabelBrain's bone shear-speed floor (≈1102.5 m/s) and the lowest `medium_properties` sound speed. PRESTUS resamples the mask to BabelBrain's native step, so resolutions need not be matched by hand. Lower `ppw` runs faster with more axial numerical dispersion.

## Installation

Initialise the submodule after cloning PRESTUS:

```bash
git submodule update --init external/BabelBrain
```

BabelBrain needs its own Python environment, because its compiled extension `BabelViscoFDTD` does not install into a generic system Python. Use the bundled conda file for your platform:

```bash
# Apple Silicon (Metal)
conda env create -f external/BabelBrain/environment_mac_arm64-314.yml
# macOS Intel
conda env create -f external/BabelBrain/environment_mac_x64-314.yml
# Linux (CUDA)
conda env create -f external/BabelBrain/environment_linux.yml
# Windows
conda env create -f external/BabelBrain/environment_win-314.yml

conda activate BabelBrain
which python   # e.g. /opt/conda/envs/BabelBrain/bin/python
```

Put the environment directory (the one containing `bin/python`) into `babelbrain.python_env`. For GPU backend choice per platform, see [GPU backend selection](doc_babelbrain_internals.md#gpu-backend-selection-and-platform-support).

## How BabelBrain compares with k-Wave

| | PRESTUS (k-Wave) | BabelBrain |
|---|---|---|
| **Solver** | Pseudo-spectral FDTD | Explicit staggered-grid FDTD |
| **Physics** | Scalar linear acoustic (pressure/density) | Viscoelastic solver; run with shear = 0 (fluid bone) in PRESTUS |
| **Skull model** | Scalar ρ, c, α | Scalar ρ, c, α from PRESTUS `medium_properties` (or BabelBrain HU mapping via `ct_mapping`) |
| **Rayleigh pre-step** | None | Rayleigh integral for the water coupling layer |
| **Narrow-beam domain** | Optional (see [limited FOV](doc_limited_fov.md)) | Off by default (the pre-focal field would not be returned); enable with `babelbrain.tight_narrow_beam: true`. Always off when `grid.acoustic_fov_length_mm` is set (mask already cropped) |
| **GPU backends** | CUDA, MATLAB GPU arrays | CUDA, OpenCL, Metal (Apple Silicon) |

See [k-Wave](doc_kwave.md) for the k-Wave backend options.

!!! warning "Tissue properties and comparability"
    On the default label-mask path, PRESTUS writes `medium_properties` to `<target>_medium_properties.json`. BabelBrain then uses PRESTUS density, sound speed and attenuation for skin, cortical bone, trabecular bone and brain (an undifferentiated `skull` fills both bone slots), with shear = 0. Remaining differences are the solver and the bone Q-correction (`q_correction`), so validate against a k-Wave run before relying on absolute values.

    PRESTUS properties do **not** apply to:

    - **water**: BabelBrain keeps 1000 kg/m³, 1500 m/s, no attenuation (PRESTUS default 994 kg/m³, 1500 m/s, 2.17e-3 dB/(cm·MHz²)). For identical water across backends, set `medium_properties.water` to 1000 kg/m³ / 1500 m/s.
    - **skull on the `ct_mapping` path**, where BabelBrain's HU models set bone (skin and brain are still injected).
    - runs where writing the JSON failed, or where the drive frequency has no entry in BabelBrain's table (5 kHz steps). These print a warning and fall back to [BabelBrain's built-in table](doc_babelbrain_internals.md#built-in-material-table).
    - **voxel-wise pseudo-CT properties** computed by PRESTUS.
    - **thermal** properties (see [Thermal tissue properties](#thermal-tissue-properties)).

## BabelBrain output folder

With the BabelBrain backend active, or with `io.export_babelbrain: true`, PRESTUS writes a `babelbrain/` folder next to the standard outputs (`nii/`, `log/`, `cache/`). It holds everything needed to re-run BabelBrain on its own:

| File | Description |
|---|---|
| `<target>_<freq>kHz_<PPW>PPW_BabelViscoInput.nii.gz` | Tissue mask in BabelBrain integer label format, including the focal voxel marker |
| `<target>_CalibWeights.h5` | Per-element calibration phases from PRESTUS calibration |
| `RunBabelBrain.py` | Ready-to-run Python script with all simulation parameters pre-filled |
| `babelbrain_config.yaml` | Parameter record for reproducibility |

To export the package while simulating with k-Wave, set `io.export_babelbrain: true` with `acoustic_backend: kwave` (or no `acoustic_backend` key). Use this to prepare a cross-validation or to run BabelBrain on another machine.

## Tissue labels and material properties

By default each tissue gets one set of properties from `medium_properties`, written as BabelBrain's integer label NIfTI (`BabelViscoInput.nii.gz`). Tissues without a BabelBrain label are written as water.

??? note "Details: label mapping"
    | PRESTUS layer name | Default config indices | BabelBrain label | BabelBrain material |
    |---|---|---|---|
    | `water` | 0, 3, 6, 9, 10 | 0 | Water |
    | `skin` | 5 | 1 | Skin |
    | `skull_cortical` | 7 | 2 | Cortical bone |
    | `skull_trabecular` | 8 | 3 | Trabecular bone |
    | `skull` (undifferentiated) | 4 | 2 | Cortical bone (fallback) |
    | `brain` | 1, 2 | 4 | Brain |
    | *(focal voxel)* | – | 5 | – |

    Label 5 marks the focus voxel. It is BabelBrain's only way to specify the target.

### CT-based skull properties (`ct_mapping`)

Set `babelbrain.ct_mapping` to let BabelBrain derive per-voxel skull properties from a pseudo-CT (HU) volume with one of its fitted models (BabelBrain's `bUseCT` mode):

```yaml
babelbrain:
  ct_mapping: Webb-Marsac   # one of the values in the table below
  ct_use_petra: false        # true only for PETRA-sequence-derived pCT
  water_amp_calibration: true # scale output to the analytic Rayleigh free-water peak
```

Use this path only to apply BabelBrain's own HU conversion to raw UTE or PETRA pseudo-CT HU values. PRESTUS writes a raw-HU `CT.nii.gz` on the same grid as `BabelViscoInput.nii.gz`, and the bridge converts it to the format BabelBrain expects.

!!! warning
    Do not set `ct_mapping` when PRESTUS's pCT pipeline (`pct.enabled: 1`) has already converted HU to acoustic properties. BabelBrain then receives only a label mask and uses its own bone values, which differ from the pCT-derived values k-Wave uses. PRESTUS's per-voxel acoustic properties cannot be passed to BabelBrain (see [Limitations](#limitations-and-unsupported-prestus-modes)).

| `ct_mapping` value | Density | SoS | Attenuation | PRESTUS equivalent |
|---|---|---|---|---|
| `Webb-Marsac` (default) | Marsac linear HU→ρ | Webb 2018 tabulated (freq-dep.) | Webb 2020 tabulated (freq-dep.) | partial: `marsac` density + no Webb |
| `Aubry` | Porosity→ρ (Aubry 2003) | Porosity→SoS | Porosity→Att | `aubry` |
| `Pichardo` | HU→ρ air–tissue linear | Density→SoS (Pichardo 2D map, freq-dep.) | Density→Att (Pichardo map) | none |
| `McDannold` | HU→ρ air–tissue linear | Density→SoS (McDannold polynomial) | Density→Att (McDannold polynomial, freq-dep.) | none |
| `Marsac-Aubry` | Marsac linear HU→ρ | Porosity→SoS (Aubry) | Porosity→Att (Aubry) | partial |
| `Pichardo-Marsac` | Marsac linear HU→ρ | Pichardo 2D map (freq-dep.) | Pichardo map | none |
| `McDannold-Marsac` | Marsac linear HU→ρ | McDannold polynomial | McDannold polynomial (freq-dep.) | none |

PRESTUS `medium_pct_*` functions support `k-plan`, `k-wave`, `marsac` and `aubry` only. For results comparable across backends, use `Webb-Marsac` (closest to the PRESTUS `marsac`/`mueller` combination) or `Aubry` (available in both).

- **PETRA density.** `ct_use_petra: true` makes `Webb-Marsac` use Björn's density formula (`SimNIBS_PETRApct_Density`). Use it for pseudo-CT from a PETRA sequence via SimNIBS's ZTE/PETRA pipeline. For UTE-derived pCT (the PRESTUS default), leave it `false`.
- **Raw UTE as CT input.** A raw UTE volume can be passed as `pseudoCT` with `ct_mapping` set. UTE intensities lack standard 120 kVp HU calibration, so `Webb-Marsac` and `McDannold` give inaccurate results. `Aubry` is the most defensible choice for raw UTE. The recommended route is PRESTUS's pCT pipeline (UTE → calibrated pseudo-CT HU → PRESTUS `marsac`/`aubry` mapping).

#### Internal air (`air_regions`)

On the pseudo-CT path, BabelBrain can treat internal air (frontal and sphenoid sinuses, mastoid cells) as a pressure-release reflector (reflection close to −1). Its solver does not update the fields inside the air mask. BabelBrain's own GUI enables this by default for CT/ZTE/PETRA input.

```yaml
babelbrain:
  ct_mapping: Webb-Marsac
  air_regions: true
  air_hu_window: [-1200, -400]   # default, as BabelBrain
  air_keep_boundary: false
```

PRESTUS writes `<target>_<f>kHz_<ppw>PPW_AirRegions.nii.gz` next to the CT input, on the same grid, and the log reports the voxel count and volume of internal air kept. Check `AirRegions.nii.gz` against the anatomy: in PRESTUS's UTE pseudo-CT, air cells inside the skull mask are captured only where their mapped HU falls in the window.

??? note "Details: how the air mask is built, and its limits"
    1. Take the pseudo-CT voxels within `air_hu_window`.
    2. Apply a 3-voxel median filter and find connected components (26-connectivity).
    3. Remove the largest component (the air around the head).
    4. With `air_keep_boundary: false`, also remove every component that touches the domain boundary. This removes exterior air split by the cropped grid, including the transducer–scalp gap (−1000 HU in the pseudo-CT but water in the simulation).

    Limitations:

    - **Pseudo-CT only.** On the label path there is no air label, so the option is ignored with a warning.
    - **Free-water runs never get the mask.** The water baseline and `water_amp_calibration` are unaffected.
    - **Boundary air is dropped.** Set `air_keep_boundary: true` if the domain edge cuts a sinus and no exterior air touches the boundary.
    - **No k-Wave equivalent.** k-Wave and BabelBrain differ wherever the beam meets internal air.

## Transducer mapping

BabelBrain only accepts named presets. Set `babelbrain.tx_system` to choose one. If it is empty, PRESTUS resolves the preset by:

1. **Serial number**: `tran.serial` from the equipment YAML.
2. **Geometry**: element count, radius of curvature, maximum outer diameter and frequency (tolerances: ±2 mm radius, ±3 mm OD, ±5% frequency).
3. **Custom CSV**: if nothing matches, it writes an element CSV in BabelBrain annular-array format and prints instructions for adding a new preset.

| PRESTUS serial / geometry | BabelBrain preset | Rings | R_curv [mm] | OD [mm] | Freq |
|---|---|---|---|---|---|
| IS_PCD15473_01001 | `R15473` | 10 | 100.0 | 65.0 | 300 kHz |
| IS_PCD15287_01001 | `R15287` | 4 | 62.9 | 64.0 | 500 kHz |
| IS_PCD15646_01001 | `R15646` | 10 | 72.0 | 72.0 | 500 kHz |
| CTX500 geometry | `CTX_500` | 7 | 62.9 | 64.0 | 500 kHz |
| CTX250 geometry | `CTX_250` | 7 | 62.9 | 64.0 | 250 kHz |
| Single-element bowl (any `elem_n: 1` annular) | `Single` | 1 | from `curv_radius_mm` | from `elem_od_mm` | any |
| IGT 64-channel phased array | `IGT64_500` | 64 | – | – | 500 kHz |

The `Single` and `BSonix` presets (`BSonix` only when set as `tx_system`) need `transducer.annular.curv_radius_mm` and `transducer.annular.elem_od_mm`.

## Calibration weights

If PRESTUS has explicit element phases for the transducer, it writes them to the input package as an `OptimizedWeightsFile`, and BabelBrain reproduces the config phases. Without calibrated phases, BabelBrain uses geometric (Rayleigh) focusing, equivalent to the PRESTUS uncalibrated mode. The file format is described in [Calibration weights file](doc_babelbrain_internals.md#calibration-weights-file).

### Free-water amplitude calibration

With `babelbrain.water_amp_calibration: true` (default), PRESTUS scales the output so BabelBrain's free-water peak matches the analytic Rayleigh peak for the same drive. This is the numerical analogue of hydrophone calibration. It requires explicit element phases and is skipped otherwise.

??? note "Details: what the calibration factor does"
    BabelBrain's saved pressure carries a residual bias of a few % from its built-in CFL/dispersion correction (observed: ≈+6% pressure / +13% intensity for R15473 at 300 kHz).

    - `k = p_Rayleigh_peak / p_BabelBrain_water_peak` scales `sensor_data.p_max_all`, and `k²` scales `freefield_isppa_wcm2`.
    - `k` is logged and stored in `sensor_data.amp_calibration`.
    - A factor more than 15% from 1 triggers a warning, because it points to a setup mismatch (geometry, phases, source amplitude) rather than solver bias.

    Why the water factor transfers to tissue runs is explained in [internals](doc_babelbrain_internals.md#free-water-calibration-transfer).

## Thermal simulation

With `thermal_backend: babelbrain`, PRESTUS calls BabelBrain's `CalculateThermalProcess` on the `*DataForSim.h5` acoustic output. BabelBrain solves the Pennes bio-heat equation on the same domain with its GPU solver. What you get:

- `results_heating` has the same structure as from the k-Wave thermal path, so the thermal analysis, NIfTI export and HTML report run unchanged.
- The cooldown starts at `ptrd`, as in the k-Wave thermal path. End-of-sonication maps come from `TempEndFUS`/`DoseEndFUS`; final maps from `FinalTemp`/`FinalDose` (after the cooldown).
- BabelBrain rescales the field so the peak ISPPA in brain equals the in-situ brain peak of the calibrated, ISPPA-scaled PRESTUS acoustic field, so both thermal backends heat from the same acoustic amplitude. The log shows it next to the free-water value and the implied transmission. Set `babelbrain.isppa_w_per_cm2` only to force a specific in-situ brain ISPPA.
- BabelBrain does not compute the ISO CEM43 fields (returned as NaN) or a focal-plane time series.
- BabelBrain's thermal step is fixed at 10 ms, so run time scales with the total simulated duration, including the cooldown.
- With `babelbrain.do_refocusing: true`, the thermal step uses the skull-corrected `/p_amp_refocus` field.

??? note "Details: timing parameter conversion"
    | PRESTUS parameter | BabelBrain parameter | Conversion | Notes |
    |---|---|---|---|
    | `timing.pd` / `timing.pri` | `DutyCycle` | `pd / pri` | |
    | `timing.pri` | `PRF` | `1 / pri` [Hz] | |
    | `timing.ptd` | `DurationUS` | direct [s] | ON duration of one pulse train |
    | `timing.ptri - timing.ptd` | `DurationOff` | direct [s] | Inter-train cooling interval; 0 if `ptri ≤ ptd` |
    | `timing.n_ptri_reps` | `Repetitions` | rounded | Number of (on + off) cycles |
    | `timing.post_ptri_dur` | `PauseBetweenGroupedSonications` | direct [s]; 0 if unset/NaN | Cooldown after the last repetition (`NumberGroupedSonications = 1`) |
    | `thermal.temp_0.brain` | `BaselineTemperature` | direct [°C] | |
    | in-brain peak of `sensor_data.p_max_all` | `Isppa` | `max p²/(2ρc)` over brain [W/cm²] | Override with `babelbrain.isppa_w_per_cm2` |

    PRESTUS reads the result with `read_babelbrain_thermal` and resamples it onto the PRESTUS thermal grid.

### Thermal tissue properties

!!! warning "Known limitation: BabelBrain thermal uses its own tissue thermal properties"
    `thermal_backend: babelbrain` ignores the thermal fields of PRESTUS `medium_properties` (`thermal_conductivity`, `specific_heat_capacity`, `perfusion`, `absorption_fraction`, `absorption_coeff`), including uncertainty-mode variants. BabelBrain has no input for them. The acoustic field amplitude is matched between the thermal backends; the tissue thermal properties are not. Use `thermal_backend: kwave` when PRESTUS thermal properties must apply.

On the label-mask path BabelBrain uses this table (`CalculateTemperatureEffects.py`):

| Tissue | Conductivity [W/m/°C] | Specific heat [J/kg/°C] | Perfusion | Absorption fraction | PRESTUS absorption fraction |
|---|---|---|---|---|---|
| Water | 0.60 | 4178 | 0 | 0 | 1 |
| Skin | 0.37 | 3391 | 106 | 0.85 | 1 |
| Cortical bone | 0.32 | 1313 | 10 | 0.16 | 0.28 |
| Trabecular bone | 0.31 | 2274 | 30 | 0.15 | 0.28 |
| Brain | 0.51 | 3630 | 559 | 0.85 | 1 |

Conductivity, specific heat and perfusion match the PRESTUS defaults, except trabecular conductivity (0.31 vs 0.32). The absorption fraction differs: BabelBrain deposits about 15% less heat in skin and brain and roughly 45% less in skull than k-Wave for the same field. Account for this when comparing temperature rise or CEM43 across thermal backends. On the CT path (`ct_mapping`), all bone voxels get the mean of the cortical and trabecular values (specific heat 1793.5, conductivity 0.315, perfusion 20, absorption fraction 0.155).

## Limitations and unsupported PRESTUS modes

| PRESTUS feature | BabelBrain support | Details |
|---|---|---|
| **Single transducer** | ✓ Full | With `grid.mode: focal_axis` (default) or `transducer_axis`. See [Geometry limits](#geometry-limits). |
| **Multi-transducer, async** (`simulation.transducer_coupling: async`) | ✗ Not supported | Errors. Use `acoustic_backend: kwave` or `jwave` with `thermal_backend: kwave`. |
| **Multi-transducer, coherent** (`simulation.transducer_coupling: coherent`) | ✗ Not supported | Errors. See [Geometry limits](#geometry-limits). |
| **`grid.mode: ras_plus`** | ✗ Not supported | The beam is generally off the grid z axis; `babelbrain_wrapper` errors. |
| **Limited acoustic FOV** (`grid.acoustic_fov_*`) | ✓ Supported | BabelBrain gets the same FOV-cropped box as k-Wave/jwave, voxel for voxel. Its narrow-beam cropping is disabled. |
| **Axisymmetric 2D** (`simulation.axisymmetric: 1`) | ✗ Not supported | BabelBrain is 3D only. Use k-Wave. |
| **`multi_isppa` mode** | ⚠ Unit-tested; end-to-end pending | All derived jobs keep `acoustic_backend`/`thermal_backend`, on `matlab` and HPC. With `thermal_backend: babelbrain`, the per-target thermal jobs run one after another on HPC, because they share one BabelBrain thermal output file. |
| **`phantom` / `water` test medium** | ✓ Supported | Free water: BabelBrain ISPPA +2.6 % vs k-Wave (kWaveArray), FLHM centre −0.5 mm, −3 dB volume −2.0 %. Bone-slab phantom: realised/specified bone attenuation 0.91 for BabelBrain (Q = 1) and k-Wave; transmission vs free field 0.128 vs 0.144. ITRUSST phantoms (`tx_system: Single`): no BabelBrain results yet. |
| **Post-hoc water simulation** | ⚠ Unit-tested; end-to-end pending | Runs BabelBrain water-only on the run's grid. With `babelbrain.cleanup_cache: false` and `io.overwrite_files: never`, it reuses the run's cached `*_Water_DataForSim.h5`. |
| **Thermal tissue properties** | ⚠ Known limitation | BabelBrain heating is lower than k-Wave's for the same field, mostly in skull. See [Thermal tissue properties](#thermal-tissue-properties). |
| **BabelBrain thermal without BabelBrain acoustic** | ✗ Not supported | `thermal_backend: babelbrain` needs the `*DataForSim.h5` of `acoustic_backend: babelbrain`. |
| **CT-based skull properties** | ✓ Supported (via `ct_mapping`) | Raw HU from UTE/PETRA only. Pass the pseudoCT HU volume as the 8th argument to `babelbrain_wrapper`. See [CT-based skull properties](#ct-based-skull-properties-ct_mapping). |
| **PRESTUS pCT acoustic properties** | ✗ Not supported | With `pct.enabled: 1`, BabelBrain uses the label path; k-Wave's per-voxel properties cannot be passed. |
| **Water (coupling medium) properties** | ⚠ Known limitation, by design | PRESTUS `medium_properties.water` is not applied (see the warning above). |
| **Uncertainty mode** | ⚠ Unit-tested; end-to-end pending | Each variant's density, sound speed, `alpha_coeff` and `alpha_power` reach BabelBrain. With `ct_mapping`, skull variants have no effect, while skin and brain variants apply. Variants never change BabelBrain **thermal** properties. No BabelBrain uncertainty run exists yet. |

### Geometry limits

BabelBrain simulates **one transducer per run, with its beam along the z axis of its domain**. PRESTUS passes its own grid without rotating it, so a beam off z would be simulated in the wrong direction while aiming at the right focus. PRESTUS refuses:

- **More than one transducer** with `acoustic_backend` or `thermal_backend` `babelbrain`, coherent or async (error `PRESTUS:babelbrain:multiTransducer`). Checked at config load, in the async pipeline and in `babelbrain_wrapper`. The Nexus flags it before launch.
- **A beam more than 2° off the grid z axis** (error `PRESTUS:babelbrain:obliqueBeam`), checked in `babelbrain_wrapper` from `trans_pos` and `focus_pos`. With one transducer, `focal_axis` and `transducer_axis` put the beam on z; `ras_plus` generally does not.

For multi-transducer runs use `acoustic_backend: kwave` or `jwave` with `thermal_backend: kwave`.

??? note "Details: version tracking"
    PRESTUS records the BabelBrain version with the results:

    ```matlab
    sensor_data.babelbrain_version.version_txt   % e.g. '0.8.1'
    sensor_data.babelbrain_version.git_describe  % e.g. '0.8.0-19-g5237cd2'
    sensor_data.babelbrain_version.git_commit    % full SHA
    sensor_data.babelbrain_version.bb_root       % path used at runtime
    ```

    It is stored in the cached `.mat` file as `acoustic_provenance.babelbrain_version`, so it survives pipeline restarts.

For the call chain, placement, output alignment and debugging, see [BabelBrain integration internals](doc_babelbrain_internals.md).

??? note "References"
    Pichardo, S. BabelBrain: An Open-Source Application for Prospective Modeling of Transcranial Focused Ultrasound for Neuromodulation Applications. *IEEE Trans Ultrason Ferroelectr Freq Control* **70**, 587–599 (2023). https://doi.org/10.1109/TUFFC.2023.3274046
