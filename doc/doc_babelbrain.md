# BabelBrain Backend

[BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain) solves the wave equation with a staggered-grid FDTD solver, a different numerical method from [k-Wave](https://github.com/ucl-bug/k-wave)'s pseudo-spectral solver, and PRESTUS gives it the same drive phases and tissue properties. Comparing a BabelBrain run of the same setup with a k-Wave run indicates how far a result depends on the solver. With `thermal_backend: babelbrain` it also computes the heating, with BabelBrain's own tissue thermal properties.

On this page:

- [Quick Start](#quick-start): selecting the backend and its Python environment.
- [Installation](#installation): the submodule and conda environment.
- [Comparison with k-Wave](#comparison-with-k-wave): solver differences and which tissue properties reach BabelBrain.
- [BabelBrain Output Folder](#babelbrain-output-folder): the input package for running BabelBrain on its own.
- [Transducer Mapping](#transducer-mapping): how PRESTUS transducers map to BabelBrain presets.
- [Calibration Weights](#calibration-weights): element phases and the free-water amplitude calibration.
- [Thermal Simulation](#thermal-simulation): BabelBrain heating, its timing, and its tissue thermal properties.
- [Limitations and Unsupported PRESTUS Modes](#limitations-and-unsupported-prestus-modes): which PRESTUS features work with BabelBrain.

!!! warning "Experimental"
    The BabelBrain backend is not yet validated against measurements or the [ITRUSST](https://itrusst.org/) benchmarks in PRESTUS. BabelBrain thermal runs use BabelBrain's own tissue thermal properties ([known limitation](#limitations-and-unsupported-prestus-modes)). Use it to cross-check k-Wave, not to replace it.

!!! note "At a glance"
    - **When you need it:** To cross-check a k-Wave result with an independent solver, or to use the Apple Silicon GPU.
    - **Before you start:** The `external/BabelBrain` submodule and a BabelBrain [conda](https://docs.conda.io/) environment ([Installation](#installation)).

PRESTUS can run [BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain) (Samuel Pichardo, University of Calgary) in place of k-Wave. Pre-processing (segmentation, grid, medium masks, placement) and post-processing (focal metrics, [CEM43](doc_glossary.md#cem43), [NIfTI](https://nifti.nimh.nih.gov/) export, HTML report) run unchanged. BabelBrain's solver, [BabelViscoFDTD](https://github.com/ProteusMRIgHIFU/BabelViscoFDTD), solves the viscoelastic wave equation and adds a Rayleigh integral step for the water coupling layer. PRESTUS runs it in fluid-bone mode by default (shear off, PRESTUS tissue properties injected), so results are comparable to k-Wave and jwave. On Apple Silicon it has a native [Metal](https://developer.apple.com/metal/) backend, where k-Wave has the slower `matlab_gpu` path and the experimental, not yet validated [`cpp_mlx`](doc_kwave.md) binary.

## Quick Start

Set three things: the backend, the [Python](https://www.python.org/) environment and (optionally) the transducer preset.

```yaml
acoustic_backend: babelbrain

babelbrain:
  python_env: /opt/conda/envs/BabelBrain   # directory containing bin/python; must be set
  tx_system: R15473                        # BabelBrain preset name; auto-detected if empty
  computing_backend: auto                  # 'auto' | 'cuda' | 'opencl' | 'metal' | 'cpu'  (or 1-4)
```

`babelbrain.python_env` must be set. Empty falls back to the system `python`, which fails with `No module named 'BabelViscoFDTD'`. PRESTUS finds BabelBrain at `external/BabelBrain` ([git submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules)); set `babelbrain.root` to override. For `tx_system`, see [Transducer mapping](#transducer-mapping).

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

??? note "All `babelbrain` keys"
    ```yaml
    babelbrain:
      tx_system: R15473            # BabelBrain preset name; auto-detected if empty
      computing_backend: auto      # 'auto' | 'cuda' | 'opencl' | 'metal' | 'cpu'  (or 1-4)
      device_name: auto            # 'auto' selects device based on computing_backend; or set explicitly e.g. 'Apple'
      ppw: 12                      # points-per-wavelength floor (default 12)
      do_refocusing: false         # skull back-propagation refocusing
      tight_narrow_beam: ~         # BabelBrain narrow-beam domain (grid.z_beyond_focus_mm past the focus); empty/false = off; ignored when grid.acoustic_fov_length_mm is set
      python_env: ''               # conda/venv environment or its Python executable; empty = system python
      root: ''                     # path to BabelBrain source; auto-detected if empty
      tx_mech_adjust_z_m: 0.0      # ADDITIONAL transducer z offset into skin [m]; a placement
                                   # correction (preset FL − |focus_pos − trans_pos|) is always added
      z_steering_m: ~              # ZSteering [m]; null = auto-derived (0 with explicit phases)
      q_correction: []             # bone Q-factor correction of the relaxation attenuation model; [] = 1.0 non-CT / 3.0 CT
      shear: 'off'                 # bone shear: 'off' (fluid bone) | 'prestus_elastic' (kwave_elastic skull properties)
      drive: prestus               # 'prestus': reproduce the PRESTUS element phases | 'native': BabelBrain steers itself
      focal_calibration: true              # geometric fallback only: steer so the achieved FWHM-centre lands at focal_distance_ep;
                                           # locate-only (report) when explicit element phases exist
      focal_calibration_target: flhm_center # quantity placed at focal_distance_ep: 'flhm_center' | 'peak'
      profile_calibration: false   # annular explicit phases only: refit BabelBrain phases to the O'Neil axial profile
      water_amp_calibration: true  # scale output so the free-water ISPPA equals target_isppa_wcm2
      use_rayleigh_for_water: false # water-only solver: false = full-wave FDTD, true = Rayleigh integral
      ct_mapping: ''               # HU-to-property mapping for CT-based skull; empty = label-based path
      air_regions: false           # pseudo-CT path only: internal air as pressure-release reflector mask
      air_hu_window: [-1200, -400] # HU range classified as air
      air_keep_boundary: false     # keep air components touching the domain boundary
      cleanup_cache: true          # delete the per-subject BabelBrain cache after reading
      dry_run: false               # prepare inputs only; stop after the acoustic stage
      # isppa_w_per_cm2: ~         # thermal only; leave empty: derived from the scaled acoustic field (in-situ brain peak)
    ```

    Placement correction, `z_steering_m`, and `focal_calibration` are explained in [Placement Correction and Drive Phases](doc_babelbrain_internals.md#placement-correction-and-drive-phases).

    **PPW and resolution.** PRESTUS uses `max(ppw, round(c_min / (freq × res_mm × 1e-3)))`, where `c_min` is the smaller of BabelBrain's bone shear-speed floor (≈1102.5 m/s) and the lowest `medium_properties` sound speed. PRESTUS resamples the mask to BabelBrain's native step, so resolutions need not be matched by hand. Lower `ppw` runs faster with more axial numerical dispersion.

**Progress and live view.** PRESTUS reads BabelBrain's acoustic and thermal time steps every 5 s and reports them to the [Nexus](doc_nexus_tabs.md#run-tab) progress bar. With [`simulation.live_view: 1`](doc_parameters_simulation.md#simulation-type-simulation) and the CUDA backend, it also saves the accumulated pressure and the temperature fields as they evolve.

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

Put the environment directory (the one containing `bin/python`), or the Python executable itself, into `babelbrain.python_env`. For GPU backend choice per platform, see [GPU backend selection](doc_babelbrain_internals.md#gpu-backend-selection-and-platform-support).

## Comparison with k-Wave

| | PRESTUS (k-Wave) | BabelBrain |
|---|---|---|
| **Solver** | Pseudo-spectral FDTD | Explicit staggered-grid FDTD |
| **Physics** | Scalar linear acoustic (pressure/density) | Viscoelastic solver; run with shear = 0 (fluid bone) in PRESTUS unless `babelbrain.shear: prestus_elastic` |
| **Skull model** | Scalar ρ, c, α | Scalar ρ, c, α from PRESTUS `medium_properties` (or BabelBrain HU mapping via `ct_mapping`) |
| **Rayleigh pre-step** | None | Rayleigh integral for the water coupling layer |
| **Narrow-beam domain** | Optional (see [partial FOV](doc_partial-fov.md)) | Off by default (the pre-focal field would not be returned); enable with `babelbrain.tight_narrow_beam: true`. Always off when `grid.acoustic_fov_length_mm` is set (mask already cropped) |
| **GPU backends** | [CUDA](https://developer.nvidia.com/cuda-toolkit), MATLAB GPU arrays | CUDA, [OpenCL](https://www.khronos.org/opencl/), Metal (Apple Silicon) |

See [k-Wave](doc_kwave.md) for the k-Wave backend options.

!!! warning "Tissue properties and comparability"
    On the default label-mask path, PRESTUS writes `medium_properties` to `<target>_medium_properties.json`. BabelBrain then uses PRESTUS density, sound speed, and attenuation for skin, cortical bone, trabecular bone, and brain (an undifferentiated `skull` fills both bone slots), with shear = 0 unless `babelbrain.shear: prestus_elastic` is set. Remaining differences are the solver and the bone Q-correction (`q_correction`), so validate against a k-Wave run before relying on absolute values.

    PRESTUS properties do **not** apply to:

    - **water**: BabelBrain keeps 1000 kg/m³, 1500 m/s, no attenuation (PRESTUS default 994 kg/m³, 1500 m/s, 2.17e-3 dB/(cm·MHz²)). For identical water across backends, set `medium_properties.water` to 1000 kg/m³ / 1500 m/s.
    - **skull on the `ct_mapping` path**, where BabelBrain's HU models set bone (skin and brain are still injected).
    - runs where writing the JSON failed, or where the drive frequency has no entry in BabelBrain's table (5 kHz steps). These print a warning and fall back to [BabelBrain's built-in table](doc_babelbrain_internals.md#built-in-material-table).
    - **voxel-wise pseudo-CT properties** computed by PRESTUS (use `babelbrain.medium: prestus`).
    - **smoothed properties** (`headmodel.smooth_properties`; use `babelbrain.medium: prestus`).
    - **thermal** properties (see [Thermal tissue properties](#thermal-tissue-properties)).

### PRESTUS Medium

With `babelbrain.medium: prestus`, BabelViscoFDTD runs on the acoustic medium k-Wave would get: the same density, sound speed, and attenuation in every voxel, including `headmodel.smooth_properties` and a pseudo-CT skull. Use it to compare the solvers on one medium. The label path differs from k-Wave whenever PRESTUS smooths the properties: the step boundaries of a label map reflect more at the skull than a smoothed medium, which changes the transmitted intensity by several percent.

PRESTUS groups voxels with equal properties into up to 1000 materials and logs the largest deviation of a voxel from its material. The medium is fluid, and BabelBrain's Q correction is 1 for every material, so the attenuation is the PRESTUS value. BabelBrain's thermal solver reads its tissue labels, so compute heating with `thermal_backend: kwave`. `ct_mapping`, `shear: prestus_elastic`, and `q_correction` act on BabelBrain's labels and are rejected.

## BabelBrain Output Folder

With the BabelBrain backend active, or with `io.export_babelbrain: true`, PRESTUS writes a `babelbrain/` folder next to the standard outputs (`nii/`, `log/`, `cache/`). It holds everything needed to re-run BabelBrain on its own:

| File | Description |
|---|---|
| `<target>_<freq>kHz_<PPW>PPW_BabelViscoInput.nii.gz` | Tissue mask in BabelBrain integer label format, including the focal voxel marker |
| `<target>_CalibWeights.h5` | Per-element calibration phases from PRESTUS calibration |
| `RunBabelBrain.py` | Ready-to-run Python script with all simulation parameters pre-filled |
| `babelbrain_config.yaml` | Parameter record for reproducibility |

To export the package while simulating with k-Wave, set `io.export_babelbrain: true` with `acoustic_backend: kwave` (or no `acoustic_backend` key). Use this to prepare a cross-validation or to run BabelBrain on another machine.

How PRESTUS tissue labels and skull properties reach BabelBrain is described in [BabelBrain Internals](doc_babelbrain_internals.md#tissue-labels-and-material-properties).

## Transducer Mapping

By default (`babelbrain.tx_system` empty) BabelBrain simulates the PRESTUS transducer itself: the rings of `transducer.annular` (`elem_id_mm`, `elem_od_mm`, `curv_radius_mm`) on BabelBrain's generic focused annular template. This is the geometry the shipped calibrations were fitted on, so their per-ring phases and amplitudes apply ring for ring, and BabelBrain, k-Wave, and jwave simulate the same source. Single-element bowls are one ring. The outer ring diameter sets BabelBrain's aperture.

Set `babelbrain.tx_system` to a BabelBrain preset to use BabelBrain's own device definition instead, for example to compare against it, or for a phased array PRESTUS does not describe as rings (`H317`, `IGT64_500`, ...). Presets are listed in BabelBrain's registry, `BabelBrain/SelFiles/transducer_list.yaml`; `REMOPD` is accepted for `TheraPhas`. The `Single` and `BSonix` presets take their aperture and focal length from `transducer.annular`.

## Calibration Weights

If PRESTUS has explicit element phases for the transducer, it writes them to the input package as an `OptimizedWeightsFile`, and BabelBrain reproduces the config phases. Without calibrated phases, BabelBrain uses geometric (Rayleigh) focusing, equivalent to the PRESTUS uncalibrated mode. The file format is described in [Calibration weights file](doc_babelbrain_internals.md#calibration-weights-file).

### Free-Water Amplitude Calibration

With `babelbrain.water_amp_calibration: true` (default), PRESTUS scales the output so the free-water [I<sub>SPPA</sub>](doc_glossary.md#isppa) of BabelBrain's own water-only run equals `transducer.target_isppa_wcm2`. This is the numerical analogue of hydrophone calibration. It is skipped, and the raw solver field kept, when there is no water run or no target.

??? note "What the calibration factor does"
    BabelBrain's saved pressure carries a residual bias of a few % from its built-in CFL/dispersion correction (observed: ≈+6% pressure / +13% intensity for R15473 at 300 kHz). Normalising to BabelBrain's own free-water run removes it together with the drive scaling.

    - `k = sqrt(target_isppa_wcm2 / free-water ISPPA)` scales `sensor_data.p_max_all` (and `p_complex`), and `freefield_isppa_wcm2` becomes the target.
    - `k` is logged and stored in `sensor_data.amp_calibration`.
    - A factor above 4 or below 0.25 triggers a warning, because it points to a setup mismatch (e.g. a misplaced or defocused free-water focus) rather than drive scaling.

    Why the water factor transfers to tissue runs is explained in [internals](doc_babelbrain_internals.md#free-water-calibration-transfer).

## Thermal Simulation

With `thermal_backend: babelbrain`, PRESTUS calls BabelBrain's `CalculateThermalProcess` on the `*DataForSim.h5` acoustic output. BabelBrain solves the Pennes bio-heat equation on the same domain with its GPU solver. What you get:

- `results_heating` has the same structure as from the k-Wave thermal path, so the thermal analysis, NIfTI export, and HTML report run unchanged.
- The cooldown starts at `ptrd`, as in the k-Wave thermal path. End-of-sonication maps come from `TempEndFUS`/`DoseEndFUS`; final maps from `FinalTemp`/`FinalDose` (after the cooldown).
- BabelBrain rescales the field so the peak I<sub>SPPA</sub> in brain equals the in-situ brain peak of the calibrated, I<sub>SPPA</sub>-scaled PRESTUS acoustic field, so both thermal backends heat from the same acoustic amplitude. The log shows it next to the free-water value and the implied transmission. Set `babelbrain.isppa_w_per_cm2` only to force a specific in-situ brain I<sub>SPPA</sub>.
- BabelBrain does not compute the ISO CEM43 fields (returned as NaN) or a focal-plane time series.
- BabelBrain's thermal step is fixed at 10 ms, so run time scales with the total simulated duration, including the cooldown.
- With `babelbrain.do_refocusing: true`, the thermal step uses the skull-corrected `/p_amp_refocus` field.

??? note "Timing parameter conversion"
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

### Thermal Tissue Properties

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

Conductivity, specific heat, and perfusion match the PRESTUS defaults, except trabecular conductivity (0.31 vs 0.32). The absorption fraction differs: BabelBrain deposits about 15% less heat in skin and brain and roughly 45% less in skull than k-Wave for the same field. Account for this when comparing temperature rise or CEM43 across thermal backends. On the CT path (`ct_mapping`), all bone voxels get the mean of the cortical and trabecular values (specific heat 1793.5, conductivity 0.315, perfusion 20, absorption fraction 0.155). BabelBrain's heat source uses the specified bone attenuation, while the wave propagates with that attenuation divided by the Q-correction (3 on the CT path; [Pichardo 2023](https://doi.org/10.1109/TUFFC.2023.3274046)). See [Bone attenuation and heat source](doc_babelbrain_internals.md#bone-attenuation-and-heat-source).

## Limitations and Unsupported PRESTUS Modes

| PRESTUS feature | BabelBrain support | Details |
|---|---|---|
| **Single transducer** | ✓ Full | With `grid.mode: focal_axis` (default) or `transducer_axis`. See [Geometry limits](#geometry-limits). |
| **Multi-transducer, async** (`simulation.transducer_coupling: async`) | ✗ Not supported | Errors. Use `acoustic_backend: kwave` or `jwave` with `thermal_backend: kwave`. |
| **Multi-transducer, coherent** (`simulation.transducer_coupling: coherent`) | ✗ Not supported | Errors. See [Geometry limits](#geometry-limits). |
| **`grid.mode: ras_plus`** | ✗ Not supported | The beam is generally off the grid z axis; `babelbrain_wrapper` errors. |
| **Partial FOV** (`grid.acoustic_fov_*`) | ✓ Supported | BabelBrain gets the same FOV-cropped box as k-Wave/jwave, voxel for voxel. Its narrow-beam cropping is disabled. |
| **Axisymmetric 2D** (`grid.axisymmetric: 1`) | ✗ Not supported | BabelBrain is 3D only. Use k-Wave. |
| **`multi_isppa` mode** | ⚠ Unit-tested; end-to-end pending | All derived jobs keep `acoustic_backend`/`thermal_backend`, on `matlab` and HPC. With `thermal_backend: babelbrain`, the per-target thermal jobs run one after another on HPC, because they share one BabelBrain thermal output file. |
| **`phantom` / `water` test medium** | ✓ Supported | Free water: BabelBrain I<sub>SPPA</sub> +2.6 % vs k-Wave (kWaveArray), FLHM centre −0.5 mm, −3 dB volume −2.0 %. Bone-slab phantom: realised/specified bone attenuation 0.91 for BabelBrain (Q = 1) and k-Wave; transmission vs free field 0.128 vs 0.144. ITRUSST phantoms (`tx_system: Single`): no BabelBrain results yet. |
| **Post-hoc water simulation** | ⚠ Unit-tested; end-to-end pending | Runs BabelBrain water-only on the run's grid. With `babelbrain.cleanup_cache: false` and `io.overwrite_files: never`, it reuses the run's cached `*_Water_DataForSim.h5`. |
| **Thermal tissue properties** | ⚠ Known limitation | BabelBrain heating is lower than k-Wave's for the same field, mostly in skull. See [Thermal tissue properties](#thermal-tissue-properties). |
| **BabelBrain thermal without BabelBrain acoustic** | ✗ Not supported | `thermal_backend: babelbrain` needs the `*DataForSim.h5` of `acoustic_backend: babelbrain`. |
| **CT-based skull properties** | ✓ Supported (via `ct_mapping`) | Raw HU from UTE/PETRA only. Pass the pseudo-CT HU volume as the 8th argument to `babelbrain_wrapper`. See [CT-based skull properties](doc_babelbrain_internals.md#ct-based-skull-properties-ct_mapping). |
| **PRESTUS pCT acoustic properties** | ✗ Not supported | With `pct.enabled: 1`, BabelBrain uses the label path; k-Wave's per-voxel properties cannot be passed. |
| **Water (coupling medium) properties** | ⚠ Known limitation, by design | PRESTUS `medium_properties.water` is not applied (see the warning above). |
| **Uncertainty mode** | ⚠ Unit-tested; end-to-end pending | Each variant's density, sound speed, `alpha_coeff`, and `alpha_power` reach BabelBrain. With `ct_mapping`, skull variants have no effect, while skin and brain variants apply. Variants never change BabelBrain **thermal** properties. No BabelBrain uncertainty run exists yet. |

### Geometry Limits

BabelBrain simulates **one transducer per run, with its beam along the z axis of its domain**. PRESTUS passes its own grid without rotating it, so a beam off z would be simulated in the wrong direction while aiming at the right focus. PRESTUS refuses:

- **More than one transducer** with `acoustic_backend` or `thermal_backend` `babelbrain`, coherent or async (error `PRESTUS:babelbrain:multiTransducer`). Checked at config load, in the async pipeline, and in `babelbrain_wrapper`. The Nexus flags it before launch.
- **A beam more than 2° off the grid z axis** (error `PRESTUS:babelbrain:obliqueBeam`), checked in `babelbrain_wrapper` from `trans_pos` and `focus_pos`. With one transducer, `focal_axis` and `transducer_axis` put the beam on z; `ras_plus` generally does not.

For multi-transducer runs use `acoustic_backend: kwave` or `jwave` with `thermal_backend: kwave`.

??? note "Version tracking"
    PRESTUS records the BabelBrain version with the results:

    ```matlab
    sensor_data.babelbrain_version.version_txt   % e.g. '0.8.8'
    sensor_data.babelbrain_version.git_describe  % e.g. 'v0.8.8'
    sensor_data.babelbrain_version.git_commit    % full SHA
    sensor_data.babelbrain_version.bb_root       % path used at runtime
    ```

    It is stored in the cached `.mat` file as `acoustic_provenance.babelbrain_version`, so it survives pipeline restarts.

For the call chain, placement, output alignment, and debugging, see [BabelBrain integration internals](doc_babelbrain_internals.md).

??? note "References"
    Pichardo, S. BabelBrain: An Open-Source Application for Prospective Modeling of Transcranial Focused Ultrasound for Neuromodulation Applications. *IEEE Trans Ultrason Ferroelectr Freq Control* **70**, 587–599 (2023). https://doi.org/10.1109/TUFFC.2023.3274046
