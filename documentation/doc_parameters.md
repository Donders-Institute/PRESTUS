## PRESTUS parameters

Parameters are organised in nested structs that map directly to YAML keys. PRESTUS uses a two-layer configuration system: `default_config.yaml` defines all parameters and their default values and should not be edited. A study-specific config (e.g. `config_study.yaml`) is loaded on top and only needs to specify values that differ from the defaults. Parameters can also be set or overridden programmatically in MATLAB before calling `prestus_pipeline_start(parameters)`. Mandatory parameters have no meaningful default and must always be provided.

---

<details>
<summary><strong>General</strong></summary>

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `subject_id` | Subject identifier used for file naming and folder management. | — | Mandatory. Set before calling `prestus_pipeline_start(parameters)`. |
| `platform` | Execution platform. | `'auto'` | `auto` / `slurm` / `qsub` / `matlab` |

</details>

---

<details>
<summary><strong><code>simulation</code> — Simulation type & execution</strong></summary>

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `medium` | Simulation medium type. | `'layered'` | `water` / `layered` / `phantom`. Mandatory. |
| `code_type` | k-Wave backend. | `'matlab_gpu'` | `matlab_cpu` / `matlab_gpu` / `cpp_cpu` / `cpp_gpu`. See [doc_backend.md](doc_backend.md). |
| `precision` | Computational precision for acoustic and thermal simulations. | `'single'` | `single` / `double` |
| `interactive` | Interactive mode with prompts and evolving plots. | `0` | `1 = yes`, `0 = no`. If `0`, see `overwrite_files`. |
| `debug` | Verbose debug mode with additional intermediate outputs. | `1` | `1 = yes`, `0 = no` |

</details>

---

<details>
<summary><strong><code>path</code> — Data paths</strong></summary>

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `anat` | Absolute path to structural input data. | — | Mandatory |
| `sim` | Absolute path to simulation outputs. | — | Mandatory |
| `seg` | Absolute path to SimNIBS segmentations (m2m folders). | — | Mandatory |
| `localite` | Path to Localite neuronavigation output folder. | — | Optional |
| `t1_pattern` | T1 image path template relative to `path.anat`. | `'sub-%1$03d_T1w.nii*'` | Supports `%03d`-style subject ID substitution |
| `t2_pattern` | T2 image path template relative to `path.anat`. | `'sub-%1$03d_T2w.nii*'` | Optional |
| `subject_subfolder` | Store outputs in subject-specific subdirectories? | `1` | `1 = yes`, `0 = no` |

</details>

---

<details>
<summary><strong><code>startup</code> — Environment & toolbox paths</strong></summary>

| **Parameter** | **Description** | **Comments** |
|---|---|---|
| `simnibs_bin_path` | Absolute path to SimNIBS binaries. | Mandatory for segmentation and MNI conversion |
| `paths_to_add` | Paths to add with `addpath()`. | [cell] e.g. `{"path/to/x"}` |
| `subpaths_to_add` | Paths to add recursively with `addpath(genpath())`. | [cell]; relative to config file location |

</details>

---

<details>
<summary><strong><code>io</code> — I/O management</strong></summary>

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `output_affix` | Optional affix for output file names. | `''` | Differentiates outputs for the same subject/transducer (e.g. different intensities or targets) |
| `overwrite_files` | File overwrite behaviour. | `'always'` | `never` / `always` / `ask`. Does NOT apply to SimNIBS segmentations. |
| `overwrite_simnibs` | Overwrite SimNIBS segmentation results? | `0` | `1 = yes`, `0 = no` |
| `save_matrices` | Save acoustic/thermal simulation outputs as `.mat`? | `0` | `1 = yes`, `0 = no`. Set to `0` for large batch runs to save disk space. |
| `save_heatingvideo` | Save a video of incremental heating? | `0` | `1 = yes`, `0 = no` |
| `adopted_heatmap` | Path to an existing ISPPA heatmap NIfTI to reuse instead of re-running acoustics. | — | Optional. Used for sequential multi-target runs. |
| `adopted_cem43` | Path to an existing CEM43 heatmap NIfTI to accumulate heating across runs. | — | Optional. Used for sequential multi-target runs. |

**Runtime-derived (set by `path_log_setup`):**

| **Field** | **Description** |
|---|---|
| `output_dir` | Resolved per-subject output directory. Set from `path.sim` (+ subject subfolder if `path.subject_subfolder = 1`). |
| `debug_dir` | Debug subfolder within `output_dir`. Created automatically. |
| `filename_output_table` | Path to the per-subject CSV results table. |
| `kwave_source_filename` | Path to cached k-Wave source `.mat` file. Set by `source_sensor_setup`. |

</details>

---

<details>
<summary><strong><code>modules</code> — Pipeline module flags</strong></summary>

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `run_grid_setup` | Setup grid and run head processing. | `1` | Mandatory for simulations. |
| `run_medium_setup` | Map medium acoustic properties. | `1` | Mandatory for simulations. |
| `run_source_setup` | Set up acoustic source. | `1` | Mandatory for simulations. |
| `run_acoustic_sims` | Run acoustic simulations. | `1` | |
| `run_acoustic_analysis` | Run acoustic analysis. | `1` | |
| `run_heating_sims` | Run thermal simulations. | `0` | Enable after acoustic results are validated. |
| `run_thermal_analysis` | Run thermal analysis. | `1` | |
| `run_nifti_creation` | Export results as NIfTI files. | `1` | |
| `run_posthoc_water_sims` | Run free-water reference simulations after head simulations. | `1` | |
| `generate_report` | Generate self-contained HTML simulation report. | `1` | |
| `segmentation_only` | Stop after segmentation; skip grid setup and all simulations. | `0` | Only has effect when `simulation.medium = 'layered'`. |

All flags: `1 = yes`, `0 = no`.

</details>

---

<details>
<summary><strong><code>transducer</code> — Transducer specification</strong></summary>

see [doc_transducer.md](doc_transducer.md)

All fields are mandatory and have no defaults — they must be set in the study config.

| **Parameter** | **Description** | **Comments** |
|---|---|---|
| `source_freq_hz` | Central frequency of the acoustic source [Hz]. | |
| `n_elements` | Number of transducer elements. | |
| `Elements_ID_mm` | Inner diameter of each element [mm]. | |
| `Elements_OD_mm` | Outer diameter of each element [mm]. | |
| `curv_radius_mm` | Radius of curvature of the transducer bowl [mm]. | |
| `dist_to_plane_mm` | Distance from geometric focus to transducer plane [mm]. | |
| `source_amp` | Pressure amplitude [Pa]. | Must be calibrated. |
| `source_phase_deg` | Source phase [degrees]. | Must be calibrated. |
| `source_phase_rad` | Source phase [radians]. | Must be calibrated. |
| `trans_pos` | Transducer bowl position (XYZ, T1 grid voxel space). | |
| `focus_pos` | Stimulation target position (XYZ, T1 grid voxel space). | |
| `expected_focal_distance_ep` | Expected distance from transducer exit plane to focus [mm]. | Alternative to specifying `trans_pos`/`focus_pos`. Either `expected_focal_distance_ep`, `expected_focal_distance_bowl`, or both pos fields must be set. |
| `expected_focal_distance_bowl` | Expected distance from transducer bowl to focus [mm]. | |

</details>

---

<details>
<summary><strong><code>placement</code> — Transducer placement</strong></summary>

see [doc_transducer.md](doc_transducer.md)

#### `placement.localite`

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `enabled` | Load transducer position from Localite neuronavigation files? | `0` | `1 = yes`, `0 = no` |
| `reference_distance_mm` | Distance from tracker to transducer exit plane [mm]. | `15` | Corrects for varying tracker-to-exit-plane distances. Only applies when `enabled=1`. |

#### `placement.heuristic`

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `save_localite_t1` | Save Localite-aligned T1 output for header correction? | `false` | See [doc_placement_heuristic.md](doc_placement_heuristic.md). |
| `dist_close` | Distance from target considered sufficiently close [mm]. | `[]` | |
| `ear_radius` | Radius of the ear exclusion zone [mm]. | `35` | |
| `left_ear_center` | Approximate coordinates for left ear [image voxels]. | `[]` | Optional |
| `right_ear_center` | Approximate coordinates for right ear [image voxels]. | `[]` | Optional |
| `criterion_intersection` | Max fraction of exit plane intersecting skin. | `0.05` | 5% |
| `criterion_skin_mean` | Mean distance of EP voxels from skin [quantile]. | `[]` | |
| `criterion_skull_mean` | Mean distance of EP voxels from skull [quantile]. | `[]` | |
| `criterion_skin_var` | Variance of EP voxel distance from skin [quantile]. | `[]` | |
| `criterion_skull_var` | Variance of EP voxel distance from skull [quantile]. | `[]` | |
| `expand_step` | Expansion step for skin intersection criterion. | `0.01` | 1% |

</details>

---

<details>
<summary><strong><code>grid</code> — Simulation grid</strong></summary>

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `resolution_mm` | Grid resolution (must be isotropic) [mm]. | `0.5` | |
| `default_dims` | Requested grid dimensions [voxels per dimension]. | `[144, 144, 400]` | Directly sets the simulation grid for `water` and `phantom` media. For `layered`, the grid is determined by head preprocessing and this value is not used. |
| `axisymmetric` | Run axisymmetric 2D simulation (`kspaceFirstOrderAS`). | `0` | `1 = yes`, `0 = no`. See [doc_simulations-acoustic.md](doc_simulations-acoustic.md). |
| `pml_size` | Perfectly Matched Layer (PML) size [voxels]. | `10` | Absorbs waves at boundaries. Recommended for 3D. See [k-Wave docs](http://www.k-wave.org/documentation/example_na_controlling_the_pml.php). |
| `source_ppw` | Points per wavelength. | `[]` | Calculated internally if not set. |
| `source_cfl` | Courant-Friedrichs-Lewy fraction. | `0.15` | |
| `source_limit_fraction` | Fraction of the stability limit to use for time step. | `0.9` | `0` = do not use stability limit |
| `max_expand` | Maximum grid expansion for prime-number FFT optimisation [voxels]. | `40` | |
| `use_kWaveArray` | Use the kWaveArray class for transducer modelling? | `1` | `1 = yes`, `0 = no` |

</details>

---

<details>
<summary><strong><code>headmodel</code> — Head model & segmentation preprocessing</strong></summary>

see [doc_preproc.md](doc_preproc.md)

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `head_pad_mm` | Symmetric padding applied to the cropped head grid prior to transducer + PML setup [mm]. | `0` | |
| `csf_expansion` | Dilation of the CSF brain mask into surrounding head regions [grid voxels]. | `40` | |
| `smooth_method` | Mask smoothing filter type. | `'gaussian'` | `gaussian` / `box` |
| `smooth_threshold_skull` | Binarisation threshold for skull mask; higher = thinner mask. | `0.5` | |
| `smooth_threshold_other` | Binarisation threshold for other masks; higher = thinner mask. | `0.5` | |
| `smooth_fwhm_mm` | FWHM of smoothing kernel [mm]. | `1` | |
| `smooth_properties` | Apply smoothing to acoustic property maps as well? | `false` | |
| `skull_fill_method` | Method for filling holes in the skull. | `'rubberwrap'` | `rubberwrap` / `imclose` |
| `skull_wrap_radius` | Rubber-wrap radius [grid voxels]. | `10` | Larger = tighter wrap ignoring bigger dents. Recommended: 2–10. |
| `skull_wrap_visualize` | Visualize rubber-wrap result? | `0` | `1 = yes`, `0 = no`. Avoid on HPC. |

</details>

---

<details>
<summary><strong><code>segmentation</code> — Segmentation settings</strong></summary>

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `use_qform` | Force qform reorientation before charm segmentation? | `0` | Set to `1` if charm reports a qform/sform mismatch error. |
| `debug` | Pass `--debug` to charm for verbose segmentation output? | `0` | `1 = yes`, `0 = no` |

</details>

---

<details>
<summary><strong><code>pct</code> — pseudo-CT skull property mapping</strong></summary>

see [doc_pseudoCT.md](doc_pseudoCT.md)

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `enabled` | Use (pseudo-)CT to inform skull medium properties? | `0` | `1 = yes`, `0 = no` |
| `mapping_density` | HU-to-density mapping algorithm. | `'k-plan'` | `k-plan` / `k-wave` / `marsac` / `aubry` / `none`. See [doc](doc_pseudoCT.md#mapping-skull-density). |
| `mapping_soundspeed` | HU-to-sound-speed mapping algorithm. | `'k-plan'` | `k-plan` / `marsac` / `aubry` / `none`. See [doc](doc_pseudoCT.md#mapping-skull-sound-speed). |
| `mapping_attenuation` | HU-to-attenuation mapping algorithm. | `'k-plan'` | `k-plan` / `mueller` / `aubry` / `none`. See [doc](doc_pseudoCT.md#mapping-skull-attenuation). |

</details>

---

<details>
<summary><strong><code>layers</code> — Simulation tissue compartments</strong></summary>

see [doc_preproc.md](doc_preproc.md)

Maps tissue compartment names to their SimNIBS charm label indices. Compartments are modelled in listed order; any voxel label not assigned to a named compartment is treated as water. Charm label values are hardcoded in `charm_seg_labels()`. Layers can be removed or added.

| **Layer** | **Default charm labels** | **Comments** |
|---|---|---|
| `water` | `[0, 3, 6, 9, 10]` | Baseline layer; all unassigned voxels are treated as water |
| `brain` | `[1, 2]` | White matter + grey matter |
| `skin` | `[5]` | |
| `skull` | `[4]` | Single skull layer; used when cortical/trabecular split is not needed |
| `skull_cortical` | `[7]` | Cortical bone; used in multi-layer skull model |
| `skull_trabecular` | `[8]` | Trabecular bone; used in multi-layer skull model |

</details>

---

<details>
<summary><strong><code>medium_properties</code> — Tissue acoustic & thermal properties</strong></summary>

see [doc_medium.md](doc_medium.md)

Each tissue compartment (`water`, `brain`, `skin`, `skull`, `skull_trabecular`, `skull_cortical`) carries the following fields:

| **Field** | **Description** | **Units / Reference** |
|---|---|---|
| `sound_speed` | Speed of sound. | m/s — ITRUSST benchmarks |
| `density` | Density. | kg/m³ — ITRUSST benchmarks |
| `alpha_coeff` | Attenuation coefficient. | dB/cm/MHz |
| `alpha_power` | Attenuation power law exponent. | — |
| `thermal_conductivity` | Thermal conductivity. | W/m/°C — Tissue Properties DB |
| `specific_heat_capacity` | Specific heat capacity. | J/kg/°C — Tissue Properties DB |
| `perfusion` | Perfusion / heat transfer rate. | mL/min/kg — Tissue Properties DB |
| `absorption_fraction` | Fraction of attenuation converted to heat. | [0–1] — Pinton et al., 2012 |

</details>

---

<details>
<summary><strong><code>timing</code> — Sonication timing protocol</strong></summary>

see [doc_simulations-thermal.md](doc_simulations-thermal.md)

Protocol duration fields must be set for thermal simulations.

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `pd` | Pulse Duration (PD) [s]. | `NaN` | Duty cycle = `pd`/`pri`. |
| `pri` | Pulse Repetition Interval (PRI) [s]. | `NaN` | PRF = 1/`pri`. |
| `ptd` | Pulse Train Duration (PTD) [s]. | `NaN` | |
| `pt_timestep` | Modelling time step within a pulse train [s]. | `0.02` | |
| `ptri` | Pulse Train Repetition Interval (PTRI) [s]. | `NaN` | OFF duration = `ptri` − `ptd`. |
| `ptrd` | Pulse Train Repetition Duration (PTRD) [s]. | `NaN` | |
| `post_ptri_dur` | Post-PTRI steady-state duration [s]. | `NaN` | |
| `post_pt_timestep` | Modelling time step following PT & PTRI [s]. | `1` | |
| `equal_step_duration` | Equal step durations for on and off cycles? | `0` | `1 = yes`, `0 = no` |

</details>

---

<details>
<summary><strong><code>thermal</code> — Thermal simulation settings</strong></summary>

see [doc_simulations-thermal.md](doc_simulations-thermal.md)

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `cem43_iso` | Calculate CEM43 per ISO norm (`1`) or kWaveDiffusion (`0`)? | `0` | |
| `temp_0.{tissue}` | Initial temperature per tissue compartment [°C]. | `37` | Tissues: `water`, `skull`, `brain`, `skin`, `skull_trabecular`, `skull_cortical` |
| `sensor_xy_halfsize` | Sensor window half-size for temperature recording [grid units]. | `100` | |
| `record_t_at_every_step` | Record temperature at every time step for the full sensor window? | `0` | Memory-intensive; disable if out-of-memory errors occur. |

</details>

---

<details>
<summary><strong><code>analysis</code> — Output analysis</strong></summary>

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `focus_area_radius` | Radius around the focus in which ISPPA is averaged for outputs [mm]. | `5` | |

</details>

---

<details>
<summary><strong><code>hpc</code> — High-performance computing</strong></summary>

see [doc_backend.md](doc_backend.md) [doc_hpc.md](doc_hpc.md)

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `gpu` | Request a specific GPU. | `''` | e.g. `"nvidia_a100-sxm4-40gb:1"` |
| `partition` | Request a dedicated queue partition. | `''` | |
| `reservation` | Request a reserved queue. | `''` | |
| `wait_for_job` | Block MATLAB until the submitted HPC job completes? | `false` | |
| `timelimit` | Job time limit. | `'04:00:00'` | |
| `memorylimit` | Memory limit [GB]. | `20` | |
| `ld_library_path` | LD_LIBRARY path for SimNIBS installation. | `''` | Set if you see `undefined symbol` errors. e.g. `/opt/gcc/7.2.0/lib64` |
| `job_prefix` | Prefix string for HPC job names. | `'PRESTUS'` | Overridden to `'TP'` for transducer positioning jobs. |
| `max_wait_checks` | Maximum number of job status checks when `wait_for_job = true`. | `540` | At 1 check/20 s, `540` ≈ 3 hours. |

</details>

---

<details>
<summary><strong><code>calibration</code> — Transducer calibration</strong></summary>

see [doc_calibration.md](doc_calibration.md)

A separate `calibration_config.yaml` applies for calibration workflows and is loaded as `parameters.calibration`.

| **Parameter** | **Description** | **Comments** |
|---|---|---|
| `path_input_axial` | Directory containing axial intensity profiles. | |
| `path_input_phase` | Directory containing phase data. | |
| `path_output` | Directory for saving free-water simulation results. | |
| `path_output_profiles` | Directory for saving optimised profile data. | |
| `filename_calibrated_CSV` | Filename of calibrated CSV data. | Mandatory only when not generated within standalone script. |
| `save_in_calibration_folder` | Save in `path_output` (`TRUE`) or `sim_path` (`FALSE`). | If `TRUE`, results are appended to existing calibration data. |
| `combinations` | Equipment combinations (must refer to entries in `equipment_config.yaml`). | Multiple combinations can be specified. |
| `focal_depths_wrt_exit_plane` | List of focal depths to characterise [mm]. | |
| `desired_intensities` | Desired free-water intensities [W/cm²]. | |
| `add_FDO` | Append Focal Distance Offset (bowl-to-exit-plane distance)? | Set to `1` if zero point in profiles reflects exit plane rather than bowl. |
| `axisymmetric2D` | Use axisymmetric 2D instead of default 3D free-water simulations? | `1 = yes`, `0 = no` |
| `force_kwavearray` | Force use of kWaveArray for free-water simulations? | If `0`, uses the setting in the default/study config. |
| `opt_method` | Optimisation method. | `FEXminimize` (open source) / `GlobalSearch` (MATLAB Global Optimization Toolbox) |
| `opt_limits` | Distance limits for optimisation [mm]. | |
| `opt_weights` | Weighting of original profile during fitting. | `1` = equal; `>1` = Gaussian (narrower with larger values) |
| `opt_seed` | Random seed for optimisation. | |
| `skip_front_peak_mm` | Distance from profile start to ignore [mm]. | Avoids near-field artefacts in peak/FWHM calculations. |

</details>
