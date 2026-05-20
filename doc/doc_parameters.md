## PRESTUS parameters

Parameters are organised in nested structs that map directly to YAML keys. PRESTUS uses a two-layer configuration system: `config_default.yaml` defines all parameters and their default values and should not be edited. A study-specific config (e.g. `config_study.yaml`) is loaded on top and only needs to specify values that differ from the defaults. Parameters can also be set or overridden programmatically in MATLAB before calling `prestus_pipeline_start(parameters)`. Mandatory parameters have no meaningful default and must always be provided.

---

### General

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `subject_id` | Subject identifier used for file naming and folder management. | — | Mandatory. Set before calling `prestus_pipeline_start(parameters)`. |
| `platform` | Execution platform. | `'auto'` | `auto` / `slurm` / `qsub` / `matlab` |

---

### Simulation type (`simulation`)

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `medium` | Simulation medium type. | `'layered'` | `water` / `layered` / `phantom`. Mandatory. |
| `code_type` | k-Wave backend. | `'matlab_gpu'` | `matlab_cpu` / `matlab_gpu` / `cpp_cpu` / `cpp_gpu`. See [doc_backend.md](doc_backend.md). |
| `precision` | Computational precision for acoustic and thermal simulations. | `'single'` | `single` / `double` |
| `interactive` | Interactive mode with prompts and evolving plots. | `0` | `1 = yes`, `0 = no`. If `0`, see `overwrite_files`. |
| `debug` | Verbose debug mode with additional intermediate outputs. | `1` | `1 = yes`, `0 = no` |

---

### Data paths (`path`)

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `anat` | Absolute path to structural input data. | — | Mandatory |
| `sim` | Absolute path to simulation outputs. | — | Mandatory |
| `seg` | Absolute path to SimNIBS segmentations (m2m folders). | — | Mandatory |
| `localite` | Path to Localite neuronavigation output folder. | — | Optional |
| `t1_pattern` | T1 image path template relative to `path.anat`. | `'sub-%1$03d_T1w.nii*'` | Supports `%03d`-style subject ID substitution |
| `t2_pattern` | T2 image path template relative to `path.anat`. | `'sub-%1$03d_T2w.nii*'` | Optional |

---

### Environment & toolbox paths (`startup`)

| **Parameter** | **Description** | **Comments** |
|---|---|---|
| `simnibs_bin_path` | Absolute path to SimNIBS binaries. | Mandatory for segmentation and MNI conversion |
| `matlab_bin_path` | Absolute path to MATLAB binary. | Required for pseudoCT generation on HPC (e.g. `/opt/matlab/R2024a/bin/matlab`). Local execution falls back to `matlabroot`. |
| `fsl_bin_path` | Absolute path to FSL bin directory. | Required for pseudoCT generation (e.g. `/opt/fsl/bin`). Prepended to PATH at runtime. |
| `ants_bin_path` | Absolute path to ANTs bin directory. | Required for pseudoCT generation (e.g. `/opt/ants/bin`). Prepended to PATH at runtime. |
| `paths_to_add` | Paths to add with `addpath()`. | [cell] e.g. `{"path/to/x"}` |
| `subpaths_to_add` | Paths to add recursively with `addpath(genpath())`. | [cell]; relative to config file location |

---

### I/O management (`io`)

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `output_affix` | Optional affix for output file names. | `''` | Differentiates outputs for the same subject/transducer (e.g. different intensities or targets) |
| `overwrite_files` | File overwrite behaviour. | `'always'` | `never` / `always` / `ask`. Does NOT apply to SimNIBS segmentations. |
| `overwrite_simnibs` | Overwrite SimNIBS segmentation results? | `0` | `1 = yes`, `0 = no` |
| `save_matrices` | Save acoustic/thermal simulation outputs as `.mat`? | `0` | `1 = yes`, `0 = no`. Set to `0` for large batch runs to save disk space. |
| `save_source_matrices` | |
| `save_acoustic_matrices` | |
| `save_thermal_matrices` | |
| `save_source_matrices` | |
| `save_property_maps` | Save per-property NIfTIs to nii/properties/ (layered medium only)? | 0 |
| `save_MNI` | Save NIfTI outputs in MNI space in addition to native T1w space (requires SimNIBS m2m folder)? |
| `save_heatingvideo` | Save a video of incremental heating? | `0` | `1 = yes`, `0 = no` |
| `adopted_heatmap` | Path to an existing **temperature** heatmap NIfTI (`heating_end.nii.gz`) to use as the thermal starting point for a follow-up simulation. Set automatically by the sequential dispatcher; can be overridden manually. | — | Optional. Used for sequential multi-target runs. |
| `adopted_cem43` | Path to an existing CEM43 heatmap NIfTI to accumulate heating across runs. Set automatically by the sequential dispatcher. | — | Optional. Used for sequential multi-target runs. |
| `adopted_cem43_iso` | Path to an existing ISO-variant CEM43 heatmap NIfTI (`CEM43_iso_end.nii.gz`) to accumulate across runs. Set automatically alongside `adopted_cem43` by the sequential dispatcher. | — | Optional. Used for sequential multi-target runs. |

**Runtime-derived (set by `path_log_setup`):**

| **Field** | **Description** |
|---|---|
| `dir_output` | Resolved per-subject output directory. Set from `path.sim` (+ subject subfolder). |
| `dir_nii` | NIfTI output subfolder (`nii/`). Created automatically. |
| `dir_nii_T1w` | Alias for `dir_nii`; space encoded in filename. |
| `dir_nii_MNI` | Alias for `dir_nii`; space encoded in filename. |
| `dir_img` | Image/figure output subfolder (`img/`). Used for all acoustic, thermal, and preprocessing figures. Created automatically. |
| `dir_tabular` | Directory for per-subject CSV results tables (same as `dir_output`). |
| `dir_reports` | Directory for simulation reports (same as `dir_output`). |
| `dir_logs` | Log subfolder (`log/`). Created automatically. |
| `dir_cache` | Cache subfolder (`cache/`). Created automatically. |
| `dir_debug` | Debug subfolder (`debug/`). Created when `simulation.debug = 1`. |
| `dir_debug_preproc` | Debug subfolder for preprocessing (`debug/preproc/`). |
| `dir_debug_medium` | Debug subfolder for medium setup (`debug/medium/`). |
| `dir_debug_source` | Debug subfolder for source setup (`debug/source/`). |
| `dir_pct` | Resolved subject pCT directory. Set from `path.seg` or `path.pct`. |
| `filename_table` | Path to the per-subject CSV results table. |
| `filename_kwave_source` | Path to cached k-Wave source `.mat` file. Set by `source_sensor_setup`. |

---

### Pipeline modules (`modules`)

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

---

### Transducer specification (`transducer`)

See [doc_transducer.md](doc_transducer.md).

All fields are mandatory and have no defaults — they must be set in the study config.

| **Parameter** | **Description** | **Comments** |
|---|---|---|
| `type` | Type of transducer array. | `annular` / `matrix`. Mandatory. |
| `freq_hz` | Fundamental frequency [Hz]. | Shared across all elements. |
| `trans_pos` | Transducer position (XYZ, T1 grid). | |
| `focus_pos` | Focus (target) position (XYZ, T1 grid). | |
| `focal_distance_ep` | Expected focal distance from transducer exit plane [mm]. | Alternative to `trans_pos`/`focus_pos`. Either `focal_distance_ep`, `focal_distance_bowl`, or both must be set. |
| `focal_distance_bowl` | Expected focal distance from transducer bowl [mm]. | |
| `focal_distance_offset` | Offset between transducer bowl and exit plane [mm]. | **Derived** from `curv_radius_mm − dist_geom_ep_mm`. Not set by user. |

#### Annular Array Definition (`annular`)

| **Parameter** | **Description** | **Comments** |
|---|---|---|
| `elem_amp` | Pressure amplitude [Pa]. | Must be calibrated. |
| `elem_phase_deg` | Source phase per element [degrees]. | Must be calibrated. |
| `elem_n` | Number of transducer elements. | |
| `elem_id_mm` | Inner diameter of each element [mm]. | |
| `elem_od_mm` | Outer diameter of each element [mm]. | |
| `curv_radius_mm` | Radius of curvature of the transducer bowl [mm]. | |
| `dist_geom_ep_mm` | Distance from geometric focus to transducer plane [mm]. | Calculated automatically from `curv_radius_mm` if not provided. |
| `depth_mm` | Transducer depth [mm]. | Visualization only. |

#### Matrix Array Definition (`matrix`)

| **Parameter** | **Description** | **Comments** |
|---|---|---|
| `elem_amp` | Pressure amplitude [Pa]. | Must be calibrated. |
| `depth_mm` | Transducer depth [mm]. | Visualization only. |
| `steering` | Steering mode of the transducer. | `1D` = axial only, `3D` = volumetric steering. |
| `elem_shape` | Shape of individual elements. | Options: `rect`, `disc`, `bowl`. Element area is defined using the rectangular dimensions and projected onto the selected shape. |
| `elem_height_mm` | Element height [mm]. | Used to define equivalent area. |
| `elem_width_mm` | Element width [mm]. | Used to define equivalent area. |
| `outer_diameter_mm` | Outer diameter of the transducer [mm]. | Defines active aperture boundary. |
| `is_curved` | Whether the array has a curved surface. | |
| `curv_radius_mm` | Radius of curvature (ROC) of the transducer bowl [mm]. | Defines natural focus. Depends on `is_curved`. |
| `dist_geom_ep_mm` | Distance from geometric focus to transducer plane [mm]. | Calculated automatically from `curv_radius_mm` and `outer_diameter_mm` if not provided. Depends on `is_curved`. |
| `is_clover_setup` | Enables Clover (multi-array) configuration. | Replicates array into multiple leaves. |
| `matrix_shape.type` | Method used to define element positions. | Options: `define_here`, `extract_from_file`. |

##### Matrix: Clover (`clover`)

| **Parameter** | **Description** | **Comments** |
|---|---|---|
| `n_leaves` | Number of Clover leaves. | Maximum: 3. |
| `ROC_parent` | Radius of curvature of the combined Clover setup [mm]. | Can differ from individual array ROC. |

##### Matrix: Grid distribution (`define_here`)

Warning: The `define_here` options are experimental and may contain bugs due to limited testing time.

| **Parameter** | **Description** | **Comments** |
|---|---|---|
| `grid_shape.type` | Grid distribution type. | Options: `rect`, `fibonacci`, `fermat`. |

###### Rectangular Grid

| **Parameter** | **Description** | **Comments** |
|---|---|---|
| `elem_n_row` | Number of element rows. | |
| `elem_n_col` | Number of element columns. | |
| `elem_spacing_height_mm` | Spacing between elements in height direction [mm]. | Edge-to-edge spacing. |
| `elem_spacing_width_mm` | Spacing between elements in width direction [mm]. | Edge-to-edge spacing. |
| `sparsity_factor` | Fraction of active elements. | Range: 0.1–1.0. |

###### Fibonacci Grid

| **Parameter** | **Description** | **Comments** |
|---|---|---|
| `elem_n` | Total number of elements. | |
| `kerf_mm` | Minimum spacing between elements [mm]. | |

###### Fermat Grid

| **Parameter** | **Description** | **Comments** |
|---|---|---|
| `elem_n` | Total number of elements. | Uses spiral distribution. |

##### Matrix: File extraction (`extract_from_file`)

| **Parameter** | **Description** | **Comments** |
|---|---|---|
| `file_path` | Path to coordinate file. | Must contain (x, y, z). |
| `start_row` | Row index where data starts. | MATLAB may skip header row. |
| `start_col` | Column index where data starts. | |
| `elem_n` | Number of elements to extract. | Can exceed final active count. |
| `select_random_subset` | Enable random subset selection. | Useful for sparse arrays. |
| `subset.random_seed` | Controls reproducibility of subset. | `true` = new subset each run. |
| `subset.subset_n_elements` | Number of elements in subset. | Must be ≤ total elements. |
| `project_on_new_ROC` | Project elements onto new curvature. | May introduce unrealistic layouts. |
| `ROC_projection.new_ROC_mm` | New radius of curvature [mm]. | Used when projection enabled. |

### Transducer placement (`placement`)

See [doc_placement.md](doc_placement.md).

#### `placement.localite`

See [doc_placement_neuronav.md](doc_placement_neuronav.md).

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `enabled` | Load transducer position from Localite neuronavigation files? | `0` | `1 = yes`, `0 = no` |
| `reference_distance_mm` | Distance from tracker to transducer exit plane [mm]. | `15` | Corrects for varying tracker-to-exit-plane distances. Only applies when `enabled=1`. |

#### `placement.heuristic`

See [doc_placement_heuristic.md](doc_placement_heuristic.md).

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

---

### Simulation grid (`grid`)

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `resolution_mm` | Grid resolution (must be isotropic) [mm]. | `0.5` | |
| `default_dims` | Requested grid dimensions [voxels per dimension]. | `[144, 144, 400]` | Directly sets the simulation grid for `water` and `phantom` media. For `layered`, the grid is determined by head preprocessing and this value is not used. |
| `axisymmetric` | Run axisymmetric 2D simulation (`kspaceFirstOrderAS`). | `0` | `1 = yes`, `0 = no`. See [doc_simulations-acoustic.md](doc_simulations-acoustic.md). |
| `pml_size` | Perfectly Matched Layer (PML) size [voxels]. | `10` | Absorbs waves at boundaries. Recommended for 3D. See [k-Wave docs](http://www.k-wave.org/documentation/example_na_controlling_the_pml.php). |
| `source_ppw` | Points per wavelength (PPW) override. | `[]` | Calculated internally from `resolution_mm` and the medium's maximum sound speed if not set. Setting this overrides the internal calculation and fixes the number of spatial samples per wavelength used to derive the time step. |
| `min_ppw` | Minimum acceptable PPW at the fundamental frequency. | `6` | Checked against the slowest medium sound speed (worst-case spatial sampling). A warning is raised if the actual PPW falls below this value, along with the maximum `resolution_mm` that would satisfy it. |
| `source_cfl` | Courant-Friedrichs-Lewy (CFL) fraction. | `0.15` | Controls temporal resolution independently of spatial resolution. Given a fixed spatial step `dx`, the time step is `dt = CFL × dx / c_max`, where `c_max` is the fastest medium sound speed (typically skull at ~2800 m/s). Smaller CFL = finer temporal sampling and must remain below the stability limit. Set to 0.15 (half k-Wave's default of 0.3) for additional stability margin in heterogeneous skull simulations. Can be relaxed toward 0.3 to reduce runtime, but should be validated against `checkStability` output. |
| `source_limit_fraction` | Fraction of the stability limit to use for time step. | `0.9` | After source setup, `checkStability` is called; if the CFL-derived `dt` exceeds the estimated limit, source setup reruns with `dt = source_limit_fraction × dt_stability_limit`. The final effective `dt` may therefore be smaller than `CFL × dx / c_max`. Set to `0` to disable. |
| `max_expand` | Maximum grid expansion for prime-number FFT optimisation [voxels]. | `40` | |
| `use_kWaveArray` | Use the kWaveArray class for transducer modelling? | `1` | `1 = yes`, `0 = no` |

---

### Segmentation (`segmentation`)

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `use_qform` | Force qform reorientation before charm segmentation? | `0` | Set to `1` if charm reports a qform/sform mismatch error. |
| `debug` | Pass `--debug` to charm for verbose segmentation output? | `0` | `1 = yes`, `0 = no` |

---

### Head model & segmentation preprocessing (`headmodel`)

See [doc_preproc.md](doc_preproc.md).

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `head_pad_mm` | Symmetric padding applied to the cropped head grid prior to transducer + PML setup [mm]. | `0` | |
| `csf_expansion` | Dilation of the CSF brain mask into surrounding head regions [grid voxels]. | `40` | |
| `smooth_method` | Mask smoothing filter type. | `'gaussian'` | `gaussian` / `box` / `off` |
| `smooth_threshold_skull` | Binarisation threshold for skull mask; higher = thinner mask. | `0.5` | |
| `smooth_threshold_other` | Binarisation threshold for other masks; higher = thinner mask. | `0.5` | |
| `smooth_fwhm_mm` | FWHM of smoothing kernel [mm]. | `1` | |
| `smooth_properties` | Apply smoothing to acoustic property maps as well? | `false` | |
| `skull_fill_method` | Method for filling holes in the skull. | `'rubberwrap'` | `rubberwrap` / `imclose` |
| `skull_wrap_radius` | Rubber-wrap radius [grid voxels]. | `10` | Larger = tighter wrap ignoring bigger dents. Recommended: 2–10. |
| `skull_wrap_visualize` | Visualize rubber-wrap result? | `0` | `1 = yes`, `0 = no`. Avoid on HPC. |

---

### Pseudo-CT skull property mapping (`pct`)

See [doc_pseudoCT.md](doc_pseudoCT.md).

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `enabled` | Use (pseudo-)CT to inform skull medium properties? | `0` | `1 = yes`, `0 = no`. When enabled, pseudoCT is generated automatically from the UTE image (specified via `path.t2_pattern`) if not already present in the SimNIBS folder. |
| `skull_mapping` | UTE→HU linear mapping algorithm used during pseudoCT generation. | `'kosciessa'` | `kosciessa` / `miscouridou` / `carpino` / `wiesinger` / `treeby`. See [doc](doc_pseudoCT.md#pseudoct-generation). |
| `debug` | Keep intermediate NIfTI files in the `pseudoCT/` subfolder. | `1` | `1 = keep`, `0 = delete after generation`. |
| `mapping_density` | HU-to-density mapping algorithm. | `'k-plan'` | `k-plan` / `k-wave` / `marsac` / `aubry` / `none`. See [doc](doc_pseudoCT.md#mapping-skull-density). |
| `mapping_soundspeed` | HU-to-sound-speed mapping algorithm. | `'k-plan'` | `k-plan` / `marsac` / `aubry` / `none`. See [doc](doc_pseudoCT.md#mapping-skull-sound-speed). |
| `mapping_attenuation` | HU-to-attenuation mapping algorithm. | `'k-plan'` | `k-plan` / `mueller` / `aubry` / `none`. See [doc](doc_pseudoCT.md#mapping-skull-attenuation). |

---

### Tissue layers (`layers`)

See [doc_preproc.md](doc_preproc.md).

Maps tissue compartment names to their SimNIBS charm label indices. Compartments are modelled in listed order; any voxel label not assigned to a named compartment is treated as water. Charm label values are hardcoded in `charm_seg_labels()`. Layers can be removed or added.

| **Layer** | **Default charm labels** | **Comments** |
|---|---|---|
| `water` | `[0, 3, 6, 9, 10]` | Baseline layer; all unassigned voxels are treated as water |
| `brain` | `[1, 2]` | White matter + grey matter |
| `skin` | `[5]` | |
| `skull` | `[4]` | Single skull layer; used when cortical/trabecular split is not needed |
| `skull_cortical` | `[7]` | Cortical bone; used in multi-layer skull model |
| `skull_trabecular` | `[8]` | Trabecular bone; used in multi-layer skull model |

---

### Tissue acoustic & thermal properties (`medium_properties`)

See [doc_medium.md](doc_medium.md).

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

---

### Sonication timing (`timing`)

See [doc_simulations-thermal.md](doc_simulations-thermal.md).

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

---

### Thermal settings (`thermal`)

See [doc_simulations-thermal.md](doc_simulations-thermal.md).

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `cem43_iso` | Calculate CEM43 per ISO norm (`1`) or kWaveDiffusion (`0`)? | `0` | |
| `temp_0.{tissue}` | Initial temperature per tissue compartment [°C]. | `37` | Tissues: `water`, `skull`, `brain`, `skin`, `skull_trabecular`, `skull_cortical` |
| `sensor_xy_halfsize` | Sensor window half-size for temperature recording [grid units]. | `100` | |
| `record_t_at_every_step` | Record temperature at every time step for the full sensor window? | `0` | Memory-intensive; disable if out-of-memory errors occur. |

---

### Output analysis (`analysis`)

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `focus_area_radius` | Radius around the focus in which IPA is averaged for outputs [mm]. | `5` | |

---

### High-performance computing (`hpc`)

See [doc_backend.md](doc_backend.md) and [doc_hpc.md](doc_hpc.md).

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

---

### Transducer calibration (`calibration`)

See [doc_calibration.md](doc_calibration.md).

A separate `config_calibration.yaml` applies for calibration workflows and is loaded as `parameters.calibration`.

| **Parameter** | **Description** | **Comments** |
|---|---|---|
| `path_input_axial` | Directory containing axial intensity profiles. | |
| `path_input_phase` | Directory containing phase data. | |
| `path_output` | Directory for saving free-water simulation results. | |
| `path_output_profiles` | Directory for saving optimised profile data. | |
| `filename_calibrated_CSV` | Filename of calibrated CSV data. | Mandatory only when not generated within standalone script. |
| `save_in_calibration_folder` | Save in `path_output` (`TRUE`) or `sim_path` (`FALSE`). | If `TRUE`, results are appended to existing calibration data. |
| `combinations` | Equipment combinations (must refer to entries in `config_equipment.yaml`). | Multiple combinations can be specified. |
| `focal_depths_wrt_exit_plane` | List of focal depths to characterise [mm]. | |
| `desired_intensities` | Desired free-water intensities [W/cm²]. | |
| `add_FDO` | Append Focal Distance Offset (bowl-to-exit-plane distance)? | Set to `1` if zero point in profiles reflects exit plane rather than bowl. |
| `axisymmetric2D` | Use axisymmetric 2D instead of default 3D free-water simulations? | `1 = yes`, `0 = no` |
| `force_kwavearray` | Force use of kWaveArray for free-water simulations? | If `0`, uses the setting in the default/study config. |
| `opt_method` | Optimisation method. | `FEXminimize` (open source) / `GlobalSearch` (MATLAB Global Optimization Toolbox) |
| `opt_limits` | Distance limits for optimisation [mm]. | |
| `opt_weights` | Weighting of original profile during fitting. | `1` = equal; `>1` = Gaussian (narrower with larger values) |
| `opt_seed` | Random seed for optimisation. | |
| `opt_use_initial_phases` | Use phases from `transducer.annular.elem_phase_rad` as the initial guess instead of random initialization. | `false` (default); set to `true` when manufacturer-provided or previously calibrated phases are available. |
| `opt_upper_velocity` | Upper velocity in global search | |
| `opt_phase_precession` | Constrain element phases during optimisation. | `false` (default) = independent; `'linear'` = strict ramp (`phase_start + i * phase_step`, 2 free params); `'monotonic'` = ordered phases via cumulative non-negative increments (N params, monotonically constrained). |
| `opt_amp_validation` | Run free-water simulation to validate phases and amplitudes?| 'always' (default) / 'never' / 'initial' (only for the first target intensity) / 'final' (only for the final target intensity) | |
| `skip_front_peak_mm` | Distance from profile start to ignore [mm]. | Avoids near-field artefacts in peak/FWHM calculations. |
