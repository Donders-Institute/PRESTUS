## PRESTUS functions

The following documents the functions provided in PRESTUS.

#### ACOUSTIC

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `acoustic_analysis`                     | Compute acoustic safety metrics (Isppa, MI, max pressure) from sensor data, extract tissue-specific maxima, focal distances, and generate overlaid ISPPA plots on segmentation. |
| `acoustic_convert_axisymmetry`          | Expand axisymmetric simulation fields to 3D or the mirrored 2D plane in the form chosen by `axisymmetric_output_mode`, after a simulation and on acoustic cache reload |
| `acoustic_fov_crop`                     | Resolve FOV settings and crop the acoustic grid if active. |
| `acoustic_grid_fov`                     | Crop the acoustic simulation domain to a beam-axis FOV box. |
| `acoustic_input_file`                   | Pressure map a run starts from instead of simulating ('' if none) |
| `acoustic_simulation`                   | Set up acoustic simulation |
| `acoustic_simulation_elastic`           | Run a k-Wave elastic (shear-wave) simulation via pstdElastic2D/3D |
| `acoustic_stage`                        | Run or load the acoustic pressure simulation (Stage 6). |
| `acoustic_water_baseline`               | Decide and run the free-water baseline simulation. |
| `acoustic_wrapper`                      | Set up acoustic simulation |
| `assemble_limited_fov_fields`           | Back-project and combine N limited-FOV acoustic results. |
| `assemble_pfov_fields`                  | Deprecated alias of ASSEMBLE_LIMITED_FOV_FIELDS. |
| `axisymmetric_output_mode`              | Form into which axisymmetric results are expanded |
| `cache_acoustic`                        | Save or load the slim acoustic results cache. |
| `coherent_drive_calibration`            | Set each transducer's drive to its own ISPPA target (coherent runs) |
| `compute_focal_metrics`                 | Compute focal-zone metrics (ITRUSST −3 dB/−6 dB volumes, FLHM span/centre) from a 3D intensity field; FLHM centre threshold is configurable via `focal_reference` (see [doc_calibration.md](doc_calibration.md)). |
| `fov_backproject_volume`                | Zero-pad a FOV sub-grid volume back to the full acoustic grid. |
| `fov_restore_transducers`               | Put every transducer back at its full-grid position after an FOV crop. |
| `fov_to_grid_pos`                       | Convert FOV-grid voxel indices to full sim-grid indices. |
| `is_acoustic_dry_run`                   | True when the configured acoustic backend only writes inputs |
| `load_acoustic_input`                   | Start a run from a prior pressure map instead of simulating |
| `needs_water_baseline`                  | Decide whether a free-water baseline simulation must run |
| `plot_axial_profile`                    | Plot and save pressure and intensity along the beam axis (depth from the exit plane) with tissue shading and the free-water O'Neil / k-Wave baseline intensity of the simulated transducer. |
| `predict_focal_shift`                   | Focus near the target predicted from the reverse field over an aperture |
| `prepare_thermal_domain`                | Restore full-grid state or adjust planimg for thermal FOV. |
| `resolve_acoustic_fov_params`           | Resolve 'auto' FOV settings and check readiness. |
| `reverse_point_source_field`            | Complex CW field of a point source at the target, at the sensor points |
| `superpose_complex_fields`              | Coherent sum of complex pressure fields from separate runs. |
| `targeting_analysis`                    | Targeting accuracy of one run from its focal-mask NIfTIs: simulated focus vs planned focus (axial/lateral), vs MNI target and ROI; writes `sub-XXX_<medium><affix>_targeting.csv`. |
| `time_reversal_pose_metrics`            | Target intensity and focus position of one pose-verification run |
| `time_reversal_pose_parameters`         | Run parameters that place the transducer at a pose of a pose search |
| `time_reversal_refocus`                 | Element transfers, refocusing phases and aberration metrics of one transducer |
| `time_reversal_stage`                   | Refocusing phases, aberration metrics and pose search from one reverse run per target |
| `write_acoustic_sidecar`                | Describe the T1w pressure NIfTI so a later run can start from it |

#### BABELBRAIN

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `babelbrain_check_geometry`             | Error for geometries the BabelBrain backend cannot simulate (several transducers, beam off the grid z axis) |
| `babelbrain_insitu_isppa`               | In-brain peak ISPPA of the scaled acoustic field — the `Isppa` BabelBrain's thermal solver normalises to. |
| `babelbrain_placement_correction`       | Mechanical Tx offset that puts BabelBrain's apex at trans_pos |
| `babelbrain_thermal_timing`             | Convert PRESTUS timing to BabelBrain thermal parameters |
| `babelbrain_thermal_wrapper`            | Run BabelBrain thermal simulation and read results |
| `babelbrain_transducer_map`             | Map a PRESTUS transducer to a BabelBrain preset |
| `babelbrain_wrapper`                    | Run BabelBrain acoustic simulation as PRESTUS backend |
| `export_babelbrain_inputs`              | Write all BabelBrain input files from a PRESTUS run |
| `export_babelbrain_weights`             | Pack PRESTUS per-element phases and particle velocity into a BabelBrain OptimizedWeightsFile (HDF5 complex64 phasors). |
| `import_babelbrain_weights`             | Read per-element phase weights from a BabelBrain OptimizedWeightsFile (.h5), derive hardware correction relative to geometric phases, and write the correction into the transducer equipment YAML. |
| `read_babelbrain_output`                | Read BabelBrain DataForSim.h5 into PRESTUS sensor_data |
| `read_babelbrain_thermal`               | Read BabelBrain thermal HDF5 into PRESTUS results_heating |
| `resample_babelbrain_vol`               | Resample a BabelBrain volume onto the PRESTUS acoustic grid |
| `resolve_bb_backend`                    | Resolve 'auto' BabelBrain computing_backend / device_name |
| `write_babelbrain_input`                | Write BabelViscoInput NIfTI from PRESTUS medium masks |

#### CALIBRATION

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `axial_transfer_matrix`                 | On-axis complex pressure per unit element drive |
| `build_steering_correction`             | Steering correction of a matrix array from hydrophone measurements |
| `calibrate_babelbrain_phase_profile`    | Whole-profile Rayleigh fit of annular drive phases for BabelBrain |
| `calibrate_focal_steering`              | Analytical steer that lands the FWHM/peak at a target |
| `calibrate_library_entry`               | Calibration core: profiles at requested depths (measured or interpolated), per-depth fits, refinements, optional k-Wave drive calibration, report (no file I/O) |
| `calibration_desired_profiles`          | Profiles on a desired-depth grid, for prediction and fitting |
| `calibration_evaluate`                  | Geometric, library, fitted and k-Wave profiles on one set of depths |
| `calibration_fit`                       | Interactive calibration fit: per-depth refinements |
| `calibration_fit_depths`                | One `fit_transducer_model` fit per depth of a profile set (shared by library build, refinement, GUI, `calibration_transducer`) |
| `calibration_fit_entry`                 | Library entry from the current entry and a CALIBRATION_FIT result |
| `calibration_fit_options`               | Per-depth fit presets: `rayleigh_phase_gain` (default) and `oneil_phase_only` (former global-search settings) |
| `calibration_focus_beyond`              | FLHM centre of each axial profile, ignoring the near field |
| `calibration_pipeline_start`            | Universal launcher for calibration_transducer (MATLAB / SLURM / qsub) |
| `calibration_save_entry`                | Write a fitted library entry to a user library folder |
| `calibration_setup`                     | Run CALIBRATION_TRANSDUCER for every equipment combination in the calibration config |
| `calibration_transducer`                | Config-driven calibration at chosen depths and intensities; writes a library entry and a calibration output folder |
| `check_equipment_config`                | Validate equipment configuration YAML files |
| `compute_flhm_center`                   | Compute the FLHM (Focal Length at Half Maximum) centre/span of an axial intensity profile at a configurable threshold fraction of peak (0.5 = −3 dB default, 0.25 = −6 dB). |
| `compute_phases`                        | Calculate per-element transducer phases (degrees) for electronic focusing by computing fractional wavelength delays to steer waves constructively to a target point relative to the exit plane. |
| `evaluate_library_entry`                | Model the axial profiles a study gets from a library entry, at any depths |
| `extract_real_intensity_profile`        | Extract or spline-interpolate axial intensity profiles between available measured focal depths (wrt exit plane), aligns peaks, skips near-field for max, and plots/saves results. |
| `field_transfer_matrix`                 | Complex pressure per unit element drive at points off the axis |
| `fit_ring_profile`                      | Fit a within-ring velocity profile shared by all focal depths |
| `fit_steering_correction`               | Fit measured / modelled focal peak over a matrix array's steering range |
| `fit_transducer_model`                  | Fit a static annular-array model to hydrophone axial profiles |
| `focal_reference_to_frac`               | Map a `focal_reference` config value (`flhm_3db`/`flhm_6db`) to its FLHM threshold fraction of peak; errors on `peak` (no associated fraction) or an unknown value. |
| `generate_calibration_report`           | Write an HTML summary of a library entry: provenance, fit quality, settings at the fitted depths, amplitude lookup. |
| `interpolate_axial_profiles`            | Axial intensity profiles at desired focal depths from measured ones |
| `kwave_calibrate_drive`                 | Calibrate the drive of per-depth refinements in free-water k-Wave (rescale to the simulated ISPPA) and report agreement with the measured profiles |
| `kwave_check_refinement`                | Deprecated name of KWAVE_CALIBRATE_DRIVE |
| `kwave_drive_metrics`                   | Absolute agreement of a simulated free-water profile with its target |
| `library_add_conditions`                | Add dephased drive conditions to a library entry [EXPERIMENTAL] |
| `library_add_refinements`               | Merge depth-specific phase offsets (and drives) into a library entry |
| `library_condition`                     | Phases, gains and drive of a dephased library condition [EXPERIMENTAL] |
| `library_entry_from_model`              | Write a FIT_TRANSDUCER_MODEL result into a library struct |
| `library_ingest_conditions`             | Add dephased conditions from a profile CSV + sidecar [EXPERIMENTAL] |
| `library_model_at_depth`                | Phases, gains and drive of a library entry at any focal depth |
| `library_refinement_from_fits`          | Write per-depth FIT_TRANSDUCER_MODEL results as library refinements |
| `library_set_ring_profile`              | Set or remove the within-ring velocity profile of a library entry |
| `library_settings_table`                | Per-depth transducer settings a library entry resolves to, as a table / CSV |
| `library_static_model`                  | Global model of a library entry (geometric, or a legacy static offset) in FIT_TRANSDUCER_MODEL form; the prior of per-depth fits |
| `library_transducer`                    | Transducer geometry a library entry was calibrated with |
| `load_equipment_config`                 | Load equipment configuration from individual YAML files |
| `load_transducer_from_library`          | Load calibrated transducer parameters from library |
| `matrix_element_layout`                 | Matrix-array elements in the transducer frame (no grid) |
| `matrix_focal_peak`                     | Modelled free-field focal peak of a matrix array per steering target |
| `matrix_rayleigh_pressure`              | Free-field complex pressure of a matrix array (Rayleigh integral) |
| `perform_global_search`                 | Globally optimize transducer element phases (rad) and particle velocity via FEXminimize or GlobalSearch to minimize error fitting O’Neil analytical intensity to measured axial profiles. |
| `phase_optimization_annulus_full_curve` | Optimizes the phase profile for an annular transducer by matching intensity curves. |
| `plot_calibration_report`               | Measured vs library vs calibrated axial profiles per depth (PNG) |
| `plot_ring_profile`                     | Show a fitted within-ring velocity profile and what it changes |
| `plot_transducer_library`               | Plot the ring phases, drive and ring gains a library entry resolves to across focal depths (as at run time; all library schemas). |
| `rayleigh_axial_intensity`              | On-axis intensity via the Rayleigh–Sommerfeld integral |
| `read_axial_profiles`                   | Read a hydrophone axial-profile CSV (one column per focal depth) |
| `read_steering_measurements`            | Read hydrophone focal peaks of a steered matrix array |
| `refine_transducer_library`             | Optional refinement of a library entry to your own hydrophone profiles |
| `resolve_transducer_from_serial`        | Populate transducer struct from equipment and calibration library |
| `ring_profile_drive`                    | Simulation drive of an annular transducer with a ring profile |
| `ring_profile_transducer`               | Annular transducer with a within-ring velocity profile |
| `save_elem_correction`                  | Compute and persist per-element phase correction into the equipment YAML |
| `save_library_entry`                    | Write a user-calibrated library entry with provenance |
| `save_optimized_values`                 | Append optimized transducer phases and amplitudes into a CSV table (indexed by intensity/focus) and save parameters as YAML for PRESTUS config integration. |
| `scan_and_ingest_calibrations`          | Scan a folder and ingest all calibration YAMLs |
| `set_real_phases`                       | Load or compute manufacturer-specific (Sonic Concepts/Imasonic) transducer element phases for given focal depth (wrt exit plane), interpolate/unwrap for virtual elements, and return phase degrees. |
| `steering_correction_factor`            | Measured / modelled focal peak at steering targets |
| `steering_path_lengths`                 | Ring-centre to target distances used for geometric steering |
| `steering_quad_terms`                   | Design matrix of a full quadratic surface in three coordinates |
| `update_transducer_library`             | Add depth-specific refinements to a library YAML |
| `validate_transducer_library`           | Validate calibration library against equipment configs |
| `water_baseline`                        | Measure free-water ISPPA at the current source amplitude |
| `write_axial_profiles`                  | Write axial profiles as a profile CSV plus YAML sidecar |

#### CORE

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `apply_isppa_scaling`                   | Scale the simulated pressure field so the free-water |
| `async_transducer_pipeline`             | Simulate N asynchronously firing transducers |
| `combine_async_intensity`               | Scale and combine N independent acoustic results |
| `grid_stage`                            | Build or load the head model grid cache (Stage 3). |
| `is_async_mode`                         | Return true when async multi-transducer mode is active. |
| `is_multi_isppa_mode`                   | True when any transducer carries more than one ISPPA target. |
| `is_pose_workflow_mode`                 | True when time_reversal.pose_workflow.enabled is set |
| `isppa_baseline`                        | Free-water ISPPA an acoustic field was produced at, and its origin |
| `isppa_intent`                          | One-line statement of the free-water ISPPA a run will model |
| `isppa_scale_factor`                    | Pressure scale factor that brings a field to a target free-water ISPPA |
| `load_grid_cache_file`                  | Load a PRESTUS grid cache and restore pipeline state. |
| `load_parameters`                       | Loads and merges configuration files for simulation parameters; a trailing `{key, value}` cell applies overrides before validation. |
| `load_transducer_parameters`            | Validate, normalize, and derive transducer parameters from config |
| `medium_stage`                          | Build or load the acoustic/thermal medium property cache (Stage 4). |
| `multi_isppa_pipeline`                  | Run acoustic simulation once, then thermal analysis at |
| `path_log_setup`                        | Set up internal paths and logging, filename for output table |
| `pose_workflow_pipeline`                | Transducer pose search, forward verification and final run in one call |
| `posthoc_water_parameters`              | Parameters of one post-hoc free-water re-entry run |
| `prestus_config_init`                   | Bootstrap a project-specific PRESTUS configuration folder |
| `prestus_config_upgrade`                | Refresh a project's config_default.yaml from the toolbox |
| `prestus_group_report_start`            | Entry point for the group HTML report and group MNI maps: study mode (runs from `study.manifest`) or subject mode (auto-discovered / listed subjects). |
| `prestus_group_start`                   | DEPRECATED entry point for group-level MNI-space analysis |
| `prestus_pipeline`                      | End-to-end transcranial focused ultrasound simulation pipeline |
| `prestus_pipeline_start`                | Universal PRESTUS pipeline launcher for MATLAB, SLURM, and qsub |
| `prestus_study_start`                   | Run the pipeline for every row of a study manifest (subjects × conditions with per-run settings), locally or on HPC, then build the study report; `DryRun` checks all rows without running. |
| `print_parameter_summary`               | Print a structured PRESTUS configuration summary to the console |
| `sequential_pipeline`                   | Dispatch follow-up simulations that chain thermal state |
| `uncertainty_pipeline`                  | Run the full PRESTUS uncertainty quantification workflow |

#### GROUP

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `changem_vectorized`                    | Replace multiple old values in array A with corresponding new values. |
| `combine_plots_by_suffix`               | Combines subject-specific plots into a single montage image. |
| `create_group_maps`                     | Group MNI maps per condition from the per-run MNI NIfTIs: −3 dB focal-volume overlap, mean/SD intensity, mean heating, plus overlap figures. |
| `create_group_MNI_plots`                | Deprecated (older output layout): group plots in MNI space; use `create_group_maps`. |
| `discover_group_subjects`               | Auto-discover completed subjects by scanning `path.sim` for folders with a matching per-subject CSV. |
| `study_read_manifest`                   | Read a study manifest (.tsv/.csv): one row per run with reserved columns and dotted config-key columns. |
| `study_row_parameters`                  | Load the full config of one study run (study config < row config < manifest columns < overrides). |
| `study_run_descriptors`                 | Per-run settings that differ between subjects (placement mode, target, transducer, depth, target intensity). |
| `study_run_files`                       | Output file names of one run (CSV, targeting sidecar, report, NIfTIs) for the group tools. |
| `study_run_table`                       | Resolve manifest rows to output locations (subject, affix from condition, medium, path.sim); rejects rows that would overwrite each other. |
| `study_setup`                           | Load a study config, read its manifest and resolve the runs and study output folder. |

#### GUI

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `calibration_gui`                       | Standalone GUI to fit transducer library entries to hydrophone profiles, predict profiles at any depth and add transducers. |
| `gui_run_worker`                        | Run the PRESTUS pipeline from a YAML config file |
| `neuronav_gui`                          | Stand-alone GUI for Localite marker ingestion: detect and review series, QC, name targets, batch ingestion, MNI view. |
| `nexus_bridge`                          | Persistent MATLAB session that serves the PRESTUS Nexus |
| `nexus_calibration`                     | Transducer calibration for the PRESTUS Nexus (through NEXUS_BRIDGE) |
| `nexus_neuronav`                        | Localite marker ingestion for the PRESTUS Nexus (through NEXUS_BRIDGE) |
| `nexus_ping`                            | Answer a NEXUS_BRIDGE health check |
| `nexus_progress`                        | Report progress of the current NEXUS_BRIDGE request to the Nexus |
| `prestus_gui`                           | Point-and-click GUI to review and edit a PRESTUS config and start runs. |
| `prestus_nexus`                         | Open the PRESTUS Nexus (standalone PyQt6 config builder, *experimental*; see [doc_nexus.md](doc_nexus.md)). |

#### HEAD

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `autopad_grid_for_transducers`          | Expand grid if any transducer/focus position falls outside |
| `head_smooth_and_crop`                  | Convert and smooth segmentations into medium map & crop grid for effiency |
| `mean_beam_axis`                        | Mean axis of several transducer beams, averaged as axes (opposed beams do not cancel) |
| `placement_run_dir`                     | Folder of one placement run (one mode, one target) |
| `placement_signature`                   | Settings a placement run's result depends on |
| `preproc_align_to_focal_axis`           | Rotate input image and its grid coordinates to the transducer-focal axis (+Z) |
| `preproc_align_to_mean_focal_axis`      | Rotate a 3D image so the mean focal axis aligns with z |
| `preproc_crop_eCSF`                     | Expand CSF to define outer edges of layered medium. |
| `preproc_crop_grid`                     | Crop grid to head + transducer + PML. |
| `preproc_deface_m2m`                    | Create a defaced copy of a SimNIBS m2m folder |
| `preproc_head`                          | Preprocesses structural brain data for simulations (segmentation, alignment, cropping). |
| `preproc_medium_mask`                   | Map segmentation indices onto the medium labels (in the config) [layered only] |
| `preproc_rescale_only`                  | Rescale a 3D image to grid resolution without rotation |
| `preproc_segmentation`                  | Setup SimNIBS segmentation |
| `preproc_transducer_placement`          | Resolve transducer and focus positions |
| `read_placement_record`                 | Read the run record of a placement run folder |
| `segmentation_run`                      | Submit SimNIBS segmentation |
| `skull_fill_holes`                      | Fill holes in skull segmentation and between skull and skin |
| `skull_rubber_wrap`                     | Inflate the skull layer to locally fill potential holes |
| `skull_rubber_wrap_visualize`           | Visualize results of skull rubber expansion |
| `smooth_img`                            | Apply 3D smoothing |
| `write_placement_positions`             | Write the resolved transducer positions as a JSON sidecar |
| `write_placement_record`                | Write the run record (placement_run.json) of a placement run |

#### HELPER

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `apply_config_overrides`                | Set config values by dotted key (`placement.localite.target`, `transducer(2).serial`); warns on unknown keys. |
| `bb_tissuemask_binary`                  | Build per-tissue binary masks from a BabelBrain material map |
| `bwconncomp_3d`                         | Find connected components in a 3-D binary volume |
| `cache_signature`                       | Settings that determine a cached pipeline stage's output |
| `cast_struct`                           | Recursively cast all numeric fields in a struct to a given type |
| `charm_seg_labels`                      | Return the fixed SimNIBS charm tissue label indices |
| `check_availability`                    | Check file availability. |
| `check_cache_signature`                 | Refuse to reuse a cache built from different settings |
| `check_layers`                          | Match requested layers to those available in the segmentation. For (p)CT, homogenize multi-skull layer into a single `skull `layer. |
| `check_sidecar_signature`               | Compare a sidecar's recorded settings with the current run's |
| `compute_ramp_factor`                   | Pulse-average / temporal-peak intensity ratio for a ramped pulse |
| `config_known_extra_keys`               | Config keys PRESTUS reads that config_default.yaml omits |
| `config_unknown_keys`                   | Dotted keys of a user config that PRESTUS does not know |
| `confirm_overwriting`                   | Check overwriting (manual). |
| `confirmation_dlg`                      | Present confirmation dialogue. |
| `default_log_file`                      | Timestamped pipeline log (diary) path for a run |
| `ensure_hdf5_plugin_path`               | Set HDF5_PLUGIN_PATH so MATLAB's h5read/h5write don't fail |
| `find_min_factor`                       | Find the number with the smallest maximum factor in a range. |
| `find_pseudoct_dir`                     | Folder holding this subject's existing pseudoCT.nii.gz |
| `find_sidecar`                          | JSON sidecar of a data file: same name with .json (BIDS convention) |
| `get_crop_dims`                         | Computes cropping dimensions for a 3D image with a margin. |
| `get_flhm_center_position`              | Calculates the center position of the full-length half-maximum (FLHM). |
| `get_output_dir`                        | Derive the simulation output directory from parameters |
| `get_pct_dir`                           | Per-subject pseudo-CT folder |
| `get_placement_dir`                     | Per-subject root of the placement outputs |
| `get_prestus_path`                      | GET_PRESTUSPATH  Return the absolute path to the PRESTUS toolbox root |
| `get_slice_by_label`                    | Extracts a specific slice from a 3D image based on axis label and slice number. |
| `get_xyz_mesh`                          | Generates a mesh of 3D coordinates for a given image. |
| `getidx`                                | Retrieves indices for requested tissues from a parameter structure. |
| `iff`                                   | Inline conditional: return a if cond is true, else b. |
| `imfill_holes_3d`                       | Fill interior holes in a 3-D binary volume |
| `kwave_version`                         | Display k-Wave version number and (if available) git hash. |
| `log_memory_usage`                      | Print this MATLAB process's RSS and the enclosing cgroup's usage. |
| `log_timer`                             | Start or stop logs for benchmarking time, RAM, and disk space use. |
| `manual_coordinate_space`               | Coordinate space of transducer.trans_pos / focus_pos |
| `masked_max_3d`                         | Computes the maximum intensity within a masked 3D region. |
| `matlab_escape_quotes`                  | Escapes text for a single-quoted MATLAB literal in generated job scripts (doubles `'`). |
| `mergestruct`                           | Recursively merge fields and subfields of two structures; treats `yaml.Null` as empty so YAML-loaded null values are always overwritten. |
| `nearest_fill_nan`                      | Replace NaN voxels with the value of the nearest non-NaN voxel |
| `outputs_follow_recompute`              | Let a recomputed stage replace its existing output files |
| `parse_memory_tags`                     | Extract peak memory figures from [memory]/[memory_bb] log lines. |
| `predict_outputs`                       | Predict the pipeline mode and the CSV/HTML files a run will write, before submission (read-only). |
| `prestus_version`                       | Print PRESTUS version info and optionally return the git hash. |
| `print_predicted_outputs`               | Print `predict_outputs`' plan to the console. |
| `r_quote_string`                        | Escape a string for safe embedding in a generated R string literal. |
| `read_ini_file`                         | Read an INI file into MATLAB |
| `read_sidecar`                          | Read a stage-output JSON sidecar |
| `resolve_chained_overwrite`             | Map io.overwrite_files onto the stages of a |
| `resolve_io_dirs`                       | Fill missing io path fields from path.sim + subject_id |
| `resolve_jwave_pml_size`                | Resolve grid.pml_size 'auto' to a wavelength-based scalar PML thickness for jwave |
| `resolve_pml_size`                      | Resolve grid.pml_size to a numeric array for k-Wave |
| `resolve_python_env`                    | Locate a Python executable inside a venv/conda environment. |
| `reuse_upstream_cache`                  | True when an existing upstream cache must be loaded |
| `round_if_integer`                      | Rounds values if they are sufficiently close to integers; otherwise raises an error. |
| `run_command`                           | Execute an external program via java.lang.ProcessBuilder. |
| `run_shell_script`                      | Execute a shell script string via ProcessBuilder, with only ONE shell in the loop. |
| `run_stats`                             | Collect a run's simulation size, duration and peak memory |
| `safe_addpath`                          | Add a directory tree to the MATLAB path, skipping hidden directories |
| `sample_beam_axis`                      | Interpolate a 2D/3D field along the beam axis from the bowl apex to the grid edge (distances from the apex). |
| `sanitize_output_affix`                 | Normalise `io.output_affix` to the BIDS `_desc-<label>` form. |
| `shell_quote_arg`                       | Quote a string as a single literal argument for system(). |
| `should_save_output`                    | Decide whether an output is saved: per-step `io` flag, else `io.save_matrices` (matrix outputs), else the built-in default. |
| `simnibs_version`                       | Get SimNIBS version of segmentation from HTML, print, allocate to `parameters`. |
| `subset_fields`                         | Copy designated structure field to new structure. |
| `tissuemask_binary`                     | Extract binary tissue segmentation masks (for indexing) |
| `trans_suffix`                          | Return '_trans-XX' when n_trans > 1, otherwise ''. |
| `upsample_to_grid`                      | Upsamples a 3D image to a higher resolution grid. |
| `validate_any_required_field`           | Error clearly unless at least one candidate field is set. |
| `validate_required_fields`              | Error clearly unless every required field is set. |
| `warn`                                  | Wrapper around warning() that suppresses the stack-trace printout. |
| `write_sidecar`                         | Write a JSON sidecar describing a stage output later runs can start from |
| `zip_fields`                            | Convert a structure's fields and values into a cell array. |

#### HPC

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `hpc_add_dependency`                    | Add scheduler job IDs to the jobs a submission waits for |
| `hpc_detect_system`                     | Detect SLURM or qsub HPC system. |
| `hpc_job_info`                          | Generate formatted job display information. |
| `hpc_job_name`                          | Generate standardized HPC job name. |
| `hpc_job_report`                        | Query scheduler for post-job resource usage and carbon metrics |
| `hpc_setup_temp_files`                  | Setup directories and generate temporary files. |
| `hpc_submit_job`                        | Submit HPC batch job (SLURM or qsub). |
| `hpc_submit_report_job`                 | Submit a lightweight dependent job to generate the post-run report |
| `hpc_validate_parameters`               | Validate HPC job parameters. |
| `hpc_wait_for_completion`               | Poll an HPC job until it completes or times out |
| `parse_sbatch_job_id`                   | Job ID from the output of sbatch |

#### JWAVE

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `jwave_wrapper`                         | Run acoustic simulation via jwave (JAX-based Python backend) |

#### LEGACY

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `config_deprecations`                   | Config keys renamed or removed since PRESTUS v0.5.0 |
| `upgrade_config_keys`                   | Convert config keys from PRESTUS v0.5.0 onwards to their current names |

#### MEDIUM

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `fitPowerLawParamsMulti`                | Fit power law absorption parameters for highly absorbing media. |
| `get_alpha_coeff`                       | Computes the amplitude attenuation coefficient for a given medium and frequency. |
| `medium_elastic_properties`             | Build the k-Wave elastic medium (compressional + shear) from a fluid medium |
| `medium_pct_attenuation`                | Skull: pCT-informed attenuation mapping |
| `medium_pct_density`                    | Skull: pCT-informed density mapping |
| `medium_pct_density_calibration`        | CT->density calibration table for pct.mapping_density 'k-plan' |
| `medium_pct_pichardo_load`              | Load and cache Pichardo 2D interpolant from BabelBrain |
| `medium_pct_soundspeed`                 | Skull: pCT-informed sound speed mapping |
| `medium_setup`                          | Set up medium |

#### NEURONAV

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `export_brainsight_trajectory`          | Write a Brainsight version-14 trajectory file. |
| `localite_matrix_to_positions`          | Convert a Localite 4×4 matrix to voxel positions. |
| `neuronav_apply_deformation`            | Map segmentation-space RAS coordinates through a SimNIBS coordinate field |
| `neuronav_build_deformation_cache`      | Precompute nearest-valid-voxel lookup for a SimNIBS coordinate field |
| `neuronav_build_target_map`             | Expand a coil_map spec into a flat target_map. |
| `neuronav_check_target_map`             | Check a target_map before NEURONAV_INGEST_MARKERS writes it |
| `neuronav_compute_series_statistics`    | Compute mean position and variability over stimulus train |
| `neuronav_compute_series_statistics_sequential`| NEURONAV_COMPUTE_SERIES_STATISTICS  Compute per-series position and orientation statistics |
| `neuronav_convert_MNI_to_native`        | Transform coordinates from MNI space to native subject space |
| `neuronav_convert_native_to_MNI`        | Convert native RAS coordinates to MNI space |
| `neuronav_convert_trigger_to_voxels`    | Convert Localite trigger positions to voxel (image) coordinates. |
| `neuronav_detect_series`                | Detect the recorded series of one session (QC, suggested names) without writing. |
| `neuronav_export_session_csv`           | Export per-session coordinate arrays to CSV with labeled voxel and RAS (mm) positions. |
| `neuronav_find_mni_template`            | SimNIBS' MNI152_T1_1mm.nii.gz, for head outlines in MNI views |
| `neuronav_find_planning_t1`             | Find the planning T1 loaded into Localite (`path.localite_t1`, else `path.localite_raw`). |
| `neuronav_get_group_mean_mni`           | Compute group-level average MNI coordinates (both mm and voxel). |
| `neuronav_head_surface`                 | Outer head surface of an MNI-space NIfTI, for plots |
| `neuronav_ingest_batch`                 | Ingest and export Localite targets for many subjects/sessions; writes `neuronav_ingest_summary.csv`. |
| `neuronav_ingest_markers`               | Ingest raw Localite markers and write normalised outputs. |
| `neuronav_label_map_from_scan`          | Build a `neuronav.target_labels` map from a label scan. |
| `neuronav_localite_pose`                | Resolve the recorded Localite instrument pose for placement (file mode, or ingested markers; ingests on first use). |
| `neuronav_map_label`                    | Map a Localite target label to the study target name (`neuronav.target_labels`). |
| `neuronav_match_entry_targets`          | Match Localite EntryTarget entries to a target map. |
| `neuronav_name_series`                  | Name the detected series of a session (coil table, planned target, default; shared GUI/script rule). |
| `neuronav_plot_placement`               | Visualise transducer placements on a segmented head. |
| `neuronav_read_entry_target`            | Read planned target positions from Localite EntryTarget XML. |
| `neuronav_read_label_candidates`        | Collect the named planned targets of a session (EntryTarget, GUMMarkers). |
| `neuronav_read_target_lookup`           | Look up pre-specified target coordinates from |
| `neuronav_resolve_paths`                | Resolve the Localite raw/ingested/T1/export folders from `path.localite_*` (with `path.localite` fallback). |
| `neuronav_resolve_target`               | Resolve each target_map entry's brain-target |
| `neuronav_scan_sessions`                | Inventory raw Localite data and planned-target labels. |
| `neuronav_select_localite`              | Select most recent Localite XML for a session. |
| `neuronav_series_qc_flag`               | Tracking problem of one detected series, as text |
| `neuronav_suggest_target_names`         | Heuristic target names for ingested series (planned target nearest the beam axis). |
| `neuronav_summary_to_mni`               | Add MNI coordinates to a neuronav ingest summary. |
| `position_transducer_localite`          | Determines transducer and focus positions in voxel space using Localite data and MRI header. |
| `position_transducer_plantus`           | Optimise transducer placement using PlanTUS. |
| `position_transducer_tus_entry`         | Determine transducer placement using TUS_entry. |

#### NIFTI

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `nifti_acoustic`                        | Export acoustic-stage data to NIfTI in T1 and MNI space |
| `nifti_medium`                          | Export medium-stage data to NIfTI |
| `nifti_medium_masks`                    | Write medium_masks pipeline output to NIfTI |
| `nifti_medium_properties`               | Write medium acoustic/thermal property maps to NIfTI |
| `nifti_pct`                             | Write pseudo-CT pipeline output to NIfTI |
| `nifti_thermal`                         | Export thermal-stage data to NIfTI in T1 and MNI space |
| `nifti_to_mni`                          | Convert a T1-space NIfTI to MNI space via SimNIBS |
| `nifti_to_t1w`                          | Back-transform one volume to T1 space and write NIfTI |
| `phantom_voxel_order`                   | Return a non-layered output volume to the voxel order of its input |
| `simnibs_t1_to_mni`                     | Register subject T1 to MNI space if not already done |

#### PSEUDO-CT

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `check_pct_sidecar`                     | Refuse an existing pseudo-CT made with other settings |
| `fit_pairwiselinear`                    | Perform a pairwise linear fit between HU and density values with optional plot. |
| `pct_create_pCT_run`                    | Submit a pseudoCT generation job |
| `pct_signature`                         | Settings a pseudo-CT depends on, and the resolved skull mapping |
| `pct_skullexpand`                       | Load SimNIBS NIfTIs + headers, run skull rubber wrap, save output. |
| `pct_skullmapping`                      | Computes pseudo-CT mapping for cortical and trabecular bone using UTE histograms. [deprecated, debug] |
| `pct_skullmapping_fit`                  | Subject-specific UTE->HU line from cortical / trabecular UTE histograms |
| `pct_soft_tissue_peak`                  | Identifies the soft tissue peak from UTE intensity distribution histograms. |

#### PLOT

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `deface_volume`                         | Remove face voxels in a volume using the SimNIBS segmentation |
| `nifti_slice_orientation`               | Compute display transforms for a NIfTI orthogonal slice. |
| `plot_coronal_slices`                   | (DEPRECATED) Visualizes coronal slices of a 3D image with optional legends for labeled images. |
| `plot_intensity_t1_overlay`             | Overlay back-transformed intensity on T1 image |
| `plot_median_montage`                   | (DEPRECATED) Creates a montage of central slices from a 3D T1 image. |
| `plot_overlay`                          | Visualizes map overlaid on a 2D slice of a 3D background image with advanced options. |
| `plot_overlay_2d`                       | Overlays map on a background image slice with key positions highlighted. |
| `plot_placement_t1_overlay`             | Quick placement QC figure on the native T1 |
| `plot_t1_with_transducer`               | Creates a plot of a T1 slice oriented along the transducer's axis with overlays. |
| `plot_transducer_overlay`               | Visualize curved transducer geometry with exit plane. |
| `show_3d_head`                          | Visualizes segmented brain images in 3D with transducer placement and target location. |
| `show_binary_mask_transducer`           | Visualise the binary mask of a matrix transducer on the simulation grid |
| `show_positioning_plots`                | Visualizes transducer positioning before and after preprocessing using segmented images. |

#### REPORT

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `build_config_dump`                     | Searchable HTML block of the full parameters struct. |
| `compute_tic`                           | Cranial Thermal Index (TIC) after IEC 62359 / AIUM-NEMA ODS. |
| `css_styles_base`                       | Shared design system for all PRESTUS HTML reports. |
| `generate_group_report`                 | Assemble the self-contained group HTML report (exposure dashboard, study design, condition comparison, targeting accuracy, group maps, outliers, roster, box plots, per-run cards) and a one-row-per-run CSV. |
| `generate_multi_isppa_report`           | Generate an HTML report comparing multiple ISPPA levels |
| `generate_placement_report`             | Stand-alone HTML placement report written after the placement stage (positions, T1 overlay, PlanTUS criteria evaluation); `placement.report`. |
| `generate_sequential_report`            | HTML summary report for a chain of sequential simulations |
| `generate_simulation_report`            | Generate a self-contained HTML simulation report |
| `generate_uncertainty_report`           | Generate an uncertainty-range HTML simulation report |
| `get_risk_limits`                       | Return ITRUSST consensus non-significant risk limits |
| `html_utils`                            | Static HTML utility methods for PRESTUS report generation |
| `preview_report_styles`                 | Report test mode: re-render an existing subject output folder in every skin (no re-simulation) into a sandbox folder with a side-by-side index.html. |
| `report_scripts`                        | Shared inline JavaScript for PRESTUS HTML reports. |
| `report_styles`                         | Alternative visual skins for the HTML reports (instrument / classic / compact / paper), selected via `modules.report_style`. |
| `risk_color`                            | Return a colour string based on a value relative to an ITRUSST risk limit |
| `table2html`                            | Convert a MATLAB table to an HTML table element |

#### SOURCE

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `create_clover_array`                   | Replicate a single matrix sub-array into a multi-leaf clover configuration |
| `create_matrix_karray`                  | Add matrix transducer elements to a kWaveArray object |
| `CWSourceMatrix`                        | Lazy rank-2 factored CW source signal matrix for k-Wave. |
| `grid_axisymmetry`                      | Convert grid to axisymmetry |
| `grid_tissue_setup`                     | Set up grid dimensions, preprocess head (if modeled), and place in grid |
| `grid_transducer_location`              | Position transducer (and target) in simulation grid |
| `source_create`                         | Creates ultrasound source signals and masks for k-Wave simulations based on transducer geometry. |
| `source_sensor_setup`                   | Configure a k-Wave simulation grid, transducer sources, and sensor mask based on input parameters for recording pressure fields in 2D or 3D ultrasonic neuromodulation setups. |
| `water_grid_dims`                       | Free-water grid size from the transducer aperture and focal depth |

#### TELEMETRY

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `generate_run_id`                       | Generate a random UUID v4 string for telemetry. |
| `prefdir_prestus`                       | Return (and create) the PRESTUS per-user preferences directory. |
| `telemetry_env_override`                | Read the PRESTUS_TELEMETRY environment variable. |
| `telemetry_has_display`                 | True when MATLAB can render a figure/dialog window. |
| `telemetry_is_headless_cluster`         | True for compiled apps and SLURM/PBS jobs. |
| `telemetry_is_noninteractive`           | True when no telemetry prompt may be shown. |
| `telemetry_run_type`                    | Resolve whether this run is 'technical' or 'deployment'. |
| `telemetry_set_consent`                 | Record a telemetry opt-in or opt-out decision. |
| `telemetry_setup`                       | Show the opt-in prompt for anonymous usage statistics. |
| `telemetry_setup_reset`                 | Remove the telemetry consent file so the opt-in |
| `track_usage`                           | Send anonymous usage statistics to the PRESTUS telemetry endpoint. |

#### THERMAL

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `cem43_iso_increment`                   | ISO thermal-dose (CEM43) increment for one time step |
| `thermal_absorption_model`              | Resolve how acoustic heat deposition is specified |
| `thermal_analysis`                      | Analyze output of thermal simulations |
| `thermal_grid_setup`                    | Build an independent k-Wave grid for thermal simulation |
| `thermal_parameters`                    | Checks and converts protocol timing setup for thermal estimation |
| `thermal_plot_protocol`                 | Visualize the requested protocol timing |
| `thermal_plot_sim`                      | Visualizes heating simulation results over time (temperature, rise, CEM43). Optional video |
| `thermal_simulation`                    | Simulate ultrasound-induced heating and thermal dose (CEM43) using k-Wave’s kWaveDiffusion from acoustic pressure data across pulsed train repetitions. |
| `thermal_stage`                         | Run or load the bio-heat thermal simulation (Stage 8). |
| `thermal_update_timeseries`             | Track max T/CEM43 per tissue layer. |

#### TRANSDUCER

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `annular_aperture_points`               | Sample points on the radiating surface of a curved annular array |
| `annular_bowl_depth_mm`                 | Distance from the transducer position (bowl apex) to the exit plane |
| `annular_geometric_phases`              | Ring phases that steer a curved annular array to a point, in water |
| `convert_to_element_pos`                | Generate element positions for matrix transducers |
| `extract_element_pos`                   | Read matrix transducer element positions from file |
| `focal_distance_calculation`            | Compute expected focal distances for (multi-)transducer setup. |
| `get_arc`                               | Generates the coordinates of an arc in 2D space. |
| `get_trans_pos_from_trigger_markers`    | Determines transducer and target positions using Localite trigger marker files. |
| `get_transducer_box`                    | Computes transducer box dimensions and positions for simulations. |
| `is_flat_matrix`                        | True for a flat matrix array (is_curved false or infinite radius) |
| `matrix_elem_diameter`                  | Aperture diameter of a disc/bowl matrix element (equal-area or `elem_diameter_mm`). |
| `matrix_element_aim`                    | Point each matrix element faces (its normal), curved or flat |
| `reciprocity_pose_candidates`           | Transducer poses around the planned one, for scoring on one reverse field |
| `reciprocity_score_poses`               | Score transducer poses on one reverse field (acoustic reciprocity) |
| `tp_candidate_mesh`                     | Build coordinate mesh and candidate transducer geometry. |
| `tp_evaluate_candidate_positions`       | Evaluate candidate transducer positions. |
| `tp_find_initial_candidate`             | Find candidate transducer positions on skull surface. |
| `tp_plot_candidate_positions`           | Plot transducer candidate positions on skull surface slice. |
| `tp_plot_geometry_overlay`              | Visualize transducer geometry on skin segmentation slice. |
| `tp_plot_heuristic_position`            | Create visualization of heuristic transducer placement. |
| `tp_remove_ear_locations`               | Exclude ear entries from heuristic transducer positions. |
| `tp_scalp_boundary`                     | Outer (scalp) boundary voxels of a segmented head volume |
| `tp_select_heuristic_position`          | Select optimal transducer position and export Localite coordinates. |
| `transducer_analyze_position`           | Computes geometric and statistical measures for a transducer position. |
| `transducer_analyze_position_fast`      | Analyzes the transducer position relative to the skull and skin. |
| `transducer_positioning`                | Determines heuristic transducer placement for a given target. |
| `transducer_positioning_start`          | Dispatch transducer positioning to MATLAB or HPC |
| `transducer_scalp_geometry`             | Place a transducer on the scalp with a standoff gap |
| `transducer_setup`                      | Create a transducer mask and label matrix for a computational grid. |
| `validate_annular_transducer`           | Validate configuration fields of an annular transducer |
| `validate_matrix_transducer`            | Validate and normalise configuration of a matrix transducer |

#### TRANSFORM

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `affine_apply_pts`                      | Apply a 4×4 affine to a set of 3-D points (row-vector convention) |
| `affine_bounds`                         | Compute output bounding box after applying an affine transform |
| `affine_resample_3d`                    | Resample a 3D (or 4D) volume under an affine transform |
| `canonical_affine_transform`            | Set NIfTI affine to a diagonal voxel-scaling matrix (LPI origin) |
| `convert_2d_to_axisymmetric`            | Converts 2D k-Wave simulation grid, medium properties, and source to axisymmetric form by halving the shorter dimension (right half from center). |
| `convert_axisymmetric_to_2d`            | Mirrors axisymmetric simulation data (sensor, medium, masks, source) left-right to full 2D, doubles radial dimension, transposes axes, and updates kgrid/positions. |
| `convert_axisymmetric_to_3d`            | Expands 2D axisymmetric data (via radialExpand2DTo3D) to full 3D grid by duplicating radial dimension, updates positions/grid/kgrid for cubic symmetry. |
| `convert_final_to_MNI_matlab`           | Converts an image from subject space to MNI space using MATLAB. |
| `convert_final_to_MNI_simnibs`          | Converts an image to MNI space using SimNIBS. |
| `ensure_ras_plus`                       | Reorient a NIfTI volume to RAS+ if needed |
| `grid_mm_to_voxel`                      | Grid position in mm to a 1-based voxel index |
| `grid_to_ras`                           | Converts 1-based voxel indices to RAS coordinates (inverse of `ras_to_grid`). |
| `grid_to_ras_affine`                    | Build the acoustic-grid-voxel -> RAS-mm affine |
| `mni2subject_coords_LDfix`              | Transforms a set of coordinates in MNI space to subject space. |
| `nifti_affine_from_1based`              | Turns an affine on 1-based voxel indices into a NIfTI header affine (0-based indices). |
| `radialExpand2DTo3D`                    | Radially expands 2D axisymmetric data into 3D Cartesian volume. |
| `ras_to_grid`                           | Converts RAS coordinates to 1-based voxel (grid) indices using the NIfTI header affine. |
| `sample_volume_at_world`                | Nearest-neighbour values of a NIfTI volume at world (RAS/MNI mm) positions; NaN outside. |
| `subject2mni_coords_LDfix`              | Transforms a set of coordinates in MNI space to subject space. |
| `transform_coordinates`                 | Transform input coordinates between coordinate systems (wrapper) space. |

#### UNCERTAINTY

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `cleanup_uncertainty_intermediates`     | Delete large intermediate cache files after uncertainty analysis |
