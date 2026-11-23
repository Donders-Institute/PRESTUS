# PRESTUS Functions

This page lists every MATLAB function under `functions/` with a one-line description, so that you can find the function behind a stage, a setting or a message.

On this page:

- One table per folder of `functions/`, in folder order from [ACOUSTIC](#acoustic) to [UNCERTAINTY](#uncertainty): each function's name and what it does. Files at the top of `functions/` are listed under [CORE](#core).

#### ACOUSTIC

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `acoustic_analysis`                     | Compute acoustic intensity metrics and generate output plots and tables |
| `acoustic_convert_axisymmetry`          | Expand axisymmetric simulation fields to 3D or the 2D plane |
| `acoustic_fov_crop`                     | Resolve FOV settings and crop the acoustic grid if active. |
| `acoustic_grid_fov`                     | Crop the acoustic simulation domain to a beam-axis FOV box. |
| `acoustic_input_file`                   | Pressure map a run starts from instead of simulating ('' if none) |
| `acoustic_simulation`                   | Dispatch a k-Wave acoustic simulation across supported compute backends |
| `acoustic_simulation_elastic`           | Run a k-Wave elastic (shear-wave) simulation via pstdElastic2D/3D |
| `acoustic_stage`                        | Run or load the acoustic pressure simulation. |
| `acoustic_water_baseline`               | Decide and run the free-water baseline simulation. |
| `acoustic_wrapper`                      | Run the acoustic simulation and convert axisymmetric output to 2D or 3D |
| `air_cavity_metrics`                    | How much of the beam enters internal air cavities |
| `axisymmetric_output_mode`              | Form into which axisymmetric results are expanded |
| `cache_acoustic`                        | Save or load the slim acoustic results cache. |
| `coherent_drive_calibration`            | Set each transducer's drive to its own ISPPA target (coherent runs) |
| `compute_focal_fwhm`                    | Axial and lateral focal widths at half intensity (−3 dB) |
| `compute_focal_metrics`                 | −3 dB and −6 dB focal zone metrics per ITRUSST reporting |
| `fov_backproject_volume`                | Zero-pad a FOV sub-grid volume back to the full acoustic grid. |
| `fov_restore_transducers`               | Put every transducer back at its full-grid position after an FOV crop. |
| `fov_to_grid_pos`                       | Convert FOV-grid voxel indices to full sim-grid indices. |
| `is_acoustic_dry_run`                   | True when the configured acoustic backend only writes inputs |
| `kwave_pml_settings`                    | PML size and absorption for a k-Wave run (PML outside the grid) |
| `load_acoustic_input`                   | Start a run from a prior pressure map instead of simulating |
| `needs_water_baseline`                  | Decide whether a free-water baseline simulation must run |
| `plot_axial_profile`                    | Save axial pressure/intensity profile with tissue-layer shading |
| `predict_focal_shift`                   | Focus near the target predicted from the reverse field over an aperture |
| `prepare_thermal_domain`                | Restore full-grid state or adjust planimg for thermal FOV. |
| `remap_acoustic_field`                  | Resample an acoustic volume onto another grid, by quantity |
| `remap_complex_field`                   | Resample a complex pressure field onto another grid without losing its phase |
| `resolve_acoustic_fov_params`           | Resolve 'auto' FOV settings and check readiness. |
| `reverse_point_source_field`            | Complex CW field of a point source at the target, at the sensor points |
| `simultaneous_peaks`                    | Peak pressures of tones at different frequencies firing together |
| `steady_state_complex_amplitude`        | Complex pressure amplitude at given frequencies from a recorded window |
| `superpose_complex_fields`              | Coherent sum of complex pressure fields from separate runs. |
| `targeting_analysis`                    | Targeting accuracy of one run from its focal-mask NIfTIs |
| `time_reversal_pose_metrics`            | Target intensity and focus position of one pose-verification run |
| `time_reversal_pose_parameters`         | Run parameters that place the transducer at a pose of a pose search |
| `time_reversal_refocus`                 | Element transfers, refocusing phases and aberration metrics of one transducer |
| `time_reversal_stage`                   | Refocusing phases, aberration metrics and pose search from one reverse run per target |
| `write_acoustic_sidecar`                | Describe the T1w pressure NIfTI so a later run can start from it |
| `write_p_complex_files`                 | Save complex pressure fields to their own cache files |

#### BABELBRAIN

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `babelbrain_check_geometry`             | Error for geometries the BabelBrain backend cannot simulate |
| `babelbrain_grid_step_mm`               | BabelBrain's own spatial step for a run |
| `babelbrain_insitu_isppa`               | In-brain peak ISPPA of an acoustic field, as BabelBrain's thermal Isppa input |
| `babelbrain_match_grid`                 | Build the PRESTUS grid at BabelBrain's own step (babelbrain.match_grid) |
| `babelbrain_medium_table`               | Quantise the PRESTUS acoustic medium into BabelBrain materials |
| `babelbrain_placement_correction`       | Mechanical Tx offset that puts BabelBrain's apex at trans_pos |
| `babelbrain_preset`                     | Resolve a BabelBrain preset from BabelBrain's transducer registry |
| `babelbrain_skull_elastic`              | Skull elastic properties for BabelBrain, matched to kwave_elastic |
| `babelbrain_thermal_timing`             | Convert PRESTUS timing to BabelBrain thermal parameters |
| `babelbrain_thermal_wrapper`            | Run BabelBrain thermal simulation and read results |
| `babelbrain_tx_geometry`                | Transducer geometry BabelBrain simulates for a PRESTUS transducer |
| `babelbrain_wrapper`                    | Run BabelBrain acoustic simulation as PRESTUS backend |
| `export_babelbrain_inputs`              | Write all BabelBrain input files from a PRESTUS run |
| `export_babelbrain_weights`             | Write calibration result to BabelBrain HDF5 format |
| `import_babelbrain_weights`             | Import a BabelBrain calibration into PRESTUS |
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
| `calibrate_focal_steering`              | Analytical steer that lands the FLHM centre/peak at a target |
| `calibrate_library_entry`               | Calibrate a transducer library entry to hydrophone profiles at chosen depths |
| `calibration_desired_profiles`          | Profiles on a desired-depth grid, for prediction and fitting |
| `calibration_evaluate`                  | Geometric, library, fitted and k-Wave profiles on one set of depths |
| `calibration_fit`                       | Interactive calibration fit: per-depth refinements |
| `calibration_fit_depths`                | One FIT_TRANSDUCER_MODEL fit per depth of a profile set |
| `calibration_fit_entry`                 | Library entry from the current entry and a CALIBRATION_FIT result |
| `calibration_fit_options`               | FIT_TRANSDUCER_MODEL options for a per-depth calibration preset |
| `calibration_focus_beyond`              | FLHM centre of each axial profile, ignoring the near field |
| `calibration_pipeline_start`            | Universal launcher for calibration_transducer (MATLAB / SLURM / qsub) |
| `calibration_save_entry`                | Write a fitted library entry to a user library folder |
| `calibration_setup`                     | Run CALIBRATION_TRANSDUCER for every equipment combination in the calibration config |
| `calibration_transducer`                | Calibrate a transducer to hydrophone profiles at chosen depths and intensities |
| `check_equipment_config`                | Validate equipment configuration YAML files |
| `compute_flhm_center`                   | Focal length at half maximum (FLHM) centre of an axial profile. |
| `compute_phases`                        | Compute element steering phases for a target focal depth |
| `equipment_geometry_fingerprint`        | Short hash of the transducer geometry a calibration was fitted for |
| `evaluate_library_entry`                | Model the axial profiles a study gets from a library entry, at any depths |
| `extract_real_intensity_profile`        | Extract or interpolate an empirical axial intensity profile |
| `field_transfer_matrix`                 | Complex pressure per unit element drive at points off the axis |
| `fit_ring_profile`                      | Fit a within-ring velocity profile shared by all focal depths |
| `fit_steering_correction`               | Fit measured / modelled focal peak over a matrix array's steering range |
| `fit_transducer_model`                  | Fit a static annular-array model to hydrophone axial profiles |
| `focal_reference_to_frac`               | Map a focal-depth convention name to an intensity threshold fraction |
| `generate_calibration_report`           | Write an HTML summary of a transducer calibration |
| `interpolate_axial_profiles`            | Axial intensity profiles at desired focal depths from measured ones |
| `kwave_calibrate_drive`                 | Calibrate the drive of per-depth refinements in free-water k-Wave |
| `kwave_drive_metrics`                   | Absolute agreement of a simulated free-water profile with its target |
| `kwave_drive_submit`                    | Submit a k-Wave drive check script as a scheduler job |
| `library_add_conditions`                | Add dephased drive conditions to a library entry [EXPERIMENTAL] |
| `library_add_refinements`               | Merge depth-specific phase offsets (and drives) into a library entry |
| `library_condition`                     | Phases, gains and drive of a dephased library condition [EXPERIMENTAL] |
| `library_entry_file`                    | Path of a calibration-library entry in a flat or nested library |
| `library_entry_from_model`              | Write a FIT_TRANSDUCER_MODEL result into a library struct |
| `library_ingest_conditions`             | Add dephased conditions from a profile CSV + sidecar [EXPERIMENTAL] |
| `library_model_at_depth`                | Phases, gains and drive of a library entry at any focal depth |
| `library_provenance`                    | PRESTUS version string and commit hash for library entry meta |
| `library_refinement_from_fits`          | Write per-depth FIT_TRANSDUCER_MODEL results as library refinements |
| `library_set_ring_profile`              | Set or remove the within-ring velocity profile of a library entry |
| `library_settings_table`                | Per-depth transducer settings a library entry resolves to, as a table / CSV |
| `library_static_model`                  | Static model of a library entry, in FIT_TRANSDUCER_MODEL form |
| `library_transducer`                    | Transducer geometry a library entry was calibrated with |
| `load_equipment_config`                 | Load the equipment library (transducers and driving systems) |
| `load_transducer_from_library`          | Load calibrated transducer parameters from library |
| `matrix_element_layout`                 | Matrix-array elements in the transducer frame (no grid) |
| `matrix_focal_peak`                     | Modelled free-field focal peak of a matrix array per steering target |
| `matrix_rayleigh_pressure`              | Free-field complex pressure of a matrix array (Rayleigh integral) |
| `perform_global_search`                 | Global optimisation of transducer phases and particle velocity |
| `phase_optimization_annulus_full_curve` | Compute intensity profile error for full-curve phase optimisation |
| `plot_calibration_report`               | Measured vs library vs calibrated axial profiles, one panel per depth |
| `plot_ring_profile`                     | Show a fitted within-ring velocity profile and what it changes |
| `plot_transducer_library`               | Visualise the settings a calibration library entry resolves to |
| `rayleigh_axial_intensity`              | On-axis intensity via the Rayleigh–Sommerfeld integral |
| `read_axial_profiles`                   | Read a hydrophone axial-profile CSV (one column per focal depth) |
| `read_offaxis_scans`                    | Hydrophone plane and volume scans as off-axis data for FIT_RING_PROFILE |
| `read_steering_measurements`            | Read hydrophone focal peaks of a steered matrix array |
| `refine_transducer_library`             | Optional refinement of a library entry to your own hydrophone profiles |
| `resolve_profile_file`                  | Hydrophone-profile CSV named by an equipment entry's char_data_path |
| `resolve_transducer_from_serial`        | Populate transducer struct from equipment and calibration library |
| `ring_profile_drive`                    | Simulation drive of an annular transducer with a ring profile |
| `ring_profile_transducer`               | Annular transducer with a within-ring velocity profile |
| `save_elem_correction`                  | Compute and persist per-element phase correction into the equipment YAML |
| `save_library_entry`                    | Write a user-calibrated library entry with provenance |
| `save_optimized_values`                 | Save optimised transducer phases and amplitude to a CSV file |
| `scan_and_ingest_calibrations`          | Scan a folder and ingest all calibration YAMLs |
| `set_real_phases`                       | Compute geometric steering phases for a focal depth |
| `steering_correction_factor`            | Measured / modelled focal peak at steering targets |
| `steering_path_lengths`                 | Ring-centre to target distances used for geometric steering |
| `steering_quad_terms`                   | Design matrix of a full quadratic surface in three coordinates |
| `transducer_library_folders`            | Your own transducer libraries, most specific first |
| `update_transducer_library`             | Add depth-specific refinements to a library YAML |
| `validate_transducer_library`           | Validate calibration library against equipment configs |
| `water_baseline`                        | Measure free-water ISPPA at the current source amplitude |
| `water_baseline_box`                    | Grid, source and sensor of the free-water baseline cropped to the beam |
| `write_axial_profiles`                  | Write axial profiles as a profile CSV plus YAML sidecar |

#### CORE

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `acoustic_store_key`                    | Content key of one transducer's acoustic field in an async run |
| `apply_isppa_scaling`                   | Scale the simulated pressure field so its free-water ISPPA matches the target |
| `async_transducer_pipeline`             | Simulate N asynchronously firing transducers with independent acoustic fields and a joint thermal model |
| `combine_async_intensity`               | Scale and combine N independent acoustic results for async (non-simultaneous) transducer firing |
| `firing_groups`                         | Which transducers fire together, and the heat sources they form |
| `grid_stage`                            | Build or load the head model grid cache. |
| `is_async_mode`                         | Return true when async multi-transducer mode is active. |
| `is_multi_isppa_mode`                   | True when any transducer carries more than one ISPPA target. |
| `is_pose_workflow_mode`                 | True when time_reversal.pose_workflow.enabled is set |
| `is_simultaneous_multifreq`             | True for coherent transducers firing together at different frequencies |
| `isppa_baseline`                        | Free-water ISPPA an acoustic field was produced at, and its origin |
| `isppa_intent`                          | One-line statement of the free-water ISPPA a run will model |
| `isppa_scale_factor`                    | Pressure scale factor that brings a field to a target free-water ISPPA |
| `load_grid_cache_file`                  | Load a PRESTUS grid cache and restore pipeline state. |
| `load_parameters`                       | Load default PRESTUS config and merge with one or more extra YAML files or structs |
| `load_transducer_parameters`            | Validate, normalize, and derive transducer parameters from config |
| `medium_stage`                          | Build or load the acoustic/thermal medium property cache. |
| `multi_isppa_pipeline`                  | Run the acoustic simulation once, then thermal analysis at several free-water target intensities in parallel |
| `path_log_setup`                        | Add PRESTUS paths, create output directories, and start diary logging |
| `pose_workflow_pipeline`                | Transducer pose search, forward verification and final run in one call |
| `posthoc_water_parameters`              | Parameters of one post-hoc free-water re-entry run |
| `prestus_config_init`                   | Bootstrap a project-specific PRESTUS configuration folder |
| `prestus_config_upgrade`                | Refresh a project's copy of the defaults file (config_default.yaml) from the toolbox |
| `prestus_group_report_start`            | Manual entry point for the group HTML report (and group MNI maps) |
| `prestus_group_start`                   | Entry point for group-level analyses in MNI space |
| `prestus_pipeline`                      | End-to-end transcranial focused ultrasound simulation pipeline |
| `prestus_pipeline_start`                | Universal PRESTUS pipeline launcher for MATLAB, SLURM, and qsub |
| `prestus_study_start`                   | Run the PRESTUS pipeline for every row of a study manifest |
| `print_parameter_summary`               | Print a structured PRESTUS configuration summary to the console |
| `resolve_config_paths`                  | Make relative file and folder settings absolute |
| `sequential_pipeline`                   | Dispatch follow-up simulations that chain thermal state |
| `uncertainty_pipeline`                  | Run the full PRESTUS uncertainty quantification workflow |

#### GROUP

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `create_group_maps`                     | Group maps in MNI space from the per-run MNI NIfTIs, per condition |
| `discover_group_subjects`               | Auto-scan path.sim for completed PRESTUS subjects |
| `study_read_manifest`                   | Read a study manifest: one row per simulation run |
| `study_row_parameters`                  | Load the full config of one study run (one manifest row) |
| `study_run_descriptors`                 | Settings that tell one run of a study apart from the others |
| `study_run_files`                       | Output file names of one run (subject + medium + affix) |
| `study_run_table`                       | Resolve each manifest row to its output location (subject, affix, medium) |
| `study_setup`                           | Load a study config and resolve its manifest into runs |

#### GUI

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `calibration_gui`                       | Standalone GUI: fit transducer library entries to hydrophone profiles |
| `neuronav_gui`                          | Stand-alone GUI for Localite marker ingestion |
| `nexus_bridge`                          | Persistent MATLAB session that serves the PRESTUS Nexus |
| `nexus_calibration`                     | Transducer calibration for the PRESTUS Nexus (through NEXUS_BRIDGE) |
| `nexus_equipment`                       | Equipment-editor checks for the PRESTUS Nexus (through NEXUS_BRIDGE) |
| `nexus_neuronav`                        | Localite marker ingestion for the PRESTUS Nexus (through NEXUS_BRIDGE) |
| `nexus_pct`                             | Pseudo-CT steps of the PRESTUS Nexus, served by NEXUS_BRIDGE |
| `nexus_ping`                            | Answer a NEXUS_BRIDGE health check |
| `nexus_progress`                        | Report progress of the current NEXUS_BRIDGE request to the Nexus |
| `prestus_nexus`                         | Open the PRESTUS Nexus (standalone config builder) |

#### HEAD

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `air_cavity_charm`                      | Internal air from the SimNIBS charm atlas labelling, on the simulation grid |
| `air_cavity_mask`                       | Internal air cavities (sinuses, mastoid air cells) on the simulation grid |
| `autopad_grid_for_transducers`          | Expand grid if any transducer/focus position falls outside |
| `export_placement`                      | Write the placed transducer poses for other tools (placement.export_formats) |
| `head_model_load_or_build`              | Pose-independent head model of one subject, built once and cached |
| `head_model_place`                      | Simulation grid of one transducer placement from the shared head model |
| `head_smooth_and_crop`                  | Convert segmentation to medium masks and crop to simulation grid |
| `mean_beam_axis`                        | Mean axis of several transducer beams, taken as axes |
| `placement_run_dir`                     | Folder of one placement run (one mode, one target) |
| `placement_signature`                   | Settings a placement run's result depends on |
| `preproc_align_to_focal_axis`           | Rotate a 3D image so the focal axis aligns with the z-axis |
| `preproc_align_to_mean_focal_axis`      | Rotate a 3D image so the mean focal axis aligns with z |
| `preproc_crop_eCSF`                     | Expand CSF mask to exclude distant bone from the simulation grid |
| `preproc_crop_grid`                     | Crop 3-D head model to the simulation-relevant region |
| `preproc_deface_m2m`                    | Create a defaced copy of a SimNIBS m2m folder |
| `preproc_head`                          | Preprocess structural MRI for k-Wave simulation |
| `preproc_medium_mask`                   | Convert SimNIBS segmentation to smoothed medium label indices |
| `preproc_rescale_only`                  | Rescale a 3D image to grid resolution without rotation |
| `preproc_segmentation`                  | Segment structural MRI using SimNIBS |
| `preproc_transducer_placement`          | Resolve transducer and focus positions |
| `read_placement_record`                 | Read the run record of a placement run folder |
| `seg_qc`                                | Record head-model postprocessing steps for the segmentation report |
| `segmentation_check_existing`           | How an existing segmentation differs from the selected SimNIBS |
| `segmentation_run`                      | Submit a SimNIBS charm segmentation job |
| `skull_bone_source`                     | Where this run's skull layer comes from |
| `skull_fill_holes`                      | Expand skull mask and fill gaps to maintain continuity |
| `skull_from_hu`                         | Define the skull layer by a Hounsfield threshold on the simulation grid |
| `skull_reconcile`                       | Make a bone mask the skull layer, reconciled with the segmentation |
| `skull_rubber_wrap`                     | Fill skull holes where brain meets skin, as a rubber sheet around the skull would |
| `skull_rubber_wrap_visualize`           | 3-D visualisation of skull balloon expansion |
| `smooth_img`                            | Smooth an image with a Gaussian or box kernel defined by FWHM |
| `tissue_from_hu`                        | Soft tissue against water by a CT threshold on the simulation grid |
| `write_placement_positions`             | Write the resolved transducer positions as a JSON sidecar |
| `write_placement_record`                | Write the run record (placement_run.json) of a placement run |

#### HELPER

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `apply_config_overrides`                | Set config values by dotted key, e.g. from a study manifest row |
| `bb_tissuemask_binary`                  | Build per-tissue binary masks from a BabelBrain material map |
| `bwconncomp_3d`                         | Find connected components in a 3-D binary volume |
| `cache_signature`                       | Settings that determine a cached pipeline stage's output |
| `cast_struct`                           | Recursively cast all numeric fields in a struct to a given type |
| `charm_seg_labels`                      | Return the fixed SimNIBS charm tissue label indices |
| `check_availability`                    | Verify that a list of files or glob patterns exist on disk |
| `check_cache_signature`                 | Refuse to reuse a cache built from different settings |
| `check_layers`                          | Validate requested medium layers against available segmentation values |
| `check_sidecar_signature`               | Compare a sidecar's recorded settings with the current run's |
| `check_tool_availability`               | Check that the external tools the requested stages need exist |
| `compute_ramp_factor`                   | Pulse-average / temporal-peak intensity ratio for a ramped pulse |
| `conda_env_dirs`                        | Folders of the conda environments on this computer |
| `config_known_extra_keys`               | Config keys PRESTUS reads that config_default.yaml omits |
| `config_unknown_keys`                   | Dotted keys of a user config that PRESTUS does not know |
| `confirm_overwriting`                   | Decide whether an existing file should be overwritten |
| `confirmation_dlg`                      | Display a two-button modal dialog and return the user's choice |
| `debug_on`                              | Whether a stage writes a kind of debug output |
| `debug_register`                        | Record a debug output in the run's debug manifest |
| `debug_save_figure`                     | Save a debug figure, register it, and close it |
| `debug_stages`                          | Stages selected by simulation.debug_stages |
| `default_log_file`                      | Timestamped pipeline log (diary) path for a run |
| `distribution_stats`                    | Summary statistics and histogram of a numeric vector |
| `ensure_hdf5_plugin_path`               | Set HDF5_PLUGIN_PATH so MATLAB's h5read/h5write don't fail |
| `find_min_factor`                       | Find the integer in a range whose largest prime factor is smallest |
| `find_nexus_python`                     | A Python that can start the PRESTUS Nexus |
| `find_pseudoct_dir`                     | Folder holding this subject's existing pseudoCT.nii.gz |
| `find_sidecar`                          | JSON sidecar of a data file: same name with .json (BIDS convention) |
| `get_crop_dims`                         | Compute bounding-box crop dimensions for a 3D binary/label image |
| `get_field`                             | Value of a (nested) struct field, or a default when it is missing |
| `get_flhm_center_position`              | Compute the centre of the focal length at half maximum (FLHM) |
| `get_nonempty_field`                    | Value of a (nested) struct field, or a default when it is missing or empty |
| `get_output_dir`                        | Derive the simulation output directory from parameters |
| `get_pct_dir`                           | Per-subject pseudo-CT folder |
| `get_placement_dir`                     | Per-subject root of the placement outputs |
| `get_prestus_path`                      | Return the absolute path to the PRESTUS toolbox root |
| `get_slice_by_label`                    | Extract a 2D slice from a 3D volume along a named axis |
| `get_xyz_mesh`                          | Generate an [Nx3] coordinate mesh for all voxels in a 3D image |
| `getidx`                                | Return tissue label indices matching a name string from a layers struct |
| `iff`                                   | Inline conditional: return a if cond is true, else b. |
| `imfill_holes_3d`                       | Fill interior holes in a 3-D binary volume |
| `kwave_version`                         | Print k-Wave toolbox version, git commit, and CUDA binary version |
| `live_view`                             | Save the evolving fields of a simulation as frames for the Nexus |
| `log_memory_usage`                      | Print this MATLAB process's RSS and the enclosing cgroup's usage. |
| `log_timer`                             | Track elapsed time, peak RAM, and drive space for labelled code blocks |
| `manual_coordinate_space`               | Coordinate space of transducer.trans_pos / focus_pos |
| `masked_max_3d`                         | Find the maximum value and its 3D subscript within a masked region |
| `matlab_escape_quotes`                  | Escape text for a single-quoted MATLAB char literal |
| `mergestruct`                           | Recursively merge the fields of structure B into structure A |
| `morph_close_sphere`                    | Morphological closing of a 3-D mask by a sphere of radius r |
| `nearest_fill_nan`                      | Replace NaN voxels with the value of the nearest non-NaN voxel |
| `outputs_follow_recompute`              | Let a recomputed stage replace its existing output files |
| `parse_memory_tags`                     | Extract peak memory figures from [memory]/[memory_bb] log lines. |
| `predict_outputs`                       | Predict which files a PRESTUS run will write, before submission. |
| `prestus_check_install`                 | Report what a PRESTUS installation has and what it lacks |
| `prestus_progress`                      | Print a progress line that the Nexus reads from the run log |
| `prestus_version`                       | Print PRESTUS version info and optionally return the git hash. |
| `print_predicted_outputs`               | Print PREDICT_OUTPUTS' plan to the console. |
| `r_quote_string`                        | Escape a string for safe embedding in a generated R string literal. |
| `read_ini_file`                         | Parse an INI file into a MATLAB struct |
| `read_sidecar`                          | Read a stage-output JSON sidecar |
| `resolve_chained_overwrite`             | Map io.overwrite_files onto the stages of a chained pipeline (async, multi-ISPPA, uncertainty) |
| `resolve_io_dirs`                       | Fill missing io path fields from path.sim + subject_id |
| `resolve_jwave_pml_size`                | Resolve grid.pml_size 'auto' to a wavelength-based scalar PML thickness for jwave |
| `resolve_pml_size`                      | Resolve grid.pml_size to a numeric array for k-Wave |
| `resolve_python_env`                    | Locate a Python executable inside a venv/conda environment. |
| `reuse_upstream_cache`                  | True when an existing upstream cache must be loaded regardless of io.overwrite_files |
| `round_if_integer`                      | Round a value to the nearest integer or error if not close enough |
| `run_command`                           | Execute an external program via java.lang.ProcessBuilder. |
| `run_shell_script`                      | Execute a shell script string via ProcessBuilder, with only ONE shell in the loop. |
| `run_stats`                             | Collect a run's simulation size, duration and peak memory |
| `safe_addpath`                          | Add a directory tree to the MATLAB path, skipping hidden directories |
| `sample_beam_axis`                      | Sample a field along the beam axis from the transducer bowl to the grid edge |
| `sanitize_output_affix`                 | Normalise io.output_affix to the BIDS desc entity |
| `segmentation_file`                     | The SimNIBS segmentation of this subject (final_tissues.nii.gz) |
| `segmentation_fingerprint`              | Content identity of the SimNIBS segmentation a stage reads |
| `shell_quote_arg`                       | Quote a string as a single literal argument for a shell command line |
| `should_save_output`                    | Decide whether a particular pipeline output should be saved |
| `simnibs_env_overrides`                 | Environment for a SimNIBS tool run with RUN_COMMAND |
| `simnibs_install_info`                  | Version and default charm settings of a SimNIBS install |
| `simnibs_ld_library_path`               | LD_LIBRARY_PATH prefix for SimNIBS tools (charm, mni2subject, ...) |
| `simnibs_version`                       | Extract SimNIBS version from charm_log.html and store in parameters |
| `struct_md5`                            | MD5 of a value's content, independent of struct field order |
| `subset_fields`                         | Copy a named subset of fields from one struct to a new struct |
| `tissuemask_binary`                     | Build per-tissue binary masks from a medium label map |
| `trans_suffix`                          | Return '_trans-XX' when n_trans > 1, otherwise ''. |
| `upsample_to_grid`                      | Resample a 3D image and update grid positions to a higher-resolution grid |
| `validate_any_required_field`           | Error clearly unless at least one candidate field is set. |
| `validate_required_fields`              | Error clearly unless every required field is set. |
| `warn`                                  | Issue a warning without the stack-trace printout |
| `write_sidecar`                         | Write a JSON sidecar describing a stage output later runs can start from |
| `zip_fields`                            | Convert a struct into an alternating {field, value, field, value, ...} cell array |

#### HPC

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `hpc_add_dependency`                    | Add scheduler job IDs to the jobs a submission waits for |
| `hpc_detect_system`                     | Detect the available HPC job scheduler |
| `hpc_job_active`                        | True while a scheduler job is queued or running |
| `hpc_job_info`                          | Build a formatted job-info struct and optionally print it |
| `hpc_job_name`                          | Generate a standardised HPC job name for a PRESTUS pipeline run |
| `hpc_job_report`                        | Query scheduler for post-job resource usage and carbon metrics |
| `hpc_scratch_finalize`                  | Remove a job's node-local scratch folder |
| `hpc_scratch_setup`                     | Put the k-Wave .h5 files of a SLURM job on node-local scratch |
| `hpc_setup_temp_files`                  | Create output/log directories and timestamped temp files |
| `hpc_submit_job`                        | Write a scheduler script and submit a MATLAB job to SLURM or qsub |
| `hpc_submit_report_job`                 | Submit a lightweight dependent job to generate the post-run report |
| `hpc_validate_parameters`               | Validate and enforce settings required for HPC batch submission |
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
| `legacy_bridge_entry`                   | Accept pipeline calls and parameter structs written for earlier PRESTUS versions |
| `legacy_bridge_layer`                   | Drop old config keys that the same layer makes redundant |
| `legacy_bridge_struct`                  | Convert old key names set on a parameters struct after LOAD_PARAMETERS |
| `legacy_hpc_submit`                     | Shared body of the v0.5.0 HPC entry points (SINGLE_SUBJECT_PIPELINE_WITH_SLURM / _QSUB) |
| `single_subject_pipeline`               | Entry point of PRESTUS v0.5.0 and earlier; runs PRESTUS_PIPELINE_START in this MATLAB |
| `single_subject_pipeline_with_qsub`     | qsub entry point of PRESTUS v0.5.0 and earlier; submits via PRESTUS_PIPELINE_START |
| `single_subject_pipeline_with_slurm`    | SLURM entry point of PRESTUS v0.5.0 and earlier; submits via PRESTUS_PIPELINE_START |
| `upgrade_config_keys`                   | Rename old config keys (PRESTUS v0.5.0 onwards) in one config layer to their current names |

#### MEDIUM

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `db_cm_mhz_to_np_m`                     | Convert an attenuation or absorption coefficient to nepers per metre |
| `fitPowerLawParamsMulti`                | Fit power-law absorption prefactors for the k-Wave fractional Laplacian solver |
| `get_alpha_coeff`                       | Compute amplitude attenuation coefficient for a given medium and frequency |
| `medium_absorption_layer`               | Keep one absorption field per tissue when merging a config layer |
| `medium_elastic_properties`             | Build the k-Wave elastic medium (compressional + shear) from a fluid medium |
| `medium_pct_attenuation`                | Map pseudo-CT Hounsfield values to skull attenuation coefficients |
| `medium_pct_density`                    | Map pseudo-CT Hounsfield values to skull bone density |
| `medium_pct_density_calibration`        | CT->density calibration table for pct.mapping_density 'k-plan' |
| `medium_pct_hu_max`                     | Densest-bone HU for the max-normalised pCT mappings |
| `medium_pct_pichardo_load`              | Load and cache Pichardo 2D interpolant from BabelBrain |
| `medium_pct_porosity`                   | Skull porosity from HU for the porosity mixture mappings |
| `medium_pct_soundspeed`                 | Map pseudo-CT Hounsfield values to skull bone sound speed |
| `medium_setup`                          | Build the k-Wave medium property struct from tissue masks and config |
| `skull_property_stats`                  | Modelled skull acoustic properties, in all skull voxels and on the beam path |

#### NEURONAV

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `export_brainsight_trajectory`          | Write a Brainsight version-14 trajectory file. |
| `export_kplan_plan`                     | Write a transducer pose and its images for k-Plan |
| `export_slicer_transform`               | Write the transducer pose as a 3D Slicer ITK transform (.tfm) |
| `grid_to_localite_mm`                   | Voxel indices to coordinates in the Localite planning image |
| `localite_matrix_to_positions`          | Convert a Localite 4×4 matrix to voxel positions. |
| `neuronav_apply_deformation`            | Map segmentation-space RAS coordinates through a SimNIBS coordinate field |
| `neuronav_build_deformation_cache`      | Precompute nearest-valid-voxel lookup for a SimNIBS coordinate field |
| `neuronav_build_target_map`             | Expand a coil_map spec into a flat target_map. |
| `neuronav_check_target_map`             | Check a target_map before NEURONAV_INGEST_MARKERS writes it |
| `neuronav_compute_series_statistics`    | Compute per-series position and orientation statistics |
| `neuronav_compute_series_statistics_sequential`| Per-series position and orientation statistics with a settable series gap |
| `neuronav_convert_MNI_to_native`        | Transform coordinates from MNI space to native subject space |
| `neuronav_convert_native_to_MNI`        | Convert native RAS coordinates to MNI space |
| `neuronav_convert_trigger_to_voxels`    | Convert Localite trigger positions to voxel coordinates |
| `neuronav_detect_series`                | Detect the recorded series of one session without writing. |
| `neuronav_export_session_csv`           | Export per-session coordinates to a labelled CSV file |
| `neuronav_find_mni_template`            | SimNIBS' MNI152_T1_1mm.nii.gz, for head outlines in MNI views |
| `neuronav_find_planning_t1`             | Locate the neuronavigation planning T1 of a subject. |
| `neuronav_get_group_mean_mni`           | Compute group-level average MNI coordinates |
| `neuronav_head_surface`                 | Outer head surface of an MNI-space NIfTI, for plots |
| `neuronav_ingest_batch`                 | Ingest and export Localite targets for many subjects/sessions. |
| `neuronav_ingest_markers`               | Ingest raw Localite markers and write normalised outputs. |
| `neuronav_label_map_from_scan`          | Build a neuronav.target_labels map from a label scan. |
| `neuronav_localite_file_series`         | The marker series of one Localite XML export. |
| `neuronav_localite_pose`                | Resolve the Localite 4×4 instrument pose for placement. |
| `neuronav_map_label`                    | Map a Localite target label to the study's target name. |
| `neuronav_match_entry_targets`          | Match Localite EntryTarget entries to a target map. |
| `neuronav_name_series`                  | Name the detected series of a session (shared GUI/script rule). |
| `neuronav_plot_placement`               | Visualise transducer placements on a segmented head. |
| `neuronav_read_entry_target`            | Read planned target positions from Localite EntryTarget XML. |
| `neuronav_read_label_candidates`        | Collect named planned targets for a session. |
| `neuronav_read_target_lookup`           | Look up pre-specified target coordinates from a study lookup table, keyed by subject and target name |
| `neuronav_resolve_paths`                | Resolve the Localite folders of a PRESTUS config. |
| `neuronav_resolve_target`               | Resolve each target_map entry's brain-target coordinate from the best available source |
| `neuronav_scan_sessions`                | Inventory raw Localite data and planned-target labels. |
| `neuronav_select_localite`              | Read the normalised InstrumentMarker XML for a session. |
| `neuronav_series_qc_flag`               | Tracking problem of one detected series, as text |
| `neuronav_suggest_target_names`         | Heuristic target names for ingested series. |
| `neuronav_summary_to_mni`               | Add MNI coordinates to a neuronav ingest summary. |
| `position_transducer_localite`          | Determine transducer/focus positions from a Localite XML file. |
| `position_transducer_plantus`           | Optimise transducer placement using PlanTUS. |
| `position_transducer_tus_entry`         | Determine transducer placement using TUS_entry. |

#### NIFTI

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `get_t1_mni_file`                       | Path of a subject's T1-in-MNI image, outside the SimNIBS m2m folder |
| `nifti_acoustic`                        | Export acoustic-stage data to NIfTI in T1 and MNI space |
| `nifti_medium`                          | Export medium-stage data to NIfTI |
| `nifti_medium_masks`                    | Write medium_masks pipeline output to NIfTI |
| `nifti_medium_properties`               | Write medium acoustic/thermal property maps to NIfTI |
| `nifti_pct`                             | Write pseudo-CT pipeline output to NIfTI |
| `nifti_thermal`                         | Export thermal-stage data to NIfTI in T1 and MNI space |
| `nifti_to_mni`                          | Convert a T1-space NIfTI to MNI space via SimNIBS |
| `nifti_to_t1w`                          | Back-transform one volume to T1 space and write NIfTI |
| `nifti_transducer`                      | Export a model of the planned transducer in T1 space |
| `phantom_voxel_order`                   | Return a non-layered output volume to the voxel order of its input |
| `simnibs_t1_to_mni`                     | Register subject T1 to MNI space if not already done |
| `write_transducer_nifti`                | Write a transducer model on the exact voxel grid of a T1 |

#### PSEUDO-CT

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `check_pct_sidecar`                     | Refuse an existing pseudo-CT made with other settings |
| `fit_pairwiselinear`                    | Piecewise linear mapping from Hounsfield Units to mass density |
| `pct_apply_limits`                      | Regularise CT HU or a CT-derived skull property to configured limits |
| `pct_bone_mask_file`                    | This subject's bone mask file (pct.bone_source 'mask') |
| `pct_create_pCT_run`                    | Submit a pseudoCT generation job |
| `pct_ct_file`                           | This subject's measured CT, or '' when the pseudo-CT is used |
| `pct_hu_reference`                      | HU that max-normalised skull mappings take as the densest bone |
| `pct_load_bone_mask`                    | The bone mask the skull is built from, on the segmentation grid |
| `pct_load_ct`                           | The CT image the skull mapping reads, with its tissue labels |
| `pct_map_properties`                    | Skull density, sound speed and attenuation from HU |
| `pct_mapping_curves`                    | Skull properties as functions of HU, for the configured mappings |
| `pct_predict_skull_properties`          | Skull properties a pseudo-CT run will model, without simulating |
| `pct_signature`                         | Settings a pseudo-CT depends on, and the resolved skull mapping |
| `pct_skull_mapped`                      | True when the skull layer's properties come from the pseudo-CT |
| `pct_skullexpand`                       | Load SimNIBS NIfTIs, run skull rubber wrap, and save output |
| `pct_skullmapping`                      | Compute pseudo-CT mapping for cortical and trabecular bone from UTE histograms |
| `pct_skullmapping_fit`                  | Subject-specific UTE->HU line from cortical / trabecular UTE histograms |
| `pct_soft_tissue_peak`                  | Identify the soft tissue peak from a UTE image intensity distribution |
| `pct_tissues_mask`                      | Tissue labels for a CT-mapped skull, from the SimNIBS segmentation |
| `pct_uses_hu_reference`                 | Whether a pCT mapping of this run scales by the densest-bone HU |

#### PLOT

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `deface_volume`                         | Remove face voxels in a volume using the SimNIBS segmentation |
| `nifti_slice_orientation`               | Compute display transforms for a NIfTI orthogonal slice. |
| `plot_intensity_t1_overlay`             | Overlay back-transformed intensity on T1 image |
| `plot_overlay`                          | Overlay a metric map on a 2D slice of a 3D background image |
| `plot_overlay_2d`                       | Overlay a metric map on a 2D background image |
| `plot_placement_t1_overlay`             | Quick placement QC figure on the native T1 |
| `plot_t1_with_transducer`               | Plot a T1 slice with transducer mask overlaid |
| `plot_transducer_overlay`               | Draw transducer geometry with labelled rectangles on the current axes |
| `show_3d_head`                          | Visualise a segmented 3D head with transducer and target positions |
| `show_binary_mask_transducer`           | Visualise the binary mask of a matrix transducer on the simulation grid |
| `show_positioning_plots`                | Visualise transducer positioning before and after preprocessing |

#### REPORT

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `build_config_dump`                     | Searchable HTML block of the full parameters struct. |
| `compute_tic`                           | Cranial Thermal Index (TIC) after IEC 62359 / AIUM-NEMA ODS. |
| `css_styles_base`                       | Shared design system for all PRESTUS HTML reports. |
| `generate_debug_report`                 | One HTML page indexing all debug outputs of a run |
| `generate_group_report`                 | Generate a self-contained group-level HTML report |
| `generate_multi_isppa_report`           | Generate an HTML report comparing multiple ISPPA levels |
| `generate_placement_report`             | Stand-alone HTML report of a transducer placement |
| `generate_segmentation_report`          | HTML report of the head-model postprocessing |
| `generate_sequential_report`            | HTML summary report for a chain of sequential simulations |
| `generate_simulation_report`            | Generate a self-contained HTML simulation report |
| `generate_uncertainty_report`           | Generate an uncertainty-range HTML simulation report |
| `get_risk_limits`                       | Return ITRUSST consensus non-significant risk limits |
| `html_utils`                            | Static HTML utility methods for PRESTUS report generation |
| `preview_report_styles`                 | Re-render an existing simulation report in several |
| `read_metrics_table`                    | Read a per-run metrics CSV with current column names |
| `report_hidden_columns`                 | Metrics-CSV columns that the HTML reports never show |
| `report_link`                           | href from a report to a file or folder |
| `report_scripts`                        | Shared inline JavaScript for PRESTUS HTML reports. |
| `report_styles`                         | Alternative visual skins for the PRESTUS HTML reports |
| `risk_color`                            | Return a colour string based on a value relative to an ITRUSST risk limit |
| `table2html`                            | Convert a MATLAB table to an HTML table element |

#### SOURCE

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `create_clover_array`                   | Replicate a single matrix sub-array into a multi-leaf clover configuration |
| `create_matrix_karray`                  | Add matrix transducer elements to a kWaveArray object |
| `cw_ramp_periods`                       | Start-ramp length for CREATECWSIGNALS that fits the simulated time |
| `CWSourceMatrix`                        | Lazy rank-2 factored CW source signal matrix for k-Wave. |
| `grid_axisymmetry`                      | Adapt grid and tissue masks for k-Wave axisymmetric simulation |
| `grid_tissue_setup`                     | Set up computational grid and tissue masks for TUS simulation |
| `grid_transducer_location`              | Map transducer and focus positions into the simulation grid |
| `source_create`                         | Build a k-Wave source struct from PRESTUS transducer parameters |
| `source_sensor_setup`                   | Create kWaveGrid, CW source, and full-grid sensor for a TUS simulation |
| `water_grid_dims`                       | Free-water grid size from the transducer aperture and focal depth |

#### TELEMETRY

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `generate_run_id`                       | A random UUID v4 string for telemetry (run ID and install ID). |
| `prefdir_prestus`                       | Return (and create) the PRESTUS per-user preferences directory. |
| `telemetry_env_override`                | Read the PRESTUS_TELEMETRY environment variable. |
| `telemetry_has_display`                 | True when MATLAB can render a figure/dialog window. |
| `telemetry_is_headless_cluster`         | True for compiled apps and SLURM/PBS jobs. |
| `telemetry_is_noninteractive`           | True when no telemetry prompt may be shown. |
| `telemetry_level_info`                  | The usage-statistics levels a user can share, in prompt order. |
| `telemetry_levels`                      | The usage-statistics levels this process may send. |
| `telemetry_notice_version`              | Version of the telemetry notice that consent is given to |
| `telemetry_run_type`                    | Resolve whether this run is 'technical', 'development' or 'deployment'. |
| `telemetry_set_consent`                 | Record which usage statistics this user shares. |
| `telemetry_setup`                       | Show the opt-in prompt for anonymous usage statistics. |
| `telemetry_setup_reset`                 | Remove the telemetry consent file so the opt-in prompt is shown again on the next run |
| `telemetry_write_consent`               | Write the consent record ~/.prestus/telemetry.json. |
| `track_usage`                           | Send anonymous usage statistics to the PRESTUS telemetry endpoint. |

#### THERMAL

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `cem43_iso_increment`                   | ISO thermal-dose (CEM43) increment for one time step |
| `thermal_absorption_model`              | Resolve how acoustic heat deposition is specified for a tissue |
| `thermal_absorption_properties`         | Absorption model, coefficient and fraction of one tissue |
| `thermal_analysis`                      | Summarise and visualise results of a thermal ultrasound simulation |
| `thermal_grid_setup`                    | Build an independent k-Wave grid for thermal simulation |
| `thermal_overlap_terms`                 | Interference of heat sources whose pulses partly overlap |
| `thermal_parameters`                    | Compute time-stepping parameters for a multi-level thermal sonication protocol |
| `thermal_plane_frames`                  | Timepoints whose focal-plane frames the thermal simulation keeps. |
| `thermal_plot_protocol`                 | Plot the sonication protocol of a thermal simulation |
| `thermal_plot_sim`                      | Visualise heating simulation results over time and save figures |
| `thermal_ramp_factors`                  | Heat factor of a ramped pulse for each heat source |
| `thermal_schedule`                      | Firing schedule of one or more heat sources as thermal time segments |
| `thermal_simulation`                    | Simulate thermal effects of transcranial ultrasound using k-Wave |
| `thermal_stage`                         | Run or load the bio-heat thermal simulation. |
| `thermal_update_timeseries`             | Accumulate per-tissue-layer thermal metric timeseries |

#### TRANSDUCER

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `annular_aperture_points`               | Sample points on the radiating surface of a curved annular array |
| `annular_bowl_depth_mm`                 | Distance from the transducer position (bowl apex) to the exit plane |
| `annular_geometric_phases`              | Ring phases that steer a curved annular array to a point, in water |
| `convert_to_element_pos`                | Generate element positions for matrix transducers |
| `exit_plane_depth_mm`                   | Distance from the bowl apex to the exit plane used for placement |
| `extract_element_pos`                   | Read matrix transducer element positions from file |
| `fill_annular_drive`                    | Steering phases and drive of an annular transducer from its depth setting |
| `focal_distance_calculation`            | Compute expected focal distances for all transducers |
| `get_arc`                               | Generate x/y coordinates of a 2-D arc |
| `get_trans_pos_from_trigger_markers`    | Compute transducer/target positions from a Localite XML file |
| `get_transducer_box`                    | Compute transducer bounding box, exit plane, and geometric focus |
| `is_flat_matrix`                        | True for a flat matrix array (is_curved false or infinite radius) |
| `matrix_elem_diameter`                  | Aperture diameter of a circular (disc/bowl) matrix element |
| `matrix_element_aim`                    | Point each matrix element faces (its normal), curved or flat |
| `reciprocity_pose_candidates`           | Transducer poses around the planned one, for scoring on one reverse field |
| `reciprocity_score_poses`               | Score transducer poses on one reverse field (acoustic reciprocity) |
| `tp_candidate_mesh`                     | Build coordinate mesh and candidate transducer geometry |
| `tp_evaluate_candidate_positions`       | Evaluate candidate transducer positions on skull surface |
| `tp_find_initial_candidate`             | Find an initial candidate transducer position on the skull surface |
| `tp_plot_candidate_positions`           | Plot transducer candidate positions on a skull surface slice |
| `tp_plot_geometry_overlay`              | Visualise transducer geometry on a skin segmentation slice |
| `tp_plot_heuristic_position`            | Create visualisation of heuristic transducer placement |
| `tp_remove_ear_locations`               | Exclude candidate positions near ear centres |
| `tp_scalp_boundary`                     | Outer (scalp) boundary voxels of a segmented head volume |
| `tp_select_heuristic_position`          | Select optimal transducer position and export Localite coordinates |
| `transducer_analyze_position`           | Compute geometric and statistical measures for a transducer position |
| `transducer_analyze_position_fast`      | GPU-accelerated transducer position analysis |
| `transducer_positioning`                | Heuristic transducer placement for MNI targets |
| `transducer_positioning_start`          | Dispatch transducer positioning to MATLAB or HPC |
| `transducer_scalp_geometry`             | Place a transducer on the scalp with a standoff gap |
| `transducer_setup`                      | Create a transducer mask and label matrix for a k-Wave grid |
| `validate_annular_transducer`           | Validate configuration fields of an annular transducer |
| `validate_matrix_transducer`            | Validate and normalise configuration of a matrix transducer |

#### TRANSFORM

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `affine_apply_pts`                      | Apply a 4×4 affine to a set of 3-D points (row-vector convention) |
| `affine_bounds`                         | Compute output bounding box after applying an affine transform |
| `affine_inverse`                        | Inverse of a 4x4 affine, in the convention it was given |
| `affine_resample_3d`                    | Resample a 3D (or 4D) volume under an affine transform |
| `canonical_affine_transform`            | Set NIfTI affine to a diagonal voxel-scaling matrix (LPI origin) |
| `convert_2d_to_axisymmetric`            | Extract the right half-plane of a 2D medium for axisymmetric simulation |
| `convert_axisymmetric_to_2d`            | Mirror axisymmetric simulation fields into a full 2D cross-section |
| `convert_axisymmetric_to_3d`            | Expand all axisymmetric simulation fields from 2D to 3D |
| `convert_final_to_MNI_simnibs`          | Transform one or more images to MNI space |
| `ensure_ras_plus`                       | Reorient a NIfTI volume to RAS+ if needed |
| `grid_mm_to_voxel`                      | Grid position in mm to a 1-based voxel index |
| `grid_remap_affine`                     | Voxel mapping from one simulation grid (or FOV box) to another |
| `grid_to_ras`                           | Convert 1-based voxel indices to RAS world coordinates |
| `grid_to_ras_affine`                    | Build the acoustic-grid-voxel -> RAS-mm affine |
| `grid_world_affine`                     | Simulation-grid voxel -> world (RAS mm) affine of a run |
| `mni2subject_coords_LDfix`              | Transform MNI coordinates to subject space via SimNIBS CLI (LD_LIBRARY_PATH hotfix) |
| `nifti_affine_from_1based`              | NIfTI header affine from a 1-based-index affine |
| `radialExpand2DTo3D`                    | Radially expand a 2D axisymmetric half-plane into a 3D Cartesian volume |
| `ras_to_grid`                           | Convert RAS world coordinates to 1-based voxel indices |
| `ras_to_voxel`                          | RAS world coordinates to 1-based voxel indices, fractional or rounded |
| `sample_volume_at_world`                | Nearest-neighbour values of a NIfTI volume at world (mm) positions |
| `subject2mni_coords_LDfix`              | Transform subject RAS+ coordinates to MNI space via SimNIBS CLI (LD_LIBRARY_PATH hotfix) |
| `transform_coordinates`                 | Convert a 3D position between MNI, RAS+, and simulation grid spaces |

#### UNCERTAINTY

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `cleanup_uncertainty_intermediates`     | Delete preprocessing cache files after uncertainty analysis |
| `uncertainty_fov_baseline`              | Baseline temperature of an absolute-temperature metric |
| `uncertainty_fov_metric_kind`           | How a metric is rescaled from the partial FOV to the whole grid |
| `uncertainty_fov_mismatch`              | Whole-grid vs. partial-FOV factors of the default medium |
| `uncertainty_fov_rescale`               | Map a partial-FOV metrics table onto the whole grid |
| `uncertainty_medium_layer`              | Medium-property layer of one uncertainty variant |
| `uncertainty_variants`                  | The simulation variants of the uncertainty pipeline |
