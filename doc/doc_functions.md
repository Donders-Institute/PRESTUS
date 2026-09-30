## PRESTUS functions

The following documents the functions provided in PRESTUS.

#### ACOUSTIC

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `acoustic_analysis`                     | Compute acoustic safety metrics (Isppa, MI, max pressure) from sensor data, extract tissue-specific maxima, focal distances, and generate overlaid ISPPA plots on segmentation. |
| `targeting_analysis`                    | Targeting accuracy of one run from its focal-mask NIfTIs: simulated focus vs planned focus (axial/lateral), vs MNI target and ROI; writes `sub-XXX_<medium><affix>_targeting.csv`. |
| `acoustic_convert_axisymmetry`          |  Expand k-Wave axisymmetric (2D rotational) simulation outputs to full 3D (via  convert_axisymmetric_to_3d ) for heating sims or 2D Cartesian (via  convert_axisymmetric_to_2d ) otherwise, updating sensor data, parameters, masks, medium, grid, source, and labels. |
| `acoustic_simulation`                   | Set up acoustic simulation                                                                      |
| `acoustic_wrapper`                      | Set up acoustic simulation                                                                      |
| `compute_focal_metrics`                 | Compute focal-zone metrics (ITRUSST −3 dB/−6 dB volumes, FLHM span/centre) from a 3D intensity field; FLHM centre threshold is configurable via `focal_reference` (see [doc_calibration.md](doc_calibration.md)). |

#### CALIBRATION

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `calc_virtual_elem`                     | (STANDALONE) Calculate additional virtual elements                                              |
| `calibration_transducer`                | Run water simulations to optimize transducer source amplitude and element phases matching a target intensity profile via O’Neil analytical solution, global search, and re-simulation. |
| `compute_flhm_center`                   | Compute the FLHM (Focal Length at Half Maximum) centre/span of an axial intensity profile at a configurable threshold fraction of peak (0.5 = −3 dB default, 0.25 = −6 dB). |
| `compute_oneil_solution`                | Compute analytical O’Neil pressure along the beam axis for focused annular transducers in water, derives particle velocity and grid adjustment factor, and plots comparisons with simulated/desired intensities. |
| `compute_phases`                        | Calculate per-element transducer phases (degrees) for electronic focusing by computing fractional wavelength delays to steer waves constructively to a target point relative to the exit plane. |
| `babelbrain_insitu_isppa`               | In-brain peak ISPPA of the scaled acoustic field — the `Isppa` BabelBrain's thermal solver normalises to. |
| `export_babelbrain_weights`             | Pack PRESTUS per-element phases and particle velocity into a BabelBrain OptimizedWeightsFile (HDF5 complex64 phasors). |
| `extract_analytical_profile`           | Compute an axial intensity profile analytically (O'Neil or Rayleigh forward model) without running a k-Wave simulation; returns the same struct layout as `extract_simulated_profile`. |
| `extract_real_intensity_profile`        | Extract or spline-interpolate axial intensity profiles between available measured focal depths (wrt exit plane), aligns peaks, skips near-field for max, and plots/saves results. |
| `extract_simulated_profile`             | Extract and visualize simulated acoustic pressure data.                                         |
| `focal_reference_to_frac`               | Map a `focal_reference` config value (`flhm_3db`/`flhm_6db`) to its FLHM threshold fraction of peak; errors on `peak` (no associated fraction) or an unknown value. |
| `import_babelbrain_weights`             | Read per-element phase weights from a BabelBrain OptimizedWeightsFile (.h5), derive hardware correction relative to geometric phases, and write the correction into the transducer equipment YAML. |
| `perform_global_search`                 | Globally optimize transducer element phases (rad) and particle velocity via FEXminimize or GlobalSearch to minimize error fitting O’Neil analytical intensity to measured axial profiles. |
| `phase_optimization_annulus_full_curve` | Optimizes the phase profile for an annular transducer by matching intensity curves.             |
| `phase_optimization_annulus`            | Optimizes the phase profile by calculating focal distance error for an annular transducer.      |
| `plot_opt_sim_results`                  | Load optimized k-Wave results, generate 2D focal plane intensity maps and axial profile comparisons (simulated vs. analytical O’Neil vs. desired), and save annotated plots.  |
| `recompute_oneil_solution`              | Recompute O’Neil analytical pressure profile using optimized phases/velocity, convert to intensity, plot comparisons with original/desired profiles (skipping near-field), and report focal metrics.  |
| `save_optimized_values`                 | Append optimized transducer phases and amplitudes into a CSV table (indexed by intensity/focus) and save parameters as YAML for PRESTUS config integration. |
| `scale_real_intensity_profile`          | Linearly scale measured focal intensity profile and update transducer source amplitude (via acoustic impedance) to match desired peak Isppa (W/cm²). |
| `set_real_phases`                       | Load or compute manufacturer-specific (Sonic Concepts/Imasonic) transducer element phases for given focal depth (wrt exit plane), interpolate/unwrap for virtual elements, and return phase degrees.  |

#### CORE

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `load_parameters`                       | Loads and merges configuration files for simulation parameters; a trailing `{key, value}` cell applies overrides before validation. |
| `prestus_study_start`                   | Run the pipeline for every row of a study manifest (subjects × conditions with per-run settings), locally or on HPC, then build the study report; `DryRun` checks all rows without running. |
| `path_log_setup`                        | Set up internal paths and logging, filename for output table                                    |
| `simulation_nifti`                      | Save key outputs as 3D Niftis                                                                   |

#### GROUP

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `changem_vectorized`                    | Replace multiple old values in array A with corresponding new values.                           |
| `combine_plots_by_suffix`               | Combines subject-specific plots into a single montage image.                                    |
| `create_group_MNI_plots`                | Deprecated (older output layout): group plots in MNI space; use `create_group_maps`.            |
| `create_group_maps`                     | Group MNI maps per condition from the per-run MNI NIfTIs: −3 dB focal-volume overlap, mean/SD intensity, mean heating, plus overlap figures. |
| `study_read_manifest`                   | Read a study manifest (.tsv/.csv): one row per run with reserved columns and dotted config-key columns. |
| `study_run_table`                       | Resolve manifest rows to output locations (subject, affix from condition, medium, path.sim); rejects rows that would overwrite each other. |
| `study_row_parameters`                  | Load the full config of one study run (study config < row config < manifest columns < overrides). |
| `study_setup`                           | Load a study config, read its manifest and resolve the runs and study output folder. |
| `study_run_descriptors`                 | Per-run settings that differ between subjects (placement mode, target, transducer, depth, target intensity). |
| `study_run_files`                       | Output file names of one run (CSV, targeting sidecar, report, NIfTIs) for the group tools. |
| `prestus_group_report_start`            | Entry point for the group HTML report and group MNI maps: study mode (runs from `study.manifest`) or subject mode (auto-discovered / listed subjects). |
| `generate_group_report`                 | Assemble the self-contained group HTML report (exposure dashboard, study design, condition comparison, targeting accuracy, group maps, outliers, roster, box plots, per-run cards) and a one-row-per-run CSV. |
| `report_styles`                         | Alternative visual skins for the HTML reports (instrument / classic / compact / paper), selected via `modules.report_style`. |
| `preview_report_styles`                 | Report test mode: re-render an existing subject output folder in every skin (no re-simulation) into a sandbox folder with a side-by-side index.html. |
| `discover_group_subjects`               | Auto-discover completed subjects by scanning `path.sim` for folders with a matching per-subject CSV. |

#### GUI

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `calibration_gui`                       | Standalone GUI to fit transducer library entries to hydrophone profiles, predict profiles at any depth and add transducers. |
| `neuronav_gui`                          | Stand-alone GUI for Localite marker ingestion: detect and review series, QC, name targets, batch ingestion, MNI view. |
| `prestus_gui`                           | Point-and-click GUI to review and edit a PRESTUS config and start runs.                         |

#### HEAD

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `head_smooth_and_crop`                  | Convert and smooth segmentations into medium map & crop grid for effiency                       |
| `preproc_medium_mask`                   | Map segmentation indices onto the medium labels (in the config) [layered only]                  |
| `preproc_align_to_focal_axis`           | Rotate input image and its grid coordinates to the transducer-focal axis (+Z)                   |
| `preproc_crop_eCSF`                     | Expand CSF to define outer edges of layered medium.                                             |
| `preproc_crop_grid`                     | Crop grid to head + transducer + PML.                                                           |
| `preproc_head`                          | Preprocesses structural brain data for simulations (segmentation, alignment, cropping).         |
| `preproc_segmentation`                  | Setup SimNIBS segmentation                                                                      |
| `segmentation_run`                      | Submit SimNIBS segmentation                                                                     |
| `skull_fill_holes`                      | Fill holes in skull segmentation and between skull and skin                                     |
| `skull_rubber_wrap_visualize`           | Visualize results of skull rubber expansion                                                     |
| `skull_rubber_wrap`                     | Inflate the skull layer to locally fill potential holes                                         |
| `smooth_img`                            | Apply 3D smoothing                                                                              |

#### HELPER

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `check_availability`                    | Check file availability.                                                                        |
| `check_layers`                          | Match requested layers to those available in the segmentation. For (p)CT, homogenize multi-skull layer into a single `skull `layer.    |
| `confirm_overwriting`                   | Check overwriting (manual).                                                                     |
| `confirmation_dlg`                      | Present confirmation dialogue.                                                                  |
| `find_min_factor`                       | Find the number with the smallest maximum factor in a range.                                    |
| `get_crop_dims`                         | Computes cropping dimensions for a 3D image with a margin.                                      |
| `get_flhm_center_position`              | Calculates the center position of the full-length half-maximum (FLHM).                          |
| `get_slice_by_label`                    | Extracts a specific slice from a 3D image based on axis label and slice number.                 |
| `get_xyz_mesh`                          | Generates a mesh of 3D coordinates for a given image.                                           |
| `getidx`                                | Retrieves indices for requested tissues from a parameter structure.                             |
| `kwave_version`                         | Display k-Wave version number and (if available) git hash.                                      |
| `apply_config_overrides`                | Set config values by dotted key (`placement.localite.target`, `transducer(2).serial`); warns on unknown keys. |
| `sanitize_output_affix`                 | Normalise `io.output_affix` to the BIDS `_desc-<label>` form. |
| `mergestruct`                           | Recursively merge fields and subfields of two structures; treats `yaml.Null` as empty so YAML-loaded null values are always overwritten. |
| `log_timer`                             | Start or stop logs for benchmarking time, RAM, and disk space use.                              |
| `masked_max_3d`                         | Computes the maximum intensity within a masked 3D region.                                       |
| `read_ini_file`                         | Read an INI file into MATLAB                                                                    |
| `round_if_integer`                      | Rounds values if they are sufficiently close to integers; otherwise raises an error.            |
| `simnibs_version`                       | Get SimNIBS version of segmentation from HTML, print, allocate to `parameters`.                 |
| `subset_fields`                         | Copy designated structure field to new structure.                                               |
| `tissuemask_binary`                     | Extract binary tissue segmentation masks (for indexing)                                         |
| `upsample_to_grid`                      | Upsamples a 3D image to a higher resolution grid.                                               |
| `zip_fields`                            | Convert a structure's fields and values into a cell array.                                      |

#### HPC

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `hpc_detect_system`                     | Detect SLURM or qsub HPC system.                                                                |
| `hpc_job_info`                          | Generate formatted job display information.                                                     |
| `hpc_job_name`                          | Generate standardized HPC job name.                                                             |
| `hpc_setup_temp_files`                  | Setup directories and generate temporary files.                                                 |
| `hpc_submit_job`                        | Submit HPC batch job (SLURM or qsub).                                                           |
| `hpc_validate_parameters`               | Validate HPC job parameters.                                                                    |
| `hpc_wait_for_job`                      | Monitor HPC job until completion.                                                               |

#### MEDIUM

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `fitPowerLawParamsMulti`                | Fit power law absorption parameters for highly absorbing media.                                 |
| `get_alpha_coeff`                       | Computes the amplitude attenuation coefficient for a given medium and frequency.                |
| `medium_properties_nifti`               | Save NifTi image of the specified medium property map.                                          |
| `medium_setup`                          | Set up medium                                                                                   |
| `medium_pct_density`                    | Skull: pCT-informed density mapping                                                             |
| `medium_pct_soundspeed`                 | Skull: pCT-informed sound speed mapping                                                         |
| `medium_pct_attenuation`                | Skull: pCT-informed attenuation mapping                                                         |

#### NEURONAV

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `neuronav_compute_series_statistics`    | Compute mean position and variability over stimulus train                                       |
| `neuronav_convert_MNI_to_native`        | Transform coordinates from MNI space to native subject space                                    |
| `neuronav_convert_native_to_MNI`        | Convert native RAS coordinates to MNI space                                                     |
| `neuronav_convert_trigger_to_voxels`    | Convert Localite trigger positions to voxel (image) coordinates.                                |
| `neuronav_detect_series`                | Detect the recorded series of one session (QC, suggested names) without writing.                |
| `neuronav_export_session_csv`           | Export per-session coordinate arrays to CSV with labeled voxel and RAS (mm) positions.          |
| `neuronav_find_planning_t1`             | Find the planning T1 loaded into Localite (`path.localite_t1`, else `path.localite_raw`).       |
| `neuronav_get_group_mean_mni`           | Compute group-level average MNI coordinates (both mm and voxel).                                |
| `neuronav_ingest_batch`                 | Ingest and export Localite targets for many subjects/sessions; writes `neuronav_ingest_summary.csv`. |
| `neuronav_label_map_from_scan`          | Build a `neuronav.target_labels` map from a label scan.                                         |
| `neuronav_localite_pose`                | Resolve the recorded Localite instrument pose for placement (file mode, or ingested markers; ingests on first use). |
| `neuronav_map_label`                    | Map a Localite target label to the study target name (`neuronav.target_labels`).                |
| `neuronav_name_series`                  | Name the detected series of a session (coil table, planned target, default; shared GUI/script rule). |
| `neuronav_read_label_candidates`        | Collect the named planned targets of a session (EntryTarget, GUMMarkers).                       |
| `neuronav_resolve_paths`                | Resolve the Localite raw/ingested/T1/export folders from `path.localite_*` (with `path.localite` fallback). |
| `neuronav_scan_sessions`                | Inventory raw Localite data and planned-target labels.                                          |
| `neuronav_select_localite`              | Select most recent Localite XML for a session.                                                  |
| `neuronav_suggest_target_names`         | Heuristic target names for ingested series (planned target nearest the beam axis).              |
| `neuronav_summary_to_mni`               | Add MNI coordinates to a neuronav ingest summary.                                               |
| `position_transducer_localite`          | Determines transducer and focus positions in voxel space using Localite data and MRI header.    |

#### PSEUDO-CT

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `fit_pairwiselinear`                    | Perform a pairwise linear fit between HU and density values with optional plot.                 |
| `pct_create_pseudoCT`                   | Generate Hounsfield pseudoCT from SimNIBS PETRA-UTE (replacing T2) via N4 bias correction, linear skull mapping, partial volume correction, smoothing, and tissue masks. |
| `pct_skullexpand`                       | Load SimNIBS NIfTIs + headers, run skull rubber wrap, save output.                              |
| `pct_skullmapping`                      | Computes pseudo-CT mapping for cortical and trabecular bone using UTE histograms. [deprecated, debug]  |
| `pct_soft_tissue_peak`                  | Identifies the soft tissue peak from UTE intensity distribution histograms.                     |

#### PLOT

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `plot_coronal_slices`                   | (DEPRECATED) Visualizes coronal slices of a 3D image with optional legends for labeled images.  |
| `plot_median_montage`                   | (DEPRECATED) Creates a montage of central slices from a 3D T1 image.                            |
| `plot_overlay_2d`                       | Overlays map on a background image slice with key positions highlighted.                        |
| `plot_overlay`                          | Visualizes map overlaid on a 2D slice of a 3D background image with advanced options.           |
| `plot_t1_with_transducer`               | Creates a plot of a T1 slice oriented along the transducer's axis with overlays.                |
| `plot_transducer_overlay`               | Visualize curved transducer geometry with exit plane.                                           |
| `show_3d_head`                          | Visualizes segmented brain images in 3D with transducer placement and target location.          |
| `show_positioning_plots`                | Visualizes transducer positioning before and after preprocessing using segmented images.        |

#### SOURCE

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `grid_axisymmetry`                      | Convert grid to axisymmetry                                                                     |
| `grid_tissue_setup`                     | Set up grid dimensions, preprocess head (if modeled), and place in grid                         |
| `grid_transducer_location`              | Position transducer (and target) in simulation grid                                             |
| `source_create`                         | Creates ultrasound source signals and masks for k-Wave simulations based on transducer geometry.|
| `source_sensor_setup`                   | Configure a k-Wave simulation grid, transducer sources, and sensor mask based on input parameters for recording pressure fields in 2D or 3D ultrasonic neuromodulation setups. |

#### THERMAL

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `thermal_analysis`                      | Analyze output of thermal simulations                                                           |
| `thermal_parameters`                    | Checks and converts protocol timing setup for thermal estimation                                |
| `thermal_plot_protocol`                 | Visualize the requested protocol timing                                                         |
| `thermal_plot_sim`                      | Visualizes heating simulation results over time (temperature, rise, CEM43). Optional video      |
| `thermal_simulation`                    | Simulate ultrasound-induced heating and thermal dose (CEM43) using k-Wave’s kWaveDiffusion from acoustic pressure data across pulsed train repetitions.  |
| `thermal_update_timeseries`             | Track max T/CEM43 per tissue layer.                                                             |

#### TRANSDUCER

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `focal_distance_calculation`            | Compute expected focal distances for (multi-)transducer setup.                                  |
| `get_arc`                               | Generates the coordinates of an arc in 2D space.                                                |
| `get_trans_pos_from_trigger_markers`    | Determines transducer and target positions using Localite trigger marker files.                 |
| `get_transducer_box`                    | Computes transducer box dimensions and positions for simulations.                               |
| `transducer_analyze_position_fast`      | Analyzes the transducer position relative to the skull and skin.                                |
| `transducer_analyze_position`           | Computes geometric and statistical measures for a transducer position.                          |
| `transducer_positioning`                | Determines heuristic transducer placement for a given target.                                   |
|   `tp_candidate_mesh`                   | Build coordinate mesh and candidate transducer geometry.                                        |
|   `tp_evaluate_candidate_positions`     | Evaluate candidate transducer positions.                                                        |
|   `tp_find_initial_candidate`           | Find candidate transducer positions on skull surface.                                           |
|   `tp_plot_candidate_positions`         | Plot transducer candidate positions on skull surface slice.                                     |
|   `tp_plot_geometry_overlay`            | Visualize transducer geometry on skin segmentation slice.                                       |
|   `tp_plot_heuristic_position`          | Create visualization of heuristic transducer placement.                                         |
|   `tp_remove_ear_locations`             | Exclude ear entries from heuristic transducer positions.                                        |
|   `tp_select_heuristic_position`        | Select optimal transducer position and export Localite coordinates.                             |
| `transducer_setup`                      | Create a transducer mask and label matrix for a computational grid.                             |

#### TRANSFORM

| **Function Name**                       | **Description**                                                                                 |
|-----------------------------------------|-------------------------------------------------------------------------------------------------|
| `sample_volume_at_world`                | Nearest-neighbour values of a NIfTI volume at world (RAS/MNI mm) positions; NaN outside. |
| `convert_2d_to_axisymmetric`            | Converts 2D k-Wave simulation grid, medium properties, and source to axisymmetric form by halving the shorter dimension (right half from center). |
| `convert_axisymmetric_to_2d`            | Mirrors axisymmetric simulation data (sensor, medium, masks, source) left-right to full 2D, doubles radial dimension, transposes axes, and updates kgrid/positions. |
| `convert_axisymmetric_to_3d`            | Expands 2D axisymmetric data (via radialExpand2DTo3D) to full 3D grid by duplicating radial dimension, updates positions/grid/kgrid for cubic symmetry. |
| `convert_final_to_MNI_matlab`           | Converts an image from subject space to MNI space using MATLAB.                                 |
| `convert_final_to_MNI_simnibs`          | Converts an image to MNI space using SimNIBS.                                                   |
| `mni2subject_coords_LDfix`              | Transforms a set of coordinates in MNI space to subject space.                                  |
| `radialExpand2DTo3D`                    | Radially expands 2D axisymmetric data into 3D Cartesian volume.                                 |
| `grid_to_ras`                           | Converts 1-based voxel indices to RAS coordinates (inverse of `ras_to_grid`).                   |
| `nifti_affine_from_1based`              | Turns an affine on 1-based voxel indices into a NIfTI header affine (0-based indices).          |
| `ras_to_grid`                           | Converts RAS coordinates to 1-based voxel (grid) indices using the NIfTI header affine.         |
| `subject2mni_coords_LDfix`              | Transforms a set of coordinates in MNI space to subject space.                                  |
| `transform_coordinates`                 | Transform input coordinates between coordinate systems (wrapper) space.                         |