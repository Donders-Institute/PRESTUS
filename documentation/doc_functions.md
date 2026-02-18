## PRESTUS function documentation

The following documents the functions provided in PRESTUS.

| **Function Name**                       | **Category**    | **Description**                                                                                |
|-----------------------------------------|-----------------|------------------------------------------------------------------------------------------------|
| `acoustic_analysis`                     | ACOUSTIC       | Compute acoustic safety metrics (Isppa, MI, max pressure) from sensor data, extract tissue-specific maxima, focal distances, and generate overlaid ISPPA plots on segmentation                                                                                           |
| `acoustic_convert_axisymmetry`          | ACOUSTIC       |  Expand k-Wave axisymmetric (2D rotational) simulation outputs to full 3D (via  convert_axisymmetric_to_3d ) for heating sims or 2D Cartesian (via  convert_axisymmetric_to_2d ) otherwise, updating sensor data, parameters, masks, medium, grid, source, and labels.    |
| `acoustic_simulation`                   | ACOUSTIC       | Set up acoustic simulation                                                                      |
| `acoustic_wrapper`                      | ACOUSTIC       | Set up acoustic simulation                                                                      |
| `calc_virtual_elem`                     | CALIBRATION    | (STANDALONE) Calculate additional virtual elements                                              |
| `calibration_transducer`                | CALIBRATION    | Run water simulations to optimize transducer source amplitude and element phases matching a target intensity profile via O’Neil analytical solution, global search, and re-simulation.                                                                                    |
| `compute_oneil_solution`                | CALIBRATION    | Compute analytical O’Neil pressure along the beam axis for focused annular transducers in water, derives particle velocity and grid adjustment factor, and plots comparisons with simulated/desired intensities.                                                      |
| `compute_phases`                        | CALIBRATION    | Calculate per-element transducer phases (degrees) for electronic focusing by computing fractional wavelength delays to steer waves constructively to a target point relative to the exit plane.                                                                           |
| `extract_real_intensity_profile`        | CALIBRATION    | Extract or spline-interpolate axial intensity profiles between available measured focal depths (wrt exit plane), aligns peaks, skips near-field for max, and plots/saves results.                                                                                     |
| `perform_global_search`                 | CALIBRATION    | Globally optimize transducer element phases (rad) and particle velocity via FEXminimize or GlobalSearch to minimize error fitting O’Neil analytical intensity to measured axial profiles.                                                                               |
| `phase_optimization_annulus_full_curve` | CALIBRATION    | Optimizes the phase profile for an annular transducer by matching intensity curves.             |
| `phase_optimization_annulus`            | CALIBRATION    | Optimizes the phase profile by calculating focal distance error for an annular transducer.      |
| `plot_opt_sim_results`                  | CALIBRATION    | Load optimized k-Wave results, generate 2D focal plane intensity maps and axial profile comparisons (simulated vs. analytical O’Neil vs. desired), and save annotated plots.                                                                                     |
| `recalculate_analytical_sol`            | CALIBRATION    | Recompute O’Neil analytical pressure profile using optimized phases/velocity, convert to intensity, plot comparisons with original/desired profiles (skipping near-field), and report focal metrics.                                                                  |
| `save_optimized_values`                 | CALIBRATION    | Append optimized transducer phases and amplitudes into a CSV table (indexed by intensity/focus) and save parameters as YAML for PRESTUS config integration.                                                                                                           |
| `scale_real_intensity_profile`          | CALIBRATION    | Linearly scale measured focal intensity profile and update transducer source amplitude (via acoustic impedance) to match desired peak Isppa (W/cm²).                                                                                                              |
| `set_real_phases`                       | CALIBRATION    | Load or compute manufacturer-specific (Sonic Concepts/Imasonic) transducer element phases for given focal depth (wrt exit plane), interpolate/unwrap for virtual elements, and return phase degrees.                                                                   |
| `load_parameters`                       | CORE           | Loads and merges configuration files for simulation parameters.                                 |
| `path_log_setup`                        | CORE           | Set up internal paths and logging, filename for output table                                    |
| `simulation_nifti`                      | CORE           | Save key outputs as 3D Niftis                                                                   |
| `changem_vectorized`                    | GROUP          | Replace multiple old values in array A with corresponding new values.                           |
| `combine_plots_by_suffix`               | GROUP          | Combines subject-specific plots into a single montage image.                                    |
| `create_group_MNI_plots`                | GROUP          | Generate group-level plots in MNI space for multiple subjects.                                  |
| `head_smooth_and_crop`                  | HEAD           | Convert and smooth segmentations into medium map & crop grid for effiency                       |
| `preproc_medium_mask`                   | HEAD           | Map segmentation indices onto the medium labels (in the config) [layered only]                  |
| `preproc_align_to_focal_axis`           | HEAD           | Rotate input image and its grid coordinates to the transducer-focal axis (+Z)                   |
| `preproc_crop_eCSF`                     | HEAD           | Expand CSF to define outer edges of layered medium.                                             |
| `preproc_crop_grid`                     | HEAD           | Crop grid to head + transducer + PML.                                                           |
| `preproc_head`                          | HEAD           | Preprocesses structural brain data for simulations (segmentation, alignment, cropping).         |
| `preproc_segmentation`                  | HEAD           | Setup SimNIBS segmentation                                                                      |
| `segmentation_run`                      | HEAD           | Submit SimNIBS segmentation                                                                     |
| `skull_fill_holes`                      | HEAD           | Fill holes in skull segmentation and between skull and skin                                     |
| `skull_rubber_wrap_visualize`           | HEAD           | Visualize results of skull rubber expansion                                                     |
| `skull_rubber_wrap`                     | HEAD           | Inflate the skull layer to locally fill potential holes                                         |
| `smooth_img`                            | HEAD           | Apply 3D smoothing                                                                              |
| `check_availability`                    | HELPER         | Check file availability.                                                                        |
| `confirm_overwriting`                   | HELPER         | Check overwriting (manual).                                                                     |
| `confirmation_dlg`                      | HELPER         | Present confirmation dialogue.                                                                  |
| `find_min_factor`                       | HELPER         | Find the number with the smallest maximum factor in a range.                                    |
| `get_crop_dims`                         | HELPER         | Computes cropping dimensions for a 3D image with a margin.                                      |
| `get_flhm_center_position`              | HELPER         | Calculates the center position of the full-length half-maximum (FLHM).                          |
| `get_slice_by_label`                    | HELPER         | Extracts a specific slice from a 3D image based on axis label and slice number.                 |
| `get_xyz_mesh`                          | HELPER         | Generates a mesh of 3D coordinates for a given image.                                           |
| `getidx`                                | HELPER         | Retrieves indices for requested tissues from a parameter structure.                             |
| `kwave_version`                         | HELPER         | Display k-Wave version number and (if available) git hash.                                      |
| `log_timer`                             | HELPER         | Start or stop logs for benchmarking time, RAM, and disk space use.                              |
| `masked_max_3d`                         | HELPER         | Computes the maximum intensity within a masked 3D region.                                       |
| `mergeStructure`                        | HELPER         | Merges multiple scalar structures into one.                                                     |
| `read_ini_file`                         | HELPER         | Read an INI file into MATLAB                                                                    |
| `round_if_integer`                      | HELPER         | Rounds values if they are sufficiently close to integers; otherwise raises an error.            |
| `subset_fields`                         | HELPER         | Copy designated structure field to new structure.                                               |
| `tissuemask_binary`                     | HELPER         | Extract binary tissue segmentation masks (for indexing)                                         |
| `upsample_to_grid`                      | HELPER         | Upsamples a 3D image to a higher resolution grid.                                               |
| `zip_fields`                            | HELPER         | Convert a structure's fields and values into a cell array.                                      |
| `single_subject_pipeline_with_qsub`     | HPC            | HPC call of pipeline with QSUB.                                                                 |
| `single_subject_pipeline_with_slurm`    | HPC            | HPC call of pipeline with SLURM.                                                                |
| `transducer_positioning_with_qsub`      | HPC            | Submits transducer positioning jobs to Qsub cluster using batch scripts.                        |
| `transducer_positioning_with_slurm`     | HPC            | Submits transducer positioning jobs to SLURM cluster using batch scripts.                       |
| `fitPowerLawParamsMulti`                | MEDIUM         | Fit power law absorption parameters for highly absorbing media.                                 |
| `get_alpha_coeff`                       | MEDIUM         | Computes the amplitude attenuation coefficient for a given medium and frequency.                |
| `medium_properties_nifti`               | MEDIUM         | Save NifTi image of the specified medium property map.                                          |
| `medium_setup`                          | MEDIUM         | Set up medium                                                                                   |
| `neuronav_compute_series_statistics`    | NEURONAV       | Compute mean position and variability over stimulus train                                       |
| `neuronav_convert_MNI_to_native`        | NEURONAV       | Transform coordinates from MNI space to native subject space                                    |
| `neuronav_convert_native_to_MNI`        | NEURONAV       | Convert native RAS coordinates to MNI space                                                     |
| `neuronav_convert_trigger_to_voxels`    | NEURONAV       | Convert Localite trigger positions to voxel (image) coordinates.                                |
| `neuronav_create_marker_average`        | NEURONAV       | Build an averaged Localite-trigger structure for export.                                        |
| `neuronav_export_session_csv`           | NEURONAV       | Export per-session coordinate arrays to CSV with labeled voxel and RAS (mm) positions.          |
| `neuronav_get_group_mean_mni`           | NEURONAV       | Compute group-level average MNI coordinates (both mm and voxel).                                |
| `neuronav_select_and_average_localite`  | NEURONAV       | Select most recent Localite XML for a session.                                                  |
| `position_transducer_localite`          | NEURONAV       | Determines transducer and focus positions in voxel space using Localite data and MRI header.    |
| `pct_create_pseudoCT`                   | PSEUDO-CT      | Generate Hounsfield pseudoCT from SimNIBS PETRA-UTE (replacing T2) via N4 bias correction, soft/skull linear mapping (5 algorithms), partial volume correction, smoothing, and tissue masks.                                                                              |
| `pct_skullmapping`                      | PSEUDO-CT      | Computes pseudo-CT mapping for cortical and trabecular bone using UTE histograms.               |
| `pct_soft_tissue_peak`                  | PSEUDO-CT      | Identifies the soft tissue peak from UTE intensity distribution histograms.                     |
| `fit_pairwiselinear`                    | PSEUDO-CT      | Perform a pairwise linear fit between HU and density values with optional plot.                 |
| `plot_coronal_slices`                   | PLOT           | (DEPRECATED) Visualizes coronal slices of a 3D image with optional legends for labeled images.  |
| `plot_median_montage`                   | PLOT           | (DEPRECATED) Creates a montage of central slices from a 3D T1 image.                            |
| `plot_overlay_2d`                       | PLOT           | Overlays map on a background image slice with key positions highlighted.                        |
| `plot_overlay`                          | PLOT           | Visualizes map overlaid on a 2D slice of a 3D background image with advanced options.           |
| `plot_t1_with_transducer`               | PLOT           | Creates a plot of a T1 slice oriented along the transducer's axis with overlays.                |
| `show_3d_head`                          | PLOT           | Visualizes segmented brain images in 3D with transducer placement and target location.          |
| `show_positioning_plots`                | PLOT           | Visualizes transducer positioning before and after preprocessing using segmented images.        |
| `grid_axisymmetry`                      | SOURCE         | Convert grid to axisymmetry                                                                     |
| `grid_tissue_setup`                     | SOURCE         | Set up grid dimensions, preprocess head (if modeled), and place in grid                         |
| `grid_transducer_location`              | SOURCE         | Position transducer (and target) in simulation grid                                             |
| `source_create`                         | SOURCE         | Creates ultrasound source signals and masks for k-Wave simulations based on transducer geometry.|
| `source_sensor_setup`                   | SOURCE         | Configure a k-Wave simulation grid, transducer sources, and sensor mask based on input parameters for recording pressure fields in 2D or 3D ultrasonic neuromodulation setups.                                                                                     |
| `thermal_analysis`                      | THERMAL        | Analyze output of thermal simulations                                                           |
| `thermal_parameters`                    | THERMAL        | Checks and converts protocol timing setup for thermal estimation                                |
| `thermal_plot_protocol`                 | THERMAL        | Visualize the requested protocol timing                                                         |
| `thermal_plot_sim`                      | THERMAL        | Visualizes heating simulation results over time (temperature, rise, CEM43). Optional video      |
| `thermal_simulation`                    | THERMAL        | Simulate ultrasound-induced heating and thermal dose (CEM43) using k-Wave’s kWaveDiffusion from acoustic pressure data across pulsed train repetitions.                                                                                                               |
| `focal_distance_calculation`            | TRANSDUCER     | Compute expected focal distances for (multi-)transducer setup.                                  |
| `get_arc`                               | TRANSDUCER     | Generates the coordinates of an arc in 2D space.                                                |
| `get_trans_pos_from_trigger_markers`    | TRANSDUCER     | Determines transducer and target positions using Localite trigger marker files.                 |
| `get_transducer_box`                    | TRANSDUCER     | Computes transducer box dimensions and positions for simulations.                               |
| `transducer_analyze_position_fast`      | TRANSDUCER     | Analyzes the transducer position relative to the skull and skin.                                |
| `transducer_analyze_position`           | TRANSDUCER     | Computes geometric and statistical measures for a transducer position.                          |
| `transducer_positioning`                | TRANSDUCER     |  Determines heuristic transducer placement for a given target.                                  |
| `transducer_setup`                      | TRANSDUCER     | Create a transducer mask and label matrix for a computational grid.                             |
| `convert_2d_to_axisymmetric`            | TRANSFORM      | Converts 2D k-Wave simulation grid, medium properties, and source to axisymmetric form by halving the shorter dimension (right half from center).                                                                                                                  |
| `convert_axisymmetric_to_2d`            | TRANSFORM      | Mirrors axisymmetric simulation data (sensor, medium, masks, source) left-right to full 2D, doubles radial dimension, transposes axes, and updates kgrid/positions.                                                                                                     |
| `convert_axisymmetric_to_3d`            | TRANSFORM      | Expands 2D axisymmetric data (via radialExpand2DTo3D) to full 3D grid by duplicating radial dimension, updates positions/grid/kgrid for cubic symmetry.                                                                                                             |
| `convert_final_to_MNI_matlab`           | TRANSFORM      | Converts an image from subject space to MNI space using MATLAB.                                 |
| `convert_final_to_MNI_simnibs`          | TRANSFORM      | Converts an image to MNI space using SimNIBS.                                                   |
| `mni2subject_coords_LDfix`              | TRANSFORM      | Transforms a set of coordinates in MNI space to subject space.                                  |
| `radialExpand2DTo3D`                    | TRANSFORM      | Radially expands 2D axisymmetric data into 3D Cartesian volume.                                 |
| `ras_to_grid`                           | TRANSFORM      | Converts RAS coordinates to voxel (grid) coordinates using NIfTI header transformation matrix.  |
| `subject2mni_coords_LDfix`              | TRANSFORM      | Transforms a set of coordinates in MNI space to subject space.                                  |