## PRESTUS function documentation

The following documents the functions provided in PRESTUS.

| **Function Name**                       | **Category**    | **Description**                                                                                |
|-----------------------------------------|-----------------|------------------------------------------------------------------------------------------------|
| `acoustic_analysis`                     | ACOUSTIC       |                                                                                                 |
| `acoustic_convert_axisymmetry`          | ACOUSTIC       |                                                                                                 |
| `acoustic_simulation`                   | ACOUSTIC       | Set up acoustic simulation                                                                      |
| `acoustic_wrapper`                      | ACOUSTIC       | Set up acoustic simulation                                                                      |
| `calc_virtual_elem`                     | CALIBRATION    | (STANDALONE) Calculate additional virtual elements                                              |
| `calibration_transducer`                | CALIBRATION    |                                                                                                 |
| `compute_oneil_solution`                | CALIBRATION    |                                                                                                 |
| `compute_phases`                        | CALIBRATION    |                                                                                                 |
| `extract_real_intensity_profile`        | CALIBRATION    |                                                                                                 |
| `perform_global_search`                 | CALIBRATION    |                                                                                                 |
| `perform_global_search`                 | CALIBRATION    |                                                                                                 |
| `phase_optimization_annulus_full_curve` | CALIBRATION    | Optimizes the phase profile for an annular transducer by matching intensity curves.             |
| `phase_optimization_annulus`            | CALIBRATION    | Optimizes the phase profile by calculating focal distance error for an annular transducer.      |
| `plot_opt_sim_results`                  | CALIBRATION    |                                                                                                 |
| `recalculate_analytical_sol`            | CALIBRATION    |                                                                                                 |
| `save_optimized_values`                 | CALIBRATION    |                                                                                                 |
| `scale_real_intensity_profile`          | CALIBRATION    |                                                                                                 |
| `set_real_phases`                       | CALIBRATION    |                                                                                                 |
| `load_parameters`                       | CORE           | Loads and merges configuration files for simulation parameters.                                 |
| `path_log_setup`                        | CORE           | Set up internal paths and logging, filename for output table                                    |
| `simulation_nifti`                      | CORE           | Save key outputs as 3D Niftis                                                                   |
| `changem_vectorized`                    | GROUP          |                                                                                                 |
| `combine_plots_by_suffix`               | GROUP          |                                                                                                 |
| `create_group_MNI_plots`                | GROUP          |                                                                                                 |
| `check_availability`                    | HELPER         |                                                                                                 |
| `confirm_overwriting`                   | HELPER         |                                                                                                 |
| `confirmation_dlg`                      | HELPER         |                                                                                                 |
| `find_min_factor`                       | HELPER         |                                                                                                 |
| `get_crop_dims`                         | HELPER         | Computes cropping dimensions for a 3D image with a margin.                                      |
| `get_flhm_center_position`              | HELPER         | Calculates the center position of the full-length half-maximum (FLHM).                          |
| `get_slice_by_label`                    | HELPER         | Extracts a specific slice from a 3D image based on axis label and slice number.                 |
| `get_xyz_mesh`                          | HELPER         | Generates a mesh of 3D coordinates for a given image.                                           |
| `getidx`                                | HELPER         | Retrieves indices for requested tissues from a parameter structure.                             |
| `masked_max_3d`                         | HELPER         | Computes the maximum intensity within a masked 3D region.                                       |
| `mergeStructure`                        | HELPER         | Merges multiple scalar structures into one.                                                     |
| `read_ini_file`                         | HELPER         |                                                                                                 |
| `round_if_integer`                      | HELPER         | Rounds values if they are sufficiently close to integers; otherwise raises an error.            |
| `subset_fields`                         | HELPER         |                                                                                                 |
| `tissuemask_binary`                     | HELPER         | Extract binary tissue segmentation masks (for indexing)                                         |
| `upsample_to_grid`                      | HELPER         |                                                                                                 |
| `zip_fields`                            | HELPER         |                                                                                                 |
| `single_subject_pipeline_with_qsub`     | HPC            |                                                                                                 |
| `single_subject_pipeline_with_slurm`    | HPC            |                                                                                                 |
| `transducer_positioning_with_qsub`      | HPC            | Submits transducer positioning jobs to Qsub cluster using batch scripts.                        |
| `transducer_positioning_with_slurm`     | HPC            | Submits transducer positioning jobs to SLURM cluster using batch scripts.                       |
| `fitPowerLawParamsMulti`                | MEDIUM         |                                                                                                 |
| `get_alpha_coeff`                       | MEDIUM         |                                                                                                 |
| `medium_mask_create`                    | MEDIUM         | Map segmentation indices onto the medium labels in the config                                   |
| `medium_properties_nifti`               | MEDIUM         | Save NifTi image of the specified medium property map.                                          |
| `medium_setup`                          | MEDIUM         | Set up medium                                                                                   |
| `neuronav_compute_series_statistics`    | NEURONAV       |                                                                                                 |
| `neuronav_convert_MNI_to_native`        | NEURONAV       |                                                                                                 |
| `neuronav_convert_native_to_MNI`        | NEURONAV       |                                                                                                 |
| `neuronav_convert_trigger_to_voxels`    | NEURONAV       |                                                                                                 |
| `neuronav_create_marker_average`        | NEURONAV       |                                                                                                 |
| `neuronav_export_session_csv`           | NEURONAV       |                                                                                                 |
| `neuronav_get_group_mean_mni`           | NEURONAV       |                                                                                                 |
| `neuronav_select_and_average_localite`  | NEURONAV       |                                                                                                 |
| `position_transducer_localite`          | NEURONAV       | Determines transducer and focus positions in voxel space using Localite data and MRI header.    |
| `pct_create_pseudoCT`                   | PSEUDO-CT      |                                                                                                 |
| `pct_skullmapping`                      | PSEUDO-CT      | Computes pseudo-CT mapping for cortical and trabecular bone using UTE histograms.               |
| `pct_soft_tissue_peak`                  | PSEUDO-CT      | Identifies the soft tissue peak from UTE intensity distribution histograms.                     |
| `plot_coronal_slices`                   | PLOT           | (DEPRECATED) Visualizes coronal slices of a 3D image with optional legends for labeled images.  |
| `plot_median_montage`                   | PLOT           | (DEPRECATED) Creates a montage of central slices from a 3D T1 image.                            |
| `plot_overlay_2d`                       | PLOT           | Overlays map on a background image slice with key positions highlighted.                        |
| `plot_overlay`                          | PLOT           | Visualizes map overlaid on a 2D slice of a 3D background image with advanced options.           |
| `plot_t1_with_transducer`               | PLOT           | Creates a plot of a T1 slice oriented along the transducer's axis with overlays.                |
| `show_3d_head`                          | PLOT           | Visualizes segmented brain images in 3D with transducer placement and target location.          |
| `show_positioning_plots`                | PLOT           | Visualizes transducer positioning before and after preprocessing using segmented images.        |
| `smooth_img`                            | PLOT           |                                                                                                 |
| `preproc_align_to_focal_axis`           | PREPROC        |                                                                                                 |
| `preproc_head`                          | PREPROC        | Preprocesses structural brain data for simulations (segmentation, alignment, cropping).         |
| `preproc_segmentation`                  | PREPROC        |                                                                                                 |
| `preproc_smooth_and_crop`               | PREPROC        |                                                                                                 |
| `segmentation_run`                      | PREPROC        |                                                                                                 |
| `grid_axisymmetry`                      | SOURCE         | Convert grid to axisymmetry                                                                     |
| `grid_tissue_setup`                     | SOURCE         | Set up grid dimensions, preprocess head (if modeled), and place in grid                         |
| `grid_transducer_location`              | SOURCE         | Position transducer (and target) in simulation grid                                             |
| `source_create`                         | SOURCE         | Creates ultrasound source signals and masks for k-Wave simulations based on transducer geometry.|
| `source_sensor_setup`                   | SOURCE         |                                                                                                 |
| `thermal_analysis`                      | THERMAL        | Analyze output of thermal simulations                                                           |
| `thermal_parameters`                    | THERMAL        | Checks and converts protocol timing setup for thermal estimation                                |
| `thermal_plot_protocol`                 | THERMAL        | Visualize the requested protocol timing                                                         |
| `thermal_plot_sim`                      | THERMAL        | Visualizes heating simulation results over time (temperature, rise, CEM43). Optional video      |
| `thermal_simulation`                    | THERMAL        |                                                                                                 |
| `focal_distance_calculation`            | TRANSDUCER     |                                                                                                 |
| `get_trans_pos_from_trigger_markers`    | TRANSDUCER     | Determines transducer and target positions using Localite trigger marker files.                 |
| `get_transducer_box`                    | TRANSDUCER     | Computes transducer box dimensions and positions for simulations.                               |
| `get_arc`                               | TRANSDUCER     |                                                                                                 |
| `transducer_analyze_position_fast`      | TRANSDUCER     |                                                                                                 |
| `transducer_analyze_position`           | TRANSDUCER     |                                                                                                 |
| `transducer_positioning`                | TRANSDUCER     |                                                                                                 |
| `transducer_setup`                      | TRANSDUCER     |                                                                                                 |
| `convert_2d_to_axisymmetric`            | TRANSFORM      |                                                                                                 |
| `convert_axisymmetric_to_2d`            | TRANSFORM      |                                                                                                 |
| `convert_axisymmetric_to_3d`            | TRANSFORM      |                                                                                                 |
| `convert_final_to_MNI_matlab`           | TRANSFORM      |                                                                                                 |
| `convert_final_to_MNI_simnibs`          | TRANSFORM      |                                                                                                 |
| `mni2subject_coords_LDfix`              | TRANSFORM      |                                                                                                 |
| `radialExpand2DTo3D`                    | TRANSFORM      |                                                                                                 |
| `ras_to_grid`                           | TRANSFORM      | Converts RAS coordinates to voxel (grid) coordinates using NIfTI header transformation matrix.  |
| `subject2mni_coords_LDfix`              | TRANSFORM      |                                                                                                 |
