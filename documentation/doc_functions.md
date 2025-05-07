## PRESTUS function documentation

The following documents the functions provided in PRESTUS.

| **Function Name**                      | **Category**   | **Description**                                                                                   |
|----------------------------------------|----------------|---------------------------------------------------------------------------------------------------|
| `align_to_focus_axis_and_scale`         | Other          |                        |
| `analyze_transducer_position_fast`      | Other          |                        |
| `analyze_transducer_position`           | Other          |                        |
| `bounded_interp1`                       | Other          |                        |
| `changem_vectorized`                    | Other          |                        |
| `check_thermal_parameters`              | Other          |                        |
| `combine_plots_by_suffix`               | Other          |                        |
| `confirm_overwriting`                   | Other          |                        |
| `confirmation_dlg`                      | Other          |                        |
| `convert_final_to_MNI_matlab`           | Other          |                        |
| `convert_final_to_MNI_simnibs`          | Other          |                        |
| `create_group_MNI_plots`                | Other          |                        |
| `create_pseudoCT`                       | Other          |                        |
| `find_min_factor`                       | Other          |                        |
| `fitPowerLawParamsMulti`                | Other          |                        |
| `get_alpha_coeff`                       | Other          |                        |
| `get_crop_dims`                         | Helper         | Computes cropping dimensions for a 3D image with a margin.                                       |
| `get_flhm_center_position`              | Helper         | Calculates the center position of the full-length half-maximum (FLHM).                          |
| `get_slice_by_label`                    | Helper         | Extracts a specific slice from a 3D image based on axis label and slice number.                 |
| `get_trans_pos_from_trigger_markers`    | LOCALITE       | Determines transducer and target positions using Localite trigger marker files.                 |
| `get_transducer_box`                    | TRANSDUCER     | Computes transducer box dimensions and positions for simulations.                               |
| `get_xyz_mesh`                          | Helper         | Generates a mesh of 3D coordinates for a given image.                                           |
| `getArc`                                | TRANSDUCER     |                        |
| `getidx`                                | Helper         | Retrieves indices for requested tissues from a parameter structure.                             |
| `load_parameters`                       | Helper         | Loads and merges configuration files for simulation parameters.                                 |
| `masked_max_3d`                         | Helper         | Computes the maximum intensity within a masked 3D region.                                       |
| `median_montage`                        | Plotting       | Creates a montage of central slices from a 3D T1 image.                                         |
| `medium_mask_create`                    | Helper       | Map segmentation indices onto the medium labels in the config |
| `mergeStructure`                        | Helper         | Merges multiple scalar structures into one.                                                     |
| `mni2subject_coords_LDfix`              | Helper         |                        |
| `pct_skullmapping`                      | PSEUDO-CT      | Computes pseudo-CT mapping for cortical and trabecular bone using UTE histograms.              |
| `pct_soft_tissue_peak`                  | PSEUDO-CT      | Identifies the soft tissue peak from UTE intensity distribution histograms.                     |
| `phase_optimization_annulus_full_curve` | CALIBRATION    | Optimizes the phase profile for an annular transducer by matching intensity curves.             |
| `phase_optimization_annulus`            | CALIBRATION    | Optimizes the phase profile by calculating focal distance error for an annular transducer.      |
| `plot_coronal_slices`                   | Plotting       | Visualizes coronal slices of a 3D image with optional legends for labeled images.               |
| `plot_heating_sims`                     | Plotting       | Visualizes heating simulation results over time (temperature, rise, CEM43).                    |
| `plot_isppa_over_image_2d`              | Plotting       | Overlays ISppa map on a background image slice with key positions highlighted.                  |
| `plot_isppa_over_image`                 | Plotting       | Visualizes ISppa map overlaid on a 2D slice of a 3D background image with advanced options.     |
| `plot_t1_with_transducer`               | Plotting       | Creates a plot of a T1 slice oriented along the transducer's axis with overlays.                |
| `position_transducer_localite`          | LOCALITE       | Determines transducer and focus positions in voxel space using Localite data and MRI header.    |
| `preprocess_brain`                      | MAPPING        | Preprocesses structural brain data for simulations (segmentation, alignment, cropping).         |
| `ras_to_grid`                           | Helper         | Converts RAS coordinates to voxel (grid) coordinates using NIfTI header transformation matrix.  |
| `round_if_integer`                      | Helper         | Rounds values if they are sufficiently close to integers; otherwise raises an error.            |
| `run_heating_simulations`               | HEATING        |                        |
| `run_segmentation`                      | SEGMENTATION   |                        |
| `run_simulations`                       | ACOUSTIC       |                        |
| `setup_grid_source_sensor`              | Other          |                        |
| `setup_medium`                          | Other          |                        |
| `setup_source`                          | Other          | Creates ultrasound source signals and masks for k-Wave simulations based on transducer geometry.|
| `show_3d_head`                          | Plotting       | Visualizes segmented brain images in 3D with transducer placement and target location.          |
| `show_positioning_plots`                | Plotting       | Visualizes transducer positioning before and after preprocessing using segmented images.         |
| `single_subject_pipeline_with_qsub`     | Other          |                        |
| `single_subject_pipeline_with_slurm`    | Other          |                        |
| `smooth_and_crop`                       | Other          |                        |
| `smooth_img`                            | Other          |                        |
| `subset_fields`                         | Other          |                        |
| `transducer_calibration`                | CALIBRATION    |                        |
| `transducer_positioning_with_qsub`      | POSITIONING    | Submits transducer positioning jobs to Qsub cluster using batch scripts.                        |
| `transducer_positioning_with_slurm`     | POSITIONING    | Submits transducer positioning jobs to SLURM cluster using batch scripts.                       |
| `transducer_positioning`                | POSITIONING    |                        |
| `transducer_setup`                      | Other          |                        |
| `upsample_to_grid`                      | Other          |                        |
| `zip_fields`                            | Helper         |                        |