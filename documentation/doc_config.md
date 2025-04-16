## PRESTUS configuration documentation

The following documents the parameters used in PRESTUS. The specification philosophy is the following: `default_config.yaml` provides a list of all parameters with default settings and will be read in first by the function `load_parameters`. This default configuration file should not be changed in standard applications to ensure that necessary fields are provided. 

To set up a specific application, an additional `config_xxx.yaml` should be provided. This file should contain exclusively the fields where defaults should be overwritten (e.g., to provide specific transducer settings). Alternatively, parameters can be specified prior to calling the `single_subject_pipeline` in MATLAB. This allows dynamic iterations over parameters of interest.

### I/O management

| **Parameter**                     | **Description**                                                                                                      | **Comments**                                                                 |
|-----------------------------------|----------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------|
| `data_path`                       | Absolute path to input data.                          | [string] Mandatory |
| `seg_path`                        | Absolute path to SimNIBS segmentations.               | [string] Mandatory |
| `data_path`                       | Absolute path to the data location.                   | [string] Mandatory |
| `sim_path`                        | Absolute path to the simulation output.               | [string] Mandatory |
| `simnibs_bin_path`                | Absolute path to SimNIBS binaries.                    | [string] Mandatory |
| `paths_to_add`                    | Toolbox paths to add with addpath().                  | [cell] e.g., `{"path/to/x", "path/to/Y"}` |
| `subpaths_to_add`                 | Toolbox paths to add with addpath(genpath()).         | [cell] e.g., `{"path/to/x", "path/to/Y"}` |
| `subject_subfolder`               | Manage outputs in subject-specific subdirectories?    | (`1 = yes [default], 0 = no`) |
| `results_filename_affix`          | Affix for result file names                           | [string] Can be used to differentiate simulation outputs for the same subject-transducer combination(e.g., different intensities and/or targets.)  |
| `overwrite_simnibs`               | Overwrite SimNIBS segmentation results?               | (`1 = yes, 0 = no` [default]) |

### Simulation type

| **Parameter**                     | **Description**                                                                                                      | **Comments**                                                                  |
|-----------------------------------|----------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------|
| `simulation_medium`               | Medium setup for simulation (`water` or `layered`).                                   | Mandatory.   |
| `layer_labels`                    | Labels for layered simulation, defining mask indices for different tissue types.      | Mandatory.   |
| `seg_labels`                      | Labels for segmentations, specifying indices for CSF, bone mask, and eye regions.     | Mandatory.  |
| `run_source_setup`                | Set up acoustic sources.                                                              | (`1 = yes, 0 = no`)  |
| `run_acoustic_sims`               | Run acoustic simulations.                                                             | (`1 = yes, 0 = no`)  |
| `run_heating_sims`                | Run heating simulations.                                                              | (`1 = yes, 0 = no`)  |
| `run_posthoc_water_sims`          | Run water simulations following head simulations.                                     | (`1 = yes, 0 = no`)  |
| `interactive`                     | Interactive mode (`1 = yes, 0 = no`).                                                 | (`1 = yes, 0 = no`) Asks prior to overwriting or starting long computations. |
| `overwrite_files`                 | File overwrite behavior (`ask`, `never`, or `always`).                                | (`ask`, `never`, or `always`) This parameter does NOT apply to SimNIBS segmentations.  |
| `use_kWaveArray`                  | Use the kWaveArray class for simulations.                                             | (`1 = yes, 0 = no`) see k-Wave documentation.    |
| `savemat`                         | Save outputs of acoustic and/or heating simulations as .mat files?                    | (`1 = yes, 0 = no`) For many parallel simulations, setting this to 0 saves disk space.    |


### Segmentation/Preprocessing

| **Parameter**                     | **Description**                                                                                                      | **Comments**                                                                  |
|-----------------------------------|----------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------|
| `segmentation_software`           | Segmentation software used (`headreco` or `charm`).                                                       | Use of the former `headreco` may result in errors (e.g., when creating pseudoCTs) due to the assumption of charm-based tissue labels in parts of the codebase.    |
| `csf_mask_expansion_factor`       | Expansion factor for cerebrospinal fluid (CSF) brain mask; controls mask dilation based on voxel size.    |   |
| `skull_smooth_threshold`          | Threshold for smoothing the skull mask; higher values result in thinner masks.                            |   |
| `other_smooth_threshold`          | Threshold for smoothing other masks; higher values result in thinner masks.                               |   |


### Transducer modeling

| **Parameter**                     | **Description**                                                                                                      | **Comments**                                                                  |
|-----------------------------------|----------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------|
| `transducer.source_freq_hz`       | Central frequency of the acoustic source (in Hz).                     |    |
| `transducer.n_elements`           | Number of elements in the transducer.                                 |    |
| `transducer.Elements_ID_mm`       | Inner diameter of each transducer element (in mm).                    |    |
| `transducer.Elements_OD_mm`       | Outer diameter of each transducer element (in mm).                    |    |
| `transducer.curv_radius_mm`       | Radius of curvature of the transducer bowl (in mm).                   |    |
| `transducer.dist_to_plane_mm`     | Distance from the geometric focus to the transducer plane (in mm).    |    |
| `transducer.source_amp`           | Amplitude of the acoustic source (in Pa).                             | Must be calibrated.   |
| `transducer.source_phase_deg`     | Phase of the acoustic source (in degrees).                            | Must be calibrated.   |
| `expected_focal_distance_mm`      | Expected distance to the stimulation focus (in mm).                   | Transducer depth setting |
| `transducer_from_localite`        | Load transducer position from Localite files?.                        | (`1 = yes, 0 = no` [default])   |
| `reference_transducer_distance_mm`  | Distance from tracker to transducer exit plane (in mm).             | Allows to correct for varying distances between the infrared trackers attached to the transducer and the exit plane. Only applies when `transducer_from_localite=1`. |


### Simulation grid

| **Parameter**                     | **Description**                                                                                                      | **Comments**                                                                  |
|-----------------------------------|----------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------|
| `grid_step_mm`                    | Resolution of the computational grid (must be isotropic, in mm).                                                     | |
| `default_grid_size`               | Default size of the simulation grid per dimension (number of points).                                                | |
| `default_grid_dims`               | Default dimensions of the simulation grid `[Nx, Ny, Nz]`.                                                            | |
| `focus_pos_t1_grid`               | Stimulation target position on T1 grid space.                                                                        | |
| `focus_area_radius`               | Radius of the target area around the focus where ISPPA is averaged (in mm).                                          | |
| `pml_size`                        | Size of the Perfectly Matched Layer (PML) used to absorb waves at the grid boundaries (default is 10 for 3D grids).  | see k-Wave documentation.|
| `prime_factor_max_grid_expansion` | Maximum expansion factor for computational grid to optimize prime numbers and speed up computations.                 | |


### Medium properties

| **Parameter**                     | **Description**                                                                                                      | **Comments**                                                                  |
|-----------------------------------|----------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------|
| `medium.water.sound_speed`                        | Speed of sound in water (in m/s).                                                 | ITRUSST benchmarks.  |
| `medium.water.density`                            | Density of water (in kg/m³).                                                      | Tissue Properties DB.|
| `medium.water.alpha_0_true`                       | Attenuation coefficient for water (in dB/cm/MHz).                                 | k-Plan documentation.|
| `medium.water.alpha_power_true`                   | Exponent for attenuation coefficient power law for water.                         | Tissue Properties DB.|
| `medium.skull.sound_speed`                        | Speed of sound in skull bone (in m/s).                                            | ITRUSST benchmarks.  |
| `medium.skull.density`                            | Density of skull bone (in kg/m³).                                                 | ITRUSST benchmarks.  |
| `medium.skull.alpha_0_true`                       | Attenuation coefficient at for skull bone (in dB/cm/MHz).                         | Pinton et al., 2011. |
| `medium.skull.alpha_power_true`                   | Exponent for attenuation coefficient power law for skull bone.                    | k-Plan documentation.|
| `medium.skull.thermal_conductivity`               | Thermal conductivity of skull bone (in W/m/°C).                                   | Tissue Properties DB.| 
| `medium.skull.specific_heat_capacity`             | Specific heat capacity of skull bone (in J/kg/°C).                                | Tissue Properties DB.| 
| `medium.brain.sound_speed`                        | Speed of sound in brain tissue (in m/s).                                          | Tissue Properties DB.|
| `medium.brain.density`                            | Density of brain tissue (in kg/m³).                                               | Tissue Properties DB.|
| `medium.brain.alpha_0_true`                       | Attenuation coefficient at for brain tissue (in dB/cm/MHz).                       | k-Plan documentation.|
| `medium.brain.alpha_power_true`                   | Exponent for attenuation coefficient power law for brain tissue.                  | k-Plan documentation.|
| `medium.brain.thermal_conductivity`               | Thermal conductivity of brain tissue (in W/m/°C).                                 | Tissue Properties DB.| 
| `medium.brain.specific_heat_capacity`             | Specific heat capacity of brain tissue (in J/kg/°C).                              | Tissue Properties DB.| 
| `medium.skin.sound_speed`                         | Speed of sound in skin tissue (in m/s).                                           | ITRUSST benchmarks.  |
| `medium.skin.density`                             | Density of skin tissue (in kg/m³).                                                | ITRUSST benchmarks.  |
| `medium.skin.alpha_0_true`                        | Attenuation coefficient for skin tissue (in dB/cm/MHz).                           | ITRUSST benchmarks.  |
| `medium.skin.alpha_power_true`                    | Exponent for attenuation coefficient power law for skin tissue.                   | ITRUSST benchmarks.  |
| `medium.skin.thermal_conductivity`                | Thermal conductivity of skin tissue (in W/m/°C).                                  | Tissue Properties DB.| 
| `medium.skin.specific_heat_capacity`              | Specific heat capacity of skin tissue (in J/kg/°C).                               | Tissue Properties DB.| 
| `medium.skull_trabecular.sound_speed`             | Speed of sound in trabecular bone (in m/s).                                       | ITRUSST benchmarks.  |
| `medium.skull_trabecular.density`                 | Density of trabecular bone (in kg/m³).                                            | ITRUSST benchmarks.  |
| `medium.skull_trabecular.alpha_0_true`            | Attenuation coefficient for trabecular bone (in dB/cm/MHz).                       | Pinton et al., 2011. |
| `medium.skull_trabecular.alpha_power_true`        | Exponent for attenuation coefficient power law for trabecular bone.               | k-Plan documentation.|
| `medium.skull_trabecular.thermal_conductivity`    | Thermal conductivity of trabecular bone (in W/m/°C).                              | Tissue Properties DB.| 
| `medium.skull_trabecular.specific_heat_capacity`  | Specific heat capacity of trabecular bone (in J/kg/°C).                           | Tissue Properties DB.| 
| `medium.skull_cortical.sound_speed`               | Speed of sound in cortical bone (in m/s).                                         | ITRUSST benchmarks.  |
| `medium.skull_cortical.density`                   | Density of cortical bone (in kg/m³).                                              | ITRUSST benchmarks.  |
| `medium.skull_cortical.alpha_0_true`              | Attenuation coefficient for cortical bone (in dB/cm/MHz).                         | Pinton et al., 2011. |
| `medium.skull_cortical.alpha_power_true`          | Exponent for attenuation coefficient power law for cortical bone.                 | k-Plan documentation.|
| `medium.skull_cortical.thermal_conductivity`      | Thermal conductivity of cortical bone (in W/m/°C).                                | Tissue Properties DB.| 
| `medium.skull_cortical.specific_heat_capacity`    | Specific heat capacity of cortical bone (in J/kg/°C).                             | Tissue Properties DB.| 


### pseudoCT mapping to skull properties
see doc_pseudoCT.md

| **Parameter**                     | **Description**                                                                                                      | **Comments**                                                                  |
|-----------------------------------|----------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------|
| `usepseudoCT`                     | Use (pseudo-)CT based mapping?  (`1 = yes, 0 = no`)               |  Maps (pseudo-)HU to tissue properties in the trabecular and cortical skull layer. |
| `pseudoCT_variant`                | Mapping algorithm (`yaakub`/`carpino`/`k-plan`/`marquet`)         | see doc_pseudoCT.md  |


### Sequence timing and baseline temperature for heating simulations
see doc_thermal-simulations.md

| **Parameter**                     | **Description**                                                                                                      | **Comments**                                                                  |
|-----------------------------------|----------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------|
| `thermal.duty_cycle`              | Fraction of stimulation duration during which the stimulation is active (`0 to 1`).                           |   |
| `thermal.iti`                     | Interval between trials, measured from the start of one trial to the start of another (in seconds).           |   |
| `thermal.n_trials`                | Number of trials simulated, determining total simulation duration (`n_trials * iti`).                         |   |
| `thermal.stim_duration`           | Duration of stimulation within a trial (in seconds).                                                          |   |
| `thermal.sim_time_steps`          | Simulation time steps during the stimulation period (in seconds).                                             |   |
| `thermal.post_stim_time_step_dur` | Time step duration during post-stimulation period (inter-trial interval, in seconds).                         |   |
| `thermal.on_off_step_duration`    | Duration of the on/off cycle during stimulation, computed based on duty cycle and time steps (in seconds).    |   |
| `thermal.equal_steps`             | Whether simulation step durations are equal for on and off cycles (`1 = yes, 0 = no`).                        |   |
| `thermal.temp_0.water`            | Initial temperature of water medium before simulation (in °C).                                                |   |
| `thermal.temp_0.skull`            | Initial temperature of skull medium before simulation (in °C).                                                |   |
| `thermal.temp_0.brain`            | Initial temperature of brain medium before simulation (in °C).                                                |   |
| `thermal.temp_0.skin`             | Initial temperature of skin medium before simulation (in °C).                                                 |   |
| `thermal.temp_0.skull_trabecular` | Initial temperature of trabecular skull medium before simulation (in °C).                                     |   |
| `thermal.temp_0.skull_cortical`   | Initial temperature of cortical skull medium before simulation (in °C).                                       |   |
| `thermal.sensor_xy_halfsize`      | Maximum size of the sensor window for temperature recording (in grid units).                                  |   |
| `thermal.record_t_at_every_step`  | Whether to record temperature at every time step for the whole sensor window (`1 = yes, 0 = no`).             |   |
| `thermal.record_t_at_every_step`  | Whether to record temperature at every time step for the whole sensor window (`1 = yes, 0 = no`).             |   |
  `thermal.continuous_protocol`     | Model a continuous pulsed protocol?                                                                           |   |
| `heatingvideo`                    | Save a video of incremental heating? (`1 = yes, 0 = no`)                                                      |   |


### GPU/HPC options
see doc_hpc.md

| **Parameter**                     | **Description**                                                                                                      | **Comments**                                                                  |
|-----------------------------------|----------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------|
| `code_type`                       | Type of k-Wave code to run (`matlab_cpu`, `matlab_gpu`, `cpp_interactive`, `cpp_noninteractive`, or `cuda`).         |   |
| `using_donders_hpc`               | Run simulations on Donders HPC (`1 = yes, 0 = no`).                                                                  |(`1 = yes, 0 = no`)   |
| `torque.ld_library_path`          | Path to LD_LIBRARY used during SimNIBS installation on Torque systems.                                     | [Optional] If you experience an `undefined symbol` error in `create_mesh_surf.cpython-39-x86_64-linux-gnu.so` define the LD_LIBRARY location (e.g.,: ```/opt/gcc/7.2.0/lib64```) |
| `slurm.ld_library_path`           | Path to LD_LIBRARY used during SimNIBS installation on Slurm systems. [Optional]                                     |  e.g.,: ```/home/'group'/'user'/.conda/envs/simnibs_env/lib/python3.9/site-packages/simnibs/mesh_tools/cgal/../../external/lib/linux``` |
| `hcp_gpu`                         | Request a specific GPU. [Optional]                                                                                   |  Not recommended by default, rely on automatic GPU detection instead. May be useful when benchmarking specific GPUs. E.g.,```"nvidia_a100-sxm4-40gb:1"```. ```scontrol show nodes \| egrep -o gres/gpu:.*=[0-9] \| egrep -o 'nvidia_.*=' \| sort \| uniq \| sed 's/=//'``` lists available GPU types. |
| `hcp_partition`                   | Request a dedicated GPU partition. [Optional]                                                                        | The Donders HCP provides a ```gpu40g``` partition that consists of nodes with GPU with vRAM > 40 GB. This is the recommended default for thermal simulations of longer protocols. |
| `hpc_reservation`                  | Request a reserved cue. [Optional]                                                                                  | |