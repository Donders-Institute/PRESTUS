# Simulation outputs

PRESTUS provides multiple outputs. Some of these outputs are optional (or can be deactivated upon request to save space.)

#### Summary table

Filename: sub-XXX_<simulation_medium>_output_table<affix>.csv

This output table provides an overview of key metrics.

**Acoustic simulation**

| Parameter                     | Description                                         |
|-------------------------------|-----------------------------------------------------|
| subject_id                    | Unique identifier for the subject                   |
| max_Isppa                     | Maximum intensity [W/cm²]                           |
| max_Isppa_after_exit_plane    | Maximum intensity after the exit plane [W/cm²]      |
| real_focal_distance           | Empirical focal distance (distance between transducer and max. intensity in brain medium (if modelled) or in any medium beyond the exit plane) [mm] |
| max_Isppa_skin                | Maximum intensity (skin medium) [W/cm²]             |
| max_Isppa_skull               | Maximum intensity (skull medium) [W/cm²]            |
| max_Isppa_brain               | Maximum intensity (brain medium) [W/cm²]            |
| max_pressure_skin             | Maximum pressure (skin medium) [Pa]                 |
| max_pressure_skull            | Maximum pressure (skull medium) [Pa]                |
| max_pressure_brain            | Maximum pressure (brain medium) [Pa]                |
| max_MI_skin                   | Maximum Mechanical Index (skin medium)              |
| max_MI_skull                  | Maximum Mechanical Index (skull medium)             |
| max_MI_brain                  | Maximum Mechanical Index (brain medium)             |
| Ix_brain                      | X-coordinate of maximum intensity (brain medium)    |
| Iy_brain                      | Y-coordinate of maximum intensity (brain medium)    |
| Iz_brain                      | Z-coordinate of maximum intensity (brain medium)    |
| trans_pos_final_1             | X-coordinate of final transducer position           |
| trans_pos_final_2             | Y-coordinate of final transducer position           |
| trans_pos_final_3             | Z-coordinate of final transducer position           |
| focus_pos_final_1             | X-coordinate of final focus position                |
| focus_pos_final_2             | Y-coordinate of final focus position                |
| focus_pos_final_3             | Z-coordinate of final focus position                |
| isppa_at_target               | Intensity value at the target location [W/cm²]      |
| avg_isppa_around_target       | Average intensity around the target [W/cm²]         |
| half_max_ISPPA_volume_brain   | Volume of brain with half maximum intensity         |

**Thermal simulation**

| Parameter                     | Description                                         |
|-------------------------------|-----------------------------------------------------|
| maxT                          | Maximum temperature (global medium) [°C]            |
| maxCEM43                      | Maximum CEM43 thermal dose (global medium)          |
| maxT_brain                    | Maximum temperature (brain medium) [°C]             |
| maxT_skull                    | Maximum temperature (skull medium) [°C]             |
| maxT_skin                     | Maximum temperature (skin medium) [°C]              |
| riseT_brain                   | Temperature rise (brain medium) [°C]                |
| riseT_skull                   | Temperature rise (skull medium) [°C]                |
| riseT_skin                    | Temperature rise (skin medium) [°C]                 |
| CEM43_brain                   | CEM43 Thermal dose (brain medium) [mins]            |
| CEM43_skull                   | CEM43 Thermal dose (skull medium) [mins]            |
| CEM43_skin                    | CEM43 Thermal dose (skin medium) [mins]             |

#### 2D/3D NIFTI images

PRESTUS outputs 2D or 3D images (depending on the input) of the following metrics, which can be used to calculat additional metrics post-hoc. 
Images are provided in subject-space (```_orig_coord_```) and in MNI-152 space (```_MNI```).

- medium_masks
- isppa [acoustic]
- MI [acoustic]
- pressure [acoustic]
- heating [heating]
- heatrise [heating]
- CEM43 [heating]

#### MATLAB structures

PRESTUS saves an overview of the parameters used to run the simulation and by default saves structures from acoustic and heating simulations. If these are detected in the results folder (and overwriting is deactivated), they will be loaded instead of performing the calculation.

- sub-XXX_<simulation_medium>_parameters<affix>.mat     | simulation parameters
- sub-XXX_<simulation_medium>_kwave_source<affix>.mat   | k-Wave source | parameters, kgrid, trans_pos_final, focus_pos_final
- sub-XXX_<simulation_medium>_results<affix>.mat        | acoustic simulation outputs [can be deactivated via ```savemat``` = 0]
- sub-XXX_<simulation_medium>_heating_res<affix>.mat    | thermal simulation outputs [can be deactivated via ```savemat``` = 0]

#### Figures

PRESTUS provides multiple figures for quick visual inspection and debugging.

TBD