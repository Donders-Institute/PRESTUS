## Transducer Calibration

To emulate transducers, PRESTUS aims to optimise the velocity and phase settings of virtual transducers such that they produce focal axis profiles similar to those measured with real driving system - transducer setups in free water (measured either in-house via hydrophones or provided by device manufacturers). This calibrated emulation will vary between different transducers, for different focal depth settings, and free-water pressures.

Transducer calibration relies on an additional config (`calibration_config`) that should be loaded as `parameters.calibration` (see the example `calibration_standalone`).

### Standalone calibration setup

**Script: `examples/calibration_standalone`**

#### Use cases

- Create a library of calibration settings for one (or multiple) transducer-depth settings
- Incorporate manufacturer information
- Template to set up and generate calibrated profiles for transducer equipment at the Donders Institute. 

#### Prerequisites

- Measured/estimated axial profiles (`path_input_axial`)
- Manufacturer-provided phase tables (`path_input_phase`)
- Entries in the equipment configuration (`PRESTUS/configs/equipment/equipment_config.yaml`)

For the Donders, the empirical profiles (steering tables) can be found on the [FUSInitiative OneDrive](https://radbouduniversiteit.sharepoint.com/:f:/r/sites/FUSInitiative-SHAREDResearchinformation/Shared%20Documents/Software/PRESTUS/Acoustic%20profiling/)

Requested calibrations for unique TPO-transducer & focal depth & intensity combinations can be specified in `calibration_config.yaml` via `combinations`, `focal_depths_wrt_exit_plane`, and `desired_intensities`.

Each line in the configuration corresponds to a unique TPO-transducer setup. In the following example specification, an IGT transducer with ID `PCD15287_01001` would be emulated for two focal depths of 40 and 50 mm from the exit plane, each for free-water intensities of 30 and 60 W/cm2 respectively. For a second transducer (`PCD15473_01001`), emulated phases and amplitudes would be provided for a depth of 40 mm at an intensity of 30 W/cm2. 

> This example fits a 32-channel transducer with 10 artificial channels. The number of emulated channels can impact the stability of the fitting solution. It is governed by the setup in `equipment_config.yaml`.

Transducer-TPO setups to be characterized:
```
  combinations:
    - IS_PCD15287_01001_IGT_32_ch_comb_10_ch
    - IS_PCD15473_01001_IGT_32_ch_comb_10_ch
```
List of focal depths (in mm) to be characterized:
```
  focal_depths_wrt_exit_plane:
    - [40, 50]
    - [40]
```
List of intensities (in free-water W/cm2) to be characterized:

```
  desired_intensities:
    - [30, 60]
    - [30]
```

#### Steps

- Define and initialize the simulation environment by setting paths and loading configuration files with equipment and user calibration data.
- Automatically load Donders-specific transducer information: Identify specific equipment combinations and extract associated transducer and driving system parameters, setting default initial amplitudes and phases.
- Set initial simulation to manufacturer data: Load measured characterization data of the actual transducer’s acoustic field including axial intensity profiles and phase information provided by manufacturers.
- Interpolate requested distance from multiple empirically measured distances
- Add Focal Distance Offset (FDO; via `add_FDO`): Translate measured focal depths relative to the transducer exit plane into simulation-relevant coordinates centered on the transducer’s mid-bowl. Missing distances are zero-interpolated.
- Select or interpolate the axial intensity profiles for the specified focal depth.
- Call `calibration_transducer`

###  Calibrate phase and amplitude settings

**Function: `calibration_transducer`**

#### Use cases

- Flexibility: only need to specify the desired axial profile
- Dynamic integration into end-to-end simulation loops (e.g.,  iterate across a amplitude-distance parameter space)

#### Prerequisites

- `profile_empirical.axial_intensity`   
    Desired intensity profile along focal beam axis
- `profile_empirical.axial_distance_bowl`   
    Distance (mm from transducer bowl)
- `desired_focal_distance_ep`   
    Requested focal distance (mm from transducer exit plane)

#### Steps

1. Scale the requested profile to the desired intensity
2. Run free-water simulation

    How simulations will be run will largely be determined by the main `parameters` (e.g., `parameters.code_type`). 
    However, additional settings in `parameters.calibration` can overwrite default behaviour:

    - `axisymmetric2D`  
    Overwrite default 3D simulation to perform axisymmetric 2D water simulations (`1` = yes, `0` = no (default)).
    - `force_kwavearray`    
    Force run free-water simulations with kwavearray (recommended)?  If set to `0`, simulations use the setting in the default or study-specific config.

3. Extract simulated intensity along the focal axis
4. Compute analytical O'Neil solution

    Also computes `simulated_analytical_scaling` factor from analytical to simulated intensity

5. Optimize the virtual transducer’s element velocity and phases

    Goal: match settings such that analytical profile corresponds to scaled empirical profile

    Multiple parameters configure the calibration:

    - `opt_method`  
    Optimization backend to use: `FEXminimize` (open source subtoolbox, default) | `GlobalSearch` (MATLAB's Global Optimization Toolbox)
    - `opt_weights  
    Weighting of the original profile during fitting (0 = uniform weighting, >1 increasingly narrow Gaussian FWHM). 
        - Uniform (`opt_weights == 0`): Equal weighting across entire profile—optimizes global shape.
        - Gaussian (`opt_weights >= 1`): FWHM-centered Gaussian peaking at focal maximum to emphasize near-focus optimization. Higher weights yield narrower Gaussians (sigma = focus_pos / weights).
    - `opt_limits`  
    Distance limits for optimization [mm]
    - `opt_seed`    
    Random seed for optimization. Specifying a seed increases reproducibility.
    - `opt_upper_velocity`  
    Upper velocity to use in global search.
    
    <br>

    > Note: `skip_front_peak_mm` specifies the distance to ignore from the start of axial profile (mm) to avoid near-field peak artifacts when calculating the peak distance and FWHM. It does not impact fitting. If the fit is intended over a narower range, define `opt_limits`.

    ![calibration_fitting](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/calibration_fitting.png)
    The above figure shows an example profile fit. Here, uniform `opt_weights` are used.
    <br>

6. Recalculate analytical solution with optimized phases and velocity
7. Calculate optimized source amplitude: 

    $$ \mathrm{amplitude\_optimized} = \left( \frac{\mathrm{velocity\_optimized}}{\mathrm{velocity\_original}} \right) \cdot \left( \frac{\mathrm{amplitude\_original}}{\mathrm{simulated\_analytical\_scaling}} \right) $$

8. Rerun water simulation with optimized phases and source amplitude
9. Extract simulated optimized intensity along the focal axis

    | Initial simulation | Optimized simulation |
    |--------------------|----------------------|
    | ![calibration_initial_intensity](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/calibration_initial_intensity.png) | ![calibration_opt_intensity](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/calibration_opt_intensity.png) |
    | Visualized in Step 3. | Visualized in Step 9. |

    Note: The black line indicates the position of the **transducer bowl**, the red line indicated the position of the **transducer exit plane**, the white line indicates the maximum estimated **(focal) intensity** (see [distance definitions](doc_transducer.md#target-distance-parameters)).
    <br>

10. Plot comparison between original and optimized results (analytical and simulated)

    ![calibration_optimized_analytical](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/calibration_optimized_analytical.png)