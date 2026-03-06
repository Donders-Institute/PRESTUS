## Transducer definition

PRESTUS does not provide default transducer calibrations. Instead, it encourages users to set up transducers and optimize settings to match their empirical setup.

See the full [parameter documentation](doc_parameters.md#transducer-specification).

#### Static parameters

`transducer.source_freq_hz`  
Central frequency of the acoustic source (in Hz).

`transducer.n_elements`  
Number of transducer elements.

`transducer.Elements_ID_mm`  
Inner diameter of each transducer element (in mm).

`transducer.Elements_OD_mm`  
Outer diameter of each transducer element (in mm).

`transducer.curv_radius_mm`  
Radius of curvature of the transducer bowl (in mm).

`transducer.dist_to_plane_mm`   
Distance from the geometric focus to the transducer plane (in mm).

#### Dynamic parameters: amplitude and focus

These parameters are traditionally calibrated in free-water simulations (see below).

`transducer.source_amp`     
Amplitude of the acoustic source (in Pa). [**CALIBRATED**]

`transducer.source_phase_deg`   
Phase of the acoustic source (in degrees). [**CALIBRATED**]

`transducer.source_phase_rad`   
Phase of the acoustic source (in radians). [**CALIBRATED**]

#### Placement parameters

See also [transducer placement](doc_placement.md). In `water` simulations, the placement will be automatically determined based on the edge of the PML layer.

`transducer.trans_pos`  
Position of transducer bowl (XYZ, T1 grid voxel space).

`transducer.focus_pos`  
Position of stimulation target (XYZ, T1 grid voxel space).

#### Target distance parameters

![PRESTUS transducer distance definitions](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/transducer_distances.png)

The figure shows PRESTUS definitions for annular transducer distances.

PRESTUS expects one either of the following to be specified (or it will attempt a geometric distance calculation based on `trans_pos` and `focus_pos`, if specified).

`expected_focal_distance_ep`    
Expected distance from the transducer exit plane to the stimulation focus (mm). This often corresponds to the focal depth setting of the driving system - transducer calibration.

`expected_focal_distance_bowl`  
Expected distance from the transducer bowl to the stimulation focus (mm). Will be internally calculated based on `expected_focal_distance_ep`, `dist_to_plane_mm`, and `curv_radius_mm`.

## Transducer modeling

The `setup_source` function creates realistic transducer sources for k-Wave simulations, supporting both simple uniform elements and advanced curved arrays that produce focused beams. The function produces binary masks (active source locations) and time-varying pressure signals. Two setups are supported:

1.	Realistic kWaveArray Elements (**recommended**)
The  kWaveArray  class (`use_kWaveArray`) in the k-Wave toolbox enables realistic simulation of multi-element phased arrays by defining transducers via physical parameters (radius, element position, apodization) rather than discrete grid points. Elements render smoothly regardless of grid resolution. Independent delays/weights per element enable dynamic focusing and element overlap is automatically handled through weighted element contributions. Supports annular rings, axisymmetric 2D simulations, and GPU acceleration. This is the recommended setup type.
<br>
2.	Simple Custom Elements
Creates single- or multi-element bowl-shaped (3D) or arc-shaped (2D) pressure source regions for each transducer element using geometric functions (k-Wave function: `makeBowl`). Each element gets identical continuous wave (CW) signals. Best for uniform, non-overlapping arrays where you want simple control. The discretization of grid points can introduce  staircasing artifacts in the resulting `source.p`  masks.

## Calibration of multi-element transducers

To emulate transducers, PRESTUS aims to optimise the velocity and phase settings of virtual transducers such that they produce focal axis profiles similar to those measured with real driving system - transducer setups in free water (measured either in-house via hydrophones or provided by device manufacturers). This calibrated emulation will vary between different transducers, for different focal depth settings, and free-water pressures.

Transducer calibration relies on an additional config (`calibration_config`) that should be loaded as `parameters.calibration` (see the example `calibration_standalone`).

### Standalone calibration setup

**Script: `examples/calibration_standalone`**

#### Use cases

- Create a library of calibration settings for one (or multiple) transducer-depth settings
- Incorporate manufacturer information
- Template to set up and generate calibrated profiles for transducer equipment at the Donders Institute. 

#### Prerequisites

- Manufacturer-provided phase tables
- Measured/estimated axial profiles
    
For the Donders, the empirical profiles (steering tables) can be found on the [FUSInitiative OneDrive](https://radbouduniversiteit.sharepoint.com/:f:/r/sites/FUSInitiative-SHAREDResearchinformation/Shared%20Documents/Software/PRESTUS/Acoustic%20profiling/)
    
#### Steps

- Define and initialize the simulation environment by setting paths and loading configuration files with equipment and user calibration data.
- Automatically load Donders-specific transducer information: Identify specific equipment combinations and extract associated transducer and driving system parameters, setting default initial amplitudes and phases.
- Set initial simulation to manufacturer data: Load measured characterization data of the actual transducer’s acoustic field including axial intensity profiles and phase information provided by manufacturers.
- Interpolate requested distance from multiple empirically measured distances
- Add Focal Distance Offset (FDO; via `addEPdistance`): Translate measured focal depths relative to the transducer exit plane into simulation-relevant coordinates centered on the transducer’s mid-bowl. Missing distances are zero-interpolated.
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

    - `submit_medium`   
    Simulation submit mode: `slurm` (recommended), `matlab`, `qsub`
    - `axisymmetric2D`  
    [*Experimental*] Overwrite default 3D simulation to perform axisymmetric 2D water simulations (`1` = yes, `0` = no (default)).
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

    Note: The black line indicates the position of the **transducer bowl**, the red line indicated the position of the **transducer exit plane**, the white line indicates the maximum estimated **(focal) intensity**.
    <br>

10. Plot comparison between original and optimized results (analytical and simulated)

    ![calibration_optimized_analytical](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/calibration_optimized_analytical.png)