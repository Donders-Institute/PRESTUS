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

### Workflow 1: Find calibration settings for one (or multiple) transducer-depth settings [standalone]

Script: `calibration_standalone`

The above script sets up and generate calibrated profiles when using the transducer equipment at the Donders Institute. 

- Prerequisites: 
    - Manufacturer-provided phase tables
    - Measured/estimated axial profiles
    
    For the Donders, the empirical profiles (steering tables) can be found on the [FUSInitiative OneDrive](https://radbouduniversiteit.sharepoint.com/:f:/r/sites/FUSInitiative-SHAREDResearchinformation/Shared%20Documents/Software/PRESTUS/Acoustic%20profiling/)
    
- Features: 
    - automatically load Donders-specific transducer information
    - set initial simulation to manufacturer data
    - interpolate requested distance from multiple empirically measured distances
    - add Focal Distance Offset (FDO) to profile 

### Workflow 2: Dynamically update calibration settings (e.g., when iterating across settings)

Function: `calibration_transducer`

- Use cases: 
    - More limited information is available (e.g., only the axial profile is known)
    - Calibration should be performed dynamically as part of a simulation loop (e.g., for sweeping an amplitude-distance parameter space)

- Prerequisites:
    - `profile_empirical.profile_focus`        - Desired intensity profile along focal beam axis
    - `profile_empirical.dist_from_tran`       - Distance (mm from transducer bowl)
    - `profile_empirical.focus_wrt_exit_plane` - Requested focal distance (mm from transducer exit plane)

### Steps:

**Script (`examples/calibration_standalone`):**
- Define and initialize the simulation environment by setting paths and loading configuration files with equipment and user calibration data.
- Identify specific equipment combinations and extract associated transducer and driving system parameters, setting default initial amplitudes and phases.
- Load measured characterization data of the actual transducer’s acoustic field including axial intensity profiles and phase information provided by manufacturers.
- Translate experimentally measured focal depths relative to the transducer exit plane into simulation-relevant coordinates centered on the transducer’s mid-bowl to align measurement and simulation references.
- Select or interpolate the axial intensity profiles for the specified focal depth.

**Function (`calibration_transducer`):**
- Run free-water simulation.
- Optimize the virtual transducer’s element amplitudes and phases so that its simulated acoustic field matches the measured and scaled transducer profile.
- Confirm via free-water simulation with optimized virtual transducer settings.