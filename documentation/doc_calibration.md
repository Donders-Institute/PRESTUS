## Calibration of multi-element transducers

To emulate the actual transducers for the use in simulations, we initially have to optimise the virtual transducers such that they match the output we obtain in free water. The acoustic profile will vary between different transducers, as well as for different focal depth settings for any given transducer. Additionally, the amplitude of stimulation that yields the desired output intensity needs to be estimated. As such, calibration of amplitude will also be specific to the desired free water ISPPA.

The calibration relies on additional config (`calibration_config`) that should be loaded as `parameters.calibration` (see the example script below).

### Workflow 1: Find calibration settings for one (or multiple) transducer-depth settings [standalone]

Script: `acoustic_profiling_standalone`

The above script provides a one-click solution to generate calibrated profiles when using the transducer equipment at the Donders Institute. 

- Prerequisites: 
    - Manufacturer-provided phase tables
    - Measured/estimated axial profiles
    
    For the Donders, the empirical profiles (steering tables) can be found on the [FUSInitiative OneDrive](https://radbouduniversiteit.sharepoint.com/:f:/r/sites/FUSInitiative-SHAREDResearchinformation/Shared%20Documents/Software/PRESTUS/Acoustic%20profiling/)
    
- Features: 
    - automatically load Donders-specific transducer information
    - set initial simulation to manufacturer data
    - interpolate requested distance from multiple empirically measured distances

### Workflow 2: Dynamically update calibration settings (e.g., when iterating across settings)

Function: `acoustic_profiling`

- Use cases: 
    - More limited information is available (e.g., only the axial profile is known)
    - Calibration should be performed dynamically as part of a simulation loop (e.g., for sweeping an amplitude-distance parameter space)

- Prerequisites:
    - `profile_empirical.profile_focus`        - Measured or theoretical intensity profile along beam axis
    - `profile_empirical.dist_from_tran`       - Distance (mm from transducer reference point)
    - `profile_empirical.focus_wrt_exit_plane` - Corresponding focal distance (mm from transducer exit plane)

### Steps:

**Script (`acoustic_profiling_standalone`):**
- Define and initialize the simulation environment by setting paths and loading configuration files with equipment and user calibration data.
- Identify specific equipment combinations and extract associated transducer and driving system parameters, setting default initial amplitudes and phases.
- Load measured characterization data of the actual transducer’s acoustic field including axial intensity profiles and phase information provided by manufacturers.
- Translate experimentally measured focal depths relative to the transducer exit plane into simulation-relevant coordinates centered on the transducer’s mid-bowl to align measurement and simulation references.
- Select or interpolate the axial intensity profiles for the specified focal depth.

**Function (`acoustic_profiling`):**
- Run free-water simulation.
- Optimize the virtual transducer’s element amplitudes and phases so that its simulated acoustic field matches the measured and scaled transducer profile.
- Confirm via free-water simulation with optimized virtual transducer settings.