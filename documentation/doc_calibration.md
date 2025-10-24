## Calibration of multi-element transducers

To emulate the transducers we use, we initially have to optimise the virtual transducers such that they match the output we obtain in free water. The acoustic profile will vary between different transducers, as well as for different focal depth settings for any given transducer. Additionally, the amplitude of stimulation that yields the desired output intensity needs to be estimated. As such, calibration of amplitude will also be specific to the desired free water ISPPA.

The calibration relies on additional config (`calibration_config`) that should be loaded as `parameters.calibration` (see the example script below).

#### Workflow 1: Find calibration settings for one (or multiple) transducer-depth settings [standalone]

Script: `acoustic_profiling_standalone`

The above script provides a one-click solution to generate calibrated profiles when using the transducer equipment at the Donders Institute. 

- Prerequisite: Manufacturer-provided phase tables, axial profiles
    
    For the Donders, the empirical profiles (steering tables) can be found on the [FUSInitiative OneDrive](https://radbouduniversiteit.sharepoint.com/:f:/r/sites/FUSInitiative-SHAREDResearchinformation/Shared%20Documents/Software/PRESTUS/Acoustic%20profiling/Phase_tables)
    
- Features: 
    - automatically load Donders-specific transducer information
    - set initial simulation to manufacturer data
    - interpolate requested distance from multiple empirically measured distances

#### Workflow 2: Dynamically update calibration settings (e.g., when iterating across settings)

Function: `acoustic_profiling`

- Use cases: 
    - More limited information is available (e.g., only the axial profile is known)
    - Calibration should be performed dynamically as part of a simulation loop (e.g., for sweeping an amplitude-distance parameter space)