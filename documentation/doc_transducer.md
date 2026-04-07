## Transducer definition

PRESTUS does not provide default transducer calibrations. Instead, it encourages users to set up transducers and optimize settings to match their empirical setup.

See the full [parameter documentation](doc_parameters.md#transducer-specification).

#### Static parameters

The transducer supports two configurations:

- **Annular array** (`type` = `annular`)
- **Matrix array** (`type` = `matrix`)

Each type is configured in its own sub-struct (`annular` or `matrix`). 

For the setup of annular transducers see [the annular definition](doc_parameters.md#annular-array-definition-annular).

For the setup of matrix transducers see [the matrix definition](doc_parameters.md#matrix-array-definition-matrix).

Shared fields such as `trans_pos`, `focus_pos`, and `exp_FD_ep`/`exp_FD_bowl` sit at the `transducer` level.

#### Dynamic parameters: amplitude and focus

These parameters are traditionally calibrated in free-water simulations (see below).

**Note:** Calibration is currently implemented only for **annular arrays**.

`transducer.annular.source_freq_hz`  
Central frequency of the acoustic source (in Hz).

`transducer.annular.source_amp`  
Amplitude of the acoustic source (in Pa). [**CALIBRATED**]

`transducer.annular.source_phase_deg`  
Phase of the acoustic source per element (in degrees). [**CALIBRATED**]

`transducer.annular.source_phase_rad`  
Phase of the acoustic source per element (in radians). [**CALIBRATED**]

For **matrix arrays**, phase delays are currently **not provided as input**.  
Instead, steering is performed automatically by computing the phase for each element based on `transducer.position.trans_pos` and `transducer.position.focus_pos`.

See [transducer calibration](doc_calibration.md).

#### Transducer placement

See also [transducer placement](doc_placement.md). 

In `water` simulations in combination with annular arrays, the placement will be automatically determined based on the edge of the PML layer.

`transducer.position.trans_pos`  
Position of transducer bowl (XYZ, T1 grid voxel space).
Used for placement and **phase calculation in matrix arrays**.

`transducer.position.focus_pos`  
Position of stimulation target (XYZ, T1 grid voxel space).
Defines the focal point and is used for **automatic steering in matrix arrays**.

#### Target distance parameters

![PRESTUS transducer distance definitions](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/transducer_distances.png)

The figure shows PRESTUS definitions for annular transducer distances.

PRESTUS expects one either of the following to be specified (or it will attempt a geometric distance calculation based on `position.trans_pos` and `position.focus_pos`, if specified).

`position.exp_FD_ep`    
Expected distance from the transducer exit plane to the stimulation focus (mm). This often corresponds to the focal depth setting of the driving system - transducer calibration.

`position.exp_FD_bowl`  
Expected distance from the transducer bowl to the stimulation focus (mm). Will be internally calculated based on `exp_FD_ep`, `dist_to_plane_mm`, and `curv_radius_mm`.

## Transducer modeling

The `setup_source` function creates realistic transducer sources for k-Wave simulations, supporting both simple uniform elements and advanced curved arrays that produce focused beams. The function produces binary masks (active source locations) and time-varying pressure signals. Two setups are supported:

1.	Realistic kWaveArray Elements (**recommended**)
The  kWaveArray  class (`use_kWaveArray`) in the k-Wave toolbox enables realistic simulation of multi-element phased arrays by defining transducers via physical parameters (radius, element position, apodization) rather than discrete grid points. Elements render smoothly regardless of grid resolution. Independent delays/weights per element enable dynamic focusing and element overlap is automatically handled through weighted element contributions. Supports annular rings, axisymmetric 2D simulations, matrix arrays, and GPU acceleration. This is the recommended setup type.
<br>
2.	Simple Custom Elements
Creates single- or multi-element bowl-shaped (3D) or arc-shaped (2D) pressure source regions for each transducer element using geometric functions (k-Wave function: `makeBowl`). Each element gets identical continuous wave (CW) signals. Best for uniform, non-overlapping arrays where you want simple control. The discretization of grid points can introduce  staircasing artifacts in the resulting `source.p`  masks.

