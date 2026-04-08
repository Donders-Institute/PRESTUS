## Transducer modeling

PRESTUS does not provide default transducer calibrations. Instead, it encourages users to set up transducers and optimize settings to match their empirical setup.

See the full [parameter documentation](doc_parameters.md#transducer-specification).

#### Static parameters

The transducer supports two alternative configurations:

- **Annular array** (`type` = `annular`)
- **Matrix array** (`type` = `matrix`)

Each type is configured in its own sub-struct (`annular` or `matrix`). Shared fields (`trans_pos`, `focus_pos`, `focal_distance_ep`/`focal_distance_bowl`, `freq_hz`) sit at the `transducer` level.

For the setup of annular transducers see [the annular definition](doc_parameters.md#annular-array-definition-annular).

For the setup of matrix transducers see [the matrix definition](doc_parameters.md#matrix-array-definition-matrix).

#### Dynamic parameters: amplitude and focus

These parameters are traditionally calibrated in free-water simulations (see below).

**Note:** Calibration is currently implemented only for **annular arrays**.

`transducer.annular.elem_amp`  
Amplitude of the acoustic source (in Pa). [**CALIBRATED**]

`transducer.annular.elem_phase_deg`  
Phase of the acoustic source per element (in degrees). [**CALIBRATED**]

For **matrix arrays**, phase delays are currently **not provided as input**.  
Instead, steering is performed automatically by computing the phase for each element based on `transducer.trans_pos` and `transducer.focus_pos`.

See [transducer calibration](doc_calibration.md).

#### Transducer placement

See also [transducer placement](doc_placement.md). 

In `water` simulations in combination with annular arrays, the placement will be automatically determined based on the edge of the PML layer.

`transducer.trans_pos`  
Position of transducer bowl (XYZ, T1 grid voxel space).
Used for placement and **phase calculation in matrix arrays**.

`transducer.focus_pos`  
Position of stimulation target (XYZ, T1 grid voxel space).
Defines the focal point and is used for **automatic steering in matrix arrays**.

#### Transducer-target distance specification

![PRESTUS transducer distance definitions](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/transducer_distances.png)

The figure shows PRESTUS definitions for annular transducer distances.

PRESTUS expects one either of the following to be specified (or it will attempt a geometric distance calculation based on `trans_pos` and `focus_pos`, if specified).

`focal_distance_ep`  
Expected distance from the transducer exit plane to the stimulation focus (mm). This often corresponds to the focal depth setting of the driving system / transducer calibration.

`focal_distance_bowl`  
Expected distance from the transducer bowl to the stimulation focus (mm). Will be internally calculated based on `focal_distance_ep`, `dist_geom_ep_mm`, and `curv_radius_mm`.

---

## Transducer geometry: annular vs. matrix array

### Annular array

An annular array consists of concentric ring-shaped elements arranged on a spherically curved surface. Because all elements share the same acoustic axis, focusing is inherently one-dimensional (axial). Lateral steering is not possible; the transducer must be physically repositioned to change the target location.

**Modeling in PRESTUS:**  
Element geometry is defined by inner/outer diameters (`Elements_ID_mm`, `Elements_OD_mm`) and bowl curvature (`curv_radius_mm`). Each element is assigned an independent pressure amplitude (`elem_amp`) and phase (`elem_phase_deg`), which must be calibrated empirically. PRESTUS includes a calibration pipeline for annular arrays to match simulated pressure to free-water measurements.

**When to use:**  
Use annular arrays when your hardware is a single-focus or multi-ring bowl transducer  with pre-calibrated per-element amplitudes and phases. This is the simpler and better-validated configuration.

---

### Matrix array

A matrix array consists of many small elements distributed across a curved or flat surface. Because each element can be driven with an independent phase delay, the acoustic focus can be steered electronically in three dimensions without moving the transducer.

**Modeling in PRESTUS:**  
Element positions are either defined analytically (Fibonacci, Fermat spiral, or rectangular grid) or loaded from a coordinate file. Element shape can be rectangular (`rect`), circular (`disc`), or curved (`bowl`). Phase delays are computed automatically from the geometric relationship between each element, `trans_pos`, and `focus_pos` — no manual phase input is required. For curved arrays, `curv_radius_mm` and `dist_geom_ep_mm` define the bowl geometry used for visualization and offset calculations. An optional Clover configuration replicates the array into a multi-aperture (up to three-leaf) arrangement.

**When to use:**  
Use matrix arrays when your hardware supports electronic steering (e.g., a phased array system) or when you want to simulate volumetric targeting without repositioning the transducer. Note that amplitude calibration is not yet implemented for matrix arrays.

---

## Transducer modeling: kWaveArray vs. simple elements

The `setup_source` function creates the acoustic source for k-Wave simulations. Two implementations are supported:

1. **Realistic kWaveArray elements** (`use_kWaveArray = 1`, recommended)  
   Uses the k-Wave `kWaveArray` class to define transducer elements via physical parameters (position, radius, apodization) rather than discrete grid voxels. Elements render smoothly at any grid resolution. Independent per-element delays and weights enable dynamic focusing. Element overlap is handled automatically through weighted contributions. Supports annular rings, axisymmetric 2D simulations, matrix arrays, and GPU acceleration.

2. **Simple custom elements** (`use_kWaveArray = 0`)  
   Creates bowl-shaped (3D) or arc-shaped (2D) pressure source regions per element using the k-Wave `makeBowl` function. Each element receives an identical continuous wave (CW) signal. Best for uniform single-element setups where simplicity is preferred. Grid discretization can introduce staircasing artifacts in the source mask.
