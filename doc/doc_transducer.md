## Transducer modeling

!!! note "At a glance"
    - **What:** The fields that define a transducer: type, frequency, element geometry, drive, phases and position.
    - **When you need it:** When your device is not in the [transducer library](doc_transducer_library.md), or to understand the fields the library fills in. With `transducer.serial` you rarely set these yourself.
    - **Before you start:** The geometry of your device (from the manufacturer).

A transducer can be taken from the shipped [transducer library](doc_transducer_library.md) (`transducer.serial`: geometry, phases and drive are filled in) or defined in full in the study config. This page describes the fields; the library and calibration pages describe where calibrated values come from.

See the full [parameter documentation](doc_parameters.md#transducer-specification-transducer).

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

`transducer.annular.elem_p_amp`  
Amplitude of the acoustic source (in Pa). [**CALIBRATED**]

`transducer.annular.elem_phase_deg`  
Phase of the acoustic source per element (in degrees). [**CALIBRATED**]

For **matrix arrays**, steering is by default performed automatically by computing the phase for each element based on `transducer.trans_pos` and `transducer.focus_pos`. Calibrated phases (one per element) can replace the geometric phases via `transducer.matrix.elem_phase_override_deg`.

See [transducer calibration](doc_calibration.md).

#### Transducer placement

See also [transducer placement](doc_placement.md). 

In `water` simulations in combination with annular arrays, the placement will be automatically determined based on the edge of the PML layer.

`transducer.trans_pos`  
Position of transducer bowl (XYZ; T1 voxels by default, or RAS/MNI mm via `placement.manual.coordinate_space`).
Used for placement and **phase calculation in matrix arrays**.

`transducer.focus_pos`  
Position of stimulation target (XYZ; same coordinate space as `trans_pos`).
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
Element geometry is defined by inner/outer diameters (`elem_id_mm`, `elem_od_mm`) and bowl curvature (`curv_radius_mm`). Each element is assigned an independent pressure amplitude (`elem_p_amp`) and phase (`elem_phase_deg`), which must be calibrated empirically. PRESTUS includes a calibration pipeline for annular arrays to match simulated pressure to free-water measurements.

**When to use:**  
Use annular arrays when your hardware is a single-focus or multi-ring bowl transducer  with pre-calibrated per-element amplitudes and phases. This is the simpler and better-validated configuration.

#### Flat elements

`curv_radius_mm: .inf` makes the elements flat, e.g. a plane piston (the ITRUSST benchmark's 20 mm source). The face lies at `trans_pos`, normal to the `trans_pos` → `focus_pos` direction, and is also the exit plane (bowl depth 0); `focus_pos` only sets the beam direction. Supported with the kWaveArray source (`grid.use_kWaveArray: 1`):

| Grid | Element | Rings (`elem_id_mm` > 0) |
|---|---|---|
| 2D axisymmetric | line across the axis (`addLineElement`) | yes |
| 2D | line normal to the beam | no |
| 3D | disc (`addDiscElement`) | no (single disc) |

Not supported: the makeBowl source (`use_kWaveArray: 0`), BabelBrain, transducer calibration and the placement planners, which assume a bowl.

---

### Matrix array

A matrix array consists of many small elements distributed across a curved or flat surface. Because each element can be driven with an independent phase delay, the acoustic focus can be steered electronically in three dimensions without moving the transducer.

**Modeling in PRESTUS:**  
Element positions are either defined analytically (Fibonacci, Fermat spiral, or rectangular grid) or loaded from a coordinate file. Element shape can be rectangular (`rect`), circular (`disc`), or curved (`bowl`). Phase delays are computed automatically from the geometric relationship between each element, `trans_pos`, and `focus_pos` — no manual phase input is required; `elem_phase_override_deg` replaces them with calibrated phases. The aperture of `disc`/`bowl` elements is equal-area to `elem_width_mm` × `elem_height_mm` unless `elem_diameter_mm` is set. For curved arrays, `curv_radius_mm` and `dist_geom_ep_mm` define the bowl geometry used for visualization and offset calculations. An optional Clover configuration replicates the array into a multi-aperture arrangement of `n_leaves` leaves on a parent bowl (`ROC_parent`), `az_sep_deg` apart (default `360 / n_leaves`) and tilted by `elevation_deg` (default: the smallest non-overlapping tilt).

**When to use:**  
Use matrix arrays when your hardware supports electronic steering (e.g., a phased array system) or when you want to simulate volumetric targeting without repositioning the transducer. Note that amplitude calibration is not yet implemented for matrix arrays.

<a id="beam-axis-constraint"></a>**Beam axis constraint:**  
Element orientation and phase delay computation (`create_matrix_karray`) assumes the transducer's natural beam axis is aligned with the grid +Z direction. This holds for the standard single-transducer placement (transducer at one face of the grid, firing toward the target). In multi-transducer configurations where a matrix transducer fires along a different axis (e.g. a second transducer opposing the first along −Z), element normals and phase delays will be computed incorrectly. Annular transducers are not affected by this constraint. If a matrix bilateral setup is required, the beam direction would need to be derived from `(focus_pos − trans_pos)` rather than hardcoded to `[0, 0, 1]` in `create_matrix_karray`.

---

## Transducer modeling: kWaveArray vs. simple elements

The `source_create` function creates the acoustic source for k-Wave simulations. Two implementations are supported:

1. **Realistic kWaveArray elements** (`use_kWaveArray = 1`, recommended)  
   Uses the k-Wave `kWaveArray` class to define transducer elements via physical parameters (position, radius, apodization) rather than discrete grid voxels. Elements render smoothly at any grid resolution. Independent per-element delays and weights enable dynamic focusing. Element overlap is handled automatically through weighted contributions. Supports annular rings, axisymmetric 2D simulations, matrix arrays, and GPU acceleration.

2. **Simple custom elements** (`use_kWaveArray = 0`)  
   Creates bowl-shaped (3D) or arc-shaped (2D) pressure source regions per element using the k-Wave `makeBowl` function. Each element receives an identical continuous wave (CW) signal. Best for uniform single-element setups where simplicity is preferred. Grid discretization can introduce staircasing artifacts in the source mask.
