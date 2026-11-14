# Transducer Definition

The transducer definition sets the source of the simulated pressure field: its frequency, element geometry, drive amplitude and phases, and its position relative to the target. The simulated focus position, focal size and intensity follow from these fields. For a device in the [transducer library](doc_transducer_library.md), `transducer.serial` fills them in; this page explains each field and how to define a device the library does not hold.

On this page:

- [Transducer Fields](#transducer-fields): the fields of a transducer and a full inline definition.
- [Amplitude and Phases](#amplitude-and-phases): the drive, set by the focal depth setting and target ISPPA, and the calibrated phases and amplitudes behind it.
- [Transducer and Focus Positions](#transducer-and-focus-positions): the transducer and focus positions.
- [Transducer-Target Distance](#transducer-target-distance): exit-plane and bowl distances to the focus.
- [Annular Array](#annular-array): ring geometry, including flat elements.
- [Matrix Array](#matrix-array): element layouts, steering and its beam-axis constraint.

!!! note "At a glance"
    - **When you need it:** When your device is not in the [transducer library](doc_transducer_library.md), or to understand the fields the library fills in. With `transducer.serial` you rarely set these yourself.
    - **Before you start:** The geometry of your device (from the manufacturer).

## Transducer Fields

A transducer comes from the shipped [transducer library](doc_transducer_library.md) (`transducer.serial`: geometry, phases and drive are filled in) or is defined in full in the study config. Every field is listed in the [parameter documentation](doc_parameters.md#transducer-specification-transducer).

A full definition of an annular array, bypassing the library:

```yaml
transducer:
  serial: IS_PCD15473_01001
  type: annular
  freq_hz: 300000.0
  focal_distance_ep: 60.0
  annular:
    elem_n: 10
    elem_id_mm: [10.0, 22.1, ...]
    elem_od_mm: [21.1, 28.8, ...]
    curv_radius_mm: 100.0
    dist_geom_ep_mm: 92.7
```

The fields fall into three groups:

| Group | Fields | Where described |
|---|---|---|
| Type and geometry | `type` (`annular` or `matrix`), `freq_hz`, and the sub-struct of that type (`annular` or `matrix`) | [Annular definition](doc_parameters.md#annular-array-definition-annular), [matrix definition](doc_parameters.md#matrix-array-definition-matrix) |
| Drive and steering | `focal_distance_ep` and `target_isppa_wcm2`; calibrated `annular.elem_phase_deg` and `annular.elem_p_amp` (advanced); matrix arrays steer automatically | [Amplitude and phases](#amplitude-and-phases) |
| Position and target | `trans_pos`, `focus_pos`, `focal_distance_ep` or `focal_distance_bowl` | [Transducer and focus positions](#transducer-and-focus-positions), [Transducer-target distance](#transducer-target-distance) |

## Amplitude and Phases

An annular array is driven by two settings: the focal depth setting (`transducer.focal_distance_ep`) steers it, and the target free-water ISPPA (`transducer.target_isppa_wcm2`) sets the drive level. PRESTUS derives the element phases and amplitudes from them:

- With a library `serial`, from the device's calibration entry at that depth (the [transducer library](doc_transducer_library.md)).
- Without one (custom geometry), geometric steering phases to that depth (`compute_phases`, as a library device's geometric baseline) and a nominal uniform drive, which ISPPA targeting then scales. Flat rings (`curv_radius_mm: .inf`) need element phases.

The element values are advanced settings, for a drive found by [calibration](doc_calibration.md) against free-water hydrophone measurements:

- `transducer.annular.elem_p_amp`: amplitude of the acoustic source per element (Pa). Without a target ISPPA (`target_isppa_wcm2: []`) it sets the drive level and is required.
- `transducer.annular.elem_phase_deg`: phase of the acoustic source per element (degrees). Replaces the steering to the focal depth setting.

!!! note
    Per-element calibration of phases, gains and drive against axial hydrophone profiles (the [transducer library](doc_transducer_library.md)) is implemented only for annular arrays. Matrix arrays can instead get a [steering correction](doc_calibration_matrix-steering.md#steering-correction-of-matrix-arrays), which scales the element amplitudes by the measured focal-peak factor at each steering target. For matrix arrays, PRESTUS computes the phase of each element from `transducer.trans_pos` and `transducer.focus_pos`. Calibrated phases (one per element) can replace the geometric phases via `transducer.matrix.elem_phase_override_deg`.

## Transducer and Focus Positions

See [transducer placement](doc_placement.md) for how a position is chosen.

- `transducer.trans_pos`: position of the transducer bowl (XYZ; T1 voxels by default, or RAS/MNI mm via `placement.manual.coordinate_space`). Used for placement and for the phase calculation of matrix arrays.
- `transducer.focus_pos`: position of the stimulation target (XYZ; same coordinate space as `trans_pos`). Defines the focal point and the automatic steering of matrix arrays.

In `water` simulations with annular arrays, the placement is determined from the edge of the PML layer.

??? note "Positions in mm for phantom and water simulations"
    For phantom and water simulations the transducer and focus positions can be given in millimetres rather than voxel indices, with two optional config fields:

    | Field | Description |
    |---|---|
    | `trans_pos_mm` | Transducer position in mm from the centre of grid voxel `[1, 1, 1]` (`[x, y, z]`; 0 mm = voxel 1). Converted to voxel indices at runtime as `round(mm / grid.resolution_mm) + 1`. |
    | `focus_pos_mm` | Focus position in mm from the centre of grid voxel `[1, 1, 1]` (`[x, y, z]`). Same conversion applies. |

    These fields are resolution-independent: the same config works unchanged when `grid.resolution_mm` changes. When both voxel-index and mm fields are present, the mm fields take precedence. Phantom NIfTI outputs reuse the input header, so world-space orientation and voxel dimensions stay consistent across all pipeline stages.

## Transducer-Target Distance

![PRESTUS transducer distance definitions](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/transducer_distances.png)

The figure shows the PRESTUS definitions for annular transducer distances. Specify one of the two fields. Without either, PRESTUS computes a geometric distance from `trans_pos` and `focus_pos`, if both are set.

- `focal_distance_ep`: expected distance from the transducer exit plane to the stimulation focus (mm). This often corresponds to the focal depth setting of the driving system or the transducer calibration ([BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain)'s TPO distance). The exit plane lies `curv_radius_mm − dist_geom_ep_mm` in front of the bowl apex; `dist_geom_ep_mm` is BabelBrain's `NaturalOutPlaneDistance`.
- `focal_distance_bowl`: expected distance from the transducer bowl to the stimulation focus (mm). Computed internally from `focal_distance_ep`, `dist_geom_ep_mm` and `curv_radius_mm`.

## Annular Array

An annular array consists of concentric ring-shaped elements on a spherically curved surface. All elements share one acoustic axis, so focusing is axial only. Lateral steering is not possible; the transducer must be repositioned to change the target.

Use an annular array when your hardware is a single-focus or multi-ring bowl transducer with calibrated per-element amplitudes and phases. This is the simpler and better-validated configuration.

Element geometry is defined by inner and outer diameters (`elem_id_mm`, `elem_od_mm`) and bowl curvature (`curv_radius_mm`). Each element has an independent pressure amplitude (`elem_p_amp`) and phase (`elem_phase_deg`), calibrated empirically. PRESTUS includes a calibration pipeline for annular arrays that matches the simulated pressure to free-water measurements.

### Flat Elements

`curv_radius_mm: .inf` makes the elements flat, e.g. a plane piston (the [ITRUSST](https://itrusst.org/) benchmark's 20 mm source). The face lies at `trans_pos`, normal to the `trans_pos` → `focus_pos` direction, and is also the exit plane (bowl depth 0); `focus_pos` only sets the beam direction. Supported with the kWaveArray source (`grid.use_kWaveArray: 1`):

| Grid | Element | Rings (`elem_id_mm` > 0) |
|---|---|---|
| 2D axisymmetric | line across the axis (`addLineElement`) | yes |
| 2D | line normal to the beam | no |
| 3D | disc (`addDiscElement`) | no (single disc) |

Not supported: the makeBowl source (`use_kWaveArray: 0`), BabelBrain, transducer calibration and the placement planners, which assume a bowl.

## Matrix Array

A matrix array consists of many small elements on a curved or flat surface. Each element is driven with an independent phase delay, so the focus can be steered electronically in three dimensions without moving the transducer.

Use a matrix array when your hardware supports electronic steering (e.g. a phased array system) or to simulate volumetric targeting without repositioning the transducer. The focal pressure a steered array delivers can be calibrated against hydrophone measurements over the steering range: see [Steering correction of matrix arrays](doc_calibration_matrix-steering.md#steering-correction-of-matrix-arrays).

Phase delays are computed from the geometry of each element, `trans_pos` and `focus_pos`; no manual phase input is required. `elem_phase_override_deg` replaces them with calibrated phases.

??? note "Matrix element layout and shape"
    - Element positions are defined analytically (Fibonacci, Fermat spiral or rectangular grid) or loaded from a coordinate file.
    - Element shape is rectangular (`rect`), circular (`disc`) or curved (`bowl`). The aperture of `disc`/`bowl` elements is equal-area to `elem_width_mm` × `elem_height_mm` unless `elem_diameter_mm` is set.
    - For curved arrays, `curv_radius_mm` and `dist_geom_ep_mm` define the bowl geometry used for visualisation and offset calculations.
    - Flat arrays (`is_curved: false`) have every element facing +z and need `rect` or `disc` elements; their coordinate file is in the transducer frame (origin at the array centre).
    - Rectangular elements laid out in more than one orientation take a per-element rotation from the coordinate file (`extract_from_file.rot_col`).
    - An optional Clover configuration replicates the array into a multi-aperture arrangement of `n_leaves` leaves on a parent bowl (`ROC_parent`), `az_sep_deg` apart (default `360 / n_leaves`) and tilted by `elevation_deg` (default: the smallest non-overlapping tilt).

### Beam Axis Constraint

Element orientation and phase delay computation (`create_matrix_karray`) assume that the transducer's natural beam axis is aligned with the grid +Z direction. This holds for the standard single-transducer placement (transducer at one face of the grid, firing toward the target). In multi-transducer configurations where a matrix transducer fires along a different axis (e.g. a second transducer opposing the first along −Z), element normals and phase delays are computed incorrectly. Annular transducers are not affected. A matrix bilateral setup would need the beam direction derived from `(focus_pos − trans_pos)` rather than hardcoded to `[0, 0, 1]` in `create_matrix_karray`.

??? note "KWaveArray vs. simple elements"
    `source_create` builds the acoustic source for [k-Wave](https://github.com/ucl-bug/k-wave) simulations in one of two ways:

    1. **Realistic kWaveArray elements** (`use_kWaveArray = 1`, recommended). Uses the k-Wave `kWaveArray` class to define elements via physical parameters (position, radius, apodisation) rather than discrete grid voxels. Elements render smoothly at any grid resolution. Independent per-element delays and weights enable dynamic focusing. Element overlap is handled through weighted contributions. Supports annular rings, axisymmetric 2D simulations, matrix arrays and GPU acceleration.
    2. **Simple custom elements** (`use_kWaveArray = 0`). Creates bowl-shaped (3D) or arc-shaped (2D) pressure source regions per element with the k-Wave `makeBowl` function. Each element receives an identical continuous wave (CW) signal. Suited to uniform single-element setups. Grid discretisation can introduce staircasing artefacts in the source mask.
