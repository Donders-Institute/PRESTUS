# Neuronavigation read-in

Extracting neuronavigation coordinates for post-hoc ultrasound simulations is critical to verify and predict the focus of the realized stimulation sites during neuromodulation experiments. Neuronavigation systems provide spatial tracking of the ultrasound transducer position relative to the individual subject’s anatomy in real-time. 

PRESTUS provides functions that aim to facilitate the read-in and preprocessing of positions recorded with the Localite system. 
To get started, see an example demo (without provided data) in `/examples/demo_localite.m`.

Multiple Localite file types are supported. Both are parsed via the same pipeline (`neuronav_compute_series_statistics`), which handles single markers, repeated trigger trains, and multi-position recordings uniformly.

### `TriggerMarkers`

Multiple triggers can be sent during stimulation to acquire updates to the transducer location over time. Files are named `TMSTrigger/TriggerMarkers_Coil0_<timestamp>`. Recorded locations are averaged within each position series (identified by time-gap segmentation). Multiple positions per recording session are supported.

### `GUMMarkers`

These files contain one or more static marker positions. They are parsed by the same `neuronav_compute_series_statistics` pipeline as TriggerMarkers.

### Pipeline entry point

Set `parameters.placement.localite.enabled = 1` to use Localite positioning in the PRESTUS pipeline (`prestus_pipeline` → `preproc_head`). Two modes are supported:

#### Simple mode — single XML file

Set `placement.localite.file` to the path of a specific Localite XML file. The file is parsed directly via `position_transducer_localite`.

```yaml
placement:
  localite:
    enabled: 1
    file: '/path/to/GUMMarkers_sub-001.xml'
```

#### Full neuronav mode — automatic file selection

Leave `placement.localite.file` empty and set `path.localite` to the Localite data folder. The pipeline automatically selects the most recent valid file for the subject and session using `neuronav_select_localite`, averages markers across triggers using `neuronav_compute_series_statistics`, and resolves positions via `localite_matrix_to_positions`.

```yaml
path:
  localite: '/path/to/localite/data/'

placement:
  localite:
    enabled: 1
    session: 1              # Session number or 'ses-01'
    markertype: 'TriggerMarkers'  # or 'GUMMarkers'
    position: 1             # Series index when multiple positions are recorded
```

This mode is appropriate when TriggerMarkers files contain repeated acquisitions per position that need to be averaged, or when automatic session-based file deduplication is required.

#### Standalone / batch use

For multi-position or multi-session workflows outside the pipeline, use `neuronav_select_localite` → `neuronav_compute_series_statistics` → `neuronav_convert_trigger_to_voxels` directly (see `examples/demo_localite.m`).

All paths share the same XML parser and position math (`localite_matrix_to_positions`). The reference distance from the IR tracker attachment point to the transducer exit plane is derived automatically from transducer geometry (`curv_radius_mm − dist_geom_ep_mm`) when defined, with fallback to `parameters.placement.localite.reference_distance_mm`.

Average transducer and target locations can be provided in multiple coordinate spaces (native image space being the most relevant for PRESTUS). The functions also facilitate standard-space definitions (experimental), e.g., to plot average transducer positions in standard space.

### Neuronavigation Coordinate Transforms

See the [coordinate system documentation](doc_coordinate_systems.md).

Recorded localite coordinates are relative to a neuronavigation planning image [1] (e.g., `T1_forneuronav`) that was loaded in the neuronavigation software. 

SimNIBS' nonlinear transform matrices to MNI are relative to the `final_tissues.nii.gz` segmentation image however, so PRESTUS initially applies a transform from the planning to the segmentation image to the coordinates before they can be ported to MNI space.

1.	Trigger to Voxel (`neuronav_convert_trigger_to_voxels`): → subject voxel indices
2.	Voxel → mm (`neuronav_grid_to_mm_batch`): → subject mm
3.	Native mm → MNI mm (`neuronav_apply_deformation`): → MNI mm

[Optional]: Interpolate average positions across subjects

4.	Aggregate means in MNI mm across subjects (`neuronav_get_group_mean_mni`)
5.	Group mean MNI mm → subject mm (`neuronav_apply_inverse_deformation`)
6.	Subject mm → subject voxel (`neuronav_mm_to_voxel`)

#### Coordinate Spaces in neuronav functions

- `neuronav_convert_trigger_to_voxels`
	- Localite/trigger data (RAS mm) → subject (native) voxel space
- `neuronav_grid_to_mm_batch`
	- Input: Subject voxel indices (i,j,k)
	- Operation: Uses subject’s T1 NIfTI affine transform
	- Output: Subject native mm coordinates (in scanner/world space)
- `neuronav_apply_deformation` (forward warp)
	- Input: Native mm locations
	- Operation: Queries subject’s Conform2MNI_nonl.nii.gz (subject-native → MNI (mm))
	- Output: MNI mm coordinates (for reporting/plotting/averaging)
- `neuronav_get_group_mean_mni`
	- Input: CSVs containing columns like Mtrans_pos_MNI_x/y/z (all from previous function outputs)
	- Operation: Aggregates these MNI mm locations across subjects/sessions, computes mean (in mm)
	- Output: MNI mm group-mean locations (still in millimeters, not voxels)
- `neuronav_apply_inverse_deformation` (inverse warp)
	- Input: MNI mm group-mean locations (x, y, z in mm)
	- Operation: Uses subject’s MNI2Conform_nonl.nii.gz nonlinear field to map standard points back into that subject’s native space.
	- Output: Native mm location (in subject’s world/scanner space, for that subject)
- `neuronav_mm_to_voxel` (helper)
	- Input: Native mm (x, y, z) location
	- Operation: Uses subject T1 affine to map mm → voxel (i, j, k)
	- Output: Subject voxel indices for subsequent image work or storage
- `neuronav_export_session_csv`
	- Columns include:
	- Native (subject) voxel indices for target/transducer: (i, j, k from subject’s grid)
	- Native mm coordinates (if desired): (x, y, z in subject’s scanner space)
	- MNI mm coordinates for each (from deformation step): (x, y, z in standard space)
	- MNI voxel indices (may be inaccurate with limited FOVs)

---

## Localite matrix origin: physical interpretation and calibration

### Background

The Localite 4×4 transformation matrix encodes the pose of an IR-tracked instrument. Column 4 (the translation vector) gives the 3D position in RAS mm; column 1 (the first axis) gives the beam direction (from the coil/transducer face toward the target). These are the only two columns used by PRESTUS.

**What PRESTUS models as the transducer position** is the *bowl rear centre* — the geometric axis point at the back of the bowl housing. This is `bowl_pos` in k-Wave's `makeBowl`, from which the acoustic focus is offset by `focal_distance_bowl` along the beam direction. The acoustic focus position is therefore self-consistent with this convention.

**The physical meaning of the Localite matrix origin (column 4) is not fixed** — it depends on how the instrument was registered in the Localite instrument definition at the time of recording. Possible interpretations and their implications:

| Possible matrix origin | Value of `tracker_to_bowl_mm` |
|---|---|
| Exit plane centre | `-(curv_radius_mm − dist_geom_ep_mm)` (negative; bowl rear is behind exit plane) |
| Bowl rear centre | `0` |
| IR tracker cluster centroid | Measured distance from cluster to bowl rear centre (likely negative) |
| Housing reference point | Measured distance from reference to bowl rear centre |

The default geometry-derived formula (`-(curv_radius_mm − dist_geom_ep_mm)`) **assumes the matrix origin is at the exit plane**, placing the bowl rear centre behind it by the bowl depth. This assumption holds only if the Localite instrument was registered with the exit plane as the origin — which is a common convention, but is not guaranteed.

An incorrect assumption introduces a systematic spatial error along the beam axis equal to the difference between the assumed and actual origin.

### How to determine the correct offset

**Option 1 — Localite instrument definition file (recommended)**

The Localite instrument definition is an XML file typically stored in the Localite software installation under `Instruments/` or provided by the instrument manufacturer. Open it and locate the `<ReferencePoint>` or `<Origin>` tag. This describes what physical point the matrix origin corresponds to and its distance from key transducer landmarks in mm.

**Option 2 — Empirical measurement (ruler / phantom)**

1. Place the transducer on a flat phantom surface and record a Localite position.
2. Measure the distance from the phantom surface (= exit plane) to the IR tracker attachment on the housing.
3. The bowl rear centre is approximately `curv_radius_mm − dist_geom_ep_mm` mm behind the exit plane.
4. Set `tracker_to_bowl_mm` = (bowl rear centre position) − (matrix origin position), in mm along the beam axis.

Alternatively, acquire an MR image of a gel phantom with the transducer in place and compare the PRESTUS-predicted bowl voxel position against the visible gel/transducer interface.

**Option 3 — Contact Localite support**

The instrument registration specification can be confirmed by Localite (https://www.localite.de). Provide the instrument type and firmware version used during the recording.

### Setting the calibrated offset in PRESTUS

Once the offset is known, set it explicitly in your config:

```yaml
placement:
  localite:
    tracker_to_bowl_mm: -15   # mm from the Localite matrix origin to the bowl rear centre
                              # negative = bowl rear is on the housing side (away from head)
```

`tracker_to_bowl_mm` takes priority over the geometry-derived distance and over `reference_distance_mm`. The geometry-derived fallback (`-(curv_radius_mm − dist_geom_ep_mm)`) remains available when transducer geometry is fully specified and the matrix origin is known to be at the exit plane.