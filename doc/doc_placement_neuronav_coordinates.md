# Localite Files, Offset and Coordinates

A Localite pose is a 4×4 matrix in the coordinates of the planning image, and its origin is a point on the tracked instrument that need not be the transducer. A wrong origin offset shifts the whole simulation along the beam axis, and converting with the wrong image header displaces positions by tens of millimetres. This page gives the folders PRESTUS reads and writes, how the offset from the tracked origin to the transducer is set, and how positions are converted to voxels and MNI space.

On this page:

- [Localite File Layout](#localite-file-layout): the folder layout PRESTUS reads from and writes to, and the path keys.
- [Transducer Offset (`tracker_to_bowl_mm`)](#transducer-offset-tracker_to_bowl_mm): from the tracked matrix origin to the bowl apex.
- [Coordinate Transforms](#coordinate-transforms): planning-image space, native voxels and MNI.

!!! note "At a glance"
    - **When you need it:** You set up the Localite paths, check where the tracked origin sits on the transducer, or convert Localite positions to voxels or MNI.
    - **Before you start:** The planning T1 that was loaded into Localite and the transducer geometry (`curv_radius_mm`, `dist_geom_ep_mm`).

---

## Localite File Layout

### Input Layout

```
<path.localite_raw>/
└── sub-XXX/                          subject folder (pipeline: sub-%03d of subject_id)
    └── ses-XX/                       session (a number n in the config means ses-0n)
        └── <placement.localite.loc_folder>/     default 'loc'; '' = no such level
            └── …any depth…           e.g. Session_<ts>/ or <patient>/Sessions/Session_<ts>/
                │                     (optionally restricted by placement.localite.session_folder_pattern)
                ├── TMSTrigger/TriggerMarkers_Coil<N>_<yyyyMMddHHmmssSSS>.xml
                ├── GUMMarkers<yyyyMMddHHmmssSSS>.xml       any depth
                └── EntryTarget<yyyyMMddHHmmssSSS>.xml      any depth; planned targets, used for naming only
```

This is Localite's own export structure placed below a `sub-XXX/ses-XX/` pair; nothing inside the export needs renaming. Selection rules:

| Rule | Detail |
|---|---|
| Marker type | TriggerMarkers first; GUMMarkers only if no usable TriggerMarkers file exists |
| Newest wins | Files are ordered by the 17-digit timestamp in the file name (file dates change on copying) |
| Size | TriggerMarkers files ≤ 10 kB (no or few pulses) are skipped |
| De-duplication | A file byte-identical to one in the previous session (`ses-(XX-1)`) is skipped; Localite re-exports older sessions |
| Series | TriggerMarkers pulses are split into series at time gaps > 3× the median interval; series whose pulse count differs from `expected_segment_length` (default 80) by more than 20 are dropped (`[]` keeps all) |
| Coils | Without `neuronav.coil_map`: Coil0 only. With it: one `TriggerMarkers_Coil<N>` file per entry; each coil contributes one series per configured sequence; untracked blocks (all positions at 0) are skipped |

The planning T1 that Localite coordinates refer to (see [Coordinate transforms](#coordinate-transforms)) is looked up as `<path.localite_t1>/[sub-XXX/]<path.localite_t1_pattern>` (default pattern `sub-XXX_*T1*.nii*`; `path.localite_t1` defaults to `path.localite_raw`).

### Output Layout

```
<path.localite_post>/sub-XXX/ses-XX/
├── InstrumentMarker_sub-XXX_ses-XX.xml    InstrumentMarkerList, one entry per target:
│                                          Marker description = target name, Matrix4D = mean pose of the series
└── InstrumentMarker_sub-XXX_ses-XX.json   QC sidecar, one object per target (key = target name)
```

| JSON field | Meaning |
|---|---|
| `series_index` | Index of the detected series (after the pulse-count filter; as used by `placement.localite.position`) |
| `transducer_id` | Coil index |
| `markertype` | `TriggerMarkers` or `GUMMarkers` |
| `source_file` | Raw file the series came from |
| `N_pulses` | Pulses averaged |
| `position_dev_mm.mean/std` | Distance of the pulses from the series mean position [mm] |
| `rotation_dev_rad.mean/std` | Angle of the pulses from the mean orientation [rad] |

With `io.overwrite_files: 'never'`, an existing pair is kept.

### Path Keys

| Key | Meaning |
|---|---|
| `path.localite_raw` | Raw Localite exports (read-only) |
| `path.localite_post` | Ingested outputs (written by ingestion, read by placement) |
| `path.localite_t1` | Planning T1s; empty = `path.localite_raw` |
| `path.localite_t1_pattern` | File-name template of the planning T1 (`%s` = `sub-XXX`) |
| `path.localite_export` | Outputs meant for import into Localite (heuristic placement: position table, plot, Localite-ready T1), per subject in `sub-XXX/` |
| `path.localite` | **Deprecated.** Read as `path.localite_raw` and as `path.localite_export` when those are empty (with a warning for the former) |

## Transducer Offset (`tracker_to_bowl_mm`)

PRESTUS models the transducer position as the bowl rear centre, the [bowl apex](doc_glossary.md) (`trans_pos`; `bowl_pos` in [k-Wave](https://github.com/ucl-bug/k-wave)'s `makeBowl`); the acoustic focus lies `focal_distance_bowl` beyond it along the beam direction. What the Localite matrix origin (column 4) refers to depends on how the instrument was registered in Localite. The offset from that origin to the bowl rear centre is taken, in this order, from `placement.localite.tracker_to_bowl_mm` (explicit), the transducer geometry `-(curv_radius_mm − dist_geom_ep_mm)` (assumes the origin is at the exit plane) or `placement.localite.reference_distance_mm` (fallback). All three share one sign convention along the beam direction (column 1, towards the head): negative puts the bowl rear centre behind the matrix origin, positive in front of it. A wrong assumption shifts the whole simulation along the beam axis by the difference between the assumed and actual origin.

```yaml
placement:
  localite:
    tracker_to_bowl_mm: -15   # example: mm from the Localite matrix origin to the bowl rear centre
                              # negative = bowl rear is on the housing side (away from head)
```

The shipped fallback `reference_distance_mm: 15` is positive: it places the bowl rear centre 15 mm in front of the matrix origin, which fits only an origin behind the bowl rear. PRESTUS does not document which physical origin this value assumes. The geometry-derived offset is used whenever `curv_radius_mm` and `dist_geom_ep_mm` are set (a flat element gives 0); with empty geometry PRESTUS falls back to `reference_distance_mm`. Set `tracker_to_bowl_mm` when the origin is not at the exit plane.

??? note "Determining the offset"
    Column 4 of the Localite 4×4 matrix is the position in RAS mm; column 1 is the beam direction (from the transducer face toward the target). These are the only two columns PRESTUS uses.

    | Possible matrix origin | Value of `tracker_to_bowl_mm` |
    |---|---|
    | Exit plane centre | `-(curv_radius_mm − dist_geom_ep_mm)` (negative; bowl rear is behind exit plane) |
    | Bowl rear centre | `0` |
    | IR tracker cluster centroid | Measured distance from cluster to bowl rear centre (positive when the cluster sits behind the bowl rear) |
    | Housing reference point | Measured distance from reference to bowl rear centre |

    **Option 1: Localite instrument definition file (recommended).** The instrument definition is an XML file, stored in the Localite installation under `Instruments/` or provided by the manufacturer. Its `<ReferencePoint>` or `<Origin>` tag describes which physical point the matrix origin corresponds to and its distance from transducer landmarks in mm.

    **Option 2: Empirical measurement (ruler / phantom).** Place the transducer on a flat phantom surface and record a Localite position. Measure the distance from the phantom surface (= exit plane) to the IR tracker attachment on the housing. The bowl rear centre is about `curv_radius_mm − dist_geom_ep_mm` mm behind the exit plane. Set `tracker_to_bowl_mm` = (bowl rear centre position) − (matrix origin position), in mm along the beam axis. Alternatively, acquire an MR image of a gel phantom with the transducer in place and compare the PRESTUS-predicted bowl voxel position against the visible gel/transducer interface.

    **Option 3: Contact Localite support** (https://www.localite.de) with the instrument type and firmware version used during the recording.

    `tracker_to_bowl_mm` takes priority over the geometry-derived distance and over `reference_distance_mm`. The geometry-derived value is correct only when the matrix origin is at the exit plane.

## Coordinate Transforms

See the [coordinate system documentation](doc_coordinate_systems.md).

Recorded Localite coordinates are relative to the planning image (e.g. `T1_forneuronav`) loaded in the neuronavigation software. That image has a **zero-origin, diagonal affine** (`world = voxel × voxel_size_mm`), not a scanner-isocentre world-space affine. Localite trigger matrices are therefore in **planning-image space**, not scanner world-space RAS.

The `InstrumentMarker` XML written by `neuronav_ingest_markers` keeps these coordinates as they are and labels the space `"RAS"`: the axes are RAS-oriented, but the origin remains at the image corner. To convert to voxels, use `ras_to_grid()` with the **planning T1 header**, not the [SimNIBS](https://simnibs.github.io/simnibs/) m2m T1 header. The m2m T1 has a full world-space affine; applying its inverse to planning-image coordinates displaces positions by tens of millimetres.

When SimNIBS is run on the planning T1, the resulting `m2m_{sub}/T1.nii.gz` shares the same voxel grid (SimNIBS keeps voxel dimensions and layout and updates only the affine). Voxel indices from `ras_to_grid(localite_mm, planning_t1_header)` are therefore valid for PRESTUS simulations that load the m2m T1; no reprojection is needed. See [Coordinate Reference Systems](doc_coordinate_systems.md#localite-planning-image-header).

SimNIBS' nonlinear transform matrices to MNI are relative to the `final_tissues.nii.gz` segmentation image, so PRESTUS applies a transform from the planning to the segmentation image before porting coordinates to MNI space. The chain is:

1. Trigger to voxel (`neuronav_convert_trigger_to_voxels`): subject voxel indices
2. Voxel to mm (`neuronav_grid_to_mm_batch`): subject mm
3. Native mm to MNI mm (`neuronav_apply_deformation`): MNI mm

Optional, to average positions across subjects:

4. Aggregate means in MNI mm across subjects (`neuronav_get_group_mean_mni`)
5. Group mean MNI mm to subject mm (`neuronav_apply_inverse_deformation`)
6. Subject mm to subject voxel (`neuronav_mm_to_voxel`)

??? note "Coordinate spaces of the neuronav functions"
    | Function | Input | Operation | Output |
    |---|---|---|---|
    | `neuronav_convert_trigger_to_voxels` | Localite/trigger data (planning-image space mm, zero-origin RAS) | Planning T1 header | Subject (native) voxel indices |
    | `neuronav_grid_to_mm_batch` | Subject voxel indices (i, j, k) | Subject T1 [NIfTI](https://nifti.nimh.nih.gov/) affine | Subject native mm (scanner/world space) |
    | `neuronav_apply_deformation` (forward warp) | Native mm locations | Queries the subject's `Conform2MNI_nonl.nii.gz` | MNI mm (for reporting, plotting, averaging) |
    | `neuronav_get_group_mean_mni` | CSVs with columns like `Mtrans_pos_MNI_x/y/z` | Mean across subjects/sessions | MNI mm group-mean locations (mm, not voxels) |
    | `neuronav_apply_inverse_deformation` (inverse warp) | MNI mm group-mean locations | Subject's `MNI2Conform_nonl.nii.gz` nonlinear field | Native mm location in that subject's world/scanner space |
    | `neuronav_mm_to_voxel` (helper) | Native mm (x, y, z) | Subject T1 affine | Subject voxel indices |
    | `neuronav_export_session_csv` | | | Native voxel indices for target and transducer, native mm, MNI mm (from the deformation step) and MNI voxel indices (may be inaccurate with limited FOVs) |
