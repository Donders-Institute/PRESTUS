# Neuronavigation Read-In

!!! note "At a glance"
    - **What:** Simulating the transducer pose that was recorded with Localite neuronavigation during the experiment.
    - **When you need it:** After a session, to simulate where the transducer actually was.
    - **Before you start:** Localite exports of the session and the planning T1 that was loaded into Localite.

PRESTUS reads transducer poses recorded with the Localite neuronavigation system and simulates them post hoc, to verify the focus of the realised stimulation. An example demo (without data) is in `examples/demo_localite.m`.

Localite data pass through two stages. **Ingestion** (`neuronav_gui`, or `neuronav_ingest_markers` in a script; see [GUI and script workflows](doc_placement_neuronav_workflows.md)) reads the raw exports of one subject and session, averages each recorded series and writes one normalised file pair. **Placement** (`placement.mode: 'localite'`) reads that pair, or in file mode a single XML you point to.

## Placing the transducer from a recorded pose

Set `placement.mode: 'localite'`. The pose's column 4 (position) and column 1 (beam direction) become the transducer and focus positions.

!!! warning "Focal distance required"
    `transducer.focal_distance_ep` (or `focal_distance_bowl`) must be set for Localite placement; `config_default.yaml` ships it as `[]`. Other placement modes derive it from `trans_pos`/`focus_pos`, but here those positions are what the step produces. If neither is set, `localite_matrix_to_positions` raises `PRESTUS:localite:missingConfig` naming the missing field. The same check covers `placement.mode = 'plantus'`, which resolves its final position through the same function.

### Automatic mode: ingested markers

Leave `placement.localite.file` empty. The ingested pair under `path.localite_post` is read. If it does not exist yet and `path.localite_raw` is set, ingestion runs first (with `neuronav.coil_map` if configured, otherwise for the requested series only, named `series_<position>` or `target`).

```yaml
path:
  localite_raw:  '/data/localite_raw'
  localite_post: '/data/derivatives/localite'
placement:
  mode: 'localite'
  localite:
    session: 2          # -> ses-02
    target: 'L_MD'      # ingested target name, or leave '' and use:
    position: 2         # detected-series index, matched against series_index in the JSON sidecar
```

To choose series and names interactively, ingest with `neuronav_gui` first; placement then only reads the result. A raw series that was not ingested is an error (`PRESTUS:localite:noSeries`), and an untracked series (position at the origin) is refused.

### File mode: a single XML file

Set `placement.localite.file`. TriggerMarkers, GUMMarkers and InstrumentMarker files are accepted; the type is detected from the content.

```yaml
placement:
  mode: 'localite'
  localite:
    file: '/path/to/TriggerMarkers_Coil0_20240820122011842.xml'
    position: 2        # detected series in the file (1-based, after the pulse-count filter), or:
    target: ''         # description of an InstrumentMarker / GUMMarkers entry
```

### Several transducers (bilateral setups)

`target` and `position` take one entry per transducer; transducer *i* uses entry *i*, and each gets its own pose. A single value applies to every transducer. A list of the wrong length is an error (`PRESTUS:localite:listLength`).

```yaml
placement:
  mode: 'localite'
  localite:
    session: 2
    target: [L_MD, R_MD]   # transducer(1) -> L_MD, transducer(2) -> R_MD
```

### Many subjects and sessions

For a study, the Nexus's Study tab turns the `neuronav_gui` outputs into manifest rows (**Import from Localite…**). It reads `neuronav_ingest_summary.csv` and the per-session `InstrumentMarker_*.json` files under `path.localite_post`, and makes either one run per target, or one run per session with the targets assigned to transducers 1, 2, … (left first, alphabetical, or a given order). The runs use `placement.mode: localite`, so PRESTUS reads the ingested markers itself. See [Group analysis](doc_group.md#from-localite-ingestion).

??? note "Details: entry points and standalone use"
    Placement runs in `preproc_transducer_placement`; the older `placement.localite.enabled = 1` switch in `preproc_head` uses the same code. Both go through `neuronav_localite_pose`, and the pose's column 4 (position) and column 1 (beam direction) go to `localite_matrix_to_positions`. All paths share the same XML parser and position math.

    For multi-position or multi-session workflows outside the pipeline, call `neuronav_ingest_markers` → `neuronav_select_localite` → `neuronav_compute_series_statistics` → `neuronav_convert_trigger_to_voxels` directly (see `examples/demo_localite.m` and the [multi-target toolkit](doc_placement_neuronav_multitarget.md)). Average transducer and target locations can be provided in several coordinate spaces; native image space is the one PRESTUS uses. Standard-space definitions (experimental) allow e.g. plotting average transducer positions in standard space.

## Localite I/O

### Input layout

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

### Output layout

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

### Path keys

| Key | Meaning |
|---|---|
| `path.localite_raw` | Raw Localite exports (read-only) |
| `path.localite_post` | Ingested outputs (written by ingestion, read by placement) |
| `path.localite_t1` | Planning T1s; empty = `path.localite_raw` |
| `path.localite_t1_pattern` | File-name template of the planning T1 (`%s` = `sub-XXX`) |
| `path.localite_export` | Outputs meant for import into Localite (heuristic placement: position table, plot, Localite-ready T1), per subject in `sub-XXX/` |
| `path.localite` | **Deprecated.** Read as `path.localite_raw` and as `path.localite_export` when those are empty (with a warning for the former) |

## Which file to use: realised angle or planned trajectory

The Localite file types answer different questions. The wrong one silently mismatches the reported angle and the reported target.

- **The realised angle at the moment of stimulation** is only available from `TriggerMarkers` (the tracked pose during the pulses, averaged per series) or `InstrumentMarker` (a static saved pose). Both store a measured rotation.
- **`GUMMarkers` (`EntryTargetPair`) carries no measured angle.** Its 4×4 matrix is synthesised from the static `Target` point and a `RotationReference` vector, which is often unset in the source file (`<Rotation angle="0.0" used="false">`). With `used="false"`, the angle `neuronav_compute_series_statistics` produces is a placeholder. Treat it as a real angle only if the source XML has `Rotation used="true"`.
- **Static planned Entry/Target coordinates are never read by placement.** This includes the `Target` field of a `GUMMarkers` `EntryTargetPair` and a standalone `EntryTarget.xml`. For single-target placement they are optional: TriggerMarkers or InstrumentMarker alone is sufficient. They matter only for the [multi-target toolkit](doc_placement_neuronav_multitarget.md), which uses them to identify which planned target a coil's series corresponds to (`neuronav_match_entry_targets`) and as a cross-check against pre-hoc planned coordinates (`neuronav_resolve_target`).
- **The reported focus position is axis-projected, not a stored coordinate.** `localite_matrix_to_positions` computes `focus_pos = trans_pos + focal_distance_bowl * ref_vec`, where `trans_pos = ref_pos + tracker_to_bowl_mm * ref_vec` is the bowl rear centre (with the default exit-plane origin this equals `ref_pos + focal_distance_ep * ref_vec`): it projects forward from the tracked position along the tracked beam axis by the calibrated focal distance. This keeps the reported target self-consistent with the reported angle. A separately clicked static Target point need not lie on the tracked beam axis.

??? note "Details: Localite XML file formats"
    All formats are parsed by `neuronav_compute_series_statistics`, which handles single markers, repeated trigger trains and multi-position recordings alike.

    **TriggerMarkers.** Multiple triggers can be sent during stimulation to record the transducer location over time. Files are named `TMSTrigger/TriggerMarkers_Coil0_<timestamp>`. Recorded locations are averaged within each position series (identified by time-gap segmentation). Multiple positions per recording session are supported.

    **GUMMarkers.** One or more static marker positions. Each `Element` is expected to be of type `EntryTargetPair` (with a `Target/Marker/ColVec3D` position and a `Rotation/RotationReference/ColVec3D` beam direction), from which a synthetic 4×4 matrix is constructed. Some Localite versions write a `GUMMarkerList` file that mixes `EntryTargetPair` and `InstrumentMarker` elements (`type="InstrumentMarker"` on the `Element` tag). Each element becomes its own series: `InstrumentMarker` elements keep their recorded 4×4 pose, `EntryTargetPair` elements the synthesised matrix. Only the former are instrument poses; placement refuses an `EntryTargetPair` series (`PRESTUS:localite:notAPose`).

    **InstrumentMarker.** Static instrument positions stored as full 4×4 matrices in an `InstrumentMarkerList` file. Each `InstrumentMarker` entry contains a `Marker/Matrix4D` block with the same `data{row}{col}` field layout as TriggerMarkers. The ingested output of `neuronav_ingest_markers` uses this format.

    **EntryTarget** (not for position or angle). Some Localite versions export an `EntryTargetList` XML with flat `Target/Marker/ColVec3D` and `Rotation/RotationReference/ColVec3D` entries at the top level (not wrapped in `GUMMarkerList/Element`). `neuronav_compute_series_statistics` does not read it, and placement never uses it. `neuronav_read_entry_target` reads it for one purpose: matching an operator-confirmed brain-target label against a `target_map` in multi-target studies (see the [multi-target toolkit](doc_placement_neuronav_multitarget.md)).

## Transducer offset (`tracker_to_bowl_mm`)

PRESTUS models the transducer position as the bowl rear centre (`bowl_pos` in k-Wave's `makeBowl`); the acoustic focus lies `focal_distance_bowl` beyond it along the beam direction. What the Localite matrix origin (column 4) refers to depends on how the instrument was registered in Localite. The offset from that origin to the bowl rear centre is taken, in this order, from `placement.localite.tracker_to_bowl_mm` (explicit), the transducer geometry `-(curv_radius_mm − dist_geom_ep_mm)` (assumes the origin is at the exit plane) or `placement.localite.reference_distance_mm` (fallback). A wrong assumption shifts the whole simulation along the beam axis by the difference between the assumed and actual origin.

```yaml
placement:
  localite:
    tracker_to_bowl_mm: -15   # mm from the Localite matrix origin to the bowl rear centre
                              # negative = bowl rear is on the housing side (away from head)
```

??? note "Details: determining the offset"
    Column 4 of the Localite 4×4 matrix is the position in RAS mm; column 1 is the beam direction (from the transducer face toward the target). These are the only two columns PRESTUS uses.

    | Possible matrix origin | Value of `tracker_to_bowl_mm` |
    |---|---|
    | Exit plane centre | `-(curv_radius_mm − dist_geom_ep_mm)` (negative; bowl rear is behind exit plane) |
    | Bowl rear centre | `0` |
    | IR tracker cluster centroid | Measured distance from cluster to bowl rear centre (likely negative) |
    | Housing reference point | Measured distance from reference to bowl rear centre |

    **Option 1: Localite instrument definition file (recommended).** The instrument definition is an XML file, stored in the Localite installation under `Instruments/` or provided by the manufacturer. Its `<ReferencePoint>` or `<Origin>` tag describes which physical point the matrix origin corresponds to and its distance from transducer landmarks in mm.

    **Option 2: Empirical measurement (ruler / phantom).** Place the transducer on a flat phantom surface and record a Localite position. Measure the distance from the phantom surface (= exit plane) to the IR tracker attachment on the housing. The bowl rear centre is about `curv_radius_mm − dist_geom_ep_mm` mm behind the exit plane. Set `tracker_to_bowl_mm` = (bowl rear centre position) − (matrix origin position), in mm along the beam axis. Alternatively, acquire an MR image of a gel phantom with the transducer in place and compare the PRESTUS-predicted bowl voxel position against the visible gel/transducer interface.

    **Option 3: Contact Localite support** (https://www.localite.de) with the instrument type and firmware version used during the recording.

    `tracker_to_bowl_mm` takes priority over the geometry-derived distance and over `reference_distance_mm`. The geometry-derived fallback applies when the transducer geometry is fully specified and the matrix origin is known to be at the exit plane.

## Coordinate transforms

See the [coordinate system documentation](doc_coordinate_systems.md).

Recorded Localite coordinates are relative to the planning image (e.g. `T1_forneuronav`) loaded in the neuronavigation software. That image has a **zero-origin, diagonal affine** (`world = voxel × voxel_size_mm`), not a scanner-isocentre world-space affine. Localite trigger matrices are therefore in **planning-image space**, not scanner world-space RAS.

The `InstrumentMarker` XML written by `neuronav_ingest_markers` keeps these coordinates as they are and labels the space `"RAS"`: the axes are RAS-oriented, but the origin remains at the image corner. To convert to voxels, use `ras_to_grid()` with the **planning T1 header**, not the SimNIBS m2m T1 header. The m2m T1 has a full world-space affine; applying its inverse to planning-image coordinates displaces positions by tens of millimetres.

When SimNIBS is run on the planning T1, the resulting `m2m_{sub}/T1.nii.gz` shares the same voxel grid (SimNIBS keeps voxel dimensions and layout and updates only the affine). Voxel indices from `ras_to_grid(localite_mm, planning_t1_header)` are therefore valid for PRESTUS simulations that load the m2m T1; no reprojection is needed. See [Coordinate Reference Systems](doc_coordinate_systems.md#localite-planning-image-zero-origin-affine-and-simnibs-voxel-grid).

SimNIBS' nonlinear transform matrices to MNI are relative to the `final_tissues.nii.gz` segmentation image, so PRESTUS applies a transform from the planning to the segmentation image before porting coordinates to MNI space. The chain is:

1. Trigger to voxel (`neuronav_convert_trigger_to_voxels`): subject voxel indices
2. Voxel to mm (`neuronav_grid_to_mm_batch`): subject mm
3. Native mm to MNI mm (`neuronav_apply_deformation`): MNI mm

Optional, to average positions across subjects:

4. Aggregate means in MNI mm across subjects (`neuronav_get_group_mean_mni`)
5. Group mean MNI mm to subject mm (`neuronav_apply_inverse_deformation`)
6. Subject mm to subject voxel (`neuronav_mm_to_voxel`)

??? note "Details: coordinate spaces of the neuronav functions"
    | Function | Input | Operation | Output |
    |---|---|---|---|
    | `neuronav_convert_trigger_to_voxels` | Localite/trigger data (planning-image space mm, zero-origin RAS) | Planning T1 header | Subject (native) voxel indices |
    | `neuronav_grid_to_mm_batch` | Subject voxel indices (i, j, k) | Subject T1 NIfTI affine | Subject native mm (scanner/world space) |
    | `neuronav_apply_deformation` (forward warp) | Native mm locations | Queries the subject's `Conform2MNI_nonl.nii.gz` | MNI mm (for reporting, plotting, averaging) |
    | `neuronav_get_group_mean_mni` | CSVs with columns like `Mtrans_pos_MNI_x/y/z` | Mean across subjects/sessions | MNI mm group-mean locations (mm, not voxels) |
    | `neuronav_apply_inverse_deformation` (inverse warp) | MNI mm group-mean locations | Subject's `MNI2Conform_nonl.nii.gz` nonlinear field | Native mm location in that subject's world/scanner space |
    | `neuronav_mm_to_voxel` (helper) | Native mm (x, y, z) | Subject T1 affine | Subject voxel indices |
    | `neuronav_export_session_csv` | | | Native voxel indices for target and transducer, native mm, MNI mm (from the deformation step) and MNI voxel indices (may be inaccurate with limited FOVs) |
