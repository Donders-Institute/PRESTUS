# Neuronavigation read-in

Extracting neuronavigation coordinates for post-hoc ultrasound simulations is critical to verify and predict the focus of the realized stimulation sites during neuromodulation experiments. Neuronavigation systems provide spatial tracking of the ultrasound transducer position relative to the individual subject’s anatomy in real-time. 

PRESTUS provides functions that aim to facilitate the read-in and preprocessing of positions recorded with the Localite system. 
To get started, see an example demo (without provided data) in `/examples/demo_localite.m`.

Three Localite XML file types are supported. All are parsed via the same pipeline (`neuronav_compute_series_statistics`), which handles single markers, repeated trigger trains, and multi-position recordings uniformly.

<details markdown>
<summary>Localite XML file formats — TriggerMarkers, GUMMarkers, InstrumentMarker</summary>

### `TriggerMarkers`

Multiple triggers can be sent during stimulation to acquire updates to the transducer location over time. Files are named `TMSTrigger/TriggerMarkers_Coil0_<timestamp>`. Recorded locations are averaged within each position series (identified by time-gap segmentation). Multiple positions per recording session are supported.

### `GUMMarkers`

These files contain one or more static marker positions. Each `Element` in the file is expected to be of type `EntryTargetPair` (with a `Target/Marker/ColVec3D` position and a `Rotation/RotationReference/ColVec3D` beam direction), from which a synthetic 4×4 matrix is constructed.

> **Mixed-type GUMMarkers files.** Some Localite versions write a `GUMMarkerList` file that mixes `EntryTargetPair` and `InstrumentMarker` elements (`type="InstrumentMarker"` on the `Element` tag). Each element becomes its own series: `InstrumentMarker` elements keep their recorded 4×4 pose, `EntryTargetPair` elements the synthesized matrix. Only the former are instrument poses; placement refuses an `EntryTargetPair` series (`PRESTUS:localite:notAPose`).

### `InstrumentMarker`

Static instrument positions stored as full 4×4 matrices in an `InstrumentMarkerList` file. Each `InstrumentMarker` entry contains a `Marker/Matrix4D` block with the same `data{row}{col}` field layout as TriggerMarkers. The ingested output of `neuronav_ingest_markers` uses this format.

### `EntryTarget` (not for position/angle -- see the multi-target toolkit)

Some Localite versions export an `EntryTargetList` XML containing flat `Target/Marker/ColVec3D` and `Rotation/RotationReference/ColVec3D` entries at the top level (not wrapped in `GUMMarkerList/Element`). This format is **not** handled by `neuronav_compute_series_statistics` and is never used for position/angle computation in the placement pipeline below -- that restriction is deliberate, not a gap; see "Choosing a file" below for why.

It **is** read by a separate function, `neuronav_read_entry_target`, for exactly one purpose: matching a static, operator-confirmed brain-target label against a `target_map` in multi-target studies (e.g. "does this subject's EntryTarget export confirm a 'PUL' target, and if so, where?"). If you need that -- naming multiple targets per session and resolving each one's confirmed coordinate with source-tracked fallbacks -- see [doc_placement_neuronav_multitarget.md](doc_placement_neuronav_multitarget.md), a toolkit layered on top of everything in this document rather than a replacement for it.

</details>

### Choosing a file: realized angle vs. planned trajectory

The file types above answer different questions. Picking the wrong one will silently mismatch your reported angle and your reported target — read this before choosing which file to use.

- **The realized angle at the moment of stimulation** is only available from `TriggerMarkers` (the tracked pose during the actual pulses, averaged per series) or `InstrumentMarker` (a static saved pose). Both store a genuine, measured rotation. Use one of these if you need "the angle at which we stimulated."
- **`GUMMarkers` (`EntryTargetPair`) does not carry a measured angle.** Its 4×4 matrix is *synthesized* from the static `Target` point plus a `RotationReference` vector — and that reference is frequently left unset in the source file (`<Rotation angle="0.0" used="false">`). When `used="false"`, the "angle" `neuronav_compute_series_statistics` produces from it is a placeholder, not captured stimulation data. Only treat it as a real angle if the source XML has `Rotation used="true"`.
- **Static planned Entry/Target coordinates are never read by the placement pipeline below.** This includes the `Target` field inside a `GUMMarkers` `EntryTargetPair`, and a standalone `EntryTarget.xml` read via `neuronav_read_entry_target`. `position_transducer_localite` → `localite_matrix_to_positions` never reads them, and never will — see the axis-projection point below for why. **If you are only running single-target placement, these static coordinates are entirely optional** — TriggerMarkers or InstrumentMarker alone is sufficient; no GUMMarkers or EntryTarget export is required. They matter only for the separate multi-target extraction toolkit ([doc_placement_neuronav_multitarget.md](doc_placement_neuronav_multitarget.md)), which uses them for two things: identifying which planned target a given coil's series corresponds to (`neuronav_match_entry_targets`), and, for studies that supply pre-hoc planned coordinates, as an independent cross-check against those (`neuronav_resolve_target`).
- **The reported focus position is axis-projected, not a stored coordinate.** `localite_matrix_to_positions` computes it as `focus_pos = trans_pos + focal_distance_bowl * ref_vec`, where `trans_pos = ref_pos + tracker_to_bowl_mm * ref_vec` is the bowl rear centre (with the default exit-plane origin this equals `ref_pos + focal_distance_ep * ref_vec`) — projecting forward from the *tracked* position along the *tracked* beam axis by the transducer's calibrated focal distance. This keeps the reported target self-consistent with the reported angle. A separately-clicked static Target point would not necessarily lie on the tracked beam axis at all, and could silently disagree with the reported angle if used instead.

### Localite I/O

Localite data pass through two stages (run them from the GUI or a script — see [GUI and script workflows](doc_placement_neuronav_workflows.md)). **Ingestion** (`neuronav_ingest_markers`, or the `neuronav_gui` front-end) reads the raw exports of one subject/session, averages each recorded series and writes one normalized file pair. **Placement** (`placement.mode: 'localite'`) and the multi-target toolkit read that pair — or, in file mode, a single XML you point to.

#### Input layout

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
| De-duplication | A file byte-identical to one in the previous session (`ses-(XX-1)`) is skipped — Localite re-exports older sessions |
| Series | TriggerMarkers pulses are split into series at time gaps > 3× the median interval; series whose pulse count differs from `expected_segment_length` (default 80) by more than 20 are dropped (`[]` keeps all) |
| Coils | Without `neuronav.coil_map`: Coil0 only. With it: one `TriggerMarkers_Coil<N>` file per entry; each coil contributes one series per configured sequence; untracked blocks (all positions at 0) are skipped |

The planning T1 that Localite coordinates refer to (see below) is looked up as `<path.localite_t1>/[sub-XXX/]<path.localite_t1_pattern>` (default pattern `sub-XXX_*T1*.nii*`; `path.localite_t1` defaults to `path.localite_raw`).

#### Output layout

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

#### Path keys

| Key | Meaning |
|---|---|
| `path.localite_raw` | Raw Localite exports (read-only) |
| `path.localite_post` | Ingested outputs (written by ingestion, read by placement) |
| `path.localite_t1` | Planning T1s; empty = `path.localite_raw` |
| `path.localite_t1_pattern` | File-name template of the planning T1 (`%s` = `sub-XXX`) |
| `path.localite_export` | Outputs meant for import into Localite (heuristic placement: position table, plot, Localite-ready T1), per subject in `sub-XXX/` |
| `path.localite` | **Deprecated.** Read as `path.localite_raw` and as `path.localite_export` when those are empty (with a warning for the former) |

### Pipeline entry point

Set `placement.mode: 'localite'` to place the transducer from a recorded Localite pose (`preproc_transducer_placement`; the older `placement.localite.enabled = 1` switch in `preproc_head` uses the same code). Both modes go through `neuronav_localite_pose`, and the pose's column 4 (position) and column 1 (beam direction) go to `localite_matrix_to_positions`.

#### File mode — single XML file

Set `placement.localite.file`. TriggerMarkers, GUMMarkers and InstrumentMarker files are accepted; the type is detected from the content.

```yaml
placement:
  mode: 'localite'
  localite:
    file: '/path/to/TriggerMarkers_Coil0_20240820122011842.xml'
    position: 2        # detected series in the file (1-based, after the pulse-count filter), or:
    target: ''         # description of an InstrumentMarker / GUMMarkers entry
```

#### Automatic mode — ingested markers

Leave `placement.localite.file` empty. The ingested pair under `path.localite_post` is read; if it does not exist yet and `path.localite_raw` is set, ingestion runs first (with `neuronav.coil_map` if configured, otherwise for the requested series only, named `series_<position>` or `target`).

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

To choose series and names interactively, ingest with `neuronav_gui` first; placement then only reads the result. A raw series that was not ingested is an error (`PRESTUS:localite:noSeries`) rather than a silent fallback, and an untracked series (position at the origin) is refused.

#### Several transducers (bilateral setups)

`target` and `position` take one entry per transducer; transducer *i* uses entry *i*, and each gets its own pose. A single value still applies to every transducer, and a list of the wrong length is an error (`PRESTUS:localite:listLength`).

```yaml
placement:
  mode: 'localite'
  localite:
    session: 2
    target: [L_MD, R_MD]   # transducer(1) -> L_MD, transducer(2) -> R_MD
```


#### Standalone / batch use

For multi-position or multi-session workflows outside the pipeline, use `neuronav_ingest_markers` → `neuronav_select_localite` → `neuronav_compute_series_statistics` → `neuronav_convert_trigger_to_voxels` directly (see `examples/demo_localite.m` and [doc_placement_neuronav_multitarget.md](doc_placement_neuronav_multitarget.md)).

All paths share the same XML parser and position math (`localite_matrix_to_positions`). The reference distance from the IR tracker attachment point to the transducer exit plane is derived automatically from transducer geometry (`curv_radius_mm − dist_geom_ep_mm`) when defined, with fallback to `parameters.placement.localite.reference_distance_mm`.

> **`transducer.focal_distance_ep` (or `focal_distance_bowl`) must be set explicitly for Localite placement — there is no default.** `config_default.yaml` ships it as `[]`, and unlike other placement modes it cannot be derived automatically here: the usual fallback computes it from `trans_pos`/`focus_pos`, but in Localite mode those positions are exactly what this step is trying to produce, so the fallback never has anything to work from. `localite_matrix_to_positions` validates this upfront and raises a clear `PRESTUS:localite:missingConfig` error naming the missing field if neither is set — this same check also covers `placement.mode = 'plantus'`, since PlanTUS placement resolves its final position through the same function.

Average transducer and target locations can be provided in multiple coordinate spaces (native image space being the most relevant for PRESTUS). The functions also facilitate standard-space definitions (experimental), e.g., to plot average transducer positions in standard space.

### Neuronavigation Coordinate Transforms

See the [coordinate system documentation](doc_coordinate_systems.md).

Recorded Localite coordinates are relative to a neuronavigation planning image (e.g., `T1_forneuronav`) that was loaded in the neuronavigation software. Neuronavigation software writes this planning image with a **zero-origin, diagonal affine** (`world = voxel × voxel_size_mm`) rather than a scanner-isocenter world-space affine. Localite trigger matrices are therefore in **planning-image space**, not scanner world-space RAS.

The `InstrumentMarker` XML written by `neuronav_ingest_markers` preserves these coordinates as-is and labels the coordinate space `"RAS"` — meaning the axes are RAS-oriented, but the origin remains at the image corner (zero-origin). To convert to voxels, use `ras_to_grid()` with the **planning T1 header**, not the SimNIBS m2m T1 header. The m2m T1 has a full world-space affine; applying its inverse to planning-image coordinates displaces positions by tens of millimetres.

When SimNIBS is run on the planning T1, the resulting `m2m_{sub}/T1.nii.gz` shares the same voxel grid (SimNIBS preserves voxel dimensions and layout, updating only the affine). Voxel indices computed via `ras_to_grid(localite_mm, planning_t1_header)` are therefore directly valid for PRESTUS simulations that load the m2m T1 — no reprojection is needed. See [Coordinate Reference Systems](doc_coordinate_systems.md#localite-planning-image-zero-origin-affine-and-simnibs-voxel-grid).

SimNIBS' nonlinear transform matrices to MNI are relative to the `final_tissues.nii.gz` segmentation image, so PRESTUS applies a transform from the planning to the segmentation image before porting coordinates to MNI space.

1.	Trigger to Voxel (`neuronav_convert_trigger_to_voxels`): → subject voxel indices
2.	Voxel → mm (`neuronav_grid_to_mm_batch`): → subject mm
3.	Native mm → MNI mm (`neuronav_apply_deformation`): → MNI mm

[Optional]: Interpolate average positions across subjects

4.	Aggregate means in MNI mm across subjects (`neuronav_get_group_mean_mni`)
5.	Group mean MNI mm → subject mm (`neuronav_apply_inverse_deformation`)
6.	Subject mm → subject voxel (`neuronav_mm_to_voxel`)

#### Coordinate Spaces in neuronav functions

- `neuronav_convert_trigger_to_voxels`
	- Localite/trigger data (planning-image space mm, zero-origin RAS) → subject (native) voxel space via planning T1 header
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

<details markdown>
<summary>Localite matrix origin: physical interpretation and calibration</summary>

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

</details>