# Localite ingestion: GUI and script workflows

!!! note "At a glance"
    - **What:** Turning raw Localite exports into ingested, named targets that `placement.mode: localite` reads, with the GUI or a script.
    - **When you need it:** Once per session (or in batch for a whole study) before simulating recorded poses.
    - **Before you start:** The raw Localite exports (`path.localite_raw`).

Ingestion turns recorded Localite positions into named, averaged targets per subject and session (folder layout and output files: [Localite I/O](doc_placement_neuronav.md#localite-io)). Run it in `neuronav_gui`, on the PRESTUS Nexus's *Placement ▸ Localite* page (settings from the study config; see [Localite page](doc_nexus.md#localite-page)), or scripted with `neuronav_ingest_batch`. All three run the same functions on the same settings:

| | GUI (`neuronav_gui`) | Script (`neuronav_ingest_batch`) |
|---|---|---|
| Best for | Checking a session, choosing series by eye, setting up the label map | Many subjects and sessions, reproducible re-runs, HPC |
| Session selection | Single session, or tickboxes in the **Batch** window (All / None / Not yet ingested) | A table or `{sub, ses}` list, e.g. rows of `neuronav_scan_sessions` |
| Settings | Form fields, **⬆ Load settings…** / **⬇ Save settings…** (YAML) | The same YAML via `load_parameters` |
| Output | `InstrumentMarker_<sub>_<ses>.xml/.json` + `neuronav_ingest_summary.csv` | identical |

## One settings file for both

Every setting that affects ingestion and naming is a config key. Keep them in the study YAML (or a separate `neuronav_settings.yaml` saved from the GUI) and use that file in both workflows:

```yaml
path:
  localite_raw:  '/data/localite_raw'          # <raw>/sub-XXX/ses-XX/[loc/]…
  localite_post: '/data/derivatives/localite'  # ingested outputs + summary CSV
  localite_t1:   '/data/localite_raw'          # planning T1s (for voxel coordinates); default = localite_raw

placement:
  localite:
    loc_folder: 'loc'              # '' if the export sits directly under ses-XX/
    expected_segment_length: 80    # pulses per TriggerMarkers series ([] = keep all)

io:
  overwrite_files: 'always'        # 'never' keeps existing ingested files

neuronav:
  coil_map: []                     # multi-coil setups, see doc_placement_neuronav_multitarget.md
  suggest_names: true              # name series after the planned target on their beam axis
  axis_tolerance_mm: 10
  target_depth_mm: [0, 150]
  target_labels:                   # planned-target labels as recorded -> study target names
    - target: 'L_PUL'
      labels: ['LEFT PUL', 'L_PUL*']
    - target: 'R_PUL'
      labels: ['RIGHT PUL', 'R_PUL*']
  batch_write: 'all'               # or 'named'
```

A file saved with **⬇ Save settings…** contains exactly these keys and is a valid PRESTUS config: `load_parameters('neuronav_settings.yaml')` merges it with the defaults. Settings the GUI has no field for (for example `transducer` geometry, needed for focus coordinates) come from a config loaded with **⬆ Load settings…**; load the study config first, then adjust.

## How series are named

Each detected series gets a name by one rule (`neuronav_name_series`), the same in the GUI and the batch:

1. **Coil table / `neuronav.coil_map`**: a target name given there for that coil/block.
2. **Planned target on the beam axis** (if `suggest_names` is true). The session's named planned targets (the latest EntryTarget export and the `EntryTargetPair`s of its GUMMarkers exports, each label at its most recent position) are collected by `neuronav_read_label_candidates`. A recorded series takes the target closest to its beam axis if it lies within `axis_tolerance_mm` and at a depth inside `target_depth_mm`. A GUMMarkers `EntryTargetPair` series uses its own label. The label is then passed through `target_labels`.
3. **Default**: `series_<k>`.

Repeated names get `_rep2`, `_rep3`, …; a planned pair that shares its label with a recorded pose in the same file becomes `<label>_target`.

### The label map (`neuronav.target_labels`)

Operators label planned targets inconsistently across sessions (`LEFT PUL` in one session, `L_PUL_1` in the next). `target_labels` maps recorded labels to the names your study uses. Each entry lists label patterns (`*` and `?` wildcards, case-insensitive); the first matching entry wins; a label that matches no entry is kept as recorded. `neuronav_map_label('L_PUL_1', map)` shows the result for a single label. **When `target_labels` is empty, labels are used as found.**

To see which labels the data contain, scan all sessions:

```matlab
parameters = load_parameters('study.yaml');
[sessions, labels] = neuronav_scan_sessions(parameters);
labels(:, {'label', 'stem', 'target', 'n_subjects', 'n_sessions'})
```

`labels.target` is the name each label gets with the current map; `labels.stem` groups numbered variants (the label without a trailing `_<number>`, e.g. `L_caudateV2_1`, `L_caudateV2_2` → `L_caudateV2`). Turn either into a map and store it:

```matlab
map = neuronav_label_map_from_scan(labels, 'stem');   % or edit labels.target, then 'target'
parameters.neuronav.target_labels = map;
yaml.dumpFile('target_labels.yaml', struct('neuronav', struct('target_labels', map)));
```

In the GUI, the **Batch** window shows the same label table: **Group labels by stem** fills the target column of unmapped labels from their stems (existing mappings are kept), **Use labels as found** resets all of them, and any cell can be edited. Only labels whose target differs from the label end up in the saved map. The map applies to single-session naming too.

Check the name sources afterwards: the GUI's **Name from** column and the summary's `rule`, `label`, `axis_dist_mm` and `depth_mm` columns show where each name came from. A match far from the usual focal depth usually means the intended target is missing from the export; narrow `target_depth_mm` to reject such matches.

## Script workflow

```matlab
parameters = load_parameters('study.yaml');   % paths, ingestion, naming, transducer

% 1. Inventory: which sessions exist, which labels they contain
[sessions, labels] = neuronav_scan_sessions(parameters);

% 2. Optional: set or refine the label map (see above)
% parameters.neuronav.target_labels = neuronav_label_map_from_scan(labels, 'stem');

% 3. Ingest all sessions with marker files, or a subset
todo = sessions(sessions.n_trigger + sessions.n_gum > 0, :);
% todo = todo(ismember(todo.sub, {'sub-001', 'sub-004'}), :);
% todo = todo(~todo.ingested, :);          % only sessions not ingested yet
summary = neuronav_ingest_batch(parameters, todo);
```

`neuronav_ingest_batch(parameters)` without a selection processes every session found. For a single session, `neuronav_ingest_batch(parameters, {'sub-001', 'ses-02'})`. Per session, the batch detects the series (nothing is written yet), names them, writes every tracked series (`batch_write: 'all'`) or only the named ones (`'named'`), and adds one summary row per written target. Untracked series and coil dropouts are never written. Sessions that fail or have nothing to write get one row whose `status` says why; the batch continues with the next session.

### Summary CSV

`<path.localite_post>/neuronav_ingest_summary.csv` collects all sessions. Re-running a subset replaces that subset's rows and keeps the others. Columns:

| Column | Meaning |
|---|---|
| `sub`, `ses`, `target` | Session and target name (empty when nothing was written) |
| `series_index`, `coil` | Detected series (as in `placement.localite.position`) and coil index |
| `rule` | `coil_map`, `heuristic` or `default` (see naming) |
| `label`, `label_source` | Planned-target label as recorded, and the file it came from |
| `axis_dist_mm`, `depth_mm` | Distance of that target from the beam axis, and its depth along it |
| `n_pulses`, `pos_dev_mm`, `rot_dev_rad` | Series QC (as in the JSON sidecar) |
| `markertype`, `source_file` | Raw file the series came from |
| `pose_x/y/z`, `beam_x/y/z` | Mean instrument position and beam direction (Localite planning-image mm) |
| `trans_ras_*`, `focus_ras_*` | Transducer (bowl rear centre) and focus [mm], via `localite_matrix_to_positions` |
| `trans_vox_*`, `focus_vox_*` | The same in planning-T1 voxels |
| `xml`, `status` | Ingested file and outcome |

Transducer, focus and voxel coordinates need the transducer focal distance (`transducer.focal_distance_ep` or `_bowl`) and the planning T1 (`path.localite_t1`, pattern `path.localite_t1_pattern`); when either is missing they are `NaN` and `status` says why. The pose columns are always filled.

The [multi-target toolkit](doc_placement_neuronav_multitarget.md) reads the ingested files for target resolution and MNI export; `placement.mode: 'localite'` reads them for simulation.

## GUI workflow

```matlab
neuronav_gui                   % empty form
neuronav_gui('study.yaml')     % pre-fill from a study config
neuronav_gui(parameters)       % pre-fill from a parameter struct
```

### Single session

1. Set the raw root (containing `<sub>/<ses>/<loc folder>/…`; leave **Localite folder** empty if the export sits directly under `<ses>/`) and the output root. The full layout is in [Localite I/O](doc_placement_neuronav.md#localite-io). The subject and session dropdowns list the `sub-*` / `ses-*` folders found; they can also be typed (a bare number `2` means `ses-02`).
2. **▶ Detect series** detects the series into a temporary folder; nothing is written yet. Check the series table and the plots.
3. Tick the series to keep, check each target name (and optionally set a `#rrggbb` colour), then **✓ Write outputs**. Changing any input after detection disables writing until series are detected again.

??? note "Details: settings read from a config, the series table and the pulse filter"
    Settings read from a config (all optional): `path.localite_raw` (falls back to `path.localite`), `path.localite_post`, `placement.localite.loc_folder` / `voxel_size` / `expected_segment_length` / `session`, `io.overwrite_files`, `neuronav.coil_map`.

    The series table lists every detected series with pulse count, position deviation (mm from the series mean position), rotation deviation (rad from the mean orientation) and a flag for untracked series (all positions at the origin). The tabs show a 3-D plot of every pulse position (small dot = one pulse, ring = series mean, arrow = beam direction; no arrow for GUMMarkers), per-pulse deviation traces (one dot per pulse), all raw marker files found with the one used marked, and the captured log/warnings.

    The **Pulses per series** filter (`expected_segment_length`, default 80) drops TriggerMarkers series whose pulse count differs from it by more than 20; leave it blank to keep all series.

### Target names

With **Suggest from planned targets** ticked (default), detection also names each series after the planned target it was aimed at. The **Name from** column shows the source and the off-axis distance and depth of each match (or the nearest miss), and the 3-D plot marks the candidates (×) and the axis to each matched target. Names typed into the coil table are kept. **↻ Re-apply to series** re-runs the naming after changing the tolerance or depth window.

A pose aimed at a target missing from the export takes the nearest remaining one, so check depth and distance: a match far outside the usual focal depth is suspect. Narrow **Target depth** around the focal distance to reject it. On the example data, all recorded poses whose planned target was in the export matched it, 0.7–4 mm off-axis.

??? note "Details: how series are matched to planned targets"
    The naming runs `neuronav_read_label_candidates` and `neuronav_suggest_target_names`:

    - Candidates are the named, set targets of the session's latest `EntryTarget` export and of the `EntryTargetPair`s in its GUMMarkers exports. Each label is kept once, at its most recent position.
    - A recorded pose (TriggerMarkers series, or a GUMMarkers `InstrumentMarker`) takes the candidate closest to its beam axis (the line from the mean position along matrix column 1), provided it lies within **Axis tolerance** (default 10 mm) and at a depth inside **Target depth** (default 0–150 mm).
    - A GUMMarkers `EntryTargetPair` series uses its own label; when a pose in the same file matched that label, the pair becomes `<label>_target` and the pose keeps the bare name. Repeated blocks on one target get `_rep2`, `_rep3`, ….

??? note "Details: single-coil and multi-coil"
    With the coil table empty, ingest uses the newest non-duplicate `TriggerMarkers_Coil0` file (falling back to `GUMMarkers`) and every detected series becomes a row. With one row per coil, ingest uses `neuronav.coil_map`: it reads one file per coil and keeps one series per target listed for that coil, skipping untracked blocks. For a coil that recorded several blocks, enter comma-separated target names (e.g. `L_MD, L_MD` → `L_MD_1`, `L_MD_2`).

### Many sessions

**▦ Batch…** opens a window that scans the raw root and lists every session with its marker files, planned labels and whether it is already ingested.

1. Tick sessions (**All**, **None**, **Not yet ingested**, or individually).
2. Check the label table on the right; group or rename targets.
3. Choose **Write** `all` or `named`, then **▶ Run batch**. Each session's status updates as it is processed; written targets are listed below and in the summary CSV.

**⬇ Save settings…** writes the settings (including the label map) as YAML for scripts or later GUI sessions.

The Batch window calls `neuronav_ingest_batch` with the main window's settings, so a saved settings file reproduces a GUI batch run in a script:

```matlab
p = load_parameters('neuronav_settings.yaml');
neuronav_ingest_batch(p, {'sub-001', 'ses-01'; 'sub-002', 'ses-01'});
```

### MNI view

**◎ MNI view…** (in the Batch window) maps the summary CSV to MNI space and plots transducer, focus and beam axis of every subject and session over an MNI head outline.

??? note "Details: loading, filtering and saving the MNI view"
    **MNI view** reads the summary CSV (default `<output root>/neuronav_ingest_summary.csv`, i.e. every session ingested so far) and, on **▶ Load & map to MNI**, maps each target's transducer, focus and Localite pose to MNI mm. Set **SimNIBS root** (`path.seg`); the planning T1 is found as in the main window (`path.localite_t1`). Subjects without a warp or planning T1 are listed in the status line and left out. The 3-D plot shows transducers (●), foci (◆), the beam axis between them and optionally the tracker pose (▲) over a head outline. Tick or untick subjects, sessions and targets in the three tables, colour by target, subject or session, and hover a point for its subject, session and target. **⬇ Save CSV…** writes the summary with the MNI columns.

??? note "Details: warp and head outline"
    The mapping uses `neuronav_summary_to_mni`: planning T1 → `m2m_<sub>/final_tissues.nii.gz` → `m2m_<sub>/toMNI/Conform2MNI_nonl.nii.gz`. The outline is taken from **Head template**: SimNIBS's `MNI152_T1_1mm.nii.gz` when a SimNIBS installation is found (`startup.simnibs_bin_path`, or the default install folder); left empty, it is taken from the first mapped subject's `final_tissues_MNI.nii.gz`.

??? note "Details: functions"
    | Function | Role |
    |---|---|
    | `neuronav_scan_sessions` | Inventory: sessions, marker files, planned labels, label table with stems |
    | `neuronav_label_map_from_scan` | Label table → `neuronav.target_labels` |
    | `neuronav_map_label` | Apply the label map to one label |
    | `neuronav_detect_series` | Detect the series of one session without writing |
    | `neuronav_name_series` | Naming rule (coil table > planned target > default) |
    | `neuronav_ingest_batch` | Loop: detect, name, write, summary CSV |
    | `neuronav_ingest_markers` | Write one session's ingested XML/JSON for a given target map |
    | `neuronav_summary_to_mni` | Add MNI coordinates (pose, transducer, focus) to the summary via the SimNIBS nonlinear warp |
