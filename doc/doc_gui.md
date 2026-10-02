# PRESTUS GUI

The PRESTUS graphical user interface (`prestus_gui`) provides a point-and-click front-end to the pipeline. Every field maps directly to a YAML parameter — see [doc_parameters.md](doc_parameters.md) for the full reference.

> A standalone config builder with a placement viewer, the **PRESTUS Nexus**, is available as an experimental alternative: see [doc_nexus.md](doc_nexus.md).

![PRESTUS GUI](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/prestus_gui.png)

## Starting the GUI

```matlab
prestus_gui
```

Use **⬆ Load** to populate fields from an existing YAML, **⬇ Save** to export the current state, and **↺ Defaults** to reset to `config_default.yaml`.

### Configuration loading

The GUI and the pipeline share the same two-layer loading logic:

1. `config_default.yaml` is always loaded first and provides values for every parameter.
2. A study-specific YAML is merged on top via `MergeStruct` — only keys present in the study file overwrite the defaults; everything else retains its default value.

**In the GUI**, clicking **⬆ Load** calls `apply_params_to_gui` on the loaded file directly (no prior reset to defaults). This means:

- Fields present in the loaded YAML are updated.
- Fields **not** present in the loaded YAML keep whatever value is currently shown.

To get a clean state before loading, click **↺ Defaults** first, then **⬆ Load**.

**When running**, the GUI writes the current field values to a temporary YAML and passes it to `load_parameters`, which merges it on top of `config_default.yaml`. Fields left blank in the GUI are omitted from the saved YAML and therefore fall back to the default. The exception is `path.t2_pattern` — leaving it blank disables T2 entirely (`charm` is run with T1 only).

---

## macOS: Full Disk Access required

On macOS, MATLAB launched from the Finder or Dock runs under stricter sandboxing than MATLAB launched from a terminal. If full disk access is required, open **System Settings → Privacy & Security → Full Disk Access** and add `MATLAB.app`. A restart of MATLAB is required for the change to take effect. This does not affect terminal-launched MATLAB, where the error does not occur.

---

## Run

Press **▶ Run Simulation** to execute the pipeline with the current settings. The log area captures all terminal output produced during the run.

## Results

After a simulation completes, the **Results** tab shows orthogonal cross-sections of the output maps.

| Sub-tab | Content |
|---|---|
| **Acoustic** | X/Y/Z slices of the acoustic pressure / ISPPA map |
| **Thermal** | X/Y/Z slices of the temperature / CEM43 map |
| **Report** | Embedded HTML summary report |

Use **Refresh** to reload results, **Open HTML Report** to open the report in a browser, and **Open Output Folder** to reveal the output directory.

---

## Localite ingest GUI

`neuronav_gui` is a stand-alone window for Stage 1 of the neuronav pipeline (`neuronav_ingest_markers`, see [doc_placement_neuronav_multitarget.md](doc_placement_neuronav_multitarget.md)). It turns a raw Localite session folder into the normalised `InstrumentMarker_<sub>_<ses>.xml` and its QC JSON sidecar, with a review step before anything is written.

```matlab
neuronav_gui                % empty form
neuronav_gui('my_config.yaml') % pre-fill from a study config
neuronav_gui(parameters)     % pre-fill from a parameter struct
```

Settings read from a config (all optional): `path.localite_raw` (falls back to `path.localite`), `path.localite_post`, `placement.localite.loc_folder` / `voxel_size` / `expected_segment_length` / `session`, `io.overwrite_files`, `neuronav.coil_map`.

**Workflow**

1. Set the raw root (containing `<sub>/<ses>/<loc folder>/…`; leave **Localite folder** empty if the export sits directly under `<ses>/`) and the output root — the full layout is in [Localite I/O](doc_placement_neuronav.md#localite-io). The subject and session dropdowns list the `sub-*` / `ses-*` folders found; they can also be typed (a bare number `2` means `ses-02`).
2. **▶ Detect series** runs the ingest file selection and series statistics into a temporary folder — nothing is written to the output root. The series table then lists every detected series with pulse count, position deviation (mm from the series mean position), rotation deviation (rad from the mean orientation) and a flag for untracked series (all positions at the origin). The tabs show a 3-D plot of every pulse position (small dot = one pulse, ring = series mean, arrow = beam direction; no arrow for GUMMarkers), per-pulse deviation traces (one dot per pulse), all raw marker files found with the one actually used marked, and the captured log/warnings.
3. Tick the series to keep, check each target name (and optionally set a `#rrggbb` colour), then **✓ Write outputs**. Changing any input after detection disables writing until series are detected again.

**Target names.** With **Suggest from planned targets** ticked (default), detection also names each series after the planned target it was aimed at, via `neuronav_read_label_candidates` and `neuronav_suggest_target_names`:

- Candidates are the named, set targets of the session's latest `EntryTarget` export and of the `EntryTargetPair`s in its GUMMarkers exports. Each label is kept once, at its most recent position.
- A recorded pose (TriggerMarkers series, or a GUMMarkers `InstrumentMarker`) takes the candidate closest to its beam axis — the line from the mean position along matrix column 1 — provided it lies within **Axis tolerance** (default 10 mm) and at a depth inside **Target depth** (default 0–150 mm).
- A GUMMarkers `EntryTargetPair` series uses its own label; when a pose in the same file matched that label, the pair becomes `<label>_target` and the pose keeps the bare name. Repeated blocks on one target get `_rep2`, `_rep3`, ….
- The **Name from** column shows the source and the off-axis distance and depth of each match (or the nearest miss), and the 3-D plot marks the candidates (×) and the axis to each matched target. Names typed into the coil table are kept. **↻ Re-apply to series** re-runs the naming after changing the tolerance or depth window.

On the example data, all recorded poses whose planned target was in the export matched it, 0.7–4 mm off-axis. A pose aimed at a target missing from the export takes the nearest remaining one, so check depth and distance: a match far outside the usual focal depth is suspect — narrow **Target depth** around the focal distance to reject it.

**Single-coil vs. multi-coil.** With the coil table empty, ingest uses the newest non-duplicate `TriggerMarkers_Coil0` file (falling back to `GUMMarkers`) and every detected series becomes a row. With one row per coil, ingest uses `neuronav.coil_map`: it reads one file per coil and keeps one series per target listed for that coil, skipping untracked blocks. For a coil that recorded several blocks, enter comma-separated target names (e.g. `L_MD, L_MD` → `L_MD_1`, `L_MD_2`).

**Settings and batch.** **⬆ Load settings…** / **⬇ Save settings…** read and write the settings as a PRESTUS YAML config (paths, ingestion, coil map, naming, label map), so scripts can use the same file. **▦ Batch…** opens a window to scan all subjects and sessions, tick the ones to ingest (**All**, **None**, **Not yet ingested**), edit the map from recorded planned-target labels to study target names, and run `neuronav_ingest_batch` on them. Its **◎ MNI view…** shows the written targets of all subjects and sessions in MNI space (see below). See [GUI and script workflows](doc_placement_neuronav_workflows.md).

**MNI view.** Reads the summary CSV (default `<output root>/neuronav_ingest_summary.csv`, i.e. every session ingested so far) and, on **▶ Load & map to MNI**, maps each target's transducer, focus and Localite pose to MNI mm with `neuronav_summary_to_mni`: planning T1 → `m2m_<sub>/final_tissues.nii.gz` → `m2m_<sub>/toMNI/Conform2MNI_nonl.nii.gz`. Set **SimNIBS root** (`path.seg`); the planning T1 is found as in the main window (`path.localite_t1`). Subjects without a warp or planning T1 are listed in the status line and left out. The 3-D plot shows transducers (●), foci (◆), the beam axis between them and optionally the tracker pose (▲) over a head outline. The outline is taken from **Head template**: SimNIBS's `MNI152_T1_1mm.nii.gz` when a SimNIBS installation is found (`startup.simnibs_bin_path`, or the default install folder); left empty, it is taken from the first mapped subject's `final_tissues_MNI.nii.gz`. Tick or untick subjects, sessions and targets in the three tables, colour by target, subject or session, and hover a point for its subject, session and target. **⬇ Save CSV…** writes the summary with the MNI columns.

The **Pulses per series** filter (`expected_segment_length`, default 80) drops TriggerMarkers series whose pulse count differs from it by more than 20; leave it blank to keep all series.
