# Multi-Target Extraction Toolkit

!!! note "At a glance"
    - **What:** Tools to extract several planned targets from Localite exports in one go.
    - **When you need it:** Advanced: for studies with several targets per session. Most users need only [Localite ingestion](doc_placement_neuronav_workflows.md).
    - **Before you start:** Localite exports with planned targets.

This is a toolkit, not a pipeline mode: a set of functions your own extraction script calls to resolve **named, per-session brain-target coordinates** (e.g. this subject's `L_PUL` on `ses-03`) from whichever raw Localite export the operator recorded, with a logged account of how confident each resolved coordinate is. It sits on top of [Neuronavigation read-in](doc_placement_neuronav.md) (read that first for the Localite file formats and coordinate spaces) and adds to the single-target placement pipeline rather than replacing it. Use it when a study needs **more than one named target per subject**, or a documented answer to where each reported coordinate came from.

## What you need to provide

1. **Raw entry-point recordings** per subject and session, in any supported format (`TriggerMarkers`, `GUMMarkers` or `InstrumentMarker`). This supplies the transducer/coil position, always read from the current session's own recording.
2. **A `target_map`** naming your targets, built by `neuronav_build_target_map` from a `coil_map` you construct per subject and session:

    ```matlab
    coil_map(1).coil   = 0;              % hardware coil index (0-based)
    coil_map(1).target = 'L_PUL';        % target_name -- your convention
    coil_map(2).coil   = 1;
    coil_map(2).target = 'R_PUL';
    target_map = neuronav_build_target_map(coil_map);
    parameters.neuronav.coil_map = coil_map;   % ingest reads one file per coil
    ```

    A coil that records several position blocks in one session lists them as `coil_map(i).sequences(j).target`; ingest then takes that many series from the coil's file, in order (untracked blocks with every pulse at the origin are skipped), matching the `series_index` that `neuronav_build_target_map` assigns.

The following inputs are optional. Each one you provide moves more targets up the confidence ranking (see [Resolution order](#resolution-order)).

- **`EntryTarget` exports**, for at least some sessions, if your operator uses Localite's target-confirmation feature. `target_map(i).localite_key` is matched (case-insensitive substring) against each `EntryTarget` description. Set it whenever your description convention differs from `target_name`, and pass a **cell array of candidates** if more than one convention is in use across sessions:

    ```matlab
    coil_map(1).localite_key = {'LEFT PUL', 'L_PUL', 'LEFT P', 'L_P'};
    ```

- **A pre-hoc target lookup table**, if your study establishes targets before the session (e.g. a placement or optimisation step run ahead of time). This is the highest-confidence source when present and is cross-checked against `EntryTarget` wherever both exist. `neuronav_read_target_lookup.m` documents the CSV schema (`subject_id, target_name, space, x, y, z`). Build this file with your own per-study script from whatever your planning step produces.

??? note "Details: why several sources"
    Across one study, the same target was confirmed in different ways for different subjects and sessions: an operator-confirmed `EntryTarget` export for the exact session; an `EntryTarget` export only from an earlier session; no `EntryTarget` export at all (workflows that only used `GUMMarkers`, which has no separate target confirmation); two label conventions in the `EntryTarget` description across sessions (`"LEFT P"` and `"L_P"`); and a pre-hoc planned coordinate from individual anatomical segmentation, independent of Localite. The toolkit resolves each target from the best source available for that target and reports which one it used.

## Resolution order

`neuronav_resolve_target` tries these in order per target and returns which one supplied each row (`target_source`). **Check this before treating a coordinate as ground truth.**

| Tier | `target_source` value | What it is | Confidence |
|---|---|---|---|
| 1 | `target_lookup` | Your pre-hoc lookup table (see above) | Highest: established before the session |
| 2 | `entrytarget_current` | This session's own `EntryTarget` export, label-matched | High: operator-confirmed, same day |
| 3 | `entrytarget:<session>` | Another session's `EntryTarget` export (nearest prior session tried first) | Medium: operator-confirmed, different day |
| 4 | `projected` | Transducer-axis geometric projection from the tracked entry pose | Estimate only: no operator confirmation exists for this target |

Because of tier 4, the resolved coordinate is **never NaN**: every target gets a value, but tier-4 rows are not confirmed. In one cross-check, tier-4 estimates differed from the tier-1 target by 3-5 voxels for a subject without any `EntryTarget` export, while tier-1 and tier-2 values matched to within 0-1 voxel. If a study needs sub-voxel target accuracy and some subjects stay on tier 4, provide a lookup table (tier 1) for those subjects.

## Minimal worked example

```matlab
for i_sub = 1:numel(all_subjects)
    sub_id = all_subjects{i_sub};
    for i_ses = 1:numel(sessions)
        session = sessions{i_ses};

        coil_map = build_my_coil_map(sub_id, session);   % your own logic
        target_map = neuronav_build_target_map(coil_map);
        parameters.neuronav.coil_map = coil_map;

        [xml_path, ~, results] = neuronav_ingest_markers(parameters, sub_id, session, target_map);
        active_sides = find(~cellfun(@isempty, results));

        [trans_ras, trans_pos, target_ras, target_pos_projected] = ...
            neuronav_convert_trigger_to_voxels(sub_id, active_sides, results, parameters);

        t1_header = niftiinfo(my_planning_t1_path(sub_id));   % see doc_placement_neuronav.md

        [target_vox, target_mm, target_source] = neuronav_resolve_target( ...
            parameters, sub_id, session, sessions, target_map, t1_header, ...
            target_pos_projected, my_target_lookup_csv_path);

        % target_vox is what you feed into PRESTUS as focus_pos.
        % target_source tells you, per row, whether that's confirmed or estimated.

        neuronav_export_session_csv(sub_id, session, {target_map.target_name}, ...
            outputs_folder, trans_pos, target_vox, trans_ras, target_ras, ...
            [], [], [], [], 0, target_mm, target_source);
    end
end
```

## Output: the exported CSV

`neuronav_export_session_csv` writes one row per target with, among other columns:

- `trans_native_vox_x/y/z`: transducer/coil position, always from the current session's own recording.
- `target_native_vox_x/y/z`: the resolved target (`target_vox` above). **A simulation script reads this as `focus_pos`.**
- `target_localite_mm_x/y/z`: the `EntryTarget`-matched mm coordinate, tiers 2-3 only. NaN for tier 1 and tier 4 rows, so it stays a reliable operator-confirmed signal independent of `target_source`.
- `target_projected_mm_x/y/z`: the tier-4 geometric estimate, always computed whichever tier was used, kept for comparison.
- `target_source`: see the resolution-order table above.

??? note "Details: function reference"
    | Function | Role |
    |---|---|
    | `neuronav_ingest_markers` | Stage 1: reads raw markers (any format), writes normalised `InstrumentMarker` XML + a QC JSON sidecar (pulse count, position deviation in mm, rotation deviation in rad, source file) per subject/session |
    | `neuronav_gui` | Stand-alone GUI for the ingest step: preview detected series and QC, name targets, write the XML/JSON (see [doc_gui.md](doc_gui.md#localite-ingest-gui)) |
    | `neuronav_ingest_batch` / `neuronav_scan_sessions` | Batch ingestion and inventory of many subjects/sessions, with a summary CSV of coordinates; see [GUI and script workflows](doc_placement_neuronav_workflows.md) |
    | `neuronav_read_label_candidates` / `neuronav_suggest_target_names` | Heuristic naming: collects named planned targets (EntryTarget, GUMMarkers pairs) and names each ingested series after the one nearest its beam axis, with distance/depth for QC |
    | `neuronav_build_target_map` | Expands your `coil_map` into the flat `target_map` every other function here consumes |
    | `neuronav_convert_trigger_to_voxels` | Converts ingested marker positions to native voxel/mm space; also computes the tier-4 axis-projected target |
    | `neuronav_read_entry_target` | Reads the latest `EntryTarget*.xml` snapshot for one subject/session |
    | `neuronav_match_entry_targets` | Label-matches `target_map` entries against `EntryTarget` descriptions (tiers 2-3's matching step) |
    | `neuronav_read_target_lookup` | Reads your pre-hoc lookup CSV (tier 1) |
    | **`neuronav_resolve_target`** | **The orchestrator**: runs all four tiers in order for one subject/session and returns the resolved coordinate + provenance per target. Call this one function rather than the matching/lookup functions individually. |
    | `neuronav_export_session_csv` | Writes the per-session summary CSV, including a `target_source` column |
    | `neuronav_plot_placement` | QC plot of every target_map entry's entry/target pair on the segmented head |
    | `neuronav_get_group_mean_mni` / `neuronav_apply_inverse_deformation` | Group-mean-MNI fallback for subjects missing trigger data entirely (a more severe gap than a missing target confirmation; see the coordinate-transform section of [Neuronavigation read-in](doc_placement_neuronav.md)) |
