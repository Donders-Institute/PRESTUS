# Named Target Coordinates

When a session has several named targets, each simulated focus needs a coordinate and a record of how far that coordinate can be trusted. `neuronav_resolve_target` takes each target's coordinate from the best source available, from a pre-hoc lookup table down to a projection along the tracked transducer axis, and reports which source it used.

On this page:

- [Target Coordinates from Several Sources](#target-coordinates-from-several-sources): what the functions resolve and when to use them.
- [Inputs You Provide](#inputs-you-provide): the recordings, the `target_map` and the optional sources.
- [Resolution Order](#resolution-order): the four tiers and their confidence.
- [Worked Example](#worked-example): a loop over subjects and sessions.
- [Exported Session CSV](#exported-session-csv): the per-session output and the function reference.

!!! note "At a glance"
    - **When you need it:** A study has more than one named target per subject and session, or needs a documented answer to where each target coordinate came from.
    - **Before you start:** The raw Localite exports and the planning T1 per subject; ingestion as in [Localite Ingestion](doc_placement_neuronav_workflows.md).

---

## Target Coordinates from Several Sources

When a session has several named targets, your own extraction script resolves **named, per-session brain-target coordinates** (e.g. this subject's `L_V1` on `ses-03`) from whichever raw Localite export the operator recorded, with a logged account of how confident each resolved coordinate is. The resolved coordinate is what a simulation uses as `focus_pos`. This is not a pipeline mode: the functions read the ingested files and add to the single-target placement pipeline rather than replacing it. Use them when a study needs **more than one named target per subject**, or a documented answer to where each reported coordinate came from. Localite file formats are in [Neuronavigation read-in](doc_placement_neuronav.md#realised-angle-vs-planned-trajectory) and coordinate spaces in [Localite Files, Offset and Coordinates](doc_placement_neuronav_coordinates.md#coordinate-transforms).

## Inputs You Provide

1. **Raw entry-point recordings** per subject and session, in any supported format (`TriggerMarkers`, `GUMMarkers` or `InstrumentMarker`). This supplies the transducer/coil position, always read from the current session's own recording.
2. **A `target_map`** naming your targets, built by `neuronav_build_target_map` from the `neuronav.coil_map` you construct per subject and session:

    ```matlab
    coil_map(1).coil   = 0;              % hardware coil index (0-based)
    coil_map(1).target = 'L_V1';        % target_name -- your convention
    coil_map(2).coil   = 1;
    coil_map(2).target = 'R_V1';
    target_map = neuronav_build_target_map(coil_map);
    parameters.neuronav.coil_map = coil_map;   % ingest reads one file per coil
    ```

    A coil that records several position blocks in one session lists them as `coil_map(i).sequences(j).target`; ingest then takes that many series from the coil's file, in order (untracked blocks with every pulse at the origin are skipped), matching the `series_index` that `neuronav_build_target_map` assigns.

The following inputs are optional. Each one you provide moves more targets up the confidence ranking (see [Resolution order](#resolution-order)).

- **`EntryTarget` exports**, for at least some sessions, if your operator uses Localite's target-confirmation feature. `target_map(i).localite_key` is matched (case-insensitive substring) against each `EntryTarget` description. Set it whenever your description convention differs from `target_name`, and pass a **cell array of candidates** if more than one convention is in use across sessions:

    ```matlab
    coil_map(1).localite_key = {'LEFT V1', 'L_V1', 'LEFT V', 'L_V'};
    ```

- **A pre-hoc target lookup table**, if your study establishes targets before the session (e.g. a placement or optimisation step run ahead of time). This is the highest-confidence source when present and is cross-checked against `EntryTarget` wherever both exist. `neuronav_read_target_lookup.m` documents the CSV schema (`subject_id, target_name, space, x, y, z`). Build this file with your own per-study script from whatever your planning step produces.

??? note "Why several sources"
    Across one study, the same target was confirmed in different ways for different subjects and sessions: an operator-confirmed `EntryTarget` export for the exact session; an `EntryTarget` export only from an earlier session; no `EntryTarget` export at all (workflows that only used `GUMMarkers`, which has no separate target confirmation); two label conventions in the `EntryTarget` description across sessions (`"LEFT V"` and `"L_V"`); and a pre-hoc planned coordinate from individual anatomical segmentation, independent of Localite. `neuronav_resolve_target` resolves each target from the best source available for that target and reports which one it used.

## Resolution Order

`neuronav_resolve_target` tries these in order per target and returns which one supplied each row (`target_source`). **Check this before treating a coordinate as ground truth.**

| Tier | `target_source` value | What it is | Confidence |
|---|---|---|---|
| 1 | `target_lookup` | Your pre-hoc lookup table (see above) | Highest: established before the session |
| 2 | `entrytarget_current` | This session's own `EntryTarget` export, label-matched | High: operator-confirmed, same day |
| 3 | `entrytarget:<session>` | Another session's `EntryTarget` export (nearest prior session tried first) | Medium: operator-confirmed, different day |
| 4 | `projected` | Transducer-axis geometric projection from the tracked entry pose | Estimate only: no operator confirmation exists for this target |

Because of tier 4, the resolved coordinate is **never NaN**: every target gets a value, but tier-4 rows are not confirmed. In one cross-check, tier-4 estimates differed from the tier-1 target by 3-5 voxels for a subject without any `EntryTarget` export, while tier-1 and tier-2 values matched to within 0-1 voxel. If a study needs sub-voxel target accuracy and some subjects stay on tier 4, provide a lookup table (tier 1) for those subjects.

## Worked Example

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

        t1_header = niftiinfo(my_planning_t1_path(sub_id));   % see doc_placement_neuronav_coordinates.md

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

## Exported Session CSV

`neuronav_export_session_csv` writes one row per target with, among other columns:

- `trans_native_vox_x/y/z`: transducer/coil position, always from the current session's own recording.
- `target_native_vox_x/y/z`: the resolved target (`target_vox` above). **A simulation script reads this as `focus_pos`.**
- `target_localite_mm_x/y/z`: the `EntryTarget`-matched mm coordinate, tiers 2-3 only. NaN for tier 1 and tier 4 rows, so it stays a reliable operator-confirmed signal independent of `target_source`.
- `target_projected_mm_x/y/z`: the tier-4 geometric estimate, always computed whichever tier was used, kept for comparison.
- `target_source`: see the resolution-order table above.

??? developer "Function reference"
    The ingestion functions (`neuronav_ingest_markers`, `neuronav_ingest_batch`, `neuronav_scan_sessions`) are listed under [GUI workflow](doc_placement_neuronav_workflows.md#gui-workflow); the naming functions under [Series naming](doc_placement_neuronav_workflows.md#series-naming).

    | Function | Role |
    |---|---|
    | `neuronav_build_target_map` | Expands your `coil_map` into the flat `target_map` every other function here consumes |
    | `neuronav_convert_trigger_to_voxels` | Converts ingested marker positions to native voxel/mm space; also computes the tier-4 axis-projected target |
    | `neuronav_read_entry_target` | Reads the latest `EntryTarget*.xml` snapshot for one subject/session |
    | `neuronav_match_entry_targets` | Label-matches `target_map` entries against `EntryTarget` descriptions (tiers 2-3's matching step) |
    | `neuronav_read_target_lookup` | Reads your pre-hoc lookup CSV (tier 1) |
    | **`neuronav_resolve_target`** | **The orchestrator**: runs all four tiers in order for one subject/session and returns the resolved coordinate + provenance per target. Call this one function rather than the matching/lookup functions individually. |
    | `neuronav_export_session_csv` | Writes the per-session summary CSV, including a `target_source` column |
    | `neuronav_plot_placement` | QC plot of every target_map entry's entry/target pair on the segmented head |
    | `neuronav_get_group_mean_mni` / `neuronav_apply_inverse_deformation` | Group-mean-MNI fallback for subjects missing trigger data entirely (a more severe gap than a missing target confirmation; see [Coordinate transforms](doc_placement_neuronav_coordinates.md#coordinate-transforms)) |
