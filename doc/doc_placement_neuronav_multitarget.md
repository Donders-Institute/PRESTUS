# Multi-target extraction toolkit

This is a toolkit, not a pipeline mode: a set of functions your own
extraction script calls to resolve **named, per-session brain-target
coordinates** (e.g. "this subject's `L_PUL` on `ses-03`") from whichever
raw Localite export the operator happened to record, with an explicit,
logged account of how confident each resolved coordinate actually is.

It sits on top of everything in [doc_placement_neuronav.md](doc_placement_neuronav.md)
(read that first for the underlying Localite file formats and coordinate
spaces) and is layered *in addition to*, not instead of, the single-target
placement pipeline described there. Use this toolkit when a study needs
**more than one named target per subject**, or needs a documented,
auditable answer to "where did this coordinate actually come from" for
every target it reports.

## Why this exists

A single subject's stimulation targets are rarely all confirmed the same
way. Across one real study, all of the following turned up for different
subjects/sessions of the *same* target:

- An operator-confirmed `EntryTarget` export, present for this exact session.
- An `EntryTarget` export, but only from an earlier session (the operator
  didn't re-export it every visit).
- No `EntryTarget` export at all for any session (some Localite/operator
  workflows only ever used `GUMMarkers`, which has no separate target
  confirmation — see doc_placement_neuronav.md).
- Two different label conventions in the `EntryTarget` description text
  across sessions for the same subject (`"LEFT P"` vs `"L_P"`).
- A pre-hoc planned coordinate, established before the session via
  individual anatomical segmentation, independent of anything Localite
  recorded.

None of this is unusual — it's what "one study, multiple RAs, multiple
software versions, multiple months" looks like in practice. The toolkit's
job is to resolve each target from the best of these sources available
*for that specific target*, and tell you which one it used.

## What you need to provide

Minimum to use this toolkit at all:

1. **Raw entry-point recordings** per subject/session, in any format
   doc_placement_neuronav.md supports (`TriggerMarkers`, `GUMMarkers`, or
   `InstrumentMarker`). This is what supplies the transducer/coil position
   — always read from the current session's own recording, never
   estimated or borrowed.
2. **A `target_map`** naming your targets, built by `neuronav_build_target_map`
   from a `coil_map` you construct per subject/session:
   ```matlab
   coil_map(1).coil   = 0;              % hardware coil index (0-based)
   coil_map(1).target = 'L_PUL';        % target_name -- your convention
   coil_map(2).coil   = 1;
   coil_map(2).target = 'R_PUL';
   target_map = neuronav_build_target_map(coil_map);
   ```

Everything below is optional and additive — each one you provide moves
more targets up the confidence ranking (see "Resolution order").

3. **`EntryTarget` exports**, for at least some sessions, if your operator
   uses Localite's target-confirmation feature. `target_map(i).localite_key`
   is what gets matched (case-insensitive substring) against each
   `EntryTarget` description; set it explicitly whenever your description
   convention doesn't match `target_name` exactly, and pass a **cell array
   of candidates** if more than one convention is in use across sessions
   (this is common — see "Why this exists" above):
   ```matlab
   coil_map(1).localite_key = {'LEFT PUL', 'L_PUL', 'LEFT P', 'L_P'};
   ```
4. **A pre-hoc target lookup table**, if your study establishes targets
   before the session (e.g. via a placement/optimization step run ahead of
   time). This is the highest-confidence source when present — not a
   post-hoc reconstruction at all — and is cross-checked against
   `EntryTarget` wherever both exist. See `neuronav_read_target_lookup.m`
   for the exact CSV schema (`subject_id, target_name, space, x, y, z`).
   Build this file with your own small per-study script that converts
   whatever your planning step already produces into that schema; nothing
   in PRESTUS needs to know your planning step's own file format.

## Resolution order

`neuronav_resolve_target` tries these in order per target, and returns
which one actually supplied each row (`target_source`) — **always check
this before treating a coordinate as ground truth**:

| Tier | `target_source` value | What it is | Confidence |
|---|---|---|---|
| 1 | `target_lookup` | Your pre-hoc lookup table (§4 above) | Highest — established before the session |
| 2 | `entrytarget_current` | This session's own `EntryTarget` export, label-matched | High — operator-confirmed, same-day |
| 3 | `entrytarget:<session>` | Another session's `EntryTarget` export (nearest prior session tried first) | Medium — operator-confirmed, different day |
| 4 | `projected` | Transducer-axis geometric projection from the tracked entry pose | Estimate only — no operator confirmation exists anywhere for this target |

Tier 4 means `target_source`/the resolved coordinate is **never NaN** —
every target always gets *some* value — but tier 4 rows should not be
treated as confirmed. In one real cross-check, tier-4 estimates differed
from the true (tier-1) target by 3-5 voxels for a subject with no
`EntryTarget` export anywhere, while tier-1/tier-2 values matched to
within 0-1 voxel. If a study cares about sub-voxel target accuracy and has
any subjects permanently stuck on tier 4, that's a signal to provide a
lookup table (tier 1) for those subjects rather than trust the estimate.

## Function reference

| Function | Role |
|---|---|
| `neuronav_ingest_markers` | Stage 1: reads raw markers (any format), writes normalized `InstrumentMarker` XML + a QC JSON sidecar (pulse count, position/rotation deviation) per subject/session |
| `neuronav_build_target_map` | Expands your `coil_map` into the flat `target_map` every other function here consumes |
| `neuronav_convert_trigger_to_voxels` | Converts ingested marker positions to native voxel/mm space; also computes the tier-4 axis-projected target |
| `neuronav_read_entry_target` | Reads the latest `EntryTarget*.xml` snapshot for one subject/session |
| `neuronav_match_entry_targets` | Label-matches `target_map` entries against `EntryTarget` descriptions (tiers 2-3's matching step) |
| `neuronav_read_target_lookup` | Reads your pre-hoc lookup CSV (tier 1) |
| **`neuronav_resolve_target`** | **The orchestrator** — runs all four tiers in order for one subject/session and returns the resolved coordinate + provenance per target. Call this one function rather than the matching/lookup functions individually. |
| `neuronav_export_session_csv` | Writes the per-session summary CSV, including a `target_source` column |
| `neuronav_plot_placement` | QC plot of every target_map entry's entry/target pair on the segmented head |
| `neuronav_get_group_mean_mni` / `neuronav_apply_inverse_deformation` | Group-mean-MNI fallback for subjects missing trigger data entirely (a different, more severe gap than a missing target confirmation — see doc_placement_neuronav.md's coordinate-transform section) |

## Minimal worked example

```matlab
for i_sub = 1:numel(all_subjects)
    sub_id = all_subjects{i_sub};
    for i_ses = 1:numel(sessions)
        session = sessions{i_ses};

        coil_map = build_my_coil_map(sub_id, session);   % your own logic
        target_map = neuronav_build_target_map(coil_map);

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

`neuronav_export_session_csv` writes one row per target with, among other
columns:

- `trans_native_vox_x/y/z` — transducer/coil position, always from the
  current session's own recording.
- `target_native_vox_x/y/z` — the resolved target (`target_vox` above).
  **This is what a simulation script should read as `focus_pos`.**
- `target_localite_mm_x/y/z` — the `EntryTarget`-matched mm coordinate,
  tiers 2-3 only. Deliberately left NaN for tier 1 and tier 4 rows, so it
  stays a reliable "was this actually operator-confirmed" signal
  independent of `target_source`.
- `target_projected_mm_x/y/z` — the tier-4 geometric estimate, always
  computed regardless of which tier ended up being used, kept as a
  reference for comparison.
- `target_source` — see the resolution-order table above.
