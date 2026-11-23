# Neuronavigation Read-In

![Pipeline steps; this page belongs to 3 Placement](assets/figures/fig_stage-strip-placement_diagram.svg){ .stage-strip }

PRESTUS reads transducer poses recorded with the [Localite](https://www.localite.de/) neuronavigation system and simulates them post hoc, to verify the focus of the realised stimulation. The recorded pose decides where the simulated transducer sits and where its beam points, so the file you read it from and the offset from the tracker to the transducer shape every focus position and intensity the simulation reports. An example demo (without data) is in `examples/demo_localite.m`.

On this page:

- [Placement from a Recorded Pose](#placement-from-a-recorded-pose): `placement.mode: 'localite'` and the focal distance it needs.
- [Realised Angle vs Planned Trajectory](#realised-angle-vs-planned-trajectory): which Localite export answers which question, and how the recorded pose becomes the simulated transducer and focus.

!!! note "At a glance"
    - **When you need it:** After a session, to simulate where the transducer actually was.
    - **Before you start:** Localite exports of the session and the planning T1 that was loaded into Localite.

Localite data pass through two stages. **Ingestion** (`neuronav_gui`, or `neuronav_ingest_markers` in a script; see [GUI and script workflows](doc_placement_neuronav_workflows.md)) reads the raw exports of one subject and session, averages each recorded series and writes one normalised file pair. **Placement** (`placement.mode: 'localite'`) reads that pair, or in file mode a single XML you point to.

Subpages:

- [Placing from a Recorded Pose](doc_placement_neuronav_recorded-pose.md): automatic and file mode, several transducers, many sessions; for anyone setting up a Localite run.
- [Localite Files, Offset, and Coordinates](doc_placement_neuronav_coordinates.md): folder layout and path keys, the tracker-to-bowl offset, and coordinate transforms; for anyone setting up paths, checking the transducer offset, or converting positions to voxels or MNI.
- [Localite Ingestion: GUI and Script Workflows](doc_placement_neuronav_workflows.md): turning raw exports into named, averaged targets in `neuronav_gui` or with `neuronav_ingest_batch`; for anyone ingesting sessions.
- [Named Target Coordinates](doc_placement_neuronav_targets.md): resolving several named targets per session from several sources in your own script; for studies with more than one target per subject.

---

## Placement from a Recorded Pose

Set `placement.mode: 'localite'`. The pose's column 4 (position) and column 1 (beam direction) become the transducer and focus positions.

!!! warning "Focal distance required"
    `transducer.focal_distance_ep` (or `focal_distance_bowl`) must be set for Localite placement; `config_default.yaml` ships it as `[]`. Other placement modes derive it from `trans_pos`/`focus_pos`, but here those positions are what the step produces. If neither is set, `localite_matrix_to_positions` raises `PRESTUS:localite:missingConfig` naming the missing field. The same check covers `placement.mode = 'plantus'`, which resolves its final position through the same function.

How to select the recorded pose (automatic or file mode, several transducers, many sessions) is described in [Placing from a Recorded Pose](doc_placement_neuronav_recorded-pose.md).

## Realised Angle vs Planned Trajectory

The Localite file types answer different questions. The wrong one silently mismatches the reported angle and the reported target.

- **The realised angle at the moment of stimulation** is only available from `TriggerMarkers` (the tracked pose during the pulses, averaged per series) or `InstrumentMarker` (a static saved pose). Both store a measured rotation.
- **`GUMMarkers` (`EntryTargetPair`) carries no measured angle.** Its 4×4 matrix is synthesised from the static `Target` point and a `RotationReference` vector, which is often unset in the source file (`<Rotation angle="0.0" used="false">`). With `used="false"`, the angle `neuronav_compute_series_statistics` produces is a placeholder. Treat it as a real angle only if the source XML has `Rotation used="true"`.
- **Static planned Entry/Target coordinates are never read by placement.** This includes the `Target` field of a `GUMMarkers` `EntryTargetPair` and a standalone `EntryTarget.xml`. For single-target placement they are optional: TriggerMarkers or InstrumentMarker alone is sufficient. They matter only for [target resolution](doc_placement_neuronav_targets.md#target-coordinates-from-several-sources), which uses them to identify which planned target a coil's series corresponds to (`neuronav_match_entry_targets`) and as a cross-check against pre-hoc planned coordinates (`neuronav_resolve_target`).
- **The reported focus position is axis-projected, not a stored coordinate.** `localite_matrix_to_positions` computes `focus_pos = trans_pos + focal_distance_bowl * ref_vec`, where `trans_pos = ref_pos + d * ref_vec` is the bowl rear centre and `d` is the resolved [transducer offset](doc_placement_neuronav_coordinates.md#transducer-offset-tracker_to_bowl_mm) (with the default exit-plane origin this equals `ref_pos + focal_distance_ep * ref_vec`): it projects forward from the tracked position along the tracked beam axis by the calibrated focal distance. This keeps the reported target self-consistent with the reported angle. A separately clicked static Target point need not lie on the tracked beam axis.

??? note "Localite XML file formats"
    All formats are parsed by `neuronav_compute_series_statistics`, which handles single markers, repeated trigger trains, and multi-position recordings alike.

    **TriggerMarkers.** Multiple triggers can be sent during stimulation to record the transducer location over time. Files are named `TMSTrigger/TriggerMarkers_Coil0_<timestamp>`. Recorded locations are averaged within each position series (identified by time-gap segmentation). Multiple positions per recording session are supported.

    **GUMMarkers.** One or more static marker positions. Each `Element` is expected to be of type `EntryTargetPair` (with a `Target/Marker/ColVec3D` position and a `Rotation/RotationReference/ColVec3D` beam direction), from which a synthetic 4×4 matrix is constructed. Some Localite versions write a `GUMMarkerList` file that mixes `EntryTargetPair` and `InstrumentMarker` elements (`type="InstrumentMarker"` on the `Element` tag). Each element becomes its own series: `InstrumentMarker` elements keep their recorded 4×4 pose, `EntryTargetPair` elements the synthesised matrix. Only the former are instrument poses; placement refuses an `EntryTargetPair` series (`PRESTUS:localite:notAPose`).

    **InstrumentMarker.** Static instrument positions stored as full 4×4 matrices in an `InstrumentMarkerList` file. Each `InstrumentMarker` entry contains a `Marker/Matrix4D` block with the same `data{row}{col}` field layout as TriggerMarkers. The ingested output of `neuronav_ingest_markers` uses this format.

    **EntryTarget** (not for position or angle). Some Localite versions export an `EntryTargetList` XML with flat `Target/Marker/ColVec3D` and `Rotation/RotationReference/ColVec3D` entries at the top level (not wrapped in `GUMMarkerList/Element`). `neuronav_compute_series_statistics` does not read it, and placement never uses it. `neuronav_read_entry_target` reads it for one purpose: matching an operator-confirmed brain-target label against a `target_map` in multi-target studies (see [Target coordinates from several sources](doc_placement_neuronav_targets.md#target-coordinates-from-several-sources)).
