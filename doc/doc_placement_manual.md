# Manual Transducer Placement

![Pipeline steps; this page belongs to 3 Placement](assets/figures/fig_stage-strip-placement_diagram.svg){ .stage-strip }

When the transducer and target positions are already known, for example from an fMRI hotspot, an anatomical marker, or a previous session, no search is needed. `placement.mode: 'manual'` takes both points from the config and uses them exactly as given.

On this page:

- [Entering Coordinates](#entering-coordinates): positions in voxels, RAS mm, or MNI mm, and how they are converted.

!!! note "At a glance"
    - **When you need it:** You already know where the transducer sits and where it aims.
    - **Before you start:** The T1w planning image; for MNI coordinates also the SimNIBS `m2m_sub-XXX` folder.

## Entering Coordinates

By default the coordinates are 1-based T1 voxel indices. Set `placement.manual.coordinate_space` to enter millimetres instead:

```yaml
placement:
  mode: 'manual'
  manual:
    coordinate_space: 'ras'   # 'voxel' (default) | 'ras' | 'mni'
transducer:
  trans_pos: [-62.1, 12.4, 48.0]   # bowl apex, subject RAS+ mm
  focus_pos: [-21.5, 10.2, 2.3]
```

- `ras`: subject RAS+ mm of the T1 (read from a viewer or neuronavigation), converted with `ras_to_grid()` via the T1 affine.
- `mni`: MNI mm. Both points are converted with charm's nonlinear warp, as in [MNI placement](doc_placement_mni.md), but the transducer is **not** snapped to the scalp, and the warp is less exact near the scalp: check the placement figure. The points are used exactly as given. Requires the SimNIBS `m2m_sub-XXX` folder.

Conversion happens once per transducer in the placement stage, so `modules.run_transducer_placement` must be `1`.

After conversion, the positions become voxel indices and `coordinate_space` is set to `'voxel'`. The original inputs are kept in `placement.manual.input_trans_pos` / `input_focus_pos`. If neither `focal_distance_ep` nor `focal_distance_bowl` is set, the distance fallback from the positions is computed after conversion.
