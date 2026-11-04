# MNI Transducer Placement

!!! note "At a glance"
    - **What:** Placement from a transducer position and a target given in MNI coordinates.
    - **When you need it:** When you know both points in template space, e.g. a group target and a standard scalp position.
    - **Before you start:** A SimNIBS segmentation of the subject (for the MNI-to-subject transforms).

`placement.mode: 'mni'` places the transducer at a position **you** give in MNI (MNI152, mm) coordinates, together with the focus. No search is performed; the result is deterministic. The [heuristic](doc_placement_heuristic.md) mode takes only an MNI target and searches the scalp for a transducer position. The `mni` mode applies the same scalp standoff geometry, so the housing clears the skin.

```yaml
placement:
  mode: 'mni'
  mni:
    trans_pos_mm: [-52,  -8,  58]   # desired transducer (scalp entry) in MNI mm
    focus_pos_mm: [-29,  -6,   7]   # focus / target in MNI mm
    skin_gap_mm:  5                 # standoff from scalp along the focal axis [mm]
```

A SimNIBS `m2m_sub-XXX` folder (under `path.seg`) is required: the conversion relies on `final_tissues.nii.gz` and the SimNIBS coordinate transforms.

## How a position is resolved

1. **Coordinate conversion.** Both points are mapped MNI → subject RAS+ → grid voxel, each with the SimNIBS transform suited to where it sits:

    | Point | Transform | Why |
    |---|---|---|
    | `focus_pos_mm` | **nonlinear** (`nonl`) | A deep in-brain target; the nonlinear warp is most accurate here, identical to the heuristic target. |
    | `trans_pos_mm` | **linear** (`12dof` affine) | The transducer is outside the brain (scalp), where the nonlinear deformation field is extrapolated and unreliable. The affine is defined everywhere. |

2. **Snap to scalp.** The converted transducer voxel is snapped to the nearest scalp voxel, so a slightly off MNI entry coordinate still lands on the head surface.

3. **Standoff.** The transducer is offset outward from the scalp along the focal axis (target → scalp) by `skin_gap_mm + dist_tp_to_ep`, where `dist_tp_to_ep` is the bowl-to-exit-plane distance derived from the transducer geometry. The default `skin_gap_mm = 5` reproduces the heuristic's gap.

The resolved `trans_pos` / `focus_pos` (T1 grid voxels) are written into every configured transducer, and a QC overlay is saved. Downstream simulation is unchanged.

??? note "Details: functions involved"
    - Conversion: `transform_coordinates` and `mni2subject_coords_LDfix`, shared with the heuristic.
    - Scalp snap: `tp_scalp_boundary` (nearest scalp / outer-boundary voxel).
    - Standoff: `transducer_scalp_geometry`, shared with `tp_candidate_mesh`; the geometry the heuristic applies to every candidate.
    - QC overlay: `plot_placement_t1_overlay`.

## Limitations

- The standoff geometry assumes annular-style transducer fields (`curv_radius_mm`, `elem_od_mm`), as the heuristic does. The mode is intended for annular (CTX-style) transducers.
- `trans_pos_mm` is a desired scalp entry. It is snapped to the scalp, so the converted point need not lie exactly on the surface.
- Check the QC overlay. Converting a scalp coordinate through any registration is approximate.

See also: [Transducer Placement overview](doc_placement.md), [Heuristic placement](doc_placement_heuristic.md), [Coordinate systems](doc_coordinate_systems.md).
