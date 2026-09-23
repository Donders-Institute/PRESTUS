# MNI Transducer Placement

The `mni` placement mode lets you specify **both** the transducer position and
the focus directly in MNI (MNI152, mm) coordinates. Unlike the
[heuristic](doc_placement_heuristic.md) mode — which takes only an MNI *target*
and then *searches* the scalp for a transducer position — `mni` mode places the
transducer at the location **you** provide, applying the same scalp standoff
geometry the heuristic uses so the housing clears the skin. No search is
performed: the result is deterministic.

Enable it with:

```yaml
placement:
  mode: 'mni'
  mni:
    trans_pos_mm: [-52,  -8,  58]   # desired transducer (scalp entry) in MNI mm
    focus_pos_mm: [-29,  -6,   7]   # focus / target in MNI mm
    skin_gap_mm:  5                 # standoff from scalp along the focal axis [mm]
```

A SimNIBS `m2m_sub-XXX` folder (under `path.seg`) is required — the conversion
relies on `final_tissues.nii.gz` and the SimNIBS coordinate transforms.

## How a position is resolved

1. **Coordinate conversion (per-point transform type).** Both points are mapped
   `MNI → subject RAS+ → grid voxel`, but with a different SimNIBS transform for
   each point, chosen by where the point sits:

   | Point | Transform | Why |
   |---|---|---|
   | `focus_pos_mm` | **nonlinear** (`nonl`) | A deep in-brain target; the nonlinear warp is most accurate here — identical to the heuristic target. |
   | `trans_pos_mm` | **linear** (`12dof` affine) | The transducer is *outside the brain* (scalp), where the nonlinear deformation field is unreliable/extrapolated. The affine is well-defined everywhere. |

   The conversion reuses the same code as the heuristic
   (`transform_coordinates`, `mni2subject_coords_LDfix`).

2. **Snap to scalp.** The converted transducer voxel is snapped to the nearest
   scalp / outer-boundary voxel (`tp_scalp_boundary`), so a slightly off MNI
   entry coordinate still lands on the head surface.

3. **Standoff geometry.** The transducer is offset outward from the scalp along
   the focal axis (target → scalp) by `skin_gap_mm + dist_tp_to_ep`, where
   `dist_tp_to_ep` is the bowl-to-exit-plane distance derived from the transducer
   geometry. This is the exact geometry the heuristic applies to every candidate
   (`transducer_scalp_geometry`, shared with `tp_candidate_mesh`). The default
   `skin_gap_mm = 5` reproduces the heuristic's hardcoded gap.

The resolved `trans_pos` / `focus_pos` (T1 grid voxels) are written into every
configured transducer, and a QC overlay (`plot_placement_t1_overlay`) is saved.
Downstream simulation is unchanged.

## Notes & limitations

- The standoff geometry assumes annular-style transducer fields
  (`curv_radius_mm`, `elem_od_mm`), the same assumption the heuristic makes. The
  mode is intended for annular (CTX-style) transducers.
- `trans_pos_mm` is interpreted as a *desired scalp entry*: it is snapped to the
  scalp, so the converted point need not lie exactly on the surface.
- Always visually verify the QC overlay — converting a scalp coordinate through
  any registration is approximate.

See also: [Transducer Placement overview](doc_placement.md),
[Heuristic placement](doc_placement_heuristic.md),
[Coordinate systems](doc_coordinate_systems.md).
