# MNI Transducer Placement

![Pipeline steps; this page belongs to 3 Placement](assets/figures/fig_stage-strip-placement_diagram.svg){ .stage-strip }

A transducer site and target given in MNI coordinates let one plan apply to every subject, but each subject's simulation needs both points in its own anatomy. `placement.mode: 'mni'` converts both points to the subject's T1 grid, snaps the transducer to the scalp, and sets its [exit plane](doc_glossary.md#exit-plane) `placement.skin_gap_mm` (default 5 mm) off the scalp; no search is performed, so the result is deterministic. The [heuristic](doc_placement_heuristic.md) mode, by contrast, takes only an MNI target and searches the scalp for a transducer position.

On this page:

- [How a Position Is Resolved](#how-a-position-is-resolved): coordinate conversion, scalp snap, and standoff.
- [Limitations](#limitations): the transducer types the mode assumes and what to check.

!!! note "At a glance"
    - **When you need it:** When you know both points in template space, e.g. a group target and a standard scalp position.
    - **Before you start:** A [SimNIBS](https://simnibs.github.io/simnibs/) segmentation of the subject (for the MNI-to-subject transforms).

```yaml
placement:
  mode: 'mni'
  mni:
    trans_pos_mm: [-52,  -8,  58]   # desired transducer (scalp entry) in MNI mm
    focus_pos_mm: [-29,  -6,   7]   # focus / target in MNI mm
  skin_gap_mm: 5                    # gap between scalp and exit plane [mm]
```

A SimNIBS `m2m_sub-XXX` folder (under `path.seg`) is required: the conversion relies on `final_tissues.nii.gz` and the SimNIBS coordinate transforms.

## How a Position Is Resolved

1. **Coordinate conversion.** Both points are mapped MNI → subject RAS+ → grid voxel with charm's nonlinear warp (`nonl`), the only MNI transform SimNIBS 4 writes; the focus as the heuristic target.

2. **Snap to scalp.** The warp is less exact near the scalp, so the converted transducer voxel is snapped to the nearest scalp voxel; a slightly off MNI entry coordinate still lands on the head surface.

3. **Standoff.** The transducer is offset outward from the scalp along the focal axis (target → scalp) so that its exit plane lies `placement.skin_gap_mm` (default 5 mm) off the scalp, as in the heuristic. The bowl apex sits a further `curv_radius_mm − dist_geom_ep_mm` behind the exit plane.

The resolved `trans_pos` / `focus_pos` (T1 grid voxels) are written into every configured transducer, and a QC overlay is saved. The log prints the resulting bowl-to-focus distance, and warns if it differs from a set `transducer.focal_distance_bowl` by more than 5 mm: the requested focus may then lie outside what the transducer can reach. Downstream simulation is unchanged.

??? note "Functions involved"
    - Conversion: `transform_coordinates` and `mni2subject_coords_LDfix`, shared with the heuristic.
    - Scalp snap: `tp_scalp_boundary` (nearest scalp / outer-boundary voxel).
    - Standoff: `transducer_scalp_geometry`, shared with `tp_candidate_mesh`; the geometry the heuristic applies to every candidate.
    - QC overlay: `plot_placement_t1_overlay`.

## Limitations

- The standoff geometry assumes annular-style transducer fields (`curv_radius_mm`, `elem_od_mm`), as the heuristic does. The mode is intended for annular (CTX-style) transducers.
- `trans_pos_mm` is a desired scalp entry. It is snapped to the scalp, so the converted point need not lie exactly on the surface.
- Check the QC overlay. Converting a scalp coordinate through any registration is approximate.

See also: [Transducer Placement overview](doc_placement.md), [Heuristic placement](doc_placement_heuristic.md), [Coordinate systems](doc_coordinate_systems.md).
