# Transducer Placement

![Pipeline steps; this page belongs to 3 Placement](assets/figures/fig_stage-strip-placement_diagram.svg){ .stage-strip }

The transducer position decides where the beam enters the head, which part of the skull it crosses, and where the focus lies relative to the target, so it shapes every acoustic and thermal result of a head simulation. This page covers how PRESTUS obtains the transducer and focus positions: entered by hand, or derived by one of several placement modes that share one exit-plane geometry.

On this page:

- [Choosing a Placement Mode](#choosing-a-placement-mode): the modes, what each optimises, and what each needs.
- [Importing a Trajectory](#importing-a-trajectory): positions from Brainsight, 3D Slicer, or k-Plan.
- [Exporting a Pose](#exporting-a-pose): the placed pose for Brainsight, 3D Slicer, or k-Plan.
- [Placement Report](#placement-report): the HTML report written after placement.

!!! note "At a glance"
    - **When you need it:** Every head simulation. Free-water simulations need no target.
    - **Before you start:** The T1w planning image; most modes also need a [SimNIBS](https://simnibs.github.io/simnibs/) segmentation.

PRESTUS works with the transducer (`transducer.trans_pos`) and the target (`transducer.focus_pos`) as voxels of the T1w planning image. The transducer coordinate is the bowl apex, not the [exit plane](doc_glossary.md#exit-plane) of a curved transducer. All modes put the exit plane `curv_radius_mm − dist_geom_ep_mm` in front of the apex ([BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain)'s `FocalLength − NaturalOutPlaneDistance`) and, where they place the transducer on the scalp, leave `placement.skin_gap_mm` (default 5 mm) between scalp and exit plane. You can enter both points directly (in voxels, subject RAS mm, or MNI mm), or let a placement mode derive them.

## Choosing a Placement Mode

| Mode | Primary objective | Requires |
|------|------------------|----------|
| ![](assets/icons/icon_place-manual_glyph.svg){ .choice-icon } `manual` | User-specified | Coordinates in config (T1 voxels, RAS+ mm or MNI mm) |
| ![](assets/icons/icon_place-localite_glyph.svg){ .choice-icon } `localite` | Recorded neuronavigation pose | [Localite](https://www.localite.de/en/) XML |
| ![](assets/icons/icon_place-mni_glyph.svg){ .choice-icon } `mni` | Transducer and focus given in MNI mm; transducer snapped to the scalp | SimNIBS `m2m_sub-XXX` folder |
| ![](assets/icons/icon_place-heuristic_glyph.svg){ .choice-icon } `heuristic` | Skull-surface proximity to MNI target | SimNIBS segmentation |
| ![](assets/icons/icon_place-plantus_glyph.svg){ .choice-icon } `plantus` | Beam–target overlap + skull acoustics | SimNIBS mesh + [PlanTUS](https://github.com/mlueckel/PlanTUS) install |
| ![](assets/icons/icon_place-tus-entry_glyph.svg){ .choice-icon } `tus_entry` | Transducer-scalp coupling geometry | Scalp mask + [TUS_entry](https://github.com/CyrilAtkinson/TUS_entry) [R](https://www.r-project.org/) install |

Each mode has its own page: [manual](doc_placement_manual.md), [MNI](doc_placement_mni.md), [heuristic](doc_placement_heuristic.md), [PlanTUS](doc_placement_plantus.md), [TUS_entry](doc_placement_tus_entry.md) and [Localite](doc_placement_neuronav.md). Whatever the mode, the resolved positions are written to `sub-XXX_placement<affix>.json` in the output folder. Coordinate conventions are described in [Coordinate systems](doc_coordinate_systems.md).

??? note "Heuristic or PlanTUS?"
    | | `heuristic` | `plantus` |
    |---|---|---|
    | **Input** | `final_tissues.nii.gz` (label volume) | Full SimNIBS mesh (`sub-XXX.msh`) + SimNIBS [Python](https://www.python.org/) env |
    | **Strategy** | Single-objective sphere-expansion on skull surface | Multi-objective optimisation (beam overlap, skin/skull angles, skull thickness, target distance) |
    | **Focal distance sweep** | No | Yes |
    | **External dependency** | None | PlanTUS + SimNIBS installation |

    Use `heuristic` when only a tissue-label segmentation is available or for batch jobs on an HPC cluster. Use `plantus` when the full SimNIBS mesh is available and you want multi-objective optimisation or transducer selection across focal distances.

## Importing a Trajectory

To reuse a plan from another tool, use *Import trajectory* on the Nexus [Placement tab](doc_nexus_tabs.md). It writes `focus_pos` and `trans_pos` for `manual` placement, except for a Localite export: that sets the [`localite` placement mode](doc_placement_neuronav.md) with `placement.localite.file`, so each run averages the chosen series' triggers and converts the pose again. Reading a Localite export needs MATLAB.

| Source | File | Sets |
|---|---|---|
| [Brainsight](https://www.rogue-research.com/tms/brainsight-tms/) | target export (`.txt`, NIfTI coordinates) | focus and beam direction; transducer on the axis at the focal depth setting |
| [3D Slicer](https://www.slicer.org/) | ITK transform (`.tfm`, `.txt`) | as Brainsight |
| [k-Plan](https://k-plan.io) | `.kdp` (device parameters), `PLAN-*-RESULTS.h5` | transducer at the plan's bowl apex, focus at its target, focal depth setting matched to the target |
| k-Plan | `.kps` (transducer positions) | transducer at the bowl apex; focus on the axis at the focal depth setting |
| [Localite](https://www.localite.de/en/) | TriggerMarkers, InstrumentMarker, or GUMMarkers export (`.xml`) | the series you choose: `placement.localite.target` if its description is unique in the file, else `.position` |

??? note "k-Plan trajectories"
    - **Image:** k-Plan counts its coordinates from the first voxel of its planning image. The import is right when that image had the T1's voxels: the T1 itself, or, for an oblique T1, a copy with the rotation removed from the header. With a rotated or flipped T1, the import says so.
    - **Focal depth:** set from the selected transducer's own exit plane, so it can differ from k-Plan's value for the same apex and target. k-Plan's CTX-500, for example, puts its exit plane 10.82 mm in front of the apex, PRESTUS's 10.56 mm.
    - **Not imported:** k-Plan's transducer model, drive, and fields.

## Exporting a Pose

To plan or navigate with the PRESTUS pose elsewhere, list the formats in `placement.export_formats`. After placement, PRESTUS writes them to `io.dir_output`, one set per transducer (`_txN` with several):

| Format | Output | Contents |
|---|---|---|
| `'brainsight'` | `sub-XXX_placement_brainsight.txt` | target and beam direction (NIfTI scanner coordinates) |
| `'slicer'` | `sub-XXX_placement_slicer.tfm` | the same pose as an ITK transform, in BabelBrain's convention |
| `'kplan'` | `sub-XXX_placement_kplan/` | T1 and CT or pseudo-CT with the rotation removed from the header, the transducer position (`.kps`), and an entry sheet with the target and focal distance |

Each export imports back into the Nexus with *Import trajectory*. For k-Plan, upload the exported images, not the originals: k-Plan counts its coordinates from their first voxel. The T1 must have (obliquely) R, A, S voxel axes. The entry sheet gives k-Plan's focal distance for the CTX-500 (exit plane 10.82 mm in front of the apex), and for other devices how to work it out.

## Placement Report

After placement, PRESTUS writes a stand-alone HTML report, `sub-XXX_placement_<mode><affix>.html`, to `io.dir_reports` (`placement.report`, default `true`).

??? note "Report contents"
    It holds the target, the transducer and focus positions (T1 voxels and subject RAS mm, with transducer–focus and focus–target distances), the T1 overlay figure, and the placement configuration. For `placement.mode: 'plantus'` it adds the [criteria evaluation](doc_placement_plantus_pipeline.md#5a-combine-and-evaluate-scores): diagnostics, each criterion at the placed spot, the criterion maps, and the weight sweep. A failed report is a warning, not an error. The report is skipped when no output directory is known.
