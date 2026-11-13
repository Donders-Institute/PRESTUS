# Transducer Placement

!!! note "At a glance"
    - **What:** How PRESTUS gets the transducer and focus positions of a head simulation: entered by hand, or derived by a placement mode (`placement.mode`).
    - **When you need it:** Every head simulation. Free-water simulations need no target.
    - **Before you start:** The T1w planning image; most modes also need a SimNIBS segmentation.

PRESTUS works with the transducer (`transducer.trans_pos`) and the target (`transducer.focus_pos`) as voxels of the T1w planning image. The transducer coordinate is the bowl apex, not the exit plane of a curved transducer. All modes put the exit plane `curv_radius_mm − dist_geom_ep_mm` in front of the apex ([BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain)'s `FocalLength − NaturalOutPlaneDistance`) and, where they place the transducer on the scalp, leave `placement.skin_gap_mm` (default 5 mm) between scalp and exit plane. You can enter both points directly (in voxels, subject RAS mm or MNI mm), or let a placement mode derive them.

## Choosing a placement mode

| Mode | Primary objective | Requires |
|------|------------------|----------|
| `manual` | User-specified | Coordinates in config (T1 voxels, RAS+ mm or MNI mm) |
| `localite` | Recorded neuronavigation pose | Localite XML |
| `mni` | Transducer and focus given in MNI mm; transducer snapped to the scalp | SimNIBS `m2m_sub-XXX` folder |
| `heuristic` | Skull-surface proximity to MNI target | SimNIBS segmentation |
| `plantus` | Beam–target overlap + skull acoustics | SimNIBS mesh + PlanTUS install |
| `tus_entry` | Transducer-scalp coupling geometry | Scalp mask + TUS_entry R install |

Each mode has its own page: [manual](#manual-coordinate-selection), [MNI](doc_placement_mni.md), [heuristic](doc_placement_heuristic.md), [PlanTUS](doc_placement_plantus.md), [TUS_entry](doc_placement_tus_entry.md) and [Localite](doc_placement_neuronav.md). Whatever the mode, the resolved positions are written to `sub-XXX_placement<affix>.json` in the output folder. Coordinate conventions are described in [Coordinate systems](doc_coordinate_systems.md).

??? note "Details: heuristic or PlanTUS?"
    | | `heuristic` | `plantus` |
    |---|---|---|
    | **Input** | `final_tissues.nii.gz` (label volume) | Full SimNIBS mesh (`sub-XXX.msh`) + SimNIBS Python env |
    | **Strategy** | Single-objective sphere-expansion on skull surface | Multi-objective optimisation (beam overlap, skin/skull angles, skull thickness, target distance) |
    | **Focal distance sweep** | No | Yes |
    | **External dependency** | None | PlanTUS + SimNIBS installation |

    Use `heuristic` when only a tissue-label segmentation is available or for batch jobs on an HPC cluster. Use `plantus` when the full SimNIBS mesh is available and you want multi-objective optimisation or transducer selection across focal distances.

## Placement report

After placement, PRESTUS writes a stand-alone HTML report, `sub-XXX_placement_<mode><affix>.html`, to `io.dir_reports` (`placement.report`, default `true`).

??? note "Details: report contents"
    It holds the target, the transducer and focus positions (T1 voxels and subject RAS mm, with transducer–focus and focus–target distances), the T1 overlay figure and the placement configuration. For `placement.mode: 'plantus'` it adds the [criteria evaluation](doc_placement_plantus.md#5a-composite-scoring-both-modes-headless-selects-the-optimum): diagnostics, each criterion at the placed spot, the criterion maps and the weight sweep. A failed report is a warning, not an error. The report is skipped when no output directory is known.

## Manual coordinate selection

Identify the coordinates in your imaging software, for example from an fMRI hotspot or an anatomical marker. By default they are 1-based T1 voxel indices. Set `placement.manual.coordinate_space` to enter millimetres instead:

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
- `mni`: MNI mm. The transducer point is converted with the linear (12dof) transform and the focus with the nonlinear warp, as in [MNI placement](doc_placement_mni.md), but the transducer is **not** snapped to the scalp. The points are used exactly as given. Requires the SimNIBS `m2m_sub-XXX` folder.

Conversion happens once per transducer in the placement stage, so `modules.run_transducer_placement` must be `1`.

??? note "Details: what happens after conversion"
    The positions become voxel indices and `coordinate_space` is set to `'voxel'`. The original inputs are kept in `placement.manual.input_trans_pos` / `input_focus_pos`. If neither `focal_distance_ep` nor `focal_distance_bowl` is set, the distance fallback from the positions is computed after conversion.
