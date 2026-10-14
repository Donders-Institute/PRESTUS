# Transducer Placement

Simulations in free water will place the transducer in a homogeneous water medium without the need to specify a target. For head simulations, any transducer needs to be positioned close to the scalp (possibly at some distance to allow for coupling), and a target needs to be specified.

Internally, PRESTUS works with the transducer (`transducer.trans_pos`) and target (`transducer.focus_pos`) as voxels of the T1w planning image. The transducer coordinate describes the bowl apex of the transducer (i.e., not the exit plane in the case of curved transducers). You can enter them directly (in voxels, subject RAS mm or MNI mm), or let a placement mode derive them (`placement.mode`).

## Choosing a placement mode

| Mode | Primary objective | Requires |
|------|------------------|----------|
| `manual` | User-specified | Coordinates in config (T1 voxels, RAS+ mm or MNI mm) |
| `localite` | Recorded neuronavigation pose | Localite XML |
| `mni` | Transducer and focus given in MNI mm; transducer snapped to the scalp | SimNIBS `m2m_sub-XXX` folder |
| `heuristic` | Skull-surface proximity to MNI target | SimNIBS segmentation |
| `plantus` | Beam–target overlap + skull acoustics | SimNIBS mesh + PlanTUS install |
| `tus_entry` | Transducer-scalp coupling geometry | Scalp mask + TUS_entry R install |


Each mode has its own page: [manual](#manual-coordinate-selection), [MNI](doc_placement_mni.md), [heuristic](doc_placement_heuristic.md), [PlanTUS](doc_placement_plantus.md), [TUS_entry](doc_placement_tus_entry.md) and [Localite](doc_placement_neuronav.md). Whatever the mode, the resolved positions are written to `sub-XXX_placement<affix>.json` in the output folder, and coordinate conventions are described in [Coordinate systems](doc_coordinate_systems.md).

## Placement report

After placement, `preproc_transducer_placement` writes a stand-alone HTML
report, `sub-XXX_placement_<mode><affix>.html` in `io.dir_reports`
(`placement.report`, default `true`). It holds the target, the transducer and
focus positions (T1 voxels and subject RAS mm, with transducer–focus and
focus–target distances), the T1 overlay figure (drawn for every mode) and the
placement configuration. For `placement.mode: 'plantus'` it adds the
[criteria evaluation](doc_placement_plantus.md#5a-composite-scoring-both-modes-headless-selects-the-optimum):
diagnostics, each criterion at the placed spot, the criterion maps and the
weight sweep. The report never stops the pipeline: failures are warnings, and
it is skipped when no output directory is known.

## Manual coordinate selection

Suitable coordinates (x/y/z) may be identified in preferred imaging software based on e.g. an fMRI hotspot, anatomical marker.

By default the coordinates are 1-based T1 voxel indices. Set `placement.manual.coordinate_space` to enter them in millimetres instead:

```yaml
placement:
  mode: 'manual'
  manual:
    coordinate_space: 'ras'   # 'voxel' (default) | 'ras' | 'mni'
transducer:
  trans_pos: [-62.1, 12.4, 48.0]   # bowl apex, subject RAS+ mm
  focus_pos: [-21.5, 10.2, 2.3]
```

- `ras` — subject RAS+ mm of the T1 (e.g. read from a viewer or neuronavigation), converted with `ras_to_grid()` via the T1 affine.
- `mni` — MNI mm. The transducer point is converted with the linear (12dof) transform and the focus with the nonlinear warp, as in [MNI placement](doc_placement_mni.md), but the transducer is **not** snapped to the scalp: the points are used exactly as given. Requires the SimNIBS `m2m_sub-XXX` folder.

Conversion happens once, per transducer, in the placement stage (so `modules.run_transducer_placement` must be `1`). Afterwards the positions are voxel indices, `coordinate_space` is set to `'voxel'`, and the original inputs are kept in `placement.manual.input_trans_pos` / `input_focus_pos`. If neither `focal_distance_ep` nor `focal_distance_bowl` is set, the distance fallback from the positions is computed after conversion.

## MNI coordinate selection

Both the transducer position and the focus can be specified directly in MNI (mm) coordinates via `placement.mode = 'mni'`. The focus is converted with the nonlinear SimNIBS warp (as the heuristic target), while the transducer point — which lies outside the brain — is converted with the linear (12dof) transform, snapped to the scalp, and offset by a configurable gap using the heuristic's standoff geometry. See [MNI Transducer Placement](doc_placement_mni.md).

## Heuristic coordinate selection

PRESTUS can identify heuristic locations for transducer placement (see [Heuristic Transducer Placement](doc_placement_heuristic.md)). This benefits iterative approaches without manual intervention.

PRESTUS also integrates [PlanTUS](doc_placement_plantus.md) (Lueckel et al., Mainz) as a native placement mode. The two automated modes differ in their inputs and optimisation strategy:

| | `heuristic` | `plantus` |
|---|---|---|
| **Input** | `final_tissues.nii.gz` (label volume) | Full SimNIBS mesh (`sub-XXX.msh`) + SimNIBS Python env |
| **Strategy** | Single-objective sphere-expansion on skull surface | Multi-objective optimisation (beam overlap, skin/skull angles, skull thickness, target distance) |
| **Focal distance sweep** | No | Yes |
| **External dependency** | None | PlanTUS + SimNIBS installation |

Use `heuristic` when only a tissue-label segmentation is available or when running batch jobs on an HPC cluster. Use `plantus` when the full SimNIBS mesh is available and multi-objective placement optimisation or transducer selection across focal distances is needed.

## Neuronavigation coordinate selection

PRESTUS provides helper functions to read coordinates acquired with neuronavigation systems. Currently, Localite is supported.

See [Neuronavigation read-in](doc_placement_neuronav.md) for full details. In brief:

1. Localite writes transducer position as a 4×4 transformation matrix (RAS mm) in one of several XML file formats (TriggerMarkers, GUMMarkers, InstrumentMarker).
2. `neuronav_compute_series_statistics` parses any of these formats and returns a mean 4×4 matrix per position series.
3. `localite_matrix_to_positions` converts the matrix to transducer bowl and focus positions in both RAS mm and T1 voxel space.

Ingestion (`neuronav_gui` or the script workflow) writes one mean pose per target to `<path.localite_post>/sub-XXX/ses-XX/InstrumentMarker_sub-XXX_ses-XX.xml`, with a QC sidecar JSON; `placement.mode: localite` reads these and converts them to grid voxels via `ras_to_grid`. See [Neuronavigation read-in § Output layout](doc_placement_neuronav.md#output-layout).

An example workflow is provided in `examples/demo_localite.m`.

## TUS_entry coordinate selection

`placement.mode = 'tus_entry'` runs the [TUS_entry](https://github.com/CyrilAtkinson/TUS_entry) R toolbox, which selects the scalp entry point with the most even transducer–scalp contact (mechanical fit), complementary to PlanTUS (acoustic objectives). See [TUS_entry placement](doc_placement_tus_entry.md).

---

## Phantom / water simulations

For phantom and water simulations the transducer and focus positions can be specified in millimetres rather than voxel indices using two optional config fields:

| Field | Description |
|---|---|
| `trans_pos_mm` | Transducer position in mm from the grid origin (`[x, y, z]`). Converted to voxel indices at runtime using `grid.resolution_mm`. |
| `focus_pos_mm` | Focus position in mm from the grid origin (`[x, y, z]`). Same conversion applies. |

These fields are resolution-independent: the same config works unchanged when `grid.resolution_mm` is changed. When both voxel-index and mm fields are present, the mm fields take precedence. Phantom NIfTI outputs reuse the input header so that world-space orientation and voxel dimensions are consistent across all pipeline stages.

