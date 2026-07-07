# Transducer Placement

Simulations in free water will place the transducer in a homogeneous water medium without the need to specify a target. For head smulations, any transducer needs to be positioned close to the scalp (possibly at some distance to allow for coupling), and a target needs to be specified.

PRESTUS expects coordinates of both the transducer (`transducer.trans_pos`) and target (`transducer.focus_pos`) to be specified as voxels in the T1w grid. The transducer coordinate describes the bowl of the transducer (i.e., not the exit plane in the case of curved transducers). Both coordinates need to be reported in the space of the planning image.

## Manual coordinate selection

Suitable coordinates (x/y/z) may be identified in preferred imaging software based on e.g. an fMRI hotspot, anatomical marker. 

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

The recommended preprocessing output is a small JSON file with `trans_pos_ras` and `focus_pos_ras` (RAS mm), generated once per session and fed back into PRESTUS via `ras_to_grid`. See [Neuronavigation read-in § Recommended simple output format](doc_placement_neuronav.md#recommended-simple-output-format).

An example workflow (`examples/demo_localite.m`) and a multi-format debugging script (`code/debug_localite_inputs.m`) are provided.

## Phantom / water simulations

For phantom and water simulations the transducer and focus positions can be specified in millimetres rather than voxel indices using two optional config fields:

| Field | Description |
|---|---|
| `trans_pos_mm` | Transducer position in mm from the grid origin (`[x, y, z]`). Converted to voxel indices at runtime using `grid.resolution_mm`. |
| `focus_pos_mm` | Focus position in mm from the grid origin (`[x, y, z]`). Same conversion applies. |

These fields are resolution-independent: the same config works unchanged when `grid.resolution_mm` is changed. When both voxel-index and mm fields are present, the mm fields take precedence. Phantom NIfTI outputs reuse the input header so that world-space orientation and voxel dimensions are consistent across all pipeline stages.

