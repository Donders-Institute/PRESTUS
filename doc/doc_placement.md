# Transducer Placement

Simulations in free water will place the transducer in a homogeneous water medium without the need to specify a target. For head simulations, any transducer needs to be positioned close to the scalp (possibly at some distance to allow for coupling), and a target needs to be specified.

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

## TUS_entry coordinate selection

`placement.mode = 'tus_entry'` integrates the [TUS_entry](https://github.com/CyrilAtkinson/TUS_entry) R toolbox as an optional PRESTUS transducer placement backend. It selects a scalp entry point that minimises coupling geometry error — specifically, the standard deviation of transducer-disk-edge to scalp distances along the chosen trajectory — subject to an optional depth constraint.

This is complementary to PlanTUS: TUS_entry optimises **mechanical fit** (how perpendicular the transducer disk sits on the scalp), whereas PlanTUS optimises **acoustic objectives** (beam–target overlap, skull thickness, etc.).

### Algorithm

TUS_entry operates in four stages:

1. **Distance filter** — generates candidate entry points on a sphere at radius `maximal_distance` around the target centroid. Retains only trajectories that cross the scalp surface and avoid any exclusion zone.
2. **Spatial decimation** — thins candidates by `resolution` to reduce redundant computation.
3. **Angle scoring** — for each candidate, constructs a circular transducer footprint (diameter = `transducer_size`) perpendicular to the trajectory and scores it by the SD of disk-edge-to-scalp distances. Lower SD = better-coupled placement.
4. **Final selection** — iterates through angle-ranked candidates, fine-tunes the entry point position, and returns the first that satisfies the depth constraint (`minimal_distance` ≤ depth ≤ `maximal_distance`).

### Installation

TUS_entry requires R and four CRAN packages. The recommended approach is a dedicated conda environment:

```bash
conda env create -f environment_tus_entry.yml
conda activate tus_entry
```

The `environment_tus_entry.yml` is provided in the PRESTUS root. Then clone TUS_entry:

```bash
git clone https://github.com/CyrilAtkinson/TUS_entry /path/to/TUS_entry
```

Point `placement.tus_entry.script_path` to `TUS_entry.R` and `placement.tus_entry.env_path` to the conda env prefix (e.g. `~/.conda/envs/tus_entry`).

### Required configuration

| Field | Type | Description |
|-------|------|-------------|
| `placement.tus_entry.script_path` | string | Full path to `TUS_entry.R` |
| `placement.tus_entry.scalp_mask` | string | Path to binary NIfTI scalp mask (same space as T1) |
| `placement.tus_entry.env_path` | string | Path to conda env containing `Rscript` |

The scalp mask is typically the `skin.nii.gz` output from SimNIBS CHARM, which PRESTUS already generates during segmentation. No additional preprocessing is needed.

### Optional configuration

| Field | Default | Description |
|-------|---------|-------------|
| `placement.tus_entry.target_mask` | *(auto)* | Path to binary NIfTI target mask. When empty, a sphere is generated from `focus_pos` or `mni_target_mm`. |
| `placement.tus_entry.target_mask_radius_mm` | `5` | Radius [mm] for auto-generated spherical target mask |
| `placement.tus_entry.mni_target_mm` | *(from focus_pos)* | [1×3] MNI target coordinate [mm]; converted to subject RAS via SimNIBS |
| `placement.tus_entry.exclu_mask` | `''` | Path to binary NIfTI exclusion zone (e.g. ventricles, vessels) |
| `placement.tus_entry.minimal_distance_cm` | `4.24` | Minimum scalp-to-target depth [cm] |
| `placement.tus_entry.maximal_distance_cm` | `7.46` | Maximum scalp-to-target depth [cm] |
| `placement.tus_entry.transducer_size_cm` | `6.4` | Transducer outer diameter [cm] |
| `placement.tus_entry.voxel_size_cm` | `0.1` | Voxel size [cm/voxel] — must match input mask resolution |
| `placement.tus_entry.resolution_cm` | `0.33` | Spatial decimation step [cm] |
| `placement.tus_entry.return_candidates` | `false` | Save full ranked candidate list as RDS (see below) |
| `placement.tus_entry.target_name` | `'target'` | Label for output filenames |
| `placement.tus_entry.work_dir` | *(auto)* | Override output directory (default: `path.seg/TUS_entry/sub-<ID>/<target_name>/`) |

> **Note on units:** TUS_entry uses centimetres internally. The PRESTUS bridge converts distance outputs to millimetres for consistency with all other PRESTUS coordinates.

### Example configuration

```yaml
placement:
  mode: tus_entry
  tus_entry:
    script_path: /path/to/TUS_entry/TUS_entry.R
    scalp_mask:  /path/to/simnibs/m2m_sub-001/skin.nii.gz
    env_path:    /home/user/.conda/envs/tus_entry
    mni_target_mm: [-40, 20, 50]     # DLPFC example
    minimal_distance_cm: 4.24
    maximal_distance_cm: 7.46
    transducer_size_cm:  6.4
    target_name: DLPFC
```

### Candidate list

By default, TUS_entry returns a single selected entry point. When `return_candidates = true`, the R runner script saves the full angle-ranked candidate list to `tus_entry_candidates.rds` in the working directory. This file can be inspected in R:

```r
cands <- readRDS("tus_entry_candidates.rds")
head(cands)
```

MATLAB cannot read RDS files directly. If you need ranked candidates for downstream PRESTUS filtering (e.g. to pick the entry with the lowest skull thickness among the top-10 coupling candidates), load the RDS in a separate R step and export to CSV, then load in MATLAB manually.

> **Note:** The `candidates` output argument of `position_transducer_tus_entry` is reserved for a future native implementation. It currently returns an empty struct array.

### Coordinate conventions

| Space | What TUS_entry provides | PRESTUS convention |
|-------|------------------------|-------------------|
| Voxel | 1-based array indices (same as MATLAB) | `trans_pos`, `focus_pos` — 1-based |
| World | NIfTI affine applied to voxel indices | `trans_pos_ras`, `focus_pos_ras` — RAS mm |

TUS_entry operates in the native voxel space of the input masks. **Use the same T1 as the reference for both SimNIBS segmentation and PRESTUS placement** to ensure coordinate consistency.

`trans_pos` is derived from the `Scalp_coordinates` returned by TUS_entry by stepping outward along the trajectory direction by the transducer standoff (`curv_radius_mm − dist_geom_ep_mm`), the same offset applied in Localite and PlanTUS modes.

`focus_pos` is derived by stepping inward from the scalp entry by `focal_distance_mm` (mean of calibrated list if available, else geometric focal distance).

### Output files

All outputs are written to `path.seg/TUS_entry/sub-<ID>/<target_name>/`:

| File | Content |
|------|---------|
| `run_tus_entry.R` | Auto-generated runner script (for reproducibility / debugging) |
| `tus_entry_report.csv` | Selected entry: voxel coords, distance [mm], angle score |
| `tus_entry_candidates.rds` | Full ranked candidate list (only if `return_candidates = true`) |
| `sub-<ID>_<target_name>_target_mask.nii.gz` | Auto-generated spherical target mask (if no target mask was provided) |

### Comparison with other placement modes

| Mode | Primary objective | Requires |
|------|------------------|----------|
| `manual` | User-specified | Coordinates in config |
| `localite` | Recorded neuronavigation pose | Localite XML |
| `heuristic` | Skull-surface proximity to MNI target | SimNIBS segmentation |
| `plantus` | Beam–target overlap + skull acoustics | SimNIBS mesh + PlanTUS install |
| `tus_entry` | Transducer-scalp coupling geometry | Scalp mask + TUS_entry R install |

When using TUS_entry placement, cite:

> Atkinson C, et al. *TUS_entry: A tool for optimising transducer entry point selection in transcranial ultrasound stimulation.* GitHub: https://github.com/CyrilAtkinson/TUS_entry

---

## Phantom / water simulations

For phantom and water simulations the transducer and focus positions can be specified in millimetres rather than voxel indices using two optional config fields:

| Field | Description |
|---|---|
| `trans_pos_mm` | Transducer position in mm from the grid origin (`[x, y, z]`). Converted to voxel indices at runtime using `grid.resolution_mm`. |
| `focus_pos_mm` | Focus position in mm from the grid origin (`[x, y, z]`). Same conversion applies. |

These fields are resolution-independent: the same config works unchanged when `grid.resolution_mm` is changed. When both voxel-index and mm fields are present, the mm fields take precedence. Phantom NIfTI outputs reuse the input header so that world-space orientation and voxel dimensions are consistent across all pipeline stages.

