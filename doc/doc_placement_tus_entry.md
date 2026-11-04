# TUS_entry Transducer Placement

!!! note "At a glance"
    - **What:** Placement with the external TUS_entry tool, which picks the entry point where the transducer sits most evenly on the scalp.
    - **When you need it:** When good mechanical contact of the transducer is the main concern.
    - **Before you start:** R with TUS_entry in a conda environment, and a SimNIBS segmentation (or another head mask).

`placement.mode: 'tus_entry'` runs the [TUS_entry](https://github.com/CyrilAtkinson/TUS_entry) R toolbox as an optional placement backend. It selects the scalp entry point with the smallest coupling geometry error: the standard deviation of the transducer-disk-edge to scalp distances along the trajectory, subject to an optional depth constraint. TUS_entry optimises **mechanical fit** (how perpendicular the transducer disk sits on the scalp); PlanTUS optimises **acoustic objectives** (beam–target overlap, skull thickness).

```yaml
placement:
  mode: tus_entry
  tus_entry:
    script_path: /path/to/TUS_entry/TUS_entry.R
    env_path:    /home/user/.conda/envs/tus_entry
    mni_target_mm: [-40, 20, 50]     # DLPFC example
    minimal_distance_cm: 4.24
    maximal_distance_cm: 7.46
    transducer_size_cm:  6.4
    target_name: DLPFC
```

??? note "Details: how TUS_entry selects the entry point"
    1. **Distance filter**: candidate entry points on a sphere of radius `maximal_distance` around the target centroid. Only trajectories that cross the scalp surface and avoid the exclusion zone are kept.
    2. **Spatial decimation**: candidates are thinned by `resolution`.
    3. **Angle scoring**: for each candidate, a circular transducer footprint (diameter `transducer_size`) perpendicular to the trajectory is scored by the SD of disk-edge-to-scalp distances. A lower SD is a better-coupled placement.
    4. **Final selection**: the angle-ranked candidates are fine-tuned in turn; the first that satisfies `minimal_distance` ≤ depth ≤ `maximal_distance` is returned.

## Installation

TUS_entry requires R and four CRAN packages. Use a dedicated conda environment:

```bash
conda env create -f examples/environment_tus_entry.yml
conda activate tus_entry
```

`environment_tus_entry.yml` is in `examples/` (run the command there or give its path). Then clone TUS_entry:

```bash
git clone https://github.com/CyrilAtkinson/TUS_entry /path/to/TUS_entry
```

Point `placement.tus_entry.script_path` to `TUS_entry.R` and `placement.tus_entry.env_path` to the conda env prefix (e.g. `~/.conda/envs/tus_entry`).

## Required configuration

| Field | Type | Description |
|-------|------|-------------|
| `placement.tus_entry.script_path` | string | Full path to `TUS_entry.R` |
| `placement.tus_entry.env_path` | string | Path to conda env containing `Rscript` |

TUS_entry needs a **solid** head mask. It walks each candidate ray from outside the head towards the target and takes the last voxel outside the mask as the scalp entry. With a hollow skin shell the head interior counts as outside, every candidate collapses onto the target, and TUS_entry fails (`replacement has 1 row, data has 0`). PRESTUS therefore builds the head mask itself: it reads `placement.tus_entry.scalp_mask` (default: the SimNIBS `m2m_sub-<ID>/final_tissues.nii.gz`), treats every nonzero voxel as head, fills it (3D and slice-wise, which also closes a skin shell open at the neck) and writes the result to the work directory. The source file is only read.

## Optional configuration

| Field | Default | Description |
|-------|---------|-------------|
| `placement.tus_entry.scalp_mask` | *(SimNIBS `final_tissues.nii.gz`)* | NIfTI head or scalp mask in the T1 grid (nonzero = head); filled into a solid head mask in the work directory |
| `placement.tus_entry.target_mask` | *(auto)* | Path to binary NIfTI target mask. When empty, a sphere is generated from `focus_pos` or `mni_target_mm`. |
| `placement.tus_entry.target_mask_radius_mm` | `5` | Radius [mm] for the auto-generated spherical target mask |
| `placement.tus_entry.mni_target_mm` | *(from focus_pos)* | [1×3] MNI target coordinate [mm]; converted to subject RAS via SimNIBS |
| `placement.tus_entry.exclu_mask` | `''` | Path to binary NIfTI exclusion zone (e.g. ventricles, vessels) |
| `placement.tus_entry.minimal_distance_cm` | `4.24` | Minimum scalp-to-target depth [cm] |
| `placement.tus_entry.maximal_distance_cm` | `7.46` | Maximum scalp-to-target depth [cm] |
| `placement.tus_entry.transducer_size_cm` | `6.4` | Transducer outer diameter [cm] |
| `placement.tus_entry.voxel_size_cm` | *(T1 voxel size)* | Voxel size [cm/voxel]; by default the mean T1 voxel size (the masks are on the T1 grid). TUS_entry assumes isotropic voxels |
| `placement.tus_entry.resolution_cm` | `0.33` | Spatial decimation step [cm] |
| `placement.tus_entry.return_candidates` | `false` | Save the full ranked candidate list as RDS (see below) |
| `placement.tus_entry.target_name` | `'target'` | Label for output filenames |
| `placement.tus_entry.work_dir` | *(auto)* | Override output directory (default: `<path.placement>/sub-<ID>/tus_entry/<target_name>/`) |

!!! note "Units"
    TUS_entry uses centimetres internally. PRESTUS converts distance outputs to millimetres, as for all other PRESTUS coordinates.

All keys are also listed in [Parameters](doc_parameters.md#placementtus_entry).

## Candidate list

By default, TUS_entry returns a single selected entry point. With `return_candidates = true`, the runner script saves the full angle-ranked candidate list to `tus_entry_candidates.rds` in the working directory. Inspect it in R:

```r
cands <- readRDS("tus_entry_candidates.rds")
head(cands)
```

MATLAB cannot read RDS files directly. To filter ranked candidates in PRESTUS (e.g. the entry with the lowest skull thickness among the top-10 coupling candidates), load the RDS in R, export to CSV and load that in MATLAB. The `candidates` output argument of `position_transducer_tus_entry` is reserved for a future native implementation and returns an empty struct array.

## Coordinate conventions

| Space | What TUS_entry provides | PRESTUS convention |
|-------|------------------------|-------------------|
| Voxel | 1-based array indices (same as MATLAB) | `trans_pos`, `focus_pos`: 1-based |
| World | NIfTI affine applied to voxel indices | `trans_pos_ras`, `focus_pos_ras`: RAS mm |

TUS_entry operates in the native voxel space of the input masks. **Use the same T1 as the reference for both SimNIBS segmentation and PRESTUS placement.**

`trans_pos` is derived from the `Scalp_coordinates` returned by TUS_entry by stepping outward along the trajectory by the transducer standoff (`curv_radius_mm − dist_geom_ep_mm`), the same offset applied in Localite and PlanTUS modes. `focus_pos` is derived by stepping inward from the scalp entry by `focal_distance_mm` (mean of the calibrated list if available, else the geometric focal distance).

## Output files

All outputs are written to the [placement run folder](doc_outputs.md#where-outputs-live), `<path.placement>/sub-<ID>/tus_entry/<target_name>/` (override with `placement.tus_entry.work_dir`):

| File | Content |
|------|---------|
| `run_tus_entry.R` | Auto-generated runner script (for reproducibility / debugging) |
| `tus_entry_report.csv` | Selected entry: voxel coords, distance [mm], angle score |
| `tus_entry_candidates.rds` | Full ranked candidate list (only if `return_candidates = true`) |
| `sub-<ID>_tus_entry_head_mask.nii.gz` | Solid head mask passed to TUS_entry as `scalp` |
| `placement_run.json` | Run record: resolved positions and the settings they depend on (reused per `placement.reuse`) |
| `sub-<ID>_<target_name>_trajectory.txt` | Brainsight trajectory export |
| `sub-<ID>_<target_name>_target_mask.nii.gz` | Auto-generated spherical target mask (if no target mask was provided) |

??? note "References"
    When using TUS_entry placement, cite: Atkinson C, et al. *TUS_entry: A tool for optimising transducer entry point selection in transcranial ultrasound stimulation.* GitHub: https://github.com/CyrilAtkinson/TUS_entry
