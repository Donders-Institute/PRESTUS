# TUS_entry Transducer Placement

![Pipeline steps; this page belongs to 3 Placement](assets/figures/fig_stage-strip-placement_diagram.svg){ .stage-strip }

A simulated transducer position is useful only if the real transducer can be placed there with good contact on the scalp. `placement.mode: 'tus_entry'` runs the [TUS_entry](https://github.com/CyrilAtkinson/TUS_entry) [R](https://www.r-project.org/) toolbox, which selects the scalp entry point with the smallest coupling geometry error (the standard deviation of the transducer-disk-edge to scalp distances), subject to an optional depth constraint. TUS_entry optimises **mechanical fit** (how perpendicular the transducer disk sits on the scalp); [PlanTUS](https://github.com/mlueckel/PlanTUS) optimises **acoustic objectives** (beam–target overlap, skull thickness).

On this page:

- [Installation](#installation): R, its packages, and TUS_entry in a [conda](https://docs.conda.io/) environment.
- [Required Configuration](#required-configuration): the two paths PRESTUS needs and the head mask it builds.
- [Optional Configuration](#optional-configuration): masks, depth window, and other settings.
- [Candidate List](#candidate-list): the full ranked candidate list.
- [Coordinate Conventions](#coordinate-conventions): voxel and world frames, and how the exit plane and focus are set.
- [Output Files](#output-files): the files written per run.

!!! note "At a glance"
    - **When you need it:** When good mechanical contact of the transducer is the main concern.
    - **Before you start:** R with TUS_entry in a conda environment, and a [SimNIBS](https://simnibs.github.io/simnibs/) segmentation (or another head mask).

!!! warning "Experimental feature"
    TUS_entry placement is under active development and is not yet covered by the PRESTUS test suite. Check the [placement report](doc_placement.md#placement-report) before simulating.

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

??? note "How TUS_entry selects the entry point"
    1. **Distance filter**: candidate entry points on a sphere of radius `maximal_distance` around the target centroid. Only trajectories that cross the scalp surface and avoid the exclusion zone are kept.
    2. **Spatial decimation**: candidates are thinned by `resolution`.
    3. **Angle scoring**: for each candidate, a circular transducer footprint (diameter `transducer_size`) perpendicular to the trajectory is scored by the SD of disk-edge-to-scalp distances. A lower SD is a better-coupled placement.
    4. **Final selection**: the angle-ranked candidates are fine-tuned in turn; the first that satisfies `minimal_distance` ≤ depth ≤ `maximal_distance` is returned.

## Installation

TUS_entry requires R and four [CRAN](https://cran.r-project.org/) packages. Use a dedicated conda environment:

```bash
conda env create -f examples/environment_tus_entry.yml
conda activate tus_entry
```

`environment_tus_entry.yml` is in `examples/` (run the command there or give its path). Then clone TUS_entry:

```bash
git clone https://github.com/CyrilAtkinson/TUS_entry /path/to/TUS_entry
```

Point `placement.tus_entry.script_path` to `TUS_entry.R` and `placement.tus_entry.env_path` to the conda env prefix (e.g. `~/.conda/envs/tus_entry`).

## Required Configuration

Both `placement.tus_entry.script_path` (path to `TUS_entry.R`) and `placement.tus_entry.env_path` (conda environment with `Rscript`) are required.

TUS_entry needs a **solid** head mask. It walks each candidate ray from outside the head towards the target and takes the last voxel outside the mask as the scalp entry. With a hollow skin shell the head interior counts as outside, every candidate collapses onto the target, and TUS_entry fails (`replacement has 1 row, data has 0`). PRESTUS therefore builds the head mask itself: it reads `placement.tus_entry.scalp_mask` (default: the SimNIBS `m2m_sub-<ID>/final_tissues.nii.gz`), treats every nonzero voxel as head, fills it (3D and slice-wise, which also closes a skin shell open at the neck) and writes the result to the work directory. The source file is only read.

## Optional Configuration

The optional keys set the target and exclusion masks, the depth range, the transducer size, the sampling resolution, and the output folder. They are listed, with their defaults, under [TUS_entry placement](doc_parameters_placement.md#tus_entry-placement-placementtus_entry).

!!! note "Units"
    TUS_entry uses centimetres internally. PRESTUS converts distance outputs to millimetres, as for all other PRESTUS coordinates.

All keys are also listed in [Parameters](doc_parameters_placement.md#tus_entry-placement-placementtus_entry).

## Candidate List

By default, TUS_entry returns a single selected entry point. With `return_candidates = true`, the runner script saves the full angle-ranked candidate list to `tus_entry_candidates.rds` in the working directory. Inspect it in R:

```r
cands <- readRDS("tus_entry_candidates.rds")
head(cands)
```

[MATLAB](https://www.mathworks.com/products/matlab.html) cannot read RDS files directly. To filter ranked candidates in PRESTUS (e.g. the entry with the lowest skull thickness among the top-10 coupling candidates), load the RDS in R, export to CSV, and load that in MATLAB. The `candidates` output argument of `position_transducer_tus_entry` is reserved for a future native implementation and returns an empty struct array.

## Coordinate Conventions

| Space | What TUS_entry provides | PRESTUS convention |
|-------|------------------------|-------------------|
| Voxel | 1-based array indices (same as MATLAB) | `trans_pos`, `focus_pos`: 1-based |
| World | NIfTI affine applied to voxel indices | `trans_pos_ras`, `focus_pos_ras`: RAS mm |

TUS_entry operates in the native voxel space of the input masks. **Use the same T1 as the reference for both SimNIBS segmentation and PRESTUS placement.**

The [exit plane](doc_glossary.md#exit-plane) lies `placement.skin_gap_mm` (default 5 mm) outward from the `Scalp_coordinates` returned by TUS_entry, along the trajectory. `trans_pos` lies a further `curv_radius_mm − dist_geom_ep_mm` outward, the offset used by all placement modes. `focus_pos` is derived by stepping inward from the exit plane by `transducer.focal_distance_ep` (or `focal_distance_bowl` minus the apex-to-exit-plane offset, as `focal_distance_calculation` converts it; if neither is set, the geometric focus `dist_geom_ep_mm`).

## Output Files

All outputs are written to the [placement run folder](doc_outputs.md#where-outputs-live), `<path.placement>/sub-<ID>/tus_entry/<target_name>/` (override with `placement.tus_entry.work_dir`):

| File | Content |
|------|---------|
| `run_tus_entry.R` | Auto-generated runner script (for reproducibility / debugging) |
| `tus_entry_report.csv` | Selected entry: voxel coords, distance [mm], angle score |
| `tus_entry_candidates.rds` | Full ranked candidate list (only if `return_candidates = true`) |
| `sub-<ID>_tus_entry_head_mask.nii.gz` | Solid head mask passed to TUS_entry as `scalp` |
| `placement_run.json` | Run record: resolved positions and the settings they depend on (reused per `placement.reuse`) |
| `sub-<ID>_<target_name>_trajectory.txt` | [Brainsight](https://www.rogue-research.com/tms/brainsight-tms/) trajectory export |
| `sub-<ID>_<target_name>_target_mask.nii.gz` | Auto-generated spherical target mask (if no target mask was provided) |

??? note "References"
    When using TUS_entry placement, cite: Atkinson C, et al. *TUS_entry: A tool for optimising transducer entry point selection in transcranial ultrasound stimulation.* GitHub: https://github.com/CyrilAtkinson/TUS_entry
