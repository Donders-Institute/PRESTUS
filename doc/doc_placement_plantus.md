# PlanTUS placement pipeline – step-by-step

> ⚠️ **Experimental.** PlanTUS placement is under active development. Results should be verified manually before use in a study.

This document describes the sequence of operations performed when
`placement.mode = 'plantus'` is set in a PRESTUS configuration.

## Overview

Given a T1-weighted image, a SimNIBS head mesh, a target-ROI mask, and a
transducer configuration YAML, PlanTUS ranks every skin surface vertex as a
candidate transducer position.

Two modes are available, controlled by `parameters.simulation.interactive`:

| Mode | Script | Vertex selection |
|------|--------|-----------------|
| **Headless** (`interactive = 0`, default) | `prestus_plantus_launcher.py` | Automatic — weighted composite score, no display required |
| **Interactive** (`interactive = 1`) | `prestus_plantus_gui_launcher.py` | Manual — five scoring maps displayed in Connectome Workbench; user clicks a vertex; the last selected vertex is placed when `wb_view` is closed |

Both modes run the same code for everything except vertex selection
(`plantus_common.py`): they compute the same scoring maps and composite score
(steps 2–5a below) and write the same output files, including the metric maps,
the Workbench scenes and `plantus_metrics.mat`.  Headless mode writes the
scenes but never opens `wb_view`.  Interactive mode additionally requires
`pynput` and Connectome Workbench `wb_view`; headless mode requires neither.

The original `PlanTUS_wrapper.py` (Mainz group, available separately at
`placement.plantus.script_path`) is *not* called directly because it has
hard-coded paths and unconditionally imports `pynput`.  The two PRESTUS scripts
above replicate the identical computation while accepting CLI arguments.

---

## Step-by-step pipeline

### 1  Parse inputs

| Input | Description |
|-------|-------------|
| `t1_filepath` | Subject T1.nii.gz (SimNIBS output) |
| `simnibs_mesh_filepath` | Head mesh `sub-XXX.msh` (CHARM output) |
| `target_mask_filepath` | Binary NIfTI sphere centred on the target |
| `tx_config_yaml` | Transducer + placement parameters (written by PRESTUS) |

The YAML is a flat key–value file.  
Key fields: `max_distance`, `min_distance`, `optimal_distance`,
`transducer_diameter`, `plane_offset`, `additional_offset`,
`focal_distance_list`, `flhm_list`, five `weight_*` scalars.

The distance fields are derived by `position_transducer_plantus.m` from the
transducer calibration (`focal_distance_list`, `flhm_list`), unless set
explicitly with `placement.plantus.min_distance_mm`, `max_distance_mm` and
`optimal_distance_mm`:

| Field | Value written by PRESTUS (when not set) |
|-------|--------------------------|
| `min_distance` | `max(0, min(focal_distance_list) − max(flhm_list))` |
| `max_distance` | `max(focal_distance_list) + max(flhm_list)` |
| `optimal_distance` | `mean(focal_distance_list)` — the desired skin–target distance |

If `max(flhm_list)` is NaN or ≤ 0, a 10 mm margin is used instead.  For a
steerable array, list the steering depths in `focal_distance_list` and set
the window to the steering range: the default margin (the largest FLHM,
often > 50 mm) otherwise admits skin spots the focus cannot reach.  If
`optimal_distance` is missing from the YAML, the launcher falls back to
`(min_distance + max_distance) / 2`.

PlanTUS derives the FLHM at the chosen focal distance from a cubic fit to
`(focal_distance_list, flhm_list)`, which needs at least 4 points. Both
launchers patch this (`plantus_compat.py`): non-finite pairs are dropped, and
with fewer than 4 points (e.g. the single O'Neill estimate used when no axial
profile CSV is given) FLHM is interpolated linearly and held constant outside
the calibrated range.

---

### 2  Convert SimNIBS mesh to surfaces

`PlanTUS.convert_simnibs_mesh_to_surface` extracts triangulated `.surf.gii`
files from the volumetric `.msh` mesh by tissue tag:

| Surface | SimNIBS tags | File |
|---------|-------------|------|
| Skin | 1005 | `skin.surf.gii` |
| Skull | 1007, 1008 | `skull.surf.gii` |

SimNIBS also produces `skull.stl` (used in the skull-entry-angle metric below).

---

### 3  Avoidance mask  (`create_avoidance_mask`)

Identifies skin vertices that **must not** be used as transducer positions:

1. **Air cavities / sinuses** – binarise `final_tissues.nii.gz`, fill holes
   with `wb_command -volume-fill-holes`, subtract to isolate air; tessellate
   to STL with FreeSurfer `mri_tessellate`; any skin vertex whose outward
   normal ray intersects an air cavity is masked out.
2. **Eyes** – extract eye mesh (tag 1006), find left/right centres via k-means
   on x-coordinates; mask out all skin vertices within 30 mm of each eye
   centre.
3. **Ears (tragus)** – read LPA/RPA from `eeg_positions/Fiducials.csv`; shift
   1.5 cm posterior; mask out all vertices within 15 mm.
4. **Erosion** – the binary avoidance metric is eroded by `transducer_radius`
   mm to leave a safety margin.

Vertices with `avoidance == 0` are excluded from the composite score.

> **External tools required**: `fslmaths`, `wb_command`, `mri_tessellate`,
> `mris_convert`.  If `final_tissues.nii.gz` is absent (e.g. CHARM not run)
> the avoidance step is skipped with a warning and all vertices remain
> eligible.

---

### 4  Per-vertex scoring metrics

For every skin vertex *i* five metrics are computed:

#### 4a  Skin–target distance  `d_i`

Euclidean distance from skin vertex *i* to the centre of mass of the target
mask (from `PlanTUS.roi_center_of_gravity`).

The score does **not** use `d_i` directly.  It uses the absolute deviation from
the desired focal distance, `|d_i − optimal_distance|`.  Smaller deviation is
better: the transducer's focal range should land on the target, not merely be
as close to it as possible.  Vertices outside
`[min_distance, max_distance]` are excluded (see 5a).

#### 4b  Skin–target angle  `α_i`

Angle (degrees) between the skin surface normal at *i* and the vector from *i*
to the target centre of mass.  
Lower is better (normal points towards target → good sonication path).

#### 4c  Skin–target intersection  `x_i`

The outward-pointing skin normal at *i* is cast as a ray of length 200 mm into
the volume.  The ray is intersected with the 3-D surface of the target mask STL
(generated by `PlanTUS.stl_from_nii`).  The metric is the total path length of
the ray *inside* the target (chord length for one entry–exit pair, sum of
chords for multiple pairs).  
Higher is better.

#### 4d  Skin–skull angle  `β_i`

The skin normal ray is intersected with the skull STL (ray length 40 mm).  The
closest skull surface vertex to the intersection point is found; the angle
between the skin normal and the skull normal at that vertex is computed.  
Lower is better (perpendicular skull entry → lower aberration).  
Where the ray meets no skull, `β_i` is NaN and scores as worst.

#### 4e  Skull thickness  `t_i`

Bone path length (mm) along the same inward skin-normal ray as 4d.  The CHARM
label volume `final_tissues.nii.gz` (next to the SimNIBS mesh) is sampled
every 0.25 mm (nearest voxel) up to 50 mm from the skin, and bone samples
(labels 4, 7, 8: bone, compact, spongy) are counted until the first
intracranial sample (labels 1–3: WM, GM, CSF).  
Lower is better (less attenuation and aberration).  
`t_i` is NaN, and scores as worst, where the ray meets no bone or does not
reach intracranial tissue within 50 mm (face, neck, rays grazing the skull
base).  The skull surface mesh is not used, because it does not include the
inner skull boundary.  Without `final_tissues.nii.gz`, `t_i` is NaN everywhere,
a warning is logged and `weight_skull_thickness` is ignored.

> **Known artefacts.**  A few vertices get implausibly thin values (≤ 1 mm):
> neck rays that pass ~45 mm of soft tissue into the spinal canal, and skin
> next to the ear where CHARM places bone or CSF directly under the skin.
> These vertices also have large skin–target angles and are rarely chosen
> unless thickness dominates the weights.

---

### 5a  Composite scoring (both modes; headless selects the optimum)

Each metric is normalised to **[0, 1]** (min–max over all finite values; NaN
or infinite values are set to 1 = worst, also for the inverted intersection
metric).  Because the range is set by the extremes, a few outlying vertices
(e.g. 40 mm skull-base rays for `t_i`) compress the spread of typical values and
make that criterion weaker than its nominal weight.  The intersection metric is
*inverted* (higher intersection → lower normalised cost).

Composite cost for vertex *i*:

```
score_i = ( w_dist  × norm(|d_i − optimal_distance|)
          + w_ang   × norm(α_i)
          + w_inter × (1 − norm(x_i))
          + w_skull × norm(β_i)
          + w_thick × norm(t_i) )
        ÷ (w_dist + w_ang + w_inter + w_skull + w_thick)
```

When `t_i` is unavailable (all NaN, or a `plantus_metrics.mat` written before
the criterion existed) `w_thick` is dropped from both numerator and
denominator.

Vertices with `avoidance == 0`, `d_i < min_distance` or `d_i > max_distance`
receive `score = inf`.

The vertex with the **minimum** composite score is selected as the transducer
position.  The interactive mode computes the same score; it reports the
optimum before `wb_view` opens, and the score and rank of each placed vertex.

**`plantus_metrics.mat`** (both modes; `placement.plantus.save_metrics`,
default `true`).  All per-vertex criteria are saved to
`<output_dir>/plantus_metrics.mat`: skin coordinates, normals and faces
(0-based), target centre, `d_i`, `α_i`, `x_i`, `β_i`, `t_i` (`skull_thickness`), the avoidance mask, the
composite score, `best_vertex` (composite optimum), `selected_vertex` (the
vertex actually placed), `selection_mode` (`auto` headless, `gui`
interactive, `manual` after `--do_only_trajectory`), `min_distance`,
`max_distance`, `optimal_distance` and the weights.  Vertex indices are
0-based.  This lets placements be re-scored under different weights without
repeating the ray casts (see `compute_composite_score` in `plantus_common.py`).
A `plantus_metrics.mat` from an earlier run in the same directory is deleted
at the start of each run (also when `save_metrics` is `false`), so it never
describes a different placement than the current one.

**Criteria evaluation** (both modes; `placement.plantus.evaluate`, default
`true`; needs `save_metrics`).  After saving the metrics, the launchers run
`plantus_evaluation.py` on `plantus_metrics.mat` (about 10 s).  It re-scores the
metrics with `compute_composite_score` under 12 weight sets (configured, equal,
one criterion at 0.55 with the rest sharing 0.45, and each criterion alone)
and writes, next to the metrics file:

| File | Contents |
|------|----------|
| `plantus_evaluation.json` | Placed vertex (rank, each criterion's value and % of eligible vertices that are better), per-criterion ranges over eligible vertices, the weight sweep and diagnostics; read by the [placement report](doc_placement.md#placement-report) |
| `plantus_weight_sensitivity.csv` | One row per weight set: chosen vertex, its criteria, distance to the placed vertex |
| `plantus_criteria_maps.png` | Each criterion and the composite score on the skin (4 views) |
| `plantus_weight_sensitivity.png` | Criteria at the placed spot and at each weight set's choice; all chosen spots on the skin |

Diagnostics flag settings under which a criterion cannot discriminate between
candidate positions:

| Code | Condition |
|------|-----------|
| `distance_range_unconstrained` | `[min_distance, max_distance]` excludes < 5% of the skin outside the avoidance mask (e.g. a large `flhm_list` value) |
| `intersection_sparse` | Beam–target intersection non-zero at < 2% of eligible vertices (small target mask) |
| `<criterion>_unavailable` / `_nan` | Criterion undefined at all / > 25% of eligible vertices |
| `placed_ineligible` | Placed vertex excluded under the configured scoring (interactive selection) |
| `placed_not_optimum` | Placed vertex is not the composite optimum (interactive selection) |
| `tradeoff_<criterion>` | At the placed vertex a weighted criterion is worse than at > 50% of eligible vertices |
| `weight_sensitive` | A single-criterion-heavy weight set moves the chosen spot > 20 mm |

Warnings are also printed to the placement log.  Evaluation failures are
logged and never abort a placement.  For an existing metrics file run
`python functions/neuronav/plantus_evaluation.py <plantus_metrics.mat>`.
The evaluation files from an earlier run are deleted with the metrics file.

### 5b  Interactive mode — GUI vertex selection

When `parameters.simulation.interactive = 1`, a planning scene is built from
the four normalised metric maps 4a–4d and loaded into Connectome Workbench (`wb_view`).
The skull-thickness map is written as `skull_thickness_skin.func.gii` but is not
part of the scene template; load it in `wb_view` by hand to inspect it.
The user inspects the maps on the skin surface and clicks a vertex; each click
replaces the current selection, and the vertex number is echoed to the MATLAB
command window as it is selected. **Closing `wb_view` confirms** the last
selected vertex, and its placement (step 6) is generated. Closing `wb_view`
without clicking a vertex aborts with an error, and no `*Localite.mat` is written.

The launcher's output is streamed to the MATLAB command window while `wb_view`
is open. MATLAB does not connect the launcher's stdin, so there is no console
yes/no prompt in this path.

When `prestus_plantus_gui_launcher.py` is run directly from a terminal, it
defaults to the original PlanTUS behaviour instead: after each click it asks
`Generate transducer placement for vertex N? (yes/no)`, and several vertices
can be confirmed in one session. The mode is set by `--selection`
(`auto` (default): `prompt` when stdin is a terminal, otherwise `last_click`).
If no vertex is confirmed, the launcher exits non-zero.

The interactive mode requires:
- A graphical display (desktop or VNC session; not a batch/HPC job)
- `pynput` Python package (mouse listener for reading wb_view output); on
  Linux it needs an X display, on macOS Accessibility permission for the
  terminal/MATLAB
- Connectome Workbench `wb_view` — either on `PATH` or specified via
  `placement.plantus.connectome_wb_path`

PRESTUS reads the most recently written `*Localite.mat`.

---

### 6  Placement geometry  (`prepare_acoustic_simulation`)

For the best vertex *v*:

1. Determine **beam direction**: if the skin normal at *v* intersects the
   target (`x_v > 0`), use the negated skin normal; otherwise use the negated
   skin-to-target vector.  Both are unit vectors pointing into the head.

2. Compute **transducer centre**:

   ```
   transducer_centre = skin_coordinate_v
                     − (plane_offset + additional_offset) × beam_direction
   ```

   `plane_offset` = distance from radiating surface to transducer exit plane.  
   `additional_offset` = gel pad / standoff thickness.

3. Build a **4 × 4 Localite position matrix** (RAS mm) with
   `PlanTUS.create_Localite_position_matrix`:
   - Column 3 (translation) = `transducer_centre`
   - Column 0 (x-axis) = `beam_direction` (unit vector)
   - Columns 1–2 = orthonormal frame computed via Gram–Schmidt

4. Save outputs to `<output_dir>/vtx<N>/`, where `<output_dir>` is the
   per-target work directory PRESTUS passes to both launchers,
   `<path.seg>/PlanTUS/sub-XXX/<target_name>/`:

   | File | Contents |
   |------|----------|
   | `*_Localite.mat` | `position_matrix` (4×4, RAS mm) — **used by PRESTUS** |
   | `*_Localite.txt` | Plain-text copy |
   | `*_kPlan.mat` | Same matrix in k-Plan convention (m, axes permuted) |
   | `*_Localite_XML.txt` | Fake Localite XML snippet |
   | `focus_position_matrix_*.txt` | Focus point transform |
   | `transducer_*.surf.gii` | Transducer surface model at position |
   | `scene.scene` | Workbench placement scene (both modes; opened only in interactive mode) |

   `*_Localite_XML.txt` is rewritten on each run for that vertex (PlanTUS
   itself appends to it).

   The per-target directory itself holds the skin/skull surfaces, the metric
   maps (`distances_skin`, `distances_skin_thresholded`, `angles_skin`,
   `target_intersection_skin`, `skin_skull_angles_skin`, `skull_thickness_skin`,
   `avoidance_skin`;
   `.func.gii`), the planning `scene.scene` and `plantus_metrics.mat`.

---

### 7  PRESTUS reads the output

`position_transducer_plantus.m` searches for the most recent `*Localite.mat`
under this target's work directory (`<path.seg>/PlanTUS/sub-XXX/<target_name>/`),
considering only files written during the current run: `*Localite.mat` files
that already existed before the launch (e.g. from an earlier run or another
vertex) are ignored unless they were modified since. If no new file is found it
errors (`PRESTUS:plantus:noOutput`). It then loads `position_matrix` and passes
it to `localite_matrix_to_positions` to convert the 4×4 RAS matrix into
`trans_pos` and `focus_pos` in voxel coordinates.

After the placement coordinates are resolved, a T1 overlay plot (`plot_placement_t1_overlay`) is generated automatically and written to the subject output folder. It shows the final transducer position and target overlaid on orthogonal T1 slices.

---

## Coordinate conventions

| Frame | Units | Description |
|-------|-------|-------------|
| Voxel / grid | integer indices | PRESTUS simulation grid |
| RAS+ (subject) | mm | T1 image world space; used by Localite matrix |
| MNI | mm | Template space; target specified here, mapped to RAS via SimNIBS |
| k-Plan | m, permuted axes | Not used by PRESTUS |

---

## External dependencies

| Tool | Headless | Interactive | Used in step |
|------|:--------:|:-----------:|-------------|
| SimNIBS Python (nibabel, scipy, nilearn) — see [installation guide](doc_installation.md#simnibs-installation) | ✓ | ✓ | All steps |
| PlanTUS Python library (`PlanTUS.py`) | ✓ | ✓ | All steps |
| PyYAML (or built-in fallback) | ✓ | ✓ | YAML config loading |
| FSL `fslmaths` | ✓ | ✓ | Avoidance mask – tissue binarisation |
| Connectome Workbench `wb_command` | ✓ | ✓ | Avoidance mask – hole filling |
| FreeSurfer `mri_tessellate` / `mris_convert` | ✓ | ✓ | Avoidance mask – air mesh |
| Connectome Workbench `wb_view` | — | ✓ | Interactive scene display |
| `pynput` Python package | — | ✓ | Mouse click listener for wb_view |

If the avoidance-mask external tools are unavailable, the launcher proceeds
without masking and prints a warning.

### Installing `pynput` for interactive mode

```bash
# In the SimNIBS conda environment:
conda activate simnibs_v4.6.0
pip install pynput
```
