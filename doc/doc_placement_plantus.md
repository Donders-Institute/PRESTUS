# PlanTUS Placement

The entry point decides which part of the skull the beam crosses and whether the focus can reach the target. Given a T1-weighted image, a [SimNIBS](https://simnibs.github.io/simnibs/) head mesh, a target-ROI mask and a transducer configuration, [PlanTUS](https://github.com/mlueckel/PlanTUS) (`placement.mode: 'plantus'`) scores every skin vertex as a candidate transducer position on skin–target distance, skin and skull angles, skull thickness and beam–target overlap. PRESTUS places the transducer at the vertex with the best weighted score, or at a vertex you choose on the maps on the Nexus PlanTUS page (or in [Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench)).

On this page:

- [Running PlanTUS](#running-plantus): the steps of a placement run.
- [Choosing the Position](#choosing-the-position): automatic, on the Nexus PlanTUS page, or in Workbench.
- [Step-by-Step Pipeline](#step-by-step-pipeline): inputs, criteria, composite score, evaluation and placement geometry.
- [Coordinate Conventions](#coordinate-conventions): the frames PlanTUS and PRESTUS exchange.
- [External Dependencies](#external-dependencies): the tools each mode needs.

!!! note "At a glance"
    - **When you need it:** When you want an optimised, documented choice of entry point, or to compare focal depths.
    - **Before you start:** The full SimNIBS mesh and a PlanTUS installation, 2.0 or later ([Installation](doc_installation.md#plantus-dependencies)). PlanTUS 1.x also works but needs Connectome Workbench, FSL and FreeSurfer.

!!! warning "Experimental"
    PlanTUS placement is under active development. Verify results manually before use in a study.

!!! warning "Limitations"
    - The Nexus PlanTUS page needs the maps of an earlier run (`save_metrics: true`) and draws them where the Nexus runs; the run itself can be headless.
    - Choosing in Workbench needs a display (desktop or VNC), Connectome Workbench `wb_view` and `pynput`. On macOS, the terminal or [MATLAB](https://www.mathworks.com/products/matlab.html) needs Accessibility permission for `pynput`. Without a display (Linux), PlanTUS runs headless.
    - The Workbench review dialog needs [PyQt5](https://www.riverbankcomputing.com/software/pyqt/) in the SimNIBS [Python](https://www.python.org/). It ships with SimNIBS; without it, the last clicked vertex is placed without review.
    - With PlanTUS 1.x, the avoidance mask needs [FSL](https://fsl.fmrib.ox.ac.uk/fsl/), [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/) and Workbench tools. Without them, all skin vertices are eligible (step 3).
    - To see errors from the launcher directly, run `prestus_plantus_gui_launcher.py` from a terminal (step 5b).

## Running PlanTUS

1. Set `placement.mode: plantus`, the `placement.plantus` paths and target, and the calibrated `focal_distance_list` / `flhm_list`.
2. In the Nexus, click *Compute maps* on the Placement card (or on Placement ▸ PlanTUS). PlanTUS computes the criterion maps (a few minutes) and places the transducer at the composite-score optimum; no simulation runs. In MATLAB, run the pipeline with `placement.plantus.vertex: []`.
3. Open Placement ▸ PlanTUS. The skin is coloured by the composite score; choose a criterion under *Map* to see one map. Click a skin position: the panel lists its criteria, composite score and rank, and the view shows the beam axis and the transducer there, in the housing Placement ▸ 3D draws.
4. Click *Place here*. The Nexus sets `placement.plantus.vertex` and places the transducer at that vertex on the same maps, without recomputing them. The config keeps the choice, so a later simulation run uses it.
5. Check the [placement report](doc_placement.md#placement-report) and the Viewer. Section [5a](#5a-combine-and-evaluate-scores) explains the evaluation it contains.

## Choosing the Position

`placement.plantus.vertex` and `placement.plantus.interactive` (`[]`, the default, follows `simulation.interactive`) select how the vertex is chosen:

| Setting | Script | Vertex selection |
|------|--------|-----------------|
| **Automatic** (`vertex: []`, `interactive: 0`) | `prestus_plantus_launcher.py` | Weighted composite score; no display required |
| **Chosen vertex** (`vertex: N`) | `prestus_plantus_launcher.py --do_only_trajectory N` | Vertex `N` (0-based) of the maps of an earlier run in the same placement run folder, chosen on the Nexus PlanTUS page; the maps are not recomputed |
| **Workbench** (`vertex: []`, `interactive: 1`) | `prestus_plantus_gui_launcher.py` | The five criterion maps in Connectome Workbench; you click a vertex, close `wb_view`, see the transducer there and confirm it in a review dialog |

Automatic and Workbench runs compute the same maps and composite score (steps 2 to 5a) and write the same files. A chosen vertex is placed only while the target is where the maps were computed for: if the target moved more than 1 mm, the run stops and asks you to clear `placement.plantus.vertex` and compute the maps again. A chosen vertex belongs to one subject's maps; for a study, give it per run in a manifest column. Only Workbench runs need `pynput`, `wb_view` and PyQt5.

Both modes share all code except vertex selection (`plantus_common.py`). The scripts replicate `PlanTUS_wrapper.py` (Mainz group, at `placement.plantus.script_path`) but take CLI arguments instead of hard-coded paths.

## Step-by-Step Pipeline

### 1  Parse Inputs

PRESTUS writes a flat transducer YAML for PlanTUS.

??? note "The skin–target distance window"
    `position_transducer_plantus.m` derives the skin–target distance window from the transducer calibration, unless `placement.plantus.min_distance_mm`, `max_distance_mm` and `optimal_distance_mm` are set:

    | Field | Value written by PRESTUS (when not set) |
    |-------|--------------------------|
    | `min_distance` | `max(0, min(focal_distance_list) − max(flhm_list))` |
    | `max_distance` | `max(focal_distance_list) + max(flhm_list)` |
    | `optimal_distance` | `mean(focal_distance_list)`, the desired skin–target distance |

    If `max(flhm_list)` is NaN or ≤ 0, a 10 mm margin is used instead. For a steerable array, list the steering depths in `focal_distance_list` and set the window to the steering range. The default margin (the largest FLHM, often > 50 mm) otherwise admits skin spots the focus cannot reach.

??? note "Launcher inputs and the FLHM fit"
    | Input | Description |
    |-------|-------------|
    | `t1_filepath` | Subject T1.nii.gz (SimNIBS output) |
    | `simnibs_mesh_filepath` | Head mesh `sub-XXX.msh` ([CHARM](https://simnibs.github.io/simnibs/build/html/documentation/command_line/charm.html) output) |
    | `target_mask_filepath` | Binary [NIfTI](https://nifti.nimh.nih.gov/) sphere centred on the target |
    | `tx_config_yaml` | Transducer + placement parameters (written by PRESTUS) |

    The YAML is flat. Key fields: `max_distance`, `min_distance`, `optimal_distance`, `transducer_diameter`, `plane_offset`, `additional_offset`, `focal_distance_list`, `flhm_list`, five `weight_*` scalars. Without `optimal_distance` in the YAML, the launcher uses `(min_distance + max_distance) / 2`.

    PlanTUS fits FLHM against focal distance with a cubic, which needs at least 4 points. With fewer (e.g. the single O'Neil estimate when no axial profile CSV is given), the launchers (`plantus_compat.py`) interpolate linearly and hold FLHM constant outside the calibrated range. Non-finite pairs are dropped.

### 2  Convert SimNIBS Mesh to Surfaces

`PlanTUS.convert_simnibs_mesh_to_surfaces` (`convert_simnibs_mesh_to_surface` in PlanTUS 1.x) extracts `.surf.gii` surfaces from the `.msh` mesh. With PlanTUS 2.0 or later, the skin surface is smoothed with `placement.plantus.skin_smoothing_iterations` Taubin passes (default 2, as in PlanTUS) before any criterion is computed; this reduces mesh noise in the normals and so in the angle maps.

??? note "Surfaces extracted"
    | Surface | SimNIBS tags | File |
    |---------|-------------|------|
    | Skin | 1005 | `skin.surf.gii` |
    | Skull | 1007, 1008 | `skull.surf.gii` |

    SimNIBS also produces `skull.stl`, used for the skull-entry angle (4d).

### 3  Build the Avoidance Mask

Skin vertices over air cavities and near the eyes and ears are excluded from the composite score.

??? note "Excluded vertices and required tools"
    Skin vertices over air cavities and sinuses, within 30 mm of each eye centre and within 15 mm of the tragus are excluded from the composite score (`avoidance == 0`), and the mask is eroded by `transducer_radius` mm as a safety margin. The step needs `final_tissues.nii.gz` and `eeg_positions/Fiducials.csv` from the m2m folder; PlanTUS 2.0 builds the mask in Python, PlanTUS 1.x also needs `fslmaths`, `wb_command`, `mri_tessellate` and `mris_convert`. If any is missing, the step is skipped with a warning and all vertices remain eligible.

??? note "How the mask is built (`create_avoidance_mask`)"
    1. **Air cavities / sinuses**: vertices whose outward normal ray hits air in `final_tissues.nii.gz` (holes filled and the air surface extracted in Python in PlanTUS 2.0; with `wb_command -volume-fill-holes` and `mri_tessellate` in 1.x).
    2. **Eyes**: vertices within 30 mm of each eye centre (tag 1006, centres by k-means on x).
    3. **Ears (tragus)**: vertices within 15 mm of LPA/RPA from `eeg_positions/Fiducials.csv`, shifted 1.5 cm posterior.
    4. **Erosion**: the mask is eroded by `transducer_radius` mm.

### 4  Score Each Vertex

PlanTUS scores every skin vertex on five metrics: 4a to 4d here, 4e below.

??? note "Metrics 4a to 4d"
    | Metric | Definition | Better |
    |---|---|---|
    | 4a Skin–target distance `d_i` | Euclidean distance from vertex *i* to the target mask's centre of mass (`PlanTUS.roi_center_of_gravity`). The score uses the absolute deviation from `optimal_distance`; vertices outside `[min_distance, max_distance]` are excluded (see 5a). | Lower deviation |
    | 4b Skin–target angle `α_i` | Angle (degrees) between the skin normal and the vector to the target centre. | Lower |
    | 4c Skin–target intersection `x_i` | Path length of the 200 mm outward-normal ray inside the target mask STL (`PlanTUS.stl_from_nii`); chords are summed for several entry–exit pairs. | Higher |
    | 4d Skin–skull angle `β_i` | Angle between the skin normal and the skull normal where the 40 mm ray meets the skull. NaN (worst) where it meets no skull. | Lower |

#### 4e  Skull Thickness

Bone path length (mm) along the 4d ray.

??? note "Sampling and missing values"
    `final_tissues.nii.gz` (next to the SimNIBS mesh) is sampled every 0.25 mm (nearest voxel) up to 50 mm from the skin. Bone samples (labels 4, 7, 8: bone, compact, spongy) are counted until the first intracranial sample (labels 1–3: WM, GM, CSF). Lower is better.

    `t_i` is NaN (worst) where the ray meets no bone or reaches no intracranial tissue within 50 mm (face, neck, skull base). Without `final_tissues.nii.gz`, `t_i` is NaN everywhere, a warning is logged and `weight_skull_thickness` is ignored.

A few vertices get implausibly thin values (≤ 1 mm): neck rays reaching the spinal canal, and skin next to the ear where CHARM places bone or CSF directly under the skin. These vertices also have large skin–target angles and are rarely chosen unless thickness dominates the weights.

### 5a  Combine and Evaluate Scores

PlanTUS combines the normalised metrics into a weighted composite score; headless mode places the vertex with the minimum score.

??? note "The composite score"
    Each metric is min–max normalised to **[0, 1]** over its finite values. NaN or infinite values are set to 1 (worst), also for the inverted intersection metric.

    ```
    score_i = ( w_dist  × norm(|d_i − optimal_distance|)
              + w_ang   × norm(α_i)
              + w_inter × (1 − norm(x_i))
              + w_skull × norm(β_i)
              + w_thick × norm(t_i) )
            ÷ (w_dist + w_ang + w_inter + w_skull + w_thick)
    ```

    If `t_i` is unavailable (all NaN, or an older `plantus_metrics.mat` without it), `w_thick` is dropped from numerator and denominator. Vertices with `avoidance == 0`, `d_i < min_distance` or `d_i > max_distance` get `score = inf`. Headless mode places the vertex with the **minimum** score. Interactive mode reports that optimum before `wb_view` opens, and the score and rank of each placed vertex. A few outlying vertices (e.g. 40 mm skull-base rays for `t_i`) compress the spread and weaken that criterion below its nominal weight.

**Criteria evaluation** (`placement.plantus.evaluate`) re-scores the placement under 12 weight sets and flags settings under which a criterion cannot discriminate.

??? note "Criteria evaluation"
    **Criteria evaluation** (`placement.plantus.evaluate`, default `true`; needs `save_metrics`) re-scores the placement under 12 weight sets (configured, equal, one criterion at 0.55 with the rest sharing 0.45, and each criterion alone) and writes, next to the metrics file:

    | File | Contents |
    |------|----------|
    | `plantus_evaluation.json` | Placed vertex (rank, each criterion's value and % of eligible vertices that are better), per-criterion ranges over eligible vertices, the weight sweep and diagnostics; read by the [placement report](doc_placement.md#placement-report) |
    | `plantus_weight_sensitivity.csv` | One row per weight set: chosen vertex, its criteria, distance to the placed vertex |
    | `plantus_criteria_maps.png` | Each criterion and the composite score on the skin (4 views) |
    | `plantus_weight_sensitivity.png` | Criteria at the placed spot and at each weight set's choice; all chosen spots on the skin |

    Diagnostics flag settings under which a criterion cannot discriminate between positions:

    | Code | Condition |
    |------|-----------|
    | `distance_range_unconstrained` | `[min_distance, max_distance]` excludes < 5% of the skin outside the avoidance mask (e.g. a large `flhm_list` value) |
    | `intersection_sparse` | Beam–target intersection non-zero at < 2% of eligible vertices (small target mask) |
    | `<criterion>_unavailable` / `_nan` | Criterion undefined at all / > 25% of eligible vertices |
    | `placed_ineligible` | Placed vertex excluded under the configured scoring (interactive selection) |
    | `placed_not_optimum` | Placed vertex is not the composite optimum (interactive selection) |
    | `tradeoff_<criterion>` | At the placed vertex a weighted criterion is worse than at > 50% of eligible vertices |
    | `weight_sensitive` | A single-criterion-heavy weight set moves the chosen spot > 20 mm |

    Warnings also go to the placement log. Evaluation failures never abort a placement.

??? note "`plantus_metrics.mat` and offline re-scoring"
    `plantus_metrics.mat` (`placement.plantus.save_metrics`, default `true`) is written to `<output_dir>/plantus_metrics.mat`. It holds skin coordinates, normals and faces, the target centre, `d_i`, `α_i`, `x_i`, `β_i`, `t_i` (`skull_thickness`), the avoidance mask, the composite score, `best_vertex` (composite optimum), `selected_vertex` (vertex placed), `selection_mode` (`auto` headless, `gui` interactive, `manual` after `--do_only_trajectory`), `min_distance`, `max_distance`, `optimal_distance` and the weights. Vertex indices and faces are 0-based. Each run deletes an existing `plantus_metrics.mat` first, even with `save_metrics: false`.

    Use `compute_composite_score` in `plantus_common.py` to re-score under other weights without new ray casts. The evaluation runs `plantus_evaluation.py` on the metrics; to evaluate an existing metrics file, run `python functions/neuronav/plantus_evaluation.py <plantus_metrics.mat>`.

### 5b  Select a Vertex Interactively

With `placement.plantus.interactive = 1`, `wb_view` opens a scene with the maps 4a–4d, one per tab. The skull-thickness map (`skull_thickness_skin.func.gii`) is loaded too and can be chosen as an overlay in any tab. Each click replaces the selection and prints the vertex's criteria, composite score and rank (and the composite optimum) to the log.

??? note "Review dialog and requirements"
    Closing `wb_view` places the last selected vertex (step 6) and opens its placement scene (transducer and focus on skin and T1) in a second `wb_view`. A review dialog then offers *Use this position*, *Choose another vertex* (the planning scene reopens; metrics are not recomputed) or *Cancel* (the launcher exits non-zero and PRESTUS stops). Closing the planning scene without a click reviews the previous vertex, or aborts if there is none.

    Interactive mode requires:

    - A graphical display (desktop or VNC session; not a batch/HPC job). The MATLAB desktop is not needed: `placement.plantus.interactive = 1` also works under `matlab -batch` and from the Nexus
    - `pynput` (mouse listener); on Linux it needs an X display, on macOS Accessibility permission for the terminal/MATLAB
    - Connectome Workbench `wb_view`, on `PATH` or via `placement.plantus.connectome_wb_path`
    - PyQt5 for the review dialog (part of SimNIBS)

??? note "The `--selection` launcher option"
    `--selection` sets how a vertex is confirmed: `review` (the dialog above; PRESTUS uses it), `prompt` (`Generate transducer placement for vertex N? (yes/no)` after each click, for a terminal), `last_click` (no review) and `auto` (default: `prompt` when stdin is a terminal, else `review` if PyQt5 is available, else `last_click`). If no vertex is confirmed, the launcher exits non-zero. If `wb_view` exits at once without a click, the planning scene is reopened without the skull-thickness map.

### 6  Compute the Placement Geometry

For the selected vertex, PlanTUS computes the beam direction and the transducer centre (`prepare_acoustic_simulation`) and writes a 4 × 4 [Localite](https://www.localite.de/en/) position matrix (RAS mm).

??? note "Output files"
    Outputs go to `<output_dir>/vtx<N>/`, where `<output_dir>` is the per-target [placement run folder](doc_outputs.md#where-outputs-live) `<path.placement>/sub-XXX/plantus/<target_name>/`:

    | File | Contents |
    |------|----------|
    | `*_Localite.mat` | `position_matrix` (4×4, RAS mm), **used by PRESTUS** |
    | `*_Localite.txt` | Plain-text copy |
    | `*_kPlan.mat` | Same matrix in [k-Plan](https://brainbox-neuro.com/products/k-plan) convention (m, axes permuted) |
    | `*_Localite_XML.txt`, `*_Localite_dummyXML.txt` (2.x) | Fake Localite XML snippet (rewritten on each run) |
    | `*_Trajectory_Brainsight.txt`, `*_Trajectory_BabelBrain.txt` | Trajectory for Brainsight and BabelBrain (PlanTUS 2.0) |
    | `*Focus*`, `*TransducerModel*.surf.gii` | Focus and transducer surface models at the position |
    | `scene.scene` | Workbench placement scene (all modes; opened only by the Workbench picker) |

    PlanTUS 2.0 names the files `<mask>_vtx<N>_…`, PlanTUS 1.x `position_matrix_<mask>_vtx<N>_…`.

    `<output_dir>` itself holds the skin/skull surfaces, the metric maps (`distances_skin`, `distances_skin_thresholded`, `angles_skin`, `target_intersection_skin`, `skin_skull_angles_skin`, `skull_thickness_skin`, `avoidance_skin`; `.func.gii`), the planning `scene.scene` and `plantus_metrics.mat`.

??? note "Beam direction and transducer centre"
    For the selected vertex *v*:

    1. **Beam direction**: PlanTUS 2.0 aims from the skin vertex at the target centre (its default `target_centered` orientation); the PlanTUS page draws this axis. PlanTUS 1.x uses the negated skin normal if it intersects the target (`x_v > 0`), otherwise the negated skin-to-target vector. All are unit vectors pointing into the head.
    2. **Transducer centre**: `transducer_centre = skin_coordinate_v − (plane_offset + additional_offset) × beam_direction`. `plane_offset` is the distance from the bowl apex to the exit plane (`curv_radius_mm − dist_geom_ep_mm`); `additional_offset` is `placement.skin_gap_mm`, the gap between skin and exit plane.
    3. **Localite position matrix** (`PlanTUS.create_Localite_position_matrix`): column 3 = `transducer_centre`, column 0 = `beam_direction`, columns 1–2 by Gram–Schmidt.

### 7  Read the Output into PRESTUS

`position_transducer_plantus.m` loads the most recent `*Localite.mat` under `<path.placement>/sub-XXX/plantus/<target_name>/` written during the current run. If there is none, it errors (`PRESTUS:plantus:noOutput`). `localite_matrix_to_positions` converts `position_matrix` into `trans_pos` and `focus_pos` in voxel coordinates. A T1 overlay plot (`plot_placement_t1_overlay`) of the transducer and target is written to the subject output folder.

## Coordinate Conventions

| Frame | Units | Description |
|-------|-------|-------------|
| Voxel / grid | integer indices | PRESTUS simulation grid |
| RAS+ (subject) | mm | T1 image world space; used by Localite matrix |
| MNI | mm | Template space; target specified here, mapped to RAS via SimNIBS |
| k-Plan | m, permuted axes | Not used by PRESTUS |

## External Dependencies

| Tool | Automatic / chosen vertex | Workbench | Used in step |
|------|:--------:|:-----------:|-------------|
| SimNIBS Python ([nibabel](https://nipy.org/nibabel/), [scipy](https://scipy.org/), [nilearn](https://nilearn.github.io/stable/index.html), vtk, h5py), see [installation guide](doc_installation.md#install-simnibs) | ✓ | ✓ | All steps |
| PlanTUS Python library (`PlanTUS.py`) | ✓ | ✓ | All steps |
| [trimesh](https://trimesh.org/), pandas (PlanTUS 2.0) | ✓ | ✓ | Surfaces, avoidance mask, target model |
| [PyYAML](https://pyyaml.org/) (or built-in fallback) | ✓ | ✓ | YAML config loading |
| FSL `fslmaths` (PlanTUS 1.x) | ✓ | ✓ | Avoidance mask: tissue binarisation |
| Connectome Workbench `wb_command` (PlanTUS 1.x) | ✓ | ✓ | Surface metrics, avoidance mask |
| FreeSurfer `mri_tessellate` / `mris_convert` (PlanTUS 1.x) | ✓ | ✓ | Avoidance mask: air mesh |
| Connectome Workbench `wb_view` | – | ✓ | Interactive scene display |
| `pynput` Python package | – | ✓ | Mouse click listener for wb_view |
| PyQt5 (part of SimNIBS) | – | ✓ | Review dialog |

The launchers detect the PlanTUS version from `code/PlanTUS.py` and check for the tools that version needs before they start.

Install trimesh (PlanTUS 2.0) and `pynput` (Workbench picker) in the SimNIBS [conda](https://docs.conda.io/) environment:

```bash
conda activate simnibs_v4.6.0
pip install trimesh pynput
```
