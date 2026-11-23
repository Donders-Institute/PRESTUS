# PlanTUS Placement

![Pipeline steps; this page belongs to 3 Placement](assets/figures/fig_stage-strip-placement_diagram.svg){ .stage-strip }

The entry point decides which part of the skull the beam crosses and whether the focus can reach the target. Given a T1-weighted image, a [SimNIBS](https://simnibs.github.io/simnibs/) head mesh, a target-ROI mask, and a transducer configuration, [PlanTUS](https://github.com/mlueckel/PlanTUS) (`placement.mode: 'plantus'`) scores every skin vertex as a candidate transducer position on skin–target distance, skin and skull angles, skull thickness, and beam–target overlap. PRESTUS places the transducer at the vertex with the best weighted score, or at a vertex you choose on the maps on the Nexus PlanTUS page (or in [Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench)).

On this page:

- [Running PlanTUS](#running-plantus): the steps of a placement run.
- [Choosing the Position](#choosing-the-position): automatic, on the Nexus PlanTUS page, or in Workbench.
- [External Dependencies](#external-dependencies): the tools each mode needs.

!!! note "At a glance"
    - **When you need it:** When you want an optimised, documented choice of entry point, or to compare focal depths.
    - **Before you start:** The full SimNIBS mesh and a PlanTUS installation, 2.0 or later ([Installation](doc_installation.md#plantus-dependencies)).

!!! warning "Experimental"
    PlanTUS placement is under active development. Verify results manually before use in a study.

!!! warning "Limitations"
    - The Nexus PlanTUS page needs the maps of an earlier run (`save_metrics: true`) and draws them where the Nexus runs; the run itself can be headless.
    - Choosing in Workbench needs a display (desktop or VNC), Connectome Workbench `wb_view`, and `pynput`. On macOS, the terminal or [MATLAB](https://www.mathworks.com/products/matlab.html) needs Accessibility permission for `pynput`. Without a display (Linux), PlanTUS runs headless.
    - The Workbench review dialog needs [PyQt5](https://www.riverbankcomputing.com/software/pyqt/) in the SimNIBS [Python](https://www.python.org/). It ships with SimNIBS; without it, the last clicked vertex is placed without review.
    - To see errors from the launcher directly, run `prestus_plantus_gui_launcher.py` from a terminal (step 5b).

## Running PlanTUS

1. Set `placement.mode: plantus`, the `placement.plantus` paths and target, and the calibrated `focal_distance_list` / `flhm_list`.

    ```yaml
    placement:
      mode: plantus
      plantus:
        script_path: /path/to/PlanTUS     # contains code/PlanTUS.py
        target_name: target               # names the output folder
        mni_target_mm: [-29, -6, 7]       # example; empty = transducer.focus_pos
    ```

2. In the Nexus, click *Compute maps* on the Placement card (or on Placement ▸ PlanTUS). PlanTUS computes the criterion maps (a few minutes) and places the transducer at the composite-score optimum; no simulation runs. In MATLAB, run the pipeline with `placement.plantus.vertex: []`.
3. Open Placement ▸ PlanTUS. The skin is coloured by the composite score; choose a criterion under *Map* to see one map. Click a skin position: the panel lists its criteria, composite score, and rank, and the view shows the beam axis and the transducer there, in the housing Placement ▸ 3D draws.
4. Click *Place here*. The Nexus sets `placement.plantus.vertex` and places the transducer at that vertex on the same maps, without recomputing them. The config keeps the choice, so a later simulation run uses it.
5. Check the [placement report](doc_placement.md#placement-report) and the Viewer. Section [5a](doc_placement_plantus_pipeline.md#5a-combine-and-evaluate-scores) explains the evaluation it contains.

## Choosing the Position

`placement.plantus.vertex` and `placement.plantus.interactive` (`[]`, the default, follows `simulation.interactive`) select how the vertex is chosen:

| Setting | Script | Vertex selection |
|------|--------|-----------------|
| **Automatic** (`vertex: []`, `interactive: 0`) | `prestus_plantus_launcher.py` | Weighted composite score; no display required |
| **Chosen vertex** (`vertex: N`) | `prestus_plantus_launcher.py --do_only_trajectory N` | Vertex `N` (0-based) of the maps of an earlier run in the same placement run folder, chosen on the Nexus PlanTUS page; the maps are not recomputed |
| **Workbench** (`vertex: []`, `interactive: 1`) | `prestus_plantus_gui_launcher.py` | The five criterion maps in Connectome Workbench; you click a vertex, close `wb_view`, see the transducer there, and confirm it in a review dialog |

Automatic and Workbench runs compute the same maps and composite score (steps 2 to 5a) and write the same files. A chosen vertex is placed only while the target is where the maps were computed for: if the target moved more than 1 mm, the run stops and asks you to clear `placement.plantus.vertex` and compute the maps again. A chosen vertex belongs to one subject's maps; for a study, give it per run in a manifest column. Only Workbench runs need `pynput`, `wb_view`, and PyQt5.

Both modes share all code except vertex selection (`plantus_common.py`). The scripts replicate `PlanTUS_wrapper.py` (Mainz group, at `placement.plantus.script_path`) but take CLI arguments instead of hard-coded paths.

How a run computes its criteria, composite score, and placement is described in [PlanTUS: How It Works](doc_placement_plantus_pipeline.md).

## External Dependencies

| Tool | Automatic / chosen vertex | Workbench | Used in step |
|------|:--------:|:-----------:|-------------|
| SimNIBS Python ([nibabel](https://nipy.org/nibabel/), [scipy](https://scipy.org/), [nilearn](https://nilearn.github.io/stable/index.html), vtk, h5py), see [installation guide](doc_installation.md#install-simnibs) | ✓ | ✓ | All steps |
| PlanTUS Python library (`PlanTUS.py`) | ✓ | ✓ | All steps |
| [trimesh](https://trimesh.org/), pandas | ✓ | ✓ | Surfaces, avoidance mask, target model |
| [PyYAML](https://pyyaml.org/) (or built-in fallback) | ✓ | ✓ | YAML config loading |
| Connectome Workbench `wb_view` | – | ✓ | Interactive scene display |
| `pynput` Python package | – | ✓ | Mouse click listener for wb_view |
| PyQt5 (part of SimNIBS) | – | ✓ | Review dialog |

The launchers check for these tools before they start.

Install trimesh, `potpourri3d`, and `pynput` (Workbench picker) in the SimNIBS [conda](https://docs.conda.io/) environment. PlanTUS uses `potpourri3d` for exact geodesic distances when it shrinks the avoidance mask; without it, PlanTUS falls back to an approximation, so the mask can differ slightly between installs:

```bash
conda activate simnibs_v4.6.0
pip install trimesh potpourri3d pynput
```
