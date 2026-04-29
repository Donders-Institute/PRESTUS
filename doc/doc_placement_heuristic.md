# Heuristic Transducer Placement

The function `transducer_positioning` identifies heuristic locations for transducer placement. It can be called with `transducer_positioning_start`.

The function targets brain regions specified by MNI coordinates, converting them to subject-native space using SimNIBS tools. It aims to identify candidate transducer positions where the geometric focus aligns with a target within a user-defined focal distance range. For each position, it computes the intersection proportion between the transducer and skin, mean/variance distances to skin and skull, and geometric focus/exit plane positions. It (optionally removes) areas including ears and selects a heuristic position based on user-defined criteria. Besides target MNI coordinates, the function does not require manual intervention, enabling automatic end-to-end workflows.

For parameters, see the [overview](doc_parameters.md#heuristic-transducer-placement).

### Candidate positions

- Convert target from MNI (mm) to subject grid space (voxels)
- Find candidate transducer positions on skull (expanding sphere)
- Calculate criteria for each location: target distance, intersection with skin, mean distance to skull, variance in distance to skin & skull

#### Table: Candidate Coordinates

The `tpars_sub-XXX_target.csv` table contains transducer position candidates with the following columns. This information can be used to select a suitable candidate (e.g., by selecting a position with low intersection of the transducer and head tissue (prop_intersect < 0.05) and minimal Euclidean distance to the target; often supported by visual inspection).

| Metric          | Description                                    | Interpretation  |
|-----------------|------------------------------------------------|-----------------|
|`idx`            | index                                          | |
|`trans_x`        | Transducer x-coordinate (voxels)               | |
|`trans_y`        | Transducer y-coordinate (voxels)               | |
|`trans_z`        | Transducer z-coordinate (voxels)               | |
|`targ_x`         | Target x-coordinate (voxels)                   | |
|`targ_y`         | Target y-coordinate (voxels)                   | |
|`targ_z`         | Target z-coordinate (voxels)                   | |
|`dist_to_target` | Euclidean distance transducer-to-target (voxels)| |
|`prop_intersect` | Fraction of voxels in the transducer’s orthogonal plane that intersect with the skull mask. | Low values (≈ 0) indicate minimal skull obstruction; high values (> 0.3) indicate significant bone interference. |
|`mean_dist_skin` | Mean distance from **non‑intersecting** voxels to the skin boundary. | Indicates average proximity of the plane to the skin surface. |
|`var_dist_skin`  | Variance of distances from non‑intersecting voxels to the skin boundary. | Higher variance → surface curvature or oblique intersection. |
|`mean_dist_skull`| Mean distance from **all** voxels in the plane to the nearest skull boundary. | Shows average clearance between the plane and skull. |
|`var_dist_skull` | Variance of distances to the skull boundary. | Quantifies irregularity of skull spacing within the beam footprint. |


### Remove ear locations

**[Optional]**  

For practical reasons, transducer locations overlapping with the ears are not desired. Possible locations that overlap with ear positions can optionally be removed within an `tp_ear_radius` (Radius of the nogo zone, mm) around `tp_left_ear_center` (approximate coordinates for left ear, voxels) and `tp_right_ear_center`.

### Select heuristic transducer position

The desired criterion for intersection with the skin can be defined via `tp_criterion_intersection`. If no placement is found within this criterion, the intersection will be iteratively expanded by `tp_expand_step` (default: 1%) until a match is identified.

| Criterion | Default | Explanation |
|:----------|---------------|-------------|
| `tp_criterion_intersect` | < 0.05 | [Fraction] Reduce skull traversal |
| `tp_criterion_skin_mean` | NaN  | [Quantile] Short EP distance to skin |
| `tp_criterion_skull_mean` | NaN | [Quantile] Short EP distance to skull |
| `tp_criterion_skin_var` | NaN | [Quantile] Coherent transducer exit plane toward skin (e.g., little curvature) |
| `tp_criterion_skull_var` | NaN | [Quantile] Coherent transducer exit plane toward skull surface (Irregular bone → phase distortion → incoherent wavefront) |

Amongst locations fulfilling the above criteria, the location with a minimum distance to the target is selected.

An an output, this function will create a table that contains the coordinates in **voxel (grid space)**, **mm (Localite space)**, and **mm (RAS+ space)** (based on the image header of the segmentation image).

> IMPORTANT: Localite mm coordinates assume that Localite planning image has `canonical_affine_transform` applied (see [below](#save-localite-planning-image))! The Localite coordinates are simply rescaled by the specified voxel sizes in the header assuming no affine transformations. Always visually verify coordinates ...

**Example coordinate output table**:
<pre style="font-size: 0.75em; line-height: 1;">
trans_x	161     
trans_y	138     
trans_z	246     
targ_x	117     
targ_y	139     
targ_z	161     

<...>       

LOCALITE_Transducer_mm	Inf       
loc_trans_x_mm	144.9       
loc_trans_y_mm	123.09      
loc_trans_z_mm	221.02      

LOCALITE_Target_mm	Inf     
loc_targ_x_mm	105.3
loc_targ_y_mm	124.88
loc_targ_z_mm	144.65

RAS_Transducer_mm	Inf
ras_trans_x_mm	43.35
ras_trans_y_mm	-16.17
ras_trans_z_mm	82.3

RAS_Target_mm	Inf
ras_targ_x_mm	5.45
ras_targ_y_mm	-20.71
ras_targ_z_mm	5.19
</pre>

### Plot heuristic transducer position

PRESTUS generates an overview of the selected transducer placement:

![ex_heuristic_placement](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/ex_heuristic_placement.png)

If `parameters.localite_path` is specified, it will also deposit a copy of the plot there.

### Save Localite planning image

**[Optional]**  

 Localite can struggle with images that contain affine matrices in the image header. PRESTUS offers the optional step (`tp_save_localiteT1`) of depositing a header-adjusted T1 planning image for Localite (based on the `T1.nii.gz` in the SimNIBS output m2m directory). For this, `parameters.localite_path` has to be specified.
