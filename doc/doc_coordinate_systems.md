# Coordinate Reference Systems

A position given in the wrong coordinate space puts the transducer or the focus in the wrong part of the head, in the case of a wrong image header by tens of millimetres. PRESTUS converts every externally given position once to T1 voxel indices and maps results back to T1 and MNI space. This page defines each space and the functions that convert between them.

On this page:

- [Canonical Coordinate Convention](#canonical-coordinate-convention): the entry-point space and where conversion happens.
- [Coordinate Spaces](#coordinate-spaces): voxel, RAS, Localite, and MNI spaces.
- [Coordinate Flow for Head Models](#coordinate-flow-for-head-models): from neuronavigation to simulation and back.
- [Coordinate Flow for Phantoms and Water](#coordinate-flow-for-phantoms-and-water): positions and an optional RAS affine for phantom grids.
- [MNI Voxel Order and Laterality](#mni-voxel-order-and-laterality): why MNI laterality needs no manual flip.
- [Localite Planning Image Header](#localite-planning-image-header): which header to use for Localite coordinates.
- [SimNIBS Warps and Transforms](#simnibs-warps-and-transforms): the SimNIBS warps and point transforms.
- [NIfTI Affine in MATLAB](#nifti-affine-in-matlab): 1-based voxels and the NIfTI affine.

!!! note "At a glance"
    - **When you need it:** When positions end up in the wrong place, or before changing any coordinate transform in the code.
    - **Before you start:** Nothing; for the basic idea see [Glossary](doc_glossary.md#head-model).

## Canonical Coordinate Convention

RAS (Right-Anterior-Superior) mm is the entry-point space in PRESTUS. [Localite](https://www.localite.de/) markers and any externally specified positions arrive in RAS mm. They are converted once to T1 voxel indices via `ras_to_grid()` and stored in `parameters.transducer.trans_pos` / `focus_pos`. Everything downstream of that call works in voxel indices. The simulation grid passed to [k-Wave](https://github.com/ucl-bug/k-wave) uses 1-based [MATLAB](https://www.mathworks.com/products/matlab.html) voxel indices and is anatomy-agnostic, because k-Wave operates in grid space, not world space.

## Coordinate Spaces

| Name | Unit | Label | Description |
|------|------|-------|-------------|
| Simulation grid | voxels (1-based) | `grid` | k-Wave internal space: the T1 voxel grid rescaled to `grid.resolution_mm`, rotated per `grid.mode` (`focal_axis`, the default, and `transducer_axis` turn the beam axis onto z; `ras_plus` keeps the T1 axes) and cropped. `planimg.transf` maps T1 voxels to this grid. |
| Subject / scanner space | mm | `ras_plus` | RAS mm in the subject's native scanner frame. PRESTUS's entry-point space: Localite markers and external coordinates arrive in RAS mm and are converted to T1 voxels via `ras_to_grid()`. |
| Localite planning image | mm | `ras_plus` (adjusted header) | Subject space expressed via the diagonal affine written by `canonical_affine_transform`. Voxel-to-world mapping is `world = (voxel − 1) × voxel_size_mm` for 1-based PRESTUS voxels (no rotation). See [Save Localite planning image](doc_placement_heuristic.md#localite-planning-image). |
| MNI space | mm | `mni` | MNI152 standard space of FSL's 1 mm template (`MNI152NLin6Asym`, the template SimNIBS ships), standard world-mm convention (negative X = left, positive X = right, the same laterality sense as RAS). The template's voxel *array* is stored in LAS order, but going through the [NIfTI](https://nifti.nimh.nih.gov/) affine (as `ras_to_grid()`/`transform_coordinates()` do) accounts for that. See [MNI voxel order and laterality](#mni-voxel-order-and-laterality). |
| MNI voxel | voxels (1-based) | *(unnamed)* | Integer indices into an MNI template volume (e.g. 182×218×182 for [SimNIBS](https://simnibs.github.io/simnibs/)'s `toMNI/final_tissues_MNI.nii.gz`). Convert via the template's NIfTI affine. |

## Coordinate Flow for Head Models

```
Localite XML
  → neuronav_compute_series_statistics   [RAS mm, native]
  → localite_matrix_to_positions         [RAS mm]  → ras_to_grid() → [T1 voxels]
                                          stored in parameters.transducer.trans_pos / focus_pos
  → grid_transducer_location             [T1 voxels] → planimg.transf → [sim-grid voxels]
  → k-Wave simulation                    [sim-grid voxels]
  → nifti_to_t1w                         [T1 voxels]    ← planimg.inv_transf
  → nifti_to_mni                         [MNI voxels]   ← SimNIBS warp
```

**Key rule:** `ras_to_grid()` / `grid_to_ras()` are the sole boundary between RAS mm and T1 voxel indices (1-based). Downstream of that call, everything is voxels. Do not replicate affine inversion elsewhere: for several points or a sub-voxel (unrounded) index, use `ras_to_voxel()`, which `ras_to_grid()` itself calls. Other 4×4 affines (grid transforms, SimNIBS matrices) are inverted with `affine_inverse()`, which reads both conventions and keeps the fixed row or column exact.

## Coordinate Flow for Phantoms and Water

Phantom and water grids have no patient-space anchor, so grid indices are arbitrary. Transducer and focus positions are either specified directly as grid voxel indices via `parameters.transducer.trans_pos` / `focus_pos`, or placed automatically (transducer at the near face, focus at the expected focal distance).

To register phantom outputs spatially (for example to overlay a simulated focal profile on a planning image, or to compare with [k-Plan](https://brainbox-neuro.com/products/k-plan)), set:

```yaml
grid:
  origin_ras_mm: [x, y, z]   # RAS mm of voxel [1,1,1]
```

PRESTUS then writes output NIfTIs with a diagonal RAS affine built from `origin_ras_mm` and `grid.resolution_mm`, so they open in [ITK-SNAP](http://www.itksnap.org/), [3D Slicer](https://www.slicer.org/) and k-Plan with correct spatial positioning. The simulation itself is unaffected.

??? note "The phantom affine"
    ```
    T = [ res   0    0   ox ]
        [  0   res   0   oy ]
        [  0    0   res  oz ]
        [  0    0    0    1 ]
    ```

## MNI Voxel Order and Laterality

The MNI152 template's voxel array is stored in **LAS** order (its sform has a negative X scale factor). This is a storage-order detail only: world-mm coordinates derived through the NIfTI affine (as `ras_to_grid()` and `transform_coordinates()` always do) come out in standard MNI convention, **negative X = left hemisphere, positive X = right**, with no manual flip.

??? note "The `side_order` clamp in `neuronav_convert_native_to_MNI`"
    `neuronav_convert_native_to_MNI` has a separate, unrelated `side_order`-conditional block that pins the *transducer's* (not the target's) MNI voxel X-coordinate to the FOV edge (`min_voxel_mni`/`max_voxel_mni`). This is a visualisation/robustness workaround for the cropped out-of-skull MNI window, not a laterality correction. It only runs when the optional `side_order` argument is supplied, and it never touches the target coordinate. Laterality is correct by construction once you go through the affine.

## Localite Planning Image Header

Neuronavigation software (e.g. Localite TMS Navigator) writes a planning T1 with a **diagonal, zero-origin affine**:

```
affine = diag([voxel_size_mm, voxel_size_mm, voxel_size_mm, 1])
```

The "RAS mm" values stored in Localite XML trigger matrices are therefore **planning-image space coordinates** (`world = voxel_0 × voxel_size_mm`, with `voxel_0` the 0-based NIfTI index), not scanner-isocenter-relative RAS mm. `ras_to_grid()` called with the planning T1 header inverts this back to the 1-based PRESTUS voxel (`voxel = round(localite_mm / voxel_size_mm) + 1`).

!!! warning "Do not re-project Localite coordinates using the SimNIBS m2m T1 header"
    The m2m T1 (`m2m_{sub_id}/T1.nii.gz`) has a full world-space affine (non-zero origin, possible small rotation). Calling `ras_to_grid(localite_mm, simnibs_t1_header)` applies the wrong inverse affine and shifts positions by tens of millimetres.

When the same T1 was used for both neuronavigation and SimNIBS input, the m2m T1 and the planning T1 share the same voxel grid: SimNIBS preserves the voxel dimensions and layout and updates only the affine. Voxel indices from the planning T1 (`ras_to_grid(localite_mm, planning_t1_header)`) are then directly valid for PRESTUS simulations that load the m2m T1, with no reprojection. If the planning T1 and the SimNIBS input T1 differ (different acquisitions, different resolution), a registration step is required before the voxel indices can be used with the m2m T1.

## SimNIBS Warps and Transforms

PRESTUS reuses nonlinear warp fields produced by SimNIBS/[charm](https://simnibs.github.io/simnibs/build/html/documentation/command_line/charm.html):

- `Conform2MNI_nonl.nii.gz`: subject native → MNI (displacement field, mm)
- `MNI2Conform_nonl.nii.gz`: MNI → subject native (inverse displacement field, mm)

The segmentation (`final_tissues.nii.gz`) is the spatial anchor for these warps. If the planning T1 differs from the segmentation T1, the affine mismatch is handled in `neuronav_convert_native_to_MNI`.

For point (coordinate) transforms, `mni2subject_coords_LDfix` / `subject2mni_coords_LDfix` wrap the SimNIBS CLI with charm's nonlinear warp (`nonl`), the only MNI transform SimNIBS 4 writes. The warp is most accurate **inside the brain** and less exact near the scalp, so `mni` transducer placement (see [MNI placement](doc_placement_mni.md)) converts both points via `mni → nonl warp → RAS+ → grid` and then snaps the transducer point to the scalp.

## NIfTI Affine in MATLAB

`niftiinfo()` stores the NIfTI voxel-to-world affine in `Transform.T` (transposed, row-vector convention). Like the NIfTI standard, it takes **0-based** voxel indices, whereas PRESTUS voxels are **1-based** MATLAB indices into the `niftiread` array:

- World from PRESTUS voxel: `world = T' * [vox - 1; 1]`  (`grid_to_ras`)
- PRESTUS voxel from world: `vox = round(T' \ [world; 1]) + 1`  (`ras_to_grid`; `ras_to_voxel` without the rounding)

PRESTUS reorders the T1 and segmentation voxels to RAS+ when it loads them (axis permutations and flips, no resampling). The affine keeps any oblique rotation, so every voxel keeps its world position.

Always use `ras_to_grid()` / `grid_to_ras()` rather than applying `Transform.T` to voxel indices directly. Affines composed on 1-based grids (e.g. `grid_to_ras_affine`) must include the 1-based → 0-based shift before `Transform.T`. Such an affine becomes a NIfTI header via `nifti_affine_from_1based`.
