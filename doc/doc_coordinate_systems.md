# Coordinate Reference Systems

## Canonical convention

**RAS (Right-Anterior-Superior) mm is the canonical *entry-point* space in PRESTUS.**

Localite markers and any externally specified positions are expressed in RAS mm. These are immediately converted to T1 voxel indices via `ras_to_grid()` and stored in `parameters.transducer.trans_pos` / `focus_pos`. All internal processing downstream of that conversion operates in voxel indices.

The internal simulation grid (passed to k-Wave) uses 1-based MATLAB voxel indices and is anatomy-agnostic. This is intentional: k-Wave operates in grid space, not world space.

---

## Named spaces

| Name | Unit | Label | Description |
|------|------|-------|-------------|
| Simulation grid | voxels (1-based) | `grid` | k-Wave internal space. Axes align with the T1 after preprocessing (rotation to focal axis etc.). Not directly tied to anatomy unless a planning image is loaded. |
| Subject / scanner space | mm | `ras_plus` | RAS mm in the subject's native scanner frame. PRESTUS's entry-point space: Localite markers and external coordinates arrive in RAS mm and are immediately converted to T1 voxels via `ras_to_grid()`. |
| Localite planning image | mm | `ras_plus` (adjusted header) | Subject space expressed via the diagonal affine written by `canonical_affine_transform`. Voxel-to-world mapping is `world = voxel * voxel_size_mm` (no rotation). See [Localite placement](doc_placement_heuristic.md#save-localite-planning-image). |
| MNI space | mm | `mni` | MNI152 standard space, standard world-mm convention (negative X = left, positive X = right — same laterality sense as RAS). The template's voxel *array* happens to be stored in LAS order, but going through the NIfTI affine (as `ras_to_grid()`/`transform_coordinates()` do) already accounts for that. See [MNI voxel storage order vs. world-mm laterality](#mni-voxel-storage-order-vs-world-mm-laterality). |
| MNI voxel | voxels (1-based) | *(unnamed)* | Integer indices into an MNI template volume (e.g. 181×217×181). Convert via the template's NIfTI affine. |

---

## Coordinate flow for patient (layered) simulations

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

**Key rule:** `ras_to_grid()` is the sole boundary between RAS mm and T1 voxel indices. Downstream of that call, everything is voxels. Do not replicate affine inversion elsewhere.

---

## Coordinate flow for phantom / free-water simulations

Phantom and water grids have no intrinsic patient-space anchor. Grid indices are arbitrary. Transducer and focus positions are either:

- specified directly as grid voxel indices via `parameters.transducer.trans_pos` / `focus_pos`, or
- placed automatically (transducer at near face, focus at expected focal distance).

**Optional: anchoring a phantom/water grid to RAS space**

If you need phantom outputs to be spatially registered (e.g., to overlay a simulated focal profile on a planning image, or to compare with k-Plan), set:

```yaml
grid:
  origin_ras_mm: [x, y, z]   # RAS mm of voxel [1,1,1]
```

When set, PRESTUS writes output NIfTIs with a diagonal RAS affine built from `origin_ras_mm` and `grid.resolution_mm`. The affine is:

```
T = [ res   0    0   ox ]
    [  0   res   0   oy ]
    [  0    0   res  oz ]
    [  0    0    0    1 ]
```

This makes the outputs directly openable in ITK-SNAP, 3D Slicer, and k-Plan with correct spatial positioning. It does not affect the simulation itself.

---

## MNI voxel storage order vs. world-mm laterality

The MNI152 template's voxel array is stored in **LAS** order internally (its sform has a negative
X scale factor), but this is purely a storage-order detail: world-mm coordinates derived through the
NIfTI affine (as `ras_to_grid()` and `transform_coordinates()` always do) already come out in
standard MNI convention — **negative X = left hemisphere, positive X = right** — with no extra
manual flip required. This was confirmed empirically: a labeled left-hemisphere target's
peak-intensity voxel in `sub-601_layered_MNI_left_dcPUT_..._intensity.nii.gz` transforms to X =
-69.0 mm via the plain sform, correctly on the left.

`neuronav_convert_native_to_MNI` (~lines 159-171) has a separate, unrelated `side_order`-conditional
block that forcibly pins the *transducer's* (not the target's) MNI voxel X-coordinate to the FOV
edge (`min_voxel_mni`/`max_voxel_mni`). This is a visualization/robustness workaround for the
relatively cropped out-of-skull MNI window, not a laterality correction — it only runs when the
optional `side_order` argument is supplied, and it never touches the target coordinate. Do not
treat it as "the" fix for MNI laterality; laterality is already correct by construction once you go
through the affine.

---

## Localite planning image: zero-origin affine and SimNIBS voxel grid

Neuronavigation software (e.g. Localite TMS Navigator) writes a planning T1 with a **diagonal, zero-origin affine**:

```
affine = diag([voxel_size_mm, voxel_size_mm, voxel_size_mm, 1])
```

This means the "RAS mm" values stored in Localite XML trigger matrices are actually **planning-image space coordinates** (`world = voxel × voxel_size_mm`), not scanner-isocenter-relative RAS mm. When `ras_to_grid()` is called with the planning T1 header, it correctly inverts this back to voxel indices:

```
voxel = localite_mm / voxel_size_mm
```

**Do not re-project Localite coordinates using the SimNIBS m2m T1 header.**  
The m2m T1 (`m2m_{sub_id}/T1.nii.gz`) has a full world-space affine (non-zero origin, possible small rotation). Calling `ras_to_grid(localite_mm, simnibs_t1_header)` applies the wrong inverse-affine and shifts positions by tens of millimetres.

**When the same T1 was used for both neuronavigation and SimNIBS input**, the m2m T1 and the planning T1 share the same voxel grid — SimNIBS preserves the voxel dimensions and layout, updating only the affine to describe world-space position. In this case, voxel indices from the planning T1 (`ras_to_grid(localite_mm, planning_t1_header)`) are directly valid for PRESTUS simulations that load the m2m T1. No coordinate reprojection is needed.

If the planning T1 and the SimNIBS input T1 differ (e.g. different acquisitions, different resolution), a registration step is required before the voxel indices can be used with the m2m T1.

---

## Available transforms (SimNIBS)

PRESTUS reuses nonlinear warp fields produced by SimNIBS/CHARM:

- `Conform2MNI_nonl.nii.gz` — subject native → MNI (displacement field, mm)
- `MNI2Conform_nonl.nii.gz` — MNI → subject native (inverse displacement field, mm)

The segmentation (`final_tissues.nii.gz`) is the spatial anchor for these warps. If the planning T1 differs from the segmentation T1, the affine mismatch is handled in `neuronav_convert_native_to_MNI`.

For point (coordinate) transforms, `mni2subject_coords_LDfix` / `subject2mni_coords_LDfix` wrap the SimNIBS CLI and accept a transform type: `nonl` (nonlinear warp), `12dof` (linear affine), or `6dof` (rigid). The nonlinear warp is most accurate **inside the brain** but is unreliable/extrapolated **outside** it. Accordingly, `mni` transducer placement (see [MNI placement](doc_placement_mni.md)) converts the in-brain focus via `mni → nonl warp → RAS+ → grid` (as the heuristic), and the out-of-brain transducer point via `mni → 12dof affine → RAS+ → grid`.

---

## NIfTI affine convention (MATLAB)

`niftiinfo()` stores the voxel-to-world affine in `Transform.T` (column-major).

- World from voxel: `world = T' * [vox; 1]`
- Voxel from world: `vox   = T' \ [world; 1]`  (used in `ras_to_grid`)

Always use `ras_to_grid()` rather than calling `T' \` directly.
