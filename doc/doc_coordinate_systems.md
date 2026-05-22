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
| MNI space | mm | `mni` | MNI152 standard space. **Note: MNI uses LAS orientation (smaller X = right hemisphere), which is the opposite of RAS.** The axis flip is handled explicitly in `neuronav_convert_native_to_MNI`. |
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

## MNI axis-flip warning

MNI space uses **LAS** orientation: the X axis increases toward the *left* hemisphere (smaller X = right). Subject RAS space has the opposite convention (smaller X = left). This flip must be accounted for when interpreting laterality in MNI-space outputs. See `neuronav_convert_native_to_MNI` lines ~149–156 for the explicit correction.

---

## Available transforms (SimNIBS)

PRESTUS reuses nonlinear warp fields produced by SimNIBS/CHARM:

- `Conform2MNI_nonl.nii.gz` — subject native → MNI (displacement field, mm)
- `MNI2Conform_nonl.nii.gz` — MNI → subject native (inverse displacement field, mm)

The segmentation (`final_tissues.nii.gz`) is the spatial anchor for these warps. If the planning T1 differs from the segmentation T1, the affine mismatch is handled in `neuronav_convert_native_to_MNI`.

---

## NIfTI affine convention (MATLAB)

`niftiinfo()` stores the voxel-to-world affine in `Transform.T` (column-major).

- World from voxel: `world = T' * [vox; 1]`
- Voxel from world: `vox   = T' \ [world; 1]`  (used in `ras_to_grid`)

Always use `ras_to_grid()` rather than calling `T' \` directly.
