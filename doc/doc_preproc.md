# Head preprocessing

!!! note "At a glance"
    - **What:** How the SimNIBS segmentation becomes the simulation grid: rotation, rescaling, smoothing, skull hole filling and cropping.
    - **When you need it:** Rarely: the defaults suit most studies. Read it when the skull looks wrong in the report, or to understand why the grid has the size it has.
    - **Before you start:** A `layered` simulation with a SimNIBS segmentation.

PRESTUS preprocesses `layered` segmentations to set up simulation grids.

![preprocessing](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/preprocessing.png)

### Orientation

PRESTUS aligns the simulation grid before cropping. The alignment strategy is controlled by `grid.mode`:

| Mode | Description |
|---|---|
| `focal_axis` (default) | Rotates the volume so the mean focal unit vector of all transducers aligns with z. For N=1 this equals `transducer_axis`; for bilateral setups it equals the inter-transducer axis (same result as the old `principal_axis`); for angled multi-transducer setups z correctly points toward the beam convergence zone. Falls back to transducer(1) focal axis if `norm(mean focal vector) < grid.focal_axis_min_norm` (default 0.5). Implemented in `preproc_align_to_mean_focal_axis`. |
| `transducer_axis` | Rotates the volume so the focal axis of transducer(1) aligns with z. Useful when transducer(1) is the primary beam and other transducers are secondary. Implemented in `preproc_align_to_focal_axis`. |
| `ras_plus` | No rotation; only isotropic rescaling to the requested grid resolution. Preserves world-space (RAS+) orientation. |
| `principal_axis` | Deprecated alias for `focal_axis`; emits a warning. Previously used PCA on transducer positions, which only equals the beam axis for bilateral setups. |

In all modes the returned `scale_rotate_recenter_matrix` encodes the full T1-space → rotated/rescaled-grid transform, which is used downstream to map all transducer positions into the simulation space. 

### Scaling

The 3D planning image is then resliced to the requested grid resolution by isotropic scaling and centered regridding. These steps are applied to the planning image, the segmentation, and the bone mask (or (p)CT, if specified). For masks, nearest neighbor interpolation is used, for continuous values (e.g., from a (pseudo-)CT, linear interpolation is applied).

### Medium Mapping (incl. segmentation smoothing)

Segmentations are mapped onto the requested layer medium masks (to later assign acoustic properties). The standard (10-tissue) SimNIBS 4 segmentation will be converted into the (by default 5) specified layers (controlled by `layers`). If a segmentation for the requested layer is not available, requested layers will be removed from the specification. All segmentation IDs that are not explicitly assigned in dedicated layers are automatically included in a baseline water layer. Medium maps contain layer indices corresponding to the respective layer label position in the medium properties (`medium_properties`). This allows for flexibility in requested layers with stable medium ids (unless `medium_properties` is edited).

Each layer is smoothed using one of the algorithms below. The size of smoothing kernels refer to the grid size of the medium mask, NOT the original dimensions of the planning image. If a multi-layer skull model is requested, both cortical and trabecular layers are initially combined into a skull segment and smoothed. Subsequently, the trabecular bone layer is smoothed and added.

Multiple smoothing algorithms are available:  
- `gaussian` (default) uses `imgaussfilt` / `imgaussfilt3` for isotropic blurring, yielding smooth gradients and moderate edge erosion (default: 1 mm FWHM kernel, converted to voxels at runtime)
- `box` applies a uniform convolution, which is fastest, but also blocky and prone to skull holes at low binarization thresholds (default: 1 mm kernel)
- `off` disables smoothing entirely, passing the raw segmentation mask through without any filtering or binarization

Different binarization thresholds can be defined for skull (default: 50%) and other tissues (default: 50%). 

PRESTUS attempts to enforce skull continuity (i.e., absence of holes) either by applying a 'rubberwrap' skull inflation algorithm (default) or 'imclose' to the skull layer. This is controlled by `headmodel.skull_fill_method`. Potential holes are converted to (cortical) skull. 

### Cropping

To increase computational efficiency, the grid is tightly cropped around the head. The cropping factor is determined by the edges of the layered medium and the location of the transducer(s). [Prior to adding the transducer, the grid can optionally be expanded via `parameters.headmodel.head_pad_mm`, which applies symmetric padding and can avoid that the transducer ends out of grid bounds. This should not impact later computational efficiency of the simulation, as regions beyond the transducer + PML will afterwards be cropped]. To guide the bounds of the head area of interest in the simulation, the CSF segmentation (if available) will be expanded by `parameters.headmodel.csf_expansion` [in mm, converted to voxels at the current grid resolution] to guide bounds of the layered medium. Areas outside of the crop mask will be assigned to the water layer (and are eligible for cropping). Tight crop dimensions are then defined by enclosing the layered medium mask, the transducer geometry, plus a PML buffer. This is followed by FFT-optimized resizing via `find_min_factor`, which expands dimensions up to `parameters.grid.max_expand` voxels for faster FFT convolution. The medium is anchored at the original origin (transducer/target unchanged relative to head), with symmetric expansion/padding around the content to fit FFT requirements while maintaining anatomical alignment for k-Wave propagation.

**Multi-transducer cropping.** When multiple transducers are configured, the crop mask includes the bowl geometry of *every* transducer, not only the first. This is important for bilateral setups (e.g., two transducers on opposite temporal sides targeting the thalamus): without including the contralateral transducer bowl, the front of the grid would be bounded only by the CSF-expanded head tissue, which extends to the far scalp plus `csf_expansion` mm — leaving unnecessary space. With all bowls included, the crop is tight to the outermost water-bath region of each transducer. Grid alignment and rotation are still driven by the first transducer's focal axis; the positions of additional transducers are mapped into the same rotated/rescaled space before the crop bounds are computed.

All image transformations are concatenated to update transducer and target positions. 

### Acoustic property mapping

See [acoustic property mapping](doc_medium.md).

### Smoothing of acoustic property maps

Sharp transitions in the acoustic properties of adjacent media can result in substantial wave reflections. To smoothen tissue transitions, acoustic property maps can be smoothed prior to the simulation. This is governed by `parameters.headmodel.smooth_properties`. Identical smoothing will be applied to all acoustic property maps. This smoothing is deactivated by default, in part to allow users to set the desired smoothing kernel at the intended grid size. If set to `true`, this will apply the same smoothing parameters as used during the creation of the medium maps (`headmodel.smooth_method` and `headmodel.smooth_fwhm_mm`, default: 1 mm FWHM gaussian).

> pCT. If pCTs are generated using PRESTUS, they automatically have smoothing applied at the skull edges ([see the pCT documentation](doc_pseudoCT.md)). If acoustic properties are mapped from pCTs, smoothing of the resulting acoustic property maps may not be intended. Recommended settings are an active area of investigation. 