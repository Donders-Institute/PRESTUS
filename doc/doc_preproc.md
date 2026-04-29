# Head preprocessing

PRESTUS preprocesses `layered` segmentations to set up simulation grids.

![preprocessing](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/preprocessing.png)

### Orientation

PRESTUS initially creates a grid by orienting the segmentation to the focal plane of the first transducer. In multi-transducer setups, alignment is done to the first specified transducer. Specifically, `preproc_align_to_focal_axis` performs a two-step rotation to align an arbitrary focal axis (defined by transducer and focus positions) with the z-axis of the coordinate system. This prepares subject-specific head models for ultrasound neuromodulation simulations by standardizing the acoustic propagation direction along z. 

### Scaling

The 3D planning image is then resliced to the requested grid resolution by isotropic scaling and centered regridding. These steps are applied to the planning image, the segmentation, and the bone mask (or (p)CT, if specified). For masks, nearest neighbor interpolation is used, for continuous values (e.g., from a (pseudo-)CT, linear interpolation is applied).

### Medium Mapping (incl. segmentation smoothing)

Segmentations are mapped onto the requested layer medium masks (to later assign acoustic properties). The standard (10-tissue) SimNIBS 4 segmentation will be converted into the (by default 5) specified layers (controlled by `parameters.layers`). If a segmentation for the requested layer is not available, requested layers will be removed from the specification. All segmentation IDs that are not explicitly assigned in dedicated layers are automatically included in a baseline water layer. Medium maps contain layer indices corresponding to the respective layer label position in the medium properties (`parameters.medium`). This allows for flexibility in requested layers with stable medium ids (unless `parameters.medium` is edited).

Each layer is smoothed using one of the algorithms below. The size of smoothing kernels refer to the grid size of the medium mask, NOT the original dimensions of the planning image. If a multi-layer skull model is requested, both cortical and trabecular layers are initially combined into a skull segment and smoothed. Subsequently, the trabecular bone layer is smoothed and added.

Multiple smoothing algorithms are available:  
- `gaussian` (default) uses smooth3 for isotropic blurring, yielding smooth gradients and moderate edge erosion (default: 2-voxel FWHM kernel)
- `box` applies a uniform convolution, which is fastest, but also blocky and prone to skull holes at low binarization thresholds (default: 2-voxel kernel)

Different binarization thresholds can be defined for skull (default: 50%) and other tissues (default: 50%). 

PRESTUS attempts to enforce skull continuity (i.e., absence of holes) either by applying a 'rubberwrap' skull inflation algorithm (default) or 'imclose' to the skull layer. This is controlled by `headmodel.skull_fill_method`. Potential holes are converted to (cortical) skull. 

### Cropping

To increase computational efficiency, the grid is tightly cropped around the head. The cropping factor is determined by the edges of the layered medium and the location of the transducer. [Prior to adding the transducer, the grid can optionally be expanded via `pad_mm`, which applies symmetric padding and can avoid that the transducer ends out of grid bounds. This should not impact later computational efficiency of the simulation, as regions beyond the transducer + PML will afterwards be cropped]. To guide the bounds of the head area of interest in the simulation, the CSF segmentation (if available) will be expanded by `csf_mask_expansion_factor` [in voxels] to guide bounds of the layered medium. Areas outside of the crop mask will be assigned to the water layer (and are elibible for cropping). Tight crop dimensions are then defined by enclosing the layered medium mask, the transducer geometry, plus a PML buffer. This is followed by FFT-optimized resizing via `find_min_factor`, which expands dimensions up to  `prime_factor_max_grid_expansion` for faster FFT convolution. The medium is anchored at the original origin (transducer/target unchanged relative to head), with symmetric expansion/padding around the content to fit FFT requirements while maintaining anatomical alignment for k-Wave propagation.

All image transformations are concatenated to update transducer and target positions. 

### Acoustic property mapping

See [acoustic property mapping](doc_medium.md).

### Smoothing of acoustic property maps

Sharp transitions in the acoustic properties of adjacent media can result in substantial wave reflections. To smoothen tissue transitions, acoustic property maps can be smoothed prior to the simulation. This is governed by `smooth_properties`. Identical smoothing will be applied to all acoustic property maps. This smoothing is deactivated by default, in part to allow users to set the desired smoothing kernel at the intended  grid size. If set to `true`, this will apply the same smoothing parameters as used during the creation of the medium maps (by default 2 voxel FWHM).

> pCT. If pCTs are generated using PRESTUS, they automatically have smoothing applied at the skull edges ([see the pCT documentation](doc_pseudoCT.md)). If acoustic properties are mapped from pCTs, smoothing of the resulting acoustic property maps may not be intended. Recommended settings are an active area of investigation. 