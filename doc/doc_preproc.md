# Head Preprocessing

!!! note "At a glance"
    - **What:** How the SimNIBS segmentation becomes the simulation grid: rotation, rescaling, smoothing, skull hole filling and cropping.
    - **When you need it:** Rarely: the defaults suit most studies. Read it when the skull looks wrong in the report, or to understand why the grid has the size it has.
    - **Before you start:** A `layered` simulation with a SimNIBS segmentation.

PRESTUS turns a `layered` segmentation into a simulation grid in five steps. Each has one or two settings; the defaults are listed with them.

![preprocessing](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/preprocessing.png)

## Orientation

`grid.mode` sets how the head is rotated before cropping, so that the beam runs along the grid's z axis:

| Mode | What it does |
|---|---|
| `focal_axis` (default) | Rotates the volume so that the mean beam axis of all transducers aligns with z. With one transducer this is its beam axis. |
| `transducer_axis` | Rotates so that the beam axis of transducer(1) aligns with z. Use it when other transducers are secondary. |
| `ras_plus` | No rotation; only isotropic rescaling to the grid resolution. Keeps the world-space (RAS+) orientation. |
| `principal_axis` | Removed; errors. Use `focal_axis`. |

??? note "Details: how the mean beam axis is formed"
    Beams are averaged as axes, not directions (`mean_beam_axis`): a beam pointing against transducer(1)'s is reversed first. Beams converging from the same side give their bisector; opposed (bilateral) beams give the bisector of the two axes, so each beam is tilted from z by about half the angle between their axes (two beams 41° off exactly opposed: 21° each). z points along transducer(1)'s beam. If the mean axis vector is shorter than `grid.focal_axis_min_norm` (default 0.5), transducer(1)'s axis is used instead; for two transducers this cannot happen (minimum 0.71). Implemented in `preproc_align_to_mean_focal_axis` (`focal_axis`) and `preproc_align_to_focal_axis` (`transducer_axis`). In all modes the returned `scale_rotate_recenter_matrix` encodes the full T1-space to grid transform, used to map every transducer and target position into the simulation space.

## Scaling

The planning image, the segmentation and the bone mask (or (pseudo-)CT) are resliced to `grid.resolution_mm` by isotropic scaling and centred regridding. Masks use nearest-neighbour interpolation; continuous images such as a (pseudo-)CT use linear interpolation.

## Layer mapping and smoothing

The SimNIBS segmentation (10 tissues) is mapped onto the layers in `layers` (5 by default). Labels not assigned to a layer become water. A requested layer without a segmentation label is dropped. Each layer index corresponds to its position in `medium_properties`, so layers can be added or removed without changing the ids of the others.

Each layer mask is smoothed and binarised. Settings:

| Setting | Default | Options |
|---|---|---|
| `headmodel.smooth_method` | `gaussian` | `gaussian` (isotropic blur, smooth edges, moderate erosion), `box` (uniform kernel, fastest, blocky, prone to skull holes), `off` (raw mask, no filtering or binarisation) |
| `headmodel.smooth_fwhm_mm` | 1 mm | Kernel size in mm, converted to voxels at the grid resolution (not the T1 resolution) |
| `headmodel.smooth_threshold_skull`, `headmodel.smooth_threshold_other` | 0.5 | Binarisation threshold after smoothing; higher gives a thinner mask |
| `headmodel.skull_fill_method` | `rubberwrap` | `rubberwrap` (skull inflation) or `imclose`; holes become cortical skull |

With a two-layer skull, cortical and trabecular bone are first combined, smoothed and hole-filled as one skull; the trabecular layer is then smoothed and added inside it.

## Cropping

The grid is cropped tightly around the head, every transducer bowl and a PML buffer, which keeps the simulation as small as possible. Two settings change the crop:

- `headmodel.csf_expansion` [mm]: how far beyond the CSF segmentation the head region extends. Voxels outside become water and are eligible for cropping.
- `headmodel.head_pad_mm`: symmetric padding added before the transducer is placed, so that a transducer far from the scalp does not fall outside the grid. Regions beyond the transducer and the PML are cropped again afterwards, so this does not slow the simulation.

??? note "Details: crop bounds and FFT-friendly sizes"
    The crop box encloses the layered medium mask, the bowl geometry of every transducer, and a PML buffer. With several transducers all bowls count, which matters for bilateral set-ups: without the contralateral bowl, that side of the grid would be bounded only by the CSF-expanded head and would carry unused space. Alignment and rotation still follow the first transducer's axis; the other transducers are mapped into the same space before the bounds are computed. The box is then enlarged by up to `grid.max_expand` voxels per dimension to sizes with small prime factors (`find_min_factor`), which makes k-Wave's FFTs faster. The medium stays anchored at its origin, so transducer and target do not move relative to the head. All transformations are concatenated and applied to the transducer and target positions.

## Smoothing of the property maps

Sharp transitions between tissues reflect sound. `headmodel.smooth_properties: true` smooths the acoustic property maps themselves, with the same `smooth_method` and `smooth_fwhm_mm` as above. It is off by default. Pseudo-CTs made with PRESTUS already smooth the skull edges ([pseudo-CT](doc_pseudoCT.md)), so property smoothing is usually not wanted with them. Recommended settings are an active area of investigation.

Which property each tissue gets is described in [Medium properties](doc_medium.md).
