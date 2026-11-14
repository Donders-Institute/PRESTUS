# Segmentation and Head Preprocessing

The skull reflects, attenuates and refracts ultrasound, so how the head is represented in the simulation shapes the predicted focus and intensity. This step turns an MRI into that representation: [SimNIBS](https://simnibs.github.io/simnibs/build/html/index.html) labels the tissues, and PRESTUS orients, resamples, smooths and crops the labels into the simulation grid.

On this page:

- [Segmentation](#segmentation): running SimNIBS `charm` and reusing its result.
- [From Segmentation to Grid](#from-segmentation-to-grid): orientation, scaling, layer mapping, cropping and smoothing of the property maps.

!!! note "At a glance"
    - **When you need it:** Rarely: the defaults suit most studies. Read it when the skull looks wrong in the report, or to understand why the grid has the size it has.
    - **Before you start:** A `layered` simulation with a SimNIBS segmentation.

## Segmentation

For a `layered` medium, the first run segments the T1 (and the T2, when `path.t2_pattern` is set) with SimNIBS `charm` into `<path.seg>/m2m_sub-XXX/`. This takes an hour or more. Later runs reuse `final_tissues.nii.gz` from that folder unless `io.overwrite_simnibs: 1`.

- `modules.segmentation_only: 1` segments and stops, e.g. to segment all subjects before simulating them.
- `segmentation.use_qform: 1` when `charm` reports a qform/sform mismatch; `segmentation.debug: 1` passes `--debug` to `charm`.
- `segmentation.noneck: 1` passes `--noneck` for a scan without a neck.

### SimNIBS Version and charm Settings

`startup.simnibs_bin_path` selects the SimNIBS install, and with it the version; PRESTUS needs SimNIBS 4 or later. While they are empty, `hpc.ld_library_path` (the install's libraries for SimNIBS tools) and `placement.plantus.env_path` (the Python PlanTUS runs in) follow it; a value set there overrides it. To change charm's settings, copy the install's `charm.ini` (in its `simnibs` package folder), edit the copy and set `segmentation.settings_file` to it; keep it outside the m2m folder. Report the settings you changed.

An existing segmentation is reused even when another version or other settings are selected. PRESTUS then warns, comparing the version in `charm_log.html` and the settings in `settings.ini` of the m2m folder with the selection. Set `io.overwrite_simnibs: 1` to segment again, which deletes the m2m folder first, or point `path.seg` to a separate folder per version to keep both. What was built from the old segmentation is detected on reuse: grid, medium, acoustic and head-model caches, placement runs and pseudo-CTs record the segmentation's checksum ([Cache reuse checks](doc_outputs_reuse.md#cache-reuse-checks)), and a defaced copy is made again.

In the Nexus, the Segmentation card (Medium tab) lists the SimNIBS installs on this computer, shows the version and the existing segmentation, opens `charm_report.html`, and has the re-run switch. When runs are submitted to a cluster over SSH (Run tab), *Scan cluster* lists the installs there: conda environments, `~/SimNIBS*`, `/opt/simnibs*` and the bin folders of the `simnibs` environment modules. The Nexus keeps the result per host, so it can check a cluster path that does not exist on this computer.


Installing SimNIBS: [Installation](doc_installation.md#install-simnibs). Errors: [Troubleshooting](doc_troubleshooting.md#segmentation-simnibs-charm).

## From Segmentation to Grid

PRESTUS turns a `layered` segmentation into a simulation grid in three steps: a single resample that rotates and rescales each image ([Orientation](#orientation) and [Scaling](#scaling)), layer mapping with mask smoothing and skull hole filling, and cropping with padding. Optional smoothing of the property maps runs later, when the medium is set up. Each step has one or two settings; the defaults are listed with them.

![preprocessing](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/preprocessing.png)

### Orientation

`grid.mode` sets how the head is rotated before cropping, so that the beam runs along the grid's z axis:

| Mode | What it does |
|---|---|
| `focal_axis` (default) | Rotates the volume so that the mean beam axis of all transducers aligns with z. With one transducer this is its beam axis. |
| `transducer_axis` | Rotates so that the beam axis of transducer(1) aligns with z. Use it when other transducers are secondary. |
| `ras_plus` | No rotation; only isotropic rescaling to the grid resolution. Keeps the world-space (RAS+) orientation. |
| `principal_axis` | Removed; errors. Use `focal_axis`. |

??? note "How the mean beam axis is formed"
    Beams are averaged as axes, not directions (`mean_beam_axis`): a beam pointing against transducer(1)'s is reversed first. Beams converging from the same side give their bisector; opposed (bilateral) beams give the bisector of the two axes, so each beam is tilted from z by about half the angle between their axes (two beams 41° off exactly opposed: 21° each). z points along transducer(1)'s beam. If the mean axis vector is shorter than `grid.focal_axis_min_norm` (default 0.5), transducer(1)'s axis is used instead; for two transducers this cannot happen (minimum 0.71). Implemented in `preproc_align_to_mean_focal_axis` (`focal_axis`) and `preproc_align_to_focal_axis` (`transducer_axis`). In all modes the returned `scale_rotate_recenter_matrix` encodes the full T1-space to grid transform, used to map every transducer and target position into the simulation space.

### Scaling

The planning image, the segmentation and the bone mask (or (pseudo-)CT) are resliced to `grid.resolution_mm` by isotropic scaling and centred regridding, in the same resample as the rotation. Images with fewer than 20 distinct values, such as masks and segmentations, use nearest-neighbour interpolation; all others, such as the planning image or a (pseudo-)CT, use linear interpolation.

### Layer Mapping and Smoothing

The SimNIBS segmentation (10 tissues) is mapped onto the layers in `layers`. The default lists six: `water`, `brain`, `skin`, `skull` (label 4), `skull_cortical` and `skull_trabecular`. Labels not assigned to another layer become water, which by default includes CSF (label 3), eyes, blood and muscle. A requested layer without a segmentation label is dropped. A charm segmentation has no label 4, so the default `skull` layer is dropped and the head consists of brain, skin, cortical and trabecular bone in water; with a pseudo-CT the skull layers are merged into one `skull` layer. Each layer index corresponds to its position in `medium_properties`, so layers can be added or removed without changing the ids of the others.

Each layer mask is smoothed and binarised. Settings:

| Setting | Default | Options |
|---|---|---|
| `headmodel.smooth_method` | `gaussian` | `gaussian` (isotropic blur, smooth edges, moderate erosion), `box` (uniform kernel, fastest, blocky, prone to skull holes), `off` (raw mask, no filtering or binarisation) |
| `headmodel.smooth_fwhm_mm` | 1 mm | Kernel size in mm, converted to voxels at the grid resolution (not the T1 resolution) |
| `headmodel.smooth_threshold_skull`, `headmodel.smooth_threshold_other` | 0.5 | Binarisation threshold after smoothing; higher gives a thinner mask |
| `headmodel.skull_fill_method` | `rubberwrap` | `rubberwrap` (skull inflation) or `imclose`; holes become cortical skull |

With a two-layer skull, cortical and trabecular bone are first combined and smoothed as one skull, which becomes cortical bone; the trabecular layer is then smoothed and added inside it. Holes in the combined skull are filled afterwards and become cortical bone.

### Cropping

The grid is cropped tightly around the head and every transducer bowl, with a 1-voxel margin; [k-Wave](https://github.com/ucl-bug/k-wave) adds the PML outside the grid. This keeps the simulation as small as possible. Two settings change the crop:

- `headmodel.csf_expansion` [mm]: how far beyond the CSF segmentation the head region extends. Voxels outside become water and are eligible for cropping.
- `headmodel.head_pad_mm`: symmetric padding added to the rotated volume before the CSF expansion and the crop. The crop box already contains every bowl, and the volume is padded automatically (with a warning) if a bowl centre falls outside it, so this setting is rarely needed. Padding beyond the head and bowls is cropped away, so it does not slow the simulation.

??? note "Crop bounds and FFT-friendly sizes"
    The crop box encloses the layered medium mask and the bowl geometry of every transducer, plus a 1-voxel margin; the PML lies outside the grid. With several transducers all bowls count, which matters for bilateral set-ups: without the contralateral bowl, that side of the grid would be bounded only by the CSF-expanded head and would carry unused space. Rotation follows the mean beam axis in the default `focal_axis` mode and transducer(1)'s axis in `transducer_axis` mode; all transducers are mapped into the same space before the bounds are computed. The box is then enlarged by up to `grid.max_expand` voxels per dimension to sizes with small prime factors (`find_min_factor`), which makes k-Wave's FFTs faster. The medium stays anchored at its origin, so transducer and target do not move relative to the head. All transformations are concatenated and applied to the transducer and target positions.

### Property Map Smoothing

Sharp transitions between tissues reflect sound. `headmodel.smooth_properties: true` smooths the acoustic and thermal property maps themselves during medium setup, with the same `smooth_method` and `smooth_fwhm_mm` as above. It is off by default. Pseudo-CTs made with PRESTUS already smooth the skull edges ([pseudo-CT](doc_pseudoCT.md)), so property smoothing is usually not wanted with them. Recommended settings are an active area of investigation.

Which property each tissue gets is described in [Medium properties](doc_medium.md).
