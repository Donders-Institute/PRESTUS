# CT-Based Skull Mapping

By default each skull layer of a `layered` simulation takes one value per property, although skull bone properties can vary within and between individuals. A CT gives the skull subject-specific properties instead. With `pct.enabled: 1`, PRESTUS maps the Hounsfield units (HU) of a CT to density, sound speed and attenuation within the skull layer. The CT is one of two images:

- a **measured CT** of the subject (`path.ct_pattern`), or
- a **pseudo-CT** that PRESTUS generates from a UTE/PETRA MR image, used when `path.ct_pattern` is empty.

Both feed the same HU-to-property mapping; only the pseudo-CT needs the UTE-to-HU step first.

On this page:

- [Setting Up a Pseudo-CT Run](#setting-up-a-pseudo-ct-run): the scan, the tools, the minimum config and the outputs.
- [Using a Measured CT](#using-a-measured-ct): the CT's requirements and the minimum config.

!!! note "At a glance"
    - **When you need it:** When you want skull properties that follow the individual skull, and have a CT or a UTE/PETRA scan.
    - **Before you start:** A CT co-registered to the T1w, or a UTE/PETRA image per subject plus [FSL](https://fsl.fmrib.ox.ac.uk/fsl/docs/) and [ANTs](https://github.com/ANTsX/ANTs).

Subpages:

- [From UTE to Pseudo-HU](doc_pseudoCT_generation.md): for choosing `pct.skull_mapping` and the matching `pct.density_calibration`; read it whenever you change either from the default.
- [Skull Property Mappings](doc_pseudoCT_property-mappings.md): for choosing how pHU become skull density, sound speed and attenuation, with the equations and references.

![pct_ute_density](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/ute_pct_density.png)

!!! warning "Experimental"
    (pseudo-)HU mapping is an experimental feature in active development. The most suitable mapping remains an active area of research; treat all mappings as experimental. See [this issue](https://github.com/Donders-Institute/PRESTUS/issues/43).

## Setting Up a Pseudo-CT Run

1. Acquire an MR sequence with detailed bone contrast (PETRA UTE; example parameters below).
2. Run [SimNIBS](https://simnibs.github.io/simnibs/build/html/index.html) with the T1w and the PETRA UTE scan (instead of the T2w) as inputs, either in PRESTUS (**default**) or in the SimNIBS GUI. The UTE image is the one in `path.t2_pattern`, which SimNIBS registers to the T1 as `m2m_sub-XXX/T2_reg.nii.gz` (`UTE_reg.nii.gz` is also read).
3. Set `pct.enabled: 1`. PRESTUS runs `pct_create_pseudoCT.sh` after segmentation if `pseudoCT.nii.gz` does not exist yet. This needs **FSL** (`fslmaths`), **ANTs** (`N4BiasFieldCorrection`, `SmoothImage`) and **[MATLAB](https://www.mathworks.com/products/matlab.html)** on `PATH`; set `startup.fsl_bin_path` and `startup.ants_bin_path` so PRESTUS can prepend them. On HPC, also set `startup.matlab_bin_path`.

Minimum config:

```yaml
path:
  t2_pattern: 'sub-%1$03d_UTE.nii.gz'  # UTE image (used for both SimNIBS and pseudoCT)

startup:
  fsl_bin_path: '/opt/fsl/bin'           # prepended to PATH for fslmaths
  ants_bin_path: '/opt/ants/bin'          # prepended to PATH for N4BiasFieldCorrection, SmoothImage
  matlab_bin_path: '/opt/matlab/R2024a/bin/matlab'  # required on HPC; local falls back to matlabroot

pct:
  enabled: 1
  skull_mapping: 'kosciessa' # UTE→HU linear mapping algorithm
  debug: 1                   # keep intermediate images in the debug_pct/ subfolder
  mapping_density: 'k-plan'
  density_calibration: ''    # HU->density table for k-plan; '' = k-Plan default (see Matching the CT Calibration)
  mapping_soundspeed: 'k-plan'
  mapping_attenuation: 'k-plan'
```

Outputs go to `<path.pct>/sub-XXX/`, or to `<parent of path.seg>/pct/sub-XXX/` when `path.pct` is empty:

| File / folder | Content |
|---|---|
| `pseudoCT.nii.gz`, `tissues_mask.nii.gz` | pseudo-CT and tissue mask |
| `pseudoCT.json` | sidecar: skull mapping, segmentation folder, PRESTUS version |
| `debug_pct/` | intermediate images and the log |
| `log_hpc/` | cluster job logs |

**Reuse.** An existing `pseudoCT.nii.gz` is used and not regenerated. Its `pseudoCT.json` must record the current `pct.skull_mapping`; otherwise the run stops with `prestus:cache:signatureMismatch` (see `io.cache_signature_check`). To regenerate, move the folder aside or set `path.pct` to another folder. A pseudo-CT without a sidecar is used with the warning `prestus:cache:unverified`. A pseudo-CT inside `m2m_sub-XXX/` is used with the warning `PRESTUS:pct:legacyLocation`; move it to the pseudo-CT folder. See [Starting from earlier stage outputs](doc_outputs_reuse.md#starting-from-earlier-stage-outputs).

**Manual.** You can also run `pct_create_pseudoCT.sh` in bash before starting PRESTUS. An example wrapper is at `examples/createPseudoCT.sh`. Required modules: SimNIBS, FSL, ANTs.

??? note "Example PETRA UTE parameters"
    From Carpino et al. (2024):

    - TR: 3.32 ms
    - TE: 0.07 ms
    - voxel size: 0.8 mm3
    - 352 sagittal slices
    - flip angle: 2°
    - FoV: 294 mm

## Using a Measured CT

Set `path.ct_pattern` to the subject's CT and `pct.enabled: 1`. The run then needs no UTE image and generates no pseudo-CT.

```yaml
path:
  ct_pattern: 'sub-%1$03d_CT.nii.gz'  # relative to path.anat, like t1_pattern

pct:
  enabled: 1
  density_calibration: ''  # this CT's HU->density table; '' = k-Plan default
  mapping_density: 'k-plan'
  mapping_soundspeed: 'k-plan'
  mapping_attenuation: 'k-plan'
```

- **Units.** The CT must hold Hounsfield units (NIfTI scaling is applied).
- **Registration.** The CT must be co-registered to the T1w: PRESTUS resamples it onto the segmentation grid through both images' world coordinates and does not register it. Voxels outside the CT are set to air (-1000 HU). The run logs the median HU of the skull voxels and warns (`PRESTUS:pct:ctMisregistered`) when it is below 300 HU.
- **Skull voxels.** The skull layer comes from the SimNIBS segmentation, built as for the pseudo-CT. No files are written to the pseudo-CT folder.
- **Calibration.** HU depend on scanner, kVp and kernel. Set `pct.density_calibration` to this CT's HU-to-density table when you have one; with the k-Plan default, PRESTUS warns (`prestus:pct:ctDefaultCalibration`).
- **Reuse.** Grid and medium caches record the CT's name, size and date, so replacing the CT rebuilds them.
