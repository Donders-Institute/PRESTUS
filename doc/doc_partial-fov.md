# Partial FOV

![Pipeline steps; this page belongs to 4 Grid & medium](assets/figures/fig_stage-strip-medium_diagram.svg){ .stage-strip }

A full acoustic simulation on a whole-head grid is slow and memory-intensive when only the focal zone and near field matter. The acoustic FOV crops the acoustic domain ([k-Wave](https://github.com/ucl-bug/k-wave), [BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain) or [jwave](https://github.com/ucl-bug/jwave)) to a beam-axis-aligned box before the simulation and back-projects the result onto the full grid afterwards. A box that misses part of the beam or the focus gives a wrong field, so its size has to enclose both.

On this page:

- [Quick Start](#quick-start): the two keys that switch the FOV on.
- [FOV Box Geometry](#fov-box-geometry): where the box starts and how `auto` sizes it.
- [Choosing FOV Dimensions](#choosing-fov-dimensions): what too small and too large boxes do, with starting values.
- [Thermal Grid](#thermal-grid): why the thermal simulation keeps the full grid by default.
- [Multi-Transducer Setups](#multi-transducer-setups): FOV boxes for coherent, async, and independent runs.

!!! note "At a glance"
    - **When you need it:** Optional: for large grids, fine resolutions, or many runs.
    - **Before you start:** A working head simulation to compare against.

## Quick Start

Add the following keys to your study [YAML](https://yaml.org/) (or set them on the `parameters` struct before calling `prestus_pipeline_start`):

```yaml
grid:
  acoustic_fov_diameter_mm: 60   # lateral width [mm], centred on the focus
  acoustic_fov_length_mm:  120   # axial depth [mm], starting at the transducer
```

Choose values that enclose the focal zone and any near-field region of interest. The pipeline prints the resulting sub-grid dimensions and offset at runtime, so you can check the crop before a full run.

## FOV Box Geometry

The FOV is a rectangular box in the beam-axis-aligned acoustic grid:

| Dimension | Extent | Anchor |
|---|---|---|
| Axial (beam axis, z in `transducer_axis` mode) | `acoustic_fov_length_mm` | starts 8 voxels behind the bowl apex (`trans_pos`) |
| Lateral (x, y) | `acoustic_fov_diameter_mm` (full width) | centred on the focal point |

The axial start sits 8 voxels behind the bowl apex, which keeps the source points spread behind it, rather than at the focus, so the whole near field between transducer and focus is captured; this matters for skull heating metrics. The axial end should extend at least a few wavelengths beyond the focus to avoid boundary artefacts from the [PML](doc_glossary.md#pml). `acoustic_fov_length_mm: 'auto'` sets the length to the apex-to-focus distance (`trans_pos` to `focus_pos`) plus `grid.z_beyond_focus_mm` (default 40 mm); `acoustic_fov_diameter_mm: 'auto'` uses the outer element diameter of `transducer(1)` plus a margin of 4 grid voxels per side, which holds the source points spread beyond the bowl rim. A numeric diameter smaller than that would crop the source; it is enlarged with a warning. The box is clamped to the full acoustic grid. If `acoustic_fov_length_mm` is too short for the focus to fall inside, a warning is raised; increase `acoustic_fov_length_mm`.

## Choosing FOV Dimensions

| Parameter | Too small | Too large |
|---|---|---|
| `acoustic_fov_length_mm` | Focus outside FOV: warning, clamped focal position, incorrect simulation | No harm; wastes grid voxels |
| `acoustic_fov_diameter_mm` | Lateral beam clipped, boundary reflections near focus | No harm; wastes grid voxels |

Recommended starting values for a single-element transducer at a ~70 mm depth:

```yaml
acoustic_fov_diameter_mm: 60    # 3–4× the −6 dB beam diameter at focus
acoustic_fov_length_mm:  120    # bowl apex–focus distance + 40 mm (grid.z_beyond_focus_mm default)
```

At 2 mm resolution this box is roughly 30 × 30 × 60 voxels (54 k), compared with a whole-head grid of 256 × 216 × 192 mm (128 × 108 × 96 voxels, 1.3 M): a ~96 % reduction in voxel count. Print the actual reduction at runtime:

```matlab
fprintf('FOV: [%s]  full: [%s]  reduction: %.1f %%\n', ...
    num2str(parameters.grid.dims), num2str(full_ac_dims), ...
    100*(1 - prod(parameters.grid.dims)/prod(full_ac_dims)));
```

## Thermal Grid

Thermal conduction has a much longer spatial scale than the acoustic focal zone. Heat deposited at the focus diffuses into surrounding tissue over the sonication duration. Cropping the thermal domain would suppress heat loss to surrounding tissue (overestimating the peak temperature rise), miss heating of non-target structures just outside the FOV (such as the skull surface), and give incorrect boundary conditions at the crop edges. By default (`grid.thermal_use_acoustic_fov: 0`) the thermal stage therefore restores the full grid and back-projects the acoustic field onto it. With `grid.thermal_use_acoustic_fov: 1`, the thermal simulation runs inside the acoustic FOV box instead; this is faster, but the cooling boundary is less accurate. For faster thermal simulations, use `grid.thermal_resolution_mm` (coarser grid) or `grid.thermal_fov_mm` (explicit thermal crop, for prototyping only; see the caveat in [PRESTUS Parameters](doc_parameters.md)).

## Multi-Transducer Setups

The FOV works in both multi-transducer couplings ([Partial Acoustic FOV](doc_multitransducer.md#partial-acoustic-fov)):

- **Coherent** (all transducers in one run): the FOV box is the union of one box per transducer, each around that transducer's own beam axis, and every transducer's position is moved into the FOV grid. `'auto'` is resolved per transducer.
- **Async** (`transducer_coupling: async`): each transducer's acoustic run has its own FOV box; `combine_async_intensity` back-projects each field onto its full grid before combining, in `shared` and `per_transducer` grid mode.

Coherent fields can also be built from per-transducer partial-FOV runs by adding their complex fields (`io.save_p_complex`, `superpose_complex_fields`); the sum is complete only where the FOV boxes overlap.

For independent partial-FOV acoustic runs outside async mode, combine them with `assemble_limited_fov_fields` before thermal:

```matlab
assemble_limited_fov_fields( ...
    {'sub-001_layered_results_t1.mat', 'sub-001_layered_results_t2.mat'}, ...
    'sub-001_layered_results_combined.mat', ...
    [target_t1_wcm2, target_t2_wcm2]);
```

The assembler back-projects each sub-grid onto the shared full acoustic grid, sums the intensities incoherently, and writes a combined cache that the thermal pipeline reads as a normal full-grid result. See [Multi-Transducer Modelling](doc_multitransducer.md) for the full multi-transducer workflow.

??? note "Pipeline stages, cache, and functions"
    | Stage | Effect |
    |---|---|
    | Medium setup | Full-grid state is saved before the crop. |
    | FOV crop | `acoustic_grid_fov` crops `kwave_medium`, `medium_masks`, `segmentation`, `trans_pos`, `focus_pos`, and `parameters.grid.dims` to the sub-grid. |
    | Source/sensor setup | `kgrid` and the source matrix are built on the sub-grid. |
    | Acoustic simulation | k-Wave runs on the sub-grid; `sensor_data.p_max_all` is sub-grid sized. |
    | Acoustic analysis and NIfTI export | Analysis runs on the sub-grid; outputs are back-projected to the full grid before [NIfTI](https://nifti.nimh.nih.gov/) writing, so the spatial location is correct. |
    | Thermal setup | Full-grid state is restored. `sensor_data.p_max_all` is back-projected to the full grid via `fov_backproject_volume`. `kgrid` is rebuilt from the full-grid parameters via `thermal_grid_setup`. With `grid.thermal_use_acoustic_fov: 1`, the FOV sub-grid is kept and only the T1 ↔ grid affines are shifted by the crop offset. |
    | Thermal simulation and analysis (Stages 8–9) | Runs on the full grid, or on the FOV sub-grid with `grid.thermal_use_acoustic_fov: 1`. |

    The acoustic cache (`.mat` file) stores the **sub-grid** pressure field together with the crop metadata (`acoustic_provenance.fov_offset_ac` and `acoustic_provenance.full_ac_dims`). This keeps the cache small and lets `assemble_limited_fov_fields` back-project and combine multiple partial-FOV results from independent transducer placements.

    | Function | Role |
    |---|---|
    | `acoustic_grid_fov` | Crops medium, masks, and positions to the FOV sub-grid; returns the crop offset. |
    | `fov_backproject_volume` | Zero-pads a FOV-sized array back to a given full-grid size. |
    | `assemble_limited_fov_fields` | Back-projects and incoherently combines N partial-FOV acoustic caches for multi-transducer thermal simulation. |
    | `thermal_grid_setup` | Rebuilds `kgrid` and resamples medium/masks for the thermal stage, on the full acoustic grid or, with `grid.thermal_use_acoustic_fov: 1`, on the FOV sub-grid. |

## Related Documentation

- [PRESTUS Parameters](doc_parameters.md): `grid.acoustic_fov_diameter_mm`, `grid.acoustic_fov_length_mm`, `grid.thermal_use_acoustic_fov`, `grid.thermal_fov_mm`, `grid.thermal_resolution_mm`
- [Multi-Transducer Modelling](doc_multitransducer.md): multi-transducer workflows
- [Uncertainty: Partial-FOV Bound Variants](doc_uncertainty.md#partial-fov-bound-variants): bound variants in the partial FOV, rescaled to the whole grid
- [Acoustic Simulations](doc_simulations-acoustic.md): acoustic simulation overview
- [Thermal Simulations](doc_simulations-thermal.md): thermal simulation overview
