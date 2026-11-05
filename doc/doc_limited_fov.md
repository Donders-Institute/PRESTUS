# Limited Acoustic FOV

!!! note "At a glance"
    - **What:** Running the acoustic simulation on a smaller box around the beam to save time and memory.
    - **When you need it:** Optional: for large grids, fine resolutions or many runs.
    - **Before you start:** A working head simulation to compare against.

A full k-Wave acoustic simulation on a whole-head grid is slow and memory-intensive when only the focal zone and near field matter. The acoustic FOV crops the acoustic domain (k-Wave, BabelBrain or jwave) to a beam-axis-aligned box before the simulation and back-projects the result onto the full grid afterwards. The thermal simulation always runs on the full head grid.

## Quick start

Add the following keys to your study YAML (or set them on the `parameters` struct before calling `prestus_pipeline_start`):

```yaml
grid:
  acoustic_fov_diameter_mm: 60   # lateral width [mm], centred on the focus
  acoustic_fov_length_mm:  120   # axial depth [mm], starting at the transducer
```

Choose values that enclose the focal zone and any near-field region of interest. The pipeline prints the resulting sub-grid dimensions and offset at runtime, so you can check the crop before a full run.

## How the FOV box is defined

The FOV is a rectangular box in the beam-axis-aligned acoustic grid:

| Dimension | Extent | Anchor |
|---|---|---|
| Axial (beam axis, z in `transducer_axis` mode) | `acoustic_fov_length_mm` | starts at the transducer centre |
| Lateral (x, y) | `acoustic_fov_diameter_mm` (full width) | centred on the focal point |

The axial start sits at the transducer centre rather than at the transducer face or the focus, so the whole near field between transducer and focus is captured; this matters for skull heating metrics. The axial end should extend at least a few wavelengths beyond the focus to avoid boundary artefacts from the PML. The box is clamped to the full acoustic grid. If `acoustic_fov_length_mm` is too short for the focus to fall inside, a warning is raised; increase `acoustic_fov_length_mm`.

## Choosing FOV dimensions

| Parameter | Too small | Too large |
|---|---|---|
| `acoustic_fov_length_mm` | Focus outside FOV: warning, clamped focal position, incorrect simulation | No harm; wastes grid voxels |
| `acoustic_fov_diameter_mm` | Lateral beam clipped, boundary reflections near focus | No harm; wastes grid voxels |

Recommended starting values for a single-element transducer at a ~70 mm depth:

```yaml
acoustic_fov_diameter_mm: 60    # 3–4× the −6 dB beam diameter at focus
acoustic_fov_length_mm:  120    # transducer–focus distance + ~30 mm margin
```

At 2 mm resolution this box is roughly 30 × 30 × 60 voxels (54 k), compared with a typical whole-head grid of ~90 × 90 × 200 voxels (1.6 M): a ~97 % reduction in voxel count. Print the actual reduction at runtime:

```matlab
fprintf('FOV: [%s]  full: [%s]  reduction: %.1f %%\n', ...
    num2str(parameters.grid.dims), num2str(full_ac_dims), ...
    100*(1 - prod(parameters.grid.dims)/prod(full_ac_dims)));
```

## Why thermal is always full-grid

Thermal conduction has a much longer spatial scale than the acoustic focal zone. Heat deposited at the focus diffuses into surrounding tissue over the sonication duration. Cropping the thermal domain would suppress heat loss to surrounding tissue (overestimating the peak temperature rise), miss heating of non-target structures just outside the FOV (such as the skull surface), and give incorrect boundary conditions at the crop edges. The acoustic FOV crop is therefore never propagated to the thermal stage. For faster thermal simulations, use `grid.thermal_resolution_mm` (coarser grid) or `grid.thermal_fov_mm` (explicit thermal crop, for prototyping only; see the caveat in [doc_parameters.md](doc_parameters.md)).

## Multi-transducer setups

The FOV works in both multi-transducer couplings ([details](doc_multitransducer.md#partial-acoustic-fov)):

- **Coherent** (all transducers in one run): the FOV box is the union of one box per transducer, each around that transducer's own beam axis, and every transducer's position is moved into the FOV grid. `'auto'` is resolved per transducer.
- **Async** (`transducer_coupling: async`): each transducer's acoustic run has its own FOV box; `combine_async_intensity` back-projects each field onto its full grid before combining, in `shared` and `per_transducer` grid mode.

Coherent fields can also be built from per-transducer limited-FOV runs by adding their complex fields (`io.save_p_complex`, `superpose_complex_fields`); the sum is complete only where the FOV boxes overlap.

For independent limited-FOV acoustic runs outside async mode, combine them with `assemble_limited_fov_fields` before thermal:

```matlab
assemble_limited_fov_fields( ...
    {'sub-001_layered_results_t1.mat', 'sub-001_layered_results_t2.mat'}, ...
    'sub-001_layered_results_combined.mat', ...
    [target_t1_wcm2, target_t2_wcm2]);
```

The assembler back-projects each sub-grid onto the shared full acoustic grid, sums the intensities incoherently, and writes a combined cache that the thermal pipeline reads as a normal full-grid result. See [doc_multitransducer.md](doc_multitransducer.md) for the full multi-transducer workflow.

??? note "Details: pipeline stages, cache and functions"
    | Stage | Effect |
    |---|---|
    | Medium setup (Stage 4) | Full-grid state is saved before the crop. |
    | FOV crop | `acoustic_grid_fov` crops `kwave_medium`, `medium_masks`, `segmentation`, `trans_pos`, `focus_pos`, and `parameters.grid.dims` to the sub-grid. |
    | Source/sensor setup (Stage 5) | `kgrid` and the source matrix are built on the sub-grid. |
    | Acoustic simulation (Stage 6) | k-Wave runs on the sub-grid; `sensor_data.p_max_all` is sub-grid sized. |
    | Acoustic analysis and NIfTI export (Stage 7) | Analysis runs on the sub-grid; outputs are back-projected to the full grid before NIfTI writing, so the spatial location is correct. |
    | Thermal setup (Stage 8) | Full-grid state is restored. `sensor_data.p_max_all` is back-projected to the full grid via `fov_backproject_volume`. `kgrid` is rebuilt from the full-grid parameters via `thermal_grid_setup`. |
    | Thermal simulation and analysis (Stages 8–9) | Runs on the full grid. |

    The acoustic cache (`.mat` file) stores the **sub-grid** pressure field together with the crop metadata (`acoustic_provenance.fov_offset_ac` and `acoustic_provenance.full_ac_dims`). This keeps the cache small and lets `assemble_limited_fov_fields` back-project and combine multiple limited-FOV results from independent transducer placements.

    | Function | Role |
    |---|---|
    | `acoustic_grid_fov` | Crops medium, masks and positions to the FOV sub-grid; returns the crop offset. |
    | `fov_backproject_volume` | Zero-pads a FOV-sized array back to a given full-grid size. |
    | `assemble_limited_fov_fields` | Back-projects and incoherently combines N limited-FOV acoustic caches for multi-transducer thermal simulation. |
    | `thermal_grid_setup` | Rebuilds `kgrid` and resamples medium/masks for the thermal stage; always operates on the full acoustic grid when `has_acoustic_fov` is true. |

## Related documentation

- [doc_parameters.md](doc_parameters.md): `grid.acoustic_fov_diameter_mm`, `grid.acoustic_fov_length_mm`, `grid.thermal_fov_mm`, `grid.thermal_resolution_mm`
- [doc_multitransducer.md](doc_multitransducer.md): multi-transducer workflows
- [doc_simulations-acoustic.md](doc_simulations-acoustic.md): acoustic simulation overview
- [doc_simulations-thermal.md](doc_simulations-thermal.md): thermal simulation overview
