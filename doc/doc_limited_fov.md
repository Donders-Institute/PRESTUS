## Limited acoustic FOV

Running the full k-Wave acoustic simulation on a whole-head grid can be slow and
memory-intensive when only the focal zone and near-field matter.  The acoustic FOV
feature crops the k-Wave domain to a beam-axis-aligned box before simulation and
back-projects the result onto the full grid afterwards.  The thermal simulation
always runs on the full head grid, regardless of any acoustic FOV setting.

---

### Quick start

Add the following keys to your study YAML (or set them on the `parameters` struct
before calling `prestus_pipeline_start`):

```yaml
grid:
  acoustic_fov_diameter_mm: 60   # lateral width [mm], centred on the focus
  acoustic_fov_length_mm:  120   # axial depth [mm], starting at the transducer
```

Choose values that comfortably enclose the focal zone and any near-field region
of interest.  The pipeline prints the resulting sub-grid dimensions and offset at
runtime so you can verify the crop before committing to a full run.

---

### How the FOV box is defined

The FOV is a rectangular box in the beam-axis-aligned acoustic grid:

| Dimension | Extent | Anchor |
|---|---|---|
| Axial (beam axis, z in `transducer_axis` mode) | `acoustic_fov_length_mm` | starts at the transducer centre |
| Lateral (x, y) | `acoustic_fov_diameter_mm` (full width) | centred on the focal point |

The axial start is placed at the transducer centre rather than at the transducer
face or the focus.  This ensures the entire near-field between the transducer and
the focus is captured, which is important for skull heating metrics.  The axial
end should extend at least a few wavelengths beyond the focus to avoid boundary
artefacts from the PML.

The box is clamped to the full acoustic grid.  If `acoustic_fov_length_mm` is too
short for the focus to fall inside, a warning is raised and `acoustic_fov_length_mm`
should be increased.

---

### Pipeline stages affected

| Stage | Effect |
|---|---|
| Medium setup (Stage 4) | Full-grid state is saved before the crop. |
| FOV crop | `acoustic_grid_fov` crops `kwave_medium`, `medium_masks`, `segmentation`, `trans_pos`, `focus_pos`, and `parameters.grid.dims` to the sub-grid. |
| Source/sensor setup (Stage 5) | `kgrid` and the source matrix are built on the sub-grid. |
| Acoustic simulation (Stage 6) | k-Wave runs on the sub-grid; `sensor_data.p_max_all` is sub-grid sized. |
| Acoustic analysis & NIfTI export (Stage 7) | Analysis runs on the sub-grid; outputs are back-projected to the full grid before NIfTI writing so spatial location is correct. |
| Thermal setup (Stage 8) | Full-grid state is fully restored.  `sensor_data.p_max_all` is back-projected to the full grid via `fov_backproject_volume`.  `kgrid` is rebuilt from the full-grid parameters via `thermal_grid_setup`. |
| Thermal simulation & analysis (Stages 8–9) | Runs on the full grid.  Results are physically correct. |

The acoustic cache (`.mat` file) stores the **sub-grid** pressure field together
with the crop metadata (`acoustic_provenance.fov_offset_ac` and
`acoustic_provenance.full_ac_dims`).  This keeps the cache small and lets
`assemble_limited_fov_fields` back-project and combine multiple limited-FOV
results from independent transducer placements.

---

### Why thermal is always full-grid

Thermal conduction has a much longer spatial scale than the acoustic focal zone.
Heat deposited at the focus diffuses into surrounding tissue over the sonication
duration; cropping the thermal domain would:

- suppress heat loss to surrounding tissue, overestimating peak temperature rise;
- miss heating of non-target structures just outside the FOV (e.g. skull surface);
- produce physically incorrect boundary conditions at the crop edges.

For these reasons the acoustic FOV crop is never propagated to the thermal stage.
If faster thermal simulations are needed, use `grid.thermal_resolution_mm` (coarser
grid) or `grid.thermal_fov_mm` (explicit thermal crop, for prototyping only — see
the caveat in [doc_parameters.md](doc_parameters.md)).

---

### Multi-transducer setups

When multiple transducers each target a different focus with independent limited-FOV
acoustic runs, use `assemble_limited_fov_fields` to combine them before thermal:

```matlab
assemble_limited_fov_fields( ...
    {'sub-001_layered_results_t1.mat', 'sub-001_layered_results_t2.mat'}, ...
    'sub-001_layered_results_combined.mat', ...
    [target_t1_wcm2, target_t2_wcm2]);
```

The assembler back-projects each sub-grid onto the shared full acoustic grid,
applies incoherent intensity summation, and writes a combined cache that the
thermal pipeline reads as a normal full-grid result.  The full-head domain is
therefore always used for thermal regardless of how many limited-FOV acoustic
runs contributed.

See [doc_multitransducer.md](doc_multitransducer.md) for the full multi-transducer
workflow.

---

### Choosing FOV dimensions

| Parameter | Too small | Too large |
|---|---|---|
| `acoustic_fov_length_mm` | Focus outside FOV → warning, clamped focal position, incorrect simulation | No harm; wastes grid voxels |
| `acoustic_fov_diameter_mm` | Lateral beam clipped, boundary reflections near focus | No harm; wastes grid voxels |

**Recommended starting values** for a single-element transducer at a ~70 mm depth:

```yaml
acoustic_fov_diameter_mm: 60    # 3–4× the −6 dB beam diameter at focus
acoustic_fov_length_mm:  120    # transducer–focus distance + ~30 mm margin
```

At 2 mm resolution the above box is roughly 30 × 30 × 60 voxels (54 k), compared
to a typical whole-head grid of ~90 × 90 × 200 voxels (1.6 M): a ~97 % reduction
in voxel count.

Print the actual reduction at runtime:

```matlab
fprintf('FOV: [%s]  full: [%s]  reduction: %.1f %%\n', ...
    num2str(parameters.grid.dims), num2str(full_ac_dims), ...
    100*(1 - prod(parameters.grid.dims)/prod(full_ac_dims)));
```

---

### Relevant functions

| Function | Role |
|---|---|
| `acoustic_grid_fov` | Crops medium, masks, and positions to the FOV sub-grid; returns crop offset. |
| `fov_backproject_volume` | Zero-pads a FOV-sized array back to a given full-grid size. |
| `assemble_limited_fov_fields` | Back-projects and incoherently combines N limited-FOV acoustic caches for multi-transducer thermal simulation. |
| `thermal_grid_setup` | Rebuilds `kgrid` and resamples medium/masks for the thermal stage; always operates on the full acoustic grid when `has_acoustic_fov` is true. |

---

### Related documentation

- [doc_parameters.md](doc_parameters.md) — `grid.acoustic_fov_diameter_mm`, `grid.acoustic_fov_length_mm`, `grid.thermal_fov_mm`, `grid.thermal_resolution_mm`
- [doc_multitransducer.md](doc_multitransducer.md) — multi-transducer workflows
- [doc_simulations-acoustic.md](doc_simulations-acoustic.md) — acoustic simulation overview
- [doc_simulations-thermal.md](doc_simulations-thermal.md) — thermal simulation overview
