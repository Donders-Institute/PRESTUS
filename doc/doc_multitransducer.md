# Multi-Transducer Modelling

When several transducers sonicate the same head, the exposure depends on how their fields combine, which a run per transducer does not capture. PRESTUS simulates them in one layered setup: transducers that fire together (`coherent`) add as complex pressures and interfere, while for transducers that fire in separate windows (`async`) acoustic safety metrics take the per-voxel maximum of the fields and heating uses the sum of their intensities.

On this page:

- [Configuration](#configuration): listing transducers and choosing how they fire.
- [Grid setup](#grid-setup): how the grid is aligned and cropped around several beams.
- [ISPPA targets](#isppa-targets): per-transducer free-water targets and their baselines.
- [Complex pressure field](#complex-pressure-field): saving fields and superposing them afterwards at any relative phase.
- [Partial acoustic FOV](#partial-acoustic-fov): limited fields of view with several transducers.
- [Known limitations](#known-limitations): what refers to the first transducer only, and unsupported setups.
- [Subpages](#subpages): async firing in separate windows.

!!! note "At a glance"
    - **When you need it:** Bilateral or dual-site stimulation, or any protocol with more than one transducer on the head.
    - **Before you start:** A layered (skull) simulation in a 3-D grid; one equipment entry and placement per transducer.

!!! note "Experimental support"
    See [this pull request](https://github.com/Donders-Institute/PRESTUS/pull/100) for examples.

## Configuration

List several transducers in a layered simulation:

```yaml
transducer:
  - name: right
    serial: IS_PCD15473_01001
    combo: {ds_serial: IGT_32_ch_comb_10_ch}
    # ... trans_pos, focus_pos, focal_distance_ep

  - name: left
    serial: IS_PCD15473_01001
    combo: {ds_serial: IGT_32_ch_comb_10_ch}
    # ...
```

`simulation.transducer_coupling` sets how the transducers fire:

| Firing | Setting |
|---|---|
| Simultaneously, same frequency, known phase | `coherent` (default): one [k-Wave](https://github.com/ucl-bug/k-wave) run with all sources active; the fields interfere. All transducers must share the fundamental frequency (a warning is issued otherwise). |
| In separate, non-overlapping windows | `async`: see [Async Multi-Transducer Mode](doc_multitransducer_async.md). |
| In sequence, or with partial overlap | Use [Sequential Simulations](doc_sequential.md): carry each heating run's temperature and CEM43 maps into the next via `adopted_heatmap` / `adopted_cem43`. This assumes the acoustic fields do not overlap in time. |

In the Nexus, *Multi-transducer coupling* selects the row above and *Target free-water Isppa* is set per transducer on the transducer page.

## Grid Setup

### Grid Alignment

The default `focal_axis` mode rotates the grid so that the mean beam axis of all transducers aligns with z. Beams are averaged as axes, not directions, so opposed (bilateral) beams do not cancel: z bisects their two axes. Each beam off z is drawn obliquely on the grid, which costs accuracy at coarse resolution (for async runs, see [`async_grid_mode: per_transducer`](doc_multitransducer_async.md#grid-geometry-mode)).

```yaml
grid:
  mode: 'focal_axis'   # default; omit to get the same behaviour
```

To align the grid to transducer(1)'s focal axis instead, set `transducer_axis`. See [Preprocessing: Orientation](doc_preproc.md#orientation) for the mode table.

### Grid Cropping

The grid is cropped to the bounding box that encloses the CSF-expanded head tissue **and the bowl geometry of every transducer**. See [Preprocessing: Cropping](doc_preproc.md#cropping).

## ISPPA Targets

`transducer(k).target_isppa_wcm2` is the free-water ISPPA of **transducer k alone**, for every coupling mode. Each transducer has its own drive (`elem_p_amp`), so the target rescales that transducer's drive, not the combined field.

```yaml
simulation:
  transducer_coupling: coherent
transducer:
  - name: left
    serial: IS_PCD15473_01001
    target_isppa_wcm2: 20     # left transducer alone: 20 W/cm² in free water
  - name: right
    serial: IS_PCD15473_01003
    target_isppa_wcm2: 40     # right transducer alone: 40 W/cm²
```

| Coupling | How the targets are applied | Sweeps (vector targets) |
|---|---|---|
| `coherent` | Before the simulation (`coherent_drive_calibration`): each transducer's `elem_p_amp` is multiplied by √(target / baseline), where baseline is the free-water ISPPA of that transducer alone at its configured drive. The combined field is simulated at these drives and not rescaled afterwards. | Not available. Errors (`prestus:isppa:coherentSweep`); run one config per point, or use `async`. |
| `async` | After the per-transducer simulations: each field is scaled to its own target in `combine_async_intensity`, then combined. | Yes: one acoustic run per transducer, one combine + thermal pair per sweep point (see [Async Multi-Transducer Mode](doc_multitransducer_async.md)). |

**Coherent baseline.** The rule is the same as for one transducer ([ISPPA targeting](doc_simulations-acoustic.md#isppa-targeting)). If the target matches the transducer's calibration ISPPA (`tx_calibration.isppa_wcm2`) within 1 %, the calibration is the baseline. Otherwise PRESTUS simulates that transducer alone in water (others at zero amplitude, same grid). With N transducers that is up to N extra free-water simulations. A transducer with `target_isppa_wcm2: []` keeps its configured drive.

??? note "Provenance and cache checks"
    The result records targets, baselines, their source and the drive scale in `acoustic_provenance.coherent_drive`. The combined free-water ISPPA (both firing) is reported but is not a target. A coherent cache without `coherent_drive` is refused when targets are set (`prestus:isppa:coherentWithoutDriveCalibration`): re-run its acoustic simulation. `path.acoustic_input` cannot be rescaled per transducer and errors when targets are set.

## Complex Pressure Field

```yaml
io:
  save_p_complex: true
```

This writes `sub-XXX_<medium>_p_complex<affix>.mat` with the complex amplitude at the fundamental on the simulation grid, at the simulated drive level, before ISPPA rescaling. Propagation in PRESTUS is linear, so fields simulated one at a time can be summed. `superpose_complex_fields` does this on one shared grid, back-projecting limited-FOV fields:

```matlab
out = superpose_complex_fields({'.../sub-001_layered_p_complex_trans-01.mat', ...
                                '.../sub-001_layered_p_complex_trans-02.mat'});
out.p_coherent    % |p1 + p2|          coherent (simultaneous) firing
out.p_incoherent  % sqrt(|p1|^2+|p2|^2) async heat source
out.p_max         % max(|p1|, |p2|)    async peak pressure
out = superpose_complex_fields(files, Weights=[1, exp(1i*pi/2)]);   % relative phase of transducer 2
```

One async run (`async_grid_mode: shared`) with `io.save_p_complex: true` therefore also gives the coherent field for any relative amplitude and phase, without another simulation. Fields on different grids (`per_transducer` mode) and fields at different frequencies are not combined. Storing `p_complex` takes a lot of memory on fine grids.

??? note "File contents"
    - `p_complex`: complex amplitude at the fundamental on the simulation grid (the acoustic FOV sub-grid with a [limited FOV](doc_limited_fov.md)), at the simulated drive level, before ISPPA rescaling. The phase refers to t = 0 of the simulation, when every source starts.
    - `p_complex_info`: `phase_reference` (`'simulation_t0'`), `freq_hz`, `record_t0_s`, `dt_s`, `drive_level` (`'simulated'`), `grid_dims`, and `fov_offset_ac` / `full_ac_dims` with a limited FOV.

    Files without `p_complex_info` are accepted only all together, with a warning (valid only for runs with the same time step and length).

## Partial Acoustic FOV

A [limited acoustic FOV](doc_limited_fov.md) (`grid.acoustic_fov_diameter_mm`, `grid.acoustic_fov_length_mm`) works with several transducers:

- **Coherent**: the FOV box is the union of one box per transducer, each enclosing a cylinder of the FOV diameter along that transducer's beam axis. `'auto'` is resolved per transducer (its outer element diameter; its focal depth + `grid.z_beyond_focus_mm`). For opposite transducers the union spans most of the head, so the saving is mostly lateral.
- **Async**: each transducer's run gets its own FOV box. Fields are back-projected onto their full grid before combining, in both grid modes. Partial-FOV and full-grid caches can be mixed.
- **Complex superposition**: each field is zero outside its own box, so the sum is complete only where boxes overlap (`out.coverage` counts the fields per voxel). Make the boxes cover the region around both foci.

The Nexus placement viewer draws every transducer's FOV box.

## Known Limitations

- Simulations without a skull or layered setup (e.g. `medium: water`) only model the first transducer.
- `kWaveArray` mode supports multiple transducers only in non-axisymmetric mode (`grid.axisymmetric = 1` is single-transducer). Matrix-type transducers in `kWaveArray` mode assume a +Z beam axis; see [Transducer geometry](doc_transducer.md#beam-axis-constraint).
- **Different carrier frequencies**: time step and PPW (`source_sensor_setup`) use the *highest* frequency; Mechanical Index (`acoustic_analysis`) uses the *lowest* frequency, which never underestimates MI = p / √f; with `io.save_p_complex: true`, one `*_p_complex_<freq>Hz.mat` is written per frequency.
- Exit-plane metrics (`real_focal_distance_ep_mm`, `real_focal_distance_bowl_mm`, `max_Isppa_after_exit_plane`) refer to the first transducer only.
- Focal-plane heating plots show only the first transducer's focal plane. Always inspect the 3D `maxT` and `CEM43` NIfTI volumes for a complete safety assessment.
- Post-hoc water-only simulations are not implemented for multiple transducers.
- The PRESTUS **[BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain) backend** (`acoustic_backend` / `thermal_backend: babelbrain`) supports one transducer per run. Multi-transducer runs need `kwave` or `jwave` with `thermal_backend: kwave`. See [BabelBrain: Geometry limits](doc_babelbrain.md#geometry-limits).

??? note "Comparison with BabelBrain"
    BabelBrain models one transducer per simulation. It treats several foci of one phased array as time-multiplexed, like PRESTUS async mode: the safety field is the per-voxel maximum of the pressure amplitudes (`CalculateFieldProcess.py`), and the thermal simulation switches between the spots' fields by duty-cycle share (`BHTEMultiplePressureFields` in `ThermalModeling/CalculateTemperatureEffects.py`). It does not add complex fields across foci.

## Subpages

- [Async Multi-Transducer Mode](doc_multitransducer_async.md): for transducers that fire in separate, non-overlapping windows; covers configuration, the physical model, grid modes, outputs and constraints.

See also: [doc_calibration.md](doc_calibration.md), [doc_multiintensity.md](doc_multiintensity.md), [doc_sequential.md](doc_sequential.md), [doc_modules.md](doc_modules.md), [doc_transducer.md](doc_transducer.md)
