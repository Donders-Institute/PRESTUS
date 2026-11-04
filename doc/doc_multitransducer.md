# Multi-Transducer Modeling

!!! note "At a glance"
    - **What:** Several transducers in one layered simulation, firing together (`coherent`) or in separate windows (`async`).
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
| Simultaneously, same frequency, known phase | `coherent` (default): one k-Wave run with all sources active; the fields interfere. All transducers must share the fundamental frequency (a warning is issued otherwise). |
| In separate, non-overlapping windows | `async`: see [Async multi-transducer mode](#async-multi-transducer-mode). |
| In sequence, or with partial overlap | Use [Sequential Simulations](doc_sequential.md): carry each heating run's temperature and CEM43 maps into the next via `adopted_heatmap` / `adopted_cem43`. This assumes the acoustic fields do not overlap in time. |

In the Nexus, *Multi-transducer coupling* selects the row above and *Target free-water Isppa* is set per transducer on the transducer page.

## Grid

### Grid alignment

The default `focal_axis` mode rotates the grid so that the mean beam axis of all transducers aligns with z. Beams are averaged as axes, not directions, so opposed (bilateral) beams do not cancel: z bisects their two axes. Each beam off z is drawn obliquely on the grid, which costs accuracy at coarse resolution (for async runs, see `async_grid_mode: per_transducer` below).

```yaml
grid:
  mode: 'focal_axis'   # default; omit to get the same behaviour
```

To align the grid to transducer(1)'s focal axis instead, set `transducer_axis`. See [Preprocessing: Orientation](doc_preproc.md#orientation) for the mode table.

### Grid cropping

The grid is cropped to the bounding box that encloses the CSF-expanded head tissue **and the bowl geometry of every transducer**. See [Preprocessing: Cropping](doc_preproc.md#cropping).

## ISPPA targets

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
| `async` | After the per-transducer simulations: each field is scaled to its own target in `combine_async_intensity`, then combined. | Yes: one acoustic run per transducer, one combine + thermal pair per sweep point (see [Async multi-transducer mode](#async-multi-transducer-mode)). |

**Coherent baseline.** The rule is the same as for one transducer ([ISPPA targeting](doc_simulations-acoustic.md#isppa-targeting)). If the target matches the transducer's calibration ISPPA (`tx_calibration.isppa_wcm2`) within 1 %, the calibration is the baseline. Otherwise PRESTUS simulates that transducer alone in water (others at zero amplitude, same grid). With N transducers that is up to N extra free-water simulations. A transducer with `target_isppa_wcm2: []` keeps its configured drive.

??? note "Details: provenance and cache checks"
    The result records targets, baselines, their source and the drive scale in `acoustic_provenance.coherent_drive`. The combined free-water ISPPA (both firing) is reported but is not a target. A coherent cache without `coherent_drive` is refused when targets are set (`prestus:isppa:coherentWithoutDriveCalibration`): re-run its acoustic simulation. `path.acoustic_input` cannot be rescaled per transducer and errors when targets are set.

## Complex pressure field

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

??? note "Details: file contents"
    - `p_complex`: complex amplitude at the fundamental on the simulation grid (the acoustic FOV sub-grid with a [limited FOV](doc_limited_fov.md)), at the simulated drive level, before ISPPA rescaling. The phase refers to t = 0 of the simulation, when every source starts.
    - `p_complex_info`: `phase_reference` (`'simulation_t0'`), `freq_hz`, `record_t0_s`, `dt_s`, `drive_level` (`'simulated'`), `grid_dims`, and `fov_offset_ac` / `full_ac_dims` with a limited FOV.

    Files without `p_complex_info` are accepted only all together, with a warning (valid only for runs with the same time step and length).

## Partial acoustic FOV

A [limited acoustic FOV](doc_limited_fov.md) (`grid.acoustic_fov_diameter_mm`, `grid.acoustic_fov_length_mm`) works with several transducers:

- **Coherent**: the FOV box is the union of one box per transducer, each enclosing a cylinder of the FOV diameter along that transducer's beam axis. `'auto'` is resolved per transducer (its outer element diameter; its focal depth + `grid.z_beyond_focus_mm`). For opposite transducers the union spans most of the head, so the saving is mostly lateral.
- **Async**: each transducer's run gets its own FOV box. Fields are back-projected onto their full grid before combining, in both grid modes. Partial-FOV and full-grid caches can be mixed.
- **Complex superposition**: each field is zero outside its own box, so the sum is complete only where boxes overlap (`out.coverage` counts the fields per voxel). Make the boxes cover the region around both foci.

The Nexus placement viewer draws every transducer's FOV box.

## Known limitations

- Simulations without a skull or layered setup (e.g. `medium: water`) only model the first transducer.
- `kWaveArray` mode supports multiple transducers only in non-axisymmetric mode (`grid.axisymmetric = 1` is single-transducer). Matrix-type transducers in `kWaveArray` mode assume a +Z beam axis; see [Transducer geometry](doc_transducer.md#beam-axis-constraint).
- **Different carrier frequencies**: time step and PPW (`source_sensor_setup`) use the *highest* frequency; Mechanical Index (`acoustic_analysis`) uses the *lowest* frequency, which never underestimates MI = p / √f; with `io.save_p_complex: true`, one `*_p_complex_<freq>Hz.mat` is written per frequency.
- Exit-plane metrics (focal distance, `max_Isppa_after_exit_plane`, `real_focal_distance`) refer to the first transducer only.
- Focal-plane heating plots show only the first transducer's focal plane. Always inspect the 3D `maxT` and `CEM43` NIfTI volumes for a complete safety assessment.
- Post-hoc water-only simulations are not implemented for multiple transducers.
- The PRESTUS **BabelBrain backend** (`acoustic_backend` / `thermal_backend: babelbrain`) supports one transducer per run. Multi-transducer runs need `kwave` or `jwave` with `thermal_backend: kwave`. See [BabelBrain: Geometry limits](doc_babelbrain.md#geometry-limits).

??? note "Details: comparison with BabelBrain"
    BabelBrain models one transducer per simulation. It treats several foci of one phased array as time-multiplexed, like PRESTUS async mode: the safety field is the per-voxel maximum of the pressure amplitudes (`CalculateFieldProcess.py`), and the thermal simulation switches between the spots' fields by duty-cycle share (`BHTEMultiplePressureFields` in `ThermalModeling/CalculateTemperatureEffects.py`). It does not add complex fields across foci.

## Async multi-transducer mode

!!! warning "Experimental"
    For transducers that fire with non-overlapping duty cycles or alternating pulses, so their fields never superpose coherently. Configure it only with a clear view of the intended exposure and the linear-regime assumptions below.

### Configuration

Async mode runs when `transducer_coupling` is `'async'` and there is more than one transducer. `prestus_pipeline_start(parameters)` dispatches to it; no other entry point is needed.

```yaml
simulation:
  transducer_coupling: async

transducer:
  - name: CTX500-A
    target_isppa_wcm2: 5.0
    freq_hz: 500000
    trans_pos: [x1, y1, z1]
    focus_pos: [fx1, fy1, fz1]
    annular:
      elem_p_amp: 100000

  - name: CTX500-B
    target_isppa_wcm2: 3.0
    freq_hz: 500000
    trans_pos: [x2, y2, z2]
    focus_pos: [fx2, fy2, fz2]
    annular:
      elem_p_amp: 100000
```

The per-transducer acoustic runs always record a free-water baseline, which the combine step needs for scaling.

**Multi-ISPPA sweep.** A vector `target_isppa_wcm2` triggers a sweep.

??? note "Details: multi-ISPPA sweep"
    Every transducer carries either a scalar (held fixed) or a vector of the same length:

    ```yaml
    transducer:
      - name: CTX500-A
        target_isppa_wcm2: [3.0, 5.0, 10.0]   # swept across 3 points

      - name: CTX500-B
        target_isppa_wcm2: 2.0                  # fixed across all 3 sweep points
    ```

    This gives three combine + thermal pairs (A→3.0, 5.0, 10.0; B→2.0 each time).

### Physical model

**Acoustic safety metrics (ISPPA, MI, peak pressure)** use the per-voxel maximum of the scaled fields. **Thermal heat deposition** uses the incoherent intensity sum, assuming equal, non-overlapping duty cycles that together fill the protocol duty cycle.

!!! warning "Limitation"
    The thermal model does not represent the firing schedule. The intensity sum is exact when both DCs are equal and non-overlapping; it over-estimates when $\mathrm{DC}_1 + \mathrm{DC}_2 < $ total protocol DC.

??? note "Details: equations"
    Acoustic safety metrics:

    $$p_\mathrm{safety} = \max_i\bigl(p_{i,\mathrm{scaled}}\bigr)$$

    where $p_{i,\mathrm{scaled}} = p_i \cdot \sqrt{I_{i,\mathrm{target}} / I_{i,\mathrm{baseline}}}$. This is `sensor_data.p_max_async` in the combined cache.

    Thermal heat deposition:

    $$Q_\mathrm{total} \propto \sum_{i=1}^{N} I_i = \sum_{i=1}^{N} \frac{p_{i,\mathrm{scaled}}^2}{2\rho c}$$

    This is `sensor_data.p_max_all`, passed to the bioheat solver as:

    $$p_\mathrm{thermal} = \sqrt{2 \cdot I_\mathrm{total} \cdot \rho c}$$

### Grid geometry mode

| `simulation.async_grid_mode` | Grid | Use when |
|---|---|---|
| `'shared'` (default) | One grid for all transducers, aligned and cropped as in [Grid alignment](#grid-alignment) and [Grid cropping](#grid-cropping). | Default. |
| `'per_transducer'` | Each transducer gets the grid a solo run would use. `combine_async_intensity` resamples every field onto transducer(1)'s grid via each run's `planimg` sidecar. | Per-transducer focal distance, ISPPA or focal volume must match a solo run, or the transducers are widely separated or steeply angled. |

??? note "Details: grid mode choice and accuracy"
    ```yaml
    simulation:
      transducer_coupling: async
      async_grid_mode: per_transducer   # default: 'shared'
    ```

    In `shared` mode an off-axis beam is voxelised at an angle, so its per-transducer results can differ from a solo run. The difference shrinks with resolution; at 2 mm (2.5 points per wavelength at 300 kHz) a beam 41° off axis lost 16% of its free-water ISPPA. In `per_transducer` mode the resampling (linear interpolation; parts outside transducer(1)'s grid are dropped) costs accuracy instead. `per_transducer` scales better with N: each grid stays sized to its own transducer, at the cost of N preproc jobs.

### Output files

Two transducers, one target each:

```
cache/sub-001_layered_results_trans-01.mat
cache/sub-001_layered_results_trans-02.mat
cache/sub-001_layered_results_desc-asyncCombined.mat
sub-001_layered_desc-asyncJoint.csv
sub-001_layered_report_desc-asyncJoint.html
```

??? note "Details: sweep and affix naming"
    A 3-point sweep writes `..._desc-asyncCombinedS01.mat` to `S03` and `sub-001_layered_desc-asyncJointS01.csv` to `S03`.

    With `io.output_affix` (e.g. `run1`), every affix starts with it: `_desc-run1_trans-01`, `_desc-run1_desc-asyncJoint`. The per-transducer runs also write their own CSV and acoustic NIfTIs at the simulated drive level; the Nexus Results tab shows them under the joint result.

**Resumption.** Re-running `prestus_pipeline_start` with the same config resumes from the last incomplete stage.

??? note "Details: resumption"
    With `io.overwrite_files: always` (default), preproc, acoustic and combine resume from their caches; the joint and water-validation stages always re-run. `force` re-runs every stage. See [Overwrite behaviour](doc_outputs.md#overwrite-behaviour).

??? note "Details: pipeline stages, HPC options and provenance"
    | Stage | Count | Affix |
    |---|---|---|
    | Preproc | 1 (`shared`) or N (`per_transducer`, parallel) | `_desc-asyncPreproc-t01`, …, `-tNN` (`per_transducer`) |
    | Acoustic, at configured `elem_p_amp`, with water baseline | N (parallel) | `_trans-01`, …, `_trans-NN` |
    | Combine: scale each field to target, write `p_max_all` and `p_max_async` | M sweep points | `_desc-asyncCombined` (M=1) or `_desc-asyncCombinedS01`… |
    | Thermal (`kWaveDiffusion`) and report | M | `_desc-asyncJoint` (M=1) or `_desc-asyncJointS01`… |

    On HPC, each acoustic job depends on its preproc job, each combine job on all N acoustic jobs (`afterok`), and each thermal job on its combine job. `prestus_pipeline_start` returns the job IDs per stage, so later runs can wait on them:

    ```matlab
    [~, jobs] = prestus_pipeline_start(parameters);   % jobs.preproc, .acoustic, .combine, .thermal, .water
    next.hpc.depend_job_ids = [jobs.preproc, jobs.acoustic];
    prestus_pipeline_start(next);
    ```

    HPC options:

    ```matlab
    options.acoustic_timelimit   = '06:00:00';   % per acoustic job
    options.combine_timelimit    = '00:30:00';   % per combine job (CPU-only, lightweight)
    options.thermal_timelimit    = '04:00:00';   % per thermal job
    options.acoustic_memorylimit = 60;           % GB
    options.thermal_memorylimit  = 40;           % GB
    options.acoustic_partition   = 'gpu';
    options.combine_partition    = 'cpu';        % combine stage needs no GPU
    ```

    The combined cache records in `acoustic_provenance`:

    | Field | Description |
    |---|---|
    | `combined_from` | Cell array of source cache file paths |
    | `combination_mode` | `'async_dual_field'`: `p_max_all` holds the thermal intensity sum; `p_max_async` holds the per-voxel pressure maximum for acoustic safety |
    | `combination_targets` | Per-transducer target ISPPAs used for scaling [W/cm²] |
    | `combination_baselines` | Per-transducer free-water baselines at the simulated drive level [W/cm²] |
    | `combination_time` | Timestamp of the combination step |

### Constraints

- **3-D grid required.** Axisymmetric (2-D) simulations are not supported.
- **Shared head model.** All transducers use the same segmentation (same `subject_id`, `simulation.medium` and output directory).
- **Linear regime.** The validity conditions of single-transducer ISPPA scaling apply (see [doc_calibration.md](doc_calibration.md)).
- **Uncertainty mode.** With `simulation.uncertainty: true`, each variant runs its own async pipeline (preproc included), and one uncertainty report is written per sweep point. See [doc_uncertainty.md](doc_uncertainty.md#combining-with-a-multi-isppa-sweep-or-async-transducers).

See also: [doc_calibration.md](doc_calibration.md), [doc_multiintensity.md](doc_multiintensity.md), [doc_sequential.md](doc_sequential.md), [doc_modules.md](doc_modules.md), [doc_transducer.md](doc_transducer.md)
