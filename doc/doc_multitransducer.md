# Multi-Transducer Modeling

> **Experimental support.** See [this pull request](https://github.com/Donders-Institute/PRESTUS/pull/100) for examples.

Multiple transducers can be specified in layered simulations:

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

### Grid alignment

The default `focal_axis` mode rotates the grid so that the mean beam axis of all transducers aligns with z. The beams are averaged as axes, not directions, so opposed (bilateral) beams do not cancel: z bisects their two axes, and each beam is tilted from z by about half the angle between them. For angled setups converging from one side (e.g. two transducers at 45° or 90°) z is their bisector. Every beam that is not on z is drawn obliquely on the grid; at coarse resolution that costs accuracy (see `async_grid_mode: per_transducer` below for async runs). No configuration is needed — this is the default:

```yaml
grid:
  mode: 'focal_axis'   # default; omit to get the same behaviour
```

If you want the grid aligned to transducer(1)'s focal axis instead (e.g. transducer(1) is the primary beam), use `transducer_axis` explicitly. See [Preprocessing — Orientation](doc_preproc.md#orientation) for the full mode table.

### Grid cropping

The simulation grid is cropped to the bounding box that encloses the CSF-expanded head tissue **and the bowl geometry of every transducer**. See [Preprocessing — Cropping](doc_preproc.md#cropping) for details.

---

## ISPPA targets

`transducer(k).target_isppa_wcm2` is the free-water ISPPA of **transducer k alone**, for every coupling mode. Each transducer has its own device settings (serial, focal-depth calibration, drive amplitude `elem_p_amp`), so a target that differs from what its drive gives needs that transducer's drive rescaled, not the combined field.

| Coupling | How the targets are applied | Sweeps (vector targets) |
|---|---|---|
| `coherent` | Before the simulation (`coherent_drive_calibration`): each transducer's `elem_p_amp` is multiplied by √(target / baseline), where baseline is the free-water ISPPA of that transducer alone at its configured drive. The combined field is simulated at these drives and not rescaled afterwards. | Not available: every sweep point needs its own simulation. Errors (`prestus:isppa:coherentSweep`); run one config per point, or use `async`. |
| `async` | After the per-transducer simulations: each field is scaled to its own target in `combine_async_intensity`, then combined. | Yes: one acoustic run per transducer, one combine + thermal pair per sweep point (see [Async multi-transducer mode](#async-multi-transducer-mode)). |

**Coherent baseline per transducer.** As for a single transducer (`modules.run_water_baseline`, [ISPPA targeting](doc_simulations-acoustic.md#isppa-targeting)): when the target matches the transducer's calibration ISPPA (`tx_calibration.isppa_wcm2`, e.g. from the equipment library) within 1 %, the calibration is the baseline and the drive is kept; otherwise PRESTUS simulates that transducer alone in water (the other transducers at zero amplitude, same grid and source geometry). With N transducers that is up to N extra free-water simulations. A transducer with `target_isppa_wcm2: []` keeps its configured drive.

The result records the per-transducer targets, baselines, their source and the drive scale in `acoustic_provenance.coherent_drive`. The free-water ISPPA of the combined source (both firing, with interference) is not a target; if the run's own free-water baseline runs, it reports that combined value. A coherent multi-transducer result cached before this behaviour existed carries no `coherent_drive` and is refused when targets are set (`prestus:isppa:coherentWithoutDriveCalibration`): re-run its acoustic simulation. `path.acoustic_input` (a prior combined pressure map) cannot be rescaled per transducer and errors when targets are set.

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

In the Nexus, the target is the per-transducer field *Target free-water Isppa* on the transducer page (switch transducers with the selector); *Multi-transducer coupling* chooses between the two rows above.

---

## Temporal relation of multi-transducer firing

In practice, multiple transducers may fire simultaneously, with some temporal offset, or fully independently with different duty cycles or pulse sequences. PRESTUS has native support only for the simultaneous coherent case. The other scenarios require workarounds described below.

### Simultaneous coherent firing (native support)

When multiple transducers are listed in the configuration, PRESTUS combines them into a single k-Wave source matrix. Their pressure fields superpose implicitly as k-Wave solves the linear acoustic PDE — this is the physically correct treatment for transducers that fire at the same frequency with a known phase relationship. All transducers must share the same fundamental frequency; a warning is issued if they differ.

This is the default mode and requires no additional configuration. Each transducer's `target_isppa_wcm2` sets its own drive before the simulation (see [ISPPA targets](#isppa-targets)).

### Asynchronous / temporally non-overlapping firing

If the transducers fire in separate, non-overlapping time windows, their pressure fields never coexist. Two distinct combination rules apply:

- **Acoustic safety metrics (ISPPA, MI, peak pressure):** the relevant quantity is the per-voxel maximum across transducers.
- **Thermal heat deposition:** time-averaged heating accumulates from all transducers; the combined heat source is the incoherent intensity sum (currently assuming identical duty cycles).

PRESTUS provides a dedicated `transducer_coupling: async` mode that handles both combination rules automatically. See [Async mode](#async-multi-transducer-mode) below for full details.

### Staggered / time-offset firing

If the transducers fire in sequence (one after another, or with partial temporal overlap), the closest approximation is to use [Sequential Simulations](doc_sequential.md). Run one transducer's heating simulation, carry over the resulting temperature and CEM43 maps via `adopted_heatmap` / `adopted_cem43`, then run the next transducer's heating simulation starting from that thermal state. This correctly captures the cumulative thermal history but assumes the acoustic fields do not overlap in time.

---

## Complex pressure field

PRESTUS extracts the complex pressure amplitude (magnitude and phase at the fundamental frequency) from the steady-state time series and makes it available as an optional cache output:

```yaml
io:
  save_p_complex: true
```

This writes `sub-XXX_<medium>_p_complex<affix>.mat` with two variables:

- `p_complex`: complex array on the simulation grid (the acoustic FOV sub-grid when a [limited FOV](doc_limited_fov.md) is active), at the simulated drive level, before any ISPPA rescaling. The magnitude is the pressure amplitude at the fundamental; the phase is referred to t = 0 of the simulation, when every source signal starts.
- `p_complex_info`: `phase_reference` (`'simulation_t0'`), `freq_hz`, `record_t0_s` and `dt_s` (start and step of the recorded window), `drive_level` (`'simulated'`), `grid_dims`, and `fov_offset_ac` / `full_ac_dims` when a limited FOV was used.

Because the phase has a common time origin, fields from separate runs can be added: acoustic propagation in PRESTUS is linear (no nonlinearity parameter is set), so the field of several transducers firing together is the sum of their fields simulated one at a time. `superpose_complex_fields` does this on one shared grid, back-projecting limited-FOV fields, and also returns the incoherent and maximum combinations used by async mode:

```matlab
out = superpose_complex_fields({'.../sub-001_layered_p_complex_trans-01.mat', ...
                                '.../sub-001_layered_p_complex_trans-02.mat'});
out.p_coherent    % |p1 + p2|          coherent (simultaneous) firing
out.p_incoherent  % sqrt(|p1|^2+|p2|^2) async heat source
out.p_max         % max(|p1|, |p2|)    async peak pressure
out = superpose_complex_fields(files, Weights=[1, exp(1i*pi/2)]);   % relative phase of transducer 2
```

The per-transducer runs of async mode (`async_grid_mode: shared`) are on the same grid as a coherent run, so one async run with `io.save_p_complex: true` gives the async result and, through `superpose_complex_fields`, the coherent field for any relative drive amplitude and phase without another simulation. Fields on different grids (`per_transducer` mode) are not added: interpolating the phase of a field sampled at a few points per wavelength is not accurate. Fields at different frequencies do not interfere on average and are rejected.

Storing `p_complex` records the pressure time series of the last periods over the whole (FOV) grid during the simulation, which takes a lot of memory on fine grids.

> Before this version, `p_complex` was saved as a column vector (k-Wave returns the recorded time series as sensor points × time) and its phase referred to the start of each run's recording window. Old files lack `p_complex_info`; `superpose_complex_fields` accepts them only all together and with a warning (valid only for runs with the same time step and length).

---

## Partial acoustic FOV

A [limited acoustic FOV](doc_limited_fov.md) (`grid.acoustic_fov_diameter_mm`, `grid.acoustic_fov_length_mm`) works with several transducers:

- **Coherent coupling**: all transducers are in one k-Wave run, so the FOV box must hold every beam. It is the union of one box per transducer, each enclosing a cylinder of the FOV diameter that runs the FOV length from that transducer's centre along its own beam axis. `'auto'` is resolved per transducer (its outer element diameter; its focal depth + `grid.z_beyond_focus_mm`). For transducers on opposite sides of the head the union spans most of the head along the inter-transducer axis, so the saving is mostly lateral.
- **Async coupling**: each transducer's acoustic run gets its own FOV box around its own beam. The fields are back-projected onto their full grid and then combined as without a FOV, in both `shared` and `per_transducer` grid mode. Partial-FOV and full-grid caches can be mixed.
- **Complex superposition**: per-transducer `p_complex` fields from partial-FOV runs can be added with `superpose_complex_fields`. Each field is zero outside its own FOV box, so the sum is only complete where the boxes overlap (`out.coverage` counts the fields covering each voxel). Make the boxes cover the region of interest around both foci.

The Nexus placement viewer draws every transducer's FOV box.

---

## Comparison with BabelBrain

BabelBrain models one transducer per simulation. Several foci of one phased array (*multi-point* focusing) are treated as time-multiplexed, like PRESTUS async mode: each focal spot is simulated separately, the displayed and safety field is the per-voxel maximum of the pressure amplitudes (`CalculateFieldProcess.py`), and the thermal simulation switches between the spots' magnitude fields with each spot on for its share of the duty cycle (`BHTEMultiplePressureFields` in `ThermalModeling/CalculateTemperatureEffects.py`). The fields are not added as complex fields. Complex fields are used inside one BabelBrain simulation, for the Rayleigh propagation of the source and for phase-conjugation refocusing.

The PRESTUS **BabelBrain backend** (`acoustic_backend` / `thermal_backend: babelbrain`) is limited to one transducer per run, coherent and async alike; multi-transducer runs need `kwave` or `jwave` with `thermal_backend: kwave`. See [BabelBrain — Geometry limits](doc_babelbrain.md#geometry-limits).

---

## Known limitations

- The BabelBrain backend supports one transducer per run (see [Comparison with BabelBrain](#comparison-with-babelbrain)).
- Simulations without a skull or layered setup (e.g. `medium: water`) only model the first transducer.
- `kWaveArray` mode supports multiple transducers in non-axisymmetric mode. Axisymmetric mode (`grid.axisymmetric = 1`) remains single-transducer only. Matrix-type transducers in `kWaveArray` mode assume a +Z beam axis; see [Transducer geometry](doc_transducer.md#beam-axis-constraint) for details.
- **Different carrier frequencies** across transducers are handled as follows:
  - **Grid time step and PPW** (`source_sensor_setup`): sized using the *highest* frequency present, ensuring all transducers are correctly sampled. A message is printed when frequencies differ.
  - **Mechanical Index** (`acoustic_analysis`): computed using the *lowest* frequency present. This is the conservative choice — MI = p / √f is largest at low frequencies, so the resulting map never underestimates MI at any voxel. A message is printed when frequencies differ.
  - **Complex pressure field** (`p_complex`, optional): when `io.save_p_complex: true`, the DFT-based extraction is run once per unique frequency; each produces a separate cache file named `*_p_complex_<freq>Hz.mat`. For a single shared frequency the original filename is preserved.
- Exit-plane metrics (focal distance, `max_Isppa_after_exit_plane`, `real_focal_distance`) are computed only with respect to the first transducer.
- Focal-plane time-course heating plots reflect only the focal plane of the first transducer and may miss hotspots near other beams. Always inspect the 3D `maxT` and `CEM43` NIfTI volumes for a complete safety assessment.
- Post-hoc water-only simulations are not implemented for multiple transducers.

---

## Async multi-transducer mode

> **Experimental.** This mode is intended for setups where two or more transducers fire asynchronously — i.e., with non-overlapping duty cycles or alternating pulses — so that their pressure fields never superpose coherently. Configuring this mode requires a clear understanding of the intended acoustic exposure and the linear-regime assumptions described below.

### When to use

| Mode | When to use |
|---|---|
| **Coherent** (`transducer_coupling: coherent`, default) | Transducers fire simultaneously; pressure fields interfere. Modelled in a single k-Wave run with all transducer sources active. |
| **Async** (`transducer_coupling: async`) | Transducers fire in separate, non-overlapping windows; their pressure fields are independent. Each is simulated separately; intensities are summed for thermal analysis and the per-voxel pressure maximum is used for acoustic safety metrics. |

### Physical model

Because only one transducer is active at any instant, the two physically relevant quantities require different combination rules.

**Acoustic safety metrics (ISPPA, MI, peak pressure)**

The instantaneous peak pressure at any voxel is the maximum over all transducer fields:

$$p_\mathrm{safety} = \max_i\bigl(p_{i,\mathrm{scaled}}\bigr)$$

where $p_{i,\mathrm{scaled}} = p_i \cdot \sqrt{I_{i,\mathrm{target}} / I_{i,\mathrm{baseline}}}$.

ISPPA, MI, and peak-pressure metrics in the output CSV and plots are derived from $p_\mathrm{safety}$ (`sensor_data.p_max_async` in the combined cache).

**Thermal heat deposition**

Time-averaged heat deposition accumulates from all transducers. Assuming equal, non-overlapping duty cycles that together fill the protocol duty cycle, the combined heat source is the incoherent intensity sum:

$$Q_\mathrm{total} \propto \sum_{i=1}^{N} I_i = \sum_{i=1}^{N} \frac{p_{i,\mathrm{scaled}}^2}{2\rho c}$$

This is stored as `sensor_data.p_max_all` in the combined cache and passed directly to the bioheat solver:

$$p_\mathrm{thermal} = \sqrt{2 \cdot I_\mathrm{total} \cdot \rho c}$$

> **Limitation.** The thermal model does not yet explicitly represent the relative firing schedule of the two transducers. The intensity-sum approximation is exact when both DCs are equal and non-overlapping; it becomes an over-estimate when $\mathrm{DC}_1 + \mathrm{DC}_2 < $ total protocol DC.

### Grid geometry mode

By default (`simulation.async_grid_mode: 'shared'`), all transducers share one grid — see
[Grid alignment](#grid-alignment) and [Grid cropping](#grid-cropping) above. This means each
transducer's own acoustic sub-run (`_trans-01`, `_trans-02`, …) is simulated on a grid oriented
along the *mean* beam axis and cropped to enclose every transducer's bowl, not just its own. A beam
that is not on the grid axis is drawn obliquely: its bowl and the skull it crosses are voxelised at
an angle. The physics is the same, but the discretisation is not, so a per-transducer result
(e.g. focal distance, ISPPA, even the free-water ISPPA) can differ from what the same transducer
shows run alone (solo/sequential mode). The difference shrinks with grid resolution; at 2 mm
(2.5 points per wavelength at 300 kHz) a beam 41° off axis lost 16% of its free-water ISPPA.

Setting `simulation.async_grid_mode: 'per_transducer'` instead gives each transducer its own
optimally aligned, minimally cropped grid — identical to what solo/sequential mode would produce
for that transducer. `combine_async_intensity` then resamples every transducer's field onto
transducer(1)'s grid (the reference frame) via each run's `planimg` sidecar (sim-voxel → world RAS
mm affine), so per-voxel max/sum combination remains correct even though each field was simulated
on a different grid. That resampling costs accuracy in turn (linear interpolation of a rotated
field; parts of a field outside transducer(1)'s grid are dropped). Use this mode when per-transducer
characterization (focal distance, ISPPA, focal volume) needs to match the solo/sequential result,
or when transducers are widely separated or steeply angled (where the shared union-bbox grid would
otherwise grow large or badly aligned).

```yaml
simulation:
  transducer_coupling: async
  async_grid_mode: per_transducer   # default: 'shared'
```

**Scalability.** `per_transducer` mode scales *better* than `shared` mode as the number of
transducers grows: each transducer's grid stays sized to its own bowl+tissue crop regardless of N,
whereas the shared grid's union bounding box grows with the spread of all N transducers. The added
cost is N independent (HPC-parallel) preproc jobs instead of one, plus one world-space affine
resample per transducer in the combine step — linear in N, and folded into the already
lightweight, CPU-only combine stage.

### Configuration

**Basic setup — two transducers, one thermal target each**

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

modules:
  run_water_baseline: 1   # required to record free-water ISPPA provenance
```

Both transducers must share the same segmentation (same `subject_id`, `simulation.medium`, and output directory). The simulation grid is determined by the union of both transducer positions.

**Multi-ISPPA sweep**

Setting `target_isppa_wcm2` to a vector triggers a sweep. All transducers must carry either a scalar target (held fixed) or a vector of the same length:

```yaml
transducer:
  - name: CTX500-A
    target_isppa_wcm2: [3.0, 5.0, 10.0]   # swept across 3 points

  - name: CTX500-B
    target_isppa_wcm2: 2.0                  # fixed across all 3 sweep points
```

This produces three combine+thermal pairs (sweep 1: A→3.0, B→2.0; sweep 2: A→5.0, B→2.0; sweep 3: A→10.0, B→2.0).

### Pipeline stages

```
[Stage 0]       Preproc grid(s).
                 async_grid_mode='shared' (default): one shared grid.
                 async_grid_mode='per_transducer': N independent grids, one
                 per transducer (parallel on HPC), affixes
                 _desc-asyncPreproc-t01, …, -tNN.

[Stage 1…N]    Acoustic simulation — one per transducer (parallel on HPC)
               Each runs at its configured elem_p_amp and records a water baseline.
               Affixes: _trans-01, _trans-02, …, _trans-NN

[Stage N+1…    Dual-field combination — one per sweep point j
 N+M]          Scales each p_i to target_i(j); writes combined acoustic cache
               containing p_max_all (intensity sum, for thermal) and
               p_max_async (per-voxel pressure max, for acoustic safety).
               In per_transducer mode, each field is first resampled onto
               transducer(1)'s grid via its planimg world-space affine.
               Affixes: _desc-asyncCombined (M=1) or _desc-asyncCombinedS01… (M>1)

[Stage N+M+1…  Thermal simulation — one per sweep point j
 N+2M]         Loads combined cache, runs kWaveDiffusion, generates report.
               Affixes: _desc-asyncJoint (M=1) or _desc-asyncJointS01… (M>1)
```

`prestus_pipeline_start` returns the submitted job IDs per stage as a second output, so later jobs can wait for an async run, e.g. runs that reuse its grid and per-transducer fields:

```matlab
[~, jobs] = prestus_pipeline_start(parameters);   % jobs.preproc, .acoustic, .combine, .thermal, .water
next.hpc.depend_job_ids = [jobs.preproc, jobs.acoustic];
prestus_pipeline_start(next);
```

On HPC (SLURM/qsub): in `shared` mode, all N acoustic stages are submitted simultaneously after the
one preproc job; in `per_transducer` mode, N preproc jobs are submitted simultaneously and each
acoustic job depends only on its own transducer's preproc job (not on the others). Each combine
stage depends on all N acoustic jobs (`afterok`); each thermal stage depends on its corresponding
combine stage.

### Output files

For two transducers, single thermal target each:

```
cache/sub-001_layered_results_trans-01.mat
cache/sub-001_layered_results_trans-02.mat
cache/sub-001_layered_results_desc-asyncCombined.mat
sub-001_layered_desc-asyncJoint.csv
sub-001_layered_report_desc-asyncJoint.html
```

With an `io.output_affix` (e.g. `run1`), every affix starts with it: `_desc-run1_trans-01`, `_desc-run1_desc-asyncJoint`. The per-transducer runs also write their own CSV and acoustic NIfTIs (at the simulated drive level), which the Nexus Results tab shows under the joint result.

For a 3-point sweep:

```
cache/sub-001_layered_results_trans-01.mat
cache/sub-001_layered_results_trans-02.mat
cache/sub-001_layered_results_desc-asyncCombinedS01.mat
cache/sub-001_layered_results_desc-asyncCombinedS02.mat
cache/sub-001_layered_results_desc-asyncCombinedS03.mat
sub-001_layered_desc-asyncJointS01.csv
sub-001_layered_desc-asyncJointS02.csv
sub-001_layered_desc-asyncJointS03.csv
```

### Activation

Async mode is activated automatically when both conditions are met:

```matlab
parameters.simulation.transducer_coupling = 'async';  % and
numel(parameters.transducer) > 1
```

No special entry point is needed:

```matlab
parameters = load_parameters('config_subject.yaml', subject_id);
prestus_pipeline_start(parameters);   % dispatches to async_transducer_pipeline
```

### Constraints and limitations

- **3-D grid required.** Axisymmetric (2-D) simulations are not supported.
- **Shared head model.** Both transducers must use the same segmented head grid.
- **Linear regime.** The same validity conditions as single-transducer ISPPA scaling apply (see [doc_calibration.md](doc_calibration.md)).
- **`modules.run_water_baseline: 1` is required.** Without it, pressure scaling is skipped.
- **Async mode composes with uncertainty mode.** With `simulation.uncertainty: true`, each uncertainty variant runs its own async pipeline (stage 0 included) with its medium, and one uncertainty report is written per sweep point. See [doc_uncertainty.md](doc_uncertainty.md#combining-with-a-multi-isppa-sweep-or-async-transducers).

### HPC options

```matlab
options.acoustic_timelimit   = '06:00:00';   % per acoustic job
options.combine_timelimit    = '00:30:00';   % per combine job (CPU-only, lightweight)
options.thermal_timelimit    = '04:00:00';   % per thermal job
options.acoustic_memorylimit = 60;           % GB
options.thermal_memorylimit  = 40;           % GB
options.acoustic_partition   = 'gpu';
options.combine_partition    = 'cpu';        % combine stage needs no GPU
```

### Resumption

Each stage checks for its sentinel output file before submitting or running. Re-running `prestus_pipeline_start` with the same config resumes from the last incomplete stage. With `io.overwrite_files: always` (default), the preproc, acoustic and combine stages resume from their caches while the joint and water-validation stages always re-run; `force` re-runs every stage. See [Overwrite behaviour](doc_outputs.md#overwrite-behaviour).

### Provenance

The combined acoustic cache records the following fields in `acoustic_provenance`:

| Field | Description |
|---|---|
| `combined_from` | Cell array of source cache file paths |
| `combination_mode` | `'async_dual_field'` — `p_max_all` holds the thermal intensity sum; `p_max_async` holds the per-voxel pressure maximum for acoustic safety |
| `combination_targets` | Per-transducer target ISPPAs used for scaling [W/cm²] |
| `combination_baselines` | Per-transducer free-water baselines at the simulated drive level [W/cm²] |
| `combination_time` | Timestamp of the combination step |

---

See also: [doc_calibration.md](doc_calibration.md), [doc_multiintensity.md](doc_multiintensity.md), [doc_sequential.md](doc_sequential.md), [doc_modules.md](doc_modules.md), [doc_transducer.md](doc_transducer.md)
