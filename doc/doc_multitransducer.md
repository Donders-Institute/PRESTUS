# Multi-Transducer Modeling

> **Experimental support.** See [this pull request](https://github.com/Donders-Institute/PRESTUS/pull/100) for examples.

Multiple transducers can be specified in layered simulations:

```yaml
transducer:
  - name: right
    elem_n: 10
    # ...

  - name: left
    elem_n: 10
    # ...
```

---

## Temporal relation of multi-transducer firing

In practice, multiple transducers may fire simultaneously, with some temporal offset, or fully independently with different duty cycles or pulse sequences. PRESTUS has native support only for the simultaneous coherent case. The other scenarios require workarounds described below.

### Simultaneous coherent firing (native support)

When multiple transducers are listed in the configuration, PRESTUS combines them into a single k-Wave source matrix. Their pressure fields superpose implicitly as k-Wave solves the linear acoustic PDE — this is the physically correct treatment for transducers that fire at the same frequency with a known phase relationship. All transducers must share the same fundamental frequency; a warning is issued if they differ.

This is the default mode and requires no additional configuration.

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

This writes `sub-XXX_<medium>_p_complex<affix>.mat` containing a single variable `p_complex` — a complex-valued array where magnitude encodes peak pressure amplitude and phase encodes the spatial phase of the wavefield at the fundamental frequency.

---

## Known limitations

- Simulations without a skull or layered setup (e.g. `medium: water`) only model the first transducer.
- `kWaveArray` mode is incompatible with multiple transducers.
- Different carrier frequencies across transducers are not supported.
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
      elem_amp: 100000

  - name: CTX500-B
    target_isppa_wcm2: 3.0
    freq_hz: 500000
    trans_pos: [x2, y2, z2]
    focus_pos: [fx2, fy2, fz2]
    annular:
      elem_amp: 100000

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
[Stage 1…N]    Acoustic simulation — one per transducer (parallel on HPC)
               Each runs at its configured elem_amp and records a water baseline.
               Affixes: _desc-t01, _desc-t02, …, _desc-tNN

[Stage N+1…    Dual-field combination — one per sweep point j
 N+M]          Scales each p_i to target_i(j); writes combined acoustic cache
               containing p_max_all (intensity sum, for thermal) and
               p_max_async (per-voxel pressure max, for acoustic safety).
               Affixes: _desc-asyncCombined (M=1) or _desc-asyncCombinedS01… (M>1)

[Stage N+M+1…  Thermal simulation — one per sweep point j
 N+2M]         Loads combined cache, runs kWaveDiffusion, generates report.
               Affixes: _desc-asyncThermal (M=1) or _desc-asyncThermalS01… (M>1)
```

On HPC (SLURM/qsub): all N acoustic stages are submitted simultaneously; each combine stage depends on all N acoustic jobs (`afterok`); each thermal stage depends on its corresponding combine stage.

### Output files

For two transducers, single thermal target each:

```
cache/sub-001_layered_results_desc-t01.mat
cache/sub-001_layered_results_desc-t02.mat
cache/sub-001_layered_results_desc-asyncCombined.mat
sub-001_layered_desc-asyncThermal.csv
sub-001_layered_desc-asyncThermal_report.html
```

For a 3-point sweep:

```
cache/sub-001_layered_results_desc-t01.mat
cache/sub-001_layered_results_desc-t02.mat
cache/sub-001_layered_results_desc-asyncCombinedS01.mat
cache/sub-001_layered_results_desc-asyncCombinedS02.mat
cache/sub-001_layered_results_desc-asyncCombinedS03.mat
sub-001_layered_desc-asyncThermalS01.csv
sub-001_layered_desc-asyncThermalS02.csv
sub-001_layered_desc-asyncThermalS03.csv
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
- **Async mode is mutually exclusive with uncertainty mode.** Use coherent mode with uncertainty if skull-property uncertainty quantification is required.

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

Each stage checks for its sentinel output file before submitting or running. Re-running `prestus_pipeline_start` with the same config resumes from the last incomplete stage.

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
