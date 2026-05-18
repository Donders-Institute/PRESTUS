# Async multi-transducer mode

> **Advanced feature.** This mode is intended for setups where two or more
> transducers fire asynchronously — i.e., with non-overlapping duty cycles
> or alternating pulses — so that their pressure fields never superpose
> coherently. Configuring this mode requires a clear understanding of the
> intended acoustic exposure and the linear-regime assumptions described below.

---

## When to use this mode

PRESTUS supports two ways to model multiple transducers:

| Mode | When to use |
|---|---|
| **Coherent** (`transducer_coupling: coherent`, default) | Transducers fire simultaneously; pressure fields interfere. Modelled in a single k-Wave run with all transducer sources active. |
| **Async** (`transducer_coupling: async`) | Transducers fire in separate, non-overlapping windows; their pressure fields are independent. Each is simulated separately; intensities are summed for thermal analysis and the per-voxel pressure maximum is used for acoustic safety metrics. |

Use async mode when the firing timing of your system means the acoustic fields from different transducers never overlap in time (e.g., alternating duty cycles, sequential pulses). If the transducers fire simultaneously or with overlapping duty cycles, use the default coherent mode.

---

## Physical model

Because only one transducer is active at any instant, the two physically relevant quantities — acoustic safety and thermal heating — require different combination rules.

### Acoustic safety metrics (ISPPA, MI, peak pressure)

The instantaneous peak pressure at any voxel is the maximum over all transducer fields:

$$p_\mathrm{safety} = \max_i\bigl(p_{i,\mathrm{scaled}}\bigr)$$

where $p_{i,\mathrm{scaled}} = p_i \cdot \sqrt{I_{i,\mathrm{target}} / I_{i,\mathrm{baseline}}}$.

ISPPA, MI, and peak-pressure metrics reported in the output CSV and plots are derived from $p_\mathrm{safety}$ (`sensor_data.p_max_async` in the combined cache). The focal ISPPA of the combined exposure therefore equals the higher of the two individual ISPPAs, not their sum.

### Thermal heat deposition

Time-averaged heat deposition accumulates from all transducers. Assuming equal, non-overlapping duty cycles that together fill the protocol duty cycle, the combined thermal heat source is the incoherent intensity sum:

$$Q_\mathrm{total} \propto \sum_{i=1}^{N} I_i = \sum_{i=1}^{N} \frac{p_{i,\mathrm{scaled}}^2}{2\rho c}$$

This is converted to a thermal-equivalent pressure stored as `sensor_data.p_max_all` in the combined cache and passed directly to the bioheat solver:

$$p_\mathrm{thermal} = \sqrt{2 \cdot I_\mathrm{total} \cdot \rho c}$$

> **Limitation.** The thermal model does not yet explicitly represent the relative firing schedule of the two transducers (e.g., T2 offset by half a PRI relative to T1). The intensity-sum approximation is exact when both DCs are equal and non-overlapping; it becomes an over-estimate when $\mathrm{DC}_1 + \mathrm{DC}_2 < $ total protocol DC. Explicit interleaved bioheat modeling is planned as a future extension.

---

## Configuration

### Basic setup — two transducers, one thermal target each

```yaml
simulation:
  transducer_coupling: async   # activate async mode

transducer:
  - name: CTX500-A
    target_isppa_wcm2: 5.0     # target free-water ISPPA for this transducer
    freq_hz: 500000
    trans_pos: [x1, y1, z1]
    focus_pos: [fx1, fy1, fz1]
    annular:
      elem_amp: 100000
      # ... other annular fields

  - name: CTX500-B
    target_isppa_wcm2: 3.0
    freq_hz: 500000
    trans_pos: [x2, y2, z2]
    focus_pos: [fx2, fy2, fz2]
    annular:
      elem_amp: 100000
      # ... other annular fields

modules:
  run_water_baseline: 1        # required to record free-water ISPPA provenance
```

Both transducers must share the same segmentation (same `subject_id`, `simulation.medium`, and output directory). The simulation grid is determined by the union of both transducer positions.

### Multi-ISPPA sweep — varying intensity across sweep points

Setting `target_isppa_wcm2` to a vector on one or more transducers triggers a sweep. Each sweep point runs an independent combine+thermal pair.

All transducers must carry either a **scalar** target (held fixed across the sweep) or a **vector of the same length** as the sweep:

```yaml
transducer:
  - name: CTX500-A
    target_isppa_wcm2: [3.0, 5.0, 10.0]   # swept across 3 points
    # ...

  - name: CTX500-B
    target_isppa_wcm2: 2.0                  # fixed across all 3 sweep points
    # ...
```

This produces three combine+thermal pairs:
- Sweep 1: A → 3.0 W/cm², B → 2.0 W/cm²
- Sweep 2: A → 5.0 W/cm², B → 2.0 W/cm²
- Sweep 3: A → 10.0 W/cm², B → 2.0 W/cm²

---

## Pipeline stages

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

On HPC (SLURM/qsub):
- All N acoustic stages are submitted simultaneously (no inter-dependency).
- Each combine stage depends on all N acoustic jobs (`afterok`).
- Each thermal stage depends on its corresponding combine stage.

---

## Output files

For two transducers, single thermal target each:

```
cache/sub-001_layered_results_desc-t01.mat         ← acoustic cache, transducer 1
cache/sub-001_layered_results_desc-t02.mat         ← acoustic cache, transducer 2
cache/sub-001_layered_results_desc-asyncCombined.mat  ← combined intensity cache
sub-001_layered_desc-asyncThermal.csv              ← thermal analysis table
sub-001_layered_desc-asyncThermal_report.html      ← report
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

---

## Activation

Async mode is activated automatically by `prestus_pipeline_start` when both conditions are met:

```matlab
parameters.simulation.transducer_coupling = 'async';  % and
numel(parameters.transducer) > 1
```

No special entry point is needed:

```matlab
parameters = load_parameters('config_subject.yaml', subject_id);
prestus_pipeline_start(parameters);   % dispatches to async_transducer_pipeline
```

---

## Constraints and limitations

- **3-D grid required.** Axisymmetric (2-D) simulations are not supported because two transducers at different positions cannot share a common axis of symmetry.
- **Shared head model.** Both transducers must use the same segmented head grid (`subject_id`, `simulation.medium`, `io.dir_output`).
- **Linear regime.** The same validity conditions as single-transducer ISPPA scaling apply (see [doc_calibration.md](doc_calibration.md)). Each transducer's pressure field is scaled independently from its own free-water baseline.
- **`modules.run_water_baseline: 1` is required.** Without it, no `freefield_isppa_wcm2` provenance is recorded and pressure scaling is skipped (a warning is printed).
- **Async mode is mutually exclusive with uncertainty mode.** The two cannot currently be combined. Use coherent mode with uncertainty if skull-property uncertainty quantification is required.

---

## HPC options

```matlab
options.acoustic_timelimit  = '06:00:00';   % per acoustic job
options.combine_timelimit   = '00:30:00';   % per combine job (CPU-only, lightweight)
options.thermal_timelimit   = '04:00:00';   % per thermal job
options.acoustic_memorylimit = 60;           % GB
options.thermal_memorylimit  = 40;           % GB
options.acoustic_partition  = 'gpu';
options.combine_partition   = 'cpu';         % combine stage needs no GPU
```

---

## Resumption

Each stage checks for its sentinel output file before submitting or running. Re-running `prestus_pipeline_start` with the same config resumes from the last incomplete stage.

---

## Provenance

The combined acoustic cache records the following fields in `acoustic_provenance`:

| Field | Description |
|---|---|
| `combined_from` | Cell array of source cache file paths |
| `combination_mode` | `'async_dual_field'` — `p_max_all` holds the thermal intensity sum; `p_max_async` holds the per-voxel pressure maximum for acoustic safety |
| `combination_targets` | Per-transducer target ISPPAs used for scaling [W/cm²] |
| `combination_baselines` | Per-transducer free-water baselines at the simulated drive level [W/cm²] |
| `combination_time` | Timestamp of the combination step |

---

See also: [doc_calibration.md](doc_calibration.md), [doc_multi_isppa.md](doc_multi_isppa.md), [doc_modules.md](doc_modules.md), [doc_transducer.md](doc_transducer.md)
