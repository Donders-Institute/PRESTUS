# Async Multi-Transducer Mode

When transducers fire in separate, non-overlapping windows, their fields never interfere, so the coherent sum of [Multi-Transducer Modelling](doc_multitransducer.md) would misstate the exposure. Async mode simulates each transducer on its own, scales each field to its own ISPPA target, and combines them: acoustic safety metrics take the per-voxel maximum and heating uses the sum of intensities.

On this page:

- [Async Configuration](#async-configuration): switching on async mode, per-transducer targets and ISPPA sweeps.
- [Physical Model](#physical-model): how the fields are combined for safety metrics and heating.
- [Grid Geometry Mode](#grid-geometry-mode): one shared grid or one grid per transducer.
- [Output Files](#output-files): per-transducer, combined and joint outputs, resumption and HPC stages.
- [Constraints](#constraints): required grid, head model and regime.

!!! note "At a glance"
    - **When you need it:** Transducers fire with non-overlapping duty cycles or alternating pulses, or you want an ISPPA sweep across several transducers.
    - **Before you start:** A multi-transducer layered setup as in [Configuration](doc_multitransducer.md#configuration); a 3-D grid.

!!! warning "Experimental"
    For transducers that fire with non-overlapping duty cycles or alternating pulses, so their fields never superpose coherently. Configure it only with a clear view of the intended exposure and the linear-regime assumptions below.

# Async Configuration

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

??? note "Multi-ISPPA sweep"
    Every transducer carries either a scalar (held fixed) or a vector of the same length:

    ```yaml
    transducer:
      - name: CTX500-A
        target_isppa_wcm2: [3.0, 5.0, 10.0]   # swept across 3 points

      - name: CTX500-B
        target_isppa_wcm2: 2.0                  # fixed across all 3 sweep points
    ```

    This gives three combine + thermal pairs (A→3.0, 5.0, 10.0; B→2.0 each time).

# Physical Model

**Acoustic safety metrics (ISPPA, MI, peak pressure)** use the per-voxel maximum of the scaled fields. **Thermal heat deposition** uses the incoherent intensity sum, assuming equal, non-overlapping duty cycles that together fill the protocol duty cycle.

!!! warning "Limitation"
    The thermal model does not represent the firing schedule. The intensity sum is exact when both DCs are equal and non-overlapping; it over-estimates when $\mathrm{DC}_1 + \mathrm{DC}_2 < $ total protocol DC.

??? note "Equations"
    Acoustic safety metrics:

    $$p_\mathrm{safety} = \max_i\bigl(p_{i,\mathrm{scaled}}\bigr)$$

    where $p_{i,\mathrm{scaled}} = p_i \cdot \sqrt{I_{i,\mathrm{target}} / I_{i,\mathrm{baseline}}}$. This is `sensor_data.p_max_async` in the combined cache.

    Thermal heat deposition:

    $$Q_\mathrm{total} \propto \sum_{i=1}^{N} I_i = \sum_{i=1}^{N} \frac{p_{i,\mathrm{scaled}}^2}{2\rho c}$$

    This is `sensor_data.p_max_all`, passed to the bioheat solver as:

    $$p_\mathrm{thermal} = \sqrt{2 \cdot I_\mathrm{total} \cdot \rho c}$$

# Grid Geometry Mode

| `simulation.async_grid_mode` | Grid | Use when |
|---|---|---|
| `'shared'` (default) | One grid for all transducers, aligned and cropped as in [Grid alignment](doc_multitransducer.md#grid-alignment) and [Grid cropping](doc_multitransducer.md#grid-cropping). | Default. |
| `'per_transducer'` | Each transducer gets the grid a solo run would use. `combine_async_intensity` resamples every field onto transducer(1)'s grid via each run's `planimg` sidecar. | Per-transducer focal distance, ISPPA or focal volume must match a solo run, or the transducers are widely separated or steeply angled. |

??? note "Grid mode choice and accuracy"
    ```yaml
    simulation:
      transducer_coupling: async
      async_grid_mode: per_transducer   # default: 'shared'
    ```

    In `shared` mode an off-axis beam is voxelised at an angle, so its per-transducer results can differ from a solo run. The difference shrinks with resolution; at 2 mm (2.5 points per wavelength at 300 kHz) a beam 41° off axis lost 16% of its free-water ISPPA. In `per_transducer` mode the resampling (linear interpolation; parts outside transducer(1)'s grid are dropped) costs accuracy instead. `per_transducer` scales better with N: each grid stays sized to its own transducer, at the cost of N preproc jobs.

# Output Files

Two transducers, one target each:

```
cache/sub-001_layered_results_trans-01.mat
cache/sub-001_layered_results_trans-02.mat
cache/sub-001_layered_results_desc-asyncCombined.mat
sub-001_layered_desc-asyncJoint.csv
sub-001_layered_report_desc-asyncJoint.html
```

??? note "Sweep and affix naming"
    A 3-point sweep writes `..._desc-asyncCombinedS01.mat` to `S03` and `sub-001_layered_desc-asyncJointS01.csv` to `S03`.

    With `io.output_affix` (e.g. `run1`), every affix starts with it: `_desc-run1_trans-01`, `_desc-run1_desc-asyncJoint`. The per-transducer runs also write their own CSV and acoustic NIfTIs at the simulated drive level; the Nexus Results tab shows them under the joint result.

**Resumption.** Re-running `prestus_pipeline_start` with the same config resumes from the last incomplete stage. With `io.overwrite_files: always` (default), preproc, acoustic and combine resume from their caches; the joint and water-validation stages always re-run. `force` re-runs every stage. See [Overwrite behaviour](doc_outputs_reuse.md#overwrite-behaviour).

??? note "Pipeline stages, HPC options and provenance"
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

# Constraints

- **3-D grid required.** Axisymmetric (2-D) simulations are not supported.
- **Shared head model.** All transducers use the same segmentation (same `subject_id`, `simulation.medium` and output directory).
- **Linear regime.** The validity conditions of single-transducer ISPPA scaling apply (see [doc_calibration.md](doc_calibration.md)).
- **Uncertainty mode.** With `simulation.uncertainty: true`, each variant runs its own async pipeline (preproc included), and one uncertainty report is written per sweep point. See [doc_uncertainty.md](doc_uncertainty.md#combining-with-multi-isppa-or-async-modes).
