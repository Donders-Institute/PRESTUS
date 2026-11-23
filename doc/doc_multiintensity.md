# Multi-I<sub>SPPA</sub> Mode

Heating and safety metrics depend on the intensity a study applies, and choosing one often means comparing several candidates. In the linear regime the pressure field scales with the square root of intensity, so this mode runs the acoustic simulation once and scales it to each target before a separate thermal simulation. The scaling holds only below the onset of nonlinear effects; at very high targets the scaled field may not be valid.

On this page:

- [Configuration](#configuration): giving several target intensities.
- [Output Files](#output-files): per-target files and the summary report.
- [Validity of Post-Hoc Scaling](#validity-of-post-hoc-scaling): when one acoustic result can serve every target.
- [Resuming a Partial Run](#resuming-a-partial-run): re-running a partial sweep, stages, and HPC options.
- [Combining with Uncertainty Mode](#combining-with-uncertainty-mode): running a sweep for each uncertainty variant.
- [Combining with Async Transducers](#combining-with-async-transducers): per-transducer targets and sweeps.

!!! note "At a glance"
    - **When you need it:** To compare heating and safety metrics across candidate intensities without repeating the acoustic simulation.
    - **Before you start:** A study config whose single-intensity run works; targets within the linear regime (below roughly 50–100 W/cm² free-water).

## Configuration

Give `transducer.target_isppa_wcm2` more than one value:

```yaml
transducer:
  target_isppa_wcm2: [10, 20, 30, 50]   # W/cm²; triggers multi-ISPPA mode
  elem_p_amp: ...
  # ... other transducer fields
```

PRESTUS then runs one acoustic simulation, caches it, and runs the thermal analysis at each requested intensity, each with its own output files and HTML report, plus a summary report. The [free-water baseline](doc_glossary.md#free-water-baseline) that the scaling needs runs automatically (`modules.run_water_baseline: auto`, the default); see [I<sub>SPPA</sub> targeting](doc_simulations-acoustic.md#isppa-targeting). In the Nexus, *Sweep…* next to the target [I<sub>SPPA</sub>](doc_glossary.md#isppa) (Transducer tab) edits the list.

!!! note
    The older form `calibration.target_isppa_wcm2` is still accepted and is moved to `transducer(1).target_isppa_wcm2` at load time with a deprecation warning.

## Output Files

Each thermal stage writes output under its own affix. For `target_isppa_wcm2: [10, 30]`:

```
sub-001_layered_desc-isppa010.csv
sub-001_layered_report_desc-isppa010.html
sub-001_layered_desc-isppa030.csv
sub-001_layered_report_desc-isppa030.html
sub-001_layered_desc-multiIsppa_report.html   ← summary
```

## Validity of Post-Hoc Scaling

The acoustic simulation ([k-Wave](https://github.com/ucl-bug/k-wave) skull propagation) is the expensive step. In the linear acoustic regime the pressure field scales as $p \propto \sqrt{I}$, so one acoustic result serves every target. Linear scaling is valid when:

- the acoustic simulation was run in the **linear regime** (target amplitudes below the onset of nonlinear harmonic generation, roughly < 50–100 W/cm² free-water for typical TUS setups);
- the **spatial pressure pattern** is the same at all target intensities (guaranteed in the linear regime);
- the thermal simulation uses the **scaled pressure** as input (the pipeline scales `sensor_data.p_max_all` before `thermal_simulation` is called).

A warning is printed if the scale factor relative to the baseline exceeds 4× or falls below 0.25×. Large scale factors mean the baseline and target are far apart in amplitude; at very high targets, nonlinear effects may invalidate the scaled field.

## Resuming a Partial Run

Each stage checks for its sentinel output file before submitting (HPC) or running ([MATLAB](https://www.mathworks.com/products/matlab.html)). Re-running `prestus_pipeline_start` with the same config resumes a partially completed run; completed stages are skipped. With `io.overwrite_files: always` (default), only the acoustic stage resumes from its cache; thermal and summary stages always re-run. `force` re-runs every stage. See [Overwrite behaviour](doc_outputs_reuse.md#overwrite-behaviour).

??? note "Pipeline stages and HPC options"
    ```
    [Stage 1]   Acoustic simulation + water baseline
                Runs once at the configured elem_p_amp.
                Saves sensor_data + acoustic_provenance.freefield_isppa_wcm2.
                No thermal simulation.

    [Stage 2…N] Thermal analysis: one job per target_isppa_wcm2
                Loads cached acoustic results.
                Scales p_max_all by sqrt(target / baseline).
                Runs thermal simulation and analysis.
                Writes output files with affix _desc-isppaNNN (NNN = integer W/cm²).
                Generates per-target HTML report.

    [Stage N+1] Summary report (optional, default: on)
                Aggregates results across all target intensities.
    ```

    On HPC ([SLURM](https://slurm.schedmd.com/)/qsub), stages 2…N are submitted in parallel with an `afterok` dependency on stage 1. The summary depends on all thermal jobs.

    Pass an `options` struct to `prestus_pipeline_start` to control per-stage resources:

    ```matlab
    options.acoustic_timelimit  = '08:00:00';   % wall time for stage 1
    options.thermal_timelimit   = '02:00:00';   % wall time per thermal job
    options.summary_timelimit   = '00:30:00';   % wall time for summary
    options.acoustic_memorylimit = 80;           % GB for stage 1
    options.thermal_memorylimit  = 40;           % GB per thermal job
    options.generate_summary_report = false;     % skip summary stage
    ```

## Combining with Uncertainty Mode

Multi-I<sub>SPPA</sub> mode and uncertainty mode compose. When both `simulation.uncertainty = true` and `transducer.target_isppa_wcm2` is a vector, uncertainty is the outer mode and each variant runs its own multi-I<sub>SPPA</sub> pipeline:

```
Stage 1    preproc / source setup                      (shared)
Variants   default, Temp-H, Temp-L, Trans-H, Trans-L, each:  (after stage 1)
             acoustic + free-water baseline
             thermal: one job per target               (_desc-temphigh_desc-isppa030 …;
                                                        none for Trans-H/L)
             multi-ISPPA summary
Stage 7    uncertainty report: one per target          (after every variant job)
```

Affixes compose: the variant affix (`_desc-temphigh`) prefixes the per-target affix, yielding `_desc-temphigh_desc-isppa030`. Each target gets an uncertainty report comparing the five variants at that intensity, e.g. `sub-001_layered_desc-uncertainty_report_desc-isppa030.html` (with the base `io.output_affix` before the target suffix when one is set). No report aggregates across intensities. See [Combining with Multi-I<sub>SPPA</sub> or Async Modes](doc_uncertainty.md#combining-with-multi-isppa-or-async-modes).

## Combining with Async Transducers

When `simulation.transducer_coupling = 'async'` and multiple transducers are defined, each transducer carries its own `target_isppa_wcm2`. The async pipeline scales each transducer's pressure field independently before summing intensities incoherently. If any transducer has a vector target, one combine+thermal pair is run per sweep point. See [Multi-Transducer Modelling](doc_multitransducer.md).

See also: [Transducer Calibration](doc_calibration.md), [Uncertainty Quantification](doc_uncertainty.md), [Pipeline Modules](doc_modules.md)
