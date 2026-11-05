# Multi-ISPPA Mode

!!! note "At a glance"
    - **What:** One acoustic simulation, scaled to several free-water target intensities, each with its own thermal simulation and report.
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

PRESTUS then runs one acoustic simulation, caches it, and runs the thermal analysis at each requested intensity, each with its own output files and HTML report, plus a summary report. The free-water baseline that the scaling needs runs automatically (`modules.run_water_baseline: auto`, the default); see [ISPPA targeting](doc_simulations-acoustic.md#isppa-targeting).

!!! note
    The older form `calibration.target_isppa_wcm2` is still accepted and is moved to `transducer(1).target_isppa_wcm2` at load time with a deprecation warning.

## Output files

Each thermal stage writes output under its own affix. For `target_isppa_wcm2: [10, 30]`:

```
sub-001_layered_desc-isppa010.csv
sub-001_layered_report_desc-isppa010.html
sub-001_layered_desc-isppa030.csv
sub-001_layered_report_desc-isppa030.html
sub-001_layered_desc-multiIsppa_report.html   ← summary
```

## Validity of post-hoc scaling

The acoustic simulation (k-Wave skull propagation) is the expensive step. In the linear acoustic regime the pressure field scales as $p \propto \sqrt{I}$, so one acoustic result serves every target. Linear scaling is valid when:

- the acoustic simulation was run in the **linear regime** (target amplitudes below the onset of nonlinear harmonic generation, roughly < 50–100 W/cm² free-water for typical TUS setups);
- the **spatial pressure pattern** is the same at all target intensities (guaranteed in the linear regime);
- the thermal simulation uses the **scaled pressure** as input (the pipeline scales `sensor_data.p_max_all` before `thermal_simulation` is called).

A warning is printed if the scale factor relative to the baseline exceeds 4× or falls below 0.25×. Large scale factors mean the baseline and target are far apart in amplitude; at very high targets, nonlinear effects may invalidate the scaled field.

## Resumption

Each stage checks for its sentinel output file before submitting (HPC) or running (MATLAB). Re-running `prestus_pipeline_start` with the same config resumes a partially completed run; completed stages are skipped. With `io.overwrite_files: always` (default), only the acoustic stage resumes from its cache; thermal and summary stages always re-run. `force` re-runs every stage. See [Overwrite behaviour](doc_outputs.md#overwrite-behaviour).

??? note "Details: pipeline stages and HPC options"
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

    On HPC (SLURM/qsub), stages 2…N are submitted in parallel with an `afterok` dependency on stage 1. The summary depends on all thermal jobs.

    Pass an `options` struct to `prestus_pipeline_start` to control per-stage resources:

    ```matlab
    options.acoustic_timelimit  = '08:00:00';   % wall time for stage 1
    options.thermal_timelimit   = '02:00:00';   % wall time per thermal job
    options.summary_timelimit   = '00:30:00';   % wall time for summary
    options.acoustic_memorylimit = 80;           % GB for stage 1
    options.thermal_memorylimit  = 40;           % GB per thermal job
    options.generate_summary_report = false;     % skip summary stage
    ```

## Relationship to uncertainty mode

Multi-ISPPA mode and uncertainty mode compose. When both `simulation.uncertainty = true` and `transducer.target_isppa_wcm2` is a vector, uncertainty is the outer mode and each variant runs its own multi-ISPPA pipeline:

```
Stage 1    preproc / source setup                      (shared)
Variants   default, Temp-H, Temp-L, Trans-H, Trans-L, each:  (after stage 1)
             acoustic + free-water baseline
             thermal: one job per target               (_desc-temphigh_desc-isppa030 …;
                                                        none for Trans-H/L)
             multi-ISPPA summary
Stage 7    uncertainty report: one per target          (after every variant job)
```

Affixes compose: the variant affix (`_desc-temphigh`) prefixes the per-target affix, yielding `_desc-temphigh_desc-isppa030`. Each target gets an uncertainty report comparing the five variants at that intensity, e.g. `sub-001_layered_desc-uncertainty_report_desc-isppa030.html` (with the base `io.output_affix` before the target suffix when one is set). No report aggregates across intensities. See [doc_uncertainty.md](doc_uncertainty.md#combining-with-a-multi-isppa-sweep-or-async-transducers).

## Relationship to async multi-transducer mode

When `simulation.transducer_coupling = 'async'` and multiple transducers are defined, each transducer carries its own `target_isppa_wcm2`. The async pipeline scales each transducer's pressure field independently before summing intensities incoherently. If any transducer has a vector target, one combine+thermal pair is run per sweep point. See [doc_multitransducer.md](doc_multitransducer.md).

See also: [doc_calibration.md](doc_calibration.md), [doc_uncertainty.md](doc_uncertainty.md), [doc_modules.md](doc_modules.md)
