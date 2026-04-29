# Multi-ISPPA mode

When `calibration.target_isppa_wcm2` is a vector with more than one value, PRESTUS automatically enters **multi-ISPPA mode**. A single acoustic simulation is run and cached; thermal analysis is then performed in parallel at each requested free-water intensity, with each job producing its own output files and HTML report.

---

## Motivation

The acoustic simulation (k-Wave skull propagation) is the computationally expensive step. Repeating it for each target intensity is unnecessary because, in the linear acoustic regime, the pressure field scales as $p \propto \sqrt{I}$. Multi-ISPPA mode exploits this: the acoustic cache is reused for every thermal job, and the pressure field is scaled post-hoc before thermal analysis.

---

## Pipeline stages

```
[Stage 1]   Acoustic simulation + water baseline
            Runs once at the configured elem_amp.
            Saves sensor_data + acoustic_provenance.freefield_isppa_wcm2.
            No thermal simulation.

[Stage 2…N] Thermal analysis — one job per target_isppa_wcm2
            Loads cached acoustic results.
            Scales p_max_all by sqrt(target / baseline).
            Runs thermal simulation and analysis.
            Writes output files with affix _isppaNNN (NNN = integer W/cm²).
            Generates per-target HTML report.

[Stage N+1] Summary report (optional, default: on)
            Aggregates results across all target intensities.
```

On HPC (SLURM/qsub), stages 2…N are submitted in parallel with an `afterok` dependency on stage 1. The summary depends on all thermal jobs.

---

## Configuration

```yaml
calibration:
  target_isppa_wcm2: [10, 20, 30, 50]   # W/cm² — triggers multi-ISPPA mode
```

**Important:** `modules.run_water_baseline` defaults to `0` — you must explicitly enable it so that the provenance required for scaling is recorded:

```yaml
modules:
  run_water_baseline: 1
```

---

## Output files

Each thermal stage writes output under its own affix. For `target_isppa_wcm2: [10, 30]`:

```
sub-001_layered_output_table_isppa010.csv
sub-001_layered_report_isppa010.html
sub-001_layered_output_table_isppa030.csv
sub-001_layered_report_isppa030.html
sub-001_layered_multi_isppa_report.html   ← summary
```

---

## Validity of post-hoc scaling

Linear scaling of the pressure field is valid when:

- The acoustic simulation was run in the **linear regime** (target amplitudes are below the onset of nonlinear harmonic generation, roughly < 50–100 W/cm² free-water for typical TUS setups)
- The **spatial pressure pattern** is the same at all target intensities (guaranteed in the linear regime)
- Downstream thermal simulation uses the **scaled pressure** as input (ensured by the pipeline — scaling is applied to `sensor_data.p_max_all` before `thermal_simulation` is called)

A warning is printed if the scale factor relative to the baseline exceeds 4× or falls below 0.25×. Large scale factors indicate the baseline and target are far apart in amplitude; at very high targets, nonlinear effects may invalidate the scaled field.

---

## HPC options

Pass an `options` struct to `prestus_pipeline_start` to control per-stage resources:

```matlab
options.acoustic_timelimit  = '08:00:00';   % wall time for stage 1
options.thermal_timelimit   = '02:00:00';   % wall time per thermal job
options.summary_timelimit   = '00:30:00';   % wall time for summary
options.acoustic_memorylimit = 80;           % GB for stage 1
options.thermal_memorylimit  = 40;           % GB per thermal job
options.generate_summary_report = false;     % skip summary stage
```

---

## Resumption

Each stage checks for its sentinel output file before submitting (HPC) or running (MATLAB). A partially completed multi-ISPPA run can be resumed by re-running `prestus_pipeline_start` with the same config — completed stages are skipped automatically.

---

## Relationship to uncertainty mode

Multi-ISPPA mode and uncertainty mode are composable. When both `simulation.uncertainty = true` and `calibration.target_isppa_wcm2` is a vector, the two modes nest automatically:

```
Stage 1   preproc / source setup         (shared)
Stage 2   acoustic — default variant     ┐
Stage 3   acoustic — liberal variant     ├ each includes water baseline
Stage 4   acoustic — conservative variant┘
Stages 5… thermal — one job per variant × target  (_liberal_isppa030, _conservative_isppa010 …)
Final     report
```

This works because `uncertainty_pipeline` passes the full `parameters` struct (including the vector `target_isppa_wcm2`) into each variant's `prestus_pipeline` call, which then dispatches to `multi_isppa_pipeline`. Affixes compose naturally: the uncertainty variant affix (`_liberal`) is already set in `io.output_affix` when `multi_isppa_pipeline` builds its per-target affixes, yielding `_liberal_isppa030` etc.

For each target intensity, an uncertainty report is generated aggregating the three variants at that intensity (e.g. `sub-001_layered_uncertainty_report_isppa030.html`). No report aggregating across intensities is produced.

See also: [doc_calibration.md](doc_calibration.md), [doc_uncertainty.md](doc_uncertainty.md), [doc_modules.md](doc_modules.md)
