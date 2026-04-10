# Uncertainty quantification

> **⚠ Experimental feature.** The uncertainty quantification workflow is under active development. The liberal and conservative medium property ranges shipped with PRESTUS are informed by published literature but should not be interpreted as formally validated bounds. They reflect plausible variation in tissue acoustic properties and are intended to support sensitivity analysis, not to replace subject-specific measurements or regulatory assessment. The ranges may not capture all sources of simulation uncertainty (e.g., segmentation errors, transducer positioning uncertainty, skull heterogeneity). Users are strongly encouraged to review and adapt the medium property ranges to their specific study population and equipment before drawing safety conclusions.

TUS simulations of wave propagation through the skull involve acoustic and thermal medium properties (e.g., speed of sound, density, attenuation) that are not precisely known for any individual subject. Published values vary considerably across studies. This uncertainty has direct implications for safety assessments: the true in-situ intensity and heating at a given location may differ from the simulated value depending on which tissue properties are assumed.

PRESTUS provides an **uncertainty mode** that makes these assumptions explicit. Rather than running a single simulation with one set of tissue properties, the uncertainty mode runs three parallel simulations — *default*, *liberal*, and *conservative* — each with a different but internally consistent set of medium properties. A dedicated postprocessing step then compares outcomes and produces an uncertainty report that summarises the range of plausible values for each safety metric.

---

## Conceptual framework

The three simulation variants are defined as:

| Variant | Interpretation | Expected effect |
|---|---|---|
| **Default** | Best-estimate medium properties | Reference simulation |
| **Liberal** | Acoustic tissue parameter combination that tends to produce higher acoustic intensity at the target and reduced heating in the skull (e.g., lower skull attenuation, lower skull absorption fraction) | Indicative upper range of in-situ exposure |
| **Conservative** | Acoustic tissue parameter combination that tends to produce lower acoustic intensity at the target and increased heating in the skull (e.g., higher skull attenuation, higher skull absorption fraction) | Indicative lower range of in-situ exposure |

The liberal/conservative distinction follows the logic that different parameter choices affect different stages of wave propagation in opposite directions. For example, higher skull attenuation reduces intracranial intensity (beneficial) but increases skull heating (adverse). The liberal and conservative parameter sets are therefore not simply "high" and "low" across all properties — they are tailored to yield the highest and lowest plausible values for the metrics of interest (intracranial intensity and heating).

> **Limitations.** The parameter value ranges used for the liberal and conservative variants remain under investigation and should not be interpreted as absolute bounds on simulation outcomes. Parameter interactions can be non-linear — both amongst tissue properties and with subject-specific head anatomy such as skull thickness and morphology. Additionally, analysis choices such as the degree of image smoothing applied to the segmentation independently affect intensity and thermal estimates. The liberal and conservative variants are therefore best understood as illustrative scenarios that bracket a plausible range, not as guaranteed worst- or best-case predictions.

---

## Medium property ranges

The following table summarises the default, liberal, and conservative medium acoustic properties currently implemented in PRESTUS. Liberal values are chosen to maximise acoustic exposure and heating at the target; conservative values minimise them.

### Brain

| Property | Conservative | Default | Liberal | Unit |
|---|---|---|---|---|
| Sound speed | 1520 | 1546 | 1568 | m/s |
| Density | 1046 | 1046 | 1046 | kg/m³ |
| Attenuation coeff. | 0.94 | 0.59 | 0.36 | dB/(cm·MHz) |
| Attenuation power | 1.0 | 1.2 | 1.0 | — |
| Thermal conductivity | 0.48 | 0.51 | 0.55 | W/(m·K) |
| Specific heat | 3583 | 3630 | 3696 | J/(kg·K) |
| Perfusion | 196 | 559 | 809 | mL/min/kg |
| Absorption fraction | 0.80 | 1.00 | 1.00 | — |

### Skull

| Property | Conservative | Default | Liberal | Unit |
|---|---|---|---|---|
| Sound speed | 2822 | 2800 | 2257 | m/s |
| Density | 1908 | 1850 | 1178 | kg/m³ |
| Attenuation coeff. | 24.4 | 13.3 | 6.9 | dB/(cm·MHz) |
| Attenuation power | 1.0 | 1.0 | 1.0 | — |
| Thermal conductivity | 0.31 | 0.32 | 0.32 | W/(m·K) |
| Specific heat | 1313 | 1793 | 2626 | J/(kg·K) |
| Perfusion | 10 | 20 | 30 | mL/min/kg |
| Absorption fraction | 1.00 | 0.28 | 0.16 | — |

> **Rationale for skull properties.** The conservative variant uses high attenuation and high acoustic impedance (high density × sound speed), which reflects a denser, more calcified skull: more acoustic energy is lost before reaching the brain, so intracranial intensity is lower. Conversely, the liberal variant uses a skull with low impedance and low attenuation, which allows more energy through to the brain. The absorption fraction modulates how much of the attenuated energy is converted to heat: a lower fraction means less skull heating per unit attenuation. In practice, the liberal variant thus produces higher intracranial intensity *and* allows more limited skull heating estimates relative to the default.

> **Multi-compartment skull models.** When a segmentation includes separate skull sub-compartments (e.g., `skull_cortical` and `skull_trabecular`), the same liberal and conservative values are currently applied identically to all skull compartments. This is a known simplification — in reality, cortical and trabecular bone have distinct acoustic properties. Users modelling detailed skull microstructure should consider customising the uncertainty configs (see [Customising medium property ranges](#customising-medium-property-ranges)) to assign compartment-specific ranges.

---

## Workflow

The uncertainty mode consists of five stages, of which stages 2–4 can run in parallel:

```
Stage 1: Preprocessing & source setup  (serial; must complete before stages 2–4)
    ↓
Stage 2: Default simulation             ┐
Stage 3: Liberal simulation             ├── parallel acoustic + thermal jobs
Stage 4: Conservative simulation        ┘
    ↓
Stage 5: Uncertainty report generation  (serial; runs once all three output files are present)
```

All five stages are managed automatically when `parameters.simulation.uncertainty = true` is set before calling `prestus_pipeline`. This is the same entry point used for standard simulations — no separate function call is required.

### Basic usage

```matlab
prestus_path = '/path/to/PRESTUS';
addpath(fullfile(prestus_path, 'functions', 'helper'));
safe_addpath(fullfile(prestus_path, 'functions'));
safe_addpath(fullfile(prestus_path, 'toolboxes'));

parameters = load_parameters('config_study.yaml');
parameters.subject_id             = 1;
parameters.simulation.medium      = 'layered';
parameters.simulation.uncertainty = true;    % activates uncertainty mode
parameters.platform               = 'auto';  % 'matlab', 'slurm', or 'qsub'
parameters.path.sim               = '/absolute/path/to/sim_outputs';  % must be set

prestus_pipeline_start(parameters);
```

> **`parameters.path.sim` must be an absolute path.** The default config sets `path.sim = ''`. If it is empty when the pipeline is called, `get_output_dir` will throw an error. Set it in your study config or explicitly in the script as shown above.

`prestus_pipeline_start` detects `simulation.uncertainty = true` on HPC and redirects to `uncertainty_pipeline` before submitting a single job. On MATLAB, `prestus_pipeline` handles the redirect after `path_log_setup` has run. The per-variant parameter structs built internally have `simulation.uncertainty` cleared, so there is no recursion.

On a local MATLAB session (`platform = 'matlab'`) the five stages run sequentially. On HPC (`platform = 'slurm'` or `'qsub'`) all five jobs are submitted immediately: stages 2–4 with a dependency on stage 1, and stage 5 with a dependency on all three simulation jobs. No MATLAB session is held open during execution.

### Customising options

Uncertainty-specific options can be passed as the second argument to `prestus_pipeline`:

```matlab
options.affixes.default       = '';
options.affixes.liberal       = '_liberal';
options.affixes.conservative  = '_conservative';

% Override built-in medium configs with study-specific files:
options.liberal_config      = '/path/to/my_medium_liberal.yaml';
options.conservative_config = '/path/to/my_medium_conservative.yaml';

% HPC wall times:
options.stage1_timelimit  = '00:30:00';
options.sim_timelimit     = '04:00:00';
options.report_timelimit  = '00:30:00';

prestus_pipeline(parameters, options);
```

You can also call `uncertainty_pipeline` directly if you need to bypass `prestus_pipeline`:

```matlab
uncertainty_pipeline(parameters, options);
```

### What happens inside

`uncertainty_pipeline` builds four parameter structs from the base parameters and dispatches them:

| Stage | Modules enabled | `io.output_affix` | Medium override |
|---|---|---|---|
| 1 | source setup only | `''` | none |
| 2 | acoustic + thermal | `''` | none (default) |
| 3 | acoustic + thermal | `'_liberal'` | `config_medium_liberal.yaml` |
| 4 | acoustic + thermal | `'_conservative'` | `config_medium_conservative.yaml` |
| 5 | uncertainty report | `''` | none |

The medium override YAMLs are merged on top of the base parameters using `MergeStruct`, so only the properties defined in the override file are changed — all other settings (transducer, grid, paths) are inherited from the base config.

> **Resuming a partial run.** On HPC, each stage is skipped automatically if its sentinel output already exists (see [Resuming a partial run](#resuming-a-partial-run) below). Only the missing stages are resubmitted, and scheduler dependencies are set only for jobs that were actually submitted.

The built-in medium override configs live in `configs/uncertainty/`:

- `config_medium_liberal.yaml` — low-impedance, low-attenuation skull; higher brain intensity
- `config_medium_conservative.yaml` — high-impedance, high-attenuation skull; lower brain intensity

---

## Output files

Each simulation variant writes its outputs into the shared `path.sim` output folder (or a subject subfolder), identified by `io.output_affix`:

| File pattern | Description | Kept when `save_matrices = 0` |
|---|---|---|
| `sub-NNN_layered_output_table<affix>.csv` | Acoustic and thermal metrics table | yes |
| `sub-NNN_layered_report<affix>.html` | Per-variant simulation report | yes |
| `sub-NNN_layered_intensity_<dim><affix>.png` | Intensity overlay image per slice dimension (x, y, z) | yes |
| `sub-NNN_layered_maxT_<dim><affix>.png` | Max temperature overlay image per slice dimension (x, y, z) | yes |
| `sub-NNN_layered_thermal_max<affix>.png` | Temperature vs. time plot | yes |
| `sub-NNN_layered_parameters<affix>.mat` | Parameter struct for the variant | yes |
| `sub-NNN_layered_heating_res<affix>.mat` | Heating simulation matrices (large) | **removed** |
| `debug/sub-NNN_after_rotating_and_scaling.mat` | Intermediate head model (debug) | **removed** |
| `debug/sub-NNN_layered_after_cropping_and_smoothing.mat` | Intermediate skull grid (debug) | **removed** |

> **Automatic cleanup.** When `parameters.io.save_matrices = 0`, the pipeline deletes the files marked **removed** above after the uncertainty report (stage 5) has been written. This is handled by `cleanup_uncertainty_intermediates`, called automatically on both the MATLAB and HPC platforms. Cleanup only runs after stage 5 completes, so a failed run can be resumed beforehand.

The uncertainty report is written as:

| File | Description |
|---|---|
| `sub-NNN_layered_uncertainty_report.html` | Unified uncertainty HTML report |

![PRESTUS GUI](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/prestus_uncertainty.png)

---

## Uncertainty report contents

The uncertainty report (`generate_uncertainty_report.m`) is a self-contained HTML file that includes:

- **Safety dashboard** — Color-coded cards for all ITRUSST consensus safety metrics (MI, CEM43, temperature rise, absolute temperature). Each card shows the **default** value with the liberal–conservative range annotated. The bar spans a metric-specific minimum (37 °C for absolute temperatures, 0 elsewhere) to `max(max_sim, limit)`. The colored fill extends to `min(max_sim, limit)`; three ticks mark the liberal (blue), default (black), and conservative (brown) values. Color is determined by the conservative value.
- **Acoustic results** — Range table showing liberal / default / conservative values for key acoustic metrics (ISPPA, focal distance, MI, half-maximum volume).
- **Thermal results** — Range table for thermal metrics (maxT, temperature rise, CEM43 per tissue). Side-by-side temperature-vs-time images (`thermal_max`) across liberal, default, and conservative variants.
- **Additional maps** — Side-by-side maximum temperature maps across variants.
- **Raw data tables** — Full CSV output tables for each variant.

---

## Resuming a partial run

The uncertainty pipeline supports automatic resumption on both platforms.

**HPC (slurm/qsub).** Before submitting each job, the pipeline checks for a sentinel file. If the sentinel exists, the job is skipped entirely and no scheduler dependency is created for it:

| Stage | Sentinel file | Location |
|---|---|---|
| 1 — Preprocessing | `sub-NNN_<medium>_after_cropping_and_smoothing.mat` | `<output_dir>/debug/` |
| 2 — Default sim | `sub-NNN_<medium>_output_table.csv` | `<output_dir>/` |
| 3 — Liberal sim | `sub-NNN_<medium>_output_table_liberal.csv` | `<output_dir>/` |
| 4 — Conservative sim | `sub-NNN_<medium>_output_table_conservative.csv` | `<output_dir>/` |
| 5 — Report | `sub-NNN_<medium>_uncertainty_report.html` | `<output_dir>/` |

To rerun a stage, delete its sentinel file before resubmitting.

**MATLAB platform.** There is no stage-level skip logic. Instead, resumability relies on `io.overwrite_files = 'never'` being set for all variant stages: each individual file write inside the pipeline is skipped if the file already exists. This means a restarted MATLAB run will fast-forward through completed steps file-by-file rather than stage-by-stage.

---

## HPC usage

When `platform` is `'slurm'` or `'qsub'`, `uncertainty_pipeline` submits all five jobs immediately using scheduler dependencies (`--dependency=afterok` for SLURM). Stages 2–4 depend on stage 1; stage 5 depends on all three simulation jobs. If stage 1 fails, all downstream jobs are automatically cancelled (`--kill-on-invalid-dep=yes`). No MATLAB session is held open during execution.

Stages 1 and 5 (preprocessing and report) do not require a GPU and are submitted to the scheduler's default CPU queue. Stages 2–4 inherit the GPU resources from the base config. To override the queue for CPU stages use `options.stage1_partition = 'gpu'` and `options.report_partition = 'gpu'`.

This requires that `hpc_submit_job` supports the `parameters.hpc.depend_job_id` (single) and `parameters.hpc.depend_job_ids` (array) fields — see [doc_hpc.md](doc_hpc.md).

---

## Customising medium property ranges

The built-in overrides are in `configs/uncertainty/`. To use study-specific ranges, copy and edit either file:

```yaml
# my_medium_liberal.yaml — only override what you want to change
medium_properties:
  skull:
    alpha_coeff: 5.0
    absorption_fraction: 0.12
```

Pass it via options:

```matlab
options.liberal_config = '/path/to/my_medium_liberal.yaml';
prestus_pipeline(parameters, options);
```

Only the properties present in the override file are changed; all other settings are inherited from the base config.

---

## References

The liberal and conservative parameter ranges are informed by published measurements and systematic reviews of transcranial ultrasound tissue properties:

- Aubry, J.-F., et al. ITRUSST consensus on safety of transcranial focused ultrasound. *Brain Stimulation* (2025). https://doi.org/10.1016/j.brs.2025.10.007
- Deffieux, T., & Konofagou, E. E. Numerical study of a simple transcranial focused ultrasound system applied to blood–brain barrier opening. *IEEE Trans Ultrason Ferroelectr Freq Control* **57**, 2637–2653 (2010).
- Pinton, G., et al. Attenuation, scattering, and absorption of ultrasound in the skull bone. *Med Phys* **39**, 299–307 (2012). https://doi.org/10.1118/1.3668316
- Marsac, L., et al. Ex vivo optimisation of a heterogeneous speed of sound model of the human skull for non-invasive transcranial focused ultrasound at 1 MHz. *Int J Hyperthermia* **33**, 635–645 (2017).
