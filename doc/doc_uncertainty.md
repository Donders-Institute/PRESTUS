# Uncertainty Quantification

The acoustic and thermal properties of skull and brain are not precisely known for an individual subject, and published values vary across studies, so the simulated in-situ intensity and heating depend on which values are assumed. Uncertainty mode makes this dependence explicit: it runs the default medium and four variants built from the [ITRUSST](https://itrusst.org/) tissue property bounds and reports the range of each safety metric.

On this page:

- [Basic Usage](#basic-usage): switching uncertainty mode on.
- [Simulation Variants](#simulation-variants): the five variants and why two pairs of bounds are needed.
- [Medium Property Ranges](#medium-property-ranges): the property values of each variant.
- [Workflow Stages](#workflow-stages): the seven stages and the options that customise them.
- [Customising Medium Property Ranges](#customising-medium-property-ranges): editing single bound values.
- [Partial-FOV Bound Variants](#partial-fov-bound-variants): running the bound variants in a cropped acoustic grid.
- [Output Files](#output-files): what each variant writes and what is cleaned up.
- [Uncertainty Report Contents](#uncertainty-report-contents): how the report presents the ranges.
- [Combining with Multi-I<sub>SPPA</sub> or Async Modes](#combining-with-multi-isppa-or-async-modes): uncertainty as the outer mode.
- [Resuming a Partial Run](#resuming-a-partial-run): skipping completed stages.
- [HPC Usage](#hpc-usage): job submission, dependencies, and queues.
- [Pseudo-CT Runs](#pseudo-ct-runs): how the bounds relate to a pseudo-CT skull.

!!! note "At a glance"
    - **When you need it:** When a safety assessment should state how much in-situ intensity and heating depend on the assumed tissue properties.
    - **Before you start:** A working study config; `simulation.uncertainty: true` and, on a cluster, a GPU budget for five runs.

!!! warning "Experimental feature"
    The uncertainty quantification workflow is under active development. The variant medium properties shipped with PRESTUS are the lower and upper bounds of the ITRUSST tissue property table. They are not formally validated bounds, and they do not capture all sources of simulation uncertainty. Each variant combines the bounds chosen one property at a time; parameter interactions can be non-linear, both amongst tissue properties and with subject-specific head anatomy such as skull thickness and morphology. Analysis choices such as the degree of image smoothing applied to the segmentation also affect intensity and thermal estimates. Treat the variants as scenarios that bracket a plausible range, not as guaranteed worst- or best-case predictions.

![Uncertainty report: safety metrics of the default variant, with the range across the Temp and Trans variants against each limit](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/prestus_uncertainty.png)

## Basic Usage

Set `simulation.uncertainty = true` (in the Nexus: the *Uncertainty analysis* switch on the Medium tab, *Uncertainty* sub-tab; the variant values are under *Values per variant* below the bounds plot) and call the usual entry point:

```matlab
prestus_path = '/path/to/PRESTUS';
addpath(fullfile(prestus_path, 'functions', 'helper'));
safe_addpath(fullfile(prestus_path, 'functions'));
safe_addpath(fullfile(prestus_path, 'external'));

parameters = load_parameters('config_study.yaml');
parameters.subject_id             = 1;
parameters.simulation.medium      = 'layered';
parameters.simulation.uncertainty = true;    % activates uncertainty mode
parameters.platform               = 'auto';  % 'matlab', 'slurm', or 'qsub'
parameters.path.sim               = '/absolute/path/to/sim_outputs';  % must be set

prestus_pipeline_start(parameters);
```

On a local [MATLAB](https://www.mathworks.com/products/matlab.html) session (`platform = 'matlab'`) the stages run one after another. On HPC (`platform = 'slurm'` or `'qsub'`) all jobs are submitted immediately with scheduler dependencies, and no MATLAB session is held open ([HPC usage](#hpc-usage)).

Uncertainty mode needs a tissue medium (`layered` or `phantom`): the variants change tissue properties only, so with `simulation.medium = 'water'` the pipeline stops with an error.

Uncertainty mode can be combined with sequential simulations; see [Sequential: Combining with uncertainty mode](doc_sequential.md#combining-with-uncertainty-mode).

## Simulation Variants

The five simulation variants are:

| Variant | Label | Medium properties | Runs heating | Used for |
|---|---|---|---|---|
| `default` | Default | PRESTUS default medium | as `modules.run_heating_sims` | Reference simulation |
| `temp_high` | Temp-H | per property, the bound that maximises brain heating | yes | Upper end of thermal metrics |
| `temp_low` | Temp-L | per property, the bound that minimises brain heating | yes | Lower end of thermal metrics |
| `trans_high` | Trans-H | per acoustic property, the bound that maximises brain intensity | no | Upper end of acoustic metrics |
| `trans_low` | Trans-L | per acoustic property, the bound that minimises brain intensity | no | Lower end of acoustic metrics |

The bounds come from the ITRUSST one-at-a-time sensitivity analysis, which varied each tissue property between its lower and upper bound with all other properties at reference. Two pairs are needed because brain sound speed, density and attenuation and skull sound speed push heating and intensity in opposite directions: higher brain attenuation, for example, raises heating in the brain but lowers the intensity that reaches deeper targets. A single high/low pair would therefore bound neither metric. The Trans variants vary acoustic properties only and run without heating.

## Medium Property Ranges

The tables list the default medium (`config/config_default.yaml`) and the four variant files in `config/uncertainty/`. Attenuation and absorption coefficients are in dB/(cm·MHz^y); the variant values are the ITRUSST 500 kHz values converted with y = 1 (`alpha_coeff` = α_500 kHz / 0.5). "–" marks a property the variant file does not list: the Trans variants run no heating, so their thermal properties are not used.

??? note "Brain and skull property tables"
    **Brain**

    | Property | Default | Temp-H | Temp-L | Trans-H | Trans-L | Unit |
    |---|---|---|---|---|---|---|
    | Sound speed | 1546 | 1573 | 1520 | 1520 | 1573 | m/s |
    | Density | 1046 | 1040 | 1050 | 1050 | 1040 | kg/m³ |
    | Attenuation coeff. | 0.59 | 0.984 | 0.308 | 0.308 | 0.984 | dB/(cm·MHz^y) |
    | Attenuation power | 1.2 | 1 | 1 | 1 | 1 | – |
    | Thermal conductivity | 0.51 | 0.48 | 0.55 | – | – | W/(m·K) |
    | Specific heat | 3630 | 3525 | 3718 | – | – | J/(kg·K) |
    | Perfusion | 559 | 196 | 809 | – | – | mL/min/kg |
    | Absorption fraction | 1 | 1 | 0.29 | – | – | – |

    **Skull**

    The variant values apply to `skull`, `skull_trabecular`, and `skull_cortical` alike. The Default column shows `skull`; the default values of the two sub-compartments are in [Medium](doc_medium.md).

    | Property | Default | Temp-H | Temp-L | Trans-H | Trans-L | Unit |
    |---|---|---|---|---|---|---|
    | Sound speed | 2800 | 2871 | 2190 | 2190 | 2871 | m/s |
    | Density | 1850 | 1080 | 1991 | 1080 | 1991 | kg/m³ |
    | Attenuation coeff. | 13.3 | 8.04 | 23.4 | 8.04 | 23.4 | dB/(cm·MHz^y) |
    | Attenuation power | 1 | 1 | 1 | 1 | 1 | – |
    | Thermal conductivity | 0.32 | 0.30 | 0.36 | – | – | W/(m·K) |
    | Specific heat | 1793 | 826 | 2524 | – | – | J/(kg·K) |
    | Perfusion | 20 | 0 | 30 | – | – | mL/min/kg |
    | Absorption | fraction 0.28 | coeff. 8.04 | coeff. 2.1 | – | – | –, dB/(cm·MHz^y) |

    **Rationale.** Temp-H combines low skull attenuation and density, which let more energy into the brain, with high brain attenuation and full brain absorption, which convert more of it to heat, and with the thermal bounds that raise temperature (low conductivity, specific heat, and perfusion). Temp-L takes the opposite bound of each. Trans-H combines low skull impedance and attenuation with low brain attenuation, so more energy reaches the target; Trans-L takes the opposite bound of each. Skull density and attenuation take the same bound in Temp-H and Trans-H; brain sound speed, density and attenuation and skull sound speed take opposite bounds.

    **Multi-compartment skull models.** With separate `skull_cortical` and `skull_trabecular` compartments, both get the single ITRUSST skull bound. To assign compartment-specific ranges, edit their values per compartment (see [Customising medium property ranges](#customising-medium-property-ranges)).

## Workflow Stages

The uncertainty mode consists of seven stages, of which stages 2–6 can run in parallel:

```
Stage 1: Preprocessing & source setup  (serial; must complete before stages 2–6)
    ↓
Stage 2: Default simulation             ┐
Stage 3: Temp-H simulation              │  acoustic + thermal
Stage 4: Temp-L simulation              ├─ parallel jobs
Stage 5: Trans-H simulation             │  acoustic only
Stage 6: Trans-L simulation             ┘
    ↓
Stage 7: Uncertainty report generation  (serial; runs once all five output files are present)
```

Variants share the stage-1 head and grid caches; each maps its own medium, cached under its own affix. Only the properties the variant's bounds set change; all other settings (transducer, grid, paths) are inherited from the base config. With `pct.enabled: 1`, stages 3–6 use a homogeneous skull layer (see [Pseudo-CT runs](#pseudo-ct-runs)).

??? note "What happens inside"
    `prestus_pipeline_start` detects `simulation.uncertainty = true` on HPC and redirects to `uncertainty_pipeline` before submitting a single job. On MATLAB, `prestus_pipeline` handles the redirect after `path_log_setup` has run. The per-variant parameter structs built internally have `simulation.uncertainty` cleared, so there is no recursion. `uncertainty_pipeline` builds one parameter struct per stage from the base parameters and dispatches them:

    | Stage | Modules enabled | `io.output_affix` | Medium override |
    |---|---|---|---|
    | 1 | source setup only | `''` | none |
    | 2 | acoustic + thermal (as `modules.run_heating_sims`) | `<base>` | none (default) |
    | 3 | acoustic + thermal | `<base>_desc-temphigh` | `config_medium_temp_high.yaml` |
    | 4 | acoustic + thermal | `<base>_desc-templow` | `config_medium_temp_low.yaml` |
    | 5 | acoustic only | `<base>_desc-transhigh` | `config_medium_trans_high.yaml` |
    | 6 | acoustic only | `<base>_desc-translow` | `config_medium_trans_low.yaml` |
    | 7 | uncertainty report | `''` | none |

    `<base>` is the base config's `io.output_affix` (usually empty).

    You can also call `uncertainty_pipeline` directly if you need to bypass `prestus_pipeline`:

    ```matlab
    uncertainty_pipeline(parameters, options);
    ```

### Customising Pipeline Options

Uncertainty-specific options can be passed as the second argument to `prestus_pipeline`:

```matlab
% Output affixes (defaults: base affix, then _desc-temphigh, _desc-templow, ...):
options.affixes.default    = '';
options.affixes.temp_high  = '_desc-temphigh';
options.affixes.temp_low   = '_desc-templow';
options.affixes.trans_high = '_desc-transhigh';
options.affixes.trans_low  = '_desc-translow';

% HPC wall times:
options.stage1_timelimit  = '00:30:00';
options.sim_timelimit     = '04:00:00';
options.report_timelimit  = '00:30:00';

prestus_pipeline(parameters, options);
```

`options.liberal_config`, `options.conservative_config`, and `options.affixes.liberal` / `.conservative` stop the run with `uncertainty_pipeline:legacyVariant`. Variant medium files are not read: `options.temp_high_config` and the other `*_config` options stop with `uncertainty_pipeline:variantFile`.

## Customising Medium Property Ranges

The variants use the ITRUSST bounds in `config/uncertainty/`. To change single values, set them in the study config under `simulation.uncertainty_overrides`, per variant, tissue, and property:

```yaml
simulation:
  uncertainty: true
  uncertainty_overrides:
    temp_high:
      skull:
        alpha_coeff: 6.0
        absorption_coeff: 6.0
```

Each edit is merged onto the variant's bounds; tissues and properties neither lists keep their `medium_properties` values. An edit that sets `absorption_coeff` for a tissue clears the bound's `absorption_fraction`, and the other way round. The Trans variants run no heating, so editing one of their thermal properties stops the run. In the Nexus, unlock the *Values per variant* table on the *Uncertainty* sub-tab and edit a cell; edited values are shown in bold red.

!!! warning "Edited bounds"
    A variant with edited values no longer brackets the ITRUSST bounds. The pipeline logs a warning and the uncertainty report lists the edited values; report them with your results.

## Partial-FOV Bound Variants

Five whole-head simulations are the main cost of an uncertainty run. With `simulation.uncertainty_fov: partial`, only the default variant runs on the whole grid. The four bound variants run in the [partial acoustic FOV](doc_partial-fov.md), and one more run repeats the default medium in that FOV:

```yaml
simulation:
  uncertainty: true
  uncertainty_fov: partial
```

| Run | Medium | Acoustic grid |
|---|---|---|
| Default | default | whole grid (`grid.acoustic_fov_*` cleared) |
| Default, partial FOV (`<base>_desc-defaultfov`) | default | partial FOV |
| Temp-H, Temp-L, Trans-H, Trans-L | bounds | partial FOV |

The partial FOV uses `grid.acoustic_fov_length_mm` and `grid.acoustic_fov_diameter_mm` where set, else `'auto'`. The thermal simulation keeps the full grid unless `grid.thermal_use_acoustic_fov: 1`.

A cropped grid changes the field, for example by removing reflections from skull outside the box. The two default runs measure that change: for each metric, the report takes the factor *whole grid / partial FOV* and applies it to the bound variants. Intensity, pressure, MI, temperature rise, and CEM43 are multiplied by their factor; absolute temperatures have their rise above the baseline (`thermal.temp_0`) multiplied. Focal position and size, volumes, and the temperature and CEM43 time courses are not rescaled. A metric whose partial-FOV default is zero or missing has no factor and stays as simulated.

The report's **Partial-FOV Rescaling** section lists each factor, and `sub-NNN_<medium>_desc-uncertainty_fovmismatch<base>.csv` holds the factors with the bound variants' values before and after rescaling. The variants' own CSVs, images, and per-variant reports keep the partial-FOV values.

!!! warning "Assumption"
    The factor comes from the default medium and is assumed to hold at the bounds. A bound that changes how much energy the cropped region returns to the beam, such as a lower skull attenuation, can break that assumption. CEM43 is exponential in temperature, so its factor is a first-order correction. Check a factor far from 1 against a whole-grid run of that variant (`uncertainty_fov: config`) before relying on the range.

Not available with a multi-I<sub>SPPA</sub> sweep, async coupling, or sequential follow-up runs; these stop with an error before any job is submitted. On HPC the repeat is stage 7 (`PRESTUS-u7-sim-default_fov_sub-XXX`) and the report stage 8.

## Output Files

Each simulation variant writes its outputs into the shared `path.sim` output folder (or a subject subfolder), identified by `io.output_affix`. The Trans variants write no thermal outputs (`maxT`, `thermal_max`, `heating_res`).

| File pattern | Description | Kept when `save_matrices` is `0` or empty |
|---|---|---|
| `sub-NNN_layered<affix>.csv` | Acoustic and thermal metrics table | yes |
| `sub-NNN_layered_report<affix>.html` | Per-variant simulation report | yes |
| `sub-NNN_layered_intensity_<dim><affix>.png` | Intensity overlay image per slice dimension (x, y, z) | yes |
| `sub-NNN_layered_maxT_<dim><affix>.png` | Max temperature overlay image per slice dimension (x, y, z) | yes |
| `sub-NNN_layered_thermal_max<affix>.png` | Temperature vs. time plot | yes |
| `cache/sub-NNN_layered_parameters_<yyMMdd_HHmm><affix>.mat` | Parameter struct for the variant | yes |
| `cache/sub-NNN_layered_heating_res<affix>.mat` | Heating simulation matrices (large) | yes |
| `cache/sub-NNN_layered_rotated_scaled.mat` | Intermediate head model | **removed** |
| `cache/sub-NNN_layered_cropped_smoothed.mat` | Intermediate skull grid | **removed** |

!!! note "Automatic cleanup"
    When `parameters.io.save_matrices` is `0` or empty (the default), the pipeline deletes the files marked **removed** above after the uncertainty report (stage 7) has been written (`cleanup_uncertainty_intermediates`, on both the MATLAB and HPC platforms). Cleanup only runs after stage 7 completes, so a failed run can be resumed beforehand. The head-model cache it removes is also the stage-1 sentinel (below), so a later re-run rebuilds the head model; set `io.save_matrices: 1` to keep it.

The uncertainty report is written as `sub-NNN_layered_desc-uncertainty_report<base>.html`. The name ends with the default variant's affix, so runs with different base affixes keep separate reports. With a nested multi-[I<sub>SPPA</sub>](doc_glossary.md#isppa) sweep or async run there is one report per sweep point, e.g. `sub-NNN_layered_desc-uncertainty_report<base>_desc-isppa030.html` (see [below](#combining-with-multi-isppa-or-async-modes)).

## Uncertainty Report Contents

The uncertainty report is a self-contained HTML file. Acoustic metrics (I<sub>SPPA</sub>, MI, pressure, focal metrics, half-maximum volume) are reported with the Trans range (Trans-L to Trans-H); thermal metrics (temperature rise, maximum temperature, [CEM43](doc_glossary.md#cem43)) with the Temp range (Temp-L to Temp-H). Every table lists the values of all five variants; the Trans variants show N/A for thermal metrics.

- **Safety dashboard**: colour-coded cards for all ITRUSST consensus safety metrics (MI, CEM43, temperature rise, absolute temperature). Each card shows the **default** value with the range of the relevant pair annotated. The bar spans a metric-specific minimum (37 °C for absolute temperatures, 0 elsewhere) to `max(max_sim, limit)`; ticks mark the variant values. The card colour follows the highest value of the default and the pair.
- **Acoustic results**: range table with the values of all variants and the Trans range for key acoustic metrics. Row colour follows the highest value of the default and the Trans pair.
- **Thermal results**: range table with the values of all variants and the Temp range for thermal metrics (maxT, temperature rise, CEM43 per tissue); row colour follows the highest value of the default and the Temp pair. Side-by-side temperature-vs-time images (`thermal_max`) of the default and Temp variants. Per-layer and global timeseries plots show temperature and CEM43 over time with the Temp-L to Temp-H band shaded around the default trajectory.
- **Pseudo-CT skull vs. bounds**: with `pct.enabled: 1`, histograms of the default run's mapped skull properties against the band of the bound variants' skull values (see [Pseudo-CT runs](#pseudo-ct-runs)).
- **Additional maps**: side-by-side maximum temperature maps across the default and Temp variants.
- **Raw data tables**: full CSV output tables for each variant.

## Combining with Multi-I<sub>SPPA</sub> or Async Modes

Uncertainty composes with the two other chained modes. When the base config is a multi-I<sub>SPPA</sub> sweep (a vector `transducer.target_isppa_wcm2`) or an async multi-transducer run (`simulation.transducer_coupling: async`), each of the five variants runs that pipeline with its own medium, and stage 7 writes one uncertainty report per sweep point. The Trans variants run no thermal jobs. Uncertainty is always the outer mode, on every platform.

| | Multi-I<sub>SPPA</sub> sweep | Async transducers |
|---|---|---|
| Stage 1 | shared preprocessing, as above | not run: each variant's async stage 0 builds its own grid and medium |
| Stages 2–6 | per variant: [`multi_isppa_pipeline`](doc_multiintensity.md) (acoustic, one thermal job per target, summary) | per variant: [`async_transducer_pipeline`](doc_multitransducer.md) (stage 0, acoustic per transducer, combine and joint per sweep point) |
| Variant outputs | `<variant affix>_desc-isppaNNN` | `<variant affix>_desc-asyncJoint` or `_desc-asyncJointSjj` |
| Stage 7 | one report per target | one report per sweep point |

On HPC, `uncertainty_pipeline` submits each variant's inner jobs itself. Job names carry the stage and variant (`PRESTUS-u3-temp_high-mi2-thermal_isppa030_sub-001`), the multi-I<sub>SPPA</sub> acoustic jobs wait for stage 1, and every report job waits for all the variants' jobs. Options for the inner pipeline (e.g. `options.thermal_timelimit`) pass through the `options` struct given to `prestus_pipeline_start`.

Not supported: sequential follow-up runs (`options.sequential_configs`) in nested mode, and a multi-I<sub>SPPA</sub> sweep with more than one coherent transducer (as without uncertainty). Both stop with an error before any job is submitted.

## Resuming a Partial Run

Re-run the same call. On HPC, a stage whose sentinel file exists is skipped and no scheduler dependency is created for it; only the missing stages are resubmitted. To rerun a stage, delete its sentinel file before resubmitting. With `io.overwrite_files: 'always'`, stages 2–7 are always resubmitted (stage 1 still resumes from its cache); with `'force'`, stage 1 is resubmitted too. See [Overwrite behaviour](doc_outputs_reuse.md#overwrite-behaviour).

On the MATLAB platform there is no stage-level skip. Resumability relies on setting `io.overwrite_files: 'never'` in the base config (the shipped default is `'always'`): each file write inside the pipeline is skipped if the file already exists, so a restarted run fast-forwards through completed steps file by file.

??? note "Sentinel files"
    | Stage | Sentinel file | Location |
    |---|---|---|
    | 1, Preprocessing | `sub-NNN_<medium>_cropped_smoothed.mat` | `<output_dir>/cache/` |
    | 2, Default sim | `sub-NNN_<medium><base>.csv` | `<output_dir>/` |
    | 3, Temp-H sim | `sub-NNN_<medium><base>_desc-temphigh.csv` | `<output_dir>/` |
    | 4, Temp-L sim | `sub-NNN_<medium><base>_desc-templow.csv` | `<output_dir>/` |
    | 5, Trans-H sim | `sub-NNN_<medium><base>_desc-transhigh.csv` | `<output_dir>/` |
    | 6, Trans-L sim | `sub-NNN_<medium><base>_desc-translow.csv` | `<output_dir>/` |
    | 7, Report | `sub-NNN_<medium>_desc-uncertainty_report<base>.html` | `<output_dir>/` |

    In nested mode (above) the variants' inner pipelines apply their own per-stage checks, and each per-point report is skipped when its file exists.

## HPC Usage

When `platform` is `'slurm'` or `'qsub'`, all seven jobs are submitted immediately using scheduler dependencies (`--dependency=afterok` for [SLURM](https://slurm.schedmd.com/)). Stages 2–6 depend on stage 1; stage 7 depends on all five simulation jobs. If stage 1 fails, all downstream jobs are cancelled (`--kill-on-invalid-dep=yes`). No MATLAB session is held open during execution.

Stages 1 and 7 (preprocessing and report) do not require a GPU and are submitted to the scheduler's default CPU queue. Stages 2–6 inherit the GPU resources from the base config. To override the queue for CPU stages use `options.stage1_partition = 'gpu'` and `options.report_partition = 'gpu'`.

??? note "Job names and dependencies"
    Job names are `PRESTUS-u1-preproc_sub-XXX`, `PRESTUS-u<k>-sim-<variant>_sub-XXX` for stages k = 2–6 (e.g. `PRESTUS-u3-sim-temp_high_sub-001`) and `PRESTUS-u7-report_sub-XXX`.

    This requires that `hpc_submit_job` supports the `parameters.hpc.depend_job_id` (single) and `parameters.hpc.depend_job_ids` (array) fields; see [Running on a Cluster (HPC)](doc_hpc.md).

## Pseudo-CT Runs

With `pct.enabled: 1`, only the default variant maps the skull from the [pseudo-CT](doc_pseudoCT.md). The four bound variants use the same head geometry (preprocessing and grid caches) but give the `skull` layer the homogeneous values of their bounds. The brain and thermal bounds apply as in other runs. The `skull_cortical` and `skull_trabecular` layers are never mapped from the pseudo-CT.

The range then compares the pseudo-CT skull with homogeneous skulls at the ITRUSST bounds. It brackets the default only where the mapped skull properties lie within those bounds, and it does not express the uncertainty of the pseudo-CT transform itself. Mapped values can fall largely outside the bounds, so check them before interpreting the range:

- The uncertainty report's **Pseudo-CT Skull vs. Bounds** section plots the default run's skull density, sound speed, and attenuation (or a prediction, when they are not cached) against the band of the bound variants' skull values and states the share of voxels outside it.
- Before simulating, `pct_predict_skull_properties` or **Predict skull properties** on the Nexus Pseudo-CT sub-tab predicts the modelled skull properties. The Nexus shows the default run's values instead when they exist.

```matlab
parameters = load_parameters('config_study.yaml');
parameters.subject_id = 1;
stats = pct_predict_skull_properties(parameters);  % writes cache/sub-001_<medium>_pct_prediction.json
stats.properties.sound_speed                        % mapped skull sound speed: mean, percentiles, histogram
stats.bounds                                        % skull values of each bound variant
```

The prediction applies the configured `pct.mapping_*` to the skull voxels of the head-preprocessing cache when it exists (`stats.source = 'grid'`, as modelled), else of `pseudoCT.nii.gz` and `tissues_mask.nii.gz` before cropping, resampling, and smoothing (`'native'`).

??? note "References"
    The variant bounds are the lower and upper values of the ITRUSST tissue property table (in preparation) and the property directions of its one-at-a-time sensitivity analysis. Sources cited in the table for the bounds include:

    - Aubry, J.-F., et al. ITRUSST consensus on biophysical safety for transcranial ultrasound stimulation. *Brain Stimulation* **18**, 1896–1905 (2025). https://doi.org/10.1016/j.brs.2025.10.007 (safety limits of the report)
    - Aubry, J.-F. et al. Benchmark problems for transcranial ultrasound simulation: Intercomparison of compressional wave models. *J Acoust Soc Am* **152**, 1003–1019 (2022). https://doi.org/10.1121/10.0013426 (skull reference values).
    - Lyons, M. E., & Parker, K. J. Absorption and attenuation in soft tissues. II. Experimental results. *IEEE Trans Ultrason Ferroelectr Freq Control* **35**, 511–521 (1988). https://doi.org/10.1109/58.4189 (lower brain attenuation and absorption).
    - White, P. J., Clement, G. T., & Hynynen, K. Longitudinal and shear mode ultrasound propagation in human skull bone. *Ultrasound Med Biol* **32**, 1085–1096 (2006). https://doi.org/10.1016/j.ultrasmedbio.2006.03.015 (lower skull attenuation, upper skull absorption).
