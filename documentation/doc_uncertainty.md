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
| **Liberal** | Parameter combination that tends to produce *higher* acoustic intensity and heating at the target | Upper bound on in-situ exposure |
| **Conservative** | Parameter combination that tends to produce *lower* acoustic intensity and heating at the target | Lower bound on in-situ exposure |

> **Safety assessment.** For the purpose of regulatory or ethical review, the **liberal** variant represents the highest plausible exposure: if safety limits are met under liberal tissue assumptions, they are expected to be met under any plausible tissue configuration. The default simulation provides the best-estimate, and the conservative variant documents the lower end of the plausible range.

The liberal/conservative distinction follows the logic that different parameter choices affect different stages of wave propagation in opposite directions. For example, higher skull attenuation reduces intracranial intensity (beneficial) but increases skull heating (adverse). The liberal and conservative parameter sets are therefore not simply "high" and "low" across all properties — they are tailored to yield the highest and lowest plausible values for the metrics of interest (intracranial intensity and heating).

---

## Medium property ranges

The following table summarises the default, liberal, and conservative medium properties currently implemented in PRESTUS. Liberal values are chosen to maximise acoustic exposure and heating at the target; conservative values minimise them.

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

> **Rationale for skull properties.** The conservative variant uses high attenuation and high acoustic impedance (high density × sound speed), which reflects a denser, more calcified skull: more acoustic energy is lost before reaching the brain, so intracranial intensity is lower. Conversely, the liberal variant uses a skull with low impedance and low attenuation, which allows more energy through to the brain. The absorption fraction modulates how much of the attenuated energy is converted to heat: a lower fraction means less skull heating per unit attenuation. In practice, the liberal variant thus produces higher intracranial intensity *and* allows more conservative skull heating estimates relative to the default.

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
addpath(genpath(fullfile(prestus_path, 'functions')));
cd(prestus_path);

parameters = load_parameters('config_study.yaml');
parameters.subject_id             = 1;
parameters.simulation.medium      = 'layered';
parameters.simulation.uncertainty = true;   % activates uncertainty mode
parameters.platform               = 'auto'; % 'matlab', 'slurm', or 'qsub'

prestus_pipeline(parameters);
```

`prestus_pipeline` detects `simulation.uncertainty = true` at startup and redirects to `uncertainty_pipeline`, which handles module flags, output affixes, medium config merging, and (on HPC) job chaining. The per-variant parameter structs built internally have `simulation.uncertainty` cleared, so there is no recursion.

On a local MATLAB session (`platform = 'matlab'`) the five stages run sequentially. On HPC (`platform = 'slurm'` or `'qsub'`) stages 2–4 are submitted in parallel with a dependency on stage 1, and stage 5 is submitted with a dependency on all three simulation jobs.

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

The built-in medium override configs live in `configs/uncertainty/`:

- `config_medium_liberal.yaml` — low-impedance, low-attenuation skull; higher brain intensity
- `config_medium_conservative.yaml` — high-impedance, high-attenuation skull; lower brain intensity

---

## Output files

Each simulation variant writes its outputs into the shared `path.sim` output folder (or a subject subfolder), identified by `io.output_affix`:

| File pattern | Description |
|---|---|
| `sub-NNN_layered_output_table<affix>.csv` | Acoustic and thermal metrics table |
| `sub-NNN_layered_heating_res<affix>.mat` | Heating simulation matrices and timeseries |
| `sub-NNN_layered_intensity<affix>.png` | Intensity overlay image |
| `sub-NNN_layered_maxT<affix>.png` | Max temperature overlay image |
| `sub-NNN_layered_report<affix>.html` | Per-variant simulation report |

The uncertainty report is written as:

| File | Description |
|---|---|
| `sub-NNN_layered_uncertainty_report.html` | Unified uncertainty HTML report |

---

## Uncertainty report contents

The uncertainty report (`generate_uncertainty_report.m`) is a self-contained HTML file that includes:

- **Safety dashboard** — Color-coded cards for all ITRUSST consensus safety metrics (MI, CEM43, temperature rise, absolute temperature), showing the worst-case (liberal) value, the conservative-to-liberal range, and a progress bar relative to the limit.
- **Acoustic results** — Range table showing conservative / default / liberal values for key acoustic metrics (ISPPA, focal distance, MI, half-maximum volume), with color coding of the worst-case value against the limit. Side-by-side intensity map comparison across all three variants.
- **Thermal results** — Range table for thermal metrics (maxT, temperature rise, CEM43 per tissue). Inline SVG timeseries plot with the default trajectory and a shaded uncertainty band spanning the conservative-to-liberal range. Side-by-side maximum temperature maps.
- **Image gallery** — Full side-by-side comparison of all output images across variants.
- **Raw data tables** — Full CSV output tables for each variant.

---

## HPC usage

When `platform` is `'slurm'` or `'qsub'`, `uncertainty_pipeline` submits all five jobs in the correct order using scheduler dependencies (`--dependency=afterok` for SLURM). Stage 1 is submitted first with `wait_for_job = true`, blocking until it completes before stages 2–4 are submitted. Stage 5 is submitted with a dependency on all three simulation job IDs and runs automatically once they have all succeeded.

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
