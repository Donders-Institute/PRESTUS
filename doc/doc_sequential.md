# Sequential Simulations

A sonication that follows another can start in tissue that is still warm from the first and adds to the thermal dose already accumulated. Run independently, each simulation starts from the baseline temperature and zero [CEM43](doc_glossary.md#cem43), so this carryover is missing. Sequential mode chains the runs so that each heating simulation starts from the end-of-run temperature and CEM43 maps of the previous one.

On this page:

- [Planning a Sequence in the Nexus](#planning-a-sequence-in-the-nexus): stances fired one after another, set up in one config.
- [Basic Usage](#basic-usage): passing follow-up configs to the pipeline.
- [Cache Reuse](#cache-reuse): which stages a follow-up run reuses and which it recomputes.
- [Sequential Simulation Parameters](#sequential-simulation-parameters): the settings that carry the thermal state forward.
- [Combining with Uncertainty Mode](#combining-with-uncertainty-mode): chaining each uncertainty variant.
- [Output Affix Handling](#output-affix-handling): how files of each run and variant are named.

!!! note "At a glance"
    - **When you need it:** When you stimulate from different coordinates in sequence and want the cumulative thermal state (for stimulation in parallel, see [Multi-Transducer Modelling](doc_multitransducer.md)). A sequence planned in advance is easiest to set up in the [Nexus](#planning-a-sequence-in-the-nexus).
    - **Before you start:** One config per sonication, loaded as structs.

## Planning a Sequence in the Nexus

The *Sonication sequence* card (Sonication tab, under the timing diagram) sets up several stances that sonicate one after another. **Add a later sonication** adds a stance as a new step, which you place in the Viewer. Each step has its own duration and pulse timing (**Pulse timing…**), and each pause is modelled as off time. Transducers that should fire at the same time go into one step (**Fire with previous step**). The strip above the table shows the whole protocol to scale.

The sequence is stored in the config with the [async `custom` schedule](doc_multitransducer_async.md#firing-schedule). Each step's start goes in `transducer(i).timing_offset_s`, and each step's timing goes in `transducer(i).timing` wherever it differs from the `timing` block, which holds the step with transducer 1. *Simulate as* chooses how the sequence runs:

| Choice | How it runs | Possible when |
|---|---|---|
| One joint run | The config as it is: one acoustic simulation per transducer and one thermal simulation for the whole sequence. The transducers of a step fire together as a [firing group](doc_multitransducer_async.md#firing-schedule): their fields are summed coherently, and steps whose pulses partly overlap interfere while both are on. | Always. |
| Chain of runs | One config per step (`<config>_step-<k>.yaml`, written when you run), chained as sequential follow-ups. Each pause becomes the previous step's `timing.post_ptri_dur`, and each step gets its own grid. A step with several transducers is one coherent run. | Steps do not overlap; each pause is a whole number of `timing.post_pt_timestep`; no uncertainty mode and no list of target I<sub>SPPA</sub>s. |
| Automatic (default) | One joint run. | |

With one joint run and `grid.async_layout: per_transducer`, heat outside transducer 1's grid is lost; use the shared grid or a chain for stances far apart. The card lists what stops a choice, and so does the Issues panel for the automatic choice. [Example Protocols](doc_multitransducer_async.md#example-protocols) shows sequences, transducers firing together, interleaved pulses, and a mixed protocol with uncertainty as the card draws them.

## Basic Usage

Each config is a struct. Given `config_1`, `config_2`, and `config_3` to run in sequence:

Without heat carryover, run each independently:

```matlab
prestus_pipeline_start(config_1)
prestus_pipeline_start(config_2)
prestus_pipeline_start(config_3)
```

With heat carryover, pass the follow-up configs via `options`. The thermal end state is handed over as T1w-space NIfTIs (`*_heating_end`, `*_CEM43_end`, `*_CEM43_iso_end`). Each follow-up resamples them into its own grid with its own transform, so the start state is correct for a follow-up with a different beam axis. The NIfTI stage is enabled for every run in a chain, and dispatch stops with `sequential_pipeline:missingThermalEnd` if a map is missing. It never starts from the water baseline silently.

```matlab
sequential_configs.config_2 = config_2;
sequential_configs.config_3 = config_3;
options.sequential_configs = sequential_configs;
prestus_pipeline_start(config_1, options);
```

Field names in `sequential_configs` must follow the pattern `config_<N>` where `<N>` is a non-negative integer (e.g. `config_2`, `config_3`, `config_11`). The dispatcher picks the field with the lowest numeric value first; only relative order matters. A field that does not match the pattern (e.g. `config_abc`) stops the run with the error `sequential_pipeline:badFieldName`.

In the Nexus, *Results ▸ Follow-up…* with *Sequential* creates a follow-up config and adds it to the chain under *Run ▸ Where to run ▸ Pipeline options*, which the run script passes as `options.sequential_configs` ([Run tab](doc_nexus_tabs.md#run-tab)).

## Cache Reuse

When a follow-up simulation is dispatched, PRESTUS decides which results are reused and which are recomputed:

| Stage | Default behaviour | How to override |
|---|---|---|
| Grid & medium setup | **Reused** from base run | Set `io.preproc_affix` in the sequential config |
| Acoustic simulation | **Re-run** with the follow-up config | Set `io.acoustic_cache_affix` to an existing affix to reuse prior acoustics |
| Thermal simulation | **Re-run** (cache file gets a unique name automatically) | Set `io.thermal_cache_affix` explicitly to control the cache filename |

The head geometry does not change between sequential targets, so preproc can always be shared. Acoustics default to a fresh run because the follow-up may target a different location or use different transducer settings.

To reuse the base run's acoustic results in the follow-up (e.g. same target, different timing protocol):

```matlab
options.sequential_configs.config_2.io.acoustic_cache_affix = config_1.io.output_affix;
```

## Sequential Simulation Parameters

| Parameter | Description |
|---|---|
| `io.adopted_heatmap` | Path to a temperature heatmap [NIfTI](https://nifti.nimh.nih.gov/) (`heating_end.nii.gz`) used as the thermal starting point. Set automatically from the previous run's output; override to supply a custom starting temperature. |
| `io.adopted_cem43` | Path to a CEM43 NIfTI (`CEM43_end.nii.gz`) carried forward for cumulative thermal dose. Set automatically; can be overridden. |
| `io.adopted_cem43_iso` | Path to the ISO-variant CEM43 NIfTI (`CEM43_iso_end.nii.gz`). Set automatically alongside `adopted_cem43`. |
| `io.preproc_affix` | Affix for grid/medium cache lookup (defaults to the base run's `output_affix`). |
| `io.acoustic_cache_affix` | Affix for acoustic cache lookup. Not set by default; acoustics are re-run unless you set this explicitly. |
| `io.thermal_cache_affix` | Affix for the thermal cache file (defaults to the run's own `io.output_affix`). |
| `options.sequential_configs` | Struct of follow-up configs to dispatch in order (see example above). |
| `options.sequential_cleanup_intermediate` | If `true`, per-run NIfTI and image files are deleted after the final summary report is successfully generated. Default: `false`. |

## Combining with Uncertainty Mode

`options.sequential_configs` can be combined with `parameters.simulation.uncertainty = true` when every run keeps the first run's transducer positions: the variants of a follow-up reuse the first run's grid. For stances that change position, use one [async run](doc_multitransducer_async.md#firing-schedule) instead, which uncertainty mode wraps. Each follow-up run is executed once per uncertainty variant (default, Temp-H, Temp-L, Trans-H, Trans-L). The default and Temp variants inherit the thermal end-state from the matched variant of the prior run; the Trans variants re-run the acoustic simulation only:

```matlab
parameters.simulation.uncertainty = true;

options.sequential_configs.config_2 = config_2;
options.sequential_configs.config_3 = config_3;

prestus_pipeline(parameters, options);
```

After all sequential runs complete, a sequential report is written for all variants. The base run appears as **run 1**; each follow-up adds a further run entry. The report shows temperature and CEM43 timeseries with the Temp-L to Temp-H band shaded around the default trajectory.

On HPC, job names are prefixed with `r1-` when sequential configs are present (e.g. `PRESTUS-r1-u2-sim-default_sub-001`) to keep scheduler job names unambiguous across chained submissions. See [Uncertainty Quantification](doc_uncertainty.md) for the uncertainty pipeline.

## Output Affix Handling

Every PRESTUS output file embeds an **affix**, a short string appended to the subject/medium prefix, that distinguishes files produced by different simulation variants or pipeline modes. PRESTUS assigns affixes from the active modes; you rarely set one yourself.

- A plain sequential follow-up appends `_seq<N>` to the affix of the run before it, with N the number of its `config_<N>` field (e.g. `_seq2`, then `_seq2_seq3`).
- An uncertainty variant gets `_desc-temphigh`, `_desc-templow`, `_desc-transhigh`, or `_desc-translow`.
- Both combined: `_desc-temphigh_seq2`.
- A custom `io.output_affix` on any sequential config is used as-is; no `_seq<N>` suffix is added.

Example filenames for subject 1, `layered` medium, two sequential runs with uncertainty:

```
sub-001_layered.csv                              ← base run, default
sub-001_layered_desc-temphigh.csv                ← base run, Temp-H
sub-001_layered_desc-templow.csv                 ← base run, Temp-L
sub-001_layered_desc-transhigh.csv               ← base run, Trans-H
sub-001_layered_desc-translow.csv                ← base run, Trans-L
sub-001_layered_seq2.csv                         ← sequential run 2, default
sub-001_layered_desc-temphigh_seq2.csv           ← sequential run 2, Temp-H
…                                                ← and so on for the other variants
```

If you see a warning that a variant output was not found but the file exists on disk, confirm that `io.dir_output` was set before calling the pipeline (i.e. `path_log_setup` has run), and that no custom `options.affixes` struct was passed with different values than those used during simulation.

??? note "Affix construction and report file lookup"
    **Affix sanitisation.** `load_parameters` sanitises the raw affix string before use: characters other than letters, digits, hyphens, and underscores are replaced with underscores; leading underscores and a leading `desc-` are stripped, and `_desc-` is prepended. So `"pilot"` becomes `"_desc-pilot"`.

    **Full affix table.**

    | Pipeline mode | `io.output_affix` value |
    |---|---|
    | Single run, no user affix | `''` (empty) |
    | Single run, user sets `io.output_affix = 'myreg'` | `'_desc-myreg'` |
    | Uncertainty, default variant | `base_affix` |
    | Uncertainty, Temp-H / Temp-L variant | `base_affix + '_desc-temphigh'` / `'_desc-templow'` |
    | Uncertainty, Trans-H / Trans-L variant | `base_affix + '_desc-transhigh'` / `'_desc-translow'` |
    | Sequential run N, no uncertainty | affix of the run before it `+ '_seq<N>'` (N from its `config_<N>` field) |
    | Sequential run N, default variant | `base_affix + '_seq<N>'` |
    | Sequential run N, other variants | `base_affix + '_desc-<variant>_seq<N>'` (e.g. `'_desc-temphigh_seq<N>'`) |
    | Sequential run N, custom affix supplied | the supplied affix, unchanged |

    `base_affix` is whatever `io.output_affix` was on the parameters struct passed to `uncertainty_pipeline` or `prestus_pipeline_start`. For most studies it is empty.

    **How report functions locate files.** `generate_uncertainty_report` reconstructs CSV paths as `<io.dir_output>/sub-NNN_<medium><variant_affix>.csv` rather than trusting `io.filename_table`, which may carry a stale value from an earlier pipeline stage. `generate_sequential_report` uses the explicit affix registry built by `sequential_pipeline` (one struct per run with `output_affix` and `thermal_cache_affix`) for all lookups; `io.filename_table` is a fallback only.

    **Thermal cache affix.** `io.thermal_cache_affix` controls the filename of the heating `.mat` cache (`sub-NNN_<medium>_heating_res<affix>.mat`). It defaults to the run's own `output_affix`. Set it explicitly only if you need to point a follow-up run at a specific pre-existing cache file.

See also: [Uncertainty Quantification](doc_uncertainty.md), [Multi-Transducer Modelling](doc_multitransducer.md), [Running on a Cluster (HPC)](doc_hpc.md)
