# Sequential Simulations

!!! note "At a glance"
    - **What:** A chain of simulations in which each heating run starts from the temperature and CEM43 maps of the previous one.
    - **When you need it:** When you stimulate from different coordinates in sequence and want the cumulative thermal state (for stimulation in parallel, see [Multi-Transducer Modelling](doc_multitransducer.md)).
    - **Before you start:** One config per sonication, loaded as structs.

## Basic usage

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

Field names in `sequential_configs` must follow the pattern `config_<N>` where `<N>` is a non-negative integer (e.g. `config_2`, `config_3`, `config_11`). The dispatcher picks the field with the lowest numeric value first; only relative order matters. Non-numeric suffixes (e.g. `config_abc`) are silently ignored.

## Cache reuse

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

## Sequential simulation parameters

| Parameter | Description |
|---|---|
| `io.adopted_heatmap` | Path to a temperature heatmap NIfTI (`heating_end.nii.gz`) used as the thermal starting point. Set automatically from the previous run's output; override to supply a custom starting temperature. |
| `io.adopted_cem43` | Path to a CEM43 NIfTI (`CEM43_end.nii.gz`) carried forward for cumulative thermal dose. Set automatically; can be overridden. |
| `io.adopted_cem43_iso` | Path to the ISO-variant CEM43 NIfTI (`CEM43_iso_end.nii.gz`). Set automatically alongside `adopted_cem43`. |
| `io.preproc_affix` | Affix for grid/medium cache lookup (defaults to the base run's `output_affix`). |
| `io.acoustic_cache_affix` | Affix for acoustic cache lookup. Not set by default; acoustics are re-run unless you set this explicitly. |
| `io.thermal_cache_affix` | Affix for the thermal cache file (defaults to the run's own `io.output_affix`). |
| `options.sequential_configs` | Struct of follow-up configs to dispatch in order (see example above). |
| `options.sequential_cleanup_intermediate` | If `true`, per-run NIfTI and image files are deleted after the final summary report is successfully generated. Default: `false`. |

## Sequential simulations with uncertainty mode

`options.sequential_configs` can be combined with `parameters.simulation.uncertainty = true`. Each follow-up run is executed once per uncertainty variant (default, Temp-H, Temp-L, Trans-H, Trans-L). The default and Temp variants inherit the thermal end-state from the matched variant of the prior run; the Trans variants re-run the acoustic simulation only:

```matlab
parameters.simulation.uncertainty = true;

options.sequential_configs.config_2 = config_2;
options.sequential_configs.config_3 = config_3;

prestus_pipeline(parameters, options);
```

After all sequential runs complete, a sequential report is written for all variants. The base run appears as **run 1**; each follow-up adds a further run entry. The report shows temperature and CEM43 timeseries with the Temp-L to Temp-H band shaded around the default trajectory.

On HPC, job names are prefixed with `r1-` when sequential configs are present (e.g. `PRESTUS-r1-u2-sim-default_sub-001`) to keep scheduler job names unambiguous across chained submissions. See [doc_uncertainty.md](doc_uncertainty.md) for the uncertainty pipeline.

## Output affix handling

Every PRESTUS output file embeds an **affix**, a short string appended to the subject/medium prefix, that distinguishes files produced by different simulation variants or pipeline modes. PRESTUS assigns affixes from the active modes; you rarely set one yourself.

- A plain sequential run gets `_seq<N>` (e.g. `_seq2` for the second run).
- An uncertainty variant gets `_desc-temphigh`, `_desc-templow`, `_desc-transhigh` or `_desc-translow`.
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

??? note "Details: affix construction and report file lookup"
    **Affix sanitisation.** `load_parameters` sanitises the raw affix string before use: characters that are not alphanumeric, hyphens, underscores, or slashes are replaced with underscores, and the prefix `_desc-` is prepended if not already present. So `"pilot"` becomes `"_desc-pilot"`.

    **Full affix table.**

    | Pipeline mode | `io.output_affix` value |
    |---|---|
    | Single run, no user affix | `''` (empty) |
    | Single run, user sets `io.output_affix = 'myreg'` | `'_desc-myreg'` |
    | Uncertainty, default variant | `base_affix` |
    | Uncertainty, Temp-H / Temp-L variant | `base_affix + '_desc-temphigh'` / `'_desc-templow'` |
    | Uncertainty, Trans-H / Trans-L variant | `base_affix + '_desc-transhigh'` / `'_desc-translow'` |
    | Sequential run N, no uncertainty | `base_affix + '_seq<N>'` |
    | Sequential run N, default variant | `base_affix + '_seq<N>'` |
    | Sequential run N, other variants | `base_affix + '_desc-<variant>_seq<N>'` (e.g. `'_desc-temphigh_seq<N>'`) |
    | Sequential run N, custom affix supplied | the supplied affix, unchanged |

    `base_affix` is whatever `io.output_affix` was on the parameters struct passed to `uncertainty_pipeline` or `prestus_pipeline_start`. For most studies it is empty.

    **How report functions locate files.** `generate_uncertainty_report` reconstructs CSV paths as `<io.dir_output>/sub-NNN_<medium><variant_affix>.csv` rather than trusting `io.filename_table`, which may carry a stale value from an earlier pipeline stage. `generate_sequential_report` uses the explicit affix registry built by `sequential_pipeline` (one struct per run with `output_affix` and `thermal_cache_affix`) for all lookups; `io.filename_table` is a fallback only.

    **Thermal cache affix.** `io.thermal_cache_affix` controls the filename of the heating `.mat` cache (`sub-NNN_<medium>_heating_res<affix>.mat`). It defaults to the run's own `output_affix`. Set it explicitly only if you need to point a follow-up run at a specific pre-existing cache file.

See also: [doc_uncertainty.md](doc_uncertainty.md), [doc_multitransducer.md](doc_multitransducer.md), [doc_hpc.md](doc_hpc.md)
