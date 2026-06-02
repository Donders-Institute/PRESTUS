# Sequential Simulations

Applicable when you stimulate from different coordinates in sequence (for stimulation in parallel, see [Multi-Transducer Modeling](doc_multitransducer.md)). Instead of starting each heating simulation with the default starting temperatures, you can start your nth stimulation with the temperature and CEM43 maps from the previous simulation. To do this, you need to feed the pipeline your other configurations.

---

## Basic usage

Each config is a struct, so multiple configs can be placed in one struct without conversion. Given `config_1`, `config_2`, and `config_3` to run in sequence:

Without heat carryover — run each independently:

```matlab
prestus_pipeline_start(config_1)
prestus_pipeline_start(config_2)
prestus_pipeline_start(config_3)
```

With heat carryover — pass subsequent configs via `options`:

```matlab
sequential_configs.config_2 = config_2;
sequential_configs.config_3 = config_3;
options.sequential_configs = sequential_configs;
prestus_pipeline_start(config_1, options);
```

Field names in `sequential_configs` must follow the pattern `config_<N>` where `<N>` is a non-negative integer (e.g. `config_2`, `config_3`, `config_11`). The dispatcher picks the field with the lowest numeric value first — only relative order matters. Non-numeric suffixes (e.g. `config_abc`) are silently ignored.

---

## Cache reuse

When a follow-up simulation is dispatched, PRESTUS automatically manages which results are reused and which are recomputed:

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

---

## Sequential simulation parameters

| Parameter | Description |
|---|---|
| `io.adopted_heatmap` | Path to a temperature heatmap NIfTI (`heating_end.nii.gz`) used as the thermal starting point. Set automatically from the previous run's output; override to supply a custom starting temperature. |
| `io.adopted_cem43` | Path to a CEM43 NIfTI (`CEM43_end.nii.gz`) carried forward for cumulative thermal dose. Set automatically; can be overridden. |
| `io.adopted_cem43_iso` | Path to the ISO-variant CEM43 NIfTI (`CEM43_iso_end.nii.gz`). Set automatically alongside `adopted_cem43`. |
| `io.preproc_affix` | Affix for grid/medium cache lookup (defaults to the base run's `output_affix`). |
| `io.acoustic_cache_affix` | Affix for acoustic cache lookup. Not set by default — acoustics are re-run unless you set this explicitly. |
| `io.thermal_cache_affix` | Affix for the thermal cache file (defaults to the run's own `io.output_affix`). |
| `options.sequential_configs` | Struct of follow-up configs to dispatch in order (see example above). |
| `options.sequential_cleanup_intermediate` | If `true`, per-run NIfTI and image files are deleted after the final summary report is successfully generated. Default: `false`. |

---

## Sequential simulations with uncertainty mode

`options.sequential_configs` can be combined with `parameters.simulation.uncertainty = true`. Each follow-up run is executed three times — once per uncertainty variant (default / liberal / conservative) — and each variant inherits the thermal end-state from the matched variant of the prior run:

```matlab
parameters.simulation.uncertainty = true;

options.sequential_configs.config_2 = config_2;
options.sequential_configs.config_3 = config_3;

prestus_pipeline(parameters, options);
```

After all sequential runs complete, `generate_sequential_report` is called with the full list of default/liberal/conservative parameter structs. The base run appears as **run 1**; each follow-up adds a further run entry. The report shows temperature and CEM43 timeseries with liberal/conservative shaded uncertainty bands overlaid on the default trajectory.

On HPC, job names are prefixed with `r1-` when sequential configs are present (e.g. `PRESTUS-r1-u2-sim-default_sub-001`) to keep scheduler job names unambiguous across chained submissions. See [doc_uncertainty.md](doc_uncertainty.md) for the full description of the uncertainty pipeline.

---

## Output affix handling

Every PRESTUS output file embeds an **affix** — a short string appended to the subject/medium prefix — that distinguishes files produced by different simulation variants or pipeline modes. You generally do not need to set this manually; PRESTUS assigns affixes automatically based on the active modes.

The affix for each run encodes which mode produced it:

- A plain sequential run gets `_seq<N>` (e.g. `_seq2` for the second run).
- An uncertainty variant gets `_desc-liberal` or `_desc-conservative`.
- Both combined: `_desc-liberal_seq2`.
- A custom `io.output_affix` on any sequential config is used as-is — no `_seq<N>` suffix is added.

Example filenames for subject 1, `layered` medium, two sequential runs with uncertainty:

```
sub-001_layered.csv                              ← base run, default
sub-001_layered_desc-liberal.csv                 ← base run, liberal
sub-001_layered_desc-conservative.csv            ← base run, conservative
sub-001_layered_seq2.csv                         ← sequential run 2, default
sub-001_layered_desc-liberal_seq2.csv            ← sequential run 2, liberal
sub-001_layered_desc-conservative_seq2.csv       ← sequential run 2, conservative
```

If you see a warning like *"Conservative variant output not found"* but the file exists on disk, confirm that `io.dir_output` was set before calling the pipeline (i.e. `path_log_setup` has run), and that no custom `options.affixes` struct was passed with different values than those used during simulation.

<details markdown>
<summary>Developer details — affix construction and report file lookup</summary>

**Affix sanitisation.** `load_parameters` sanitises the raw affix string before use: characters that are not alphanumeric, hyphens, underscores, or slashes are replaced with underscores, and the prefix `_desc-` is prepended if not already present. So `"liberal"` becomes `"_desc-liberal"`.

**Full affix table.**

| Pipeline mode | `io.output_affix` value |
|---|---|
| Single run, no user affix | `''` (empty) |
| Single run, user sets `io.output_affix = 'myreg'` | `'_desc-myreg'` |
| Uncertainty — default variant | `base_affix` |
| Uncertainty — liberal variant | `base_affix + '_desc-liberal'` |
| Uncertainty — conservative variant | `base_affix + '_desc-conservative'` |
| Sequential run N, no uncertainty | `base_affix + '_seq<N>'` |
| Sequential run N — default variant | `base_affix + '_seq<N>'` |
| Sequential run N — liberal variant | `base_affix + '_desc-liberal_seq<N>'` |
| Sequential run N — conservative variant | `base_affix + '_desc-conservative_seq<N>'` |
| Sequential run N, custom affix supplied | the supplied affix, unchanged |

`base_affix` is whatever `io.output_affix` was on the parameters struct passed to `uncertainty_pipeline` or `prestus_pipeline_start`. For most studies it is empty.

**How report functions locate files.** `generate_uncertainty_report` reconstructs CSV paths as `<io.dir_output>/sub-NNN_<medium><variant_affix>.csv` rather than trusting `io.filename_table`, which may carry a stale value from an earlier pipeline stage. `generate_sequential_report` uses the explicit affix registry built by `sequential_pipeline` (one struct per run with `output_affix` and `thermal_cache_affix`) for all lookups; `io.filename_table` is a fallback only.

**Thermal cache affix.** The `io.thermal_cache_affix` field controls the filename of the heating `.mat` cache (`sub-NNN_<medium>_heating_res<affix>.mat`). It defaults to the run's own `output_affix`. Set it explicitly only if you need to point a follow-up run at a specific pre-existing cache file.

</details>

---

See also: [doc_uncertainty.md](doc_uncertainty.md), [doc_multitransducer.md](doc_multitransducer.md), [doc_hpc.md](doc_hpc.md)
