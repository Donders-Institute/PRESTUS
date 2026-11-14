# Reusing Earlier Outputs and Caches

A run can start from results an earlier run wrote instead of recomputing them: a segmentation, pseudo-CT, placement or pressure map, or the caches in its own `cache/` folder. Reuse saves time, but a reused file built from other settings gives results that do not match the configuration, so PRESTUS records and checks the settings each file depends on.

On this page:

- [Starting from earlier stage outputs](#starting-from-earlier-stage-outputs): starting a run from a placement, pseudo-CT or pressure map.
- [Overwrite behaviour](#overwrite-behaviour): when existing files are reused or recomputed.
- [Cache reuse checks](#cache-reuse-checks): how a reused cache is checked against the current settings.

!!! note "At a glance"
    - **When you need it:** You re-run or resume a subject, run a chained pipeline, or start a run from an earlier placement or pressure map.
    - **Before you start:** An earlier run's output folder; see [Simulation Outputs](doc_outputs.md) for the folder layout.

---

## Starting from Earlier Stage Outputs

PRESTUS keeps two kinds of intermediate results apart:

- **Within-run caches** (`cache/` of a simulation folder) are reloaded only by the run that wrote them, or a re-run of it. `io.overwrite_files` decides reuse ([Overwrite behaviour](#overwrite-behaviour)).
- **Stage outputs** are complete results that a later, separate run starts from. Each has a [BIDS](https://bids.neuroimaging.io/)-style JSON sidecar whose `signature` lists the settings the output depends on. The later run compares it with its own settings.

??? note "Stage outputs and compared settings"
    | Stage output | Sidecar | Picked up by | Compared settings |
    |---|---|---|---|
    | SimNIBS segmentation, `path.seg/m2m_sub-XXX/` | none (standard SimNIBS output) | every layered run, when `final_tissues.nii.gz` exists (`io.overwrite_simnibs`) | SimNIBS version (`charm_log.html`) and charm settings (`settings.ini`): a warning, see [Segmentation](doc_preproc.md#simnibs-version-and-charm-settings) |
    | Pseudo-CT, `<pct folder>/pseudoCT.nii.gz` | `pseudoCT.json` | every run with `pct.enabled: 1`, when it exists | `pct.skull_mapping` (with `custom_slope` / `custom_intercept`) and the segmentation (checksum of `final_tissues.nii.gz`) |
    | Placement run, `<path.placement>/sub-XXX/<mode>/<target_name>/` | `placement_run.json` | heuristic, PlanTUS and TUS_entry placement, per `placement.reuse` | `placement.<mode>` (without tool install paths), the transducer without its positions, and the segmentation (checksum of `final_tissues.nii.gz`) |
    | Pressure map, `nii/sub-XXX_<medium>_T1w<affix>_pressure.nii.gz` | `…_pressure.json` | a run with `path.acoustic_input` set to the NIfTI | `simulation.medium`, `grid` (without `thermal_*`), `layers`, `headmodel`, `pct`, `medium_properties`, the transducer (geometry and grid positions), the subject and, for `layered`, the segmentation |

`io.cache_signature_check` handles a mismatch (`error` by default, `warn`, `off`). An output without a sidecar (e.g. a pressure map from another tool) is used with the warning `prestus:cache:unverified`.

**Placement (`placement.reuse`).**

| Value | Effect |
|---|---|
| `auto` (default) | Reuse a recorded run whose settings match; otherwise run the placement again and replace the record. Under `io.overwrite_files: never`, a mismatch is handled like `always` |
| `always` | Require a recorded run |
| `never` | Always run the tool |

`placement.reuse` is independent of `io.overwrite_files`: a PlanTUS or TUS_entry result is reused while the simulation stages are recomputed. `placement_run.mat` records are also read.

**Pressure map (`path.acoustic_input`).** Point `path.acoustic_input` at a T1w pressure NIfTI from an earlier run. The new run skips the acoustic simulation and the free-water baseline, puts the map onto its own grid (including an acoustic FOV), rescales it to its `transducer.target_isppa_wcm2` using the sidecar's `isppa.level_wcm2`, and then runs acoustic analysis and the thermal stage as after a simulation. It needs the T1 → grid transform (`simulation.medium` layered or phantom, a 3D grid) and the k-Wave grid of the source stage (`modules.run_source_setup: 1`).

Limits:

- The map is stored at T1 resolution, so a map from a finer simulation grid has lost some peak sharpness.
- The map has no peak negative pressure, so the mechanical index uses the peak pressure (with a warning).
- A map without a sidecar is used as it is and cannot be rescaled.

??? note "Sidecar format"
    Every sidecar starts with the same header:

    ```json
    {
      "sidecar_version": 1,
      "kind": "acoustic_pressure",
      "created": "2026-10-07T16:00:00",
      "generated_by": {"name": "PRESTUS", "version": "v0.6.1-42-gabc1234"},
      "subject_id": 1,
      "...": "kind-specific fields",
      "signature": {"...": "settings the output depends on"}
    }
    ```

---

## Overwrite Behaviour

`io.overwrite_files` controls whether existing files are reused or recomputed.

| Value | Single run | Chained pipelines: intermediate stages | Chained pipelines: output stages |
|---|---|---|---|
| `never` | Reuse every existing file | Resume from cache | Skip if the output exists |
| `always` (default) | Recompute everything | Resume from cache | Always re-run, overwrite outputs |
| `force` | Same as `always` | Always re-run | Always re-run, overwrite outputs |
| `ask` | Prompt per file (interactive only) | Resume from cache | Skip if the output exists |

Under `never` or `ask`, an acoustic or thermal stage that recomputes its results (no cache to reuse, e.g. with `io.save_matrices: 0`) still replaces its existing NIfTI outputs, so the NIfTI and CSV outputs of a folder always come from the same run. `io.overwrite_simnibs` governs SimNIBS segmentations separately.

Use `force` when an input of the intermediate caches has changed: transducer placement, transducer or source parameters, grid, segmentation, or medium properties. Grid, medium and acoustic caches detect such changes when they are reused ([Cache reuse checks](#cache-reuse-checks)); `force` rebuilds them instead of stopping the run. Alternatively, delete the relevant `cache/` files before resubmitting.

??? note "Chained pipelines"
    Chained pipelines split one run into separate jobs. Each job reuses what the previous stage wrote:

    | Pipeline | Intermediate stages | Output stages |
    |---|---|---|
    | Async multi-transducer ([doc](doc_multitransducer.md)) | Stage-0 preproc, per-transducer acoustic, combine | Joint thermal + report, water validation |
    | Multi-ISPPA ([doc](doc_multiintensity.md)) | Stage-1 acoustic | Per-ISPPA thermal + report, summary report |
    | Uncertainty ([doc](doc_uncertainty.md)) | Stage-1 preprocessing | Variant simulations + reports, uncertainty report |

    With the default `always`, resubmitting a finished chained run reuses the acoustic simulations and regenerates the CSV, NIfTI and HTML outputs. Within a chained run, sub-jobs that consume a shared cache (e.g. per-transducer acoustic jobs loading the Stage-0 grid) always load it, even under `force`. Only the stage that owns the cache rebuilds it.

---

## Cache Reuse Checks

Grid, medium and acoustic caches are named by subject, medium and affix only. Each cache therefore stores the settings it depends on (`cache_signature`). A reused cache built from different settings stops the run by default (`io.cache_signature_check`).

??? note "Mismatch handling"
    When a cache is reused (under `io.overwrite_files: never`, after a declined overwrite prompt, or in the intermediate stages of a chained pipeline under `always`) and the settings differ, `io.cache_signature_check` decides:

    | Value | Effect |
    |---|---|
    | `error` (default) | Stop with `prestus:cache:signatureMismatch`, naming the differing keys. |
    | `warn` | Warn and reuse the cache anyway. |
    | `off` | No check. |

    A cache without a signature is reused with the warning `prestus:cache:unverified`. Stage-output sidecars are checked the same way ([Starting from earlier stage outputs](#starting-from-earlier-stage-outputs)).

??? note "What each cache compares"
    | Cache | Compared settings |
    |---|---|
    | Grid (`*_grid_cache*.mat`) | `simulation.medium`, `grid`, `layers`, `headmodel`, `pct`, and the transducer (geometry and positions after placement) |
    | Medium (`*_medium_cache*.mat`) | the grid settings above, plus `medium_properties` |
    | Acoustic (`*_results*.mat`) | the medium settings above except `absorption_coeff` (thermal only), plus `subject_id` and `grid.acoustic_fov_*`; not `grid.thermal_*` |

    For `simulation.medium: phantom`, each also compares the phantom file (size and modification date of `final_tissues.nii.gz`). For `layered`, each compares the segmentation by content (checksum of `final_tissues.nii.gz`), so a re-run of charm with another SimNIBS version or other settings is detected, while a copied or moved segmentation still matches. A cache written before the checksum was recorded counts the segmentation as changed when `final_tissues.nii.gz` is newer than the cache file.

    Not compared, because they do not change these caches: paths, I/O, HPC and module flags; the drive level (`target_isppa_wcm2`, `tx_calibration`), which only scales the acoustic field; and, for grid and medium caches, the partial acoustic FOV (`grid.acoustic_fov_*`), which is cropped after the grid and medium stages. A run that reuses a cache keeps its own values of these settings.

    Sub-jobs that load a cache owned by another stage of the same chained run are not checked.
