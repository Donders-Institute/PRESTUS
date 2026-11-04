# Simulation outputs

!!! note "At a glance"
    - **What:** Every file, table column and folder a run writes, and how a later run can start from them.
    - **When you need it:** To find a specific file or column, or to reuse an earlier run. Read [Reading your results](doc_results.md) first.
    - **Before you start:** A finished run.

---

## Output directory structure

A finished run writes to `<path.sim>/sub-XXX/`:

```
sub-XXX_<medium>_report<affix>.html       ← self-contained HTML report
sub-XXX_<medium><affix>.csv               ← quantitative metrics
sub-XXX_<medium><affix>_targeting.csv     ← targeting accuracy (layered or phantom)
sub-XXX_placement<affix>.json             ← resolved transducer positions (T1 voxels; every placement mode)

nii/                                       ← NIfTI maps (medium, acoustic, thermal)

cache/                                     ← regenerable intermediates
    sub-XXX_<medium>_parameters_<yyMMdd_HHmm><affix>.mat   parameters snapshot at start
    sub-XXX_<medium>_report_parameters<affix>.mat          final parameters (report re-rendering)
    sub-XXX_<medium>_cropped_smoothed<affix>.mat           cropped/smoothed head (layered)
    sub-XXX_<medium>_grid_cache<affix>.mat                 grid stage (io.save_grid_cache)
    sub-XXX_<medium>_planimg<affix>.mat                    planning-image sidecar (always)
    sub-XXX_<medium>_medium_cache<affix>.mat               medium stage (io.save_medium_cache)
    sub-XXX_<medium>_kwave_source<affix>.mat               source (io.save_source_matrices)
    sub-XXX_<medium>_results<affix>.mat                    acoustic results + provenance
    sub-XXX_<medium>_heating_res<affix>.mat                thermal results

log/
    sub-XXX_<medium>_<timestamp><affix>.txt                MATLAB diary of the run
    sub-XXX_<medium>_parameters<affix>.yaml                resolved parameters (human-readable)

debug/                                     ← diagnostic artefacts (debug=1 only)
    preproc/        head preprocessing (rotation, cropping, skull visualisations)
    medium/         medium mapping (grid-space property matrices, pCT)
    source/         source/transducer setup (element distribution plots)
```

- **Cache files** are regenerable. With `io.overwrite_files = 'never'`, existing files are reloaded instead of recomputed ([Overwrite behaviour](#overwrite-behaviour)). `io.save_matrices = 0` writes no cache files; `1` saves all of them. Per-step flags such as `io.save_acoustic_matrices` override it when set.
- **Debug files** are written only with `parameters.simulation.debug = 1`.
- **Run statistics**: each run with an acoustic or thermal simulation adds one line to `<path.sim>/run_stats.jsonl` (one file per `path.sim`, shared by all subjects): grid size, time steps, thermal steps, simulation durations, solver, partition and peak RAM. The Nexus estimates resources from it ([Resource estimates](doc_hpc.md#resource-estimates-from-earlier-runs)).

---

## Summary table

Filename: `sub-XXX_<simulation_medium><affix>.csv`

One CSV holds all quantitative metrics: acoustic metrics first, then thermal metrics if heating was simulated. Tissue-specific columns exist for layered simulations only.

Targeting accuracy (`modules.run_targeting_analysis`, layered or phantom media) goes to `sub-XXX_<simulation_medium><affix>_targeting.csv`; see [Targeting accuracy](doc_group.md#targeting-accuracy).

### Intensity nomenclature

All intensity values are **pulse-average intensity (IPA)**, computed from the steady-state peak pressure amplitude:

```
I = p² / (2ρc)   [W/cm²]
```

`p` is the peak pressure amplitude, `ρ` density and `c` sound speed. The factor 2 averages over the sinusoid's cycle. The result is the time-average intensity *during* the pulse (IPA). It is not the temporal peak (ITP, without the factor 2), nor the average over the repetition interval (ISPTA, which also applies the duty cycle PD/PRI).

The **spatial peak** of IPA is ISPPA. IPA at a specific location (target point, radius average) is labelled IPA.

The **mechanical index** per voxel is MI = p⁻ / √f, with p⁻ the peak rarefaction (negative) pressure in MPa and f the centre frequency in MHz (IEC 62359 / FDA). PRESTUS uses the absolute value of `p_min_all` from k-Wave. For multi-transducer async combinations it uses the per-voxel most-negative pressure across transducers (`p_min_async`).

[Darmani et al. (2022, Clinical Neurophysiology)](https://doi.org/10.1016/j.clinph.2021.12.010) give an overview of the nomenclature.

![intensity_nomenclature](https://ars.els-cdn.com/content/image/1-s2.0-S1388245721008920-gr2.jpg)

??? note "References"
    Darmani, G., Bergmann, T. O., Butts Pauly, K., Caskey, C. F., de Lecea, L., Fomenko, A., Fouragnan, E., Legon, W., Murphy, K. R., Nandi, T., Phipps, M. A., Pinton, G., Ramezanpour, H., Sallet, J., Yaakub, S. N., Yoo, S. S., & Chen, R. (2022). Non‑invasive transcranial ultrasound stimulation for neuromodulation. *Clinical Neurophysiology, 133*(7), 167–177. https://doi.org/10.1016/j.clinph.2022.03.011

### Acoustic metrics

#### Global

| Column | Unit | Description |
|---|---|---|
| `subject_id` | none | Subject identifier |
| `freq_Hz` | Hz | Transducer centre frequency |
| `Isppa` | W/cm² | ISPPA: maximum pulse-average intensity across the entire simulation grid |
| `Isppa_after_exitplane` | W/cm² | Maximum IPA beyond the transducer exit plane (or within brain if modelled) |
| `real_focal_distance_mm` | mm | Distance from the transducer face to the location of maximum IPA |
| `Ipa_target` | W/cm² | IPA at the specified target coordinate |
| `Ipa_target_radius` | W/cm² | Mean IPA within a sphere of radius `analysis.focus_area_radius` mm centred on the target |
| `trans_pos_vox` | voxels | Transducer position in simulation grid coordinates |
| `focus_pos_vox` | voxels | Target focus position in simulation grid coordinates |

Water/free-field simulations also include `Psptp` (spatial peak temporal peak pressure), `Ptp_target` (temporal peak pressure at the target) and `MI` (global mechanical index from peak rarefaction pressure, unitless; ITRUSST limit 1.9, labelled "MI (free water)" in the report). They also include the ITRUSST focal metrics (`focal_vol_3dB_mm3`, `focal_vol_6dB_mm3`, `focal_centroid_3dB_dist_mm`), computed beyond the transducer exit plane.

#### Tissue-specific *(layered simulations only)*

| Column | Unit | Description |
|---|---|---|
| `Isppa_brain` | W/cm² | Maximum IPA within brain tissue (GM + WM) |
| `Isppa_skull` | W/cm² | Maximum IPA within skull bone |
| `Isppa_skin` | W/cm² | Maximum IPA within skin |
| `Psptp_brain` | Pa | Spatial peak temporal peak pressure within brain tissue |
| `Psptp_skull` | Pa | Spatial peak temporal peak pressure within skull bone |
| `Psptp_skin` | Pa | Spatial peak temporal peak pressure within skin |
| `Psptp` | Pa | Global spatial peak temporal peak pressure across the entire grid (also written for water/free-field runs). Max pressure safety limit: 2 MPa |
| `Ptp_target` | Pa | Temporal peak pressure at the target coordinate |
| `MI_brain` | none | Maximum mechanical index within brain tissue (GM + WM) |
| `MI_skull` | none | Maximum mechanical index within skull bone |
| `MI_skin` | none | Maximum mechanical index within skin |
| `MI_tc` | none | **Transcranial MI**: maximum mechanical index across all intracranial voxels (WM, GM, CSF, blood; skull bone excluded). Primary safety metric for cavitation risk at the sonication target. ITRUSST limit: 1.9 |
| `Ix_brain_vox`, `Iy_brain_vox`, `Iz_brain_vox` | voxels | Grid coordinates of peak IPA within brain |
| `halfmax_ISPPA_volume_brain_mm3` | mm³ | Volume of the −3 dB intensity region within brain (IPA ≥ 50 % of `Isppa_brain`); equals `focal_vol_3dB_mm3` |
| `focal_vol_3dB_mm3` | mm³ | Volume of the −3 dB focal zone within brain (IPA ≥ 50 % of peak; 10·log₁₀(0.5) = −3 dB) |
| `focal_vol_6dB_mm3` | mm³ | Volume of the −6 dB focal zone within brain (IPA ≥ 25 % of peak = pressure ≥ 50 % of peak; FWHM) |
| `focal_centroid_3dB_dist_mm` | mm | Distance from the transducer to the geometric centroid of the −3 dB focal region |

### Thermal metrics *(only written when heating simulations are run)*

Temperature is simulated over the full sonication protocol (on/off cycles) with the Pennes bioheat equation. All values are spatial maxima within the respective tissue.

#### Global

| Column | Unit | Description |
|---|---|---|
| `maxT` | °C | Global maximum temperature reached at any point during sonication |
| `endT` | °C | Global maximum temperature at the end of the last sonication pulse |
| `maxCEM43` | min | Global maximum CEM43 (kWave form) over the full protocol |
| `maxCEM43end` | min | Global maximum CEM43 (kWave form) at end of last pulse |
| `maxCEM43iso` | min | Global maximum CEM43 (ISO form) over the full protocol |
| `maxCEM43isoend` | min | Global maximum CEM43 (ISO form) at end of last pulse |

#### Tissue-specific

| Column | Unit | Description |
|---|---|---|
| `maxT_brain/skull/skin` | °C | Maximum temperature in each tissue during sonication. ITRUSST limit: 39 °C |
| `endT_brain/skull/skin` | °C | Maximum temperature in each tissue at end of last pulse |
| `riseT_brain/skull/skin` | °C | Maximum temperature rise above baseline (maxT − baseline). ITRUSST limit: 2 °C |
| `rise_endT_brain/skull/skin` | °C | Temperature rise above baseline at end of last pulse |
| `CEM43_brain/skull/skin` | min | Maximum CEM43 (kWave form) in each tissue. ITRUSST limits: brain 2 min, skull 16 min, skin 21 min |
| `CEM43_end_brain/skull/skin` | min | CEM43 (kWave form) in each tissue at end of last pulse |
| `CEM43iso_brain/skull/skin` | min | Maximum CEM43 (ISO form) in each tissue |
| `CEM43iso_end_brain/skull/skin` | min | CEM43 (ISO form) in each tissue at end of last pulse |

Both CEM43 formulations are always computed; see [Thermal simulations](doc_simulations-thermal.md) for the formulas. The HTML report shows the one selected by `parameters.thermal.cem43_iso`.

### Non-significant risk limits

The HTML report flags each metric against the ITRUSST consensus non-significant risk limits for transcranial ultrasound (Aubry et al., 2025, Brain Stimulation):

| Metric | Limit |
|---|---|
| MI transcranial | 1.9 |
| MI (free water) | 1.9 |
| Max pressure | 2 MPa |
| Temperature rise (a) | 2 °C |
| Maximum temperature (b) | 39 °C |
| CEM43 (brain) (c) | 2 min |
| CEM43 (skull) (c) | 16 min |
| CEM43 (skin) (c) | 21 min |

For thermal estimates, non-significant risk can be characterised by a, b, or c. These limits are for quick reference. Always check the original reference and possible consensus updates. PRESTUS does not make independent recommendations.

---

## NIfTI images

PRESTUS exports 3D volumes (2D for axisymmetric simulations) in subject-native T1 space and MNI-152 space. Each group of maps is written as soon as its stage finishes.

- Subject T1 space: `sub-XXX_<medium>_T1w<affix>_<type>.nii.gz`
- MNI space: `sub-XXX_<medium>_MNI<affix>_<type>.nii.gz`

### Medium maps *(written after medium setup)*

| Type | Resampling | Description |
|---|---|---|
| `medium_masks` | nearest | Tissue label mask (integer indices) |
| `sound_speed` | nearest | Speed of sound [m/s]; edges filled with water value |
| `density` | nearest | Mass density [kg/m³]; edges filled with water value |
| `alpha_coeff` | nearest | Absorption coefficient [dB/MHz/cm]; edges filled with water value |

### Acoustic maps *(written after acoustic analysis)*

| Type | Description |
|---|---|
| `intensity` | Pulse-average intensity (IPA) [W/cm²] |
| `MI` | Mechanical index |
| `pressure` | Peak pressure [Pa]; with a JSON sidecar so a later run can start from it ([Starting from earlier stage outputs](#starting-from-earlier-stage-outputs)) |
| `focal_mask_3dB` | −3 dB focal zone binary mask (IPA ≥ 50 % of peak; restricted to brain for layered, to beyond-exit-plane for free-field) |
| `focal_mask_6dB` | −6 dB / FWHM focal zone binary mask (IPA ≥ 25 % of peak = P ≥ Pmax/2; same region restriction) |

### Thermal maps *(written after thermal analysis)*

| Type | Description |
|---|---|
| `heating` | Maximum temperature during sonication [°C] |
| `heating_end` | Temperature at end of last pulse [°C] |
| `heatrise` | Temperature rise above baseline [°C] |
| `heatrise_end` | Temperature rise above baseline at end of last pulse [°C] |
| `CEM43` | CEM43 thermal dose, kWave form [min] |
| `CEM43_end` | CEM43 (kWave) at end of last pulse [min] |
| `CEM43_iso` | CEM43 thermal dose, ISO form [min] |
| `CEM43_iso_end` | CEM43 (ISO) at end of last pulse [min] |

---

## Figures

PRESTUS saves PNG overlays of intensity and temperature on the segmentation.

![PRESTUS_fig_example_thermal](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/thermal_fig_examples.png)

---

## HTML report

`sub-XXX_<medium>_report<affix>.html` is a self-contained report in the output directory. Sections:

| Section | Content |
|---|---|
| **Safety dashboard** | colour-coded cards for all ITRUSST-limited metrics (green / amber / red); thermal cards use the formulation selected by `parameters.thermal.cem43_iso` |
| **Simulation summary** | key acoustic and thermal values |
| **Methods** | ITRUSST-aligned methods text (CC0) with reference list: software, transducer, grid, tissue properties, pulse timing, thermal model, exposure values. *Copy methods text* / *Copy text + references* copy it as formatted text with a plain-text fallback. Check it against your study before use |
| **Configuration summary** | parameters used for the run |
| **Medium properties** | per-tissue acoustic and thermal properties (layered only) |
| **Acoustic results** | intensity overlays and the full acoustic metrics table |
| **Thermal results** | temperature, CEM43 (kWave) and CEM43 (ISO) plots (if heating was simulated) |
| **Post-hoc water simulation** | free-field reference results (if run) |
| **Log** | full pipeline execution log |

![PRESTUS_html](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/html_report.png)

Example report on a phantom: [Download example HTML](https://raw.githubusercontent.com/jkosciessa/PRESTUS_bin/main/examples/example_report.html)

---

## Where outputs live

Each location belongs to one level. A level only reads the levels above it, and each folder has one writer.

| Level | Depends on | Locations | Written by |
|---|---|---|---|
| **Inputs** | none | `path.anat` (T1/T2), `path.localite_raw` (raw Localite exports), `path.localite_t1` (planning T1s) | you / the scanner / Localite; PRESTUS only reads them |
| **Subject derivatives** | the subject's images | `path.seg` (SimNIBS `m2m_sub-XXX`), `path.pct` (pseudo-CT), `path.localite_post` (ingested Localite markers) | segmentation (layered media; `io.overwrite_simnibs`), pseudo-CT generation, `neuronav_ingest_markers` |
| **Placement** | subject derivatives, target, transducer geometry, placement settings | `path.placement` (default `<path.sim>/placement`): `sub-XXX/<mode>/<target_name>/` | the placement stage (heuristic, PlanTUS, TUS_entry) |
| **Simulation** | everything above, plus the acoustic/thermal setup | `path.sim/sub-XXX/` (this page) | the simulation pipeline |

- **Placement** does not depend on the medium, grid or drive level. Point several setups at one `path.placement` to reuse a heuristic, PlanTUS or TUS_entry result. Each run folder holds everything the run produced and a run record, `placement_run.json`. `manual`, `mni` and `localite` placements have no run folder. The simulation records the positions it used in `sub-XXX_placement<affix>.json`.
- **`path.localite_export`** is a hand-off folder, not a level. It holds copies of placement results for import into Localite (heuristic position table, Localite-ready T1).
- **SimNIBS `m2m_sub-XXX`** folders are only read after segmentation. For the pseudo-CT location, see [(pseudo-)CT](doc_pseudoCT.md).

---

## Starting from earlier stage outputs

PRESTUS keeps two kinds of intermediate results apart:

- **Within-run caches** (`cache/` of a simulation folder) are reloaded only by the run that wrote them, or a re-run of it. `io.overwrite_files` decides reuse ([Overwrite behaviour](#overwrite-behaviour)).
- **Stage outputs** are complete results that a later, separate run starts from. Each has a BIDS-style JSON sidecar whose `signature` lists the settings the output depends on. The later run compares it with its own settings.

| Stage output | Sidecar | Picked up by | Compared settings |
|---|---|---|---|
| SimNIBS segmentation, `path.seg/m2m_sub-XXX/` | none (standard SimNIBS output) | every layered run, when `final_tissues.nii.gz` exists (`io.overwrite_simnibs`) | none |
| Pseudo-CT, `<pct folder>/pseudoCT.nii.gz` | `pseudoCT.json` | every run with `pct.enabled: 1`, when it exists | `pct.skull_mapping` (with `custom_slope` / `custom_intercept`) |
| Placement run, `<path.placement>/sub-XXX/<mode>/<target_name>/` | `placement_run.json` | heuristic, PlanTUS and TUS_entry placement, per `placement.reuse` | `placement.<mode>` (without tool install paths) and the transducer without its positions |
| Pressure map, `nii/T1w/sub-XXX_<medium>_T1w<affix>_pressure.nii.gz` | `…_pressure.json` | a run with `path.acoustic_input` set to the NIfTI | `simulation.medium`, `grid` (without `thermal_*`), `layers`, `headmodel`, `pct`, `medium_properties`, the transducer (geometry and grid positions) and the subject |

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

??? note "Details: sidecar format"
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

## Overwrite behaviour

`io.overwrite_files` controls whether existing files are reused or recomputed.

| Value | Single run | Chained pipelines: intermediate stages | Chained pipelines: output stages |
|---|---|---|---|
| `never` | Reuse every existing file | Resume from cache | Skip if the output exists |
| `always` (default) | Recompute everything | Resume from cache | Always re-run, overwrite outputs |
| `force` | Same as `always` | Always re-run | Always re-run, overwrite outputs |
| `ask` | Prompt per file (interactive only) | Resume from cache | Skip if the output exists |

Under `never` or `ask`, an acoustic or thermal stage that recomputes its results (no cache to reuse, e.g. with `io.save_matrices: 0`) still replaces its existing NIfTI outputs, so the NIfTI and CSV outputs of a folder always come from the same run. `io.overwrite_simnibs` governs SimNIBS segmentations separately.

Use `force` when an input of the intermediate caches has changed: transducer placement, transducer or source parameters, grid, segmentation, or medium properties. Grid and medium caches detect such changes ([Cache reuse checks](#cache-reuse-checks)); the acoustic cache only detects a changed grid size. Alternatively, delete the relevant `cache/` files before resubmitting.

??? note "Details: chained pipelines"
    Chained pipelines split one run into separate jobs. Each job reuses what the previous stage wrote:

    | Pipeline | Intermediate stages | Output stages |
    |---|---|---|
    | Async multi-transducer ([doc](doc_multitransducer.md)) | Stage-0 preproc, per-transducer acoustic, combine | Joint thermal + report, water validation |
    | Multi-ISPPA ([doc](doc_multiintensity.md)) | Stage-1 acoustic | Per-ISPPA thermal + report, summary report |
    | Uncertainty ([doc](doc_uncertainty.md)) | Stage-1 preprocessing | Variant simulations + reports, uncertainty report |

    With the default `always`, resubmitting a finished chained run reuses the acoustic simulations and regenerates the CSV, NIfTI and HTML outputs. Within a chained run, sub-jobs that consume a shared cache (e.g. per-transducer acoustic jobs loading the Stage-0 grid) always load it, even under `force`. Only the stage that owns the cache rebuilds it.

## Cache reuse checks

Grid, medium and acoustic caches are named by subject, medium and affix only. Each cache therefore stores the settings it depends on (`cache_signature`). When a cache is reused (under `io.overwrite_files: never`, after a declined overwrite prompt, or in the intermediate stages of a chained pipeline under `always`) and the settings differ, `io.cache_signature_check` decides:

| Value | Effect |
|---|---|
| `error` (default) | Stop with `prestus:cache:signatureMismatch`, naming the differing keys. |
| `warn` | Warn and reuse the cache anyway. |
| `off` | No check. |

A cache without a signature is reused with the warning `prestus:cache:unverified`. Stage-output sidecars are checked the same way ([Starting from earlier stage outputs](#starting-from-earlier-stage-outputs)).

??? note "Details: what each cache compares"
    | Cache | Compared settings |
    |---|---|
    | Grid (`*_grid_cache*.mat`) | `simulation.medium`, `grid`, `layers`, `headmodel`, `pct`, and the transducer (geometry and positions after placement) |
    | Medium (`*_medium_cache*.mat`) | the grid settings above, plus `medium_properties` |
    | Acoustic (`*_results*.mat`) | the medium settings above, plus `subject_id` and `grid.acoustic_fov_*`; not `grid.thermal_*` |

    For `simulation.medium: phantom`, each also compares the phantom file (size and modification date of `final_tissues.nii.gz`).

    Not compared, because they do not change these caches: paths, I/O, HPC and module flags; the drive level (`target_isppa_wcm2`, `tx_calibration`), which only scales the acoustic field; and, for grid and medium caches, the partial acoustic FOV (`grid.acoustic_fov_*`), which is cropped after the grid and medium stages. A run that reuses a cache keeps its own values of these settings.

    Sub-jobs that load a cache owned by another stage of the same chained run are not checked.

---

## Group outputs

The group tools read existing per-run outputs and never simulate; see [Group analysis](doc_group.md).

**Study** (`prestus_group_report_start` on a config with `study.manifest`, or the end of `prestus_study_start`), written to `<path.sim>/group/<study.name>/`:

```
group_<name>_report.html                  ← self-contained group HTML report
group_<name>_runs.csv                     ← one row per run: settings, metrics, targeting
plots/
    group_<medium>_<metric>_<name>.png      box plots (intensity, mi, conditions, targeting, maxT, riseT, cem43)
    group_<cond>_focal_overlap.png          MNI overlap figure per condition
maps/
    group_<cond>_focal_overlap_3dB.nii.gz   fraction of runs covering each voxel
    group_<cond>_<Ipa|Itp>_mean.nii.gz      mean intensity (and _sd)
    group_<cond>_heatrise_mean.nii.gz       mean temperature rise (heating runs)
log/
    plan_<time>.tsv / launch_<time>.tsv     study dry-run plan / launch log
log_hpc/                                  ← dependent study-report job (HPC)
```

**Subjects** (no manifest: all subjects with the same `simulation.medium` and `io.output_affix`; `prestus_group_report_start`, or `modules.generate_group_report = 1` after each subject), written to `path.sim`:

```
group_<medium>_report<affix>.html         ← self-contained group HTML report
group_<medium>_runs<affix>.csv            ← one row per subject
group_plots/
    group_<medium>_<metric><affix>.png      box plots: intensity (Isppa / Ipa_target), mi (per tissue),
                                            targeting (MNI centroids, focus error), maxT, riseT, cem43
group/<medium><affix>/maps/, plots/       ← group MNI maps (manual prestus_group_report_start only)
```

Thermal plots (`maxT`, `riseT`, `cem43`) need a layered medium with heating results. See [Group HTML report](doc_group.md#group-html-report) for the report contents.
