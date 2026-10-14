# Simulation outputs

PRESTUS provides multiple outputs. Some of these outputs are optional (or can be deactivated upon request to save space.)

---

## Summary table

Filename: `sub-XXX_<simulation_medium><affix>.csv`

Targeting accuracy (`modules.run_targeting_analysis`, layered or phantom media) goes to a sidecar, `sub-XXX_<simulation_medium><affix>_targeting.csv`; see [Targeting accuracy](doc_group.md#targeting-accuracy).

All quantitative metrics are written to a single CSV file. Acoustic metrics are written first; thermal metrics are appended if heating simulations were run. Tissue-specific columns are only present for layered simulations.

### Intensity nomenclature

All intensity values are **pulse-average intensity (IPA)**, computed from the steady-state peak pressure amplitude as:

```
I = p² / (2ρc)   [W/cm²]
```

where `p` is the peak pressure amplitude, `ρ` is density, and `c` is sound speed. The cycle-averaging factor of 2 reflects that `p` is the amplitude of a sinusoid — this yields the time-average intensity *during* the pulse, i.e. IPA, not the temporal peak (ITP, which would lack the factor of 2) and not the temporal average over the full repetition interval (ISPTA, which would additionally apply the duty cycle PD/PRI).

The **spatial peak** of IPA is ISPPA. When IPA is reported at a specific location (target point, radius average) rather than at the spatial maximum, it is labelled IPA.

The **mechanical index** at each voxel is defined as MI = p⁻ / √f, where p⁻ is the peak rarefaction (negative) pressure in MPa and f is the centre frequency in MHz (IEC 62359 / FDA). PRESTUS records `p_min_all` from k-Wave and uses its absolute value for MI; for multi-transducer async combinations the per-voxel most-negative pressure across transducers (`p_min_async`) is used instead.

An overview of the nomenclature is provided in [Darmani et al. (2022, Clinical Neurophysiology)](https://doi.org/10.1016/j.clinph.2021.12.010).

![intensity_nomenclature](https://ars.els-cdn.com/content/image/1-s2.0-S1388245721008920-gr2.jpg)

<br>
<span style="font-size: 12px;">
  <strong>References:</strong><br>
  Darmani, G., Bergmann, T. O., Butts Pauly, K., Caskey, C. F., de Lecea, L., Fomenko, A., Fouragnan, E., Legon, W., Murphy, K. R., Nandi, T., Phipps, M. A., Pinton, G., Ramezanpour, H., Sallet, J., Yaakub, S. N., Yoo, S. S., & Chen, R. (2022). Non‑invasive transcranial ultrasound stimulation for neuromodulation. <em>Clinical Neurophysiology, 133</em>(7), 167–177. https://doi.org/10.1016/j.clinph.2022.03.011
</span>
<br>

### Acoustic metrics

#### Global

| Column | Unit | Description |
|---|---|---|
| `subject_id` | — | Subject identifier |
| `freq_Hz` | Hz | Transducer centre frequency |
| `Isppa` | W/cm² | ISPPA — maximum pulse-average intensity across the entire simulation grid |
| `Isppa_after_exitplane` | W/cm² | Maximum IPA beyond the transducer exit plane (or within brain if modelled) |
| `real_focal_distance_mm` | mm | Distance from the transducer face to the location of maximum IPA |
| `Ipa_target` | W/cm² | IPA at the specified target coordinate |
| `Ipa_target_radius` | W/cm² | Mean IPA within a sphere of radius `analysis.focus_area_radius` mm centred on the target |
| `trans_pos_vox` | voxels | Transducer position in simulation grid coordinates |
| `focus_pos_vox` | voxels | Target focus position in simulation grid coordinates |

Water/free-field simulations additionally include `Psptp` (spatial peak temporal peak pressure), `Ptp_target` (temporal peak pressure at the target), `MI` (global mechanical index from peak rarefaction pressure, unitless; ITRUSST limit 1.9, labelled "MI (free water)" in the report), and the full set of ITRUSST focal metrics (`focal_vol_3dB_mm3`, `focal_vol_6dB_mm3`, `focal_centroid_3dB_dist_mm`) computed over the region beyond the transducer exit plane.

#### Tissue-specific *(layered simulations only)*

| Column | Unit | Description |
|---|---|---|
| `Isppa_brain` | W/cm² | Maximum IPA within brain tissue (GM + WM) |
| `Isppa_skull` | W/cm² | Maximum IPA within skull bone |
| `Isppa_skin` | W/cm² | Maximum IPA within skin |
| `Psptp_brain` | Pa | Spatial peak temporal peak pressure within brain tissue |
| `Psptp_skull` | Pa | Spatial peak temporal peak pressure within skull bone |
| `Psptp_skin` | Pa | Spatial peak temporal peak pressure within skin |
| `Psptp` | Pa | Global spatial peak temporal peak pressure across the entire grid (same column also written for water/free-field runs). Max pressure safety limit: 2 MPa |
| `Ptp_target` | Pa | Temporal peak pressure at the target coordinate |
| `MI_brain` | — | Maximum mechanical index within brain tissue (GM + WM) |
| `MI_skull` | — | Maximum mechanical index within skull bone |
| `MI_skin` | — | Maximum mechanical index within skin |
| `MI_tc` | — | **Transcranial MI** — maximum mechanical index across all intracranial voxels (WM, GM, CSF, blood; skull bone excluded). Primary safety metric for cavitation risk at the sonication target. ITRUSST limit: 1.9 |
| `Ix_brain_vox`, `Iy_brain_vox`, `Iz_brain_vox` | voxels | Grid coordinates of peak IPA within brain |
| `halfmax_ISPPA_volume_brain_mm3` | mm³ | Volume of the −3 dB intensity region within brain (IPA ≥ 50 % of `Isppa_brain`; preserved for backward compatibility — equals `focal_vol_3dB_mm3`) |
| `focal_vol_3dB_mm3` | mm³ | Volume of the −3 dB focal zone within brain (IPA ≥ 50 % of peak; 10·log₁₀(0.5) = −3 dB) |
| `focal_vol_6dB_mm3` | mm³ | Volume of the −6 dB focal zone within brain (IPA ≥ 25 % of peak = pressure ≥ 50 % of peak; FWHM) |
| `focal_centroid_3dB_dist_mm` | mm | Distance from the transducer to the geometric centroid of the −3 dB focal region |

---

### Thermal metrics *(only written when heating simulations are run)*

Temperature is simulated over the full sonication protocol (on/off cycles) using the Pennes bioheat equation. All values are spatial maxima within the respective tissue.

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

Both CEM43 formulations are always computed. See [doc_simulations-thermal.md](doc_simulations-thermal.md) for formula details. The HTML report displays the formulation selected by `parameters.thermal.cem43_iso`.

---

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

For thermal estimates, non-significant risk can be characterized by a, b, or c. NSR limits above and in the HTML are provided for quick reference. Always check the original reference for details and check for possible consensus updates. PRESTUS does not make independent recommendations.

---

## NIfTI images

PRESTUS exports 3D volumes (or 2D for axisymmetric simulations) in NIfTI format, provided in subject-native T1 space and MNI-152 space. NIfTI export is staged: medium property maps are written after medium setup, acoustic maps after acoustic analysis, and thermal maps after thermal analysis.

Filename patterns:
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
| `pressure` | Peak pressure [Pa] |
| `focal_mask_3dB` | −3 dB focal zone binary mask (IPA ≥ 50 % of peak; restricted to brain for layered, to beyond-exit-plane for free-field) |
| `focal_mask_6dB` | −6 dB / FWHM focal zone binary mask (IPA ≥ 25 % of peak = P ≥ Pmax/2; same region restriction) |

### Thermal maps *(written after thermal analysis)*

| Type | Description |
|---|---|
| `heating` | Maximum temperature during sonication [°C] |
| `heating_end` | Temperature at end of last pulse [°C] |
| `heatrise` | Temperature rise above baseline [°C] |
| `heatrise_end` | Temperature rise above baseline at end of last pulse [°C] |
| `CEM43` | CEM43 thermal dose — kWave form [min] |
| `CEM43_end` | CEM43 (kWave) at end of last pulse [min] |
| `CEM43_iso` | CEM43 thermal dose — ISO form [min] |
| `CEM43_iso_end` | CEM43 (ISO) at end of last pulse [min] |

---

## Output directory structure

```
sub-XXX_<medium>_report<affix>.html       ← self-contained HTML report
sub-XXX_<medium><affix>.csv               ← quantitative metrics
sub-XXX_<medium><affix>_targeting.csv     ← targeting accuracy (layered or phantom)
sub-XXX_placement<affix>.json             ← resolved transducer positions (T1 voxels; every placement mode)

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

**Cache files** are regenerable: if `io.overwrite_files = 'never'`, existing files are reloaded instead of recomputed. Set `io.save_matrices = 0` to skip saving intermediate matrices (cache files will not be written), or `1` to save all of them; per-step flags such as `io.save_acoustic_matrices` override it when set. See [Overwrite behaviour](#overwrite-behaviour) below.

**Debug files** are only written when `parameters.simulation.debug = 1`.

**Run statistics**: each run that ran an acoustic or thermal simulation adds one line to `<path.sim>/run_stats.jsonl` (one file per `path.sim`, shared by all subjects). The line records the grid size, time steps, thermal steps, simulation durations, solver, partition and peak RAM. The Nexus estimates resources from this file; see [Resource estimates](doc_hpc.md#resource-estimates-from-earlier-runs).

### Overwrite behaviour

`io.overwrite_files` controls whether existing files are reused or recomputed.

| Value | Single run | Chained pipelines: intermediate stages | Chained pipelines: output stages |
|---|---|---|---|
| `never` | Reuse every existing file | Resume from cache | Skip if the output exists |
| `always` (default) | Recompute everything | Resume from cache | Always re-run, overwrite outputs |
| `force` | Same as `always` | Always re-run | Always re-run, overwrite outputs |
| `ask` | Prompt per file (interactive only) | Resume from cache | Skip if the output exists |

SimNIBS segmentations are governed separately by `io.overwrite_simnibs`.

**Chained pipelines** split one run into separate jobs, each of which reuses what the previous stage wrote:

| Pipeline | Intermediate stages | Output stages |
|---|---|---|
| Async multi-transducer ([doc](doc_multitransducer.md)) | Stage-0 preproc, per-transducer acoustic, combine | Joint thermal + report, water validation |
| Multi-ISPPA ([doc](doc_multiintensity.md)) | Stage-1 acoustic | Per-ISPPA thermal + report, summary report |
| Uncertainty ([doc](doc_uncertainty.md)) | Stage-1 preprocessing | Variant simulations + reports, uncertainty report |

So with the default `always`, resubmitting a finished chained run reuses the expensive acoustic simulations and regenerates the CSV, NIfTI and HTML outputs. Use `force` when an input that feeds the intermediate caches has changed — transducer placement, transducer or source parameters, grid, segmentation, or medium properties. Grid and medium caches detect such changes (see below); the acoustic cache only detects a changed grid size, so after changing source or medium parameters use `force` or delete the relevant `cache/` files before resubmitting.

Within a chained run, sub-jobs that consume a shared cache (e.g. per-transducer acoustic jobs loading the Stage-0 grid) always load it rather than rebuilding it, even under `force`; the stage that owns the cache is the one that rebuilds it. Internally this is signalled by `io.reuse_upstream_caches` (set by the pipelines; not a user setting).

### Cache reuse checks

Grid and medium caches are named by subject, medium and affix only. So that a cache built from different settings is not reused silently, each cache stores the settings it depends on (`cache_signature`):

| Cache | Compared settings |
|---|---|
| Grid (`*_grid_cache*.mat`) | `simulation.medium`, `grid`, `layers`, `headmodel`, `pct`, and the transducer (geometry and positions after placement) |
| Medium (`*_medium_cache*.mat`) | the grid settings above, plus `medium_properties` |

Paths, I/O, HPC and module flags are not compared (moving data or changing downstream stages does not change these caches), nor is the drive level (`target_isppa_wcm2`, `tx_calibration`), which only scales the acoustic field.

When a cache is reused because of `io.overwrite_files: never` (or a declined overwrite prompt, or the intermediate stages of a chained pipeline under `always`) and the settings differ, `io.cache_signature_check` decides:

| Value | Effect |
|---|---|
| `error` (default) | Stop with `prestus:cache:signatureMismatch`, naming the differing keys. |
| `warn` | Warn and reuse the cache anyway. |
| `off` | No check. |

Caches written before signatures existed cannot be checked; they are reused with the warning `prestus:cache:unverified`. Sub-jobs that load a cache owned by another stage of the same chained run (`io.reuse_upstream_caches`) are not checked.

---

## Figures

PRESTUS saves PNG overlays of intensity and temperature on the segmentation for quick visual inspection.

![PRESTUS_fig_example_thermal](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/thermal_fig_examples.png)

---

## HTML report

A self-contained HTML report is written to the output directory (`sub-XXX_<medium>_report<affix>.html`). It includes:

- **Safety dashboard** — color-coded cards for all ITRUSST-limited metrics (green / amber / red)
- **Simulation summary** — key acoustic and thermal values at a glance
- **Methods** — auto-generated, ITRUSST-aligned methods text (CC0) covering software, transducer, simulation grid, tissue properties, pulse timing, thermal model and the simulated exposure values, with a matching reference list. *Copy methods text* / *Copy text + references* put it on the clipboard as formatted text (sub-/superscripts kept) with a plain-text fallback, ready to paste into a manuscript. Check it against your study before use.
- **Configuration summary** — parameters used for the run
- **Medium properties** — per-tissue acoustic and thermal properties (layered only)
- **Acoustic results** — intensity overlay images and the full acoustic metrics table
- **Thermal results** — temperature, CEM43 (kWave), and CEM43 (ISO) plots (if heating simulations were run); the safety dashboard uses the formulation selected by `parameters.thermal.cem43_iso`
- **Post-hoc water simulation** — free-field reference results (if run)
- **Log** — full pipeline execution log

![PRESTUS_html](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/html_report.png)

A full example report on a phantom can be found here: [🔗 Download Example HTML](https://raw.githubusercontent.com/jkosciessa/PRESTUS_bin/main/examples/example_report.html)

---

## Group outputs

The group tools read existing per-run outputs and never simulate; see [Group analysis](doc_group.md).

**Study** (`prestus_group_report_start` on a config with `study.manifest`, or the end of `prestus_study_start`): everything goes to `<path.sim>/group/<study.name>/`:

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

**Subjects** (no manifest: all subjects under the same `simulation.medium` and `io.output_affix`; `prestus_group_report_start`, or `modules.generate_group_report = 1` after each subject), written to `path.sim`:

```
group_<medium>_report<affix>.html         ← self-contained group HTML report
group_<medium>_runs<affix>.csv            ← one row per subject
group_plots/                              ← group-level box-plot images
    group_<medium>_intensity<affix>.png     Isppa / Ipa_target across subjects
    group_<medium>_mi<affix>.png            Mechanical Index per tissue
    group_<medium>_targeting<affix>.png     MNI centroids and focus error (targeting sidecars)
    group_<medium>_maxT<affix>.png          maximum temperature  (heating runs)
    group_<medium>_riseT<affix>.png         temperature rise     (heating runs)
    group_<medium>_cem43<affix>.png         CEM43 thermal dose   (heating runs)
group/<medium><affix>/maps/, plots/       ← group MNI maps (manual prestus_group_report_start only)
```

Thermal plots need a layered medium with heating results. Box plots auto-scale to the number of plotted columns. The Safety dashboard's Mechanical grid always shows both `MI` (free water) and `MItc` tiles plus Max pressure; whichever MI variant has no data for the medium renders greyed out ("N/A"). See [Group HTML report](doc_group.md#group-html-report) for details.
