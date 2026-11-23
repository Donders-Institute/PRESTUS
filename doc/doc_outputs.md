# Simulation Outputs

![Pipeline steps; this page belongs to 7 Outputs](assets/figures/fig_stage-strip-outputs_diagram.svg){ .stage-strip }

A run writes its results as files: a metrics table, [NIfTI](https://nifti.nimh.nih.gov/) maps, figures, and an HTML report, plus caches and stage outputs that a later run can reuse instead of recomputing. Knowing which file holds a value lets you check that the numbers you report come from the settings you intended.

On this page:

- [Output Directory Structure](#output-directory-structure): the folders and files of a run.
- [Where Outputs Live](#where-outputs-live): which folder each stage writes and what it depends on.
- [NIfTI Images](#nifti-images): the medium, acoustic, and thermal maps in T1 and MNI space.
- [Figures](#figures): the PNG overlays of intensity and temperature.
- [HTML Report](#html-report): the sections of the per-run report.
- [Debug Outputs](#debug-outputs): the intermediate figures and volumes of debug mode, and the debug report.
- [Group Outputs](#group-outputs): the files of the group report and maps.

!!! note "At a glance"
    - **When you need it:** To find a specific file, or to see where a stage writes its results. Read [Reading your results](doc_results.md) first.
    - **Before you start:** A finished run.

Subpages:

- [Metrics Table Columns](doc_outputs_metrics-table.md): every column of the metrics CSV and the reference levels; for anyone reporting or checking a value.
- [Reusing Earlier Outputs and Caches](doc_outputs_reuse.md): starting a run from earlier stage outputs, overwrite behaviour, and cache checks; for anyone re-running, resuming, or chaining runs.

---

## Output Directory Structure

A finished run writes to `<path.sim>/sub-XXX/`:

```
sub-XXX_<medium>_report<affix>.html       ← self-contained HTML report
sub-XXX_layered_segmentation_report<affix>.html ← head-model postprocessing, step by step (layered; report.segmentation)
sub-XXX_<medium>_debug_report<affix>.html ← all debug outputs (simulation.debug; report.debug)
sub-XXX_<medium><affix>.csv               ← quantitative metrics
sub-XXX_<medium><affix>_targeting.csv     ← targeting accuracy (layered or phantom)
sub-XXX_placement<affix>.json             ← resolved transducer positions (T1 voxels; every placement mode)

nii/                                       ← NIfTI maps (medium, transducer model, acoustic, thermal); not created with modules.run_nifti_creation: 0
img/                                       ← PNG figures (acoustic, thermal, preprocessing)
babelbrain/                                ← BabelBrain inputs (io.export_babelbrain or acoustic_backend: babelbrain)

cache/                                     ← regenerable intermediates
    sub-XXX_<medium>_parameters_<yyMMdd_HHmm><affix>.mat   parameters snapshot at start
    sub-XXX_<medium>_report_parameters<affix>.mat          final parameters (report re-rendering)
    sub-XXX_<medium>_cropped_smoothed<affix>.mat           cropped/smoothed head (layered)
    sub-XXX_<medium>_segqc<affix>.mat                      head-model postprocessing records (segmentation report)
    sub-XXX_<medium>_grid_cache<affix>.mat                 grid stage (io.save_grid_cache)
    sub-XXX_<medium>_planimg<affix>.mat                    planning-image sidecar (always)
    sub-XXX_<medium>_medium_cache<affix>.mat               medium stage (io.save_medium_cache)
    sub-XXX_<medium>_kwave_source<affix>.mat               source (io.save_source_matrices)
    sub-XXX_<medium>_results<affix>.mat                    acoustic results + provenance
    sub-XXX_<medium>_heating_res<affix>.mat                thermal results (io.save_thermal_matrices)
    sub-XXX_<medium>_HEAT<affix>.mat                       thermal time courses (thermal analysis; not BabelBrain)

log/
    sub-XXX_<medium>_<timestamp><affix>.txt                MATLAB diary of the run
    sub-XXX_<medium>_parameters<affix>.yaml                resolved parameters (human-readable; io.prestus_version: version, commit, commit date, uncommitted changes)

debug/                                     ← diagnostic artefacts (simulation.debug: 1 only)
    manifest.jsonl  one line per debug file: stage, caption, writing function
    preproc/        head preprocessing (rotation, cropping, skull visualisations)
    medium/         medium mapping (grid-space property matrices, pCT)
    source/         source/transducer setup (element distribution plots)
    babelbrain/     BabelBrain FOV and domain-coverage diagnostics
```

- **Cache files** are regenerable. With `io.overwrite_files = 'never'`, existing files are reloaded instead of recomputed ([Overwrite behaviour](doc_outputs_reuse.md#overwrite-behaviour)). `io.save_matrices = 0` writes none of the optional cache matrices; `1` saves all of them. The parameter snapshots and the planning-image sidecar are always written. Per-step flags such as `io.save_acoustic_matrices` override it when set.
- **Debug files** are written only with `parameters.simulation.debug = 1` ([Debug outputs](#debug-outputs)).
- **Run statistics**: each run with an acoustic or thermal simulation adds one line to `<path.sim>/run_stats.jsonl` (one file per `path.sim`, shared by all subjects): grid size, time steps, thermal steps, simulation durations, solver, partition, and peak RAM. The Nexus estimates resources from it ([Resource estimates](doc_run-time.md#from-earlier-runs)).

---

## Where Outputs Live

Each location belongs to one level. A level only reads the levels above it, and each folder has one writer.

| Level | Depends on | Locations | Written by |
|---|---|---|---|
| **Inputs** | none | `path.anat` (T1/T2), `path.localite_raw` (raw [Localite](https://www.localite.de/) exports), `path.localite_t1` (planning T1s) | you / the scanner / Localite; PRESTUS only reads them |
| **Subject derivatives** | the subject's images | `path.seg` ([SimNIBS](https://simnibs.github.io/simnibs/build/html/index.html) `m2m_sub-XXX`), `path.pct` (pseudo-CT), `<parent of path.seg>/mni/sub-XXX/` (T1 in MNI space), `path.localite_post` (ingested Localite markers) | segmentation (layered media; `io.overwrite_simnibs`), pseudo-CT generation, `neuronav_ingest_markers` |
| **Placement** | subject derivatives, target, transducer geometry, placement settings | `path.placement` (default `<path.sim>/placement`): `sub-XXX/<mode>/<target_name>/` | the placement stage (heuristic, [PlanTUS](https://github.com/mlueckel/PlanTUS), [TUS_entry](https://github.com/CyrilAtkinson/TUS_entry)) |
| **Simulation** | everything above, plus the acoustic/thermal setup | `path.sim/sub-XXX/` (this page) | the simulation pipeline |

- **Placement** does not depend on the medium, grid, or drive level. Point several setups at one `path.placement` to reuse a heuristic, PlanTUS, or TUS_entry result. Each run folder holds everything the run produced and a run record, `placement_run.json`. `manual`, `mni`, and `localite` placements have no run folder. The simulation records the positions it used in `sub-XXX_placement<affix>.json`.
- **`path.localite_export`** is a hand-off folder, not a level. It holds copies of placement results for import into Localite (heuristic position table, Localite-ready T1).
- **SimNIBS `m2m_sub-XXX`** folders are only read after segmentation. For the pseudo-CT location, see [(pseudo-)CT](doc_pseudoCT.md).

---

## NIfTI Images

PRESTUS exports 3D volumes (2D for axisymmetric simulations) in subject-native T1 space and MNI-152 space. Each group of maps is written as soon as its stage finishes. MNI maps (layered media, `io.save_MNI`) are warped with the nonlinear SimNIBS charm transform (`subject2mni`) onto SimNIBS's `MNI152_T1_1mm.nii.gz`: FSL's 1 mm MNI152 template, BIDS/TemplateFlow `MNI152NLin6Asym` (identical to FSL's file in SimNIBS 4.1.0).

- Subject T1 space: `sub-XXX_<medium>_T1w<affix>_<type>.nii.gz`
- MNI space: `sub-XXX_<medium>_MNI<affix>_<type>.nii.gz`

??? note "Medium maps (written after medium setup)"
    | Type | Resampling | Description |
    |---|---|---|
    | `medium_masks` | nearest | Tissue label mask (integer indices) |

    With `io.save_property_maps: 1` (layered media), each acoustic and thermal property map is also written, to `nii/properties/<property>_T1w<affix>.nii.gz` (and `_MNI<affix>`), without the `sub-XXX_<medium>` prefix: `sound_speed`, `density`, `alpha_coeff`, `alpha_power`, `thermal_conductivity`, `specific_heat`, `perfusion_coeff`, `absorption_fraction`, and `absorption_np`.

??? note "Transducer model (written after medium setup)"
    | Type | Resampling | Description |
    |---|---|---|
    | `transducer` | none (built on the T1 grid) | Model of the planned transducer bowl(s), from the planned T1-voxel positions; the bowl's rear centre is at `trans_pos`. `1` = emitting face across the active aperture (smallest `elem_id_mm` to largest `elem_od_mm`), `0.5` = the display thickness behind it (`io.transducer_nii_back_mm`), `0` = elsewhere, including the central opening. The aperture edges follow the profile's mm limits, not whole voxels. Layered medium, annular transducers, `io.save_transducer_nii`. T1 space only. |

??? note "Acoustic maps (written after acoustic analysis)"
    | Type | Description |
    |---|---|
    | `Itp` | Temporal-peak intensity [W/cm²] |
    | `Ipa` | Pulse-average intensity (IPA) [W/cm²]; with `io.save_Ipa` and `timing.pd` |
    | `Ita` | Temporal-average intensity [W/cm²]; with `io.save_Ita` |
    | `MI` | Mechanical index |
    | `pressure` | Peak pressure [Pa]; with a JSON sidecar so a later run can start from it ([Starting from earlier stage outputs](doc_outputs_reuse.md#starting-from-earlier-stage-outputs)) |
    | `focal_mask_3dB` | −3 dB focal zone binary mask (IPA ≥ 50 % of peak; restricted to brain for layered, to beyond-exit-plane for free-field) |
    | `focal_mask_6dB` | −6 dB / FWHM focal zone binary mask (IPA ≥ 25 % of peak = P ≥ Pmax/2; same region restriction) |

??? note "Thermal maps (written after thermal analysis)"
    | Type | Description |
    |---|---|
    | `heating` | Maximum temperature during sonication [°C] |
    | `heating_end` | Temperature at the end of the simulated protocol (after the last pulse train and `timing.post_ptri_dur`) [°C] |
    | `heatrise` | Temperature rise above baseline [°C] |
    | `heatrise_end` | Temperature rise above baseline at the end of the simulated protocol (after the last pulse train and `timing.post_ptri_dur`) [°C] |
    | `CEM43` | CEM43 thermal dose, k-Wave form [min] |
    | `CEM43_end` | CEM43 (k-Wave) at the end of the simulated protocol (after the last pulse train and `timing.post_ptri_dur`) [min] |
    | `CEM43_iso` | CEM43 thermal dose, ISO form [min] |
    | `CEM43_iso_end` | CEM43 (ISO) at the end of the simulated protocol (after the last pulse train and `timing.post_ptri_dur`) [min] |

---

## Figures

PRESTUS saves PNG overlays of intensity and temperature on the segmentation.

![Thermal outputs: maximum temperature overlaid on the segmentation, and temperature over time in skull, skin, brain, and target in the focal plane](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/thermal_fig_examples.png)

---

## HTML Report

`sub-XXX_<medium>_report<affix>.html` is a self-contained report in the output directory. Sections:

| Section | Content |
|---|---|
| **Safety dashboard** | colour-coded cards for all ITRUSST-limited metrics (green / amber / red); thermal cards use the formulation selected by `parameters.thermal.cem43_iso` |
| **Simulation summary** | key acoustic and thermal values |
| **Methods** | ITRUSST-aligned methods text ([CC0](https://creativecommons.org/publicdomain/zero/1.0/)) with reference list: software, transducer, grid, tissue properties, pulse timing, thermal model, exposure values. *Copy methods text* / *Copy text + references* copy it as formatted text with a plain-text fallback. Check it against your study before use |
| **Configuration summary** | parameters used for the run |
| **Medium properties** | per-tissue acoustic and thermal properties (layered only) |
| **Acoustic results** | intensity overlays and the full acoustic metrics table |
| **Thermal results** | temperature, CEM43 (k-Wave), and CEM43 (ISO) plots (if heating was simulated) |
| **Post-hoc water simulation** | free-field reference results (if run) |
| **Diagnostics** | links to the segmentation report and the debug report |
| **Log** | full pipeline execution log |

![Simulation report: safety dashboard with MI, temperature rise, CEM43, and I<sub>SPPA</sub> per tissue against their limits, followed by the run header and summary](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/html_report.png)

Example report on a phantom: [Download example HTML](https://raw.githubusercontent.com/jkosciessa/PRESTUS_bin/main/examples/example_report.html)

---

## Debug Outputs

When a result looks wrong, the intermediate states of the run show where it went wrong. `simulation.debug: 1` writes them to `debug/`, one subfolder per stage, and `sub-XXX_<medium>_debug_report<affix>.html` collects them on one page: the figures of each stage with a caption and the function that wrote them, and a table of the volumes. If the run fails, the page is still written and opens with the error and its stack.

Finer settings choose what debug mode writes. `simulation.debug_stages` limits the outputs to some stages (`preproc`, `medium`, `source`, `babelbrain`), `simulation.debug_figures` and `simulation.debug_volumes` switch the figures (PNG) and the intermediate volumes (NIfTI) on or off, and `report.debug` the report. Each volume is about the size of the head volume, so `simulation.debug_volumes: 0` saves most of the disk space. See [simulation parameters](doc_parameters_simulation.md#simulation-type-simulation).

The report also links the [segmentation report](doc_head-model.md#segmentation-report), the SimNIBS `charm` outputs (`segmentation.debug` keeps more of them), the pseudo-CT folder (`pct.debug` keeps its intermediates), and the run logs. Figures are linked, not embedded, so keep the report next to its `debug/` folder when you move it.

??? note "Reading the debug volumes"
    The `*_final` volumes in `debug/preproc/` are back-projected to the T1 and overlay it in any NIfTI viewer. The `matrix_*` volumes in `debug/medium/`, the `balloon_mask*` volumes, and `debug/source/*_source_mask` are in simulation-grid voxels without a T1 header; open them alongside each other, not over the T1. The source mask shows where the transducers sit against the skin in `matrix_density`; with a partial acoustic FOV it is on the cropped grid and no longer overlays the medium matrices.

---

## Group Outputs

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
