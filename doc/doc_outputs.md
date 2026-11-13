# Simulation Outputs

A run writes its results as files: a metrics table, [NIfTI](https://nifti.nimh.nih.gov/) maps, figures and an HTML report, plus caches and stage outputs that a later run can reuse instead of recomputing. Knowing which file holds a value lets you check that the numbers you report come from the settings you intended.

On this page:

- [Output directory structure](#output-directory-structure): the folders and files of a run.
- [Where outputs live](#where-outputs-live): which folder each stage writes and what it depends on.
- [NIfTI images](#nifti-images): the medium, acoustic and thermal maps in T1 and MNI space.
- [Figures](#figures): the PNG overlays of intensity and temperature.
- [HTML report](#html-report): the sections of the per-run report.
- [Group outputs](#group-outputs): the files of the group report and maps.

!!! note "At a glance"
    - **When you need it:** To find a specific file, or to see where a stage writes its results. Read [Reading your results](doc_results.md) first.
    - **Before you start:** A finished run.

Subpages:

- [Metrics Table Columns](doc_outputs_metrics-table.md): every column of the metrics CSV and the reference levels; for anyone reporting or checking a value.
- [Reusing Earlier Outputs and Caches](doc_outputs_reuse.md): starting a run from earlier stage outputs, overwrite behaviour and cache checks; for anyone re-running, resuming or chaining runs.

---

## Output Directory Structure

A finished run writes to `<path.sim>/sub-XXX/`:

```
sub-XXX_<medium>_report<affix>.html       ← self-contained HTML report
sub-XXX_<medium><affix>.csv               ← quantitative metrics
sub-XXX_<medium><affix>_targeting.csv     ← targeting accuracy (layered or phantom)
sub-XXX_placement<affix>.json             ← resolved transducer positions (T1 voxels; every placement mode)

nii/                                       ← NIfTI maps (medium, acoustic, thermal); not created with modules.run_nifti_creation: 0
img/                                       ← PNG figures (acoustic, thermal, preprocessing)
babelbrain/                                ← BabelBrain inputs (io.export_babelbrain or acoustic_backend: babelbrain)

cache/                                     ← regenerable intermediates
    sub-XXX_<medium>_parameters_<yyMMdd_HHmm><affix>.mat   parameters snapshot at start
    sub-XXX_<medium>_report_parameters<affix>.mat          final parameters (report re-rendering)
    sub-XXX_<medium>_cropped_smoothed<affix>.mat           cropped/smoothed head (layered)
    sub-XXX_<medium>_grid_cache<affix>.mat                 grid stage (io.save_grid_cache)
    sub-XXX_<medium>_planimg<affix>.mat                    planning-image sidecar (always)
    sub-XXX_<medium>_medium_cache<affix>.mat               medium stage (io.save_medium_cache)
    sub-XXX_<medium>_kwave_source<affix>.mat               source (io.save_source_matrices)
    sub-XXX_<medium>_results<affix>.mat                    acoustic results + provenance
    sub-XXX_<medium>_heating_res<affix>.mat                thermal results (io.save_thermal_matrices)
    sub-XXX_<medium>_HEAT<affix>.mat                       thermal time courses (thermal analysis; not BabelBrain)

log/
    sub-XXX_<medium>_<timestamp><affix>.txt                MATLAB diary of the run
    sub-XXX_<medium>_parameters<affix>.yaml                resolved parameters (human-readable)

debug/                                     ← diagnostic artefacts (debug=1 only)
    preproc/        head preprocessing (rotation, cropping, skull visualisations)
    medium/         medium mapping (grid-space property matrices, pCT)
    source/         source/transducer setup (element distribution plots)
    babelbrain/     BabelBrain FOV and domain-coverage diagnostics
```

- **Cache files** are regenerable. With `io.overwrite_files = 'never'`, existing files are reloaded instead of recomputed ([Overwrite behaviour](doc_outputs_reuse.md#overwrite-behaviour)). `io.save_matrices = 0` writes none of the optional cache matrices; `1` saves all of them. The parameter snapshots and the planning-image sidecar are always written. Per-step flags such as `io.save_acoustic_matrices` override it when set.
- **Debug files** are written only with `parameters.simulation.debug = 1`.
- **Run statistics**: each run with an acoustic or thermal simulation adds one line to `<path.sim>/run_stats.jsonl` (one file per `path.sim`, shared by all subjects): grid size, time steps, thermal steps, simulation durations, solver, partition and peak RAM. The Nexus estimates resources from it ([Resource estimates](doc_hpc.md#resource-estimates-from-earlier-runs)).

---

## Where Outputs Live

Each location belongs to one level. A level only reads the levels above it, and each folder has one writer.

| Level | Depends on | Locations | Written by |
|---|---|---|---|
| **Inputs** | none | `path.anat` (T1/T2), `path.localite_raw` (raw [Localite](https://www.localite.de/) exports), `path.localite_t1` (planning T1s) | you / the scanner / Localite; PRESTUS only reads them |
| **Subject derivatives** | the subject's images | `path.seg` ([SimNIBS](https://simnibs.github.io/simnibs/build/html/index.html) `m2m_sub-XXX`), `path.pct` (pseudo-CT), `path.localite_post` (ingested Localite markers) | segmentation (layered media; `io.overwrite_simnibs`), pseudo-CT generation, `neuronav_ingest_markers` |
| **Placement** | subject derivatives, target, transducer geometry, placement settings | `path.placement` (default `<path.sim>/placement`): `sub-XXX/<mode>/<target_name>/` | the placement stage (heuristic, [PlanTUS](https://github.com/mlueckel/PlanTUS), [TUS_entry](https://github.com/CyrilAtkinson/TUS_entry)) |
| **Simulation** | everything above, plus the acoustic/thermal setup | `path.sim/sub-XXX/` (this page) | the simulation pipeline |

- **Placement** does not depend on the medium, grid or drive level. Point several setups at one `path.placement` to reuse a heuristic, PlanTUS or TUS_entry result. Each run folder holds everything the run produced and a run record, `placement_run.json`. `manual`, `mni` and `localite` placements have no run folder. The simulation records the positions it used in `sub-XXX_placement<affix>.json`.
- **`path.localite_export`** is a hand-off folder, not a level. It holds copies of placement results for import into Localite (heuristic position table, Localite-ready T1).
- **SimNIBS `m2m_sub-XXX`** folders are only read after segmentation. For the pseudo-CT location, see [(pseudo-)CT](doc_pseudoCT.md).

---

## NIfTI Images

PRESTUS exports 3D volumes (2D for axisymmetric simulations) in subject-native T1 space and MNI-152 space. Each group of maps is written as soon as its stage finishes.

- Subject T1 space: `sub-XXX_<medium>_T1w<affix>_<type>.nii.gz`
- MNI space: `sub-XXX_<medium>_MNI<affix>_<type>.nii.gz`

??? note "Medium maps (written after medium setup)"
    | Type | Resampling | Description |
    |---|---|---|
    | `medium_masks` | nearest | Tissue label mask (integer indices) |
    | `sound_speed` | nearest | Speed of sound [m/s]; edges filled with water value |
    | `density` | nearest | Mass density [kg/m³]; edges filled with water value |
    | `alpha_coeff` | nearest | Attenuation coefficient [dB/(cm·MHz^y)]; edges filled with water value |

??? note "Acoustic maps (written after acoustic analysis)"
    | Type | Description |
    |---|---|
    | `intensity` | Pulse-average intensity (IPA) [W/cm²] |
    | `MI` | Mechanical index |
    | `pressure` | Peak pressure [Pa]; with a JSON sidecar so a later run can start from it ([Starting from earlier stage outputs](doc_outputs_reuse.md#starting-from-earlier-stage-outputs)) |
    | `focal_mask_3dB` | −3 dB focal zone binary mask (IPA ≥ 50 % of peak; restricted to brain for layered, to beyond-exit-plane for free-field) |
    | `focal_mask_6dB` | −6 dB / FWHM focal zone binary mask (IPA ≥ 25 % of peak = P ≥ Pmax/2; same region restriction) |

??? note "Thermal maps (written after thermal analysis)"
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
| **Thermal results** | temperature, CEM43 (kWave) and CEM43 (ISO) plots (if heating was simulated) |
| **Post-hoc water simulation** | free-field reference results (if run) |
| **Log** | full pipeline execution log |

![PRESTUS_html](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/html_report.png)

Example report on a phantom: [Download example HTML](https://raw.githubusercontent.com/jkosciessa/PRESTUS_bin/main/examples/example_report.html)

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
