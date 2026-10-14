# Group analysis

!!! note "At a glance"
    - **What:** Running many subjects or conditions from one table, and summarising them in group maps and a group report.
    - **When you need it:** For studies with more than a handful of runs.
    - **Before you start:** A study config that works for one subject.

PRESTUS supports multi-subject studies end to end:

- **[Study setup](#study-setup)**: a *study manifest* lists every run (a subject, optionally under a named condition) with its own settings: placement, target, depth, intensity, transducer, … `prestus_study_start` runs all rows locally or on HPC.
- **[Targeting accuracy](#targeting-accuracy)**: each run measures how well the simulated focus hits the planned focus, an MNI target and an MNI ROI (`modules.run_targeting_analysis`).
- **[Group MNI maps](#group-mni-maps)**: focal-volume overlap, mean/SD intensity and mean heating per condition in MNI space (`create_group_maps`).
- **[Group HTML report](#group-html-report)**: one self-contained, interactive report across runs: safety dashboard, study design, condition comparison, targeting accuracy, group maps, outliers, and a tidy one-row-per-run CSV (`prestus_group_report_start`).

The report and the maps never run simulations; they read the per-run outputs that already exist on disk.

## Study setup

A study is a normal PRESTUS config (the *study config*) whose `study.manifest` names a table with one row per run. The study config holds everything the runs share. The manifest holds what differs between runs.

```yaml
# config_amythal.yaml (study config)
path:
  sim: /project/sims
  seg: /project/simnibs
  anat: /project/simnibs
  t1_pattern: 'm2m_sub-%1$03d/T1.nii.gz'
simulation:
  medium: layered
transducer:
  serial: CTX_500_026
  target_isppa_wcm2: 30
placement:
  mode: localite
study:
  manifest: study.tsv      # relative to this config's folder
  name: amythal            # outputs in <path.sim>/group/amythal/
```

```text
# study.tsv — one row per run; tab-separated (.tsv) or comma-separated (.csv)
subject_id  condition  placement.localite.session  placement.localite.target  transducer.target_isppa_wcm2  analysis.target_mni_mm
1           L_MD       ses-02                      L_MD                                                     [-8, -16, 6]
1           R_MD       ses-02                      R_MD                                                     [8, -16, 6]
2           L_MD       ses-02                      L_MD                       20                            [-8, -16, 6]
```

### Manifest columns

| Column | Meaning |
|---|---|
| `subject_id` | Required. Integer, `003` or `sub-003`. |
| `condition` | Optional run label (target, protocol, session, …). Becomes the run's output affix and the grouping variable in the group report. |
| `include` | Optional `1`/`0`, `true`/`false`, `yes`/`no`; default `1`. |
| `config` | Optional per-row YAML merged over the study config. The path is relative to the manifest. |
| `notes` | Ignored free text (also `note`, `comment`). |
| any dotted config key | A per-run value, e.g. `placement.localite.target`, `transducer.focal_distance_ep`, `grid.resolution_mm`. Struct arrays take an index: `transducer(2).trans_pos` for the second transducer of a multi-transducer setup. A plain `transducer.x` is only allowed while the config defines one transducer. |

Cells are parsed like YAML values: `30`, `[-38, 20, 52]`, `L_MD`, `true`. An **empty cell keeps the study value**, and `[]` sets an empty value. Lines starting with `#` are comments. In `.csv` files, quote cells that contain commas: `"[-38, 20, 52]"`. Unknown keys are set anyway, with a warning (a typo guard).

### Merge order

For each row, later sources win:

`config_default.yaml` < study config < row `config` YAML < manifest columns < `Overrides` of `prestus_study_start` < the run's `subject_id` and output affix

All of this goes through one `load_parameters` call. Derived values (transducer geometry, thermal timing, …) are therefore computed from the row's own settings. `load_parameters` accepts the same kind of overrides directly: `load_parameters('config.yaml', dir, {'subject_id', 3; 'placement.localite.target', 'L_MD'})`.

### Output names

Each run writes the usual per-subject outputs to `<path.sim>/sub-XXX/`, distinguished by its output affix:

- An explicit `io.output_affix` column (or row-config value) is used as is.
- Otherwise the condition is appended to the study affix with `-`. For example, study affix `liberal` with condition `M1` gives `_desc-liberal-M1`; without a study affix, condition `M1` gives `_desc-M1`.
- The result is sanitised like any `io.output_affix`.

Two included rows that would write to the same files (same subject, affix, medium and `path.sim`) are rejected before anything runs. Give them distinct conditions.

### Running a study

```matlab
prestus_study_start('config/config_amythal.yaml', 'DryRun', true)   % load + check every row, print the plan
prestus_study_start('config/config_amythal.yaml')                   % run / submit all included rows
prestus_study_start(cfg, 'Subjects', [1 3], 'Conditions', "L_MD")   % a subset
prestus_study_start(cfg, 'Rows', 1, 'Platform', 'matlab', ...
    'Overrides', {'grid.resolution_mm', 2})                         % quick debug run of one row
```

- **Config errors surface first.** All selected rows are loaded before anything runs. The dry run prints each row's placement mode, target, transducer, planned depth, target intensity and platform, and writes `<path.sim>/group/<name>/log/plan_<time>.tsv`. A real launch writes `launch_<time>.tsv`, with the job ID and status per row.
- **Platform.** Each run is started with `prestus_pipeline_start` and its own `platform`. On `matlab` the rows run one after another in the current session; a failing row is logged and the others continue (unless `'StopOnError', true`). On `slurm`/`qsub` there is one job per row, or a job graph for uncertainty, multi-ISPPA and async rows.
- **One report at the end.** The per-subject automatic group report is switched off for study runs (`modules.generate_group_report = 0`, unless the manifest sets it). The study report is built once at the end instead: directly on `matlab`, and on HPC as a dependent CPU job (`afterany`, so it also runs when some subject jobs fail).
- **Job graphs.** Rows that submit a job graph (uncertainty, multi-ISPPA, async multi-transducer) return no single final job ID. For those, run `prestus_group_report_start(cfg)` once their jobs have finished.

### From Localite ingestion

Ingest the recordings on the Nexus's *Placement ▸ Localite* page or with `neuronav_gui` (single sessions via *Write outputs*, many via *Ingest all sessions…* / the batch window), then use **Import from Localite…** on the Nexus's Study tab (or the shortcut on the Placement tab). The importer reads what ingestion writes under `path.localite_post`: `neuronav_ingest_summary.csv` and the per-session `InstrumentMarker_*.json` files. It then builds rows with `placement.mode = localite`, `placement.localite.session` and `placement.localite.target`, in one of two ways:

- **One run per target** (one transducer): the condition is the target, with the session appended when a subject has the same target in several sessions.
- **One run per session** (several transducers): the session's targets are assigned to transducers 1, 2, … as `placement.localite.target: [L_MD, R_MD]`. Order them left first (`L…` before `R…`), alphabetically, or in a given order. Sessions whose target count does not match the configured transducers are skipped and listed. The condition is the session.

PRESTUS reads the ingested markers itself at run time, so no positions are copied into the manifest. This replaces study-specific scripts that inject `trans_pos`/`focus_pos` per session.

## Targeting accuracy

With `modules.run_targeting_analysis: 1` (default, layered or phantom media, `run_nifti_creation` on), each run writes `sub-XXX_<medium><affix>_targeting.csv` next to its per-subject CSV. The analysis reads the focal-mask NIfTIs and never fails the pipeline; unavailable measures are `NaN`. The same function can be run post hoc: `targeting_analysis(parameters, 'PlannedFocus', f, 'PlannedTrans', t)` with T1 voxel positions.

| Columns | Meaning |
|---|---|
| `condition`, `placement_mode`, `target_label`, `transducer_name`, `n_transducers`, `focal_distance_ep_planned_mm`, `target_isppa_wcm2` | The run's settings, shown in the report's study design table. |
| `focal_centroid_ras_*`, `peak_ras_*` | Centroid of the −3 dB focal volume (brain-masked) and the intensity peak inside it, in T1 RAS mm. |
| `planned_focus_ras_*`, `planned_trans_ras_*` | Focus and transducer positions resolved by placement, recorded before the grid stage, in RAS mm. |
| `focus_error_mm` | Distance between the focal centroid and the planned focus. |
| `focus_error_axial_mm`, `focus_error_lateral_mm` | That error along and across the planned beam axis (transducer → planned focus); axial `> 0` means the simulated focus is deeper than planned. |
| `peak_error_mm` | Distance between the peak and the planned focus. |
| `focal_centroid_mni_*`, `peak_mni_*` | The same positions in MNI mm (needs `io.save_MNI`). |
| `target_mni_*`, `target_mni_source` | The intended MNI target: `analysis.target_mni_mm`, else the MNI target of the placement mode (`placement.heuristic` / `placement.plantus` `.mni_target_mm`, `placement.mni.focus_pos_mm`). |
| `target_error_mni_mm`, `peak_target_error_mni_mm` | Distance from the MNI focal centroid (or peak) to the MNI target. |
| `roi_volume_mm3`, `roi_overlap_mm3`, `roi_focal_fraction`, `roi_coverage`, `roi_mean_intensity`, `roi_max_intensity`, `roi_centroid_dist_mm` | With `analysis.target_roi_mni` (an MNI NIfTI mask, any grid): ROI volume, overlap with the focal volume, fraction of the focal volume inside the ROI, fraction of the ROI covered, intensity in the ROI (`intensity_quantity` = Ipa, else Itp), and the distance between the centroids. |

With several transducers, the planned focus and beam axis of `transducer(1)` are used. Targets and ROIs can differ per run: set `analysis.target_mni_mm` or `analysis.target_roi_mni` as manifest columns.

## Group MNI maps

`prestus_group_report_start` builds these by default (`'Maps', false` to skip). They come from `create_group_maps` and use each run's MNI NIfTIs. The maps are computed per condition, or across all runs when there are no conditions:

| File in `<out>/maps/` | Content |
|---|---|
| `group_<cond>_focal_overlap_3dB.nii.gz` | Fraction of runs whose −3 dB focal volume covers each voxel (0–1). |
| `group_<cond>_<Ipa\|Itp>_mean.nii.gz`, `…_sd.nii.gz` | Voxelwise mean and SD intensity (Ipa when every run has it, else Itp when every run has that; no intensity maps otherwise). |
| `group_<cond>_heatrise_mean.nii.gz` | Mean temperature rise (when every run has heating maps). |

`<out>/plots/group_<cond>_focal_overlap.png` shows three orthogonal MNI slices through the overlap peak. The background is the mean MNI-space T1 of the subjects (`<path.seg>/m2m_sub-XXX/toMNI/T1_to_MNI_post-hoc.nii.gz`, when present). The figure also shows the ROI outline, the per-run focal centroids and the targets. Runs whose MNI grid differs from the first run's are skipped with a warning. `<out>` is `<path.sim>/group/<study.name>/` for a study and `<path.sim>/group/<medium><affix>/` otherwise.

## Group HTML report

```matlab
prestus_group_report_start('config/config_amythal.yaml');           % study: runs from the manifest
prestus_group_report_start('config/config_single.yaml');            % no manifest: auto-discover subjects
prestus_group_report_start('config/config_single.yaml', [1 2 3 5]); % explicit subject list
prestus_group_report_start(cfg, [], 'Maps', false);                 % report only
```

| Mode | Runs | Written to |
|---|---|---|
| Study (`study.manifest` set) | Included manifest rows | `<path.sim>/group/<name>/group_<name>_report.html`, `group_<name>_runs.csv`, `plots/`, `maps/` |
| Subjects | All subjects in `path.sim` with `sub-XXX_<medium><affix>.csv` (or the given list) | `<path.sim>/group_<medium>_report<affix>.html`, `group_<medium>_runs<affix>.csv`, `group_plots/` |

**Automatically**, with `modules.generate_group_report: 1`, `prestus_pipeline` rebuilds the subject-mode report after each subject (best-effort, wrapped in `try`/`catch`, without group maps). Runs started by `prestus_study_start` skip this and build one study report at the end.

The report reads each run's `sub-XXX_<medium><affix>.csv` (last row) and, when present, the targeting sidecar and the per-subject PNGs. Runs without a CSV are listed as missing. Report contents:

- **Header**: study, manifest, runs, subjects, conditions, missing runs, medium, and a link to the one-row-per-run CSV.
- **Safety banner**: runs above any ITRUSST non-significant-risk reference level, with links to their rows.
- **Filter**: checkboxes per run and per-condition buttons. The dashboard, roster and cards update live; box plots and figures show the full set.
- **Group safety dashboard and other metrics**: mean ± SD, range and N per ITRUSST metric, coloured by the worst run against the NSR limit. The Mechanical sub-grid shows both MI (free water) and MItc tiles plus Max pressure; the tile that does not apply to the medium is greyed out.
- **Study design** (study, or more than one condition): runs per condition, the manifest keys and how many values each takes, and each run's placement mode, target, transducer, planned depth, target intensity and manifest settings.
- **Condition comparison** (more than one condition): mean ± SD (n) of the key metrics per condition, plus box plots. It is descriptive only: runs of one subject under several conditions are not independent.
- **Targeting accuracy** (targeting sidecars present): focus and MNI-target errors, ROI measures per condition and per run, and a figure of the MNI centroids and the axial/lateral error.
- **Group MNI maps** (when built): overlap figures and links to the NIfTIs.
- **Outliers**: runs whose value lies far from their condition's median (modified z-score 0.6745·(x − median)/MAD, |z| > 3.5; Iglewicz & Hoaglin, 1993), for conditions with at least 5 runs. They are also flagged in the roster.
- **Data analysis**: client-side distribution, correlation, and by-group views (group by subject, condition, placement mode, target or transducer when they vary).
- **Roster**: one sortable row per run with key metrics, NSR and outlier flags, and a link to the run's own report.
- **Group acoustic and thermal summaries**: box plots (Isppa per tissue, MI per tissue; maximum temperature, temperature rise and CEM43 when thermal results exist).
- **Per-run cards**: the run's settings and figures.
- **Configuration summary**: of the loaded (study) config. Per-run settings are in the study design table.

## Deprecated: `prestus_group_start` / `create_group_MNI_plots`

`create_group_MNI_plots` expects an older output layout (`parameters.temp_output_dir`, `*_final_intensity_MNI*.nii.gz`, `*_output_table*.csv`) and no longer finds the files the pipeline writes. `prestus_group_start` now warns and runs `prestus_group_report_start` (group MNI maps + report) instead. Its plot options are ignored; for an ROI, set `analysis.target_roi_mni`. The ROI statistics of the old function (intensity in the ROI, overlap, centroid distance) are now part of the [targeting analysis](#targeting-accuracy), in mm rather than voxels.
