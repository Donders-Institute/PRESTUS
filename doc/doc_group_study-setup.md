# Study Setup

A study runs many simulations, one per subject, target, or condition. A study manifest lists them with one row per run and the settings that vary, so one config and one command run them all, locally or on a cluster.

On this page:

- [Study Setup](#study-setup): the manifest, running a study, its columns, output names, and importing runs from Localite.

!!! note "At a glance"
    - **When you need it:** To run several subjects or conditions with one config.
    - **Before you start:** A working single-run study config ([Configuring a Study](doc_configuration.md)).

## Study Setup

A study is a normal PRESTUS config (the *study config*) whose `study.manifest` names a table with one row per run. The study config holds everything the runs share. The manifest holds what differs between runs.

```yaml
# config_pilot.yaml (study config)
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
  name: pilot            # outputs in <path.sim>/group/pilot/
```

```text
# study.tsv: one row per run; tab-separated (.tsv) or comma-separated (.csv)
subject_id  condition  placement.localite.session  placement.localite.target  transducer.target_isppa_wcm2  analysis.target_mni_mm
1           L_M1       ses-02                      L_M1                                                     [-38, -20, 58]
1           R_M1       ses-02                      R_M1                                                     [38, -20, 58]
2           L_M1       ses-02                      L_M1                       20                            [-38, -20, 58]
```

### Running a Study

```matlab
prestus_study_start('config/config_pilot.yaml', 'DryRun', true)   % load + check every row, print the plan
prestus_study_start('config/config_pilot.yaml')                   % run / submit all included rows
prestus_study_start(cfg, 'Subjects', [1 3], 'Conditions', "L_M1")   % a subset
prestus_study_start(cfg, 'Rows', 1, 'Platform', 'matlab', ...
    'Overrides', {'grid.resolution_mm', 2})                         % quick debug run of one row
```

??? note "Dry run, platform, and the study report"
    - **Config errors surface first.** All selected rows are loaded before anything runs. The dry run prints each row's placement mode, target, transducer, planned depth, target intensity, and platform, and writes `<path.sim>/group/<name>/log/plan_<time>.tsv`. A real launch writes `launch_<time>.tsv` with the job ID and status per row.
    - **Platform.** Each run is started with `prestus_pipeline_start` and its own `platform`. On `matlab` the rows run one after another; a failing row is logged and the others continue (unless `'StopOnError', true`). On `slurm`/`qsub` there is one job per row, or a job graph for uncertainty, multi-I<sub>SPPA</sub>, and async rows.
    - **One report at the end.** Study runs switch off the per-subject group report (`modules.generate_group_report = 0`). `prestus_study_start` applies this after the manifest columns, so a manifest column cannot switch it back on; its `'Overrides'` option can. The per-subject simulation report (`modules.generate_report`) is not affected. The study report is built once at the end: directly on `matlab`, and on HPC as a dependent CPU job (`afterany`, so it also runs when some subject jobs fail). Rows that submit a job graph return no single final job ID; for those, run `prestus_group_report_start(cfg)` once their jobs have finished.

### Manifest Columns

`subject_id` is the only required column; any dotted config key sets a per-run value.

??? note "Column reference and cell syntax"
    | Column | Meaning |
    |---|---|
    | `subject_id` | Required. Integer, `003` or `sub-003`. |
    | `condition` | Optional run label (target, protocol, session, …). Becomes the run's output affix and the grouping variable in the group report. |
    | `include` | Optional `1`/`0`, `true`/`false`, `yes`/`no`; default `1`. |
    | `config` | Optional per-row [YAML](https://yaml.org/) merged over the study config. The path is relative to the manifest. |
    | `notes` | Ignored free text (also `note`, `comment`). |
    | any dotted config key | A per-run value, e.g. `placement.localite.target`, `transducer.focal_distance_ep`, `grid.resolution_mm`. Struct arrays take an index: `transducer(2).trans_pos` for the second transducer of a multi-transducer setup. A plain `transducer.x` is only allowed while the config defines one transducer. |

    Cells are parsed like YAML values: `30`, `[-38, 20, 52]`, `L_M1`, `true`. An **empty cell keeps the study value**, and `[]` sets an empty value. Lines starting with `#` are comments. In `.csv` files, quote cells that contain commas: `"[-38, 20, 52]"`. Unknown keys are set anyway, with a warning (a typo guard).

### Output Names

Each run writes the usual per-subject outputs to `<path.sim>/sub-XXX/`, distinguished by its output affix.

??? note "Affix rules and clashes"
    - An explicit `io.output_affix` column (or row-config value) is used as is.
    - Otherwise the condition is appended to the study affix with `-`. For example, study affix `pilot` with condition `M1` gives `_desc-pilot-M1`; without a study affix, condition `M1` gives `_desc-M1`.
    - The result is sanitised like any `io.output_affix`.

    Two included rows that would write to the same files (same subject, affix, medium, and `path.sim`) are rejected before anything runs. Give them distinct conditions.

??? note "Merge order"
    For each row, later sources win:

    `config_default.yaml` < study config < row `config` YAML < manifest columns < `Overrides` of `prestus_study_start` < the run's `subject_id` and output affix

    All of this goes through one `load_parameters` call. Derived values (transducer geometry, thermal timing, …) are therefore computed from the row's own settings. `load_parameters` accepts the same kind of overrides directly: `load_parameters('config.yaml', dir, {'subject_id', 3; 'placement.localite.target', 'L_M1'})`.

### Importing Runs from Localite

Ingest the [Localite](https://www.localite.de/) recordings on the Nexus's *Placement ▸ Localite* page or with `neuronav_gui` (single sessions via *Write outputs*, many via *Ingest all sessions…* / the batch window), then use **Import from Localite…** on the Nexus's Study tab (or the shortcut on the Placement tab).

??? note "Importer files and row layouts"
    The importer reads what ingestion writes under `path.localite_post`: `neuronav_ingest_summary.csv` and the per-session `InstrumentMarker_*.json` files. It builds rows with `placement.mode = localite`, `placement.localite.session`, and `placement.localite.target`, in one of two ways:

    - **One run per target** (one transducer): the condition is the target, with the session appended when a subject has the same target in several sessions.
    - **One run per session** (several transducers): the session's targets are assigned to transducers 1, 2, … as `placement.localite.target: [L_M1, R_M1]`. Order them left first (`L…` before `R…`), alphabetically, or in a given order. Sessions whose target count does not match the configured transducers are skipped and listed. The condition is the session.

PRESTUS reads the ingested markers itself at run time, so no positions are copied into the manifest.
