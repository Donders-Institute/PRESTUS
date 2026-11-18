# PRESTUS Nexus

A standalone desktop app for putting a PRESTUS study config together before running it in MATLAB. It writes a normal override YAML for `load_parameters`, containing only the values that differ from `config_default.yaml`. MATLAB still does all loading, validation and simulation.

In place: config builder, transducer library panel, placement viewer with 3D view, MATLAB validate/run, results browser with follow-up runs, deprecation of the old tools, and the Localite and Calibration pages (MATLAB GUIs integrated through a persistent MATLAB session). User documentation is in [doc/doc_nexus.md](../doc/doc_nexus.md).

## Running

Install (Python ≥ 3.10; all packages are listed in `pyproject.toml`):

```bash
bash nexus/install.sh                      # macOS/Linux; nexus\install.bat on Windows
# or, with conda set up, in this folder:
conda env create -f environment.yml        # creates the env 'prestus-nexus' with everything
# without conda, in any Python ≥ 3.10 environment:  pip install -e .
```

The installers create or update the env `prestus-nexus` and install Miniforge first if conda is missing.

Run:

```bash
prestus-nexus                              # untitled study; File > New study… / Open…
prestus-nexus path/to/config_x.yaml
python -m prestus_nexus [config_x.yaml]    # the same, from this folder
```

From MATLAB: `prestus_nexus` or `prestus_nexus('config_x.yaml')`. The Python is found by `find_nexus_python` (`$PRESTUS_NEXUS_PYTHON`, the last one used, the env `prestus-nexus`, `python3`) and remembered; `Python='/path/to/python'` overrides it.

PRESTUS itself is found automatically, as the checkout this folder lives in. You can override that with `--prestus-root` or `$PRESTUS_ROOT`.

## Workflow

1. **Study config.** Opens the given config, or an untitled one. *File ▸ New study…* creates a config and can copy `config_default.yaml` into the study folder, as `prestus_config_init` does (this fixes the defaults the study is merged onto). *Open…* / *Open recent* switch configs. The subject inputs are on the Setup tab; the Issues panel checks them as you type (the T1 pattern must resolve with MATLAB's `sprintf` and glob semantics, and `m2m_sub-XXX` must exist or SimNIBS must be available to run charm).
2. **Tabs and sections.**
   - *QuickStart* sets up a first run in four steps (T1 image, transducer, MNI target, run); everything else keeps its default. It derives `path.anat`, `t1_pattern`, `subject_id`, the T2, `path.seg` and a results folder from the chosen T1 (`quickstart.py`).
   - Nine tabs (Setup, Transducer, Placement, Medium, Sonication, Solver, Study, Run, Results) hold one card per config section; `overlay.yaml` (`tabs`, `display`) assigns sections to tabs and gives fields readable labels, advanced flags and a home section. Advanced settings stay hidden until a card's advanced button or the tab's *Advanced* toggle is used.
   - The goal tags in `overlay.yaml` (`tune_goals`, `tuning`) are kept for a *Tune for* control that is switched off (`TUNE_FOR` in `ui/main_window.py`).
   - Each card's ↗ button, the field help and *Help ▸ Docs for this tab* (F1) open the page named by `doc:` / `docs:` in `overlay.yaml` (`docs.py`).
   - Changed values are marked ●, each field has a reset button, and hovering a field (or **?**) shows its help text, taken from the comments in `config_default.yaml`, with its key and default.
   - Mode-specific fields appear only when they apply, e.g. `placement.localite.*` only for `placement.mode: localite`.
   - Advanced sections are hidden unless they contain changes.
3. **Transducers.** The transducer section can hold several transducers. As in PRESTUS, transducer 2 and later inherit every value set on transducer 1 unless they override it.
4. **Issues panel.**
   - Lists problems PRESTUS would otherwise report only at load or run time: wrong types, unknown values, unknown or misspelt keys (with suggestions), deprecated keys, backend combinations, frequency mismatches between transducers, and missing inputs.
   - Closed by default: the header shows the error and warning counts; clicking them (or *View ▸ Issues…*, Ctrl+Shift+I) opens the list as a popup.
   - Clicking an issue jumps to the field.
5. **Review & save.** Shows the exact YAML that will be written and which `config_default.yaml` it is merged onto. Saves are atomic (`fsio.py`); *Edit ▸ Undo / Redo* step through changes (`undo.py`).
6. **Results and follow-ups.** The Results tab browses the runs in the output folder; *Follow-up…* writes a new study config that starts from the selected run (`followup.py`), and the Run tab's *Pipeline options* pass sequential follow-ups and per-stage limits to `prestus_pipeline_start` (`launcher.PipelineOptions`).

## How it is built

| File | Role |
|---|---|
| `schema.py` | Builds the key tree from `config/config_default.yaml`: key order, defaults, and help text taken from inline comments and the comment blocks above keys. A project copy of `config_default.yaml` is merged on top, as `load_parameters` does (`stale_keys`, `pinned`). |
| `overlay.yaml` | Adds only what the YAML can't express: widget kinds, allowed values, units, visibility rules, keys the code reads that the default omits (`extra`), deprecated aliases, tuning goals per key (`tuning`), readable drop-down labels for choice values (`choice_labels`; the config keeps the value), and schematic icons and figures for drop-down choices (`choice_icons`, `choice_figures`; drawn in `ui/choice_icons.py`). |
| `model.py` | Defaults plus sparse overrides, merged like `mergestruct` + `restore_struct_arrays` in `load_parameters`. |
| `yamlio.py` | YAML reading and writing that matches MATLAB's loader (see below). |
| `checks.py` | Validation. Each rule notes the MATLAB source it mirrors. |
| `equipment.py` | Equipment library (`config/equipment`) and calibration coverage (`config/calibration`). |
| `placement.py` | Voxel conventions, bowl geometry, scalp snapping, Brainsight/Slicer trajectory import. |
| `launcher.py` | MATLAB discovery, `run_<study>.m` / `submit_<study>.m` launchers, `-batch` commands with success/failure markers. |
| `jobs.py` | Cluster jobs submitted through PRESTUS: IDs from the log, `squeue`/`sacct`/`qstat` state, cancel. |
| `ui/jobs_view.py` | The *Jobs* side tab and its pop-out window: job table, and the selected job's output/error/report logs, followed live. |
| `remote.py` | Submitting over SSH: local → cluster path mapping, the `ssh … bash -lc` submit and probe commands, scheduler commands run on the host, and the shared connection (*Connect…*: one password login, reused via `ControlPath`). |
| `askpass.py` | ssh's `SSH_ASKPASS` prompt as a Nexus dialog (`prestus-nexus --askpass`). |
| `merge.py` | `mergestruct` semantics, shared by the schema (toolbox + project defaults) and the model. |
| `project_defaults.py` | Refreshing a project's `config_default.yaml` copy (*Tools ▸ Update project defaults…*). |
| `outputs.py` | Where the pipeline writes per-subject outputs; the placement record (`sub-XXX_placement<affix>.json`). |
| `ui/` | PyQt6 main window (tabs), section forms, transducer panel and cross-section, placement viewer, timing diagram, study manifest and Localite import, run page, QuickStart tab (`start_page.py`), Calibration page (`calibration_page.py`) and its library plots (`calibration_plots.py`), empty-state panels (`empty_state.py`), Localite page with batch and MNI dialogs (`localite_page.py`), MATLAB session status (`matlab_status.py`), style. |
| `timing.py` | Protocol maths and the TUS calculator settings exchange. |
| `study.py`, `neuronav_import.py` | Study manifest model; manifest rows from the Localite ingestion outputs. |
| `bridge.py` | The MATLAB session for interactive steps: file protocol (`BridgeClient`), the Qt client (`MatlabBridge`, one per app), JSON helpers. MATLAB side: `functions/gui/nexus_bridge.m`. |
| `calibration.py`, `localite.py` | Calibration and Localite page logic without Qt: depth lists, the study's focal depths, `meta.fit` ↔ fit settings; raw-root resolution, settings snapshot, target map. |
| `labels.py` | Readable field labels. |
| `quickstart.py` | QuickStart tab logic without Qt: config keys from one T1 file (pattern, subject, T2, segmentation folder), default results folder, SimNIBS discovery on macOS/Windows/Linux. |
| `fsio.py` | Atomic file writes: a temporary file in the same folder replaces the target only once complete. |
| `undo.py` | *Edit ▸ Undo / Redo*: the history of model changes; typing in one field forms one step. |
| `telemetry.py` | *Tools ▸ Usage statistics…*: reads and writes the consent record `~/.prestus/telemetry.json` that MATLAB uses; sends nothing. |
| `docs.py` | Links into the documentation (published pages or local `doc/*.md`) and the Help menu. |
| `multitransducer.py` | Grid orientation of a multi-transducer run under `grid.mode` and each beam's tilt (ports of `mean_beam_axis.m` and `preproc_head.m`). |
| `followup.py` | Follow-up runs without Qt: the five types, their config keys, validation and the `# nexus-followup:` lineage line. |
| `proctree.py` | Stopping a started program with everything it started (MATLAB's launcher and binary). |
| `ui/export.py` | The download button on every figure: SVG, PDF, PNG and copy. |
| `ui/errors.py` | Unexpected errors as a dialog instead of a crash, logged to `nexus_errors.log` in the cache folder. |
| `ui/transducers_page.py` | *Placement ▸ Transducers* (2+ transducers): per-transducer table, grid orientation, coupling and timing notes. |
| `ui/followup_dialog.py`, `ui/pose_table.py` | The *Follow-up…* dialog; the *Poses* sub-tab of a time-reversal pose search with *Use this pose…*. |
| `ui/pipeline_options.py` | *Run ▸ Where to run ▸ Pipeline options*: sequential chain, per-stage limits, summary report switch (stored in the session, not the config). |
| `ui/isppa_sweep.py` | The target ISPPA editor with *Sweep…* for a list of targets (multi-ISPPA). |

### YAML compatibility with MATLAB

PRESTUS reads YAML with `yaml.loadFile` (SnakeYAML 1.30, YAML 1.1). The Nexus therefore:

- **Reads floats with SnakeYAML's own rule.** PyYAML reads `2.17e-3` (used in `config_default.yaml`) as a string; SnakeYAML reads it as a float.
- **Quotes every string.** Unquoted `yes`/`no`/`on`/`off` would load as booleans, and `04:00:00` as a number.
- **Writes floats in a form both parsers read the same way** (`1.0e-05`).
- **Writes every element of a transducer list with the same top-level keys.** `yaml.loadFile` turns the list into a struct array, and that fails if the elements' fields differ. Values that are inherited anyway are written out explicitly, and sub-groups are written as `{}`.

## Tests

```bash
python -m pytest -q                 # everything; the MATLAB test runs if MATLAB is found
python -m pytest -q -m "not matlab" # without MATLAB (~30 s)
```

- `test_schema.py` is the drift test. It checks that every overlay rule (including `tuning`) still matches a key and every tuning goal is declared, that every `extra`/`deprecated` key is still referenced in `functions/`, that no `extra` key has since been added to the default, and that allowed values contain the default. MATLAB's `load_parameters` reads the same `extra`/`deprecated` lists (`config_known_extra_keys`) for its unknown-key warning, so this test also keeps that check accurate.
- `test_matlab_roundtrip.py` writes a study config with the Nexus, loads it with the real `load_parameters` and compares every value with the Nexus's own merge. Set `PRESTUS_MATLAB` if `matlab` is not on your `PATH`.
- `test_real_configs.py` loads and round-trips every study config in `$PRESTUS_REAL_CONFIG_DIR` (skipped when unset).
- `test_bridge.py` tests the bridge protocol and the Qt client against `fake_bridge.py` (no MATLAB), and the real `nexus_bridge` with MATLAB. `test_bridge_pages.py` tests the Calibration and Localite pages on recorded MATLAB answers (`tests/fixtures/bridge`; rerun `make_bridge_fixtures.m` in MATLAB after changing `nexus_calibration` / `nexus_neuronav`). `test_bridge_pages_matlab.py` drives both pages end to end with MATLAB.

When `config_default.yaml` gains or loses keys, run the drift test and update `overlay.yaml` if it fails.

## Features by phase

1. **Foundation.** Schema and drift test, override-only YAML, session memory, schema-driven forms, issues panel, review.
2. **Transducer library** (`equipment.py`, `ui/transducer_panel.py`).
   - Device and driving-system pickers, backed by `config/equipment`.
   - Calibration coverage, looked up by the same file names as `resolve_transducer_from_serial`, including `transducer.library_path`.
   - Checks against `min_foc`/`max_foc`.
   - Library values are shown as the baseline, and a value set in the study wins, as in PRESTUS.
3. **Placement** (`placement.py`, `ui/placement_viewer.py`).
   - Three orthogonal slices of the T1 with the segmentation overlay.
   - Click-to-place, with the bowl drawn from the device geometry.
   - Tools: snap to scalp with a gap, place at the focal distance, and match the focal distance to the placement.
   - Brainsight/Slicer trajectory import using BabelBrain's conventions.
   - Tissue at the focus, and the exit plane's distance from the scalp.
4. **Run** (`launcher.py`, `jobs.py`, `ui/run_page.py`).
   - A `run_<study>.m` launcher.
   - *Validate in MATLAB* (`load_parameters` only).
   - *Run locally* via `matlab -batch`, with a live log and success/failure markers.
   - *Submit to cluster* through PRESTUS' own HPC path, and a job list that follows the submitted jobs.
   - *Submit over SSH* from a computer without `sbatch`/`qsub` that sees the study folder on a shared drive.
   - **Tools menu:** the Localite, PlanTUS, TUS entry and Calibration sub-tabs, and `neuronav_gui` (pre-filled from the saved study config) and `calibration_gui` (saving to `transducer.library_path`) in an interactive MATLAB session.
5. **MATLAB GUIs in the Nexus** (`bridge.py`, `ui/calibration_page.py`, `ui/localite_page.py`).
   - One persistent MATLAB session (`nexus_bridge`, `matlab -batch`) serves interactive steps through JSON request/response files; it ends after 30 min idle or with the app.
   - *Transducer ▸ Calibration*: `calibration_gui`'s predict / evaluate / fit / save for the config's transducers, saving to `transducer.library_path`, predicting at the study's focal depths; k-Wave drive check as a Run-tab script. Computed by `nexus_calibration` → `calibration_fit` etc., the same functions `calibration_gui` uses.
   - *Placement ▸ Localite*: `neuronav_gui`'s detect / name / write with the study config's settings, *Use in this study*, batch ingestion (label map → `neuronav.target_labels`, `ingest_<study>.m` on the Run tab), MNI view. Computed by `nexus_neuronav`.
6. **Retirement.** `launch_web_alignment` and `functions/gui_py` have been removed. The MATLAB entry point is `functions/gui/prestus_nexus.m`. Packaging: `pyproject.toml` (all dependencies) and `environment.yml` (installs it with pip).

## Voxel convention note

PRESTUS voxels are 1-based MATLAB indices into the T1 array, while MATLAB's `niftiinfo` transform maps *0-based* indices. `ras_to_grid` / `grid_to_ras` in PRESTUS and `placement.py` in the Nexus both apply this rule, so `coordinate_space: voxel` and `ras` place the transducer at the same anatomical position.

## Not done yet

- Heuristic, MNI and Localite placements are edited as forms; the viewer shows their resolved positions read-only (from the last run, see `outputs.placement_record`) but cannot place them interactively.
- The standalone app (`packaging/prestus_nexus.spec`) is not built or signed in CI yet.
