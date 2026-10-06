# PRESTUS Nexus

A standalone desktop app for putting a PRESTUS study config together before running it in MATLAB. It writes a normal override YAML for `load_parameters`, containing only the values that differ from `config_default.yaml`. MATLAB still does all loading, validation and simulation.

Status: *experimental*. All five planned phases are in place: config builder, transducer library panel, placement viewer, MATLAB validate/run, and deprecation of the old tools. User documentation is in [doc/doc_nexus.md](../doc/doc_nexus.md).

## Running

Requires Python ≥ 3.10, PyQt6 and PyYAML. The placement viewer also needs numpy, nibabel and matplotlib. `conda env create -f environment.yml` creates a `prestus-nexus` env with all of them.

```bash
cd PRESTUS/nexus
python -m prestus_nexus                    # untitled study; File > New study… / Open…
python -m prestus_nexus path/to/config_x.yaml
# or: pip install -e ".[viewer]"   then   prestus-nexus
```

From MATLAB: `prestus_nexus` or `prestus_nexus('config_x.yaml', Python='/path/to/python')`.

PRESTUS itself is found automatically, as the checkout this folder lives in. You can override that with `--prestus-root` or `$PRESTUS_ROOT`.

## Workflow

1. **Study config.** Opens the given config, or an untitled one. *File ▸ New study…* creates a config and can copy `config_default.yaml` into the study folder, as `prestus_config_init` does (this fixes the defaults the study is merged onto). *Open…* / *Open recent* switch configs. The subject inputs are on the Setup tab; the Issues panel checks them as you type (the T1 pattern must resolve with MATLAB's `sprintf` and glob semantics, and `m2m_sub-XXX` must exist or SimNIBS must be available to run charm).
2. **Tabs and sections.**
   - Six tabs (Setup, Transducer, Placement, Simulation, Study, Run) hold one card per config section; `overlay.yaml` (`tabs`, `display`) assigns sections to tabs and gives fields readable labels, advanced flags and a home section. Advanced settings stay hidden until a card's *+ n advanced* button or the tab's *Show advanced options* is used.
   - *Tune for* (shared by all tabs) shows and marks the keys tagged with a goal (speed, accuracy, memory, disk space) in `overlay.yaml` (`tune_goals`, `tuning`), advanced ones included.
   - Changed values are marked ●, each field has a reset button, and hovering a field (or **?**) shows its help text, taken from the comments in `config_default.yaml`, with its key and default.
   - Mode-specific fields appear only when they apply, e.g. `placement.localite.*` only for `placement.mode: localite`.
   - Advanced sections are hidden unless they contain changes.
3. **Transducers.** The transducer section can hold several transducers. As in PRESTUS, transducer 2 and later inherit every value set on transducer 1 unless they override it.
4. **Issues panel.**
   - Lists problems PRESTUS would otherwise report only at load or run time: wrong types, unknown values, unknown or misspelt keys (with suggestions), deprecated keys, backend combinations, frequency mismatches between transducers, and missing inputs.
   - Clicking an issue jumps to the field.
5. **Review & save.** Shows the exact YAML that will be written and which `config_default.yaml` it is merged onto.

## How it is built

| File | Role |
|---|---|
| `schema.py` | Builds the key tree from `config/config_default.yaml`: key order, defaults, and help text taken from inline comments and the comment blocks above keys. A project copy of `config_default.yaml` is merged on top, as `load_parameters` does (`stale_keys`, `pinned`). |
| `overlay.yaml` | Adds only what the YAML can't express: widget kinds, allowed values, units, visibility rules, keys the code reads that the default omits (`extra`), deprecated aliases, and tuning goals per key (`tuning`). |
| `model.py` | Defaults plus sparse overrides, merged like `mergestruct` + `restore_struct_arrays` in `load_parameters`. |
| `yamlio.py` | YAML reading and writing that matches MATLAB's loader (see below). |
| `checks.py` | Validation. Each rule notes the MATLAB source it mirrors. |
| `equipment.py` | Equipment library (`config/equipment`) and calibration coverage (`config/calibration`). |
| `placement.py` | Voxel conventions, bowl geometry, scalp snapping, Brainsight/Slicer trajectory import. |
| `launcher.py` | MATLAB discovery, `run_<study>.m` / `submit_<study>.m` launchers, `-batch` commands with success/failure markers. |
| `jobs.py` | Cluster jobs submitted through PRESTUS: IDs from the log, `squeue`/`sacct`/`qstat` state, cancel. |
| `merge.py` | `mergestruct` semantics, shared by the schema (toolbox + project defaults) and the model. |
| `project_defaults.py` | Refreshing a project's `config_default.yaml` copy (*Tools ▸ Update project defaults…*). |
| `outputs.py` | Where the pipeline writes per-subject outputs; the placement record (`sub-XXX_placement<affix>.json`). |
| `ui/` | PyQt6 main window (tabs), section forms, transducer panel and cross-section, placement viewer, timing diagram, study manifest and Localite import, run page, style. |
| `timing.py` | Protocol maths and the TUS calculator settings exchange. |
| `study.py`, `neuronav_import.py` | Study manifest model; manifest rows from the `neuronav_gui` outputs. |
| `labels.py` | Readable field labels. |

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
   - **Tools menu:** opens `neuronav_gui` (pre-filled from the saved study config) and `calibration_gui` (saving to `transducer.library_path`) in an interactive MATLAB session.
5. **Retirement.** `launch_web_alignment` and `functions/gui_py` have been removed. The MATLAB entry point is `functions/gui/prestus_nexus.m`. Packaging: `pyproject.toml` with a `viewer` extra, and `environment.yml`.

## Voxel convention note

PRESTUS voxels are 1-based MATLAB indices into the T1 array, while MATLAB's `niftiinfo` transform maps *0-based* indices. `ras_to_grid` / `grid_to_ras` in PRESTUS and `placement.py` in the Nexus both apply this rule, so `coordinate_space: voxel` and `ras` place the transducer at the same anatomical position.

## Not done yet

- A 3-D view (NiiVue via QtWebEngine) is not included; the viewer uses 2-D orthogonal slices, as BabelBrain does.
- Non-manual placement modes are edited as forms; the viewer shows their resolved positions read-only (from the last run, see `outputs.placement_record`) but cannot place them interactively.
- The standalone app (`packaging/prestus_nexus.spec`) is not built or signed in CI yet.
