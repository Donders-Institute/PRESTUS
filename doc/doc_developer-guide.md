# Developer Guide

Many changes to PRESTUS touch more than one file: a new config key also needs a row in the parameter docs, and a change to a cached stage needs an updated cache signature. This page shows how the code is organised and what happens during a run, so that you know where a change belongs, what else to update and which checks to run before you commit.

On this page:

- [Pipeline and Design Principles](#pipeline-and-design-principles): the pipeline in one paragraph and the three ideas that hold the code together.
- [Code Layout](#code-layout): which folder holds which subsystem, and the file to start reading at.
- [Steps of a Run](#steps-of-a-run): the steps of a run in order, from loading the config to writing outputs.
- [Required Updates per Change](#required-updates-per-change): what to edit, what else to update and which check catches an omission.
- [Pre-Commit Checks](#pre-commit-checks): the checks CI runs, and which test levels a change needs.
- [Nexus Code Structure](#nexus-code-structure): how the Nexus builds its forms, and which Python file does what.

!!! note "At a glance"
    - **When you need it:** Before your first change to the pipeline, the config or the Nexus.
    - **Before you start:** A working installation ([Installation](doc_installation.md)) and one finished simulation ([Your First Simulation](doc_getting-started.md)).

## Pipeline and Design Principles

PRESTUS turns one configuration into one simulation. The configuration is a set of [YAML](https://yaml.org/) files merged into a single [MATLAB](https://www.mathworks.com/products/matlab.html) struct, `parameters`. The pipeline then runs a fixed sequence of stages: segment the head, place the transducer, build the grid and the tissue properties, simulate the sound field, simulate heating, and write reports. Each stage receives what it needs as function arguments and returns its results, so a stage can be read, tested and skipped on its own.

Three ideas hold the code together:

- **The config is the interface.** Every behaviour a user can choose is a key in `config/config_default.yaml`. Code reads its settings from `parameters`, not from globals.
- **Expensive stages are cached.** Grid, medium, acoustic and thermal results are saved under the subject's `cache/` folder. A re-run reuses them; grid, medium and acoustic caches are first checked against the settings they were made with.
- **Results are reported at a free-water intensity.** Simulations run at a nominal drive, and the field is scaled afterwards so that the same transducer in water would give the requested ISPPA.

The Nexus is a separate [Python](https://www.python.org/) application. It edits the YAML and starts MATLAB; it never simulates.

## Code Layout

| Folder | Holds | Start reading at |
|---|---|---|
| `functions/core/` | Config loading, output folders, logging | `load_parameters.m` |
| `functions/prestus_pipeline.m` | The stage sequence | this file |
| `functions/head/`, `medium/` | Segmentation, head preprocessing, tissue properties | `medium_setup.m` |
| `functions/transducer/`, `source/` | Transducer geometry, placement, [k-Wave](https://github.com/ucl-bug/k-wave) source | `source_sensor_setup.m` |
| `functions/acoustic/`, `thermal/` | Solver calls and their analysis | `acoustic_stage.m`, `thermal_stage.m` |
| `functions/transform/` | Coordinate conversions | `ras_to_grid.m` |
| `functions/calibration/` | Fitting transducers to hydrophone data | `refine_transducer_library.m` |
| `functions/hpc/` | [SLURM](https://slurm.schedmd.com/) and qsub submission | `hpc_submit_job.m` |
| `functions/report/`, `group/` | Subject and group reports | `generate_simulation_report.m` |
| `functions/helper/` | Shared utilities, cache signatures, sidecars | `cache_signature.m` |
| `config/` | Defaults, equipment library, calibrations | `config_default.yaml` |
| `nexus/` | The Nexus (Python, [PyQt6](https://www.riverbankcomputing.com/software/pyqt/)) | `prestus_nexus/app.py` |
| `tests/` | MATLAB tests, grouped in levels | `run_all_tests.m` |
| `tests/docs/` | Doc generator and checkers run in CI | |
| `doc/` | This documentation ([MkDocs](https://www.mkdocs.org/)) | `mkdocs.yaml` for the navigation |

[Functions](doc_functions.md) lists every MATLAB function with a one-line description.

## Steps of a Run

1. **Load.** `load_parameters` merges the toolbox defaults, an optional project copy of the defaults, the study YAML and any overrides, in that order. It then fills in derived values, such as the transducer geometry from the equipment library, and validates the result.
2. **Decide where to run.** `prestus_pipeline_start` runs the pipeline in the current MATLAB session or submits it as a cluster job. Some modes (several transducers fired separately, several intensities, uncertainty) split into chained jobs.
3. **Place.** The placement stage turns the chosen placement method into transducer and focus positions in T1 voxels. This is the only place where RAS millimetres become voxels.
4. **Build the head model.** The grid stage crops and resamples the head onto the simulation grid; the medium stage maps tissue labels (or a pseudo-CT) to sound speed, density and attenuation.
5. **Simulate sound.** The source stage builds the transducer on the grid. A free-water run of the same source gives the baseline intensity, and the acoustic stage calls the chosen solver (k-Wave, [jwave](https://github.com/ucl-bug/jwave) or [BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain)) and scales the field to the target ISPPA.
6. **Simulate heat.** The thermal stage runs the sonication protocol on the scaled field.
7. **Write outputs.** Analysis functions write the metrics CSV, [NIfTI](https://nifti.nimh.nih.gov/) maps and the HTML report to the subject's output folder.

The stage table in [How PRESTUS Works](doc_architecture.md#pipeline-stages) gives each stage's functions, inputs and cache files. Each stage can be switched off with a `modules.*` flag ([Modules](doc_modules.md)).

## Required Updates per Change

| To | Edit | Also update | Checked by |
|---|---|---|---|
| Add a config key | `config/config_default.yaml`, with a comment (the Nexus shows it as help) | a row in [Parameters](doc_parameters.md) | `check_param_docs.py`, Nexus `test_schema.py` |
| Read a key the defaults do not define | the reading function | an entry under `extra:` in `nexus/prestus_nexus/overlay.yaml` | Nexus `test_schema.py` |
| Rename a released key | `config_default.yaml` | an entry under `deprecated:` in `overlay.yaml`; `upgrade_config_keys` then converts old configs | Nexus `test_schema.py` |
| Change what a cached stage computes | the stage function | `cache_signature.m` for that stage, so old caches are rejected | `test_cache_signature.m` |
| Add an output a later run can start from | the writing stage | a JSON sidecar via `write_sidecar` | `check_sidecar_signature` on reuse |
| Add a function | a new file with a first help line `% NAME  Short description` | nothing: the function index is generated | `gen_function_index.py --check` |
| Add a test | `tests/test_*.m` | the test's level in `run_all_tests.m` | `test_suite_wiring` |
| Change a coordinate transform | read [Coordinate Systems](doc_coordinate_systems.md) first | | unit and demo tests |
| Change results | | a line in the changelog with ⚠️ and the action users must take | physics and golden-metric tests |

Docs describe current behaviour as tasks, and the changelog has one line per change ([Contributing](CONTRIBUTING.md)).

## Pre-Commit Checks

Run the checks CI runs before you commit:

```bash
python tests/docs/gen_function_index.py --check
python tests/docs/check_param_docs.py
python tests/docs/check_links.py
mkdocs build --strict -d /tmp/prestus-site
(cd nexus && python -m pytest -q -m "not matlab")
matlab -batch "cd tests; run_all_tests('unit')"
```

Unit tests need no data. A change to the physics needs `run_all_tests('physics')` as well; a change to the pipeline needs the demo levels, which use the demo subject. [Testing](doc_testing.md) lists every level, the data it needs and its run time.

## Nexus Code Structure

The Nexus builds its forms from `config_default.yaml`: every key becomes a field and its comment becomes the help text (`schema.py`). `overlay.yaml` adds what the YAML cannot say: widget types, choices, units, advanced flags, and the extra and deprecated keys. A new config key therefore appears in the Nexus without Python changes.

| File | Does |
|---|---|
| `model.py`, `merge.py`, `yamlio.py` | Hold the study config and merge and write it the way MATLAB reads it |
| `checks.py` | The Issues panel: checks run before MATLAB sees the config |
| `quickstart.py` | The QuickStart tab's logic |
| `launcher.py`, `remote.py`, `jobs.py` | Run locally, submit over SSH, follow cluster jobs |
| `bridge.py` | The persistent MATLAB session for the Calibration and Localite pages (`functions/gui/nexus_bridge.m`) |
| `results.py` | Finds the runs in a subject's output folder |
| `ui/` | One module per tab or page; `main_window.py` assembles them |

The Nexus tests run without MATLAB (`pytest -m "not matlab"`); the tests marked `matlab` check the round trip through MATLAB.

## Further Reading

- [How PRESTUS Works](doc_architecture.md): config layers, chained jobs, stage inputs and outputs, caches, intensity scaling.
- [Coordinate Systems](doc_coordinate_systems.md): before touching any position or transform.
- [Testing](doc_testing.md): test levels and how to add a test.
- [Contributing](CONTRIBUTING.md): issues, pull requests, doc and changelog style.
