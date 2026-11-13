# Testing

PRESTUS is tested at five levels, from checks that take seconds and run on every change to benchmarks that need a GPU or a real subject and run before a release. This page tells you which level a change needs, what data, hardware and time each level takes, and how to add a test.

On this page:

- [Level 0: Static Checks](#level-0-static-checks): the doc, parameter, link and Nexus checks that run on every push.
- [Level 1: Unit Tests](#level-1-unit-tests): data-free tests of each function's logic, grouped by area.
- [Level 2: Physics Benchmarks](#level-2-physics-benchmarks): free-water fields against the [O'Neil](https://doi.org/10.1121/1.1906542) solution and the [ITRUSST](https://itrusst.org/) benchmarks, with their tolerances.
- [Level 3: Pipeline Tests](#level-3-pipeline-tests): the whole pipeline on synthetic water, a segmented subject and the demo subject, with golden metrics.
- [Level 4: Backends and the Nexus Round Trip](#level-4-backends-and-the-nexus-round-trip): [jwave](https://github.com/ucl-bug/jwave) and [BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain) against [k-Wave](https://github.com/ucl-bug/k-wave), and Nexus configs loaded in MATLAB.
- [Release Checklist](#release-checklist): the tests to run before tagging a release.
- [Rules for New Tests](#rules-for-new-tests): where a new test goes and the rules it follows.

The [MATLAB](https://www.mathworks.com/products/matlab.html) tests use [`matlab.unittest`](https://www.mathworks.com/help/matlab/matlab-unit-test-framework.html) and live in `tests/`; the Nexus tests use [pytest](https://docs.pytest.org/) and live in `nexus/tests/`. `run_all_tests(level)` runs a MATLAB level and fails (non-zero exit under `matlab -batch`) when any test fails.

| Level | What it guarantees | Command | Needs | Time | Runs |
|---|---|---|---|---|---|
| 0 Static | Docs match the code, the site builds, the Nexus works | see [Level 0](#level-0-static-checks) | [Python](https://www.python.org/) | ~2 min | CI, every push |
| 1 Unit | Every function's logic, incl. config, placement, scaling, caches, HPC parsing, reports | `run_all_tests('unit')` | MATLAB | ~15 min (macOS, checkout on a network drive; less on a local disk) | CI, every push |
| 2 Physics | The solver path reproduces analytic and published reference fields | `run_all_tests('physics')`, `('itrusst')` | MATLAB (+ ITRUSST data) | ~5 min CPU (2D), plus the unit tests, which `run_all_tests('physics')` also runs | CI weekly; 3D on a GPU node before a release |
| 3 Pipeline | Whole pipeline on synthetic and real data, numbers stable against a reference | `run_all_tests('water')` … `('demo_golden')` | k-Wave CPU; [SimNIBS](https://simnibs.github.io/simnibs/build/html/index.html) outputs; demo data | 10 min – hours | before a release |
| 4 Backends | jwave / BabelBrain agree with k-Wave; Nexus ↔ MATLAB round trip | `runtests('tests/test_jwave.m')` etc. | [conda](https://docs.conda.io/) envs, GPU, MATLAB | 10–60 min | after backend or environment changes |

---

## Level 0: Static Checks

```bash
python tests/docs/gen_function_index.py --check   # doc_functions.md lists every function
python tests/docs/check_param_docs.py             # every config key documented; defaults match
python tests/docs/check_links.py                  # relative links and anchors resolve
mkdocs build --strict -d /tmp/prestus-site        # site builds without warnings (mkdocs < 2)
(cd nexus && python -m pytest -q -m "not matlab")
```

`gen_function_index.py` without `--check` regenerates `doc/doc_functions.md` (it keeps written descriptions and takes new ones from each file's help line). `check_param_docs.py` compares the Default column of `doc_parameters.md` with `config_default.yaml`. The Nexus tests include the schema drift test (`test_schema.py`: overlay vs. `config_default.yaml` and `functions/`), which also guards the list of extra config keys MATLAB uses for its unknown-key check.

---

## Level 1: Unit Tests

```matlab
cd tests
run_all_tests          % same as run_all_tests('unit')
```

No data, SimNIBS, k-Wave simulation or GPU is needed. A failed test prints its full diagnostic; a skipped test prints why it was skipped.

| Area | Files |
|---|---|
| Configuration | `test_load_parameters`, `test_config_layering` (toolbox → project → study layering, stale project defaults, unknown keys, non-interactive warnings, `prestus_config_upgrade`), `test_config_upgrade` (keys from v0.5.0 on converted to current names), `test_study` |
| Intensity scaling | `test_isppa_policy` (`needs_water_baseline`, `isppa_baseline`, `isppa_scale_factor`, `apply_isppa_scaling`), `test_focal_reference` |
| Caches | `test_cache_signature` |
| Transducers & calibration | `test_transducer_library`, `test_transducer_model` (tag `model`), `test_library_conditions`, `test_matrix_transducer`, `test_cw_source_matrix`, `test_calibration_session` (the GUI/Nexus calibration steps, `meta.fit`) |
| Placement & neuronavigation | `test_transducer_placement` (tag `manual`), `test_plantus_placement` (tags `validation`, `dispatch`), `test_neuronav_io`, `test_neuronav_batch`, `test_placement_report`, `test_placement_positions` |
| Coordinates | `test_transform` |
| Head model | `test_head_preprocessing`, `test_integration_headpreproc`, `test_deface_m2m` |
| Thermal | `test_thermal_parameters`, `test_thermal_dose` (ISO CEM43 against [Sapareto & Dewey](https://doi.org/10.1016/0360-3016(84)90379-1) values) |
| Backends (logic only) | `test_babelbrain_insitu_isppa`, `test_babelbrain_thermal_timing`, `test_babelbrain_air_regions` (needs Python with [numpy](https://numpy.org/)/[scipy](https://scipy.org/)/[nibabel](https://nipy.org/nibabel/); skipped otherwise) |
| HPC | `test_hpc_job_report`, `test_parse_memory_tags`, `test_hpc_submit` |
| Security | `test_shell_quote_arg`, `test_r_quote_string`, `test_run_command`, `test_command_injection` |
| Nexus MATLAB session | `test_nexus_bridge` (`nexus_bridge` protocol and allow-list, `nexus_neuronav` and `nexus_calibration` on synthetic data) |
| Helpers, reports, suite wiring | `test_helper`, `test_reports`, `test_suite_wiring` |

---

## Level 2: Physics Benchmarks

These check numbers, not just that code runs. Both use a focused bowl in lossless water (the ITRUSST BM1 source: ROC 64 mm, aperture 64 mm, 500 kHz, 60 kPa) run through `prestus_pipeline`.

**`test_physics_water`**: against the analytic O'Neil solution (k-Wave `focusedBowlONeil`), no data needed:

| Check | Tolerance (2D axisymmetric) |
|---|---|
| Focal position on axis | 0.5 mm |
| Focal peak pressure | 2 % |
| −3 dB (intensity) axial focal length | 5 % |
| Four-ring annular array in phase (CTX-500 ring diameters on the same bowl) against O'Neil summed over the rings: focal position / peak pressure / −3 dB length | 1 mm / 2 % / 5 % (2026-10-08: 0.50 mm, +0.10 %, +0.04 %; the summed peak is flat, ±0.5 mm costs 0.08 % pressure) |
| Doubling the drive doubles the pressure | 1 % |
| A target ISPPA on an uncalibrated transducer is met in water | 2 % |
| With `run_water_baseline: 0`, an unmet target stops the run | error `prestus:isppa:noBaseline` |

**`test_itrusst_benchmark`**: benchmarks PH1-BM1-SC1 (focused bowl) and PH1-BM1-SC2 (20 mm plane piston, a [flat annular element](doc_transducer.md#flat-elements)) of [Aubry et al. (2022)](https://doi.org/10.1121/10.0013426): the pressure amplitude in the central plane against the published k-Wave reference (60 points per wavelength). Download `INTERCOMPARISON-PH1-DATA-V1.0` from [Zenodo](https://doi.org/10.5281/zenodo.6020543) and set `PRESTUS_ITRUSST_DATA` to it; the test skips otherwise.

| Metric | Bowl: tolerance (2D) | Bowl: result 2026-10-02 | Piston: tolerance (2D) | Piston: result 2026-10-02 |
|---|---|---|---|---|
| Focal pressure | 2 % | +0.2 % | 2 % | 0.0 % |
| Focal position shift | 0.5 mm | 0.0 mm | 0.5 mm | 0.0 mm |
| L∞ (max abs. difference / reference max) | 4 % | 2.1 % | 3 % | 1.6 % |
| L2 (relative norm) | 3 % | 1.5 % | 1 % | 0.4 % |

For the piston, every on-axis near-field maximum has the same amplitude (2 p₀), so the focal pressure and position are taken beyond half the near-field distance (a²/2λ), which leaves only the last axial maximum.

```matlab
run_all_tests('physics')       % unit + test_physics_water (2D axisymmetric, ~5 min CPU)
setenv('PRESTUS_ITRUSST_DATA', '/path/to/INTERCOMPARISON-PH1-DATA-V1.0')
run_all_tests('itrusst')       % 2D BM1 bowl + piston, 3D BM1 bowl
run_all_tests('physics_3d')    % 3D O'Neil check
```

The 3D variants (tags `physics_3d`, `itrusst_3d`; ~6 M voxels at 0.5 mm) belong on a GPU node: set `PRESTUS_PHYSICS_CODE_TYPE` (e.g. `cpp_gpu`) and submit the MATLAB call as a job. Their tolerances are provisional (focal pressure 5 %, position 1 mm; ITRUSST L∞ 10 % and L2 6 %; O'Neil focal length 10 %) until their first run on the reference node.

**Tolerance policy.** A benchmark's tolerances are set from its first passing run on the reference machine plus a margin, recorded with that run's values in the test file, and then kept fixed. Loosen a tolerance only together with an explanation in the changelog; a benchmark that starts failing after a change is a finding, not a test to adjust.

---

## Level 3: Pipeline Tests

### Synthetic Water Pipeline

```matlab
run_all_tests('water')     % unit + slower unit checks + 3D water pipeline + multi-ISPPA
```

`test_integration_water` runs the full pipeline with `simulation.medium: water` on a small 3D grid (CPU k-Wave), `test_isppa_rescaling` the multi-ISPPA pipeline. The level also runs the slower data-free checks: `test_transducer_model` tag `fit` (optimiser) and `test_plantus_placement` tag `integration` (mocked [PlanTUS](https://github.com/mlueckel/PlanTUS) subprocess; needs Python 3). Fixture: `fixtures/make_minimal_parameters.m` (single-element bowl, drive level, no target ISPPA), `fixtures/make_synthetic_segmentation.m`.

### Subject Smoke Tests

`test_integration_pipeline` reaches checkpoints on a SimNIBS-segmented subject; each level includes the ones above it:

| Level | Coverage | Duration |
|---|---|---|
| `run_all_tests('head')` | `smoke_config`, `smoke_head` (head preprocessing to medium masks) | 1–5 min |
| `run_all_tests('acoustic')` | + `smoke_acoustic` | 10–60 min |
| `run_all_tests('all')` | + `smoke_thermal` | 10–60 min |

```bash
export PRESTUS_TEST_DATA=/path/to/demo/data      # folder containing m2m_sub-001/
export PRESTUS_DEMO_CONFIG=/path/to/config.yaml   # optional; default config/config_tutorial.yaml
```

### Demo Subject (All Placement Modes)

`test_integration_demo` uses demo subject sub-009 (SimNIBS segmentation, UTE-derived pseudo-CT, [Localite](https://www.localite.de/) data) in a folder outside the repository, `PRESTUS_DEMO_DATA` (check its layout with `setup_demo_data('/path/to/prestus_testdata')`):

```
prestus_testdata/
├── bids/sub-009/anat/           sub-009_T1w.nii.gz, sub-009_UTE.nii.gz
├── simnibs/m2m_sub-009/          SimNIBS output (incl. pseudoCT, MNI transforms)
├── localite/sub-009/ses-02/localite/Session_*/TMSTrigger/TriggerMarkers_Coil0*.xml
├── coords/sub-009_heuristic_target.json   {"target_name": "MD_thalamus", "mni_target_mm": [-6,-16,4]}
├── configs/config_demo_{manual,localite,heuristic}.yaml
├── golden/sub-009_manual_metrics.json     reference metrics (see below)
└── sim_outputs/                  written by the tests
```

| Level | Checks | Duration |
|---|---|---|
| `demo_inputs` | input files exist | ~1 s |
| `demo_config` | the three placement configs load | ~1 s |
| `demo_localite` | Localite XML → finite voxel positions | ~1 s |
| `demo_pct` | pseudo-CT readable | ~1 s |
| `demo_head` | head preprocessing (manual, pseudo-CT) | 1–5 min |
| `demo_acoustic` | full acoustic pipeline (3D) | 10–60 min (GPU) |
| `demo_thermal` | thermal on the cached acoustics | 10–60 min |
| `demo_golden` | key metrics vs. the stored reference | seconds |

**Golden metrics.** `test_demo_golden` reads the result table of the manual-placement demo run and compares `Isppa_brain` (5 %), `Psptp_brain` and `MI_brain` (3 %), the focus position `I{x,y,z}_brain_vox` (1 voxel), `maxT_brain` (0.2 °C) and `CEM43iso_brain` (20 %) with `golden/sub-009_manual_metrics.json`. Without a reference it writes `…metrics.candidate.json` and reports the test as incomplete: review the candidate (e.g. against the last release), then rename it. Update the reference only for intended changes and say so in the changelog. Run it after `demo_acoustic` and `demo_thermal`.

---

## Level 4: Backends and the Nexus Round Trip

| File | Checks | Needs |
|---|---|---|
| `smoke_test_backends.m` | BabelBrain and jwave importable, bridge scripts present (<30 s) | conda envs |
| `test_babelbrain.m` | staged: paths, input writing, transducer mapping, dry run, full FDTD (tags `quick`, `nogpu`, `full`) | `BABELBRAIN_PYTHON`, GPU for `full` |
| `test_jwave.m` | jwave (Helmholtz) vs k-Wave on the same water simulation: peak, FWHM, field correlation (tag `quick` without the timing sweep) | `JWAVE_PYTHON` |
| `nexus/tests/test_matlab_roundtrip.py` | configs written by the Nexus load in `load_parameters` with identical values | MATLAB (`PRESTUS_MATLAB`) |
| `nexus/tests/test_real_configs.py` | every study config in `$PRESTUS_REAL_CONFIG_DIR` round-trips | real configs |
| `nexus/tests/test_bridge.py` (matlab-marked part), `test_bridge_pages_matlab.py` | the real `nexus_bridge`; the Localite and Calibration pages end to end (synthetic Localite study; CTX_250_001 profiles from `local/profiles`) | MATLAB |

```matlab
runtests('tests/test_babelbrain.m', 'Tag', 'nogpu')
runtests('tests/test_jwave.m', 'Tag', 'quick')
```

```bash
cd nexus && python -m pytest -q      # includes the MATLAB round trip when MATLAB is found
```

---

## Release Checklist

1. CI green (levels 0–1, and the latest weekly physics run).
2. `run_all_tests('itrusst')` and `run_all_tests('physics_3d')` on the GPU node.
3. `run_all_tests('demo_acoustic')`, `run_all_tests('demo_thermal')`, `run_all_tests('demo_golden')`.
4. `test_babelbrain` (`full`) and `test_jwave` if a backend or its environment changed.
5. Nexus round trip (`cd nexus && python -m pytest -q`).

---

## Rules for New Tests

- **New function or fix** → a `methods (Test)` block in the `test_*.m` of the function's area; a new file only for a new area. `run_all_tests` finds every test class in `tests/` and runs a test at the levels its level tag belongs to (`prestus_test_levels.m`). Give each class the shared path fixture and one level tag, and add no `addpath` of its own:

    ```matlab
    classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_x < matlab.unittest.TestCase
    ```

    A file that spans levels tags each `methods (Test, TestTags = {...})` block instead; `standalone` marks tests no level runs. `test_suite_wiring` fails for a test without exactly one level tag.
- **Anything that changes numbers** → extend a physics benchmark, or the golden metrics, so the change is visible.
- **New pipeline stage** → a tagged method in `test_integration_pipeline.m` and, if it warrants it, a level in `prestus_test_levels.m`.
- **New demo subject or placement mode** → extend `test_integration_demo.m` with a private loader like `load_demo_config`.
- Tests write only to `tempdir` (or `PRESTUS_DEMO_DATA/sim_outputs`) and skip, not fail, when optional data or tools are missing (`assume*`).
