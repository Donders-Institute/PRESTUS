# Testing

PRESTUS uses MATLAB's built-in [`matlab.unittest`](https://www.mathworks.com/help/matlab/matlab-unit-test-framework.html) framework. Tests live in `tests/` at the repo root and are organised into two tiers.

---

## Tier 1 — Unit tests

Unit tests cover pure, deterministic functions — no external tools, no file I/O, no k-Wave. They run in seconds on any machine.

| File | Functions covered |
|---|---|
| `test_helper.m` | `round_if_integer`, `find_min_factor`, `get_crop_dims`, `masked_max_3d`, `charm_seg_labels`, `get_flhm_center_position`, `cast_struct`, `get_xyz_mesh`, `zip_fields`, `subset_fields` |
| `test_thermal_parameters.m` | `thermal_parameters` — duty cycle, pulse counts, step discretisation, validation errors |
| `test_transform.m` | `ras_to_grid`, axisymmetric round-trip size checks |
| `test_load_parameters.m` | `load_parameters` — default keys, config merging, affix sanitisation |
| `test_head_preprocessing.m` | `get_crop_dims`, `preproc_medium_mask`, `skull_fill_holes` on synthetic segmentation volumes |

Run all unit tests from MATLAB:

```matlab
run_all_tests          % equivalent to run_all_tests('unit')
```

From the command line (e.g. in CI):

```bash
matlab -batch "run_all_tests"
```

---

## Tier 2 — Integration / smoke tests

Integration tests verify that pipeline stages reach expected output checkpoints on real demo data. They are tagged by pipeline depth so only the relevant subset runs.

| Tag | Pipeline coverage | Typical duration |
|---|---|---|
| `smoke_config` | Config loading only | seconds |
| `smoke_head` | Head preprocessing up to medium masks | 1–5 min |
| `smoke_acoustic` | Full acoustic simulation | 10–60 min |
| `smoke_thermal` | Thermal simulation (reuses cached acoustic output) | 10–60 min |

Each level is a superset of the one above; `run_all_tests('acoustic')` also runs `smoke_config` and `smoke_head`.

### Prerequisites

1. Run SimNIBS `charm` segmentation for the demo subject (or set `modules.segmentation_only = 1` to do this in isolation; see [Modules](doc_modules.md)).
2. Set two environment variables before starting MATLAB:

```bash
export PRESTUS_TEST_DATA=/path/to/demo/data     # folder containing m2m_sub-001/
export PRESTUS_DEMO_CONFIG=/path/to/config.yaml  # optional; defaults to tutorial_config.yaml
```

### Running a specific level

```matlab
run_all_tests('head')      % unit + smoke_config + smoke_head
run_all_tests('acoustic')  % ... + smoke_acoustic
run_all_tests('all')       % full suite
```

Or target a single tag directly:

```matlab
runtests('tests/test_integration_pipeline.m', 'Tag', 'smoke_head')
```

Tests that require demo data skip gracefully (rather than failing) when `PRESTUS_TEST_DATA` is not set, so the unit test level is always safe to run in CI without data.

---

## What is not tested

The following are intentionally out of scope for automated tests and are validated through the demo/tutorial runs instead:

- Acoustic and thermal k-Wave simulations (require licensed k-Wave and GPU/CPU compute)
- SimNIBS segmentation calls
- HPC job submission
- Plot and visualisation functions
- NIfTI read/write as primary logic

---

## Adding new tests

- **New pure function** → add a `methods (Test)` block to the relevant `test_*.m` file, or create a new class file following the same `matlab.unittest.TestCase` pattern.
- **New pipeline stage** → add a tagged method to `test_integration_pipeline.m` and a corresponding case in `run_all_tests.m` if the stage warrants its own level.
