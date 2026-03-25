# PRESTUS Test Suite

Tests are organised into two tiers.

## Tier 1 — Unit tests (no data required)

| File | What it covers |
|---|---|
| `test_helper.m` | Pure helper functions: `round_if_integer`, `find_min_factor`, `get_crop_dims`, `masked_max_3d`, `charm_seg_labels`, `get_flhm_center_position`, `cast_struct`, `get_xyz_mesh`, `zip_fields`, `subset_fields` |
| `test_thermal_parameters.m` | Timing arithmetic in `thermal_parameters`: duty cycle, pulse counts, step sums, validation errors |
| `test_transform.m` | Coordinate transforms: `ras_to_grid` with known affines, axisymmetric round-trip size checks |
| `test_load_parameters.m` | Config merging in `load_parameters`: default keys, override behaviour, affix sanitisation |
| `test_head_preprocessing.m` | Head preprocessing on synthetic segmentation volumes: `get_crop_dims`, `preproc_medium_mask`, `skull_fill_holes` |

Run all unit tests from MATLAB:

```matlab
run_all_tests          % or run_all_tests('unit')
```

From the command line (CI):

```bash
matlab -batch "run_all_tests"
```

---

## Tier 2 — Water pipeline test (no data required)

File: `test_integration_water.m`

Runs the full pipeline end-to-end with `simulation.medium = 'water'` on a small synthetic grid using CPU k-Wave. No MRI, no SimNIBS, no GPU needed. Typical runtime: 2–5 min.

```matlab
run_all_tests('water')
```

Fixtures used:
- `fixtures/make_minimal_parameters.m` — builds a valid parameter struct for a single-element bowl transducer on a small grid
- `fixtures/make_synthetic_segmentation.m` — generates a concentric-shell segmentation volume for head preprocessing tests

---

## Tier 3 — Integration / smoke tests (demo data required)

File: `test_integration_pipeline.m`

Tests are tagged by pipeline depth. Each level is a superset of the one above.

| Tag | What it runs | Typical duration | Data needed |
|---|---|---|
| `smoke_config` | Config loading only | seconds | none |
| `smoke_head` | Head preprocessing up to medium masks | 1–5 min | SimNIBS m2m output |
| `smoke_acoustic` | Full acoustic simulation | 10–60 min | T1/T2 + m2m |
| `smoke_thermal` | Thermal simulation (reuses acoustic outputs) | 10–60 min | completed acoustic run |

### Prerequisites

1. Run SimNIBS `charm` segmentation for the demo subject first (or use `modules.segmentation_only = 1`).
2. Set environment variables before starting MATLAB:

```bash
export PRESTUS_TEST_DATA=/path/to/demo/data    # folder containing m2m_sub-001/
export PRESTUS_DEMO_CONFIG=/path/to/tutorial_config.yaml   # optional; defaults to configs/tutorial_config.yaml
```

### Running a specific level

```matlab
run_all_tests('head')      % unit + smoke_config + smoke_head
run_all_tests('acoustic')  % ... + smoke_acoustic
run_all_tests('all')       % full suite
```

Or run a single tag directly:

```matlab
runtests('tests/test_integration_pipeline.m', 'Tag', 'smoke_head')
```

---

## Adding new tests

- **New pure function** → add a `methods (Test)` block to the relevant `test_*.m` file, or create a new one following the same `matlab.unittest.TestCase` pattern.
- **New pipeline stage** → add a tagged method to `test_integration_pipeline.m` and a case in `run_all_tests.m` if it warrants its own level.
