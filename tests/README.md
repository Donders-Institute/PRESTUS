# PRESTUS Test Suite

Tests are organised into four tiers.

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
export PRESTUS_DEMO_CONFIG=/path/to/config_tutorial.yaml   # optional; defaults to config/config_tutorial.yaml
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

## Tier 4 — Demo-data tests (real subject, all placement modes)

File: `test_integration_demo.m`

Uses a real demo subject (sub-009) with pre-computed SimNIBS segmentation, UTE-derived
pCT, and Localite neuronavigation data stored in a standardised folder outside the repo
(`prestus_testdata/`). Three placement modes are covered — manual, Localite, and
heuristic — and pCT-based skull mapping is verified as a separate tag.

| Tag | What it checks | Typical duration |
|---|---|---|
| `demo_inputs` | All required input files exist | ~1 s |
| `demo_config` | All three placement-mode configs load without error | ~1 s |
| `demo_localite` | Localite XML parses to finite voxel positions | ~1 s |
| `demo_pct` | pCT NIfTIs are present and readable | ~1 s |
| `demo_head` | Head preprocessing end-to-end (manual + pCT variants) | 1–5 min |
| `demo_acoustic` | Full acoustic pipeline | 10–60 min |
| `demo_thermal` | Thermal pipeline on cached acoustic outputs | 10–60 min |

### Demo data layout (`prestus_testdata/`)

```
prestus_testdata/
├── bids/sub-009/anat/
│   ├── sub-009_T1w.nii.gz       ← symlink → simnibs/m2m_sub-009/T1.nii.gz
│   └── sub-009_UTE.nii.gz       ← symlink → simnibs/m2m_sub-009/UTE_reg.nii.gz
├── simnibs/
│   └── m2m_sub-009/             ← symlink → SimNIBS output (incl. pCT, MNI transforms)
├── localite/
│   └── sub-009/ses-02/localite/
│       └── Session_*/TMSTrigger/TriggerMarkers_Coil0*.xml
├── coords/
│   └── sub-009_heuristic_target.json   ← {"target_name": "MD_thalamus", "mni_target_mm": [-6,-16,4]}
├── configs/
│   ├── config_demo_manual.yaml
│   ├── config_demo_localite.yaml
│   └── config_demo_heuristic.yaml
└── sim_outputs/                 ← written by the pipeline
```

### Prerequisites

1. Validate the layout using the diagnostic script:

```matlab
setup_demo_data('/path/to/prestus_testdata')
```

2. Set the environment variable before running tests:

```bash
export PRESTUS_DEMO_DATA=/path/to/prestus_testdata
```

### Running

```matlab
run_all_tests('demo_inputs')    % file checks only (~1 s)
run_all_tests('demo_localite')  % inputs + config + XML parsing
run_all_tests('demo_head')      % all fast checks + head preprocessing
run_all_tests('demo_acoustic')  % full pipeline up to acoustics
run_all_tests('demo_thermal')   % thermal on cached acoustics
```

Or run a single tag directly:

```matlab
runtests('tests/test_integration_demo.m', 'Tag', 'demo_localite')
```

Tests skip gracefully when `PRESTUS_DEMO_DATA` is not set.

---

## Adding new tests

- **New pure function** → add a `methods (Test)` block to the relevant `test_*.m` file, or create a new one following the same `matlab.unittest.TestCase` pattern.
- **New pipeline stage** → add a tagged method to `test_integration_pipeline.m` and a case in `run_all_tests.m` if it warrants its own level.
- **New demo subject or placement mode** → extend `test_integration_demo.m`; add a helper method in the `Access = private` section to load and configure the new variant.
