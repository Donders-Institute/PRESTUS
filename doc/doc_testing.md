# Testing

PRESTUS uses MATLAB's built-in [`matlab.unittest`](https://www.mathworks.com/help/matlab/matlab-unit-test-framework.html) framework. Tests live in `tests/` at the repo root and are organised into four tiers.

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
| `test_deface_m2m.m` | `preproc_deface_m2m` — folder creation, sentinel, idempotency, skip paths (uses a stub binary; no real pydeface required) |

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
export PRESTUS_DEMO_CONFIG=/path/to/config.yaml  # optional; defaults to config_tutorial.yaml
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

## Tier 3 — Demo-data tests (real subject, all placement modes)

File: `test_integration_demo.m`

Uses a real demo subject (sub-009) with pre-computed SimNIBS segmentation, UTE-derived pCT,
and Localite TriggerMarkers data, stored in a standardised folder outside the repository
(`prestus_testdata/`). Three coordinate-input modes are covered in dedicated configs —
manual voxel coordinates, Localite XML, and heuristic MNI target — and pCT-based skull
mapping is verified as a separate tag.

| Tag | What it checks | Typical duration |
|---|---|---|
| `demo_inputs` | All required input files are present | ~1 s |
| `demo_config` | All three placement-mode configs load and have the expected structure | ~1 s |
| `demo_localite` | Localite TriggerMarkers XML parses to finite voxel positions | ~1 s |
| `demo_pct` | Pre-computed pCT NIfTIs are present and contain a 3-D volume | ~1 s |
| `demo_head` | Head preprocessing end-to-end: manual placement + pCT variant | 1–5 min |
| `demo_acoustic` | Full acoustic pipeline on sub-009 | 10–60 min |
| `demo_thermal` | Thermal pipeline using cached acoustic outputs | 10–60 min |

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
│   └── sub-009_heuristic_target.json   ← {"target_name":"MD_thalamus","mni_target_mm":[-6,-16,4]}
├── configs/
│   ├── config_demo_manual.yaml         ← manual voxel coordinates
│   ├── config_demo_localite.yaml       ← Localite mode (ses-02, TriggerMarkers)
│   └── config_demo_heuristic.yaml      ← heuristic mode (mni_target_mm injected from JSON)
└── sim_outputs/                        ← written by the pipeline; may be empty initially
```

The heuristic config has `placement.heuristic.mni_target_mm: []` as a placeholder. The test
fixture reads `coords/sub-009_heuristic_target.json` and injects the coordinate at runtime,
keeping the coordinate file as the single source of truth.

### Prerequisites

1. Run the diagnostic script to confirm all inputs are in place:

```matlab
setup_demo_data('/path/to/prestus_testdata')
```

2. Set the environment variable before running tests:

```bash
export PRESTUS_DEMO_DATA=/path/to/prestus_testdata
```

### Running

```matlab
run_all_tests('demo_inputs')    % file presence checks only (~1 s)
run_all_tests('demo_localite')  % inputs + config loading + XML parsing
run_all_tests('demo_head')      % all fast checks + head preprocessing (~5 min)
run_all_tests('demo_acoustic')  % full pipeline up to acoustics
run_all_tests('demo_thermal')   % thermal on cached acoustics
```

Or target a single tag:

```matlab
runtests('tests/test_integration_demo.m', 'Tag', 'demo_pct')
```

All demo tests skip gracefully when `PRESTUS_DEMO_DATA` is not set, so they never block CI.

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
- **New demo subject or placement mode** → extend `test_integration_demo.m`; add a private helper method following the pattern of `load_demo_config` / `load_localite_config` / `load_heuristic_config` to configure the new variant.
