# Transducer Library

PRESTUS ships with a two-tier library for transducer definitions:

| Folder | Content | Purpose |
|---|---|---|
| `config/equipment/` | One YAML per physical device | Geometry, frequency, hardware corrections |
| `config/calibration/` | One YAML per calibration key | Calibrated phases and amplitude lookup table |

When a study config contains `transducer.serial`, PRESTUS resolves geometry and calibrated parameters automatically at load time. Fields already set in the study config always override library values.

---

## Library keys and fallback chain

| Filename | When used |
|---|---|
| `{tran_serial}_{ds_serial}.yaml` | DS-specific empirical calibration (highest accuracy) |
| `{tran_serial}.yaml` | DS-agnostic baseline — geometric phases from manufacturer report, or multi-DS joint fit |

PRESTUS resolves in order:
1. Try `{tran_serial}_{ds_serial}.yaml` — use if found
2. Fall back to `{tran_serial}.yaml` — issue a warning, continue
3. If neither exists — error with instructions

Specifying `combo.ds_serial` in the study config enables step 1; omitting it goes straight to step 2.

---

## Study config: minimal to fully specified

**Minimal — let the library fill in phases and amplitude:**

```yaml
transducer:
  serial: IS_PCD15473_01001
  target_isppa_wcm2: 32.0
  focal_distance_ep: 60.0
  # resolves config/calibration/IS_PCD15473_01001.yaml (DS-agnostic)
```

**With a specific driving system** (preferred for best accuracy):

```yaml
transducer:
  serial: IS_PCD15473_01001
  combo:
    ds_serial: IGT_32_ch_comb_10_ch
  target_isppa_wcm2: 32.0
  focal_distance_ep: 60.0
  # resolves config/calibration/IS_PCD15473_01001_IGT_32_ch_comb_10_ch.yaml
```

**Fully specified — bypass the library entirely:**

```yaml
transducer:
  serial: IS_PCD15473_01001
  type: annular
  freq_hz: 300000.0
  focal_distance_ep: 60.0
  annular:
    elem_n: 10
    elem_id_mm: [10.0, 22.1, ...]
    elem_od_mm: [21.1, 28.8, ...]
    curv_radius_mm: 100.0
    dist_geom_ep_mm: 92.7
    elem_phase_deg: [0, 45, 90, ...]
    elem_amp: [150000, 150000, ...]
```

---

## Shipped combos

<!-- Run validate_transducer_library and calib_library_postfit_benchmark to update this table -->

| Combo | Manufacturer | Depths (mm) | Calibration type | Pre-cal R² | Post-cal R² |
|---|---|---|---|---|---|
| IS_PCD15287_01001_IGT_32_ch_comb_10_ch | Imasonic | TBD | empirical | 0.985 | — |
| IS_PCD15287_01002_IGT_32_ch_comb_10_ch | Imasonic | TBD | empirical | 0.983 | — |
| IS_PCD15473_01001_IGT_32_ch_comb_10_ch | Imasonic | TBD | empirical | 0.973 | — |
| CTX_500_006_SC_203_035 | Sonic Concepts | TBD | empirical | 0.862 | — |

Run `calib_library_summary` then `calib_library_postfit_benchmark` (both in `v05/code/`) to populate the post-calibration column with values for your specific library version.

**Accuracy notes:**

- **Imasonic**: high geometric model accuracy (R² > 0.97) — empirical calibration adds incremental improvement by correcting DS-introduced hardware phase offset.
- **Sonic Concepts CTX-500**: moderate geometric accuracy (R² ~0.86) — empirical calibration improves focal placement.
- **Sonic Concepts CTX-250**: geometric model degrades significantly when steering more than ~30% below the natural focal length. Library calibration substantially improves accuracy in that regime. Depths outside the characterised range use linear extrapolation — a warning is emitted.

---

## How the global model works

Each library YAML stores a `global_model` section built by `build_transducer_library`:

```yaml
meta:
  tran_serial: IS_PCD15473_01001
  calibration_type: empirical      # or: geometric
  created_at: "2026-06-09"

global_model:
  depths_ep_mm: [40, 50, 60, 70, 80]   # characterised depths [mm from exit plane]
  elem_phase_deg:                        # n_depths × n_elem matrix
    - [0.0, 45.2, 90.1, ...]            # row for depth 40 mm
    - [0.0, 48.7, 96.3, ...]            # row for depth 50 mm
    - ...
  amplitude_scaling:
    i32p5: 150430    # ISPPA 32.5 W/cm² → elem_amp [Pa equivalent]
    i65:   212780
```

At runtime `load_transducer_from_library` **linearly interpolates** phases between characterised depths. Amplitude is **sqrt-scaled** when the requested intensity is not in the table (a warning is printed for deviations > 5%).

Per-depth refinement entries (marked `is_refinement: true`) override the global model when the requested depth is within 1 mm.

---

## Calibration levels

### Level 0 — Use as-is (most users)

No action needed. Set `serial`, `focal_distance_ep`, and `target_isppa_wcm2`. The library handles the rest.

### Level 1 — Manufacturer-report baseline (no empirical data required)

Use this when you have a transducer not in the shipped library, or a new driving system, and want a working baseline immediately.

```matlab
equip_path   = fullfile(prestus_root, 'config', 'equipment');
library_path = fullfile(prestus_root, 'config', 'transducer');

% DS-agnostic entry (uses first available phase table):
build_geometric_library_entry('MY_SERIAL', equip_path, library_path)

% DS-specific entry:
build_geometric_library_entry('MY_SERIAL', equip_path, library_path, 'Combo', 'MY_DS_SERIAL')
```

Phases are computed from ring geometry (elem_id_mm, elem_od_mm, curv_radius_mm) at 10 auto-selected depths spanning the focal range. No manufacturer phase table file or water-tank data required. Accuracy equals the geometric model baseline.

### Level 2 — Single-point refinement

Use this when you have a single hydrophone measurement at one focal depth in your own setup.

1. Run the standard single-depth calibration workflow (see [`doc_calibration.md`](doc_calibration.md)).
2. The output YAML is saved to `path_output_profiles/`.
3. Ingest into the library:
   ```matlab
   update_transducer_library(combo_name, path_output_profiles, library_path)
   ```
   The new depth is added as a `is_refinement: true` entry that overrides global-model interpolation at that depth.

### Level 3 — Full multi-depth recalibration (best accuracy)

Use this when you have hydrophone measurements at ≥3 depths in your own water tank.

1. Organise axial profile CSVs in the standard format (see [Appendix: Input file formats](#appendix-input-file-formats-for-level-3)).
2. Add `char_data_path` to your combo entry in `config/equipment/{serial}.yaml`:
   ```yaml
   combos:
     MY_DS_SERIAL:
       ds:
         serial: MY_DS_SERIAL
         name: My Driving System
         manufact: MyCompany
       char_data_path: Axial_profiles_MY_SERIAL~MY_DS.csv
       phase_table: MY_SERIAL_phases.ini   # optional — omit to use geometry-only phases
   ```
3. Edit `v05/code/build_transducer_library.m` to set `COMBOS_TO_RUN` and `DESIRED_INTENSITIES`, then run it.
4. Validate and benchmark:
   ```matlab
   validate_transducer_library(equip_path, library_path)
   % then run v05/code/calib_library_postfit_benchmark.m
   ```

### Level 4 — Contribute back

Open a PR with:
- Pre/post R² table from `calib_library_postfit_benchmark`
- `validate_transducer_library` passes with no errors
- `library_version.yaml` patch-bumped
- Characterisation conditions in the combo YAML `notes` field (water temp, hydrophone model, tank geometry)

---

## Adding a new transducer or driving system

### New driving system for an existing transducer

1. Add a combo entry to `config/equipment/{tran_serial}.yaml`:
   ```yaml
   combos:
     MY_NEW_DS:
       ds:
         serial: MY_NEW_DS
         name: My Driving System
         manufact: MyCompany
       # phase_table: same_as_other_combos.ini   # optional; omit to use geometry-only phases
       # char_data_path: add when empirical data is available
   ```
2. Generate a geometric baseline:
   ```matlab
   build_geometric_library_entry('TRAN_SERIAL', equip_path, library_path, 'Combo', 'MY_NEW_DS')
   ```
3. Use in study configs with `combo.ds_serial: MY_NEW_DS`.
4. Upgrade to Level 3 when empirical data is available.

### New transducer (not in equipment config)

1. Create `config/equipment/{MY_SERIAL}.yaml` — copy an existing file as template.

   Required fields for annular arrays:
   ```yaml
   type: transducer
   serial: MY_SERIAL
   name: Human-readable name
   manufact: Imasonic    # or: Sonic Concepts
   min_foc: 30.0         # minimum steerable focal distance from exit plane [mm]
   max_foc: 90.0
   transducer:
     annular:
       elem_n: 10
       elem_id_mm: [0.0, ...]
       elem_od_mm: [21.0, ...]
       curv_radius_mm: 100.0
       dist_geom_ep_mm: 92.7
     freq_hz: 300000.0
   ```
2. Run `check_equipment_config(equip_path)` to validate geometry and physical plausibility.
3. Proceed from Level 1.

---

## Library limitations

| Limitation | Effect | Mitigation |
|---|---|---|
| Linear phase interpolation | Accuracy degrades between characterised depths for non-linear steering relationships | Add more characterisation depths |
| Amplitude sqrt-scaling fallback | ~5% intensity error per 2× intensity step | Calibrate at each desired intensity |
| Water-tank basis | Does not account for tissue acoustic properties | Combined with PRESTUS forward simulation in tissue medium |
| DS-agnostic fallback | Ignores DS hardware phase offset (typically 0–20°) | Use DS-specific entry for publication results |
| Extrapolation beyond characterised range | Unknown accuracy, warning emitted | Keep `focal_distance_ep` within characterised range |

---

## Appendix: Input file formats for Level 3

One input file is required per transducer–driving-system combo when running `build_transducer_library`; the phase table is optional:

| File | Purpose | Format | Required? |
|---|---|---|---|
| Axial profile CSV | Empirical intensity vs. depth at multiple focal settings | See below | Yes |
| Phase table | Manufacturer-supplied steering phases (improves element-position accuracy) | Manufacturer-specific (`.ini` or `.csv`) | No — geometry-only used when absent |

---

### Axial profile CSV

#### Layout

```
,Focal depths (mm),,,,,,,
Distance WRT Exit plane (mm),40.0,50.0,60.0,70.0,80.0
0.5,0.12,0.31,0.04,0.01,0.0
1.0,0.18,0.44,0.08,0.02,0.0
...
```

**Row 1 (title row):** First cell is a label string (ignored). Remaining cells label each column as focal depth in mm from the transducer exit plane.

**Row 2 onward (data rows):** First column is axial distance from exit plane in mm. Remaining columns are measured intensity in W/cm² at each focal depth.

#### Rules

- Separator: comma
- Encoding: UTF-8 (BOM is tolerated — `readmatrix` strips it automatically)
- Focal depths: any number of columns ≥ 1; more depths improve global model accuracy
- Axial distance: must start at or below the face of the transducer (0 mm is fine; negative values are ignored)
- Intensity values: must be non-negative; `readmatrix` pads ragged CSVs with NaN — trailing NaN columns are trimmed automatically
- No units row: intensity is always interpreted as W/cm²

#### Reading in MATLAB

```matlab
raw = readmatrix(csv_path, 'NumHeaderLines', 1);
focal_depths_ep = round(raw(1, 2:end), 2);   % row 1, cols 2+
dist_ep_mm      = raw(2:end, 1);             % col 1, rows 2+
intens_data     = raw(2:end, 2:end);         % data block

% Drop NaN columns (ragged padding)
valid = isfinite(focal_depths_ep);
focal_depths_ep = focal_depths_ep(valid);
intens_data     = intens_data(:, valid);
```

#### File naming convention

```
Axial_profiles_{TRAN_SERIAL}~{DS_SERIAL}.csv
```

Examples:
- `Axial_profiles_IS_PCD15473_01001~IGT_32_ch_comb_10_ch.csv`
- `Axial_profiles_CTX-250-001~TPO-203-035.csv`

The tilde (`~`) is the legacy separator in filenames; the library YAML keys use underscore (`_`). The filename separator does not affect PRESTUS — only the `char_data_path` field in the equipment YAML matters.

#### Depth coverage recommendations

| Transducer type | Minimum depths | Recommended |
|---|---|---|
| Imasonic (10-element annular) | 3 | 5–7 spaced 10 mm apart |
| Sonic Concepts CTX (4-element) | 4 | 6–8, including natural focal depth |
| Any | — | Cover full intended depth range; include the natural focal depth |

---

### Phase tables

#### Imasonic (`.ini` format)

Imasonic supplies a `.ini` file per transducer (not per driving system). All driving systems using the same physical transducer share the same phase table.

```ini
[Header]
NumberOfChannels=10
Frequency=300000

[Focus_40]
Channel_1=45.2
Channel_2=90.1
...

[Focus_50]
Channel_1=48.7
...
```

Read with:
```matlab
phase_table = read_ini_file(phase_table_path);
```

The file is passed to `set_real_phases(phase_table, tran, depth, params)` which interpolates to the requested focal depth.

#### Sonic Concepts (`.csv` format)

Sonic Concepts supplies a steer table CSV per transducer–DS combo:

```
Focus [mm],Ch1,Ch2,Ch3,Ch4
14.0,0,90,180,270
16.0,0,85,170,255
...
```

Read with:
```matlab
phase_table = readtable(phase_table_path);
```

`set_real_phases` selects the row nearest the requested focal depth by linear interpolation.

#### File placement

Phase tables must be accessible at the path specified in `phase_table` under the combo entry in the equipment YAML. PRESTUS resolves this path from `config_calibration.yaml`:

```yaml
calibration:
  path_input_phase: /path/to/phase_tables/
```

The `phase_table` field in the equipment YAML is a filename relative to that folder:

```yaml
combos:
  IGT_32_ch_comb_10_ch:
    phase_table: transducer_15473_10_300kHz.ini
```

---

### Checklist before running build_transducer_library

- [ ] Axial profile CSV is readable by `readmatrix(path, 'NumHeaderLines', 1)`
- [ ] Row 1 contains focal depth values (not strings)
- [ ] Column 1 contains axial distances in ascending order
- [ ] Intensity values are non-negative
- [ ] Phase table file exists at `path_input_phase/{phase_table}`
- [ ] Equipment YAML validates: `check_equipment_config(equip_path)` passes
- [ ] Focal depths in CSV are within `[min_foc, max_foc]` for the transducer
