# Transducer Library

PRESTUS ships with a two-tier library for transducer definitions:

| Folder | Content | Purpose |
|---|---|---|
| `config/equipment/` | One YAML per physical device | Geometry, frequency, hardware corrections |
| `config/transducer/` | One YAML per calibration key | Calibrated phases and amplitude lookup table |

When a study config contains `transducer.serial`, PRESTUS resolves geometry and (if available) calibrated parameters automatically at load time.  Inline fields in the study config always override library values.

### Calibration library keys

Calibrations can be generic (transducer only) or specific to a transducer–driving-system pair:

| Scenario | Library filename | When to use |
|---|---|---|
| Generic / DS-agnostic | `{serial}.yaml` | Calibrated against characterisation data without a specific DS, or DS doesn't matter |
| DS-specific | `{serial}_{ds_serial}.yaml` | Calibration was performed with a specific driving system and phases/amplitude differ per DS |

Specifying `combo.ds_serial` in the study config selects the DS-specific key; omitting it uses the generic key.  Both are optional — if neither library file exists and phases/amplitude are absent, the pipeline errors with instructions.

---

## Study config: minimal vs. fully specified

**Minimal — let the library fill in the rest:**

```yaml
transducer:
  serial: IS_PCD15473_01001
  name: left                      # optional human label
  target_isppa_wcm2: 32.0
  focal_distance_ep: 60.0
  # looks up config/transducer/IS_PCD15473_01001.yaml
```

**With a specific driving system** (uses `{serial}_{ds_serial}.yaml`):

```yaml
transducer:
  serial: IS_PCD15473_01001
  combo:
    ds_serial: IGT_32_ch_comb_10_ch
  target_isppa_wcm2: 32.0
  focal_distance_ep: 60.0
  # looks up config/transducer/IS_PCD15473_01001_IGT_32_ch_comb_10_ch.yaml
```

**Fully specified — bypass the library entirely:**

```yaml
transducer:
  serial: IS_PCD15473_01001       # still useful for provenance
  type: annular
  freq_hz: 300000.0
  focal_distance_ep: 60.0
  annular:
    elem_n: 10
    elem_id_mm: [10.0, 22.1, ...]
    elem_od_mm: [21.1, 28.8, ...]
    curv_radius_mm: 100.0
    dist_geom_ep_mm: 92.7
    elem_phase_deg: [159.2, 159.2, ...]
    elem_amp: [150238.0, ...]
```

---

## Adding a new transducer geometry

1. Create `config/equipment/{SERIAL}.yaml` (replace `{SERIAL}` with the device serial number):

```yaml
type: transducer
serial: MY_TRANSDUCER_001
name: My Transducer 1 ch
manufact: MyManufacturer
n_elem: 4
min_foc: 20
max_foc: 80
freq_hz: 500000.0

transducer:
  freq_hz: 500000.0
  annular:
    elem_n: 4
    elem_id_mm: [0, 34, 53, 70]
    elem_od_mm: [32, 51, 67, 80]
    curv_radius_mm: 63.2
    dist_geom_ep_mm: 62.0

combos:
  MY_DRIVER_001:
    ds:
      serial: MY_DRIVER_001
      name: My Driver
      manufact: MyManufacturer
      available_ch: 4
    char_data_path: Axial_profiles_MY_TRANSDUCER_001~MY_DRIVER_001.csv
    phase_table: phase_table_MY_TRANSDUCER_001.ini
```

2. Verify PRESTUS picks it up:

```matlab
eq = load_equipment_config();
disp(fieldnames(eq.trans))   % should include 'MY_TRANSDUCER_001'
```

The transducer is now available in the GUI serial dropdown and can be referenced by `serial` in study configs.  No calibration is required to run simulations — you can provide phases and amplitude inline.

---

## Adding a calibration (depositing results)

Calibration produces a YAML file in the `calibration.path_output_profiles` folder.  To make it available for automatic resolution:

1. Run the calibration workflow for each required focal depth and intensity:

```matlab
parameters = load_parameters('config_calibration.yaml', '/path/to/config');
calibration_standalone(parameters);
```

   This writes `{combo_name}-F{focal}mm-I{intensity}wpercm2.yaml` to the output folder.

2. Run `update_transducer_library` to convert the per-run YAMLs into the parametric library format:

```matlab
update_transducer_library('MY_TRANSDUCER_001_MY_DRIVER_001', ...
    '/path/to/calibration/output', ...
    fullfile(get_prestus_path(), 'config', 'transducer'));
```

   This creates or updates `config/transducer/MY_TRANSDUCER_001_MY_DRIVER_001.yaml`.

3. Verify:

```matlab
lib = load_transducer_from_library('MY_TRANSDUCER_001_MY_DRIVER_001', 60, 32, ...
    load_equipment_config());
disp(lib.transducer.annular.elem_phase_deg)
```

After step 2, any study config referencing `serial: MY_TRANSDUCER_001` with `target_isppa_wcm2` and `focal_distance_ep` will have phases and amplitude resolved automatically.

---

## Hardware corrections (`elem_phase_correction`)

If element-to-element phase offsets have been measured (e.g., via hydrophone mapping), they can be stored in the equipment YAML:

```yaml
# inside config/equipment/MY_TRANSDUCER_001.yaml
elem_phase_correction:
  ref_depth_ep_mm: 60.0
  deg: [0.5, -1.2, 0.8, -0.3]
```

`load_transducer_from_library` applies these corrections automatically on top of the calibrated phases.  Use `save_elem_correction` after running the correction measurement pipeline to write this block.

---

## Reproducibility

Every pipeline run writes a resolved-parameters YAML to `{output}/log/sub-NNN_{medium}_parameters{affix}.yaml`.  This file contains the fully expanded transducer struct (geometry + calibrated phases + amplitude) as it was at runtime, so results can be reproduced even if the library is updated later.
