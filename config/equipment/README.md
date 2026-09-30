# Equipment Configuration Files

This folder contains transducer YAML files and supporting data used by the calibration pipeline.

---

## Transducer YAML (`<SERIAL>.yaml`)

One file per physical transducer unit. Loaded by `load_equipment_config()`.

```yaml
type: transducer
serial: CTX_250_001              # unique identifier; must match the filename stem
name: NeuroFUS 4 ch. CTX_250_001 # human-readable label
manufact: Sonic Concepts          # informational label, e.g. 'Sonic Concepts' | 'Imasonic'
n_elem: 4                         # number of physical elements
min_foc: 13.79                    # minimum steerable focal depth from exit plane [mm]
max_foc: 61.48                    # maximum steerable focal depth from exit plane [mm]

transducer:
  annular:
    elem_n: 4                          # virtual element count used in simulation
    elem_id_mm: [0.0, 32.9, 46.1, 56.0]  # inner diameters of each ring [mm]
    elem_od_mm: [32.4, 45.6, 55.5, 64.0] # outer diameters of each ring [mm]
    curv_radius_mm: 63.2               # radius of curvature of the bowl [mm]
    dist_geom_ep_mm: 52.38             # distance from exit plane to geometric focus (centre of curvature) [mm]
  freq_hz: 250000.0                    # operating frequency [Hz]

combos:
  <combo_name>:                        # arbitrary key used in calibration config
    ds:
      serial: SC_203_035
      name: TPO junior SC_203_035
      manufact: Sonic Concepts
      available_ch: 4
    char_data_path: <TRAN>/tran-<TRAN>_ds-<DS>_axial.csv  # measured axial profiles, relative to local/profiles/
```

### Matrix transducers

Instead of `transducer.annular`, a `transducer.matrix` block uses the fields of `transducer.matrix` in `config_default.yaml`; `calibration_gui` (**New…**) writes this form. `n_elem` is the element count. Matrix entries are geometry only in the calibration pipeline.

```yaml
transducer:
  matrix:
    steering: 3D                 # 1D | 3D
    elem_shape: rect             # rect | disc | bowl
    elem_width_mm: 3.0
    elem_height_mm: 3.0
    outer_diameter_mm: 60.0
    is_curved: true
    curv_radius_mm: 62.94
    dist_geom_ep_mm: 52.38
    matrix_shape:
      type: define_here
      define_here:
        grid_shape:
          type: fibonacci        # rect | fibonacci | fermat
          fibonacci: {elem_n: 256, kerf_mm: 0.5}
  freq_hz: 250000.0
```

### Steering phases

Steering phases are computed from the ring geometry in this YAML (`elem_id_mm`, `elem_od_mm`, `curv_radius_mm`, `dist_geom_ep_mm`) by `set_real_phases` / `compute_phases`, for every manufacturer. No manufacturer phase table is read, and combos have no `phase_table` field. Any systematic difference between geometric and hardware phases is absorbed by the fitted per-element correction (see below and `doc/doc_transducer_library.md`).

### Optional: per-element hardware correction

Written by `save_elem_correction` after a reference-depth calibration. Enables geometric-steering mode for subsequent depths.

```yaml
elem_phase_correction:
  ref_depth_ep_mm: 40.0
  deg: [0.0, -12.3, 5.7, 2.1]   # one value per physical element
```

When present, `calibration_setup` uses `compute_phases(depth) + correction` instead of a global search.

---

## Axial Profile CSV (`local/profiles/<TRAN>/tran-<TRAN>_ds-<DS>[_acq-<LABEL>]_axial.csv`)

Hydrophone measurements of on-axis intensity for one transducer–driving-system combination. They are site data, so they live in the git-ignored `local/profiles/` folder, not here. Each CSV has a YAML sidecar with the same stem.

```
,Focal depths (mm),,,,,        # title row (skipped)
Distance WRT Exit plane (mm),10,15,20,...   # focal depths from exit plane [mm], cols 2…N
0,4.38,1.45,0.22,...           # distance from exit plane [mm] | intensity [W/cm²] per depth
0.5,4.74,2.17,0.18,...
...
```

Read by `read_axial_profiles`. The folder layout, sidecar fields, focal-label conventions and format rules are in `doc/doc_transducer_library.md` (Appendix).

---

## Legacy phase-table examples

`example_phase_table_SC.csv` (Sonic Concepts steering table) and `example_phase_table_IS.ini` (Imasonic element positions) are kept for reference only. PRESTUS no longer reads either format.

---

## Adding a new transducer

1. Create `<SERIAL>.yaml` following the template above.
2. Check that `elem_id_mm`, `elem_od_mm`, `curv_radius_mm`, and `dist_geom_ep_mm` are correct; steering phases depend on them.
3. Add an axial-profile CSV (+ sidecar) per driving-system combination under `local/profiles/<SERIAL>/`, and set `char_data_path` to it.
4. Register the combo in the calibration config (`parameters.calibration.combinations`).
