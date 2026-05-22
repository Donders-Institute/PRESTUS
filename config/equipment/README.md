# Equipment Configuration Files

This folder contains transducer YAML files and supporting data used by the calibration pipeline.

---

## Transducer YAML (`<SERIAL>.yaml`)

One file per physical transducer unit. Loaded by `load_equipment_config()`.

```yaml
type: transducer
serial: CTX_250_001              # unique identifier; must match the filename stem
name: NeuroFUS 4 ch. CTX_250_001 # human-readable label
manufact: Sonic Concepts          # 'Sonic Concepts' | 'Imasonic' | 'generic' (or any string)
n_elem: 4                         # number of physical elements
min_foc: 13.79                    # minimum steerable focal depth from exit plane [mm]
max_foc: 61.48                    # maximum steerable focal depth from exit plane [mm]

transducer:
  annular:
    elem_n: 4                          # virtual element count used in simulation
    elem_id_mm: [0.0, 32.9, 46.1, 56.0]  # inner diameters of each ring [mm]
    elem_od_mm: [32.4, 45.6, 55.5, 64.0] # outer diameters of each ring [mm]
    curv_radius_mm: 63.2               # radius of curvature of the bowl [mm]
    dist_geom_ep_mm: 52.38             # axial distance from bowl apex to exit plane [mm]
  freq_hz: 250000.0                    # operating frequency [Hz]

combos:
  <combo_name>:                        # arbitrary key used in calibration config
    ds:
      serial: SC_203_035
      name: TPO junior SC_203_035
      manufact: Sonic Concepts
      available_ch: 4
    char_data_path: Axial_profiles_<TRAN>~<DS>.csv  # empirical axial intensity profiles
    phase_table: <file>                # see manufacturer-specific rules below
```

### `manufact` and `phase_table` rules

| `manufact`      | `phase_table` value       | File format       | How phases are computed              |
|-----------------|---------------------------|-------------------|--------------------------------------|
| `Sonic Concepts`| path to CSV (required)    | see SC section    | looked up by focal depth             |
| `Imasonic`      | path to INI (required)    | see IS section    | `compute_phases` using INI positions |
| anything else   | `''` or omit field        | not used          | `generate_tran_ini_from_geometry` derives positions from ring geometry in this YAML, then calls `compute_phases` |

### Optional: per-element hardware correction

Written by `save_elem_correction` after a reference-depth calibration. Enables geometric-steering mode for subsequent depths.

```yaml
elem_phase_correction:
  ref_depth_ep_mm: 40.0
  deg: [0.0, -12.3, 5.7, 2.1]   # one value per physical element
```

When present, `calibration_setup` uses `compute_phases(depth) + correction` instead of a global search.

---

## Axial Profile CSV (`Axial_profiles_<TRAN>~<DS>.csv`)

Empirical hydrophone measurements of on-axis intensity for one transducer–driving-system combination.

```
% optional comment rows (any number; readmatrix skips leading text rows)
% second comment row
0,0,0,...                    # filler row (parsed as charac_data(1,:), ignored)
0,40.0,55.0,70.0,...         # focal depths from exit plane [mm] (charac_data(2,:))
5.0,0.12,0.05,0.01,...       # distance from exit plane [mm] | intensity columns
10.0,0.45,0.20,0.04,...
...
```

- Column 1: axial distance from exit plane [mm]
- Columns 2…N: measured intensity at each focal depth (arbitrary units; scaled later)
- The focal-depth row (`charac_data(2, 2:end)`) must match the depths listed in the calibration config
- Read by `readmatrix`; leading rows starting with `%` or non-numeric text are skipped automatically

---

## Sonic Concepts Phase Table CSV

Required when `manufact: Sonic Concepts`. Maps focal depth to manufacturer-specified driving phases.

```
Distance,Phase_El1,Phase_El2,Phase_El3,Phase_El4,Gain
40.0,0.0,330.0,300.0,275.0,1.0
55.0,0.0,345.0,323.0,308.0,1.0
70.0,0.0,0.0,0.0,0.0,1.0
```

- Column 1: focal depth from exit plane [mm] — must exactly match the requested depth (no interpolation)
- Columns 2…N−1: phases in degrees [0–360) for each element
- Last column: ignored (gain or notes)
- Code extracts phases as `phase_table(:, 2:end-1)`

---

## Imasonic Element-Position INI

Required when `manufact: Imasonic`. Provides 3-D positions of element representative points for `compute_phases`.

```ini
[elements]
x1=5.0|0.0|-69.82
x2=12.125|0.0|-68.94
x3=15.75|0.0|-68.20
x4=18.75|0.0|-67.44
```

**Coordinate system:**
- Origin: natural focus (centre of curvature of the bowl)
- Z-axis: axial, positive toward the target
- All elements have Z < 0 (they sit behind the focal plane)

**Formula for a spherical bowl** with radius of curvature R [mm] and element mid-radius r_mid [mm]:

```
r_mid = (elem_id_mm/2 + elem_od_mm/2) / 2
z     = -sqrt(R² - r_mid²)
```

For `generic` transducers these positions are derived automatically from the ring geometry in the YAML by `generate_tran_ini_from_geometry`; no INI file is needed.

---

## Adding a new transducer

1. Create `<SERIAL>.yaml` following the template above.
2. Set `manufact` to `Sonic Concepts`, `Imasonic`, or `generic`.
3. For `Sonic Concepts`: add a phase-table CSV and reference it under `phase_table`.
4. For `Imasonic`: add an element-position INI and reference it under `phase_table`.
5. For `generic`: leave `phase_table: ''`; ensure `elem_id_mm`, `elem_od_mm`, and `curv_radius_mm` are correct.
6. Add an axial-profile CSV measured for each driving-system combination.
7. Register the combo in the calibration config (`parameters.calibration.combinations`).
