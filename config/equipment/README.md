# Equipment Configuration Files

This folder contains transducer YAML files and supporting data used by the calibration pipeline.

---

## Naming

| Entry | Serial | Example | Geometry from |
|---|---|---|---|
| Device model | `<MODEL>` | `CTX_500`, `IS_PCD15473`, `H246` | the model's specification |
| Measured unit | `<MODEL>_<unit id>` | `CTX_500_026`, `IS_PCD15473_01001` | the unit, with its own calibration |
| Variant of a unit | `<MODEL>_<unit id>_<variant>` | `IS_PCD15473_01001_OPM` (excluded) | the unit in another configuration |

`<MODEL>` is the model designation with `-` and spaces written as `_` (e.g. NeuroFUS CTX-250 2 ch. → `CTX_250_2ch`, BrainSonix BSonix 55 mm → `BSonix_55mm`); Imasonic models keep the `IS_PCD` prefix of their units. Use a unit's entry when it is in the library, since its calibration is specific to it, and the model entry otherwise. A study naming a unit serial that has no entry (`CTX_500_031`) uses the model entry of its longest `<MODEL>_` prefix (`CTX_500`), with a warning (`resolve_transducer_from_serial`). Every entry names its model in `model:` (absent for `CTX_500_006`, a 2-channel CTX-500 with no model entry). The units `CTX_250_009` and `CTX_250_014` predate this scheme: they are of model `CTX_250_2ch`.

## Layout

```
config/equipment/
  equipment_info.yaml
  ds/ds-<DS>_driver.yaml                          # driving systems
  model-<MODEL>/tran-<MODEL>_transducer.yaml      # device model: geometry, steering range, sources
  model-<MODEL>/tran-<SERIAL>_transducer.yaml     # measured units of the model: overrides and combos only
  housing/                                        # housing meshes (Nexus 3D view)
```

Your own entries belong in your own library folder (`transducer.library_path`, `path.transducer_library` or the `PRESTUS_TRANSDUCER_LIBRARY` environment variable; see `doc_transducer_library.md`), not here: equipment files there named `*_transducer.yaml` or `*_driver.yaml` add transducers or replace shipped ones of the same serial (`load_equipment_config(equipment_path, user_folders)`).

`load_equipment_config()` reads every YAML in the folder and its subfolders, by its `type` (`transducer` or `driving_system`); file names are for orientation only. A flat `<SERIAL>.yaml`, as in earlier versions or written by `calibration_gui`, is read as well. A unit (`model:` naming another entry) takes the fields it omits from its model: `transducer` (geometry), `min_foc`/`max_foc`, `n_elem`, `manufact` and, in the Nexus, `housing`. A combo without a `ds` block takes it from the driving-system entry with the same serial.

The calibration library (`config/calibration/`) has the same model folders: `model-<MODEL>/tran-<TRAN>[_ds-<DS>]_calib.yaml` (see `library_entry_file`).

## Transducer YAML

One file per device model or physical transducer unit.

The steering range is a property of the hardware, not of the measurements: the depths with hydrophone profiles (the measured range) are read from the calibration entry. A unit that is retired or replaced gets `status: deprecated`, an entry withheld until it is documented `status: excluded` (not listed in the Nexus or `calibration_gui`, skipped by `build_transducer_library`). Both stay in the library so earlier studies re-run unchanged, and PRESTUS and the Nexus warn when a study names them.

```yaml
type: transducer
serial: CTX_250_001              # unique identifier; must match the filename stem
name: NeuroFUS 4 ch. CTX_250_001 # human-readable label
manufact: Sonic Concepts          # informational label, e.g. 'Sonic Concepts' | 'Imasonic'
n_elem: 4                         # number of physical elements
model: CTX_250                    # serial of the device model entry
min_foc: 13.79                    # hardware steering range from the exit plane [mm]; optional on a unit,
max_foc: 61.48                    #   which then uses its model's. Give the source in a comment.
status: active                    # optional: active (default) | deprecated | excluded
superseded_by: ""                 # optional, with deprecated/excluded: serial to use instead
status_note: ''                   # optional, why
housing:                          # optional, display only (Nexus 3D view); a unit falls back to its model's
  diameter_mm: 70                 # housing diameter [mm]
  depth_mm: 20                    # from the exit plane (front face) to the back [mm]

transducer:
  annular:
    elem_n: 4                          # virtual element count used in simulation
    elem_id_mm: [0.0, 32.9, 46.1, 56.0]  # inner diameters of each ring [mm]
    elem_od_mm: [32.4, 45.6, 55.5, 64.0] # outer diameters of each ring [mm]
    curv_radius_mm: 63.2               # radius of curvature of the bowl [mm]
    dist_geom_ep_mm: 52.38             # distance from exit plane to geometric focus (centre of curvature) [mm]
  freq_hz: 250000.0                    # operating frequency [Hz]

combos:
  <DS>:                                # driving-system serial (ds/ds-<DS>_driver.yaml)
    ds:                                # optional; without it, taken from the driving-system entry
      serial: SC_203_035
      name: TPO junior SC_203_035
      manufact: Sonic Concepts
      available_ch: 4
    char_data_path: <TRAN>/tran-<TRAN>_ds-<DS>_axial.csv  # measured axial profiles, relative to local/profiles/
    drive_limit:                       # optional: largest drive of this driving system [Pa, as elem_p_amp]
      max_elem_p_amp_pa: 168000
      source: '...'                    # how the value was obtained
```

## Driving-system YAML (`ds/ds-<DS>_driver.yaml`)

```yaml
type: driving_system
serial: IGT_128_ch_comb_10_ch      # the key combos use
name: IGT 128 ch. - 10 ch.
manufact: IGT
available_ch: 128
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

1. Create `model-<MODEL>/tran-<SERIAL>_transducer.yaml` following the template above; for a unit of a model already in the library, give only `type`, `serial`, `model`, `name` and its combos.
2. Check that `elem_id_mm`, `elem_od_mm`, `curv_radius_mm`, and `dist_geom_ep_mm` are correct; steering phases depend on them.
3. Add an axial-profile CSV (+ sidecar) per driving-system combination under `local/profiles/<SERIAL>/`, and set `char_data_path` to it.
4. Register the combo in the calibration config (`parameters.calibration.combinations`).
