# Adding and Maintaining Library Entries

![Pipeline steps; this page belongs to 2 Transducer](assets/figures/fig_stage-strip-transducer_diagram.svg){ .stage-strip }

A study can name a device with `transducer.serial` only if PRESTUS has its geometry in `config/equipment/` and a calibration entry in `config/calibration/`. This page covers adding a transducer or driving system, starting from a geometric baseline entry that uses geometric steering alone, and how maintainers rebuild, validate, and contribute shipped entries.

On this page:

- [Adding a Transducer or Driving System](#adding-a-transducer-or-driving-system): equipment entries by hand or from the GUI, and a geometric baseline entry.
- [Rebuilding a Shipped Entry](#rebuilding-a-shipped-entry): rebuilding from hydrophone profiles and validating the result.
- [Contributing an Entry](#contributing-an-entry): what a pull request for a new entry contains.

!!! note "At a glance"
    - **When you need it:** Your device has no entry in `config/equipment/` or `config/calibration/`.
    - **Before you start:** The device geometry (ring diameters, curvature, exit-plane distance, focal range); for a refined entry, hydrophone axial profiles.

How studies use the library is described in [Transducer Library](doc_transducer_library.md); fitting an entry to hydrophone profiles in [Calibration](doc_calibration.md).

## Adding a Transducer or Driving System

### New Driving System for an Existing Transducer

1. Add a combo entry to the transducer's entry in `config/equipment/` (`model-<MODEL>/tran-<TRAN>_transducer.yaml`), and describe a new driving system in `config/equipment/ds/ds-<DS>_driver.yaml`:

    ```yaml
    combos:
      MY_NEW_DS:
        ds:
          serial: MY_NEW_DS
          name: My Driving System
          manufact: MyCompany
        # char_data_path: add when hydrophone data is available
    ```

2. Generate a [geometric baseline](#geometric-baseline-entry).
3. Use it in study configs with `combo.ds_serial: MY_NEW_DS`.
4. [Refine it](doc_calibration_hydrophone-profiles.md#refining-a-library-entry) when hydrophone data is available.

### New Transducer by Hand

1. Create `config/equipment/model-<MODEL>/tran-<MY_SERIAL>_transducer.yaml`, copying an existing file as a template (a flat `config/equipment/<MY_SERIAL>.yaml` is read as well). A measured unit of a model already in the library needs only `type`, `serial`, `model`, `name`, and its combos; it takes the geometry and steering range from the model entry. Required fields for a device model of an annular array:

    ```yaml
    type: transducer
    serial: MY_SERIAL
    name: Human-readable name
    manufact: Imasonic    # or: Sonic Concepts
    min_foc: 30.0         # hardware steering range from the exit plane [mm]; a unit may omit it
    max_foc: 90.0         # and use its model's
    transducer:
      annular:
        elem_n: 10
        elem_id_mm: [0.0, ...]
        elem_od_mm: [21.0, ...]
        curv_radius_mm: 100.0
        dist_geom_ep_mm: 92.7
      freq_hz: 300000.0
    ```

2. Run `check_equipment_config(equip_path)` to validate the geometry and its physical plausibility.
3. Continue with the [geometric baseline](#geometric-baseline-entry).

### New Transducer from the GUI

The Nexus [Equipment sub-tab](#new-transducer-or-driving-system-in-the-nexus) does the same. Click **New…** under the transducer list in `calibration_gui`. The form takes the type (annular or matrix), frequency, curvature radius, exit-plane-to-geometric-focus distance (`dist_geom_ep_mm`) and focal range, then the type-specific geometry:

- **Annular:** ring inner and outer diameters and the number of physical channels, plus optionally one driving system with its profile CSV.
- **Matrix:** steering (1D/3D), curved or flat, aperture diameter, element shape and size, and a `define_here` grid (`rect` rows × columns with spacing, `fibonacci` count and kerf, or `fermat` count). Coordinate-file layouts (`extract_from_file`) are not offered in the form.

**Copy from** fills the form from an existing transducer. A side view and a face view update as you type: rings coloured by channel, or matrix element footprints, with the [exit plane](doc_glossary.md#exit-plane), geometric focus, and steering range. The message box reports the ring/channel count or the element count and median pitch, and flags overlapping rings or elements. On save the GUI writes the entry into your own library (the output folder, `OutputPath`): into the entry's own file when it already has one, else `tran-<SERIAL>_transducer.yaml` in that folder, and selects the new transducer. It does not change `config/equipment/`; only maintainers (`AllowShipped`) save there (into the existing file, else `model-<MODEL>/tran-<SERIAL>_transducer.yaml` when the model folder exists, else `<SERIAL>.yaml`). For an annular transducer, the geometric prediction is shown right away. Studies resolve the new serial when the folder is one of their libraries (`transducer.library_path`, `path.transducer_library`, or `PRESTUS_TRANSDUCER_LIBRARY`); own folders are searched for `*_transducer.yaml` and `*_driver.yaml` files only.

??? note "What the form checks and stores"
    On save the GUI:

    - rejects a serial that is not a valid identifier, unequal inner/outer diameter counts, a ring whose outer diameter does not exceed its inner diameter, a ring count that is not a multiple of the physical channel count (`fit_transducer_model` groups consecutive rings into channels), a curved matrix whose curvature radius is not larger than half the aperture, and a matrix grid with no element inside the aperture;
    - runs `check_equipment_config` on the entry and asks for confirmation if it reports issues.

    The profile CSV path is stored relative to `ProfileRoot` when it lies under it, otherwise as an absolute path. For matrix entries, the top-level `n_elem` is the element count inside the aperture.

    **View geometry** draws the selected transducer the same way. Matrix transducers are geometry only in the calibration GUI: the axial forward model behind prediction and fitting (`axial_transfer_matrix`, `steering_path_lengths`) is defined for annular arrays.

### New Transducer or Driving System in the Nexus

The Nexus has the same form: *Transducer* tab, **Advanced**, **Equipment** sub-tab. It also adds driving-system combinations to a transducer, including driving systems not yet in the library, and shows any library transducer (*View geometry* on the Calibration page opens it there).

- **Start from** shows a library transducer (side view and face view) and fills the form from it: its geometry under a new serial, or, with *Edit this entry*, the entry itself with its serial, other fields (housing, status, corrections), and combinations, to add driving systems or change it.
- **Type** is annular (ring diameters and physical channels) or matrix (steering, curved or flat, aperture, element shape and size, and a `rect`, `fibonacci`, or `fermat` grid). MATLAB computes the matrix element layout with `convert_to_element_pos`, as the simulation places the elements; the face view draws the footprints and the message box gives the element count and median pitch.
- **Driving systems** lists the combinations, each a driving-system serial and an optional axial-profile CSV (stored relative to `local/profiles/` when it lies there). A serial that is not in the library is a new driving system: its name, manufacturer, and channel count are asked for, and it is written as `ds-<DS>_driver.yaml` next to the transducer.
- **Save** writes `tran-<SERIAL>_transducer.yaml` (or rewrites the entry's own file) into the chosen own library folder, after `check_equipment_config` in MATLAB, whose issues ask for confirmation. An own entry with a shipped serial replaces the shipped one. Only with *Shipped library (maintainers)* does it write into `config/equipment/`, in the shipped layout (as `calibration_gui`'s `AllowShipped`). Studies and the Calibration sub-tab see the entry when the folder is one of their libraries (`transducer.library_path`, `path.transducer_library`, or `PRESTUS_TRANSDUCER_LIBRARY`).

The form applies the same checks as `calibration_gui` (below), plus a driving system with fewer channels than the transducer and a driving system listed twice. Profile CSVs belong to annular transducers only.

### Geometric Baseline Entry

For a transducer or driving system that is not in the shipped library and has no hydrophone data:

```matlab
COMBOS_TO_RUN = {'MY_SERIAL_MY_DS_SERIAL'};              % a combo without char_data_path
LIBRARY_PATH  = '/project/my_study/transducer_library';  % omit to write config/calibration/
run(fullfile(prestus_root, 'functions', 'calibration', 'build_transducer_library.m'))
```

This writes an entry with an empty offset (pure geometric steering). Accuracy equals the geometric baseline in the [shipped combos](doc_transducer_library_shipped.md#shipped-entries-and-fit-quality) table. Add the device to `config/equipment/` first (above).

## Rebuilding a Shipped Entry

Maintainers only. Set `char_data_path` in the transducer's equipment entry and run `build_transducer_library`. `CHAR_DATA_ROOT` defaults to `local/profiles/` (see [Input file formats](doc_calibration_profile-formats.md#input-file-formats)), and `LIBRARY_PATH` to `config/calibration/`. A missing profile CSV is an error; it does not fall back silently to a geometry-only entry. Then check the entry:

```matlab
validate_transducer_library(library_path, equip_path)   % must report no 'error'
```

??? note "What the build fits and records"
    Every measured depth becomes a refinement (`PER_DEPTH_FIT_OPTIONS`, as in [Calibration](doc_calibration_hydrophone-profiles.md#refining-a-library-entry)), on a geometric global model with the measured drive per depth. `meta.quality` records:

    - the in-sample R² (median, min) and the geometric R²;
    - the **held-out interpolation** R²: each interior depth predicted from its two neighbours, without axial shift, as at run time.

    Held-out interpolation shows how much of the refinement carries over between measured depths. On CTX-250-001 + TPO-203-035 it is 0.938, against 0.970 in-sample and 0.849 for geometric steering. Space the measured depths closely where you will stimulate.

## Contributing an Entry

Open a PR with:

- the rebuilt YAML (with `global_model.forward_model` set)
- `validate_transducer_library` passing with no errors
- the `meta.quality` summary (in-sample, held-out interpolation, and geometric R²)
- characterisation conditions in the combo YAML `notes` field (water temperature, hydrophone model, tank geometry)
