# Adding and maintaining library entries

This page is for adding a transducer or driving system that is not in the shipped library, and for maintainers who rebuild or contribute shipped entries. How the library works, how studies use it and how to refine an entry to your own hydrophone profiles is described in [Transducer Library](doc_transducer_library.md).

## Rebuild a shipped entry

Maintainers only. Set `char_data_path` in `config/equipment/{serial}.yaml` and run `build_transducer_library`. `CHAR_DATA_ROOT` defaults to `local/profiles/` (see [Input file formats](doc_transducer_library.md#appendix-input-file-formats)), and `LIBRARY_PATH` to `config/calibration/`.

Every measured depth becomes a refinement (`PER_DEPTH_FIT_OPTIONS`, as in [Level 2](doc_transducer_library.md#level-2-refine-to-your-own-hydrophone-profiles-optional)), on a geometric global model with the measured drive per depth. `meta.quality` records:
- the in-sample R² (median, min) and the geometric R²;
- the **held-out interpolation** R²: each interior depth predicted from its two neighbours, without axial shift, as at run time.

Held-out interpolation shows how much of the refinement carries over between measured depths. On CTX-250-001 + TPO-203-035 it is 0.938, against 0.970 in-sample and 0.849 for geometric steering. Space the measured depths closely where you will stimulate.

A missing profile CSV is an error; it does not fall back silently to a geometry-only entry. Then check the entry:

```matlab
validate_transducer_library(library_path, equip_path)   % must report no 'error'
```

## Contribute an entry back

Open a PR with:
- the rebuilt YAML (with `global_model.forward_model` set)
- `validate_transducer_library` passing with no errors
- the `meta.quality` summary (in-sample, held-out interpolation and geometric R²)
- characterisation conditions in the combo YAML `notes` field (water temperature, hydrophone model, tank geometry)

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
       # char_data_path: add when hydrophone data is available
   ```
2. Generate a geometric baseline ([Level 1](doc_transducer_library.md#level-1-geometry-only-baseline-no-hydrophone-data)).
3. Use it in study configs with `combo.ds_serial: MY_NEW_DS`.
4. Refine ([Level 2](doc_transducer_library.md#level-2-refine-to-your-own-hydrophone-profiles-optional)) when hydrophone data is available.

### New transducer (not in equipment config)

1. Create `config/equipment/{MY_SERIAL}.yaml`, copying an existing file as a template.

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
2. Run `check_equipment_config(equip_path)` to validate the geometry and its physical plausibility.
3. Continue from [Level 1](doc_transducer_library.md#level-1-geometry-only-baseline-no-hydrophone-data).

Alternatively, click **New…** under the transducer list in `calibration_gui`. The form takes the type (annular or matrix), frequency, curvature radius, exit-plane-to-geometric-focus distance (`dist_geom_ep_mm`) and focal range, then the type-specific geometry:

- **Annular:** ring inner and outer diameters and the number of physical channels, plus optionally one driving system with its profile CSV.
- **Matrix:** steering (1D/3D), curved or flat, aperture diameter, element shape and size, and a `define_here` grid (`rect` rows × columns with spacing, `fibonacci` count and kerf, or `fermat` count). Coordinate-file layouts (`extract_from_file`) are not offered in the form.

**Copy from** fills the form from an existing transducer. A side view and a face view update as you type: rings coloured by channel, or matrix element footprints, with the exit plane, geometric focus and steering range. The message box reports the ring/channel count or the element count and median pitch, and flags overlapping rings or elements. On save the GUI:

- rejects a serial that is not a valid identifier, unequal inner/outer diameter counts, a ring whose outer diameter does not exceed its inner diameter, a ring count that is not a multiple of the physical channel count (`fit_transducer_model` groups consecutive rings into channels), a curved matrix whose curvature radius is not larger than half the aperture, and a matrix grid with no element inside the aperture;
- runs `check_equipment_config` on the entry and asks for confirmation if it reports issues;
- writes `{serial}.yaml` to the GUI's equipment folder (`EquipPath`, default `config/equipment/`) and selects the new transducer. For an annular transducer, the geometric prediction is shown right away.

The profile CSV path is stored relative to `ProfileRoot` when it lies under it, otherwise as an absolute path. For matrix entries, the top-level `n_elem` is the element count inside the aperture. Serial-based study configs (`transducer.serial`) always read `config/equipment/`, so save there if the pipeline should resolve the new serial.

**View geometry** draws the selected transducer the same way. Matrix transducers are geometry only in the calibration GUI: the axial forward model behind prediction and fitting (`axial_transfer_matrix`, `steering_path_lengths`) is defined for annular arrays.
