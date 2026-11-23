# Profile File Formats

![Pipeline steps; this page belongs to 2 Transducer](assets/figures/fig_stage-strip-transducer_diagram.svg){ .stage-strip }

The calibration fits read your hydrophone measurements from files. This page gives the layout of each file, how its focal depth is labelled, and how the files are named and placed, so the fits find them.

On this page:

- [Input File Formats](#input-file-formats): axial profile CSVs, plane and volume scans, depth labels, and folder naming.

!!! note "At a glance"
    - **When you need it:** You prepare hydrophone profiles for [Emulating Hydrophone Profiles](doc_calibration_hydrophone-profiles.md).
    - **Before you start:** Axial profiles (and optionally plane or volume scans) measured at the focal depths you use.

## Input File Formats

`refine_transducer_library` and `build_transducer_library` need one axial profile CSV per transducer–driving-system combo. Steering phases are computed from the ring geometry (`set_real_phases`), so manufacturer phase tables are not used. Plane and volume scans are optional: they constrain a [within-ring velocity profile](doc_calibration_hydrophone-profiles.md#within-ring-velocity-profile).

### Axial Profile CSV

```
,Focal depths (mm),,,,,,,
Distance WRT Exit plane (mm),40.0,50.0,60.0,70.0,80.0
0.5,0.12,0.31,0.04,0.01,0.0
1.0,0.18,0.44,0.08,0.02,0.0
...
```

| Row | Content |
|---|---|
| 1 (title) | Label strings only; skipped |
| 2 (depth) | First cell a label (ignored); then the focal depth of each column in mm from the exit plane |
| 3 onward (data) | Axial distance from the exit plane in mm; then intensity in W/cm² at each focal depth |

Rules:

- Separator: comma. Encoding: UTF-8 (a BOM is tolerated).
- Any number of focal-depth columns ≥ 1.
- Axial distance starts at or below the transducer face (0 mm is fine; negative values are ignored).
- Intensity is non-negative and always W/cm² (no units row). Trailing NaN columns are trimmed.
- Focal depths outside the transducer's `[min_foc, max_foc]` give a warning. `check_equipment_config(equip_path)` validates the equipment YAML before a rebuild.

#### Focal-Depth Labels

Column headers are nominal driver setpoints. They correspond to the [FLHM](doc_glossary.md#flhm) centre of the measured profile, at a manufacturer-dependent threshold ([focal depth convention](doc_calibration.md#focal-depth-convention-flhm)). On the Donders profiles (2026-09-30; median |header − FLHM centre|, first 10 mm skipped):

| Equipment | Convention | Median error |
|---|---|---|
| Sonic Concepts CTX + TPO | `flhm_3db` | 0.03–0.26 mm (−6 dB: 1.0–7.0 mm) |
| Imasonic + IGT 32 ch | `flhm_6db` | 0.04–0.08 mm (−3 dB: 2.1–2.6 mm) |
| Imasonic + IGT 128 ch | closest to `flhm_6db` | 0.60–0.75 mm (−3 dB: 1.2–1.7 mm) |

#### Depth Coverage Recommendations

Cover the full intended depth range and include the natural focal depth.

#### Profile Folder and File Naming

Measured profiles are site data in `local/profiles/` inside the PRESTUS root (ignored by [Git](https://git-scm.com/)), one subfolder per transducer. `build_transducer_library` and `calibration_setup` read them from there.

??? note "Folder layout, file names, and sidecar"
    ```
    local/profiles/
    ├── IS_PCD15473_01001/
    │   ├── tran-IS_PCD15473_01001_ds-IGT_32_ch_comb_10_ch_axial.csv
    │   ├── tran-IS_PCD15473_01001_ds-IGT_32_ch_comb_10_ch_axial.yaml         ← sidecar
    │   ├── tran-IS_PCD15473_01001_ds-IGT_32_ch_comb_10_ch_acq-OPM_axial.csv  ← repeat measurement
    │   ├── tran-IS_PCD15473_01001_ds-IGT_128_ch_comb_10_ch_axial.csv
    │   └── sourcedata/        ← never read: vendor spreadsheets, steering tables, deprecated files
    └── CTX_250_001/
        ├── tran-CTX_250_001_ds-SC_105_010_axial.csv
        └── tran-CTX_250_001_ds-SC_203_035_axial.csv
    ```

    File names follow `tran-{TRAN_SERIAL}_ds-{DS_SERIAL}[_acq-{LABEL}]_axial.csv`, using the serials of `config/equipment/` (underscores, not hyphens). `acq-` marks a repeat or alternative measurement of the same pair. One subfolder per transducer holds all its measurements whatever driving system recorded them.

    `char_data_path` in the equipment YAML is relative to this folder, for example `CTX_250_001/tran-CTX_250_001_ds-SC_105_010_axial.csv`. Point `calibration.path_input_axial` (for `calibration_setup`) or `CHAR_DATA_ROOT` (for `build_transducer_library`) at the folder.

    **Sidecar.** Each CSV has a YAML sidecar with the same stem, recording `quantity` and `units`, `distance_reference` (`exit_plane` | `bowl`), whether headers are `nominal` setpoints or `measured` positions (`focal_label_is`), their FLHM convention (`focal_label_reference`), acquisition conditions, and provenance. PRESTUS does not read sidecars: set `calibration.add_FDO` and `calibration.focal_reference` in the calibration config to match it.

### Plane and Volume Scans

Optional. Off-axis hydrophone scans for the ring-profile fit, read by `read_offaxis_scans` (`OffAxis` in `RingProfileOptions`; *Ring profile ▸ Options…* in the GUIs; `ring_profile_offaxis` in `config_calibration.yaml`).

**What to measure.** One transverse plane a few millimetres before the focus, at a setting that also has an axial profile, made the ring profile identifiable on CTX-500-026 (see [Within-ring velocity profile](doc_calibration_hydrophone-profiles.md#within-ring-velocity-profile)). Further planes or a volume add checks (e.g. of an off-centre or tilted axis) more than constraint.

Accepted files are a Sonic Concepts FUS Metrology export or a CSV with one row per point.

??? note "Accepted files, CSV columns, and reading call"
    **Accepted files.**

    - **Sonic Concepts FUS Metrology export**: `FUS_Metrology_Export_BableBrain.mat`, its folder, or a `.zip` of the export folder with `Planes/` and `Axial/`. Each entry of the cell array `data` needs `type`, `focus` (focal setting, mm from the exit plane), `x`, `y`, `z` (mm), and `ISPPA.data` with `ISPPA.unit = 'W/cm2'`. Entries of type `Axial` are used only to rescale the scans.
    - **CSV**, one row per point:

    | Column | Meaning |
    |---|---|
    | `focus_ep_mm` | Focal setting the scan was recorded at [mm from the exit plane] |
    | `x_mm`, `y_mm` | Lateral position [mm]; the beam axis at x = y = 0 (or use `Centre = 'centroid'`) |
    | `z_ep_mm` | Distance from the exit plane [mm] |
    | `intensity_wcm2` | Intensity [W/cm²] |
    | `label` (optional) | Name of the scan; separates scans at the same setting |

    ```
    focus_ep_mm,x_mm,y_mm,z_ep_mm,intensity_wcm2,label
    40,-15.0,-15.0,35,0.012,XY 35
    40,-14.5,-15.0,35,0.014,XY 35
    ...
    ```

    ```matlab
    [oa, scans] = read_offaxis_scans('SC CTX-500-026 (4 ch.).zip', struct('Axial', meas));
    % oa: OffAxis input of fit_ring_profile; scans: maps for plotting
    ```

??? note "How scans are read and reduced"
    - **Kind from the coordinates.** Points need not lie on a grid.

        | Varying coordinates | Kind | Reduced to |
        |---|---|---|
        | x, y (constant z) | Transverse plane | Azimuthally averaged radial profile |
        | x or y with z | Axial plane | (\|x\|, z) map, both sides folded together |
        | x, y, and z | 3D volume | (r, z) map, averaged around one axis parallel to z in each z bin |

    - **One scan per setting and label.** Several planes at the same setting need distinct `label` values in a CSV; without them, their points read as one volume.
    - **Absolute level.** A scan is compared with the axial profile at its focal setting (within 1 mm) relative to that profile's peak, so it must be in W/cm² at the same drive. When a Sonic Concepts export includes its own axial profiles, each scan is rescaled by the ratio of the fitted axial peak to the export's at that setting.
    - **Axis.** The scan axis is x = y = 0, or with `Centre = 'centroid'` the intensity centroid of the top 50 % of a transverse plane or volume. A volume's axis is assumed parallel to z; tilt is not corrected.
    - **Reduction options** (`OffAxisOptions`): `Kind` (`all`, `transverse`, `axial`, `volume`, or a list), `Weight` (per scan, relative to one axial profile; default 1), `RMaxMm` (30), `BinMm` (radial bin, 0.5), `ZBinMm` (z bin of axial planes and volumes, 1; the fit's field matrix grows with the number of reduced points), `Centre`.
