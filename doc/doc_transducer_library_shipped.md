# Shipped Library Entries

![Pipeline steps; this page belongs to 2 Transducer](assets/figures/fig_stage-strip-transducer_diagram.svg){ .stage-strip }

PRESTUS ships library entries for common devices and driving systems. This page lists them with how well each fits its hydrophone profiles, so you can judge how far to trust a simulation with an entry before using it.

On this page:

- [Shipped Entries and Fit Quality](#shipped-entries-and-fit-quality): the entries, their device models and units, drive limits, and deprecated entries.

!!! note "At a glance"
    - **When you need it:** To choose a device, or to judge how well its entry matches measurements.
    - **Before you start:** Nothing.

## Shipped Entries and Fit Quality

Agreement of the shipped entries (library model v3, `library_mode: per_depth`, built 2026-09-30) with the hydrophone axial profiles they were fitted to, from `meta.quality`. R² is on peak-normalised intensity beyond 10 mm from the [exit plane](doc_glossary.md#exit-plane), as a median over measured depths.

Columns: *geo*, geometric steering only; *in-sample*, the refinement at each measured depth (minimum over depths in brackets); *held-out*, each interior depth predicted from its two neighbours, as at run time between measured depths; *worst Δ*, the largest held-out loss against geometric steering at any depth.

| Combo | Manufacturer | Depths mm (n) | R² geo | R² in-sample (min) | R² held-out | worst Δ vs geo | Flags |
|---|---|---|---|---|---|---|---|
| IS_PCD15287_01001_IGT_128_ch_comb_10_ch | [Imasonic](https://www.imasonic.com) | 7–96 (13) | 0.990 | 0.998 (0.930) | 0.998 | +0.000 | n/a |
| IS_PCD15287_01001_IGT_32_ch_comb_10_ch | Imasonic | 10–90 (10) | 0.986 | 0.998 (0.988) | 0.998 | -0.002 | n/a |
| IS_PCD15287_01002_IGT_128_ch_comb_10_ch | Imasonic | 7–94 (13) | 0.982 | 0.998 (0.889) | 0.998 | +0.002 | n/a |
| IS_PCD15287_01002_IGT_32_ch_comb_10_ch | Imasonic | 10–90 (10) | 0.988 | 0.999 (0.955) | 0.999 | +0.001 | n/a |
| IS_PCD15473_01001_IGT_128_ch_comb_10_ch | Imasonic | 11–93 (11) | 0.981 | 0.996 (0.980) | 0.995 | -0.001 | n/a |
| IS_PCD15473_01001_IGT_32_ch_comb_10_ch | Imasonic | 10–90 (10) | 0.983 | 0.998 (0.983) | 0.998 | -0.005 | n/a |
| IS_PCD15473_01003_IGT_128_ch_comb_10_ch | Imasonic | n/a | n/a | n/a | n/a | n/a | geometric only (no profiles) |
| IS_PCD15473_01003_IGT_32_ch_comb_10_ch | Imasonic | 10–90 (10) | 0.982 | 0.997 (0.981) | 0.997 | -0.012 | n/a |
| CTX_250_001_SC_105_010 | [Sonic Concepts](https://sonicconcepts.com/) | 14–61 (12) | 0.786 | 0.984 (0.916) | 0.975 | -0.069 | n/a |
| CTX_250_001_SC_203_035 | Sonic Concepts | 14–61 (12) | 0.849 | 0.970 (0.904) | 0.938 | -0.277 | nearfield_excess |
| CTX_250_009_SC_105_010 | Sonic Concepts | 16–46 (12) | 0.838 | 0.917 (0.607) | 0.926 | -0.002 | nearfield_excess |
| CTX_250_009_SC_203_035 | Sonic Concepts | 16–38 (10) | 0.845 | 0.912 (0.243) | 0.913 | -0.004 | nearfield_excess |
| CTX_250_014_SC_105_010 | Sonic Concepts | 13–44 (12) | 0.869 | 0.906 (0.754) | 0.898 | -0.003 | n/a |
| CTX_250_014_SC_203_035 | Sonic Concepts | 13–44 (12) | 0.858 | 0.902 (0.751) | 0.892 | -0.001 | n/a |
| CTX_250_026_SC_105_010 | Sonic Concepts | 22–62 (12) | 0.852 | 0.975 (0.893) | 0.975 | +0.062 | nearfield_excess |
| CTX_250_026_SC_203_035 | Sonic Concepts | 22–61 (12) | 0.858 | 0.975 (0.911) | 0.976 | +0.056 | nearfield_excess |
| CTX_500_006_SC_105_010 | Sonic Concepts | 33–79 (10) | 0.910 | 0.921 (0.450) | 0.944 | -0.020 | n/a |
| CTX_500_006_SC_203_035 | Sonic Concepts | 36–78 (7) | 0.849 | 0.873 (0.477) | 0.924 | -0.002 | nearfield_excess |
| CTX_500_026_SC_105_010 | Sonic Concepts | 40–80 (12) | 0.909 | 0.973 (0.800) | 0.960 | +0.015 | n/a |
| CTX_500_026_SC_203_035 | Sonic Concepts | 40–79 (12) | 0.896 | 0.971 (0.832) | 0.959 | +0.013 | n/a |

`nearfield_excess`: the fitted model has a hotter near field than the measured profile at some depth.

Geometric steering alone reproduces the Imasonic profiles well (R² ≈ 0.98–0.99). For the Sonic Concepts CTX-500 / CTX-250 it is less accurate (R² ≈ 0.79–0.91), so refinement matters most there.

### Device Models and Units

Each device has a **model** entry, named after the model (`CTX_500`), and the library's **measured units** of it, named model plus unit id (`CTX_500_026`). A unit's calibration is fitted to that unit's hydrophone profiles, so use the unit's serial when your device is in the library. Otherwise name the model, which gives the model's geometry and a [geometric baseline](doc_transducer_library_maintenance.md#geometric-baseline-entry) calibration: accuracy is that of the *geo* column above. A serial of a unit that is not in the library (`CTX_500_031`) falls back to its model (`CTX_500`), with a warning; to use the unit's own calibration, [add it](doc_transducer_library_maintenance.md#adding-a-transducer-or-driving-system). A model entry has no driving system; with `combo.ds_serial` set, PRESTUS falls back to the model's generic entry and says so. Each entry names its model in `model:`; the naming rules are in `config/equipment/README.md`.

The model geometry is based on the [BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain) specification, and for Imasonic on the manufacturer's specification and drawings; where the library also has measured units, the unit geometry is identical.

| Model | Device | Rings | ROC mm | Exit plane → geometric focus mm | Steering range mm | Frequency kHz | Measured units |
|---|---|---|---|---|---|---|---|
| `CTX_500` | NeuroFUS CTX-500 | 4 | 62.94 | 52.38 | 33.6–82.5 | 500 | `CTX_500_026` |
| `CTX_250` | NeuroFUS CTX-250 | 4 | 62.94 | 52.38 | 25–60 | 250 | `CTX_250_001`, `_026` |
| `CTX_250_2ch` | NeuroFUS CTX-250, 2 ch. | 2 | 62.9 | 56.9 | 20–80 | 250 | `CTX_250_009`, `_014` |
| `DPX_500` | NeuroFUS DPX-500 | 4 | 150 | 144.9 | 50–120 | 500 | |
| `DPXPC_300` | NeuroFUS DPX-PC-300 | 4 | 150 | 144.9 | 50–120 | 300 | |
| `IS_PCD15287` | Imasonic PCD15287 | 10 | 75 | 65.3 | 7.05–100.3 | 300 | `IS_PCD15287_01001`, `_01002` |
| `IS_PCD15473` | Imasonic PCD15473 | 10 | 100 | 92.7 | 7.12–102.7 | 300 | `IS_PCD15473_01001`, `_01003` |
| `H246` | Sonic Concepts H-246 (flat) | 2 | flat | n/a | 25–95 | 500 | |
| `BSonix_35mm` … `_80mm` | BrainSonix BSonix | 1 | 43 / 68.5 / 74 / 87.5 | 32.68 / 55.66 / 64.92 / 80.51 | fixed | 650 | |

The steering range (`min_foc`–`max_foc`, from the exit plane) is the hardware range of the model; a unit without its own range uses its model's. The measured range is the span of depths with hydrophone profiles in the calibration entry (column *Depths mm* above). Within the steering range but outside the measured range, PRESTUS extrapolates the fitted model and warns (`load_transducer_from_library:extrapolated`); outside the steering range it warns that the device cannot be set there. For the Imasonic models, the provenance of each value is in the equipment YAML. `CTX_500_006` (2-channel CTX-500, 500 kHz) has no model entry. `H246` has flat rings and needs the kWaveArray source ([Flat Elements](doc_transducer.md#flat-elements)).

### Drive Limits

A driving system delivers a limited drive per channel, so at deep foci it may not reach the target [I<sub>SPPA</sub>](doc_glossary.md#isppa). With a limit, PRESTUS compares the library drive for `target_isppa_wcm2` with it and warns (`resolve_transducer_from_serial:driveLimit`), giving the largest free-water I<sub>SPPA</sub> reachable at that depth. The limit is `transducer.max_elem_p_amp` in the study config, else the combo's `drive_limit.max_elem_p_amp_pa` in the equipment entry. Without either, there is no check. The drive is compared as one value; per-ring gains are not included.

The shipped limits are estimates for the IGT 128-channel driver from the Radboud characterisation (full amplitude for a 1.00 MPa focus at about 83–89 mm): `IS_PCD15287_01001` and `_01002`, and `IS_PCD15473_01001`. Their source is in the equipment entry. Your amplifier may differ (for example the optional 10 W channels), so set `transducer.max_elem_p_amp` for your setup.

### Deprecated and Excluded Entries

A unit that is retired or replaced gets `status: deprecated`; an entry withheld until it is documented gets `status: excluded` (both optionally with `superseded_by` and `status_note`). Both stay in the library, so earlier studies re-run unchanged, and PRESTUS and the Nexus warn when a study names them. Excluded entries are not listed in the Nexus or `calibration_gui` and are skipped by `build_transducer_library`. The OPM variants of `IS_PCD15473_01001` and `_01003` are excluded until their setup is documented.
