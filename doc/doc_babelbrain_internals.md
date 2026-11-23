# BabelBrain Backend: Integration Internals

PRESTUS runs [BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain) through a bridge: `babelbrain_wrapper.m` writes the inputs and calls `run_babelbrain_case.py`, which patches BabelBrain at runtime and runs it, and the readers map the results back onto the PRESTUS grid. Most problems with this backend are misalignment between the two grids. This page documents each step of the bridge, for debugging a run (start with [Output Alignment](#output-alignment) and the [zero-pressure checklist](#zero-pressure-debugging-checklist)) and for changing the bridge; for using the backend see [BabelBrain backend](doc_babelbrain.md).

On this page:

- [Call Chain](#call-chain): the functions from the pipeline to BabelBrain, and the runtime patches the bridge applies.
- [GPU Backend Selection and Platform Support](#gpu-backend-selection-and-platform-support): `computing_backend` and `device_name` per platform.
- [Grid Resolution](#grid-resolution): why the input mask is resampled to BabelBrain's step.
- [Placement Correction and Drive Phases](#placement-correction-and-drive-phases): the placement correction, steering and phases, the free-water check, and the bowl rim.
- [Calibration Weights File](#calibration-weights-file): what the weights store and the [HDF5](https://www.hdfgroup.org/solutions/hdf5/) layout.
- [Free-Water Amplitude Calibration](#free-water-amplitude-calibration): the correction of BabelBrain's source bias and the water baseline checks.
- [Inputs, Geometry, and Modes](#inputs-geometry-and-modes): CT, air, and water inputs, the single-transducer limit, and mode-specific behaviour.
- [Output Alignment](#output-alignment): how the field lands on the PRESTUS masks.
- [BabelBrain Output Format](#babelbrain-output-format): the datasets read, axis flips, interpolation, and the [NIfTI](https://nifti.nimh.nih.gov/) affine.
- [Zero-Pressure Debugging Checklist](#zero-pressure-debugging-checklist): a checklist for missing pressure.
- [Bone Attenuation and Heat Source](#bone-attenuation-and-heat-source): BabelBrain's Q-factor attenuation model, its bone Q-correction, and why bone heats from the uncorrected attenuation.
- [Built-In Material Table](#built-in-material-table): BabelBrain's own tissue values, used only if material injection fails.
- [Tissue Labels and Material Properties](#tissue-labels-and-material-properties): label mapping, CT-based skull properties, and internal air.

!!! note "At a glance"
    - **Tested with:** BabelBrain 0.8.8, the `external/BabelBrain` submodule pinned at the v0.8.8 tag (`7aeccd3`). Other versions may break the runtime patches below; each run records the version it used ([version tracking](doc_babelbrain.md)).

## Call Chain

```
prestus_pipeline.m
  └─ acoustic_stage.m            ← dispatches on parameters.acoustic_backend (run_acoustic_backend)
       └─ babelbrain_wrapper.m   ← MATLAB orchestrator
            ├─ babelbrain_tx_geometry.m      (PRESTUS rings or BabelBrain preset)
            ├─ write_babelbrain_input.m       (NIfTI label mask)
            ├─ export_babelbrain_weights.m    (calibration HDF5)
            └─ system() → run_babelbrain_case.py   ← Python bridge
                 └─ BabelBrain RUN_SIM().RunCases()
```

The bridge runs in `babelbrain.python_env`, writes `*DataForSim.h5` to `<dir_cache>/<subject_id>/` and patches BabelBrain at runtime (never `external/`). Without the skull-attenuation patch, injected bone attenuation is not realised (QCorr stays 1.0), although the material table looks correct.

**Runtime patches**

| Patch | Effect | Disable (bridge flag, not a config key) |
|---|---|---|
| `_inject_prestus_materials` | Overwrites `MatFreq` from `--medium_properties_json` (attenuation converted from dB/(cm·MHz^y) to Np/m) | runs only without `--no_skull_atten_fix` |
| `_patch_noncat_skull_attenuation` | Bone shear = 0 (fluid bone), so the specified skull attenuation is realised | `--no_skull_atten_fix` |
| `_patch_zdomain_reference` | `CalculateDomainZReference` returns 0 (see [Output alignment](#output-alignment)) | `--no_zdomain_fix` |
| `_patch_annular_surface` | Ring tessellation samples exactly the configured bands (`np.arange` endpoint overshoot added a band to some rings) | always on |
| `_patch_annular_array_weights` | Expands per-ring calibration weights to sub-elements (focused annular arrays) | with a weights file |
| `_patch_resave_normalized`, `_patch_calculate_temperature_effects` | Crash guards (negative mask index; tissue class absent from the domain) | always on |

## GPU Backend Selection and Platform Support

| `computing_backend` | Backend | Hardware | Notes |
|---|---|---|---|
| `1` | [CUDA](https://developer.nvidia.com/cuda-toolkit) | NVIDIA GPUs | Linux; request a GPU in [SLURM](https://slurm.schedmd.com/) (`--gres=gpu:1`) |
| `2` | [OpenCL](https://www.khronos.org/opencl/) | AMD GPUs, Intel GPUs, some NVIDIA | Needs [`pyopencl`](https://documen.tician.de/pyopencl/) in the BabelBrain environment |
| `3` | [Metal](https://developer.apple.com/metal/) | Apple Silicon (M1/M2/M3/M4); macOS Intel with `environment_mac_x64-314.yml` (slow) | Only GPU acoustic solver on Apple Silicon |
| `4` | CPU | Any | Fallback; much slower |

`device_name: auto` picks the first GPU (`nvidia-smi`, `pyopencl`) or `'Apple'`. Set it explicitly to choose among several GPUs. Windows is untested (Unix-style path separators in `system()`).

## Grid Resolution

BabelBrain's step is `SmallestSOS / (PPW × freq)` (`SmallestSOS` normally ≈1102.5 m/s, `DensityToSSoSPichardo`), but it takes the input array shape as its grid. `write_babelbrain_input.m` therefore resamples the mask (nearest-neighbour) to that step (`bb_native_res_mm`, logged as "BabelBrain native resolution").

## Placement Correction and Drive Phases

BabelBrain puts the apex at `z_mm ≈ 0` and the FocalSpot at preset `FocalLength` (`Babel_{tx_system}/default.yaml`); PRESTUS has the apex `|focus_pos − trans_pos|` from the target. `babelbrain_placement_correction.m` (shared with `export_babelbrain_inputs`) moves the transducer by the difference:

```
TxMechanicalAdjustmentZ = preset_FocalLength − |focus_pos − trans_pos| · res_mm   (+ any user standoff)
```

The apex-to-target distance differs from the preset (e.g. 71 mm vs 92 mm for R15473 at 85 mm EP with a 64 mm-deep target). Uncorrected, the field shifts by `|focus_pos − trans_pos| − FocalLength` along the beam (−7.7 mm for R15473 with the target at 85 mm EP; −28.75 mm at 64 mm EP).

The standoff is `babelbrain.tx_mech_adjust_z_m`. Without `trans_pos`, the fallback is `preset_FocalLength − focal_distance_bowl`. BabelBrain aims at `center_z = ZDim[FocalSpot] + ZSteering + TxMechanicalAdjustmentZ`, so transducer and aim move together.

| Drive | `ZSteering` | Behaviour |
|---|---|---|
| Explicit phases (`tr.annular.elem_phase_rad`) | 0 | Weights carry the full drive as a delta from the flat natural-focus reference; BabelBrain reproduces the [k-Wave](https://github.com/ucl-bug/k-wave)/[jwave](https://github.com/ucl-bug/jwave)/Rayleigh ring phases. `focal_calibration` only reports the predicted focus. |
| No explicit phases | `aim_bowl − preset_FocalLength` | `aim_bowl` is the steer `D*` from `calibrate_focal_steering.m` (lands `flhm_center`/`peak` on `focal_distance_ep`) with `focal_calibration: true`, else `focal_distance_bowl`. |
| `babelbrain.drive: native` | `focal_distance_bowl − preset_FocalLength` | BabelBrain computes its own ring phases for that aim, as a run from its interface with the same focal setting. No weights file and no focal calibration. Its free-water focus can differ from the PRESTUS drive's by a few millimetres. |
| `babelbrain.z_steering_m` set | pinned | Overrides both; the weight-delta reference follows `FocalLength + ZSteering`. |

`Single` and `BSonix` presets lack `TxDiam`/`FocalLength`, so the wrapper passes `--aperture` (max `elem_od_mm`) and `--focal_length` (`curv_radius_mm`).

**Free-water check and bowl rim**

**Free-water check.** With explicit phases, the BabelBrain axial profile must overlay the Rayleigh profile of the same phases. A rigid offset of near-field lobes and peak means the placement correction is missing. A lobe *shape* mismatch means the drive phases differ (e.g. weights delta referenced to the wrong aim).

**Bowl rim.** A rim inside the domain makes BabelBrain snap the transducer back in whole voxels (log `Warning: Need to apply correction to reposition Tx for …`; non-zero `source_residual_mm`). `write_babelbrain_input.m` therefore trims `ceil((d_apex + bowl_depth)/dx) + 1` transducer-side slices (bowl depth 5.4 mm for R15473; `d_apex` negative if the apex is outside the grid). It warns if that slab contains tissue; it is zero-filled on paste-back.


## Calibration Weights File

The file stores deltas from BabelBrain's ring phases at its aim `FocalLength + ZSteering`:

```
w[n] = exp(i × (elem_phase[n] − geometric_phase_at_aim[n]))
```

**File layout and import/export**

- `/CALIBRATION` is a 1-D complex64 dataset ([h5py](https://www.h5py.org/) compound `{r, i}`) as in BabelBrain's `SaveToH5py`; every root item has the `type` attribute `ReadFromH5py` requires. `convert_weights_if_needed` converts the split `/CALIBRATION/real` + `/imag` layout.
- `export_babelbrain_weights` writes one weight per ring; the bridge expands them to sub-elements.
- `import_babelbrain_weights` reads `TxCalibration` files (one weight per sub-element) given the per-ring element counts (6th argument).

<a id="free-water-calibration-transfer"></a>
## Free-Water Amplitude Calibration

`water_amp_calibration` scales the head simulation by the factor that brings BabelBrain's free-water [I<sub>SPPA</sub>](doc_glossary.md#isppa) to `transducer.target_isppa_wcm2`, which also removes the few-% pressure bias of BabelBrain's source correction. The thermal input is not rescaled, because BabelBrain renormalises it to `Isppa`.

**Why the factor transfers, and water baseline checks**

BabelBrain's `Correction = 100/(100 − ExpectedError(CFL))`, ×√2 (`babel_integration_base`) leaves a few % pressure bias at CFL ~0.5. The factor transfers to tissue because the simulation is linear and the bias sits in the shared source chain. The transcranial timestep differs (sub-602 pFOV, R15473 300 kHz 6 PPW: CFL/water 0.300 → Correction 0.290 vs 0.471 → 0.465). The correction tracks it within ~2% in pressure (0.290/0.465 = 0.624 vs CFL ratio 0.637), ~4% in intensity.

The thermal renormalisation is `PressureRatio = sqrt(Isppa·2ρc)/max(p_amp[brain])`.

**Water baseline checks** (`babelbrain_wrapper.m`):

- The `--water_only` run must get the same `--source_amp_pa` as the transcranial run. Otherwise BabelBrain drives at ~50 kPa and `freefield_isppa_wcm2` is ~(156.6/50)^2 ≈ 10x too low for R15473.
- `Water_DataForSim.h5` in `<dir_cache>/<subject_id>/` is reused if present, unless `io.overwrite_files` is `'always'`/`'force'`. With `cleanup_cache: true` it is moved to `<dir_cache>/` after the run and not reused. If the I<sub>SPPA</sub>-normalised field looks wrong, compare its mtime with the main `DataForSim.h5`.

## Inputs, Geometry, and Modes

PRESTUS does not rotate its grid to the trajectory, so the beam must lie on z and one transducer is supported.

**Water, CT, and air inputs**

- **CT.** `write_babelbrain_input.m` writes raw-HU `CT.nii.gz` on the `BabelViscoInput.nii.gz` grid. The bridge quantises it (uniform, 10-bit, BabelBrain's default) into the index map BabelBrain expects and writes `CT-cal.npz`. The wrapper passes `--use_ct`, `--ct_mapping`, and optionally `--use_petra`.
- **Air.** BabelBrain does not update the fields inside the air mask (`ReflectorMask`), so air needs no material and no tiny timestep. `bExtractAirRegions` is on only with `air_regions`, because PRESTUS otherwise supplies no air mask.
- **Uncertainty.** Variants go in via `--medium_properties_json`; on the CT path only skin and brain apply.
- **Water.** Simulated water comes from the module-level `Material` (the 1120 kHz `MatFreq` entry), so the injected `Water` line affects only `GetSmallestSOS`. Water is not injected, because it also sets the source factor `FactorConvPtoU`; the impedance difference is ~0.6%.

**Geometry rationale and mode-specific behaviour**

BabelBrain's own workflow rotates the grid instead (`BabelDatasetPreps.GetSkullMaskFromSimbNIBSSTL`). Multi-transducer support would need beam-aligned per-transducer domains (as `grid.async_layout: per_transducer` for k-Wave) and a thermal step that accepts a combined heat source.

- **`multi_isppa`.** The thermal result is named after the shared acoustic file (`<stem>-ThermalField_AllCombinations.h5`), so thermal jobs are serialised on HPC.
- **Post-hoc water.** The wrapper passes `--water_only` under the baseline's target name, which enables cache reuse.
- **Partial FOV.** The wrapper disables `bTightNarrowBeamDomain` on a cropped mask, because its lateral footprint under-shoots the beam near the skull.

<a id="output-alignment-tissue-marker-anchoring-and-the-zdomainstart-fix"></a>
## Output Alignment

Field and `MaterialMap` must land exactly on `medium_masks`.

- **`ZDomainStart` = 0.** BabelBrain adds a skin-surface heuristic (`GlobalStartSkin − LineOfSightStartSkin`) to the source position in `CalculateRayleighFieldsForward`. PRESTUS already sets `TxMechanicalAdjustmentZ`, so the patch removes it. The log line "ZDomainStart neutralized: ... mm" shows the removed amount.
- **Marker anchor.** The readers anchor the resample on `/TargetLocation`, the label-5 marker at `focus_pos` (benchmark: label agreement 83% → 99%, bone 20% → 98%). Preset `FocalLength` (`focal_length_mm`) is the fallback for files without it.

`<target>_prestus_anchor.json` records realised-minus-requested source position. Expect ≤ 1 native voxel; `read_babelbrain_output.m` warns above 1.5, which points to a missing `ZDomainStart` patch.

## BabelBrain Output Format

`sensor_data.p_max_all` is `/p_amp` (peak pressure [Pa], `[Nx, Ny, Nz]`), or `/p_amp_refocus` with `do_refocusing: true`. Axes are permuted to `[lateral X, lateral Y, axial]` (dim 3 = beam).

**Axis flips and interpolation**

| Item | Handling (change both sides together) |
|---|---|
| Acoustic arrays | BabelBrain flips them (`np.flip(axis=2)`) before saving, but not `z_vec`; the reader uses `zdim = z_vec_bb(end:-1:1)`. |
| `/TargetLocation` | Computed before the flip; index the raw `/z_vec`. |
| Input mask | `write_babelbrain_input.m` pre-flips when the beam runs along +z (`focus_pos(3) >= trans_pos(3)`), so BabelBrain's flip on load cancels it; with the beam along −z, BabelBrain's flip alone orients it. The reader resamples with the same `beam_sign`. |
| Thermal arrays | Saved in `z_vec` order; not reversed. Resampled onto `parameters_th` from `thermal_stage`, anchored on the thermal focus (`grid.thermal_focus_pos_exact` etc. from `thermal_grid_setup`, else `transducer(1).focus_pos`). Rises over `BaselineTemperature`, so outside voxels read baseline. |
| Interpolation | Fields `'linear'`; `MaterialMap` `method='mode'`, because `'nearest'` drops sub-voxel cortical bone. |

The NIfTI affine comes from `grid_to_ras_affine.m`: `planimg.inv_transf`, `planimg.t1_header`, `planimg.fov_offset_ac`, then the dim-3 flip, in that order. A diagonal-only affine misaligns whenever `grid.mode` rotates the head (the default).

## Zero-Pressure Debugging Checklist

Inspect `/MaterialMap`, `/p_amp`, `/TargetLocation`, and `/z_vec` in `*DataForSim.h5`, then:

| Check | Points to |
|---|---|
| Console "zero-filled (no BabelBrain data)" line non-zero | Domain too short; voxels outside it are zero-filled |
| `ZShrink_L` > 0, `XShrink_*`/`YShrink_*`, or fewer tissue voxels in `MaterialMap` than in `BabelViscoInput.nii.gz` | Narrow-beam mode cropped tissue |
| "BabelBrain native resolution" vs `grid.resolution_mm` | Mask not resampled to BabelBrain's step |
| No `/TargetLocation` in the file | Fallback anchor; field may be shifted along the beam |
| Skin→trabecular in a central-line label dump | Raw segmentation error |

`simulation.debug: 1` writes axis plots (`sub-<id>_babelbrain_fov_debug<output_affix>.png`/`.csv`) and `sub-<id>_babelbrain_material_map<output_affix>.nii.gz` to `<dir_output>/debug/babelbrain/`.

<a id="built-in-material-table"></a>
## Bone Attenuation and Heat Source

BabelBrain models attenuation as a quality factor, Q = ω / (2cα) ([Pichardo 2023](https://doi.org/10.1109/TUFFC.2023.3274046), eq. 3), and fits one relaxation mechanism per material to it (BabelViscoFDTD `OptimalTauForQFactor`). The correction factor q<sub>c</sub> (`QCorrection`) multiplies Q, so a material propagates with α / q<sub>c</sub>: BabelBrain's own fit gives realised/specified attenuation 1/q<sub>c</sub> (0.333 for every bone material of a CT-path run at 500 kHz with q<sub>c</sub> = 3). On the CT path q<sub>c</sub> applies to the bone materials only (`QCorrArr[6:]` in `babel_integration_base.py`); soft tissue keeps 1. The paper motivates q<sub>c</sub> by HU-to-attenuation maps that measure total loss, including scattering that a voxel-level solver produces again; Webb-Marsac with q<sub>c</sub> = 3 reproduced measured transcranial transmission and is BabelBrain's CT default. That calibration was made for BabelBrain's solver at 6 PPW, so PRESTUS applies q<sub>c</sub> only on the BabelBrain backend. A bone-slab phantom (6 mm cortical, 78.5 Np/m at 300 kHz, PPW 6) gave realised/specified attenuation 0.91 at Q = 1 (k-Wave 0.91), 0.29 at 3, 0.14 at 6, 0.087 at 10, and 0.063 at 14. BabelBrain's CT default of 3.0 gave transcranial I<sub>SPPA</sub> 10.6 vs 7.6 W/cm² for k-Wave (sub-602 pFOV), so the non-CT path uses Q = 1. `babelbrain.q_correction` > 0 reaches `BabelFTD_Simulations_BASE.__init__` via `--q_correction` → `BABEL_PYTEST_QFACTOR`. MATLAB `setenv` does not work here, because `run_command` uses Java's ProcessBuilder with the JVM start-up environment.

BabelBrain's heat source uses the specified attenuation, not α / q<sub>c</sub>: per voxel it deposits p² / (2ρc) × absorption fraction × (1 − e<sup>−2hα</sup>) / h, with h the grid step (`getQCoeff` in BabelViscoFDTD `RayleighAndBHTE.py`). For small hα this is 2α × fraction; at the highest CT-path bone attenuation (≈650 Np/m, h = 0.37 mm) it is about 20 % lower. On the CT path the bone wave therefore loses α / q<sub>c</sub> while the bone heats from α × 0.155. To reproduce this with k-Wave, set bone attenuation to α / q<sub>c</sub> and the bone absorption fraction to q<sub>c</sub> × 0.155, since PRESTUS heats from the propagation attenuation × absorption fraction.

## Built-In Material Table

Literature fits tabulated every 5 kHz from 100 to 1125 kHz (`MatFreq` in `TranscranialModeling/babel_integration/babel_integration_helpers.py`). They apply only if material injection fails.

**Built-in values and references**

| Tissue | Density [kg/m³] | Long. SoS [m/s] 300 / 500 / 650 kHz | Long. att. [Np/m] 300 / 500 / 650 kHz |
|---|---|---|---|
| Water | 1000 | 1500 / 1500 / 1500 | 0 / 0 / 0 |
| Skin | 1116 | 1537 / 1537 / 1537 | 1.4 / 2.3 / 3.0 |
| Brain | 1041 | 1562 / 1562 / 1562 | 2.1 / 3.5 / 4.5 |
| White matter | 1041 | 1537 / 1537 / 1537 | 3.1 / 5.1 / 6.6 |
| Gray matter | 1045 | 1520 / 1520 / 1520 | 1.3 / 2.2 / 2.9 |
| CSF | 1007 | 1507 / 1507 / 1507 | 0.03 / 0.05 / 0.06 |
| Cortical bone | 1897 | 2452 / 2476 / 2494 | 49 / 81 / 106 |
| Trabecular bone | 1738 | 2148 / 2205 / 2247 | 49 / 81 / 105 |

Bone longitudinal SoS is linear in frequency, fitted at 270 kHz and 836 kHz ([Pichardo 2011](https://doi.org/10.1088/0031-9155/56/1/014)). PRESTUS sets bone shear to 0:

| Tissue | Shear SoS [m/s] 300 / 500 / 650 kHz | Shear att. [Np/m] 300 / 500 / 650 kHz |
|---|---|---|
| Cortical bone | 1473 / 1542 / 1594 | 99 / 164 / 214 |
| Trabecular bone | 1270 / 1313 / 1345 | 99 / 164 / 214 |
| All soft tissues | 0 / 0 / 0 | 0 / 0 / 0 |

Versus PRESTUS defaults ([Medium properties](doc_medium.md)) at 500 kHz: cortical attenuation 81 vs ≈77 Np/m (13.3 dB/(cm·MHz), y = 1; about 6% higher), cortical SoS 2476 vs 2800 m/s (about 12% lower), brain SoS 1562 vs 1546 m/s (~1%).

<span style="font-size: 12px;">
  <strong>References for BabelBrain fitting functions:</strong><br>
  Pichardo, S. et al. Multi-frequency characterization of the speed of sound and attenuation coefficient for longitudinal transmission of freshly excised human skulls. <em>Phys Med Biol</em> <strong>56</strong>, 219–250 (2011). https://doi.org/10.1088/0031-9155/56/1/014<br>
  Pichardo, S. et al. A viscoelastic model for the prediction of transcranial ultrasound propagation: application for the estimation of shear acoustic properties in the human skull. <em>Phys Med Biol</em> <strong>62</strong>, 6938–6962 (2017). https://doi.org/10.1088/1361-6560/aa7ccc<br>
  Goss, S.A., Johnston, R.L. & Dunn, F. Comprehensive compilation of empirical ultrasonic properties of mammalian tissues. <em>J Acoust Soc Am</em> <strong>64</strong>, 423–457 (1978). https://doi.org/10.1121/1.382016<br>
  Labuda, C. et al. Two-dimensional mapping of the ultrasonic attenuation and speed of sound in brain. <em>Ultrasonics</em> <strong>124</strong>, 106742 (2022). https://doi.org/10.1016/j.ultras.2022.106742
</span>

## Tissue Labels and Material Properties

By default each tissue gets one set of properties from `medium_properties`, written as BabelBrain's integer label NIfTI (`BabelViscoInput.nii.gz`). Tissues without a BabelBrain label are written as water.

**Label mapping**

| PRESTUS layer name | Default config indices | BabelBrain label | BabelBrain material |
|---|---|---|---|
| `water` | 0, 3, 6, 9, 10 | 0 | Water |
| `skin` | 5 | 1 | Skin |
| `skull_cortical` | 7 | 2 | Cortical bone |
| `skull_trabecular` | 8 | 3 | Trabecular bone |
| `skull` (undifferentiated) | 4 | 2 | Cortical bone (fallback) |
| `brain` | 1, 2 | 4 | Brain |
| *(focal voxel)* | – | 5 | – |

Label 5 marks the focus voxel. It is BabelBrain's only way to specify the target.

### CT-Based Skull Properties (`ct_mapping`)

Set `babelbrain.ct_mapping` to let BabelBrain derive per-voxel skull properties from a pseudo-CT (HU) volume with one of its fitted models (BabelBrain's `bUseCT` mode):

```yaml
babelbrain:
  ct_mapping: Webb-Marsac   # one of the values in the table below
  ct_use_petra: false        # true only for PETRA-sequence-derived pCT
  water_amp_calibration: true # scale output to the target free-water ISPPA
```

Use this path only to apply BabelBrain's own HU conversion to raw UTE or PETRA pseudo-CT HU values. PRESTUS writes a raw-HU `CT.nii.gz` on the same grid as `BabelViscoInput.nii.gz`, and the bridge converts it to the format BabelBrain expects.

The bone on this path is the PRESTUS skull layer. BabelBrain's own bone definition (an HU threshold on its grid) is not applied by the bridge; set `pct.bone_source: hu` to define the skull that way, or `mask` to use a saved bone mask ([Skull from the CT](doc_pseudoCT.md#skull-from-the-ct)). The bridge then raises bone HU below `pct.bone_hu_threshold` to it on BabelBrain's grid, as BabelBrain does, because the resampling mixes soft-tissue HU into bone-edge voxels.

BabelBrain's Marsac density scales by the largest bone HU in its table. In a native BabelBrain run that table covers the whole head; through the bridge it would cover only the simulation domain. The bridge therefore passes the PRESTUS whole-head reference ([The HU Reference](doc_pseudoCT_property-mappings.md#the-hu-reference)) as the table's maximum. To emulate a native run, set `pct.hu_reference_percentile: 100`, `pct.interpolation: spline`, and `pct.limits.hu_max: 2100`. A native run resamples the CT onto its own grid before taking the maximum, so its value includes the spline overshoot at the skull surface, which the reference on the segmentation grid does not; `pct.hu_reference_resolution_mm` set to BabelBrain's grid step reproduces it.

!!! warning
    Do not set `ct_mapping` when PRESTUS's pCT pipeline (`pct.enabled: 1`) has already converted HU to acoustic properties. BabelBrain then receives only a label mask and uses its own bone values, which differ from the pCT-derived values k-Wave uses. PRESTUS's per-voxel acoustic properties cannot be passed to BabelBrain (see [Limitations](doc_babelbrain.md#limitations-and-unsupported-prestus-modes)).

| `ct_mapping` value | Density | SoS | Attenuation | PRESTUS equivalent |
|---|---|---|---|---|
| `Webb-Marsac` (default) | Marsac linear HU→ρ | Webb 2018 tabulated (freq-dep.) | Webb 2021 tabulated (freq-dep.) | `marsac-babelbrain` density, `webb` sound speed and attenuation (attenuation without q<sub>c</sub>) |
| `Aubry` | Porosity→ρ (Aubry 2003) | Porosity→SoS | Porosity→Att | `aubry` |
| `Pichardo` | HU→ρ air–tissue linear | Density→SoS (Pichardo 2D map, freq-dep.) | Density→Att (Pichardo map) | none |
| `McDannold` | HU→ρ air–tissue linear | Density→SoS (McDannold polynomial) | Density→Att (McDannold polynomial, freq-dep.) | none |
| `Marsac-Aubry` | Marsac linear HU→ρ | Porosity→SoS (Aubry) | Porosity→Att (Aubry) | partial |
| `Pichardo-Marsac` | Marsac linear HU→ρ | Pichardo 2D map (freq-dep.) | Pichardo map | none |
| `McDannold-Marsac` | Marsac linear HU→ρ | McDannold polynomial | McDannold polynomial (freq-dep.) | none |

PRESTUS `medium_pct_*` functions support density `k-plan`, `k-wave`, `marsac`, `marsac-babelbrain`, `aubry`, and `air-tissue`; sound speed `k-plan`, `marsac`, `aubry`, `webb`, `pichardo`, and `mcdannold`; attenuation `k-plan`, `mueller`, `aubry`, `webb`, `pichardo`, and `mcdannold` ([Property mappings](doc_pseudoCT_property-mappings.md)). For results comparable across backends, use `Webb-Marsac` (closest to the PRESTUS `marsac`/`mueller` combination) or `Aubry` (available in both).

- **PETRA density.** `ct_use_petra: true` makes `Webb-Marsac` use Björn's density formula (`SimNIBS_PETRApct_Density`). Use it for pseudo-CT from a PETRA sequence via [SimNIBS](https://simnibs.github.io/simnibs/build/html/index.html)'s ZTE/PETRA pipeline. For UTE-derived pCT (the PRESTUS default), leave it `false`.
- **Raw UTE as CT input.** A raw UTE volume can be passed as `pseudoCT` with `ct_mapping` set. UTE intensities lack standard 120 kVp HU calibration, so `Webb-Marsac` and `McDannold` give inaccurate results. `Aubry` is the most defensible choice for raw UTE. The recommended route is PRESTUS's pCT pipeline (UTE → calibrated pseudo-CT HU → PRESTUS `marsac`/`aubry` mapping).

#### Internal Air (`air_regions`)

On the pseudo-CT path, BabelBrain can treat internal air (frontal and sphenoid sinuses, mastoid cells) as a pressure-release reflector (reflection close to −1). Its solver does not update the fields inside the air mask. BabelBrain's own GUI enables this by default for CT/ZTE/PETRA input.

```yaml
babelbrain:
  ct_mapping: Webb-Marsac
  air_regions: true
  air_hu_window: [-1200, -400]   # default, as BabelBrain
  air_keep_boundary: false
```

PRESTUS writes `<target>_<f>kHz_<ppw>PPW_AirRegions.nii.gz` next to the CT input, on the same grid, and the log reports the voxel count and volume of internal air kept. Check `AirRegions.nii.gz` against the anatomy: in PRESTUS's UTE pseudo-CT, air cells inside the skull mask are captured only where their mapped HU falls in the window.

**How the air mask is built, and its limits**

1. Take the pseudo-CT voxels within `air_hu_window`.
2. Apply a 3-voxel median filter and find connected components (26-connectivity).
3. Remove the largest component (the air around the head).
4. With `air_keep_boundary: false`, also remove every component that touches the domain boundary. This removes exterior air split by the cropped grid, including the transducer–scalp gap (−1000 HU in the pseudo-CT but water in the simulation).

Limitations:

- **Pseudo-CT only.** On the label path there is no air label, so the option is ignored with a warning.
- **Free-water runs never get the mask.** The water baseline and `water_amp_calibration` are unaffected.
- **Boundary air is dropped.** Set `air_keep_boundary: true` if the domain edge cuts a sinus and no exterior air touches the boundary.
- **No k-Wave equivalent.** k-Wave and BabelBrain differ wherever the beam meets internal air.
