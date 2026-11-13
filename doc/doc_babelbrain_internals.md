# BabelBrain Backend: Integration Internals

PRESTUS runs [BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain) through a bridge: `babelbrain_wrapper.m` writes the inputs and calls `run_babelbrain_case.py`, which patches BabelBrain at runtime and runs it, and the readers map the results back onto the PRESTUS grid. Most problems with this backend are misalignment between the two grids. This page documents each step of the bridge, for debugging a run (start with [Output Alignment](#output-alignment) and the [zero-pressure checklist](#zero-pressure-debugging-checklist)) and for changing the bridge; for using the backend see [BabelBrain backend](doc_babelbrain.md).

On this page:

- [Call Chain](#call-chain): the functions from the pipeline to BabelBrain, and the runtime patches the bridge applies.
- [GPU Backend Selection and Platform Support](#gpu-backend-selection-and-platform-support): `computing_backend` and `device_name` per platform.
- [Grid Resolution](#grid-resolution): why the input mask is resampled to BabelBrain's step.
- [Placement Correction and Drive Phases](#placement-correction-and-drive-phases): the placement correction, steering and phases, and the free-water and bone checks.
- [Calibration Weights File](#calibration-weights-file): what the weights store and the [HDF5](https://www.hdfgroup.org/solutions/hdf5/) layout.
- [Free-Water Amplitude Calibration](#free-water-amplitude-calibration): the correction of BabelBrain's source bias and the water baseline checks.
- [Inputs, Geometry and Modes](#inputs-geometry-and-modes): CT, air and water inputs, the single-transducer limit and mode-specific behaviour.
- [Output Alignment](#output-alignment): how the field lands on the PRESTUS masks.
- [BabelBrain Output Format](#babelbrain-output-format): the datasets read, axis flips, interpolation and the [NIfTI](https://nifti.nimh.nih.gov/) affine.
- [Zero-Pressure Debugging Checklist](#zero-pressure-debugging-checklist): a checklist for missing pressure.
- [Built-In Material Table](#built-in-material-table): BabelBrain's own tissue values, used only if material injection fails.

!!! note "At a glance"
    - **Tested with:** BabelBrain 0.8.1, the `external/BabelBrain` submodule pinned at upstream `main` commit `36d6a3e` (after the v0.8.1 tag). Other versions may break the runtime patches below; each run records the version it used ([version tracking](doc_babelbrain.md)).

## Call Chain

```
prestus_pipeline.m
  └─ acoustic_stage.m            ← dispatches on parameters.acoustic_backend (run_acoustic_backend)
       └─ babelbrain_wrapper.m   ← MATLAB orchestrator
            ├─ babelbrain_transducer_map.m   (preset resolution)
            ├─ write_babelbrain_input.m       (NIfTI label mask)
            ├─ export_babelbrain_weights.m    (calibration HDF5)
            └─ system() → run_babelbrain_case.py   ← Python bridge
                 └─ BabelBrain RUN_SIM().RunCases()
```

The bridge runs in `babelbrain.python_env`, writes `*DataForSim.h5` to `<dir_cache>/<subject_id>/` and patches BabelBrain at runtime (never `external/`). Without the skull-attenuation patch, injected bone attenuation is not realised (QCorr stays 1.0), although the material table looks correct.

??? note "Runtime patches"
    | Patch | Effect | Disable (bridge flag, not a config key) |
    |---|---|---|
    | `_inject_prestus_materials` | Overwrites `MatFreq` from `--medium_properties_json` (attenuation converted from dB/(cm·MHz^y) to Np/m) | runs only without `--no_skull_atten_fix` |
    | `_patch_noncat_skull_attenuation` | Bone shear = 0 (fluid bone), so the specified skull attenuation is realised | `--no_skull_atten_fix` |
    | `_patch_zdomain_reference` | `CalculateDomainZReference` returns 0 (see [Output alignment](#output-alignment)) | `--no_zdomain_fix` |

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
| `babelbrain.z_steering_m` set | pinned | Overrides both; the weight-delta reference follows `FocalLength + ZSteering`. |

`Single` and `BSonix` presets lack `TxDiam`/`FocalLength`, so the wrapper passes `--aperture` (max `elem_od_mm`) and `--focal_length` (`curv_radius_mm`).

??? note "Free-water check, bowl rim and bone Q-correction"
    **Free-water check.** With explicit phases, the BabelBrain axial profile must overlay the Rayleigh profile of the same phases. A rigid offset of near-field lobes and peak means the placement correction is missing. A lobe *shape* mismatch means the drive phases differ (e.g. weights delta referenced to the wrong aim).

    **Bowl rim.** A rim inside the domain makes BabelBrain snap the transducer back in whole voxels (log `Warning: Need to apply correction to reposition Tx for …`; non-zero `source_residual_mm`). `write_babelbrain_input.m` therefore trims `ceil((d_apex + bowl_depth)/dx) + 1` transducer-side slices (bowl depth 5.4 mm for R15473; `d_apex` negative if the apex is outside the grid). It warns if that slab contains tissue; it is zero-filled on paste-back.

    **Bone Q-correction.** Realised bone attenuation scales ~1/Q. A bone-slab phantom (6 mm cortical, 78.5 Np/m at 300 kHz, PPW 6) gave realised/specified attenuation 0.91 at Q = 1 (k-Wave 0.91), 0.29 at 3, 0.14 at 6, 0.087 at 10 and 0.063 at 14. BabelBrain's CT default of 3.0 gave transcranial Isppa 10.6 vs 7.6 W/cm² for k-Wave (sub-602 pFOV), so the non-CT path uses Q = 1. `babelbrain.q_correction` > 0 reaches `BabelFTD_Simulations_BASE.__init__` via `--q_correction` → `BABEL_PYTEST_QFACTOR`. MATLAB `setenv` does not work here, because `run_command` uses Java's ProcessBuilder with the JVM start-up environment.

## Calibration Weights File

The file stores deltas from BabelBrain's ring phases at its aim `FocalLength + ZSteering`:

```
w[n] = exp(i × (elem_phase[n] − geometric_phase_at_aim[n]))
```

??? note "File layout and import/export"
    - `/CALIBRATION` is a 1-D complex64 dataset ([h5py](https://www.h5py.org/) compound `{r, i}`) as in BabelBrain's `SaveToH5py`; every root item has the `type` attribute `ReadFromH5py` requires. `convert_weights_if_needed` converts the split `/CALIBRATION/real` + `/imag` layout.
    - `export_babelbrain_weights` writes one weight per ring; the bridge expands them to sub-elements.
    - `import_babelbrain_weights` reads `TxCalibration` files (one weight per sub-element) given the per-ring element counts (6th argument).

<a id="free-water-calibration-transfer"></a>
## Free-Water Amplitude Calibration

`water_amp_calibration` removes the few-% pressure bias of BabelBrain's source correction by comparing the free-water field with `rayleigh_axial_intensity.m`, and applies the same factor to the head simulation. The thermal input is not rescaled, because BabelBrain renormalises it to `Isppa`.

??? note "Why the factor transfers, and water baseline checks"
    BabelBrain's `Correction = 100/(100 − ExpectedError(CFL))`, ×√2 (`BabelIntegrationBASE`) leaves a few % pressure bias at CFL ~0.5. The factor transfers to tissue because the simulation is linear and the bias sits in the shared source chain. The transcranial timestep differs (sub-602 pFOV, R15473 300 kHz 6 PPW: CFL/water 0.300 → Correction 0.290 vs 0.471 → 0.465). The correction tracks it within ~2% in pressure (0.290/0.465 = 0.624 vs CFL ratio 0.637), ~4% in intensity.

    The thermal renormalisation is `PressureRatio = sqrt(Isppa·2ρc)/max(p_amp[brain])`.

    **Water baseline checks** (`babelbrain_wrapper.m`):

    - The `--water_only` run must get the same `--source_amp_pa` as the transcranial run. Otherwise BabelBrain drives at ~50 kPa and `freefield_isppa_wcm2` is ~(156.6/50)^2 ≈ 10x too low for R15473.
    - `Water_DataForSim.h5` is reused if present, unless `io.overwrite_files` is `'always'`/`'force'`. If the ISPPA-normalised field looks wrong, compare its mtime with the main `DataForSim.h5`.

## Inputs, Geometry and Modes

PRESTUS does not rotate its grid to the trajectory, so the beam must lie on z and one transducer is supported.

??? note "Water, CT and air inputs"
    - **CT.** `write_babelbrain_input.m` writes raw-HU `CT.nii.gz` on the `BabelViscoInput.nii.gz` grid. The bridge quantises it (uniform, 10-bit, BabelBrain's default) into the index map BabelBrain expects and writes `CT-cal.npz`. The wrapper passes `--use_ct`, `--ct_mapping` and optionally `--use_petra`.
    - **Air.** BabelBrain does not update the fields inside the air mask (`ReflectorMask`), so air needs no material and no tiny timestep. `bExtractAirRegions` is on only with `air_regions`, because PRESTUS otherwise supplies no air mask.
    - **Uncertainty.** Variants go in via `--medium_properties_json`; on the CT path only skin and brain apply.
    - **Water.** Simulated water comes from the module-level `Material` (the 1120 kHz `MatFreq` entry), so the injected `Water` line affects only `GetSmallestSOS`. Water is not injected, because it also sets the source factor `FactorConvPtoU`; the impedance difference is ~0.6%.

??? note "Geometry rationale and mode-specific behaviour"
    BabelBrain's own workflow rotates the grid instead (`BabelDatasetPreps.GetSkullMaskFromSimbNIBSSTL`). Multi-transducer support would need beam-aligned per-transducer domains (as `simulation.async_grid_mode: per_transducer` for k-Wave) and a thermal step that accepts a combined heat source.

    - **`multi_isppa`.** The thermal result is named after the shared acoustic file (`<stem>-ThermalField_AllCombinations.h5`), so thermal jobs are serialised on HPC.
    - **Post-hoc water.** The wrapper passes `--water_only` under the baseline's target name, which enables cache reuse.
    - **Limited FOV.** The wrapper disables `bTightNarrowBeamDomain` on a cropped mask, because its lateral footprint under-shoots the beam near the skull.

<a id="output-alignment-tissue-marker-anchoring-and-the-zdomainstart-fix"></a>
## Output Alignment

Field and `MaterialMap` must land exactly on `medium_masks`.

- **`ZDomainStart` = 0.** BabelBrain adds a skin-surface heuristic (`GlobalStartSkin − LineOfSightStartSkin`) to the source position in `CalculateRayleighFieldsForward`. PRESTUS already sets `TxMechanicalAdjustmentZ`, so the patch removes it. The log line "ZDomainStart neutralized: ... mm" shows the removed amount.
- **Marker anchor.** The readers anchor the resample on `/TargetLocation`, the label-5 marker at `focus_pos` (benchmark: label agreement 83% → 99%, bone 20% → 98%). Preset `FocalLength` (`focal_length_mm`) is the fallback for files without it.

`<target>_prestus_anchor.json` records realised-minus-requested source position. Expect ≤ 1 native voxel; `read_babelbrain_output.m` warns above 1.5, which points to a missing `ZDomainStart` patch.

## BabelBrain Output Format

`sensor_data.p_max_all` is `/p_amp` (peak pressure [Pa], `[Nx, Ny, Nz]`), or `/p_amp_refocus` with `do_refocusing: true`. Axes are permuted to `[lateral X, lateral Y, axial]` (dim 3 = beam).

??? note "Axis flips and interpolation"
    | Item | Handling (change both sides together) |
    |---|---|
    | Acoustic arrays | BabelBrain flips them (`np.flip(axis=2)`) before saving, but not `z_vec`; the reader uses `zdim = z_vec_bb(end:-1:1)`. |
    | `/TargetLocation` | Computed before the flip; index the raw `/z_vec`. |
    | Input mask | `write_babelbrain_input.m` pre-flips when the beam runs along +z (`focus_pos(3) >= trans_pos(3)`), so BabelBrain's flip on load cancels it; with the beam along −z, BabelBrain's flip alone orients it. The reader resamples with the same `beam_sign`. |
    | Thermal arrays | Saved in `z_vec` order; not reversed. Resampled onto `parameters_th` from `thermal_stage`, anchored on the thermal focus (`grid.thermal_focus_pos_exact` etc. from `thermal_grid_setup`, else `transducer(1).focus_pos`). Rises over `BaselineTemperature`, so outside voxels read baseline. |
    | Interpolation | Fields `'linear'`; `MaterialMap` `method='mode'`, because `'nearest'` drops sub-voxel cortical bone. |

The NIfTI affine comes from `grid_to_ras_affine.m`: `planimg.inv_transf`, `planimg.t1_header`, `planimg.fov_offset_ac`, then the dim-3 flip, in that order. A diagonal-only affine misaligns whenever `grid.mode` rotates the head (the default).

## Zero-Pressure Debugging Checklist

Inspect `/MaterialMap`, `/p_amp`, `/TargetLocation` and `/z_vec` in `*DataForSim.h5`, then:

| Check | Points to |
|---|---|
| Console "zero-filled (no BabelBrain data)" line non-zero | Domain too short; voxels outside it are zero-filled |
| `ZShrink_L` > 0, `XShrink_*`/`YShrink_*`, or fewer tissue voxels in `MaterialMap` than in `BabelViscoInput.nii.gz` | Narrow-beam mode cropped tissue |
| "BabelBrain native resolution" vs `grid.resolution_mm` | Mask not resampled to BabelBrain's step |
| No `/TargetLocation` in the file | Fallback anchor; field may be shifted along the beam |
| Skin→trabecular in a central-line label dump | Raw segmentation error |

`simulation.debug: 1` writes axis plots (`sub-<id>_babelbrain_fov_debug<output_affix>.png`/`.csv`) and `sub-<id>_babelbrain_material_map<output_affix>.nii.gz` to `<dir_output>/debug/babelbrain/`.

<a id="built-in-material-table"></a>
## Built-In Material Table

Literature fits tabulated every 5 kHz from 100 to 1125 kHz (`TranscranialModeling/BabelIntegrationBASE.py`). They apply only if material injection fails.

??? note "Built-in values and references"
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
