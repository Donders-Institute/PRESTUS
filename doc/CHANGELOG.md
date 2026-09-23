# Changelog

Notable changes to this project are documented here. 

---

## v0.7-pre (development)

Adds defacing, an acoustic FOV crop, multi-transducer incoherent summation, GUI updates for Placement, Multi-Transducer, and Calibration tabs, and a standalone PyQt6-based Python GUI. Introduces BabelBrain and jwave as selectable acoustic solver backends; BabelBrain also contains a CT path for HU-based skull property mapping, and a thermal wrapper. Principal-axis grid alignment becomes the new default for multi-transducer setups. Three additional pCT skull-property algorithms are added (Webb, Pichardo, McDannold). The transducer calibration library moves to `config/calibration/` and gains a version manifest, geometric-baseline entry generation, and geometry-only calibration mode. HPC post-run reporting is automated via a lightweight dependent SLURM job. Several fixes address MI computation, thermal crashes, multi-frequency transducer handling, and path-rename regressions across NeuroNav and calibration workflows. Async multi-transducer simulations are stabilised with a shared grid preproc stage and ras_plus enforcement. Pulse amplitude ramping is modelled in thermal simulation. TUS_entry is added as a new transducer placement mode.

#### BabelBrain acoustic and thermal backend

BabelBrain can now be used as a drop-in acoustic solver via `acoustic_backend: babelbrain` (see [doc_babelbrain.md](doc_babelbrain.md)).

- **Solver bridge** — `babelbrain_wrapper` translates the PRESTUS grid, medium masks, and transducer geometry into BabelBrain's NIfTI input format, invokes the Python FDTD engine, and returns `sensor_data` in the same format as the k-Wave backend
- **Transducer preset resolution** — transducer serial is automatically matched to a named BabelBrain preset; novel transducers are written to a temporary CSV so BabelBrain can handle them without a matching built-in entry
- **Calibrated phase export** — optimised element phases from the transducer library are exported as an `OptimizedWeightsFile` (H5); geometric (Rayleigh) focusing is used as a fallback when no library entry is present
- **CT path** — `babelbrain_wrapper` optionally accepts a pseudo-CT HU volume and writes it as `CT.nii.gz` alongside `BabelViscoInput.nii.gz`; activated by setting `babelbrain.ct_mapping` (e.g. `'Webb-Marsac'`, `'Pichardo'`); `babelbrain.ct_use_petra` selects the PETRA-specific density formula
- **Non-CT skull attenuation fix** — on the binary-mask (non-CT) path BabelBrain built its per-material Q-correction array as `np.ones(...)`, applying *no* Q correction to bone; in the relaxation FDTD model this drove the skull's target Q to infinity → near-zero skull attenuation and inflated ISPPA (~30 W/cm² vs ~11 W/cm² for k-Wave / jwave with identical settings). `run_babelbrain_case.py` now reproduces BabelBrain's own CT-mode bone treatment on the non-CT path (bone `ShearSoS = 0`, `QCorrArr[bone] = QCorrection`) via a runtime patch of `Step1_InitializeConditions` — no edits under `external/`, no throw-away run to discover the beam-cropped domain. Disable with `--no_skull_atten_fix` to reproduce the un-patched behaviour
- **PPW consistency** — points-per-wavelength is resolved from the simulation frequency and grid resolution and passed consistently to both the BabelBrain input NIfTI filename and the write step
- **Provenance** — BabelBrain version and git describe string are stored in the acoustic provenance cache so cached results are invalidated on submodule updates
- **Simulation report** — backend, transducer preset, PPW, and GPU compute backend (CUDA / OpenCL / Metal / CPU) reported in the HTML simulation report; material properties (density, SoS, attenuation per tissue type) read from DataForSim.h5 and displayed when `acoustic_backend: babelbrain`
- **Grid registration and HDF5 output parsing overhauled** — dimension ordering corrected (MATLAB h5read Fortran-order vs BabelBrain C-order), resampling anchored on TargetLocation + focus_pos; same fix applied to all volumetric thermal output fields
- **Thermal timing mapping corrected** — `ptd`→DurationUS, `ptri−ptd`→DurationOff, `ptrd/ptri`→Repetitions, `post_ptri_dur`→PostPause (previously `ptrd` was used as DurationUS, producing a single-block simulation instead of repeated trains)
- **Preset resolution improved** — folder name mapping and element counts corrected for CTX_500, CTX_250, and R15287
- **Focal-anchor ground-truth check** — `read_babelbrain_output`/`read_babelbrain_thermal` now cross-check the geometric `focal_length_mm` resample anchor (preset `FocalLength` / `TxMechanicalAdjustmentZ` / `ZSteering` chain) against BabelBrain's own `/TargetLocation` (the index where it placed the label-5 focal marker PRESTUS wrote into the input mask), correcting the anchor only if they actually disagree. `/TargetLocation` is computed by BabelBrain before it flips `p_amp`/`MaterialMap` along the depth axis, so it indexes the *raw, unflipped* `z_vec` — not the flip-corrected coordinate array used to index the (flipped) data arrays. An initial version of this check read it against the flipped array by mistake, which manufactured a spurious ~11.5mm shift on an otherwise-correct anchor (confirmed against a real run: `z_vec[TargetLocation]` equals the preset `FocalLength` exactly once indexed correctly) — fixed to index the raw array directly. `resample_babelbrain_vol` now takes an explicit interpolation `method` (default `'linear'`; MaterialMap/label resampling must use `'nearest'`)
- **Skull label containment fix** — `preproc_medium_mask` now smooths the combined (skull+cortical+trabecular) extent once and constrains the trabecular assignment to that combined shell, instead of smoothing cortical and trabecular independently. Independent smoothing let their boundaries disagree at the smoothing edge, letting trabecular (a physically-interior/diploe layer) end up directly adjacent to skin with no cortical between them in ~15% of trabecular voxels in one benchmark run
- **BabelViscoInput.nii.gz header fix** — `write_babelbrain_input` previously built a purely diagonal (axis-aligned) affine from the T1 origin, discarding the acoustic grid's actual head-alignment rotation (`grid.mode` defaults to a beam/focal-axis alignment that rotates the head relative to the T1) and any acoustic-FOV-crop voxel offset — the written NIfTI's header did not correspond to the T1/whole-brain medium mask it was derived from, making the skull layer look sparse/scattered when overlaid in a viewer even though the underlying label volume was fine. A new shared helper (`grid_to_ras_affine.m`) composes the grid's head-alignment transform, the FOV-crop offset (now threaded through via `planimg.fov_offset_ac`), and the load-time axis flip into the correct affine; the debug MaterialMap NIfTI (see above) uses the same helper instead of `niftiwrite`'s default identity header
- **FOV/domain-coverage debug figure** — with `simulation.debug: 1`, `read_babelbrain_output` writes `debug/babelbrain/sub-<id>_babelbrain_fov_debug<output_affix>.png`/`.csv` (BabelBrain-native pressure profile + its own `/MaterialMap` tissue-index profile + the PRESTUS-resampled profile, zero-filled region shaded) and `sub-<id>_babelbrain_material_map<output_affix>.nii.gz` (BabelBrain's `/MaterialMap`, resampled with nearest-neighbor onto the PRESTUS grid for overlay against `medium_masks`)
- **`*_Sub_NORM.nii.gz` gated behind `simulation.debug`** — BabelBrain's `Step10_GetResults` unconditionally calls `ResaveNormalized` to write a brain-masked, max-normalized copy of `FullElasticSolution_Sub.nii.gz` for its own GUI viewer; PRESTUS only ever reads `DataForSim.h5`, so `run_babelbrain_case.py` now only writes this file when `simulation.debug: 1` (new `--debug` flag, wired from `babelbrain_wrapper.m`) — useful there for visualizing where the brain/skull mask boundary falls relative to the solved field. `ResaveNormalized` is still patched unconditionally for its negative-index bounds-clamp fix
- **`grid_to_ras_affine` numerical cleanup** — chaining `planimg.inv_transf` (built via `inv()`) with `planimg.t1_header.Transform.T` and the FOV-crop offset left floating-point noise (~1e-14) in the composed matrix's last column, which `affine3d` rejects outright (`Error using images.internal.affine3dImpl/set.T ... final column ... must consist of zeroes`) — `write_babelbrain_input` crashed on layered runs as soon as this path was exercised. `grid_to_ras_affine` now zeroes the last column explicitly after composing the transform

#### jwave acoustic backend

jwave (JAX-based) can be used as an acoustic solver via `acoustic_backend: jwave` (see [doc_jwave.md](doc_jwave.md)).

- Supports `time_domain` and `helmholtz` solvers (Helmholtz is the primary workflow), configurable PML size, and `single`/`double` floating-point precision
- Translates PRESTUS kgrid/medium/source structs to HDF5, invokes `run_jwave_sim.py`, and returns `sensor_data` in the same format as the k-Wave backend
- Axisymmetric grids and kWaveArray sources are not supported; both raise a clear error at pipeline entry
- `benchmark_jwave` test class compares k-Wave CPU, jwave time-domain, and jwave Helmholtz on a water phantom (wall-clock time, peak pressure, focal FWHM, field correlation)
- **Helmholtz convergence diagnostics** — spectral residual reported after each solve; optional time-domain FOV benchmark (`jwave.benchmark_fov_mm`) cross-validates GMRES convergence on a cropped subgrid; time_domain retained for cross-validation
- **Helmholtz lossless-water fix** — `alpha_coeff = 0` (truly lossless) made the Helmholtz operator near-singular, causing GMRES to stall and underestimate peak pressure by ~18% in free-water benchmarks. `run_helmholtz()` now applies an attenuation floor of 0.002 dB/(cm·MHz) before building the jwave Medium; this is the physical value for distilled water and has negligible acoustic effect (< 0.006 % over 100 mm at 500 kHz)
- **Grid crop margin** — `preproc_crop_grid` now reserves `pml_size` voxels when the jwave backend is active (jwave embeds its PML inside the declared grid N; k-Wave adds it outside)
- **Time-domain (PSTD/FD) Stokes absorption sign fix** — `_absorb_tau_from_power_law` in `run_jwave_sim.py` was missing a negative sign, making the Stokes (y=2) absorption correction energy-injecting instead of dissipative; this caused a persistent NaN on small, PML-dominated FOVs that survived several rounds of CFL/stability-margin tuning, since the instability was structural rather than a matter of insufficient headroom. Fixed; `auto_pad_for_pml`'s stencil-clearance floor and `jwave_wrapper.m`'s `jw_absorb_safety` margin (added while chasing this bug) are kept as legitimate, independent hardening. Both time-domain solvers now also track `max|p|` at every step (not just within the recording window) so a future instability's onset step is visible directly in the log

#### Acoustic simulation — memory optimisations

- Large-aperture transducers and high-resolution grids no longer cause out-of-memory errors during source construction; source weights are now computed one element at a time rather than all at once
- `matlab_cpu` simulations with CW sources use substantially less RAM for the source signal, making CPU-based runs feasible for large arrays that would otherwise require many GB of storage
- GPU out-of-memory errors during source setup are now caught and handled gracefully; simulations fall back to CPU automatically with a warning
- ⚠️ **Fixed:** the lazy CW source (`CWSourceMatrix`, `matlab_cpu`/`matlab_gpu`) now supports k-Wave's source scaling (`source.p .* s` and the per-voxel loop for heterogeneous media); previously these runs errored before the time loop. Row index 58 was also misread as `:` (its character code), and single-element transducers produced a sparse source that could not be cast to single. Results match the dense source matrix

#### Acoustic simulation
- `resolve_acoustic_fov_params` resolves `'auto'` sentinel values for `grid.acoustic_fov_length_mm` and `grid.acoustic_fov_diameter_mm` from focal geometry; `has_acoustic_fov` flag gates all downstream crop/restore logic
- `resolve_water_baseline` encapsulates the decision of whether to run a free-water baseline simulation and stores the result in `acoustic_provenance.freefield_isppa_wcm2`
- `prepare_thermal_domain` restores full-grid state (or adjusts `planimg` affines for thermal-in-FOV mode) after an acoustic FOV crop
- `assemble_pfov_fields` back-projects N per-FOV acoustic results onto the shared full grid and combines them via incoherent intensity sum, producing a standard PRESTUS cache compatible with the thermal pipeline
- `acoustic_analysis` accepts an optional `acoustic_provenance` struct; per-transducer Ipa at each focus stored as bracketed CSV strings via internal `vec2csvstr` helper
- `resolve_thermal_grid` selects between an independent thermal grid and direct acoustic-grid reuse; centralises the logic previously spread across `prestus_pipeline`
- Beam-axis-aligned partial FOV (pFOV) crop reduces memory and compute for deep targets; renamed from "limited FOV" throughout; FOV direction corrected for ras_plus mode
- pFOV crop now compatible with async multi-transducer mode: per-transducer pFOV caches are backprojected before field combination
- Intensity plot orientation corrected in ras_plus mode; matrix transducer beam direction generalised beyond +Z
- Multi-transducer incoherent field assembly for independent-element arrays
- External acoustic NIfTI can now be passed directly as input to the thermal pipeline
- Per-transducer operating frequency supported in kWaveArray multi-transducer mode
- Exit-plane mask now unions contributions from all transducers (was transducer(1) only)
- FOV length and diameter derived automatically from focal geometry; post-hoc water simulation runs independently per transducer
- `planimg` affine sidecar serialised alongside the medium cache, enabling stage-3+ resumption without re-running the head pipeline
- ⚠️ **Fixed:** MI is now computed from peak rarefaction pressure (`abs(p_min_all)`) per IEC 62359; previously `p_max_all` was used
- ⚠️ **Fixed:** simulation end time (`t_end`) now derived from the minimum sound speed (soft tissue / water, ~1500 m/s) rather than the maximum (bone, ~2800 m/s); the previous behaviour truncated the simulation ~1.9× too early, causing an artefactual pressure drop before steady state was reached

#### Pseudo-CT acoustic property mapping

Three new algorithms added for density, sound speed, and attenuation mapping (see [doc_pseudoCT.md](doc_pseudoCT.md)):

- `air-tissue` density — linear HU→ρ calibrated to air and cortical bone; matches `HUtoDensityAirTissue` in BabelBrain; required input for `pichardo` and `mcdannold` workflows
- `webb` sound speed — linear HU→c from Webb et al. (2018); scanner-specific slope and intercept configurable via `pct.webb_sos_slope` / `pct.webb_sos_intercept`
- `webb` attenuation — exponential HU model from Webb et al. (2020); scanner-specific coefficients configurable via `pct.webb_att_alpha0`, `.webb_att_beta`, `.webb_att_c`
- `pichardo` sound speed and attenuation — frequency-dependent 2D ρ–freq interpolation from the `MapPichardo.h5` look-up table in the BabelBrain submodule
- `mcdannold` sound speed and attenuation — polynomial ρ-based model (McDannold et al.)

#### Placement
- **TUS_entry placement mode** (`placement.mode: 'tus_entry'`) — integrates the TUS_entry R toolbox (Atkinson et al.) to select a scalp entry point minimising coupling geometry error (SD of transducer-disk-edge to scalp distances); complementary to PlanTUS which optimises acoustic objectives

#### Head pipeline / grid
- Pipeline grid, medium, and thermal grid stages extracted into standalone functions (`run_grid_stage`, `run_medium_stage`, `resolve_thermal_grid`)
- **Focal-axis alignment** is now the default grid orientation (`grid.mode: focal_axis`): rotates the grid so the mean focal unit vector of all transducers aligns with z; for bilateral setups this equals the inter-transducer axis (identical to the old `principal_axis` result); for angled multi-transducer setups z correctly points toward the beam convergence zone rather than the inter-transducer axis. `principal_axis` is retained as a deprecated alias and emits a warning
- Grid cropping now includes the bowl geometry of every configured transducer in the crop mask (was transducer(1) only), preventing the contralateral side from being over-extended
- Automated defacing via pydeface, with optional integration into the segmentation pipeline
- Defaced regions shown transparently in placement and 3D head views

#### Transducer library
- Scan equipment, store validated parametric models, plot and compare entries
- `calibration.focal_reference` (`flhm_3db` default, `flhm_6db`, `peak`) sets which axial-profile convention a quoted focal depth refers to; used by `compute_flhm_center`/`compute_focal_metrics` and recorded in the calibration sidecar (see [doc_calibration.md](doc_calibration.md))
- Joint-depth phase/amplitude fit for more robust free-field characterisation
- Transducers resolved automatically from a serial-number-based equipment library; DS-agnostic `{serial}.yaml` fallback used when no DS-specific entry exists, with a warning
- `build_geometric_library_entry` creates a manufacturer-report baseline entry (no empirical data required) as an immediate-use fallback before characterisation
- `config/calibration/library_version.yaml` — version manifest listing all combo YAMLs with calibration metadata; follows semantic versioning

#### Multi-transducer async (experimental)
- Shared grid+medium preproc stage (Stage 0) runs over all transducers before per-transducer acoustic jobs, ensuring a common voxel grid for field combination
- `grid.mode: ras_plus` enforced at startup with a clear error and config snippet
- pFOV crop compatible with async: per-transducer pFOV-cropped caches backprojected before field combination
- HPC: acoustic jobs depend on Stage 0 via SLURM job dependencies

#### Thermal simulation
- **Pulse amplitude ramping** — `timing.pulse_ramp_dur` and `timing.pulse_ramp_shape` (`'linear'` or `'tukey'`) model symmetric rise/fall ramps; the ramp_factor (ISPPA/ISPTP) propagates into duty cycle products and the volumetric heat source Q
- Focal-axis temperature extraction and plot orientation corrected for ras_plus mode
- Thermal FOV can optionally be cropped to match the acoustic FOV (`thermal.match_acoustic_fov`), reducing memory for off-axis targets
- `thermal_grid_setup` exposed as a standalone utility
- ⚠️ **Fixed:** three runtime crashes in `thermal_analysis` corrected
- ⚠️ **Fixed:** initial temperature (`temp_0`) now correctly loaded from cached heating files

#### NeuroNav
- `placement.mode = 'tus_entry'` integrates the TUS_entry R toolbox (Atkinson et al.) as an optional placement backend: selects a scalp entry point that minimises coupling geometry error (SD of transducer-disk-edge to scalp distances) subject to an optional depth constraint; see [doc_placement.md](doc_placement.md#tus_entry-coordinate-selection)
- Multi-coil ingestion supported via `parameters.neuronav.coil_map`
- ⚠️ **Fixed:** T1 NIfTI lookup now works reliably across all conversion functions
- ⚠️ **Fixed:** `neuronav_convert_MNI_to_native` / `neuronav_convert_native_to_MNI` updated to use `parameters.path.seg` (was `parameters.seg_path`)

#### Sequential / Reporting
- HTML simulation report shows a per-transducer Ipa table for `asyncJoint` runs, listing each transducer's Isppa at its own focus sampled from the joint combined intensity map
- Run labels derived automatically from output affixes in consolidated reports
- Post-job resource usage reported via SLURM `sacct`, including a CO₂e footprint estimate
- HTML simulation report extended with BabelBrain and jwave backend sections
- ⚠️ **Fixed:** uncertainty pipeline CSV saved to correct path (no extra `tabular/` subdirectory)

#### GUI
- New Placement, Multi-Transducer, and Calibration tabs
- Transducer library browser tab
- Lazy tab initialisation (tabs load on first selection, reducing startup time)
- Parameter values seeded into unbuilt tabs on config load
- ⚠️ **Fixed:** NaN placeholder fields no longer write `0` into `grid.source_ppw`, which caused an Inf time step
- ⚠️ **Fixed:** 2D phantom / axisymmetric coordinate handling in collect and load
- ⚠️ **Fixed:** thermal timing enable state re-synced on config load

#### Infrastructure
- `doc_limited_fov.md` — user guide for the acoustic FOV crop feature
- `config/transducer/IS_PCD15473_01001_IGT_32_ch_comb_10_ch.yaml` — transducer config for the IGT 32-channel combined 10-channel array
- `iff(cond, a, b)` — inline conditional utility for one-liner expressions
- `compute_ramp_factor(pd, ramp_dur, shape)` — returns ISPPA/ISPTP ratio for linear or Tukey amplitude ramps applied to a pulse
- Grid-setup and medium-setup blocks extracted from `prestus_pipeline` into standalone functions `run_grid_stage`, `run_medium_stage`, and `load_grid_cache_file`; both cache-hit and cache-miss paths share identical load logic
- `overwrite_files` setting now respected by `simnibs_t1_to_mni`
- Compromised polyfill.io CDN replaced with Cloudflare drop-in
- `pydeface` added to environment setup script
- ⚠️ **Fixed:** SimNIBS path guard added in `subject2mni_coords_LDfix` to prevent crashes when SimNIBS is not on the path

---

## v0.6.1
*(2026-05-22)*

Extends the calibration pipeline with a unified entry point, parametric model store, and HPC dispatch. Introduces async multi-transducer simulations, a sequential simulation dispatcher, PlanTUS-based automated placement, and RAS as the canonical coordinate space. Drops the Image Processing Toolbox requirement.

#### Core / Config
- MATLAB Image Processing Toolbox no longer required — morphological ops and resampling use base-MATLAB equivalents
- `transducer.target_isppa_wcm2` migrated from calibration struct to transducer struct
- Optional `transducer.name` field added to transducer config

#### Calibration
- Unified `fit_source_parameters` entry point with HPC dispatch, multi-intensity mode, and parametric model store
- Geometric steering mode with per-element hardware correction
- Multiple analytical forward models — O'Neil and Rayleigh
- Transducer library (draft)
- Free-water correction moved post-analytical to make optional
- Free-water validation can be deactivated or limited to first/last target ISPPA
- Pluggable analytical forward model interface
- Option to use manufacturer phases to initiate global search (`opt_use_initial_phases`)
- ⚠️ **Fixed:** amplitude scaling bug present since early PRESTUS versions — correction factor is the square root of intensity ratios, not the ratio itself

#### Placement
- PlanTUS placement mode wraps the external PlanTUS optimiser for automated transducer targeting — **experimental** (see [doc_placement_plantus.md](doc_placement_plantus.md))
- T1 overlay plot auto-generated after heuristic and PlanTUS placement
- Browser-based transducer alignment viewer (*experimental*)

#### Acoustic/Thermal
- Async multi-transducer pipeline — **experimental** (see [doc_multitransducer.md](doc_multitransducer.md))
- Sequential simulation dispatcher: runs multiple parameter configurations in series, chains thermal simulations, and produces a consolidated multi-run report (see [doc_advanced.md](doc_advanced.md))
- ⚠️ **Fixed:** sequential simulation produced incorrect outputs prior to this release due to mismatching grid->T1w transforms and back of the adopted heatmaps; the heatmaps were not interpolated with the baseline free-water value leading to interpolation artefacts
- Uncertainty pipeline produces uncertainty bands across sequential multi-run parameter sweeps
- Complex pressure field optionally saved (`save_p_complex`) (*experimental*)

#### Head pipeline
- RAS established as canonical entry-point coordinate space; RAS+ orientation enforced at NIfTI load time; `ras_plus` grid mode rescales without rotation
- `tissuemask_binary` exposes a `water` mask field
- `smooth_img` gains an `off` method (no smoothing); exposed in GUI

#### Phantom
- `trans_pos_mm` / `focus_pos_mm` config fields allow resolution-independent transducer and focus position specification
- ⚠️ **Fixed:** phantom NIfTI outputs used incorrect voxel dimensions and inconsistent world-space orientation across pipeline stages

#### NIfTI / Transform
- Batched `tformarray` with optional parallel workers for MNI transforms
- Explicit FOV mask used in MNI outputs instead of zero sentinel (unreliable for pCT Hounsfield values)
- ⚠️ **Fixed:** pCT path layout corrected; NaN padding used instead of zero fill

#### Telemetry
- System resources reported at run start: OS, CPU model/cores, RAM, GPU, CUDA version

#### Reporting
- Public usage statistics (see [doc_usage_statistics.md](doc_usage_statistics.md))
- Group HTML report — `prestus_group_report_start` / `generate_group_report` aggregate per-subject CSVs and PNGs (across all subjects sharing a `medium` + `output_affix`) into a self-contained HTML with a mean ± SD exposure dashboard, subject roster, split acoustic (Intensity, MI) and thermal (maxT, riseT, CEM43) box plots, per-subject cards, and a client-side subject filter. Subjects are auto-discovered via `discover_group_subjects`; generated iteratively after each subject finishes when `modules.generate_group_report = 1`, or on demand. See [doc_group.md](doc_group.md#group-html-report).

---

## v0.6.0 
*(2026-04-30)*

Introduces matrix transducer support, an uncertainty quantification pipeline, and a rudimentary GUI. The pipeline is further integrated with automated pCT generation, placement dispatch, multi-ISPPA targeting, and a BIDS-inspired output structure. Parameter names and output paths have changed throughout.

> ⚠️ **Migration note:** Several parameter names and all output paths have changed in this release. Existing configs and post-processing scripts will need updating.
>
> **Parameters:**
> - `parameters.submit_medium` renamed to `parameters.platform`
> - Segmentation smoothing now specified in FWHM mm (parameter name changed)
> - Transducer field suffixes renamed consistently — check transducer config fields
> - Further parameter renames throughout — consult `config_default.yaml` for current names
>
> **Output paths:**
> - `parameters.io` fields renamed: `*_dir` → `dir_*`, `*_filename` → `filename_*`
> - Preprocessing cache files follow a consistent `res_N_descriptor` naming scheme
> - `space-` prefix removed from all NIfTI output filenames
> - `output_affix` now appended at the end of all output filenames
> - HPC log directory renamed: `hpc_log/` → `log_hpc/`
> - Outputs reorganised into BIDS-inspired typed subdirectories (`acoustic/`, `thermal/`, `cache/`, `debug/`, etc.)
>
> **Directory and config naming:**
> - `toolboxes/` directory renamed to `external/`
> - Config files consistently prefixed with `config_`

#### Added

- [**feature**] Uncertainty quantification pipeline (see [doc_uncertainty.md](doc_uncertainty.md))
    - Three-variant simulation workflow (default / liberal / conservative medium properties) producing a unified HTML report with per-metric safety cards colour-coded against ITRUSST non-significant risk limits, acoustic and thermal range tables, side-by-side intensity and temperature maps, and an SVG temperature timeseries with uncertainty band
    - Supports MATLAB (sequential) and HPC (parallel, scheduler-dependent) execution with automatic resumption from partial runs
- [**feature**] Matrix transducer support (see [doc_transducer.md](doc_transducer.md), [PR #117](https://github.com/Donders-Institute/PRESTUS/pull/117))
    - Flexible element layout options: rectangular grid, Fibonacci, Fermat spiral, and Clover multi-array configuration
    - Element positions can be defined inline or loaded from file
- [**pCT**] pCT generation fully automated inside the pipeline; a T1w image, UTE image, and mapping algorithm are sufficient — pre-generation remains possible for debugging; pCT section included in HTML report when used
- [**calibration**] In-pipeline amplitude calibration to target ISPPA
- [**feature**] Multi-ISPPA targeting — run acoustics once, scale to N independent thermal targets via post-hoc pressure scaling (see [doc_multiintensity.md](doc_multiintensity.md))
- [**thermal**] Dual CEM43 output: both kWave and ISO formulations propagated through pipeline, CSV, NIfTIs, and report
- [**report**] ISO CEM43 shown in Exposure Dashboard and methods section when enabled
- [**feature**] Placement dispatch — unified entry point for manual, Localite, and heuristic transducer positioning (see [doc_placement.md](doc_placement.md))
- [**feature**] `prestus_config_init` — bootstraps a project-specific config directory and `config_default.yaml`; `load_parameters` now resolves the project config automatically so PRESTUS can be called from any working directory
- [**neuronav**] Unified Localite pipeline via `localite_matrix_to_positions`; `InstrumentMarker` support added; `tracker_to_bowl_mm` parameter
- [**feature**] PRESTUS GUI — tabbed parameter editor and simulation launcher (see [doc_gui.md](doc_gui.md))
- [**feature**] `prestus_group_start` — new entry point for group-level MNI-space analysis
- [**feature**] `segmentation_only` pipeline mode
- [**report**] MItc (mechanical index transcranial) reported across all intracranial tissues
- [**thermal**] `simulation.log_thermal_updates` — suppress verbose per-step k-Wave progress output (default: off)
- [**refactor**] `simulation_nifti` split into per-pipeline-stage calls; bone mask separated from pseudoCT NIfTIs
- [**telemetry**] Opt-in anonymous usage telemetry (disabled by default; see [doc_telemetry.md](doc_telemetry.md))

#### Changed

- [**I/O**] Subject-specific output subfolders enforced
- [**output**] Intensity output naming cleaned up: redundant `max` labels removed; naming aligned with standardised pressure and intensity metric conventions
- [**grid**] `pml_size` defaults to `'auto'`; effective PML size resolved and tracked in log
- [**config**] `simulation.debug` defaults to `0`; property-map and MNI NIfTI export flags added to `config_default.yaml`
- [**validation**] Function inputs consistently validated via MATLAB `arguments` blocks; `help <function>` exposes header information inline
- [**report**] Safety Dashboard renamed to **Exposure Dashboard**; water-medium limits now populated; pressure units corrected; dynamic pressure display added; CEM and temperature removed from water-medium report

#### Fixed

- [**paths**] `load_parameters` resolves `config_default.yaml` via an absolute path, preventing stray output folders
- [**paths**] Cross-platform SimNIBS path compatibility
- [**transducer**] Annular array element format and default thickness corrected
- [**nifti**] Property map export: corrected datatype mismatch, scalar `alpha_power` handling, and `FillValue` assignment
- [**hpc**] Transducer index suffix omitted from positioning filenames for single-transducer runs
- [**plots**] T1 intensity overlay orientation corrected using NIfTI sform; deprecated `LineSmoothing` property removed
- [**logging**] Post-hoc water simulation writes to its own independent log file; PPW warning no longer emits a backtrace
- [**debug**] Debug directory creation and plots skipped when `simulation.debug = 0`

---

## v0.5.0
*(2026-02-27)*

Major restructuring of the pipeline into discrete processing stages. Adds parallel multi-transducer support, an HTML simulation report, C++ backend, skull rubber-wrap, and heuristic transducer placement. Parameter names have changed and existing configs will need updating.

> ⚠️ **Migration note** — configs from v0.4.x will need updating:
> - `layer_labels` renamed to `layers`
> - Transducer and focus position fields renamed to `trans_pos` / `focus_pos` consistently
> - `grid_step_m` and `default_grid_size` removed
> - Grid expansion parameter renamed
> - Attenuation parameters renamed — consult `config_default.yaml`
> - Layer indexing field migrated from `layers` to `medium` throughout
> - Pipeline refactored into discrete processing stages — consult the documentation for current function entry points

#### Added

- [**feature**] C++ backend CPU & GPU support by @jkosciessa
- [**feature**] Parallel multi-transducer support for layered simulations by @dreamstimlab (see [PR #100](https://github.com/Donders-Institute/PRESTUS/pull/100))
- [**feature**] HTML simulation report — self-contained per-subject summary of acoustic and thermal outputs by @sirmrmarty
- [**feature**] Heuristic transducer placement — automated skull-surface search with configurable positioning criteria and ear exclusion (see [PR #109](https://github.com/Donders-Institute/PRESTUS/pull/109))
- [**skull**] Rubber-wrap skull mask inflation by @meijer-s (see [PR #103](https://github.com/Donders-Institute/PRESTUS/pull/103))
- [**thermal**] Flexible and standardised protocol timing specification aligned with TUSCalculator; protocol plot output
- [**thermal**] Additional thermal timeseries outputs (max per tissue, heating at focus)
- [**feature**] Per-step pipeline benchmarking: wall time, peak RAM, and disk usage logged at each stage
- [**grid**] Minimum PPW validation at source setup (`grid.min_ppw`, default 6); warns with the maximum `resolution_mm` that would satisfy the requirement
- [**calibration**] `calibration.opt_phase_precession` — optional constraint on element phases during amplitude/phase optimisation
- [**neuronav**] Updated Localite GUMMarkers parsing
- [**doc**] GitHub Pages [documentation](https://donders-institute.github.io/PRESTUS/) — parameters, functions, installation, quick start, and processing steps

#### Changed

- [**refactor**] Pipeline and core processing steps refactored; functions organised by processing stage
- [**axisymmetric**] Axisymmetric acoustic simulations now default to 3D thermal output

#### Fixed

- [**pCT**] kPlan skull density and sound speed regularised to water minimum
- [**acoustic**] GPU array correctly moved to CPU before NIfTI export
- [**calibration**] Unoptimized parameters no longer passed to Oneil recomputation
- [**savemat**] k-Wave source matrix saving can now be correctly deactivated; water simulations force-save matrix outputs

#### Deprecated/Removed

- [**feature**] Explicit `skull`-only simulation mode removed; equivalent behaviour is achieved by removing layers from the simulation setup

---

## v0.4.2 
*(2025-12-04)*

Refactors the calibration and acoustic profiling pipeline and introduces multi-element axisymmetry support.

#### Added

- [**calibration**] Acoustic profiling refactor with documented standalone calibration pipeline (see [doc_calibration.md](doc_calibration.md))
- [**calibration**] Multi-element axisymmetry support (experimental)

#### Changed

- [**default**] kWaveArray used by default

---

## v0.4.1 
*(2025-08-26)*

Post-release additions to v0.4.0: heat trace outputs, Localite neuronavigation integration, and 3D thermal as the default for axisymmetric runs.

#### Added

- [**thermal**] Heat trace saving
- [**logging**] MATLAB log written to output directory
- [**neuronav**] Localite coordinate processing

#### Changed

- [**axisymmetric**] 3D thermal simulation is now the default for axisymmetric acoustic runs

---

## v0.4.0 
*(2025-06-04)*

This release features major updates. Please consult the documentation for current parameter choices. Some features remain unstable and new bugs may pop up. Let's remain vigilant and try to quickly introduce hotfixes in those cases.

- [**Feature**] Acoustic profiling (beta) @MaCuinea
- [**Feature**] pCT->acoustic skull medium incl. multiple mapping algorithms (beta, doc) @eleonoracarpino @jkosciessa
- [**Feature**] Run consecutive simulations @KTZ228
- [**Feature**] GPU updates (doc) @jkosciessa
- [**Feature**] CEM43 ISO definition, see doc @jkosciessa
- [**Feature**] Modeling of absorption (as fraction of attenuation) and perfusion @jkosciessa. Be advised that this renders the heating estimates more liberal than in the past. See the documentation. Also see some comments in the issue.
- [**Feature**] axisymmetry simulations with 2D input grids (doc) @jkosciessa
- [**Feature**] simulations with computational phantoms (doc) @jkosciessa
- [**Documentation**] Initial draft of documentation (e.g., parameters) (see links above and here) @jkosciessa
- Many smaller bug fixes and QOL improvements (e.g., ability to disable .mat saving, more debugging outputs)
- DOC: minor typo in README.md by @mekman in #58
- Run simnibs segmentation code also supports SLURM by @MaCuinea in #62
- Prepare acoustic_profiling for merge by @MaCuinea in #69
- Semi-automated acoustic profiling by @MaCuinea in #71
- Added note to default config regarding ld_library_path in slurm by @KTZ228 in #70
- Merge pull request #70 from KTZ228/development by @MaCuinea in #77
- Fixed tutorial by @MaCuinea in #73
- Feature: consecutive simulations by @KTZ228 in #79
- Documentation, pCT, HCP updates & more, see above by @jkosciessa in #59
- The last one summarizes more than half a year of updates to code and documentation.

---

## v0.3.0

- PseudoCT feature by @eleonoracarpino and @jkosciessa
- Clean-up repository by @jkosciessa
- Incorrect allocation of attenuation values bug fix by @jkosciessa
- Improve correct labeling of tissue by @jkosciessa
- Improvement of 3D plots by @jkosciessa
- Check simulation stability and adjust time step accordingly by @MaCuinea and @jkosciessa
- Enable tissue-specific temp_0 specification by @jkosciessa
- Fix heating plot bug by @jkosciessa

---

## v0.2.0

- DOC: add GNU license by @mekman in #17
- Documentation changes and bug fixes by @KTZ228 in #14
- Disabled masking of neural tissue in skin in smooth_and_crop_layered by @KTZ228 in #19
- Adjusted create_Group_MNI_plots by @KTZ228 in #22
- Documentation changes by @KTZ228 in #24
- Minor changes to comments for heating setup by @jkosciessa in #29
- Figure saving updates, some bug fixes by @jkosciessa in #8
- Kenneth master by @KTZ228 in #33
- Group plots fixes and improvements by @KTZ228 in #26
- FIX: revert to working version by @mekman in #35
- DOC: reworked installation instructions by @mekman in #36

---

## v0.1.0

Initial in-house development version by @achetverikov

- Dedicated segmentation path option, minor bug fixes, updated simnibs call by @jkosciessa in #2
- fix small errors by @jkosciessa in #3
- adjusted skull values by @eleonoracarpino in #4
- Bug fixes and documentation changes by @KTZ228 in #6
- DOC: simplify readme.md wording regarding requirements by @mekman in #7
