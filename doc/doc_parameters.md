# PRESTUS Parameters

This page lists every configuration parameter with its default, grouped by the YAML struct it belongs to.

On this page:

- [Subject and Platform (`subject_id`, `platform`)](#subject-and-platform-subject_id-platform): subject identifier and execution platform.
- [Simulation Type (`simulation`)](#simulation-type-simulation): medium, [k-Wave](https://github.com/ucl-bug/k-wave) code type, precision, interactive and debug modes, transducer coupling.
- [Simulation Backends (`acoustic_backend`, `thermal_backend`)](#simulation-backends-acoustic_backend-thermal_backend): which solver runs the acoustic and thermal steps, and its settings.
- [Data Paths (`path`)](#data-paths-path): input, segmentation, simulation and [Localite](https://www.localite.de/) folders.
- [Environment & Toolbox Paths (`startup`)](#environment-toolbox-paths-startup): [SimNIBS](https://simnibs.github.io/simnibs/build/html/index.html), [MATLAB](https://www.mathworks.com/products/matlab.html) and other external tool paths.
- [I/O Management (`io`)](#io-management-io): output affix, overwriting, cache checks and which matrices are saved.
- [Pipeline Modules (`modules`)](#pipeline-modules-modules): which pipeline stages run.
- [Transducer Specification (`transducer`)](#transducer-specification-transducer): ISPPA target, transducer type, frequency, library entry and positions.
- [Transducer Placement (`placement`)](#transducer-placement-placement): placement mode, skin gap, reuse and the per-mode keys.
- [Simulation Grid (`grid`)](#simulation-grid-grid): resolution, dimensions, axisymmetric mode and PML.
- [Segmentation (`segmentation`)](#segmentation-segmentation): options passed to SimNIBS [charm](https://simnibs.github.io/simnibs/build/html/documentation/command_line/charm.html), and defacing.
- [Head Model & Segmentation Preprocessing (`headmodel`)](#head-model-segmentation-preprocessing-headmodel): padding and mask preprocessing of the head model.
- [Pseudo-CT Skull Property Mapping (`pct`)](#pseudo-ct-skull-property-mapping-pct): pseudo-CT generation and its mapping to skull properties.
- [Tissue Layers (`layers`)](#tissue-layers-layers): which segmentation labels form each tissue layer.
- [Tissue Acoustic & Thermal Properties (`medium_properties`)](#tissue-acoustic-thermal-properties-medium_properties): per-tissue acoustic and thermal properties.
- [Sonication Timing (`timing`)](#sonication-timing-timing): pulse and pulse-train durations and repetition intervals.
- [Thermal Settings (`thermal`)](#thermal-settings-thermal): CEM43 variant, baseline temperatures and thermal sensor and logging.
- [Time Reversal (`time_reversal`)](#time-reversal-time_reversal): the extra point-source run per target.
- [Output Analysis (`analysis`)](#output-analysis-analysis): focus-area radius and MNI target or ROI for analysis.
- [Study Setup (`study`)](#study-setup-study): study manifest and output folder name.
- [High-Performance Computing (`hpc`)](#high-performance-computing-hpc): GPU, partition, time and memory limits for cluster jobs.
- [Transducer Calibration (`calibration`)](#transducer-calibration-calibration): keys of `config_calibration.yaml`.

Parameters are organised in nested structs that map directly to [YAML](https://yaml.org/) keys. `config_default.yaml` defines all parameters and their default values and should not be edited; a study config is merged on top and only needs the values that differ (see [Configuring a study](doc_configuration.md) for the full layering, project defaults and overrides). Parameters can also be set in MATLAB before calling `prestus_pipeline_start(parameters)`. Mandatory parameters have no meaningful default and must always be provided.

---

### Subject and Platform (`subject_id`, `platform`)

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `subject_id` | Subject identifier used for file naming and folder management. | none | Mandatory. Set before calling `prestus_pipeline_start(parameters)`. |
| `platform` | Execution platform. | `'auto'` | `auto` / `slurm` / `qsub` / `matlab`. `auto` detects [SLURM](https://slurm.schedmd.com/) or qsub and errors when neither is installed; use `matlab` to run in the current session. See [HPC](doc_hpc.md#where-a-run-executes). |

---

### Simulation Type (`simulation`)

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `medium` | Simulation medium type. | `'layered'` | `water` / `layered` / `phantom`. Mandatory. |
| `code_type` | k-Wave backend. | `'matlab_gpu'` | `matlab_cpu` / `matlab_gpu` / `cpp_cpu` / `cpp_gpu` / `cpp_mlx`. See [k-Wave](doc_kwave.md). |
| `precision` | Computational precision for acoustic and thermal simulations. | `'single'` | `single` / `double` |
| `interactive` | Interactive mode with prompts and evolving plots. | `0` | `1 = yes`, `0 = no`. If `0`, see `overwrite_files`. |
| `debug` | Verbose debug mode with additional intermediate outputs. | `0` | `1 = yes`, `0 = no` |
| `transducer_coupling` | How several transducers interact acoustically. | `'coherent'` | `coherent`: pressure fields superpose in one k-Wave run (interference). `async`: each transducer is simulated separately and intensities are summed incoherently before the thermal simulation. See [doc_multitransducer.md](doc_multitransducer.md). |
| `uncertainty_temp_high` | Medium properties of the Temp-H variant (bounds that maximise brain heating) when `simulation.uncertainty` is on. | `''` | A YAML file with a `medium_properties` block, merged onto the run's config: an absolute path, or a file name in `config/uncertainty/`. `''` = `config/uncertainty/config_medium_temp_high.yaml`. The `options.temp_high_config` argument of `uncertainty_pipeline` takes precedence. See [doc_uncertainty.md](doc_uncertainty.md#customising-pipeline-options). |
| `uncertainty_temp_low` | Medium properties of the Temp-L variant (bounds that minimise brain heating). | `''` | As `uncertainty_temp_high`; `''` = `config/uncertainty/config_medium_temp_low.yaml`; `options.temp_low_config` takes precedence. |
| `uncertainty_trans_high` | Medium properties of the Trans-H variant (bounds that maximise brain intensity; acoustic only). | `''` | As `uncertainty_temp_high`; `''` = `config/uncertainty/config_medium_trans_high.yaml`; `options.trans_high_config` takes precedence. |
| `uncertainty_trans_low` | Medium properties of the Trans-L variant (bounds that minimise brain intensity; acoustic only). | `''` | As `uncertainty_temp_high`; `''` = `config/uncertainty/config_medium_trans_low.yaml`; `options.trans_low_config` takes precedence. |

---

### Simulation Backends (`acoustic_backend`, `thermal_backend`)

Top-level keys (not nested under `simulation`). Pre-processing (segmentation, medium, grid, placement) and post-processing (analysis, [NIfTI](https://nifti.nimh.nih.gov/) export, report) are shared by all backends. See [Simulation backends](doc_backend.md) for how the backends compare and when results are directly comparable.

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `acoustic_backend` | Acoustic solver. | `'kwave'` | `kwave` / `kwave_elastic` / `babelbrain` / `jwave`. Validated in `load_parameters`. `babelbrain` takes one transducer, with the beam along the grid z axis (not `grid.mode: ras_plus`); see [BabelBrain geometry limits](doc_babelbrain.md#geometry-limits). `kwave_elastic` models shear waves in the skull (see [Elastic skull](doc_kwave.md#elastic-skull-model-shear-waves)). |
| `thermal_backend` | Bio-heat solver. | `'kwave'` | `kwave` / `babelbrain`. `babelbrain` requires `acoustic_backend: babelbrain`, one transducer and no `path.acoustic_input`; invalid combinations error at config load. [jwave](https://github.com/ucl-bug/jwave) has no thermal solver; pair it with `kwave`. |

#### BabelBrain Settings (`babelbrain`)

Used when `acoustic_backend` or `thermal_backend` is `babelbrain`, or when `io.export_babelbrain: true`. Full explanation in [BabelBrain](doc_babelbrain.md).

??? note "BabelBrain settings"

    | **Parameter** | **Description** | **Default** | **Comments** |
    |---|---|---|---|
    | `python_env` | BabelBrain [conda](https://docs.conda.io/) env (directory containing `bin/python`) or [Python](https://www.python.org/) executable. | `''` | Must be set; system `python3` lacks `BabelViscoFDTD`. |
    | `root` | BabelBrain source tree. | `''` | Empty = `external/BabelBrain` submodule. |
    | `computing_backend` | GPU/CPU backend. | `'auto'` | `auto` / `cuda` / `opencl` / `metal` / `cpu` (or `1`–`4`). `auto`: Metal on macOS, else CUDA if `nvidia-smi` finds a GPU, else CPU. The thermal engine has no CPU backend. |
    | `device_name` | Device name. | `'auto'` | Set explicitly when several GPUs are present. |
    | `tx_system` | BabelBrain transducer preset. | `''` | Empty = auto-detect from serial/geometry. See [Transducer mapping](doc_babelbrain.md#transducer-mapping). |
    | `ppw` | Points-per-wavelength floor for BabelBrain's own FDTD grid. | `12` | Lower = faster, more numerical dispersion. |
    | `do_refocusing` | Apply BabelBrain's refocusing step. | `false` | Thermal then uses `/p_amp_refocus`. |
    | `z_beyond_focus_mm` | Domain extent past the focus [mm]. | `40` | |
    | `tight_narrow_beam` | BabelBrain's own lateral domain shrink. | `false` | Always off when `grid.acoustic_fov_*` is set. |
    | `focal_calibration` | Geometric-fallback steering calibration. | `true` | Locate-only when explicit element phases exist. |
    | `focal_calibration_target` | Quantity placed at `focal_distance_ep`. | `'flhm_center'` | `flhm_center` / `peak` |
    | `profile_calibration` | Refit drive phases to the continuous-bowl ideal. | `false` | BabelBrain only; moves it away from the k-Wave/jwave phases. Leave off for cross-backend comparison. |
    | `water_amp_calibration` | Scale output so its free-water peak matches the analytic Rayleigh peak. | `true` | Removes BabelBrain's CFL-correction amplitude bias; requires explicit element phases. |
    | `use_rayleigh_for_water` | Water-only solver. | `false` | `false` = FDTD (like-for-like with k-Wave/jwave). |
    | `q_correction` | Bone Q-factor correction (BabelBrain `QCorrection`) of the relaxation attenuation model. | `[]` | `[]` = 1.0 on the non-CT path (realises the specified bone `alpha_coeff`), BabelBrain's 3.0 on the CT path. Realised bone attenuation scales ~1/Q; values > 0 override both paths. |
    | `ct_mapping` | BabelBrain HU→property mapping. | `''` | Empty = label mask with PRESTUS `medium_properties`. See [CT-based skull properties](doc_babelbrain.md#ct-based-skull-properties-ct_mapping). |
    | `air_regions` | Internal air (sinuses, mastoid) as a pressure-release reflector mask. | `false` | Pseudo-CT path only (`ct_mapping` set); not applied to the free-water run. See [Internal air](doc_babelbrain.md#internal-air-air_regions). |
    | `air_hu_window` | HU range classified as air for `air_regions`. | `[-1200, -400]` | BabelBrain's default. |
    | `air_keep_boundary` | Keep air components touching the domain boundary. | `false` | By default they are removed with the exterior air (cropping can split it). |
    | `ct_use_petra` | PETRA density formula for `Webb-Marsac`. | `false` | |
    | `tx_mech_adjust_z_m` | Additional transducer standoff [m]. | `0` | Added to the automatic placement correction. |
    | `z_steering_m` | Pin BabelBrain `ZSteering` [m]. | `~` | Empty = auto-derived. |
    | `isppa_w_per_cm2` | Override the Isppa given to BabelBrain thermal [W/cm²]. | `~` | Interpreted as the **in-situ peak ISPPA in brain**. Empty = derived from the scaled acoustic field (recommended; matches k-Wave thermal). |
    | `cleanup_cache` | Delete the per-subject BabelBrain cache after reading. | `true` | Result [HDF5](https://www.hdfgroup.org/solutions/hdf5/)s are moved out first. |
    | `dry_run` | Prepare the BabelBrain inputs without simulating. | `false` | The pipeline stops after the acoustic stage; nothing is cached. |

#### jwave Settings (`jwave`)

Used when `acoustic_backend: jwave`.

??? note "jwave settings"

    PML thickness is set by `grid.pml_size`. Full reference: [jwave configuration](doc_jwave.md#configuration-reference).

    | **Parameter** | **Description** | **Default** | **Comments** |
    |---|---|---|---|
    | `python_env` | Python environment or executable with jwave. | `''` | Empty = system `python3`. |
    | `python_env_gpu` | GPU-capable Python environment. | `''` | Empty = same as `python_env`. |
    | `use_gpu` | Use GPU acceleration. | `~` | `true` / `false` / `~` (auto-detect). |
    | `solver` | jwave solver. | `'helmholtz'` | `helmholtz` (frequency-domain GMRES) / `pseudo_spectral` (time-domain, as `kspaceFirstOrder*`) / `finite_differences` (time-domain, FD stencils). |
    | `fd_accuracy` | FD stencil order. | `8` | `finite_differences` only. |
    | `gmres_restart` | GMRES restart parameter (Krylov subspace size). | `10` | jwave's default; larger values can destabilise GMRES in single precision. |
    | `gmres_maxiter` | Maximum GMRES iterations. | `1000` | jwave's default. |
    | `gmres_tol` | GMRES convergence tolerance. | `1.0e-3` | jwave's default. |
    | `convergence_threshold` | Relative Helmholtz residual above which the run errors. | `0.05` | See `ignore_convergence_threshold`. |
    | `ignore_convergence_threshold` | Warn instead of erroring above `convergence_threshold`. | `false` | The unconverged result is then returned. |
    | `helmholtz_operator` | Density term of the Helmholtz solver. | `'conservative'` | `conservative`: ρ∇·((1/ρ)∇p), source scaled by the local sound speed. `jwave`: jwave's built-in term, which underestimates transmission through density jumps (~20 % through a skull plate). See [jwave: density term](doc_jwave.md#density-term-helmholtz_operator). |
    | `precision` | Computational precision. | `'single'` | `single` / `double` |
    | `dry_run` | Write the HDF5 solver inputs only. | `false` | The pipeline stops after the acoustic stage and caches nothing. |
    | `cleanup_cache` | Delete the HDF5 solver input after a successful run. | `true` | The input is rewritten on every run and never reused. Kept on dry runs and failed runs; set `false` to keep it for debugging. |

#### Elastic Skull Settings (`elastic`)

Used when `acoustic_backend: kwave_elastic`.

??? note "Elastic skull settings"

    Applied to the skull layers (`skull`, `skull_cortical`, `skull_trabecular`). Attenuations are power laws a·f^y [dB/(cm·MHz^y)], converted to the solver's y = 2 law so the value at the source frequency is exact. Sources and rationale: [k-Wave elastic backend](doc_kwave.md#elastic-skull-model-shear-waves).

    | **Parameter** | **Description** | **Default** | **Comments** |
    |---|---|---|---|
    | `shear_mapping` | Skull shear speed mapping. | `'pichardo'` | `pichardo`: cS = a(f0)·ρ + b(f0), Pichardo et al. 2017 fits (a = 0.458, b = 545.52 at 270 kHz; a = 0.386, b = 815.51 at 836 kHz), interpolated linearly to the transducer frequency, extrapolated with a warning outside that range. Uses the density map (pseudo-CT or layer densities). `off`: cS = 0 (fluid skull through the elastic solver). |
    | `shear_density_min` | No shear support at or below this density [kg/m³]. | `1200` | |
    | `poisson_min` | Minimum Poisson ratio; caps cS at cL·√((1−2ν)/(2(1−ν))). | `0.1` | Bone: ν ≈ 0.2–0.35. |
    | `alpha_coeff_shear` | Shear attenuation [dB/(cm·MHz^y)]. | `43.63` | Power law through Pichardo et al. 2017 (57 Np/m at 270 kHz, 373 Np/m at 836 kHz). |
    | `alpha_power_shear` | Shear attenuation exponent y. | `1.662` | |
    | `alpha_coeff_long` | Skull longitudinal attenuation in elastic runs [dB/(cm·MHz^y)]. | `27.31` | Power law through the values Pichardo et al. 2017 paired with their shear fit (33 Np/m at 270 kHz, 231 Np/m at 836 kHz). Used with `long_attenuation: elastic`. |
    | `alpha_power_long` | Longitudinal attenuation exponent y. | `1.722` | |
    | `long_attenuation` | Source of skull longitudinal attenuation. | `'elastic'` | `elastic`: `alpha_coeff_long`. `medium`: keep the fluid/pCT skull `alpha_coeff` (effective values that already include shear losses). |
    | `min_ppw_shear` | Warn below this many points per shear wavelength. | `6` | |

---

### Data Paths (`path`)

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `anat` | Absolute path to structural input data. | none | Mandatory |
| `sim` | Absolute path to simulation outputs. | none | Mandatory |
| `seg` | Absolute path to SimNIBS segmentations (m2m folders). | none | Mandatory |
| `placement` | Placement outputs: one folder per run, `<placement>/sub-XXX/<mode>/<target_name>/` (heuristic, [PlanTUS](https://github.com/mlueckel/PlanTUS), [TUS_entry](https://github.com/CyrilAtkinson/TUS_entry)). | `''` | Empty = `<sim>/placement`. Share one folder across setups to reuse placements. See [Where outputs live](doc_outputs.md#where-outputs-live). |
| `localite_raw` | Raw Localite exports, `<localite_raw>/sub-XXX/ses-XX/[<loc_folder>/]…`. | `''` | See [Localite I/O](doc_placement_neuronav_coordinates.md#localite-file-layout). |
| `localite_post` | Ingested Localite markers, `<localite_post>/sub-XXX/ses-XX/InstrumentMarker_sub-XXX_ses-XX.xml/.json`. | `''` | Written by `neuronav_ingest_markers` / `neuronav_gui`; read by `placement.mode: 'localite'` (automatic mode). |
| `localite_t1` | Folder with the planning T1s loaded into Localite (`[sub-XXX/]<localite_t1_pattern>`). | `''` | Empty = `localite_raw`. |
| `localite_t1_pattern` | File-name template of the planning T1; `%s` = `sub-XXX`. | `'%s_*T1*.nii*'` | Replaces `neuronav_t1_pattern` (still read). |
| `localite_export` | Output folder for files to import into Localite (heuristic placement), per subject `sub-XXX/`. | `''` | Resolved into `io.dir_localite_export`. |
| `pct` | Pseudo-CT outputs, per subject `<pct>/sub-NNN/`. | `''` | Empty = `<parent of path.seg>/pct/sub-NNN/`, next to the SimNIBS folder (never inside it). Resolved into `io.dir_pct`. |
| `acoustic_input` | Start from a prior pressure map instead of simulating: a peak-pressure NIfTI [Pa] in this subject's T1w space, e.g. an earlier run's `nii/*_T1w*_pressure.nii.gz`. | `''` | Its JSON sidecar is checked against the current settings and used to rescale the map to `transducer.target_isppa_wcm2`; acoustic analysis and thermal simulation then run as usual. Replaces `io.external_acoustic_nifti` (still read). See [Starting from earlier stage outputs](doc_outputs_reuse.md#starting-from-earlier-stage-outputs). |
| `localite` | **Deprecated.** | `''` | Read as `localite_raw` and as `localite_export` when those are empty. |
| `t1_pattern` | T1 image path template relative to `path.anat`. | `'sub-%1$03d_T1w.nii*'` | Supports `%03d`-style subject ID substitution |
| `t2_pattern` | T2 image path template relative to `path.anat`. | `''` | Optional; empty skips the T2 (e.g. `'sub-%1$03d_T2w.nii*'`). |

---

### Environment & Toolbox Paths (`startup`)

Segmentation and MNI conversion need `startup.simnibs_bin_path`. The other paths are needed for pseudo-CT generation or to add folders to the MATLAB path.

??? note "Environment and toolbox paths"

    | **Parameter** | **Description** | **Comments** |
    |---|---|---|
    | `simnibs_bin_path` | Absolute path to SimNIBS binaries. | Mandatory for segmentation and MNI conversion |
    | `matlab_bin_path` | Absolute path to MATLAB binary. | Required for pseudoCT generation on HPC (e.g. `/opt/matlab/R2024a/bin/matlab`). Local execution falls back to `matlabroot`. |
    | `fsl_bin_path` | Absolute path to [FSL](https://fsl.fmrib.ox.ac.uk/fsl/) bin directory. | Required for pseudoCT generation (e.g. `/opt/fsl/bin`). Prepended to PATH at runtime. |
    | `ants_bin_path` | Absolute path to [ANTs](https://github.com/ANTsX/ANTs) bin directory. | Required for pseudoCT generation (e.g. `/opt/ants/bin`). Prepended to PATH at runtime. |
    | `ants_simg_path` | ANTs [Singularity](https://apptainer.org/) image (`.simg`/`.sif`). | When set, ANTs tools run via `singularity exec`. |
    | `modules_to_load` | HPC environment modules loaded before pseudoCT generation. | e.g. `['matlab/R2023b', 'fsl/6.0.7', 'ants/2.5.3']` |
    | `paths_to_add` | Paths to add with `addpath()`. | [cell] e.g. `{"path/to/x"}` |
    | `subpaths_to_add` | Paths to add recursively with `addpath(genpath())`. | [cell]; relative to the MATLAB working directory |

---

### I/O Management (`io`)

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `output_affix` | Optional affix for output file names. | `''` | Differentiates outputs for the same subject/transducer (e.g. different intensities or targets) |
| `overwrite_files` | File overwrite behaviour. | `'always'` | `never` / `always` / `force` / `ask`. Does NOT apply to SimNIBS segmentations. In chained pipelines (async, multi-ISPPA, uncertainty), `always` re-runs only the output stages (CSV, NIfTI, HTML report) and reuses cached intermediates; `force` re-runs every stage. See [Overwrite behaviour](doc_outputs_reuse.md#overwrite-behaviour). |
| `cache_signature_check` | Response when a reused grid, medium or acoustic cache, or a stage output a run starts from (placement run, pseudo-CT, pressure map), was built from different settings. | `'error'` | `error` / `warn` (reuse anyway) / `off`. Caches: applies when one is reused because of `overwrite_files: never` (or a declined overwrite). See [Cache reuse checks](doc_outputs_reuse.md#cache-reuse-checks) and [Starting from earlier stage outputs](doc_outputs_reuse.md#starting-from-earlier-stage-outputs). |
| `overwrite_simnibs` | Overwrite SimNIBS segmentation results? | `0` | `1 = yes`, `0 = no` |
| `save_matrices` | Save all intermediate matrices (`.mat` / cache files)? | `[]` | `1` = save all, `0` = save none, `[]` = use the per-step defaults. Applies to every per-step matrix flag left empty; a per-step flag set to `0`/`1` always wins. `0` also triggers cleanup of uncertainty intermediates. |
| `save_source_matrices` | Save `kwave_source.mat` (transducer source array)? | `[]` → `1` | `[]` inherits `save_matrices`, else default. Needed by multi-stage pipelines (e.g. uncertainty). |
| `save_acoustic_matrices` | Save acoustic results (`sensor_data`, `kgrid`, medium; can be very large)? | `[]` → `0` | `[]` inherits `save_matrices`, else default. Needed to reload acoustics (e.g. separate thermal runs, calibration debugging). |
| `save_thermal_matrices` | Save `heating_res.mat` (temperature time series and field maps)? | `[]` → `1` | `[]` inherits `save_matrices`, else default. |
| `save_grid_cache` / `save_medium_cache` | Save grid / medium caches? | unset → `0` | Inherit `save_matrices` when unset. |
| `save_property_maps` | Save per-property NIfTIs to nii/properties/ (layered medium only)? | `0` | Not affected by `save_matrices`. |
| `save_MNI` | Save NIfTI outputs in MNI space in addition to native T1w space (requires SimNIBS m2m folder)? | `1` | Not affected by `save_matrices`. |
| `save_heatingvideo` | Save a video of incremental heating? | `0` | `1 = yes`, `0 = no` |
| `save_Ipa` | Save the pulse-average intensity (Ipa) NIfTI map? | `1` | Suppressed automatically when `timing.pd` is absent. |
| `save_Ita` | Save the temporal-average intensity (Ita) NIfTI map? | `1` | Suppressed automatically when `timing.pd`/`timing.pri` are absent. |
| `save_p_complex` | Save the complex pressure field (`p_complex`) in the acoustic cache? | `0` | Large. |
| `deface_plots` | Remove facial anatomy from all QC figures? | `1` | `1 = yes`, `0 = no`. Uses the SimNIBS eye label to anchor the cut plane ~30 mm posterior to the eye socket (ear-canal level). Requires the SimNIBS segmentation (`path.seg`). Has no effect on saved NIfTIs; see `segmentation.deface` for that. |
| `adopted_heatmap` | Path to an existing **temperature** heatmap NIfTI (`heating_end.nii.gz`) to use as the thermal starting point for a follow-up simulation. Set automatically by the sequential dispatcher; can be overridden manually. | none | Optional. Used for sequential multi-target runs. |
| `adopted_cem43` | Path to an existing CEM43 heatmap NIfTI to accumulate heating across runs. Set automatically by the sequential dispatcher. | none | Optional. Used for sequential multi-target runs. |
| `adopted_cem43_iso` | Path to an existing ISO-variant CEM43 heatmap NIfTI (`CEM43_iso_end.nii.gz`) to accumulate across runs. Set automatically alongside `adopted_cem43` by the sequential dispatcher. | none | Optional. Used for sequential multi-target runs. |

??? note "Runtime-derived io fields"

    **Runtime-derived (set by `path_log_setup`):**

    | **Field** | **Description** |
    |---|---|
    | `dir_output` | Resolved per-subject output directory. Set from `path.sim` (+ subject subfolder). |
    | `dir_nii` | NIfTI output subfolder (`nii/`). Created automatically. |
    | `dir_nii_T1w` | Alias for `dir_nii`; space encoded in filename. |
    | `dir_nii_MNI` | Alias for `dir_nii`; space encoded in filename. |
    | `dir_img` | Image/figure output subfolder (`img/`). Used for all acoustic, thermal, and preprocessing figures. Created automatically. |
    | `dir_tabular` | Directory for per-subject CSV results tables (same as `dir_output`). |
    | `dir_reports` | Directory for simulation reports (same as `dir_output`). |
    | `dir_logs` | Log subfolder (`log/`). Created automatically. |
    | `dir_cache` | Cache subfolder (`cache/`). Created automatically. |
    | `dir_debug` | Debug subfolder (`debug/`). Created when `simulation.debug = 1`. |
    | `dir_debug_preproc` | Debug subfolder for preprocessing (`debug/preproc/`). |
    | `dir_debug_medium` | Debug subfolder for medium setup (`debug/medium/`). |
    | `dir_debug_source` | Debug subfolder for source setup (`debug/source/`). |
    | `dir_pct` | Resolved subject pCT directory (`get_pct_dir`): `<path.pct>/sub-NNN`, or `<parent of path.seg>/pct/sub-NNN`. |
    | `dir_localite_export` | Per-subject Localite export folder (`<path.localite_export>/sub-XXX/`, fallback `path.localite`); `''` if unset. |
    | `filename_table` | Path to the per-subject CSV results table. |
    | `filename_kwave_source` | Path to cached k-Wave source `.mat` file. Set by `source_sensor_setup`. |

---

### Pipeline Modules (`modules`)

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `run_transducer_placement` | Resolve `trans_pos`/`focus_pos` from `placement.mode`. | `1` | Required when `placement.manual.coordinate_space` is `ras` or `mni` (positions are converted here). |
| `run_grid_setup` | Setup grid and run head processing. | `1` | Mandatory for simulations. |
| `run_medium_setup` | Map medium acoustic properties. | `1` | Mandatory for simulations. |
| `run_source_setup` | Set up acoustic source. | `1` | Mandatory for simulations. |
| `run_water_baseline` | Run a free-water simulation that records the ISPPA of the current drive, for scaling to `transducer.target_isppa_wcm2`. | `'auto'` | `auto`: when a target is set and the transducer has no calibration ISPPA or the target differs from it. `1`: always. `0`: only when the target differs from a calibration ISPPA. A target that cannot be met stops the run (`prestus:isppa:noBaseline`). See [ISPPA targeting](doc_simulations-acoustic.md#isppa-targeting). |
| `run_acoustic_sims` | Run acoustic simulations. | `1` | |
| `run_acoustic_analysis` | Run acoustic analysis. | `1` | |
| `run_heating_sims` | Run thermal simulations. | `0` | Enable after acoustic results are validated. |
| `run_thermal_analysis` | Run thermal analysis. | `1` | |
| `run_nifti_creation` | Export results as NIfTI files. | `1` | |
| `run_posthoc_water_sims` | Run free-water reference simulations after head simulations. | `1` | |
| `generate_report` | Generate self-contained HTML simulation report. | `1` | |
| `report_style` | Visual skin of the HTML report: `instrument`, `classic` (pre-v0.6.1 look), `compact`, `paper`. Changes presentation only, not content. | `'instrument'` | Try skins on an existing output folder with `preview_report_styles(<subject dir>)`. |
| `report_style_switcher` | Embed all skins plus a dropdown to switch between them in the report. | `0` | |
| `run_targeting_analysis` | Measure targeting accuracy from the focal-mask NIfTIs: simulated focus vs planned focus (native), vs MNI target / ROI (MNI). Writes `sub-XXX_<medium><affix>_targeting.csv`. | `1` | Layered or phantom media; needs `run_nifti_creation`. See [Targeting accuracy](doc_group.md#targeting-accuracy). |
| `generate_group_report` | Aggregate completed subjects into a group HTML report at the end of the run. | `1` | Best-effort; auto-discovers subjects. Or run `prestus_group_report_start` manually. Switched off for runs started by `prestus_study_start`, which builds one study report at the end. |
| `segmentation_only` | Stop after segmentation; skip grid setup and all simulations. | `0` | Only has effect when `simulation.medium = 'layered'`. |

All flags: `1 = yes`, `0 = no`.

---

### Transducer Specification (`transducer`)

See [doc_transducer.md](doc_transducer.md).

Except `target_isppa_wcm2` (default `30`), the fields have no defaults and must be set in the study config or come from the equipment library (`serial`).

| **Parameter** | **Description** | **Comments** |
|---|---|---|
| `target_isppa_wcm2` | Target free-water ISPPA [W/cm²] of this transducer alone. | Default `30`. Scalar: the field is scaled to it (needs a free-water baseline; see `modules.run_water_baseline`). Vector: one thermal job per value ([multi-intensity](doc_multiintensity.md)). Empty `[]`: no scaling, the field is reported at `elem_p_amp`. Several transducers: each transducer's own target; coherent runs set each drive before the simulation, async runs scale each field before combining (see [Multi-transducer: ISPPA targets](doc_multitransducer.md#isppa-targets)). |
| `type` | Type of transducer array. | `annular` / `matrix`. Mandatory. |
| `freq_hz` | Fundamental frequency [Hz]. | Shared across all elements. |
| `serial` | Equipment serial; geometry, phases and drive come from the equipment config and transducer library. | See [doc_transducer_library.md](doc_transducer_library.md). Values set in the study config take precedence. |
| `combo.ds_serial` | Driving-system serial. | Selects `<serial>_<ds_serial>.yaml`; without it only `<serial>.yaml` is looked up. |
| `library_path` | Folder with your own library entries. | Used when it holds an entry for this combo, otherwise the shipped `config/calibration/` (with a warning). |
| `condition` | **Experimental, not validated.** Dephased library condition (`conditions.<label>` in the library entry). | Replaces the focal-depth lookup; `focal_distance_ep` defaults to the condition's nominal setpoint. See [Dephased conditions](doc_transducer_library.md#dephased-conditions). |
| `trans_pos` | Transducer position (XYZ, bowl apex). | T1 voxels by default; RAS+ or MNI mm via `placement.manual.coordinate_space`. |
| `focus_pos` | Focus (target) position (XYZ). | T1 voxels by default; RAS+ or MNI mm via `placement.manual.coordinate_space`. |
| `focal_distance_ep` | Focal depth setting: the depth the device is set to focus at, from the exit plane [mm]. | For annular arrays it selects the library phases and drive, so the simulated focus lands here whatever the `trans_pos`→`focus_pos` distance (the planned target depth). Empty = that distance minus `focal_distance_offset` (`focal_distance_calculation`). |
| `focal_distance_bowl` | Focal depth setting from the bowl apex [mm]. | `focal_distance_ep` + bowl-to-exit-plane offset; set either. |
| `focal_distance_offset` | Offset between transducer bowl and exit plane [mm]. | **Derived** from `curv_radius_mm − dist_geom_ep_mm`. Not set by user. |

#### Annular Array Definition (`annular`)

| **Parameter** | **Description** | **Comments** |
|---|---|---|
| `elem_p_amp` | Pressure amplitude [Pa]; scalar or one value per element. | Must be calibrated. `elem_amp` is the deprecated name (renamed on load with a warning). |
| `elem_phase_deg` | Source phase per element [degrees]. | Must be calibrated. |
| `elem_n` | Number of transducer elements. | |
| `elem_id_mm` | Inner diameter of each element [mm]. | |
| `elem_od_mm` | Outer diameter of each element [mm]. | |
| `curv_radius_mm` | Radius of curvature of the transducer bowl [mm]. | `Inf` (YAML `.inf`) = flat element, e.g. a plane piston: its face at `trans_pos` is the exit plane; needs `grid.use_kWaveArray: 1`, in 3D a single disc. See [Flat elements](doc_transducer.md#flat-elements). |
| `dist_geom_ep_mm` | Distance from geometric focus to transducer plane [mm]. | Calculated automatically from `curv_radius_mm` if not provided (`Inf` for flat elements). |
| `depth_mm` | Transducer depth [mm]. | Visualization only. |

#### Matrix Array Definition (`matrix`)

| **Parameter** | **Description** | **Comments** |
|---|---|---|
| `elem_p_amp` | Pressure amplitude [Pa]; scalar or one value per element. | Must be calibrated. `elem_amp` is the deprecated name (renamed on load with a warning). |
| `depth_mm` | Transducer depth [mm]. | Visualization only. |
| `steering` | Steering mode of the transducer. | `1D` = axial only, `3D` = volumetric steering. |
| `elem_shape` | Shape of individual elements. | Options: `rect`, `disc`, `bowl`. Element area is defined using the rectangular dimensions and projected onto the selected shape, unless `elem_diameter_mm` is set. |
| `elem_height_mm` | Element height [mm], along y. | `rect` elements: size along y; `disc`/`bowl`: sets the equal-area diameter. |
| `elem_width_mm` | Element width [mm], along x. | `rect` elements: size along x; `disc`/`bowl`: sets the equal-area diameter. |
| `elem_diameter_mm` | Aperture diameter of `disc`/`bowl` elements [mm]. | Optional. Empty = equal-area diameter from `elem_width_mm` × `elem_height_mm` (flat disc area for `disc`, spherical-cap area for `bowl`). Also sets the element size in the transducer mask. |
| `elem_phase_override_deg` | Calibrated per-element phases [deg]. | Optional, one value per element (after Clover replication). Replaces the geometric steering phases toward `focus_pos`; element orientation is unchanged. Empty = geometric steering. |
| `use_steering_correction` | Scale the element drive by the steering correction k(target). | Default `true`. Used when the transducer library (or `steering_correction`) holds one; see [Steering correction of matrix arrays](doc_calibration_matrix-steering.md#steering-correction-of-matrix-arrays). |
| `steering_correction` | Steering correction given inline. | Optional struct as written by `build_steering_correction`. Empty = from the transducer library. |
| `outer_diameter_mm` | Outer diameter of the transducer [mm]. | Defines active aperture boundary. |
| `is_curved` | Whether the array has a curved surface. | `false` = flat array: every element faces +z, `curv_radius_mm` is not used, `elem_shape` must be `rect` or `disc`. |
| `curv_radius_mm` | Radius of curvature (ROC) of the transducer bowl [mm]. | Defines natural focus. Depends on `is_curved`. |
| `dist_geom_ep_mm` | Distance from geometric focus to transducer plane [mm]. | Calculated automatically from `curv_radius_mm` and `outer_diameter_mm` if not provided. Depends on `is_curved`. |
| `is_clover_setup` | Enables Clover (multi-array) configuration. | Replicates array into multiple leaves. |
| `matrix_shape.type` | Method used to define element positions. | Options: `define_here`, `extract_from_file`. |

##### Matrix: Clover (`clover`)

??? note "Clover keys"

    | **Parameter** | **Description** | **Comments** |
    |---|---|---|
    | `n_leaves` | Number of Clover leaves. | Any positive integer; one leaf stays on the parent axis unless `elevation_deg` is set. |
    | `ROC_parent` | Radius of curvature of the combined Clover setup [mm]. | Can differ from individual array ROC; elements then point at their own leaf's natural focus. |
    | `az_sep_deg` | Azimuthal separation between adjacent leaves [deg]. | Optional. Empty = `360 / n_leaves`. `n_leaves × az_sep_deg` must not exceed 360. |
    | `elevation_deg` | Tilt of each leaf from the parent axis [deg]. | Optional. Empty = smallest tilt at which adjacent leaf centres are one `outer_diameter_mm` apart on the parent bowl. |

##### Matrix: Grid Distribution (`define_here`)

??? note "Grid distribution keys"

    Warning: The `define_here` options are experimental and may contain bugs due to limited testing time.

    | **Parameter** | **Description** | **Comments** |
    |---|---|---|
    | `grid_shape.type` | Grid distribution type. | Options: `rect`, `fibonacci`, `fermat`. |

    **Rectangular Grid**

    | **Parameter** | **Description** | **Comments** |
    |---|---|---|
    | `elem_n_row` | Number of element rows (along y). | Row pitch: `elem_height_mm + elem_spacing_height_mm`. |
    | `elem_n_col` | Number of element columns (along x). | Column pitch: `elem_width_mm + elem_spacing_width_mm`. Elements whose centre lies outside `outer_diameter_mm` are dropped. |
    | `elem_spacing_height_mm` | Spacing between elements in height direction [mm]. | Edge-to-edge spacing. |
    | `elem_spacing_width_mm` | Spacing between elements in width direction [mm]. | Edge-to-edge spacing. |
    | `sparsity_factor` | Fraction of active elements. | Range: 0.1–1.0. Validated but **not applied** yet: all elements in the grid are used. |

    **Fibonacci Grid**

    | **Parameter** | **Description** | **Comments** |
    |---|---|---|
    | `elem_n` | Total number of elements. | |
    | `kerf_mm` | Minimum spacing between elements [mm]. | |

    **Fermat Grid**

    | **Parameter** | **Description** | **Comments** |
    |---|---|---|
    | `elem_n` | Total number of elements. | Uses spiral distribution. |

##### Matrix: File Extraction (`extract_from_file`)

??? note "File extraction keys"

    | **Parameter** | **Description** | **Comments** |
    |---|---|---|
    | `file_path` | Path to coordinate file. | Must contain (x, y, z) [mm]. Curved array: origin at the centre of curvature, z towards the array. Flat array: origin at the array centre, z along the beam. |
    | `start_row` | Row index where data starts. | MATLAB may skip header row. |
    | `start_col` | Column index where data starts. | |
    | `rot_col` | Column with each element's rotation about its own normal [deg]. | Optional. For rectangular elements laid out in different orientations (e.g. `90` = turned on its side). Empty = all elements as `elem_width_mm` × `elem_height_mm` along x × y. |
    | `elem_n` | Number of elements to extract. | Can exceed final active count. |
    | `select_random_subset` | Enable random subset selection. | Useful for sparse arrays. |
    | `subset.random_seed` | Controls reproducibility of subset. | `true` = new subset each run. |
    | `subset.subset_n_elements` | Number of elements in subset. | Must be ≤ total elements. |
    | `project_on_new_ROC` | Project elements onto new curvature. | May introduce unrealistic layouts. |
    | `ROC_projection.new_ROC_mm` | New radius of curvature [mm]. | Used when projection enabled. |

### Transducer Placement (`placement`)

See [doc_placement.md](doc_placement.md).

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `mode` | Placement mode. | `'manual'` | `'manual'`, `'localite'`, `'heuristic'`, `'mni'`, `'plantus'` or `'tus_entry'`. |
| `report` | Write the stand-alone placement report `sub-XXX_placement_<mode>.html` to `io.dir_reports`. | `true` | See [placement report](doc_placement.md#placement-report). |
| `skin_gap_mm` | Gap between scalp and transducer exit plane along the beam axis [mm]. | `5` | e.g. a gel pad. Used by heuristic, MNI, PlanTUS and TUS_entry placement and by their [time-reversal pose search](doc_time_reversal.md); manual and Localite poses keep their own gap. The exit plane lies `curv_radius_mm − dist_geom_ep_mm` in front of the bowl apex. |
| `reuse` | When to reuse a recorded heuristic, PlanTUS or TUS_entry run (`placement_run.json` in its run folder). | `'auto'` | `'auto'`: reuse it when its settings match, else run again; `'always'`: require it (mismatch per `io.cache_signature_check`); `'never'`: always run. Independent of `io.overwrite_files`. See [Starting from earlier stage outputs](doc_outputs_reuse.md#starting-from-earlier-stage-outputs). |

#### Manual Placement (`placement.manual`)

See [doc_placement.md](doc_placement.md#manual-coordinate-selection).

??? note "Manual placement keys"

    | **Parameter** | **Description** | **Default** | **Comments** |
    |---|---|---|---|
    | `coordinate_space` | Space of `transducer.trans_pos` / `focus_pos` in `mode='manual'`. | `'voxel'` | `'voxel'` (1-based T1 voxels), `'ras'` (subject RAS+ mm), `'mni'` (MNI mm; transducer via 12dof, focus via nonlinear warp, no scalp snapping). Converted per transducer in the placement stage, which must be enabled. |

#### Localite Placement (`placement.localite`)

See [doc_placement_neuronav.md](doc_placement_neuronav.md).

??? note "Localite placement keys"

    | **Parameter** | **Description** | **Default** | **Comments** |
    |---|---|---|---|
    | `enabled` | Load transducer position from Localite neuronavigation files in head preprocessing? | `0` | `1 = yes`, `0 = no`. Older switch; `placement.mode: 'localite'` is the placement-stage equivalent. Both use `neuronav_localite_pose`. |
    | `file` | Single Localite XML (TriggerMarkers, GUMMarkers or InstrumentMarker; type detected). | `''` | File mode. Empty = automatic mode (ingested markers under `path.localite_post`). |
    | `session` | Session of the automatic mode. | `1` | Integer `n` → `ses-0n`, or `'ses-XX'`. |
    | `target` | Select the series by name. | `''` | Ingested target name (automatic mode) or entry description (file mode). Several transducers: one name per transducer, e.g. `[L_M1, R_M1]` (transducer *i* uses entry *i*); a single name applies to all. |
    | `position` | Otherwise: index of the detected series (1-based, after the `expected_segment_length` filter). | `1` | Automatic mode matches it against `series_index` in the JSON sidecar. Several transducers: one index per transducer, e.g. `[1, 2]`. |
    | `markertype` | **Deprecated, ignored.** | `'TriggerMarkers'` | The marker type is detected from the files. |
    | `loc_folder` | Folder between `ses-XX/` and the Localite export. | `'loc'` | `''` = the export sits directly under `ses-XX/`. |
    | `session_folder_pattern` | Wildcard that a folder between `loc_folder` and the marker files must match. | `'*'` | e.g. `'Session_2024*'`. |
    | `expected_segment_length` | Pulses per TriggerMarkers series. | `80` | Series off by more than 20 pulses are dropped; `[]` keeps all. |
    | `voxel_size` | Voxel size for `position_dev_vox` in the series statistics [mm]. | `0.9` | Does not affect positions. |
    | `reference_distance_mm` | Fallback offset from the Localite matrix origin to the bowl rear centre, along the beam axis [mm]. | `15` | Used only when `tracker_to_bowl_mm` is empty and the transducer has no scalar `curv_radius_mm` and `dist_geom_ep_mm`. Same sign convention as `tracker_to_bowl_mm`. |
    | `tracker_to_bowl_mm` | Measured offset from the Localite matrix origin to the bowl rear centre [mm]. | `[]` | Overrides geometry and `reference_distance_mm`. |

#### Neuronavigation Ingestion (`neuronav`)

See [Localite ingestion](doc_placement_neuronav_workflows.md#shared-settings-file); `coil_map` also feeds [target resolution](doc_placement_neuronav_targets.md#target-coordinates-from-several-sources).

??? note "Neuronav ingestion keys"

    | **Parameter** | **Description** | **Default** | **Comments** |
    |---|---|---|---|
    | `coil_map` | Multi-coil ingestion: one entry per coil with `.coil` (index of `TriggerMarkers_Coil<N>`) and `.target` or `.sequences(:).target`. | `[]` | `[]` = single-coil ingestion (Coil0, GUMMarkers fallback). |
    | `suggest_names` | Name series after the planned target nearest their beam axis. | `true` | See [naming](doc_placement_neuronav_workflows.md#series-naming). |
    | `axis_tolerance_mm` | Max. distance of that target from the beam axis [mm]. | `10` | |
    | `target_depth_mm` | Accepted target depth along the beam axis [mm]. | `[0, 150]` | Narrow around the focal distance to reject far matches. |
    | `target_labels` | Label map: entries `.target` + `.labels` (patterns with `*`/`?`, case-insensitive). | `[]` | `[]` = labels used as recorded; list them with `neuronav_scan_sessions`. |
    | `batch_write` | Series written by `neuronav_ingest_batch`: `'all'` tracked series or only `'named'` ones. | `'all'` | |

#### Heuristic Placement (`placement.heuristic`)

See [doc_placement_heuristic.md](doc_placement_heuristic.md).

??? note "Heuristic placement keys"

    | **Parameter** | **Description** | **Default** | **Comments** |
    |---|---|---|---|
    | `save_localite_t1` | Save Localite-aligned T1 output for header correction? | `false` | See [doc_placement_heuristic.md](doc_placement_heuristic.md). |
    | `dist_close` | Distance from target considered sufficiently close [mm]. | `[]` | |
    | `ear_radius` | Radius of the ear exclusion zone [mm]. | `35` | |
    | `left_ear_center` | Approximate coordinates for left ear [image voxels]. | `[]` | Optional |
    | `right_ear_center` | Approximate coordinates for right ear [image voxels]. | `[]` | Optional |
    | `criterion_intersection` | Max fraction of exit plane intersecting skin. | `0.05` | 5% |
    | `criterion_skin_mean` | Mean distance of EP voxels from skin [quantile]. | `[]` | |
    | `criterion_skull_mean` | Mean distance of EP voxels from skull [quantile]. | `[]` | |
    | `criterion_skin_var` | Variance of EP voxel distance from skin [quantile]. | `[]` | |
    | `criterion_skull_var` | Variance of EP voxel distance from skull [quantile]. | `[]` | |
    | `expand_step` | Expansion step for skin intersection criterion. | `0.01` | 1% |

#### MNI Placement (`placement.mni`)

See [doc_placement_mni.md](doc_placement_mni.md).

??? note "MNI placement keys"

    | **Parameter** | **Description** | **Default** | **Comments** |
    |---|---|---|---|
    | `trans_pos_mm` | Desired transducer (scalp entry) position [MNI mm]. | `[]` | Required for `mode='mni'`. Converted with the linear (12dof) transform, then snapped to scalp. |
    | `focus_pos_mm` | Focus / target position [MNI mm]. | `[]` | Required for `mode='mni'`. Converted with the nonlinear (nonl) transform, like the heuristic target. |

#### PlanTUS Placement (`placement.plantus`)

See [doc_placement_plantus.md](doc_placement_plantus.md).

??? note "PlanTUS placement keys"

    Headless or interactive selection is set by `interactive`, which defaults to `simulation.interactive`.

    | **Parameter** | **Description** | **Default** | **Comments** |
    |---|---|---|---|
    | `script_path` | PlanTUS installation directory (contains `code/PlanTUS.py`). | `''` | Required. |
    | `env_path` | SimNIBS installation or conda env whose Python runs PlanTUS. | `''` | Required. |
    | `target_name` | Target label, used in output folder and file names. | `''` | Required. Outputs go to `<path.placement>/sub-XXX/plantus/<target_name>/`. |
    | `mni_target_mm` | Target [MNI mm]. | `[]` | Empty = use `transducer.focus_pos` (voxels). |
    | `focal_distance_list` | Focal distances from the calibration [mm]. | `[]` | Required. Set the eligible and optimal skin–target distance. |
    | `flhm_list` | Focal length at half maximum per entry of `focal_distance_list` [mm]. | `[]` | Required. |
    | `max_angle_deg` | Maximum transducer tilt [deg]. | `10` | |
    | `mask_radius_mm` | Radius of the spherical target mask [mm]. | `2` | |
    | `min_distance_mm` | Smallest eligible skin–target distance [mm]. | `[]` | `[]` = `min(focal_distance_list) − max(flhm_list)`. For a steerable array, set it and `max_distance_mm` to the steering range; the default margin otherwise admits spots the focus cannot reach. |
    | `max_distance_mm` | Largest eligible skin–target distance [mm]. | `[]` | `[]` = `max(focal_distance_list) + max(flhm_list)`. |
    | `optimal_distance_mm` | Optimum of the distance criterion [mm]. | `[]` | `[]` = `mean(focal_distance_list)`. |
    | `weight_skin_target_distances` | Weight: deviation of skin–target distance from the optimal distance. | `0.2` | The five weights should sum to 1. See [composite scoring](doc_placement_plantus.md#5a-combine-and-evaluate-scores). |
    | `weight_skin_target_angles` | Weight: angle between skin normal and direction to target. | `0.2` | |
    | `weight_skin_target_intersections` | Weight: beam path length inside the target (higher is better). | `0.2` | |
    | `weight_skin_skull_angles` | Weight: angle between skin and skull normals. | `0.2` | |
    | `weight_skull_thickness` | Weight: bone path length along the skin normal. | `0.2` | From `m2m_sub-XXX/final_tissues.nii.gz`; ignored if that file is missing. See [skull thickness](doc_placement_plantus.md#4e-skull-thickness). |
    | `save_metrics` | Write per-vertex criteria, score and selected vertex to `plantus_metrics.mat`. | `true` | Both modes. Used by the criteria evaluation and offline re-scoring. |
    | `evaluate` | Evaluate the placement from `plantus_metrics.mat` (weight sweep, diagnostics, figures). | `true` | Both modes; needs `save_metrics`. Shown in the placement report. See [criteria evaluation](doc_placement_plantus.md#5a-combine-and-evaluate-scores). |
    | `interactive` | Choose the vertex in [Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench). | `[]` | `1` = interactive, `0` = headless, `[]` = follow `simulation.interactive`. `1` needs a display but not the MATLAB desktop; without a display PlanTUS runs headless. |
    | `connectome_wb_path` | Connectome Workbench `bin/` directory. | `''` | Needed unless `wb_command` is on PATH or in the env's `bin/`; interactive mode also needs `wb_view`. |

#### TUS_entry Placement (`placement.tus_entry`)

See [TUS_entry placement](doc_placement_tus_entry.md).

??? note "TUS_entry placement keys"

    Not in `config_default.yaml`; the keys are listed in `nexus/prestus_nexus/overlay.yaml`.

    | **Parameter** | **Description** | **Default** | **Comments** |
    |---|---|---|---|
    | `script_path` | Full path to `TUS_entry.R`. | none | Required. |
    | `env_path` | Conda environment containing `Rscript`. | none | Required. |
    | `scalp_mask` | Head or scalp mask in the T1 grid (nonzero = head). | SimNIBS `final_tissues.nii.gz` | Filled into a solid head mask in the work directory. |
    | `target_mask` | Binary NIfTI target mask. | *(auto)* | Empty: a sphere around `focus_pos` or `mni_target_mm`. |
    | `target_mask_radius_mm` | Radius of the auto-generated target sphere [mm]. | `5` | |
    | `mni_target_mm` | Target [MNI mm]. | *(from `focus_pos`)* | Converted to subject space with SimNIBS. |
    | `exclu_mask` | Binary NIfTI exclusion zone (e.g. ventricles, vessels). | `''` | |
    | `minimal_distance_cm` | Minimum scalp-to-target depth [cm]. | `4.24` | |
    | `maximal_distance_cm` | Maximum scalp-to-target depth [cm]. | `7.46` | |
    | `transducer_size_cm` | Transducer outer diameter [cm]. | `6.4` | |
    | `voxel_size_cm` | Voxel size [cm]. | *(T1 voxel size)* | TUS_entry assumes isotropic voxels. |
    | `resolution_cm` | Spatial decimation step [cm]. | `0.33` | |
    | `return_candidates` | Save the ranked candidate list as RDS. | `false` | |
    | `target_name` | Label for output folder and file names. | `'target'` | |
    | `work_dir` | Output directory override. | *(auto)* | Default `<path.placement>/sub-XXX/tus_entry/<target_name>/`. |

---

### Simulation Grid (`grid`)

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `resolution_mm` | Grid resolution (must be isotropic) [mm]. | `0.5` | |
| `default_dims` | Requested grid dimensions [voxels per dimension], or `'auto'`. | `'auto'` | Directly sets the simulation grid for `water` and `phantom` media (a phantom uses its own size unless it matches). For `layered`, the grid is determined by head preprocessing and this value is not used. `'auto'` (water only, `water_grid_dims`) sizes the grid to the transducer, as the acoustic FOV of a head simulation: aperture + 2 × `water_margin_mm` (default 5 mm) laterally, focal distance from the bowl + `z_beyond_focus_mm` (default 40 mm) axially; `[Nlateral, Naxial]` with `axisymmetric`. For the shipped transducers this is 15–80 % of the points of `[144, 144, 400]` (the default until 2026-10) with the same free-water ISPPA (within 0.003 %) and focus (k-Wave drive benchmark). A non-annular (matrix) transducer placed by `trans_pos`/`focus_pos`, which are grid indices, keeps `[144, 144, 400]` (`[144, 400]` axisymmetric) with a warning. |
| `axisymmetric` | Run axisymmetric 2D simulation (`kspaceFirstOrderAS`). | `0` | `1 = yes`, `0 = no`. See [doc_simulations-acoustic.md](doc_simulations-acoustic.md). |
| `axisymmetric_output` | Form of axisymmetric results. | `'auto'` | `'auto'`: expanded to 3D for `phantom` media or when heating is requested, otherwise the mirrored 2D plane. `'2d'`: always the plane (not with heating); for fine grids, where the 3D expansion would not fit in memory. `'3d'`: always expanded. |
| `mode` | Grid orientation before cropping. | `'focal_axis'` | `focal_axis`: rotate so the mean beam axis of all transducers (beams averaged as axes) is z. `transducer_axis`: focal axis of transducer 1. `ras_plus`: no rotation, only rescaling. (`principal_axis` was removed and errors.) See [Head preprocessing](doc_preproc.md#orientation). |
| `focal_axis_min_norm` | Fallback threshold of `focal_axis` mode. | `0.5` | Not in `config_default.yaml`. When the mean beam-axis vector is shorter than this (beam axes too spread; only possible with three or more transducers), transducer 1's axis is used. |
| `pml_size` | Perfectly Matched Layer (PML) size [voxels]. | `'auto'` | `'auto'` chooses a size that keeps FFT grid sizes small-prime (k-Wave `getOptimalPMLSize`); a number fixes it. Absorbs waves at boundaries. Recommended for 3D. See [k-Wave docs](http://www.k-wave.org/documentation/example_na_controlling_the_pml.php). |
| `water_baseline_grid` | Grid of the k-Wave free-water baseline (ISPPA targeting). | `'beam'` | `'beam'`: the box around the source and each beam axis from the bowl to 40 mm beyond the deeper of focus and geometric focus, +10 mm (WATER_BASELINE_BOX); in water nothing outside reflects back, so the free-field result is the same at a fraction of the cost. `'full'`: the whole simulation grid. jwave and axisymmetric runs use the full grid. |
| `pml_size_range` | Search range `[min, max]` of `pml_size: 'auto'` [grid points]. | `[10, 40]` | Passed to k-Wave's `getOptimalPMLSize`, which picks the FFT-friendliest size in the range. A thicker layer absorbs better at grazing incidence and costs grid points. |
| `pml_alpha` | PML absorption [Nepers per grid point]. | `[]` | Empty: k-Wave's default (2). Passed as k-Wave's `PMLAlpha`. |
| `source_ppw` | Points per wavelength (PPW) override. | `[]` | Calculated internally from `resolution_mm` and the medium's maximum sound speed if not set. Setting this overrides the internal calculation and fixes the number of spatial samples per wavelength used to derive the time step. |
| `min_ppw` | Minimum acceptable PPW at the fundamental frequency. | `6` | Checked against the slowest medium sound speed (worst-case spatial sampling). A warning is raised if the actual PPW falls below this value, along with the maximum `resolution_mm` that would satisfy it. |
| `source_cfl` | Courant-Friedrichs-Lewy (CFL) fraction. | `0.15` | Controls temporal resolution independently of spatial resolution. Given a fixed spatial step `dx`, the time step is `dt = CFL × dx / c_max`, where `c_max` is the fastest medium sound speed (typically skull at ~2800 m/s). Smaller CFL = finer temporal sampling and must remain below the stability limit. Set to 0.15 (half k-Wave's default of 0.3) for additional stability margin in heterogeneous skull simulations. Can be relaxed toward 0.3 to reduce runtime, but should be validated against `checkStability` output. |
| `source_limit_fraction` | Fraction of the stability limit to use for time step. | `0.9` | After source setup, `checkStability` is called; if the CFL-derived `dt` exceeds the estimated limit, source setup reruns with `dt = source_limit_fraction × dt_stability_limit`. The final effective `dt` may therefore be smaller than `CFL × dx / c_max`. Set to `0` to disable. |
| `steady_state_transits` | Number of domain-diagonal one-way transits (at min. sound speed) before recording starts. | `1` | Total simulation duration is `transits × diagonal / c_min`. Compact, highly reflective domains (e.g. a partial-FOV crop with a high skull/water impedance mismatch, ~50% reflection coefficient for typical skull properties) may need more than one transit for multiply-reflected energy to decay to true periodic steady state. jwave's own cross-validation test uses `roundtrips=20`. |
| `steady_state_num_periods` | Number of trailing wave periods recorded as the steady-state window. | `3` | `sensor.record_start_index` is set so the last `num_periods` periods of the simulation are recorded for `p_max_all`/`p_min_all`/`p_final`. Increase together with `steady_state_transits` if reported Isppa/focal depth are still drifting as either is increased; a converged setup should stop changing further. |
| `bli_tolerance` | kWaveArray band-limited interpolation (BLI) tolerance. | `[]` | Empty: 0.1 (3D), 0.01 (axisymmetric). Smaller values keep more of the sinc tails (more source grid points). |
| `bli_upsampling_rate` | kWaveArray integration points per grid point. | `[]` | Empty: 10 (3D), 100 (axisymmetric). |
| `source_ppp` | Time steps per wave period (PPP) override. | `[]` | When set, `dt = period / source_ppp`; `source_ppw` and `source_cfl` are not used. For matching a reference simulation's time step (e.g. the [ITRUSST](https://itrusst.org/) k-Wave references). The stability recheck (`source_limit_fraction`) still applies. |
| `sim_time_us` | Fixed simulation duration [µs]. | `[]` | When set, replaces `steady_state_transits` (and jwave's PML margin); the last `steady_state_num_periods` periods of this duration are recorded. |
| `max_expand` | Maximum grid expansion for prime-number FFT optimisation [voxels]. | `40` | |
| `use_kWaveArray` | Use the kWaveArray class for transducer modelling? | `1` | `1 = yes`, `0 = no` |
| `acoustic_fov_diameter_mm` | Lateral extent of the acoustic simulation FOV [mm], centred on the beam axis. When set together with `acoustic_fov_length_mm`, the acoustic simulation runs on a cropped sub-grid rather than the full head grid. Thermal simulation runs on the full head grid unless `thermal_use_acoustic_fov: 1`. | `[]` | Optional. `'auto'`: the outer element diameter of `transducer(1)`. See [doc_limited_fov.md](doc_limited_fov.md). |
| `acoustic_fov_length_mm` | Axial extent of the acoustic simulation FOV [mm], starting from the bowl apex (`trans_pos`) and extending distally. Must be large enough to enclose the focus and any near-field region of interest. | `[]` | Optional. Requires `acoustic_fov_diameter_mm`. `'auto'`: the `trans_pos`-to-`focus_pos` distance + `z_beyond_focus_mm`. See [doc_limited_fov.md](doc_limited_fov.md). |
| `thermal_use_acoustic_fov` | Run the thermal simulation within the acoustic FOV instead of the full grid. | `0` | `1` is faster, with a less accurate cooling boundary. See [doc_limited_fov.md](doc_limited_fov.md). |
| `thermal_resolution_mm` | Grid resolution for the thermal simulation [mm]. When set to a value larger than `resolution_mm`, the acoustic pressure field is resampled to a coarser grid before thermal simulation, reducing memory and runtime. | `[]` | Optional. Defaults to `resolution_mm` (same grid). |
| `thermal_fov_mm` | Spatial extent of the thermal simulation domain [mm], as a `[1×3]` vector `[x y z]`. When set, the thermal grid is cropped to a box of this size centred on the focus. Unlike the acoustic FOV crop, the thermal FOV has no physical justification for excluding regions: heat diffuses well beyond the focal zone, and cropping can underestimate temperature rise in neighbouring tissue. Use only for rapid prototyping, not for safety-critical results. | `[]` | Optional. See [doc_limited_fov.md](doc_limited_fov.md). |

---

### Segmentation (`segmentation`)

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `use_qform` | Force qform reorientation before charm segmentation? | `0` | Set to `1` if charm reports a qform/sform mismatch error. |
| `debug` | Pass `--debug` to charm for verbose segmentation output? | `0` | `1 = yes`, `0 = no` |
| `deface` | Create a defaced copy of the SimNIBS m2m folder (`m2m_sub-NNN_defaced/`) for safe sharing? | `0` | `1 = yes`, `0 = no`. Runs [`pydeface`](https://github.com/poldracklab/pydeface) (must be installed in `startup.simnibs_bin_path`) on T1, UTE, bias-corrected, and MNI-space anatomical images. The original m2m folder is never modified. Only affects saved files; see `io.deface_plots` for figure defacing. |

---

### Head Model & Segmentation Preprocessing (`headmodel`)

See [doc_preproc.md](doc_preproc.md).

??? note "Head model keys"

    | **Parameter** | **Description** | **Default** | **Comments** |
    |---|---|---|---|
    | `head_pad_mm` | Symmetric padding of the rotated head volume before CSF expansion and cropping [mm]; bowls outside the volume are auto-padded. | `0` | |
    | `csf_expansion` | Dilation of the CSF brain mask into surrounding head regions [mm]. | `40` | Converted to voxels at the grid resolution (`preproc_crop_eCSF`). |
    | `smooth_method` | Mask smoothing filter type. | `'gaussian'` | `gaussian` / `box` / `off` |
    | `smooth_threshold_skull` | Binarisation threshold for skull mask; higher = thinner mask. | `0.5` | |
    | `smooth_threshold_other` | Binarisation threshold for other masks; higher = thinner mask. | `0.5` | |
    | `smooth_fwhm_mm` | FWHM of smoothing kernel [mm]. | `1` | |
    | `smooth_properties` | Apply smoothing to acoustic property maps as well? | `false` | |
    | `skull_fill_method` | Method for filling holes in the skull. | `'rubberwrap'` | `rubberwrap` / `imclose` |
    | `skull_wrap_radius` | Rubber-wrap radius [grid voxels]. | `10` | Larger = tighter wrap ignoring bigger dents. Recommended: 2–10. |
    | `skull_wrap_visualize` | Visualize rubber-wrap result? | `0` | `1 = yes`, `0 = no`. Avoid on HPC. |
    | `reuse_across_placements` | Build the smoothed head model once per subject and only rotate and crop it for each transducer placement. | `false` | Saves the smoothing and hole filling (most of head preprocessing at fine resolution) in every run after the first of a subject: other poses, aims or affixes in the same output folder. Cached as `<dir_cache>/sub-XXX_<medium>_headmodel_<res>mm.mat` with the settings it depends on; rebuilt when they change. With `grid.mode: ras_plus` the result equals per-placement preprocessing exactly (checked on one head, two poses). With a rotated mode the smoothing precedes the rotation and the labels are resampled once more (nearest neighbour), which changes thin skull layers (skull-layer Dice 0.80–0.90 against per-placement preprocessing at 1 mm; a warning says so). Needs the planning T1 and the segmentation on one voxel grid. |

---

### Pseudo-CT Skull Property Mapping (`pct`)

See [doc_pseudoCT.md](doc_pseudoCT.md).

??? note "Pseudo-CT mapping keys"

    | **Parameter** | **Description** | **Default** | **Comments** |
    |---|---|---|---|
    | `enabled` | Use (pseudo-)CT to inform skull medium properties? | `0` | `1 = yes`, `0 = no`. When enabled, pseudoCT is generated automatically from the UTE image (specified via `path.t2_pattern`) if not already present in the SimNIBS folder. |
    | `skull_mapping` | UTE→HU linear mapping algorithm used during pseudoCT generation. | `'kosciessa'` | `kosciessa` / `miscouridou` / `wiesinger` / `treeby` / `custom`; `carpino` is deprecated. See [doc](doc_pseudoCT_generation.md#choosing-the-ute-to-hu-mapping). |
    | `custom_slope` | Slope of the `custom` UTE→HU map (HU per normalised UTE unit). | `[]` | Scalar; required when `skull_mapping` is `custom`. |
    | `custom_intercept` | Intercept of the `custom` UTE→HU map [HU]. | `[]` | Scalar; required when `skull_mapping` is `custom`. |
    | `debug` | Keep intermediate NIfTI files in the `pseudoCT/` subfolder. | `1` | `1 = keep`, `0 = delete after generation`. |
    | `mapping_density` | HU-to-density mapping algorithm. | `'k-plan'` | `k-plan` / `k-wave` / `marsac` / `aubry` / `air-tissue` / `none`. See [doc](doc_pseudoCT_property-mappings.md#mapping-skull-density). |
    | `density_calibration` | HU-to-density table for `mapping_density: k-plan`; use the table of the CT the UTE→HU line was fitted to. | `''` | `''` ([k-Plan](https://k-plan.io/) default) / `ucl_lowdose_80kvp` (for `treeby`) / `cph2025` / path to a CSV (HU, density [kg/m³]). See [doc](doc_pseudoCT_property-mappings.md#mapping-skull-density). |
    | `mapping_soundspeed` | HU-to-sound-speed mapping algorithm. | `'k-plan'` | `k-plan` / `marsac` / `aubry` / `webb` / `pichardo` / `mcdannold` / `none`. See [doc](doc_pseudoCT_property-mappings.md#mapping-skull-sound-speed). |
    | `mapping_attenuation` | HU-to-attenuation mapping algorithm. | `'k-plan'` | `k-plan` / `mueller` / `aubry` / `webb` / `pichardo` / `mcdannold` / `none`. See [doc](doc_pseudoCT_property-mappings.md#mapping-skull-attenuation). |

---

### Tissue Layers (`layers`)

See [doc_preproc.md](doc_preproc.md).

Maps tissue compartment names to their SimNIBS charm label indices. Compartments are modelled in listed order; any voxel label not assigned to a named compartment is treated as water. Charm label values are hardcoded in `charm_seg_labels()`. Layers can be removed or added.

| **Layer** | **Default charm labels** | **Comments** |
|---|---|---|
| `water` | `[0, 3, 6, 9, 10]` | Baseline layer; all unassigned voxels are treated as water |
| `brain` | `[1, 2]` | White matter + grey matter |
| `skin` | `[5]` | |
| `skull` | `[4]` | Single skull layer; used when cortical/trabecular split is not needed |
| `skull_cortical` | `[7]` | Cortical bone; used in multi-layer skull model |
| `skull_trabecular` | `[8]` | Trabecular bone; used in multi-layer skull model |

---

### Tissue Acoustic & Thermal Properties (`medium_properties`)

See [doc_medium.md](doc_medium.md).

Each tissue compartment (`water`, `brain`, `skin`, `skull`, `skull_trabecular`, `skull_cortical`) carries the following fields:

| **Field** | **Description** | **Units / Reference** |
|---|---|---|
| `sound_speed` | Speed of sound. | m/s; source: ITRUSST benchmarks |
| `density` | Density. | kg/m³; source: ITRUSST benchmarks |
| `alpha_coeff` | Attenuation coefficient. | dB/cm/MHz |
| `alpha_power` | Attenuation power law exponent. | dimensionless |
| `thermal_conductivity` | Thermal conductivity. | W/m/°C; source: [Tissue Properties DB](https://itis.swiss/virtual-population/tissue-properties/database) |
| `specific_heat_capacity` | Specific heat capacity. | J/kg/°C; source: Tissue Properties DB |
| `perfusion` | Perfusion / heat transfer rate. | mL/min/kg; source: Tissue Properties DB |
| `absorption_fraction` | Fraction of attenuation converted to heat. Set this or `absorption_coeff`; a config layer that sets one empties the other ([Heat deposition](doc_simulations-thermal.md#heat-deposition)). Default for every tissue. | [0–1]; source: [Pinton et al., 2012](https://doi.org/10.1118/1.3668316) |
| `absorption_coeff` | Absorption coefficient, with the power law `alpha_power`; empty by default. Set this or `absorption_fraction`. | dB/cm/MHz; source: Pinton et al., 2012 |

---

### Sonication Timing (`timing`)

See [doc_simulations-thermal.md](doc_simulations-thermal.md).

Protocol duration fields must be set for thermal simulations.

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `pd` | Pulse Duration (PD) [s]. | `0` | Duty cycle = `pd`/`pri`. |
| `pri` | Pulse Repetition Interval (PRI) [s]. | `0` | PRF = 1/`pri`. |
| `ptd` | Pulse Train Duration (PTD) [s]. | `0` | |
| `pt_timestep` | Modelling time step within a pulse train [s]; used with `equal_step_duration: 1`. | `0.02` | |
| `ptri` | Pulse Train Repetition Interval (PTRI) [s]. | `0` | OFF duration = `ptri` − `ptd`. |
| `ptrd` | Pulse Train Repetition Duration (PTRD) [s]. | `0` | |
| `post_ptri_dur` | Cooling period after the PTRD [s]. | `0` | |
| `post_pt_timestep` | Modelling time step between pulse trains and after the PTRD [s]. | `1` | |
| `equal_step_duration` | Equal step durations for on and off cycles? | `0` | `1 = yes`, `0 = no` |
| `pulse_ramp_dur` | Rise (= fall) duration of the amplitude ramp at each end of a pulse [s]. | *(not set)* = `0` | Not in `config_default.yaml`. Scales the pulse-average intensity: Ipa = factor × Itp (`compute_ramp_factor`). |
| `pulse_ramp_shape` | Ramp shape: `'linear'` or `'tukey'` (raised cosine). | *(not set)* = `'linear'` | Used with `pulse_ramp_dur`. |

The Nexus draws this protocol (pulse, pulse train, repetitions) and exchanges it with the [TUS calculator](https://www.socsci.ru.nl/fusinitiative/tuscalculator/); see [Nexus](doc_nexus_tabs.md#timing-diagram-and-tus-calculator).

---

### Thermal Settings (`thermal`)

See [doc_simulations-thermal.md](doc_simulations-thermal.md).

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `cem43_iso` | CEM43 form used in the HTML report (dashboard, safety verdict, methods text): ISO 14708-3 / ITRUSST (`1`) or k-Wave `kWaveDiffusion` (`0`). | `0` | Both forms are always computed and written to CSV, NIfTI and figures. With `thermal_backend: babelbrain` it is set to `0` with a warning (`prestus:cem43IsoBabelBrain`): BabelBrain does not compute the ISO form. |
| `temp_0.{tissue}` | Initial temperature per tissue compartment [°C]. | `37` | Tissues: `water`, `skull`, `brain`, `skin`, `skull_trabecular`, `skull_cortical` |
| `log_thermal_updates` | Show k-Wave's thermal output for every time step. | `0` | `0`: only the first step is shown in full. |

---

### Time Reversal (`time_reversal`)

See [Time reversal](doc_time_reversal.md).

??? note "Time reversal keys"

    One extra k-Wave run per target (point source at the target), between the free-water baseline and the acoustic simulation.

    | **Parameter** | **Description** | **Default** | **Comments** |
    |---|---|---|---|
    | `enabled` | Run the reverse simulation and write the refocusing / aberration report. | `false` | Cartesian 2D/3D grids with the k-Wave source setup; not with `acoustic_backend: babelbrain` (use `babelbrain.do_refocusing`). |
    | `apply_phases` | Use the time-reversal phases in the acoustic run. | `false` | Matrix: per element; annular: per ring; ignored for a single element. ISPPA targeting keeps the drive amplitudes of the geometric drive. |
    | `save_field` | Save the complex reverse field (sensor region) with the results. | `false` | `<dir_cache>/<prefix>_timereversal.mat`. |
    | `placement_search.enabled` | Score transducer poses around the planned one on the same reverse run. | `false` | Curved annular transducers, 3D. Suggestions only: the acoustic run keeps the planned pose. |
    | `placement_search.max_angle_deg` | Largest tilt of the beam direction from the planned one [deg]. | `20` | Poses lie on an orbit around the target, against the scalp. |
    | `placement_search.angle_step_deg` | Angular step of the candidate directions [deg]. | `2.5` | |
    | `placement_search.aim_iterations` | Re-aim at target − predicted shift this many times. | `2` | `0` = geometric aim only. |
    | `placement_search.jitter_mm` | Also score the planned pose with its apex moved by these distances [mm], aimed at the target. | `[]` | One or more, e.g. `0.5` (the local search users often run with forward simulations) or `[0.5 1 2 5]`. Rows of kind `jitter`, distance in column `jitter_mm`. |
    | `placement_search.jitter_mode` | How `jitter_mm` moves the apex. | `'grid'` | `'grid'`: ± along each grid axis, as when editing voxel coordinates (6 poses per distance; the bowl is not re-seated, so beyond ~1 mm it can cut into or lift off the scalp). `'scalp'`: slid over the scalp perpendicular to the planned beam in 8 directions and re-seated against the scalp (8 poses per distance). |
    | `only` | Stop after the time-reversal stage. | `false` | No water baseline, acoustic or thermal simulation: the report and pose table only. |
    | `pose_workflow.enabled` | Run the pose workflow (search → forward verification → final full-grid run) instead of a single simulation. | `false` | One curved annular transducer, `simulation.medium: layered`. See [Pose workflow](doc_time_reversal.md#pose-workflow). |
    | `pose_workflow.n_top` | Number of best poses (by predicted target pressure after aim correction) to verify. | `5` | The planned pose is always verified as the reference. |
    | `pose_workflow.poses` | Pose ids of the pose table to verify as well. | `[]` | |
    | `pose_workflow.aim` | Aim of the verification runs. | `'corrected'` | `'corrected'` (aberration-corrected aim point), `'geometric'` (at the target) or `'both'`. |
    | `pose_workflow.verify_fov` | Grid of the verification runs. | `'partial'` | `'partial'`: acoustic FOV around each beam (`grid.acoustic_fov_*`, `'auto'` unless set); `'full'`; `'config'`: as configured. |
    | `pose_workflow.target_radius_mm` | Radius of the target sphere over which the simulated intensity is averaged [mm]. | `2` | |
    | `pose_workflow.rank_by` | How the best pose is chosen. | `'target_intensity'` | `'target_intensity'` (highest mean intensity in the target sphere) or `'focus_error'` (smallest distance from the −3 dB focal centroid to the target). |
    | `pose_workflow.max_focus_error_mm` | Only poses whose −3 dB focal centroid lies within this distance of the target qualify [mm]. | `[]` | `[]` = no limit. |
    | `pose_workflow.final_run` | Run the best pose on the full grid with the configured modules. | `true` | Affix `<affix>_posebest`. |
    | `pose_workflow.search_memorylimit` | `hpc.memorylimit` of the search run [GB]. | `[]` | `[]` = `hpc.memorylimit`. The search records the field over all candidate poses. |

---

### Output Analysis (`analysis`)

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `focus_area_radius` | Radius around the focus in which IPA is averaged for outputs [mm]. | `5` | |
| `target_mni_mm` | Intended target `[x y z]` in MNI mm for the targeting analysis. | `[]` | `[]` = the MNI target of the placement mode (`placement.heuristic`/`plantus.mni_target_mm`, `placement.mni.focus_pos_mm`), if any. Can differ per run (study manifest column). |
| `target_roi_mni` | Target ROI mask in MNI space (NIfTI) for the targeting analysis and the group maps. | `''` | Overlap with the −3 dB focal volume, ROI coverage, intensity in the ROI; outlined on the group overlap figures. |

---

### Study Setup (`study`)

Multi-subject studies with per-run settings; see [Group analysis](doc_group.md#study-setup).

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `manifest` | Study manifest (`.tsv`/`.csv`): one row per run (`subject_id`, optional `condition`, `include`, `config`) plus any config key as a column. | `''` | Relative paths resolve against the study config's folder. Run with `prestus_study_start(<config>)`; report with `prestus_group_report_start(<config>)`. |
| `name` | Name of the study outputs folder `<path.sim>/group/<name>/`. | `''` | `''` = manifest file name without extension. |

---

### High-Performance Computing (`hpc`)

See [k-Wave](doc_kwave.md) and [HPC](doc_hpc.md).

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `gpu` | Request a specific GPU. | `''` | e.g. `"nvidia_a100-sxm4-40gb:1"` |
| `partition` | Request a dedicated queue partition. | `''` | |
| `reservation` | Request a reserved queue. | `''` | |
| `exclude` | Nodes to avoid, comma-separated (SLURM `sbatch --exclude`). | `''` | e.g. a node whose GPU is held by a stale process. |
| `wait_for_job` | Block MATLAB until the submitted HPC job completes? | `false` | |
| `timelimit` | Job time limit. | `'04:00:00'` | |
| `memorylimit` | Memory limit [GB]. | `20` | |
| `posthoc_water_timelimit` | Time limit of the post-hoc free-water jobs (`modules.run_posthoc_water_sims`). | `'05:00:00'` | |
| `posthoc_water_memorylimit` | Memory limit of the post-hoc free-water jobs [GB]. | `40` | |
| `ld_library_path` | LD_LIBRARY path for SimNIBS installation. | `''` | Set if you see `undefined symbol` errors. e.g. `/opt/gcc/7.2.0/lib64` |
| `job_prefix` | Prefix string for HPC job names. | `'PS'` | Keep short; the job name is capped at 20 characters. Overridden to `'TP'` for transducer positioning jobs. |
| `max_wait_checks` | Maximum number of job status checks when `wait_for_job = true`. | `540` | At 1 check/20 s, `540` ≈ 3 hours. |
| `depend_type` | Scheduler dependency type for `depend_job_id(s)`. | `'afterok'` | `'afterok'` or `'afterany'`. Internal (not in `config_default.yaml`); `prestus_study_start` uses `'afterany'` for the study-report job. |


#### Snellius GPU Partitions (`snellius`)

Used when `hpc.name: 'snellius'`: `hpc.partition` selects one of these blocks, which sets the job's resources.

??? note "Snellius partition keys"

    | **Parameter** | **Description** | **Default** | **Comments** |
    |---|---|---|---|
    | `gpu_a100.cores` | CPU cores requested on `gpu_a100`. | `18` | |
    | `gpu_a100.memorylimit` | Memory [GB] on `gpu_a100`. | `120` | |
    | `gpu_a100.timelimit` | Wall time on `gpu_a100`. | `'120:00:00'` | Partition maximum. |
    | `gpu_a100.n_gpu` | GPUs requested on `gpu_a100`. | `1` | |
    | `gpu_h100.cores` | CPU cores requested on `gpu_h100`. | `16` | |
    | `gpu_h100.memorylimit` | Memory [GB] on `gpu_h100`. | `180` | |
    | `gpu_h100.timelimit` | Wall time on `gpu_h100`. | `'120:00:00'` | Partition maximum. |
    | `gpu_h100.n_gpu` | GPUs requested on `gpu_h100`. | `1` | |

---

### Transducer Calibration (`calibration`)

See [doc_calibration.md](doc_calibration.md).

??? note "Calibration keys"

    A separate `config_calibration.yaml` applies for calibration workflows and is loaded as `parameters.calibration`.

    | **Parameter** | **Description** | **Comments** |
    |---|---|---|
    | `path_input_axial` | Directory containing axial intensity profiles; `char_data_path` in the equipment YAML is relative to it. | Usually `<PRESTUS>/local/profiles/` (see doc_transducer_library.md, Appendix). |
    | `path_output` | Directory for saving free-water simulation results. | |
    | `path_output_profiles` | Directory for saving optimised profile data. | |
    | `filename_calibrated_CSV` | Filename of calibrated CSV data. | Mandatory only when not generated within standalone script. |
    | `save_in_calibration_folder` | Save in `path_output` (`TRUE`) or `sim_path` (`FALSE`). | If `TRUE`, results are appended to existing calibration data. |
    | `combinations` | Equipment combinations (must refer to entries in `config_equipment.yaml`). | Multiple combinations can be specified. |
    | `focal_depths_wrt_exit_plane` | List of focal depths to characterise [mm]. | |
    | `desired_intensities` | Desired free-water intensities [W/cm²]. | |
    | `add_FDO` | Append Focal Distance Offset (bowl-to-exit-plane distance)? | Set to `1` if zero point in profiles reflects exit plane rather than bowl. |
    | `axisymmetric2D` | Use axisymmetric 2D instead of default 3D free-water simulations? | `1 = yes`, `0 = no` |
    | `force_kwavearray` | Force use of kWaveArray for free-water simulations? | If `0`, uses the setting in the default/study config. |
    | `opt_method` | Optimisation method. | [`FEXminimize`](https://github.com/rodyo/FEX-minimize) (open source) / `GlobalSearch` (MATLAB [Global Optimization Toolbox](https://www.mathworks.com/products/global-optimization.html)) |
    | `opt_limits` | Distance limits for optimisation [mm]. | |
    | `opt_weights` | Weighting of original profile during fitting. | `1` = equal; `>1` = Gaussian (narrower with larger values) |
    | `opt_seed` | Random seed for optimisation. | |
    | `opt_use_initial_phases` | Use phases from `transducer.annular.elem_phase_rad` as the initial guess instead of random initialization. | `false` (default); set to `true` when manufacturer-provided or previously calibrated phases are available. |
    | `opt_upper_velocity` | Upper velocity in global search | |
    | `opt_phase_precession` | Constrain element phases during optimisation. | `false` (default) = independent; `'linear'` = strict ramp (`phase_start + i * phase_step`, 2 free params); `'monotonic'` = ordered phases via cumulative non-negative increments (N params, monotonically constrained). |
    | `opt_amp_validation` | Run free-water simulation to validate phases and amplitudes?| 'always' (default) / 'never' / 'initial' (only for the first target intensity) / 'final' (only for the final target intensity) | |
    | `skip_front_peak_mm` | Distance from profile start to ignore [mm]. | Avoids near-field artefacts in peak/FWHM calculations. |
    | `profile_align` | Focus alignment when interpolating measured profiles to an unmeasured focal depth. | `'label'` (default) = column focal-depth headers; `'peak'` = intensity peak beyond `skip_front_peak_mm`. See `interpolate_axial_profiles`. |
    | `max_extrapolation_mm` | How far outside the measured depth range a desired depth may lie [mm]. | `0` (default) = no extrapolation (error). Otherwise the nearest measured profile, shifted to the desired depth. |
    | `focal_reference` | Convention for what a quoted focal depth means. | `flhm_3db` (default, ITRUSST, intensity ≥ peak/2) / `flhm_6db` ([Imasonic](https://www.imasonic.com/)/IGT, intensity ≥ peak/4) / `peak` (raw intensity peak). Recorded into the sidecar YAML's `meta.calibration.focal_reference`. See [doc_calibration.md](doc_calibration.md). |
