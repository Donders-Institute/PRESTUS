# Changelog

Notable changes to this project are documented here. 

---

## v0.6.0-development

Extends the calibration pipeline with a unified entry point, parametric model store, and HPC dispatch. Introduces async multi-transducer simulations, PlanTUS-based automated placement, and RAS as the canonical coordinate space.

#### Added

- [**calibration**] HPC dispatch, multi-intensity mode, and parametric calibration pipeline
- [**calibration**] Geometric steering mode with per-element hardware correction
- [**calibration**] Multiple analytical forward models — O'Neil & Rayleigh model
- [**calibration**] Transducer library draft
- [**feature**] Async multi-transducer pipeline (see [doc_async_transducer.md](doc_async_transducer.md))
- [**feature**] Sequential simulation dispatcher with multi-run report and thermal chaining (see [doc_advanced.md](doc_advanced.md))
- [**feature**] PlanTUS placement mode — wraps external PlanTUS optimiser for automated transducer targeting (see [doc_placement_plantus.md](doc_placement_plantus.md))
- [**feature**] [Draft] Browser-based real-time transducer alignment viewer

#### Changed

- [**calibration**] Variable names generalised in `perform_global_search` and `phase_optimization_annulus_full_curve`; `opt_use_initial_phases` flag introduced
- [**coord**] RAS established as canonical entry-point coordinate space; phantom and water pipelines use a RAS anchor
- [**core**] `transducer.target_isppa_wcm2` migrated from calibration struct to transducer struct
- [**telemetry**] Sequential pipeline mode detected and reported; version string (tag/VERSION) reported separately from git hash
- [**core**] `MergeStruct` renamed to `mergestruct`; external copy removed
- [**calibration**] Free-water validation can now be deactivated or run only for the first/last target ISPPA
- [**head**] `tissuemask_binary` now exposes a `water` mask field
- [**head**] `smooth_img` gains an `off` method (no smoothing); exposed in GUI
- [**head**] RAS+ orientation enforced at NIfTI load time; `ras_plus` grid mode rescales without rotation
- [**core**] Image Processing Toolbox calls replaced with base-MATLAB helpers throughout (no toolbox required for morphological ops, resampling)
- [**config**] Optional `transducer.name` field
- [**nifti**] Batched `tformarray` with optional parallel workers for MNI transforms
- [**acoustic**] Save complex pressure field; `save_p_complex` flag guards time-series recording

#### Fixed

- [**calibration**] Amplitude scaling bug present since early PRESTUS versions: correction factor must be the square root of intensity ratios, not the ratio itself
- [**pCT**] pseudoCT path layout corrected; NaN padding used instead of zero fill
- [**nifti**] Zero fill replaced with explicit FOV mask in MNI outputs; `FillValues` option added to `convert_final_to_MNI_simnibs`
- [**calibration**] Cache folder removed after calibration when `save_matrices` is disabled
- [**telemetry**] Telemetry check added for standalone calibration runs
- [**thermal**] `thermal_plot_sim` legend guarded against empty tissue masks
- [**report**] Sequential report: fallback paths and real-time time axis corrected
- [**transform**] `isequal` used for size comparison; `affine3d` replaces deprecated `affineTransform3d`

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
- [**feature**] Multi-ISPPA targeting — run acoustics once, scale to N independent thermal targets via post-hoc pressure scaling (see [doc_multi_isppa.md](doc_multi_isppa.md))
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
