# Changelog

Notable changes to this project will be documented here.

## `development` v0.5.0 [*unreleased*]

Note: Parameters, naming, and default values have changed. Please consult the [documentation](https://donders-institute.github.io/PRESTUS/) for the current parameter specification.

#### Added

- [**feature**] Uncertainty quantification pipeline by @jkosciessa
    - Three-variant simulation workflow (default / liberal / conservative medium properties) producing a unified HTML report with safety dashboard, acoustic and thermal range tables, timing summary, and side-by-side maps
    - Supports MATLAB (sequential) and HPC (parallel, scheduler-dependent) execution with automatic resumption from partial runs
- [**feature**] Matrix transducer support (see [PR #117](https://github.com/Donders-Institute/PRESTUS/pull/117))
    - Flexible element layout options (rectangular grid, Fibonacci, Fermat spiral) and Clover multi-array configuration
    - Element positions can be defined inline or loaded from file
- [**feature**] Heuristic transducer placement
    - Automated skull-surface search with configurable positioning criteria and ear exclusion
- [**feature**] C++ backend CPU & GPU support by @jkosciessa
- [**feature**] HTML simulation report by @sirmrmarty
    - Self-contained per-subject summary of acoustic and thermal outputs
- [**feature**] Parallel multi-transducer support for layered simulations by @dreamstimlab (see [PR #100](https://github.com/Donders-Institute/PRESTUS/pull/100))
- [**feature**] Per-step pipeline benchmarking: wall time, peak RAM, and disk usage logged at each stage
- [**grid**] Minimum PPW validation at source setup (`grid.min_ppw`, default 6); warns with the maximum `resolution_mm` that would satisfy the requirement
- [**skull**] Rubber-wrap skull mask inflation by @meijer-s (see [PR #103](https://github.com/Donders-Institute/PRESTUS/pull/103))
- [**thermal**] More flexible and standardised protocol timing specification with output plot
- [**thermal**] Additional thermal timeseries outputs (max per tissue, heating at focus)
- [**doc**] GitHub Pages [documentation](https://donders-institute.github.io/PRESTUS/) covering parameters, functions, installation, quick start, and processing steps

#### Changed

- [**refactor**] Pipeline and core processing steps refactored; functions organised by processing stage
- [**validation**] Function inputs are now consistently validated via MATLAB `arguments` blocks; `help <function>` calls expose header information inline. Some checks may surface new errors — please report unexpected failures
- [**calibration**] Dedicated config and unified pipeline for transducer calibration [@MaCuinea, @jkosciessa]
- [**hpc**] HPC submission refactored; platform selection unified under `parameters.platform`
- [**preproc**] Segmentation smoothing now specified in FWHM mm; grid interpolation updated for continuous data
- [**pCT**] Updated to new folder structure; skull density and sound speed regularised to water minimum
- [**pCT**] pCT generation is now fully integrated into the pipeline; pre-generating pCTs beforehand remains possible (e.g. for debugging) but is no longer required — a T1w image, a UTE image, and a mapping algorithm are sufficient
- [**refactor**] PML (Perfectly Matched Layer) is now added only at the acoustic wave simulation stage, rather than carried through from preprocessing or phantom generation; the default `"auto"` setting optimises PML size automatically
- [**I/O**] Subject-specific output subfolders are now enforced to simplify the I/O landscape
- [**axisymmetric**] Axisymmetric acoustic simulations now default to 3D thermal output
- [**parameter**] Various parameters renamed or restructured; consult `default_config.yaml` and the documentation for current names

#### Fixed

Not all (hot-)fixes are reported here.

- [**source**] Each uncertainty variant now recomputes its own k-Wave source with its own time axis; shared caching across variants with different medium sound speeds could produce incorrect pressure maps
- [**paths**] `load_parameters` resolves `default_config.yaml` via an absolute path, preventing stray output folders when called from an unexpected working directory
- [**calibration**] Potential bugs with transducer distance calculation
- [**savemat**] k-Wave source matrix saving can now be correctly deactivated; water simulations force-save matrix outputs

#### Deprecated/Removed

- [**feature**] Removed explicit ‘skull’-only simulations; equivalent behaviour is achieved by removing layers in the simulation setup

## v0.4.0

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

## v0.3.0

- PseudoCT feature by @eleonoracarpino and @jkosciessa
- Clean-up repository by @jkosciessa
- Incorrect allocation of attenuation values bug fix by @jkosciessa
- Improve correct labeling of tissue by @jkosciessa
- Improvement of 3D plots by @jkosciessa
- Check simulation stability and adjust time step accordingly by @MaCuinea and @jkosciessa
- Enable tissue-specific temp_0 specification by @jkosciessa
- Fix heating plot bug by @jkosciessa

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

## v0.1.0

Initial in-house development version by @achetverikov

- Dedicated segmentation path option, minor bug fixes, updated simnibs call by @jkosciessa in #2
- fix small errors by @jkosciessa in #3
- adjusted skull values by @eleonoracarpino in #4
- Bug fixes and documentation changes by @KTZ228 in #6
- DOC: simplify readme.md wording regarding requirements by @mekman in #7