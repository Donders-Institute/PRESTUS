# Changelog

Notable changes to this project will be documented here.

## `development` v0.5.0 [*unreleased*]

Note: While intended to be largely backwards compatible, parameters, their naming, and default values have changed. Please consult the [documentation](https://donders-institute.github.io/PRESTUS/) for the current parameter specification. 

#### Added

- [**feature**] C++ support by @jkosciessa
- [**feature**] Advanced logging & benchmarking
    - Key parameters are displayed at onset of simulation
    - PRESTUS and k-Wave commit version logged (git-only)
- [**feature**] html summary output by @sirmrmarty
    - Acoustic properties will only be printed for modeled layers @jkosciessa
- [**feature**] Parallel multi-transducer support for layered simulation by @dreamstimlab (see [PR](https://github.com/Donders-Institute/PRESTUS/pull/100))
- [**skull**] new option `parameters.rubberwrap` by @meijer-s (see [PR](https://github.com/Donders-Institute/PRESTUS/pull/103))
    - Locally inflate the skull mask. 
    - New default for the layered simulation and only option during pCT creation. The sequence of where this is implemented currently varies between layered/pCT. For pCT, the rubber expansion is implemented following segmentation. For layered, it follows medium map creation (i.e., following the smoothing of masks) to catch potential issues where smoothing and binarization may (re-)introduce holes. 
- [**doc**] GitHub pages [documentation](https://donders-institute.github.io/PRESTUS/)
    - Full parameter & function overview
    - Installation Guide
    - Quick Start Guide
    - Example demos
    - Draft of processing steps
- [**thermal**] More flexible and standardised protocol timing specification
- [**thermal**] Output plot of requested timing protocol
- [**thermal**] Additional thermal timeseries of max. in tissue and heating at focus
- [**thermal**] Heating timeseries are saved as .mat
- [**parameter**] I/O handling of intermediate outputs. `savemat` can be set to 0 to not save processing matrices (save disk space)
- [**parameter**] Dedicated `debug` parameter
    - saves plots of postprocessing etc in debug folder of sim directory
    - parameter overview when loading config (cf. pipeline onset)
    - saves nifti images of raw acoustic properties
    - currectly active by default

#### Changed

- [**refactor**] Pipeline and main processing steps refactored
    - Functions organized by relevant processing stage
    - main pipeline: 1150 lines -> 450 lines (incl. updated logging)
- [**calibration**] Dedicated config and unified pipeline for transducer calibration, new parameters. Function for calibration. [@MaCuinea, @jkosciessa]
- [**parameter**] Minor parameters changes: some were dropped, some renamed, others added, check the function doc if in doubt. default_config and documentation updated accordingly. Behaviour should be comparable to v0.4.0
- [**pCT**] pCT creation updated to new folder structure and simplified based on @meijer-s  implementation.
- [**savemat**] kwave source matrix saving now can be deactivated (as intended)
- [**feature**] updates to axisymmetry support
- [**advanced**] Updates to the sequential simulations with different transducer-target specs @Kenneth van der Zee
- [**preproc**] segmentation smoothing choices changed @jkosciessa
- [**preproc**] grid interpolation for continuous data changed @jkosciessa

#### Fixed

Not all (hot-)fixes are reported here.

- [**calibration**] potential bugs with distance calculation
- [**savemat**] kwave source matrix saving now can be deactivated (as intended)
- [**savemat**] water simulations will force-save matrix outputs

#### Deprecated/Removed

- [**feature**] Removed explicit ‘skull’-only simulations. These can be run by removing layers in the simulation setup.

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