[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15095860.svg)](https://doi.org/10.5281/zenodo.15095860)

![PRESTUS logo](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/logo_PRESTUS.png)

# PRESTUS: PREprocessing & Simulations for Transcranial Ultrasound Stimulation

PRESTUS (PREprocessing & Simulations for Transcranial Ultrasound Stimulation) is an open-source MATLAB toolbox that aims to streamline imaging-informed simulations of Transcranial Ultrasound Stimulation (TUS): from the segmentation of T1-weighted MRI head scans, [preprocessing](doc_preproc.md) of their tissue layers, and mapping of [acoustic tissue properties](doc_medium.md) (possibly informed for skull via [(pseudo-)CT](doc_pseudoCT.md) images) in a simulation grid, to the execution of [acoustic](doc_simulations-acoustic.md) and [thermal](doc_simulations-thermal.md) simulations using the widely adopted [k-Wave engine](doc_backend.md). High-performance computing ([HPC](doc_hpc.md)) support via SLURM and CentOS enables efficient parallelization, and large-scale analyses. Output 3D NifTI images are automatically mapped to a standard template (MNI) space to facilitate [group](doc_group.md) reporting.

Key features of the modular end-to-end pipeline include:

- Automated MRI segmentation (using SimNIBS 4 charm) and [preprocessing](doc_preproc.md).
- 2D / 3D grid setup.
- Multi-layer [medium property mapping](doc_medium.md) (water, skin, multi-layer skull, and brain).
- [(pseudo-)CT](doc_pseudoCT.md)-informed continuous skull mapping.
- Virtual multi-element [transducer calibration](doc_calibration.md) (free-water profile emulation).
- Estimation of entry-target [coordinates](doc_placement.md) incl. [heuristic placement](doc_placement_heuristic.md) and [neuronavigation](doc_placement_neuronav.md) read-in.
- Flexible [temporal protocol specification](doc_simulations-thermal.md) (e.g., including breaks).
- k-Wave integration for robust [acoustic](doc_simulations-acoustic.md) and [thermal](doc_simulations-thermal.md) simulations.
- Support for [high-performance computing](doc_hpc.md).
- 3D NifTI [outputs](doc_outputs.md) for major reporting metrics (in subject- & MNI-space).

These features make PRESTUS an accessible tool for researchers to plan, validate, and report transcranial ultrasound targeting. PRESTUS is intended solely for basic research use in non-clinical applications.

A recent overview poster can be found [here](https://jkosciessa.github.io/downloads/2025-FUN25-PRESTUS.pdf).

# Where to start

New to PRESTUS? Follow these in order:

1. [Concepts and glossary](doc_concepts.md): what PRESTUS does, and the terms used throughout.
2. [Installation](doc_installation.md): MATLAB, PRESTUS and k-Wave are enough to begin.
3. [Your first simulation](doc_getting-started.md): a free-water run in a few minutes, then what a head simulation needs.
4. [Reading your results](doc_results.md): what to check in the report, and how to read the safety summary.
5. [Configuring a study](doc_configuration.md): your own config files for your own data.

Looking for something specific?

| I want to… | Go to |
|---|---|
| simulate a transducer I own | [Transducer Library](doc_transducer_library.md) |
| decide where the transducer goes | [Placement](doc_placement.md) |
| simulate heating for my protocol | [Thermal](doc_simulations-thermal.md) |
| simulate a pose recorded with Localite | [Localite placement](doc_placement_neuronav.md) |
| run many subjects | [Studies & group analysis](doc_group.md) |
| run on a cluster | [Running on a cluster](doc_hpc.md) |
| look up a setting | [Parameters](doc_parameters.md) |
| fix an error | [Troubleshooting](doc_troubleshooting.md) |

# Installation

Follow the [Installation Guide](doc_installation.md), run [your first simulation](doc_getting-started.md), and set up your own study with [Configuring a study](doc_configuration.md).

### Versions

PRESTUS is under active development. The default `main` branch is considered more stable. The [`development` branch](https://github.com/Donders-Institute/PRESTUS/tree/development) may contain both new features and stability updates, but is subject to more dynamic changes. [Releases](https://github.com/Donders-Institute/PRESTUS/releases) bundle major updates and are associated with a version number. Starting in PRESTUS v0.5 logs will print the git hash when PRESTUS has been cloned as a git repository.

# License

Released under GNU General Public License v3.0 (see LICENSE).

> **Disclaimer**
> This software is currently under development and is provided “as is” without warranty of any kind, either express or implied, including but not limited to the implied warranties of merchantability, fitness for a particular purpose, and non-infringement. This tool is intended for research purposes only and is not designed, intended, or approved for medical or clinical use, diagnosis, or treatment of patients. The contributors accept no liability for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limited to, procurement of substitute goods or services; loss of use, data, or profits; or business interruption) arising in any way out of the use of this software. Use at your own risk.

# Development

If you would like to assist the community development of this tool, please consider [CONTRIBUTING](CONTRIBUTING.md).