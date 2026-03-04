[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15095860.svg)](https://doi.org/10.5281/zenodo.15095860)

![PRESTUS logo](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/logo_PRESTUS.png)

# PRESTUS: PREprocessing & Simulations for Transcranial Ultrasound Stimulation

PRESTUS (PREprocessing & Simulations for Transcranial Ultrasound Stimulation) is an open-source MATLAB toolbox that aims to streamline imaging-informed simulations of Transcanial Ultrasound Stimulation (TUS): from the segmentation of T1-weighted MRI head scans, [preprocessing](doc_preproc.md) of their tissue layers, and mapping of [acoustic tissue properties](doc_medium.md) (possibly informed for skull via [(pseudo-)CT](doc_pseudoCT.md) images) in a simulation grid, to the execution of [acoustic](doc_simulations-acoustic.md) and [thermal](doc_simulations-thermal.md) simulations using the widely adopted [k-Wave engine](doc_backend.md). High-performance computing ([HPC](doc_hpc.md)) support via SLURM and CentOS enables efficient parallelization, and large-scale analyses. Output 3D NifTI images are automatically mapped to a standard template (MNI) space to facilitate [group](doc_group.md) reporting.

Key features of the modular end-to-end pipeline include:

- Automated MRI segmentation (using SimNIBS 4 charm) and [preprocessing](doc_preproc.md).
- 2D / 3D grid setup.
- Multi-layer [medium property mapping](doc_medium.md) (water, skin, multi-layer skull, and brain).
- [(pseudo-)CT](doc_pseudoCT.md)-informed continuous skull mapping.
- Virtual multi-element [transducer calibration](doc_calibration.md) (free-water profile emulation).
- Estimation of entry-target [coordinates](doc_placement.md) incl. [neuronavigation](doc_neuronav.md) read-in.
- Flexible [temporal protocol specification](doc_simulations-thermal.md) (e.g., including breaks).
- k-Wave integration for robust [acoustic](doc_simulations-acoustic.md) and [thermal](doc_simulations-thermal.md) simulations.
- Support for [high-performance computing](doc_hpc.md).
- 3D NifTI [outputs](doc_outputs.md) for major reporting metrics (in subject- & MNI-space).

These features make PRESTUS an accessible tool for researchers to plan, validate, and report transcranial ultrasound targeting. PRESTUS is intended solely for basic research use in non-clinical applications.

A recent overview poster can be found [here](https://jkosciessa.github.io/downloads/2025-FUN25-PRESTUS.pdf).

# Installation

Follow the [Installation Guide](doc_installation.md) and continue with the [Quick Start Guide](doc_getting-started.md).

### Versions

PRESTUS is under active development. The default `main` branch is considered more stable. The [`development` branch](https://github.com/Donders-Institute/PRESTUS/tree/development) may contain both new features and stability updates, but is subject to more dynamic changes. [Releases](https://github.com/Donders-Institute/PRESTUS/releases) bundle major updates and are associated with a version number. Starting in PRESTUS v0.5 logs will print the git hash when PRESTS has been cloned as a git repository.

# License

Released under GNU General Public License v3.0 (see LICENSE).

> **Disclaimer**
> This software is currently under development and is provided “as is” without warranty of any kind, either express or implied, including but not limited to the implied warranties of merchantability, fitness for a particular purpose, and non-infringement. This tool is intended for research purposes only and is not designed, intended, or approved for medical or clinical use, diagnosis, or treatment of patients. The contributors accept no liability for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limited to, procurement of substitute goods or services; loss of use, data, or profits; or business interruption) arising in any way out of the use of this software. Use at your own risk.

# Development

If you would like to assist the community development of this tool, please consider [CONTRIBUTING](CONTRIBUTING.md).