[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15095860.svg)](https://doi.org/10.5281/zenodo.15095860)

# PRESTUS: PREprocessing & Simulations for Transcranial Ultrasound Stimulation

PRESTUS (PREprocessing & Simulations for Transcranial Ultrasound Stimulation) is an open-source MATLAB toolbox that aims to streamline imaging-informed simulations of Transcanial Ultrasound Stimulation (TUS): from the segmentation of T1-weighted MRI head scans and mapping of medium tissue properties (possibly informed for skull via (pseudo-)CT images) in a simulation grid, to the execution of acoustic and thermal simulations using the widely adopted k-Wave engine. High-performance computing (HPC) support via SLURM and CentOS enables efficient parallelization, and large-scale analyses. Output 3D NifTI images are automatically mapped to a standard template (MNI) space to facilitate group reporting.

Key features include:
- Automated MRI segmentation (using SimNIBS 4 charm) and preprocessing.
- 2D / 3D grid setup.
- Multi-layer medium property mapping (water, skin, multi-layer skull, and brain).
- (pseudo-)CT-informed continuous skull mapping.
- Virtual multi-element transducer calibration (free-water profile emulation).
- Estimation of entry-target coordinates.
- Flexible temporal protocol specification (e.g., including session breaks).
- k-Wave integration for robust acoustic and heating simulations.
- Support for high-performance computing.
- 3D NifTI outputs for major reporting metrics (in subject- & MNI-space).

These features make PRESTUS an accessible tool for researchers to plan, validate, and report transcranial ultrasound targeting. 
The ```documentation``` folder contains evolving documentation of the main workflows and parameters, alongside a tutorial on how to run the integrated pipeline.

A recent overview poster can be found [here](https://jkosciessa.github.io/downloads/2025-FUN25-PRESTUS.pdf).

If you would like to assist the community development of this tool, please consider [CONTRIBUTING](CONTRIBUTING.md).

# Installation

## Donders Institute HPC Cluster

When working on the Donders High-Performance-Computing cluster, PRESTUS and its dependencies (SimNIBS, k-Wave) are already installed. 

Type ``module load simnibs/4.0.0`` (or add the command to your .bashrc so that it is executed automatically once you login) and add ``addpath('/opt/prestus/dev')`` to your matlab path. This will use the most up-to-date version of PRESTUS (i.e., the current development branch). If you want to use an older version you can also use ``addpath('/opt/prestus/0.2.0')``, or older versions. Now you can start matlab R2022b.

For more information on HPC usage, see [doc_hpc](documentation/doc_hpc.md)

If you want to get started with simulations, you can use the PRESTUS example dataset. This command will copy the dataset to your home directory:

```
cp /opt/prestus/example_data/PRESTUS_example_data.zip ${HOME}
```

## Outside the Donders HPC

Follow the [installation guide](documentation/doc_installation.md).

# Reference

Chetverikov, A.\*, Kosciessa, J. Q.\*, Cornelissen, M., Carpino, E., van der Zee, K., & Verhagen, L. (2025). PRESTUS (0.4.0). Zenodo. https://doi.org/10.5281/zenodo.15965832

# Contributors

[![Contributors](https://img.shields.io/github/contributors/Donders-Institute/PRESTUS.svg?color=00B4D8&style=flat-square)](https://github.com/Donders-Institute/PRESTUS/graphs/contributors)
[![All Contributors](https://img.shields.io/github/all-contributors/Donders-Institute/PRESTUS?color=00B4D8&style=flat-square)](#contributors)
<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="http://juliankosciessa.eu"><img src="https://avatars.githubusercontent.com/u/40263608?v=4?s=100" width="100px;" alt="Julian Kosciessa"/><br /><sub><b>Julian Kosciessa</b></sub></a><br /><a href="#code-jkosciessa" title="Code">💻</a> <a href="#ideas-jkosciessa" title="Ideas, Planning, & Feedback">🤔</a> <a href="#tutorial-jkosciessa" title="Tutorials">✅</a> <a href="#maintenance-jkosciessa" title="Maintenance">🚧</a> <a href="#bug-jkosciessa" title="Bug reports">🐛</a> <a href="#doc-jkosciessa" title="Documentation">📖</a></td>
      <td align="center" valign="top" width="14.28%"><a href="http://andreychetverikov.org"><img src="https://avatars.githubusercontent.com/u/1465806?v=4?s=100" width="100px;" alt="Andrey Chetverikov"/><br /><sub><b>Andrey Chetverikov</b></sub></a><br /><a href="#code-achetverikov" title="Code">💻</a> <a href="#ideas-achetverikov" title="Ideas, Planning, & Feedback">🤔</a> <a href="#tutorial-achetverikov" title="Tutorials">✅</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/KTZ228"><img src="https://avatars.githubusercontent.com/u/51954604?v=4?s=100" width="100px;" alt="Kenneth van der Zee"/><br /><sub><b>Kenneth van der Zee</b></sub></a><br /><a href="#code-KTZ228" title="Code">💻</a> <a href="#ideas-KTZ228" title="Ideas, Planning, & Feedback">🤔</a> <a href="#maintenance-KTZ228" title="Maintenance">🚧</a> <a href="#bug-KTZ228" title="Bug reports">🐛</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/MaCuinea"><img src="https://avatars.githubusercontent.com/u/134381864?v=4?s=100" width="100px;" alt="Margely Cornelissen"/><br /><sub><b>Margely Cornelissen</b></sub></a><br /><a href="#code-MaCuinea" title="Code">💻</a> <a href="#maintenance-MaCuinea" title="Maintenance">🚧</a> <a href="#tutorial-MaCuinea" title="Tutorials">✅</a> <a href="#bug-MaCuinea" title="Bug reports">🐛</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/eleonoracarpino"><img src="https://avatars.githubusercontent.com/u/123380299?v=4?s=100" width="100px;" alt="Eleonora Carpino"/><br /><sub><b>Eleonora Carpino</b></sub></a><br /><a href="#code-eleonoracarpino" title="Code">💻</a> <a href="#ideas-eleonoracarpino" title="Ideas, Planning, & Feedback">🤔</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/mekman"><img src="https://avatars.githubusercontent.com/u/139282?v=4?s=100" width="100px;" alt="Matthias Ekman"/><br /><sub><b>Matthias Ekman</b></sub></a><br /><a href="#bug-mekman" title="Bug reports">🐛</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/sirmrmarty"><img src="https://avatars.githubusercontent.com/u/140894211?v=4?s=100" width="100px;" alt="Martin Wimmers"/><br /><sub><b>Martin Wimmers</b></sub></a><br /><a href="#tutorial-sirmrmarty" title="Tutorials">✅</a> <a href="#bug-sirmrmarty" title="Bug reports">🐛</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/neurodream"><img src="https://avatars.githubusercontent.com/u/117816806?v=4?s=100" width="100px;" alt="Nico Adelhöfer"/><br /><sub><b>Nico Adelhöfer</b></sub></a><br /><a href="#ideas-neurodream" title="Ideas, Planning, & Feedback">🤔</a></td>
      <td align="center" valign="top" width="14.28%"><a href="http://lennartverhagen.com"><img src="https://avatars.githubusercontent.com/u/12236166?v=4?s=100" width="100px;" alt="Lennart Verhagen"/><br /><sub><b>Lennart Verhagen</b></sub></a><br /><a href="#ideas-lennartverhagen" title="Ideas, Planning, & Feedback">🤔</a></td>
    </tr>
  </tbody>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

# License

Released under GNU General Public License v3.0 (see LICENSE).

> **Disclaimer**
> This software is currently under development and is provided “as is” without warranty of any kind, either express or implied, including but not limited to the implied warranties of merchantability, fitness for a particular purpose, and non-infringement. This tool is intended for research purposes only and is not designed, intended, or approved for medical or clinical use, diagnosis, or treatment of patients. The contributors accept no liability for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limited to, procurement of substitute goods or services; loss of use, data, or profits; or business interruption) arising in any way out of the use of this software. Use at your own risk.