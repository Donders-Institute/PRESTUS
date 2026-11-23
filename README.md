[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15095860.svg)](https://doi.org/10.5281/zenodo.15095860)

# PRESTUS: PREprocessing & Simulations for Transcranial Ultrasound Stimulation

PRESTUS (PREprocessing & Simulations for Transcranial Ultrasound Stimulation) is an open-source MATLAB toolbox that aims to streamline imaging-informed simulations of Transcranial Ultrasound Stimulation (TUS): from the segmentation of T1-weighted MRI head scans and mapping of medium tissue properties (possibly informed for skull via (pseudo-)CT images) in a simulation grid, to the execution of acoustic and thermal simulations using the widely adopted k-Wave engine. High-performance computing (HPC) support via SLURM and CentOS enables efficient parallelization, and large-scale analyses. Output 3D NifTI images are automatically mapped to a standard template (MNI) space to facilitate group reporting.

Key features include:

- Automated MRI segmentation (using SimNIBS 4 charm) and preprocessing.
- 2D / 3D grid setup.
- Multi-layer medium property mapping (water, skin, multi-layer skull, and brain).
- (pseudo-)CT-informed continuous skull mapping.
- Virtual multi-element transducer calibration (free-water profile emulation).
- Estimation of entry-target coordinates.
- Flexible temporal protocol specification (e.g., including session breaks).
- k-Wave integration for robust acoustic and heating simulations, with [BabelBrain and jwave as alternative backends](https://donders-institute.github.io/PRESTUS/doc_backend/).
- Support for high-performance computing.
- 3D NifTI outputs for major reporting metrics (in subject- & MNI-space).

These features make PRESTUS an accessible tool for researchers to plan, validate, and report transcranial ultrasound targeting. 
A recent overview poster can be found [here](https://jkosciessa.github.io/downloads/2025-FUN25-PRESTUS.pdf).

# Installation

```bash
git clone https://github.com/Donders-Institute/PRESTUS.git
```

Then, in MATLAB in the PRESTUS folder, run `prestus_setup(Save=true)`. It puts PRESTUS on the path, fetches the required submodules and reports what else is installed. For the PRESTUS Nexus, run `bash nexus/install.sh` (macOS, Linux) or `nexus\install.bat` (Windows). Head simulations also need SimNIBS 4. See [Installation](https://donders-institute.github.io/PRESTUS/doc_installation/) for details and optional components.

# Documentation

Please see the [PRESTUS documentation](https://donders-institute.github.io/PRESTUS/).

# Reference

If you use this tool for publications, please cite:
Kosciessa, J. Q., Chetverikov, A., Cornelissen, M., van der Zee, K., Carpino, E., Wimmers, M., Adelhöfer, N., Meijer, S., & Verhagen, L. (2026). PRESTUS (0.6.0). Zenodo. https://doi.org/10.5281/zenodo.15095860

# Contributors

[![Contributors](https://img.shields.io/github/contributors/Donders-Institute/PRESTUS.svg?color=00B4D8&style=flat-square)](https://github.com/Donders-Institute/PRESTUS/graphs/contributors)
<p>
  <a href="http://juliankosciessa.eu" title="Julian Kosciessa: code, ideas, tutorials, maintenance, bug reports, docs"><img src="https://avatars.githubusercontent.com/u/40263608?v=4&s=112" width="56" height="56" alt="Julian Kosciessa"/></a>
  <a href="http://andreychetverikov.org" title="Andrey Chetverikov: code, ideas, tutorials"><img src="https://avatars.githubusercontent.com/u/1465806?v=4&s=112" width="56" height="56" alt="Andrey Chetverikov"/></a>
  <a href="https://github.com/KTZ228" title="Kenneth van der Zee: code, ideas, maintenance, bug reports"><img src="https://avatars.githubusercontent.com/u/51954604?v=4&s=112" width="56" height="56" alt="Kenneth van der Zee"/></a>
  <a href="https://github.com/MaCuinea" title="Margely Cornelissen: code, maintenance, tutorials, bug reports"><img src="https://avatars.githubusercontent.com/u/134381864?v=4&s=112" width="56" height="56" alt="Margely Cornelissen"/></a>
  <a href="https://github.com/eleonoracarpino" title="Eleonora Carpino: code, ideas"><img src="https://avatars.githubusercontent.com/u/123380299?v=4&s=112" width="56" height="56" alt="Eleonora Carpino"/></a>
  <a href="https://github.com/sirmrmarty" title="Martin Wimmers: code, ideas, tutorials, maintenance, bug reports"><img src="https://avatars.githubusercontent.com/u/140894211?v=4&s=112" width="56" height="56" alt="Martin Wimmers"/></a>
  <a href="https://github.com/neurodream" title="Nico Adelhöfer: code, ideas"><img src="https://avatars.githubusercontent.com/u/117816806?v=4&s=112" width="56" height="56" alt="Nico Adelhöfer"/></a>
  <a href="https://github.com/meijer-s" title="Sjoerd Meijer: code, ideas"><img src="https://avatars.githubusercontent.com/u/243970592?v=4&s=112" width="56" height="56" alt="Sjoerd Meijer"/></a>
  <a href="https://github.com/Marwan-Engels" title="Marwan Engels: code, ideas, docs"><img src="https://avatars.githubusercontent.com/u/70513219?v=4&s=112" width="56" height="56" alt="Marwan Engels"/></a>
  <a href="https://github.com/mekman" title="Matthias Ekman: bug reports"><img src="https://avatars.githubusercontent.com/u/139282?v=4&s=112" width="56" height="56" alt="Matthias Ekman"/></a>
  <a href="http://lennartverhagen.com" title="Lennart Verhagen: ideas"><img src="https://avatars.githubusercontent.com/u/12236166?v=4&s=112" width="56" height="56" alt="Lennart Verhagen"/></a>
</p>

- [Julian Kosciessa](http://juliankosciessa.eu): code, ideas, tutorials, maintenance, bug reports, docs
- [Andrey Chetverikov](http://andreychetverikov.org): code, ideas, tutorials
- [Kenneth van der Zee](https://github.com/KTZ228): code, ideas, maintenance, bug reports
- [Margely Cornelissen](https://github.com/MaCuinea): code, maintenance, tutorials, bug reports
- [Eleonora Carpino](https://github.com/eleonoracarpino): code, ideas
- [Martin Wimmers](https://github.com/sirmrmarty): code, ideas, tutorials, maintenance, bug reports
- [Nico Adelhöfer](https://github.com/neurodream): code, ideas
- [Sjoerd Meijer](https://github.com/meijer-s): code, ideas
- [Marwan Engels](https://github.com/Marwan-Engels): code, ideas, docs
- [Matthias Ekman](https://github.com/mekman): bug reports
- [Lennart Verhagen](http://lennartverhagen.com): ideas

If you would like to assist the community development of this tool, please consider [CONTRIBUTING](CONTRIBUTING.md).

# License

Released under GNU General Public License v3.0 (see LICENSE).

> **Disclaimer**
> This software is currently under development and is provided “as is” without warranty of any kind, either express or implied, including but not limited to the implied warranties of merchantability, fitness for a particular purpose, and non-infringement. This tool is intended for research purposes only and is not designed, intended, or approved for medical or clinical use, diagnosis, or treatment of patients. The contributors accept no liability for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limited to, procurement of substitute goods or services; loss of use, data, or profits; or business interruption) arising in any way out of the use of this software. Use at your own risk.