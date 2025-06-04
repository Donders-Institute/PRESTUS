[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15095860.svg)](https://doi.org/10.5281/zenodo.15095860)

# PRESTUS: PREprocessing & Simulations for Transcranial Ultrasound Stimulation

PRESTUS (PREprocessing & Simulations for Transcranial Ultrasound Stimulation) is an open-source MATLAB toolbox that aims to streamline imaging-informed TUS simulations: from the segmentation of T1-weighted MRI head scans and mapping of medium tissue properties (possibly informed for skull via (pseudo-)CT images) in a simulation grid, to the execution of acoustic and thermal simulations using the widely adopted k-Wave engine. High-performance computing (HPC) support via SLURM and CentOS enables efficient parallelization, and large-scale analyses. Output 3D NifTI images are automatically mapped to a standard template (MNI) space to facilitate group reporting.

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

If you would like to assist the community development of this tool, please consider [CONTRIBUTING](CONTRIBUTING.md).

# Installation

## DCC/N Cluster

When working on the DCCN cluster, PRESTUS and its dependencies (SimNIBS, k-Wave) are already installed. 

Type ``module load simnibs/4.0.0`` (or add the command to your .bashrc so that it is executed automatically once you login) and add ``addpath('/opt/prestus/dev')`` to your matlab path*. Now you can start matlab R2022b.

If you want to get started with simulations, you can use the PRESTUS example dataset. This command will copy the dataset to your home directory:

```
cp /opt/prestus/example_data/PRESTUS_example_data.zip ${HOME}
```

*this will use the most up-to-date version of PRESTUS (i.e., the current development branch). If you want to use an older version you can also use ``addpath('/opt/prestus/0.2.0')``, or older versions.

## Outside the DCC/N

Download and install these tools:

- MATLAB (R2022b)
- [SimNIBS 4](https://github.com/simnibs/simnibs)
- toolboxes included as submodules if this repository is cloned (must be added on MATLAB startup)
    - [k-Wave (1.4)](https://github.com/ucl-bug/k-wave.git)
    - [export_fig](https://github.com/altmany/export_fig)
    - [FEX-minimize](https://github.com/rodyo/FEX-minimize.git)
    - [xml2struct](https://github.com/joe-of-all-trades/xml2struct)

Tested on MATLAB 2022b, set up to work on Donders HPC (can work on a local PC as well). 

Before using the package, you need to have some libraries on your path. Major dependencies are included as submodules in the toolbox folder. If you clone this repository, you can retrieve the submodules as follows:
```
git clone --recurse-submodules https://github.com/Donders-Institute/PRESTUS.git
```

If you cloned this repository in the past, and updated it, you can retrieve submodules as follows:
```
cd PRESTUS
git submodule init
git submodule update
```
You additionally need to install SimNIBS (https://simnibs.github.io/simnibs/build/html/index.html#simnibs-4).

*Note: If you do not clone this repository, you must manually download and add the toolboxes specified above.*

Example data: We are currently working on a solution to make the example dataset available for users outside the DCCN.

# Reference

Chetverikov, A., Kosciessa, J. Q., Cornelissen, M., van der Zee, K., & Verhagen, L. (2024). PRESTUS (0.3.0). Zenodo. https://doi.org/10.5281/zenodo.15095861

# Contributors

- Andrey Chetverikov (@achetverikov)
- Julian Kosciessa (@jkosciessa)
- Kenneth van der Zee (@KTZ228)
- Margely Cornelissen (@MaCuinea)
- Matthias Ekman (@mekman)
- Eleonora Carpino (@eleonoracarpino)
- Martin Wimmers (@sirmrmarty)

# License

Released under GNU General Public License v3.0 (see LICENSE).

> **Disclaimer**
> This software is currently under development and is provided “as is” without warranty of any kind, either express or implied, including but not limited to the implied warranties of merchantability, fitness for a particular purpose, and non-infringement. This tool is intended for research purposes only and is not designed, intended, or approved for medical or clinical use, diagnosis, or treatment of patients. The contributors accept no liability for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limited to, procurement of substitute goods or services; loss of use, data, or profits; or business interruption) arising in any way out of the use of this software. Use at your own risk.