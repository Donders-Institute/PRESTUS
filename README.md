[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15095860.svg)](https://doi.org/10.5281/zenodo.15095860)

# PRESTUS: PREprocessing & Simulations for Transcranial Ultrasound Stimulation

PRESTUS is an open-source toolbox that facilitates imaging-informed TUS simulations. PRESTUS provides an integrated solution for defining relevant tissue parameters based on personalized MR or CT images, and using them to run  acoustic and heating simulations (with HCP support), and summarizing outputs in standard template space. These features make PRESTUS an accessible tool for researchers to plan, validate, and report transcranial ultrasound targeting.

The examples folder contains some examples of how to run the pipeline, including a tutorial.

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

Chetverikov, A., Kosciessa, J. Q., van der Zee, K., & Verhagen, L. (2024). PRESTUS (0.3.0). Zenodo. https://doi.org/10.5281/zenodo.15095861

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
