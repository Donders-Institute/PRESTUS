# PRESTUS: PREprocessing & Simulations for Transcranial Ultrasound Stimulation package

PRESTUS provides a set of functions for ultrasonic simulations including brain preprocessing, extracting data from Localite neuronavigation files, and more. 

The examples folder contains some examples of how to run the pipeline, including a tutorial.

# Installation

## DCC/N Cluster

When working on the DCCN cluster, PRESTUS and its dependencies (SimNIBS, k-Wave) are already installed. 

Type ``module load simnibs/4.0.0`` (or add the command to your .bashrc so that it is executed atomatically once you login) and add ``addpath('/opt/prestus/dev')`` to your matlab path*. Now you can start matlab R2022b.

If you want to get started with simulations, you can use the PRESTUS example dataset. This command will copy the dataset to your home directory:

```
cp /opt/prestus/example_data/PRESTUS_example_data.zip ${HOME}
```

*this will use the most up-to-date version of PRESTUS (i.e., the current development branch). If you want to use an older version you can also use ``addpath('/opt/prestus/0.2.0')``, or older versions.

## Outside the DCC/N

Download and install these tools:

- [SimNIBS 4](https://simnibs.github.io/simnibs/build/html/index.html#simnibs-4)
- [k-Wave 1.4](http://www.k-wave.org/download.php) (should be installed in the 'toolbox' folder of PRESTUS or be added automatically on MATLAB startup if you use HPC)
- toolboxes included in this repository

Example data: We are currently working on a solution to make the example dataset available for users outside the DCCN.

# Authorship

The initial version was developed by Andrey Chetverikov (website: http://andreychetverikov.org, GitLab: @A.Chetverikov, GitHub: https://github.com/achetverikov).

# Contributors

- Kenneth van der Zee (GitLab: @kenneth.vanderzee, GitHub: https://github.com/KTZ228)
- Julian Kosciessa (GitLab: @julian-kosciessa, GitHub: https://github.com/jkosciessa)
- Matthias Ekman (GitLab: @m.ekman, GitHub: https://github.com/mekman)
- Eleonora Carpino (GitLab: @eleonora.carpino, GitHub: https://github.com/eleonoracarpino)

# License

Released under GNU General Public License v3.0 (see LICENSE).
