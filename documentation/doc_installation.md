# Installation Guide

### PRESTUS installation

Download and install these tools:

- MATLAB (R2022b*)
- [SimNIBS 4.0.0](https://github.com/simnibs/simnibs) (see instructions below)
- toolboxes | They are automatically included as submodules if this repository is recursively cloned (see below). They must be added on MATLAB startup.
    - [k-Wave (1.4)](https://github.com/ucl-bug/k-wave.git)
    - [export_fig](https://github.com/altmany/export_fig)
    - [FEX-minimize](https://github.com/rodyo/FEX-minimize.git)
    - [xml2struct](https://github.com/joe-of-all-trades/xml2struct)

\* Other versions may work as well, but R2022b was the default deployment. Especially on HPCs with GPUs, more recent MATLAB versions can cause issues.

Before using the package, you need to have some libraries on your path. 
Major dependencies are included as submodules in the toolbox folder. 
If you clone this repository, you can retrieve the submodules as follows:
```
git clone --recurse-submodules https://github.com/Donders-Institute/PRESTUS.git
```

If you want to recursively clone the `development` branch, you can use the following command:
```
git clone --recurse-submodules -b development https://github.com/Donders-Institute/PRESTUS.git
```

*Note: If you do not clone this repository, you must manually download the toolboxes specified above, and place them in the respective PRESTUS subdirectory.*

If you cloned this repository in the past, and updated it, you can retrieve submodules as follows:
```
cd PRESTUS
git submodule init
git submodule update
```

### SimNIBS installation

You additionally need to install SimNIBS (https://simnibs.github.io/simnibs/build/html/index.html#simnibs-4): e.g., `simnibs_installer/install -s -t /home/USER/SimNIBS`.
We recommend installing SimNIBS within an **anaconda environment** (especially on HPC clusters that often constrict individual user permissions).
You can follow the [official installation guide](https://simnibs.github.io/simnibs/build/html/installation/conda.html).

The following is a step-by-step guide to install SIMNIBS on the computing cluster. Some steps and filepaths may differ between computing environments.

1. Log in to the cluster
2. Open terminal
3. Start an interactive session to come into another execution node (change walltime and memory as needed): `qsub -I -X -N 'simnibs' -l walltime=01:00:00,mem=8g`
4. Download the evironment for SIMNIBS: `wget https://github.com/simnibs/simnibs/releases/latest/download/environment_linux.yml`
5. Copy the generated environment file to your home directory (i.e., /home/neuromod/USER)
6. Navigate to your home directory (i.e., `cd ~`)
7. Load the anaconda environment (e.g., `module load anaconda3/2020.07`)
8. Create the anaconda environment: `conda env create -f ~/environment_linux.yml`
9. Activate the simnibs environment: `source activate simnibs_env`
10. Install the latest version of SIMNIBS: `pip install -f https://github.com/simnibs/simnibs/releases/latest simnibs`
11. [Optional] To setup menu icons, file association, the MATLAB library, and add SimNIBS to the system path, run the postinstall_simnibs script

```
mkdir $HOME/SimNIBS
postinstall_simnibs --setup-links -d $HOME/SimNIBS
```

#### [Optional: Starting SimNIBS after first installation]

Exit out of everything and following the steps below for starting SimNIBS after first installation to check whether installation was successful. If installation was successful you should now be able to see the SimNIBS gui.

1. Start an interactive session, e.g., `qsub -I -X -N 'simnibs' -l walltime=01:00:00,mem=8g`
2. `module load anaconda3/2020.07`
3. `source activate simnibs_env`
4. `simnibs_gui`
5. To close the virtual SimNIBS environment: conda deactivate

#### Specifying SimNIBS paths in PRESTUS

In PRESTUS, specify either in the config or directly in MATLAB both the path to the SimNIBS binaries, and the path to shared libraries for GCC 7.2.0:
```
parameters.simnibs_bin_path = fullfile('/home', 'neuromod', 'USER', '.conda', 'envs', 'simnibs_env', 'bin');
parameters.ld_library_path = "/opt/gcc/7.2.0/lib64";
```