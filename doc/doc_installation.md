# Installation Guide

## PRESTUS installation

Download and install these tools:

- `MATLAB (R2023b)`. Other versions may work as well, but R2023b is the current default deployment. Especially on HPCs with GPUs, more recent MATLAB versions can cause issues. The HPC scripts currently hardcode the MATLAB R2023b module on the Donders HPC.
- [`SimNIBS 4`](https://github.com/simnibs/simnibs) (see [SimNIBS installation](#simnibs-installation))
- external | They are automatically included as submodules if this repository is recursively cloned (see below). They must be added on MATLAB startup.
    - [`k-Wave (1.4.1)`](https://github.com/ucl-bug/k-wave.git). While k-Wave 1.4 is supported, version 1.4.1 (currently [GitHub exclusive](https://github.com/ucl-bug/k-wave/releases/tag/v1.4.1)) introduced GPU support for thermal simulations with kWaveDiffusion. We therefore recommend cloning k-Wvae 1.4.1 from GitHub.
    - [`export_fig`](https://github.com/altmany/export_fig)
    - [`FEX-minimize`](https://github.com/rodyo/FEX-minimize.git)
    - [`xml2struct`](https://github.com/joe-of-all-trades/xml2struct)

Major dependencies are included as submodules in the toolbox folder. 
If you clone this repository, you can retrieve the submodules as follows:
```
git clone --recurse-submodules https://github.com/Donders-Institute/PRESTUS.git
```

If you want to recursively clone the `development` branch, you can use the following command:
```
git clone --recurse-submodules -b development https://github.com/Donders-Institute/PRESTUS.git
```

*Note: If you do not clone this repository, you must manually download the external dependencies listed above, and place them in the respective PRESTUS subdirectory.*

If you cloned this repository in the past, and updated it, you can retrieve submodules as follows:
```
cd PRESTUS
git submodule init
git submodule update
```

Ensure that the paths and subfolders are added in MATLAB. See `simple_main.m` for an example.

#### Donders Institute HPC Cluster

When working on the Donders High-Performance-Computing cluster, PRESTUS and its dependencies (SimNIBS, k-Wave) are already installed. Note that this may not be the most recent version of either the `main` or `development` branch.

Type `module load simnibs/4.0.0` (or add the command to your .bashrc so that it is executed automatically once you login) and add `addpath('/opt/prestus/dev')` to your matlab path. If you want to use an older version you can also use `addpath('/opt/prestus/0.2.0')`, or older versions. Now you can start matlab R2022b.

> `simnibs/4.0.0` and `simnibs/4.1.0` are currently available on the Donders HPC.

> `/opt/prestus/0.1.0`, `/opt/prestus/0.2.0`, and (an outdated) `/opt/prestus/dev` are currently available on the Donders HPC.

For more information on HPC usage, see the [HPC guide](doc_hpc.md)

If you want to get started with simulations, you can use the PRESTUS example dataset. This command will copy the dataset to your home directory:

```
cp /opt/prestus/example_data/PRESTUS_example_data.zip ${HOME}
```

#### [Optional: Download C++ binaries]

For the `cpp_cpu` and `cpp_gpu` modes, download the C++ binaries from the [k-Wave website](http://www.k-wave.org/download.php) and place them into `k-Wave/binaries/`; see also the instructions [here](http://www.k-wave.org/documentation/kspaceFirstOrder3DC.php).

For more information on installing, (potentially) compiling, and specifying C++ binaries, see [the backend documentation](doc_backend.md).

## SimNIBS installation

You additionally need to install [SimNIBS 4](https://simnibs.github.io/simnibs/build/html/index.html#simnibs-4). We recommend installing SimNIBS within an **anaconda environment** (especially on HPC clusters that often constrict individual user permissions). This allows for the flexible installation of multiple SimNIBS versions.
You can follow the [official installation guide](https://simnibs.github.io/simnibs/build/html/installation/conda.html).

The following is a step-by-step guide to install SIMNIBS on the Linus computing cluster. Some steps and filepaths may differ between computing environments.

1. Log in to the cluster
2. Open terminal
3. Start an interactive session to come into another execution node (change walltime and memory as needed)  
    
    QSUB: `qsub -I -X -N 'simnibs' -l walltime=01:00:00,mem=8g` 
    SLURM (Donders default): `srun --pty --x11 -J simnibs -t 1:00:00 --mem=8G bash` 

4. Download the environment for the desired SIMNIBS version to your home directory (i.e., `/home/neuromod/USER`)

    `wget https://github.com/simnibs/simnibs/releases/download/v4.6.0/environment_linux.yml -O ~/environment_linux_v4.6.0.yml`. This includes dependencies such as Python, MATLAB Runtime, and other prerequisites.

    <br>

    > **SimNIBS environment versions**
    Multiple SimNIBS 4+ environment versions are available: 
    https://github.com/simnibs/simnibs/releases/download/v4.0.0/environment_linux.yml   
    https://github.com/simnibs/simnibs/releases/download/v4.1.0/environment_linux.yml   
    https://github.com/simnibs/simnibs/releases/download/v4.5.0/environment_linux.yml   
    https://github.com/simnibs/simnibs/releases/download/v4.6.0/environment_linux.yml   
    To download the latest environment, select the following:   
    https://github.com/simnibs/simnibs/releases/latest/download/environment_linux.yml   

    <br>

5. Navigate to your home directory (i.e., `cd ~`)
6. Load the anaconda environment    

    e.g., `module load anaconda3/2020.07`   

7. Create the anaconda environment

    `conda env create -f ~/environment_linux_v4.6.0.yml -n simnibs_v4.6.0`
    
    This will create an environment in `~/.conda/envs/`.    

8. Activate the simnibs environment 

    `source activate simnibs_v4.6.0`    

9. Install the desired SimNIBS version

    `pip install https://github.com/simnibs/simnibs/releases/download/v4.6.0/simnibs-4.6.0-cp311-cp311-linux_x86_64.whl` 
    
    <br>
    
    > **SimNIBS versions**  
    Multiple SimNIBS 4+ LINUX versions are available:   
    https://github.com/simnibs/simnibs/releases/download/v4.0.0/simnibs-4.0.0-cp39-cp39-linux_x86_64.whl    
    https://github.com/simnibs/simnibs/releases/download/v4.1.0/simnibs-4.1.0-cp39-cp39-linux_x86_64.whl    
    https://github.com/simnibs/simnibs/releases/download/v4.5.0/simnibs-4.5.0-cp311-cp311-linux_x86_64.whl  
    https://github.com/simnibs/simnibs/releases/download/v4.6.0/simnibs-4.6.0-cp311-cp311-linux_x86_64.whl  
    To download the latest environment, select the following:   
    https://github.com/simnibs/simnibs/releases/latest  

    <br>

10. [Optional] To setup menu icons, file association, the MATLAB library, and add SimNIBS to the system path, run the postinstall_simnibs script

    ```
    mkdir $HOME/SimNIBS_v4.6.0
    postinstall_simnibs --setup-links -d $HOME/SimNIBS_v4.6.0
    ```

#### [Optional: Starting SimNIBS after first installation]

Exit out of everything and following the steps below for starting SimNIBS after first installation to check whether installation was successful. If installation was successful you should now be able to see the SimNIBS gui.

1. Start an interactive session (see Step 3 [above](#simnibs-installation))
2. `module load anaconda3/2020.07`
3. `source activate simnibs_v4.6.0`
4. `simnibs_gui`
5. To close the virtual SimNIBS environment: `conda deactivate`

#### Specifying SimNIBS paths in PRESTUS

In PRESTUS, specify either in the config or directly in MATLAB both the path to the SimNIBS binaries, and the path to shared libraries for GCC 7.2.0:
```
parameters.simnibs_bin_path = fullfile('/home', 'neuromod', 'USER', '.conda', 'envs', 'simnibs_v4.6.0', 'bin');
parameters.ld_library_path = "/opt/gcc/7.2.0/lib64";
```