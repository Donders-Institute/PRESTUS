# Installation

PRESTUS runs in [MATLAB](https://www.mathworks.com/products/matlab.html) and simulates with [k-Wave](https://github.com/ucl-bug/k-wave) by default; head simulations also need [SimNIBS](https://simnibs.github.io/simnibs/build/html/index.html) to segment the MRI. MATLAB, PRESTUS and k-Wave are enough for a free-water simulation, and the Python environments further down serve optional features only. This page lists what each kind of simulation needs and how to install it, so you can install only what your study uses.

On this page:

- [Required Components](#required-components): which component each kind of simulation requires.
- [Install PRESTUS](#install-prestus): the clone with its submodules and the MATLAB path.
- [Install SimNIBS](#install-simnibs): the SimNIBS environment for head segmentation, and how PRESTUS finds it.
- [Install Optional Components](#install-optional-components): environments for other backends, placement tools, the Nexus, C++ binaries and defacing.

!!! note "At a glance"
    - **When you need it:** Before your [first simulation](doc_getting-started.md).
    - **Before you start:** MATLAB R2023b (other versions may work; R2023b is the tested default). Installing SimNIBS takes the most time, so start with the free-water simulation if you only want to try PRESTUS.

## Required Components

| Component | Needed for | Where it comes from |
|---|---|---|
| MATLAB R2023b | Everything | MathWorks. More recent versions can cause issues on HPCs with GPUs. |
| PRESTUS with its submodules (k-Wave 1.4.1, export_fig, FEX-minimize, xml2struct) | Everything | The recursive clone below |
| SimNIBS 4 | Head simulations (segmentation with `charm`) | [Install SimNIBS](#install-simnibs) |
| k-Wave C++ binaries | `cpp_cpu`, `cpp_gpu` and `cpp_mlx` code types | [k-Wave C++ binaries](#k-wave-c-binaries) |
| Python environments | [BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain) or [jwave](https://github.com/ucl-bug/jwave) backends, [PlanTUS](https://github.com/mlueckel/PlanTUS) and [TUS_entry](https://github.com/CyrilAtkinson/TUS_entry) placement, the Nexus | [Install Optional Components](#install-optional-components) |

## Install PRESTUS

Clone the repository with its submodules. This fetches k-Wave and the other MATLAB dependencies into `external/`:

```bash
git clone --recurse-submodules https://github.com/Donders-Institute/PRESTUS.git
# or the development branch:
git clone --recurse-submodules -b development https://github.com/Donders-Institute/PRESTUS.git
```

Then add PRESTUS to the MATLAB path at the start of every session, as `simple_main.m` does:

```matlab
addpath(genpath('functions')); addpath(genpath('config')); addpath(genpath('external'));
```

That is all a free-water simulation needs. Continue with [Your first simulation](doc_getting-started.md) and come back here for SimNIBS when you want a head simulation.

??? note "Submodules and manual downloads"
    k-Wave 1.4 is supported, but version 1.4.1 (currently [GitHub only](https://github.com/ucl-bug/k-wave/releases/tag/v1.4.1)) adds GPU support for thermal simulations with `kWaveDiffusion`, so the submodule points there. If you updated an existing clone, fetch the submodules with:

    ```bash
    cd PRESTUS
    git submodule init
    git submodule update
    ```

    Without git, download [k-Wave](https://github.com/ucl-bug/k-wave.git), [export_fig](https://github.com/altmany/export_fig), [FEX-minimize](https://github.com/rodyo/FEX-minimize.git) and [xml2struct](https://github.com/joe-of-all-trades/xml2struct) yourself and place them in the matching `external/` subfolders.

??? note "Donders Institute HPC cluster"
    On the Donders HPC, PRESTUS and its dependencies (SimNIBS, k-Wave) are installed centrally, but these versions lag behind `main` and `development`. Load a SimNIBS module (`module avail simnibs` lists them; put the `module load` line in your `.bashrc`) and add an installed PRESTUS to your MATLAB path, e.g. `addpath('/opt/prestus/dev')` (`ls /opt/prestus` lists the versions). For current features, clone PRESTUS yourself as above. The example dataset can be copied to your home directory:

    ```bash
    cp /opt/prestus/example_data/PRESTUS_example_data.zip ${HOME}
    ```

    The HPC scripts hardcode the MATLAB R2023b module. See [Running on a cluster](doc_hpc.md).

## Install SimNIBS

Head simulations need [SimNIBS 4](https://simnibs.github.io/simnibs/build/html/index.html#simnibs-4) for the `charm` segmentation. Install it in its own [conda](https://docs.conda.io/) environment; this also works on clusters that restrict user permissions and lets several SimNIBS versions coexist.

The quickest route is the setup script in `examples/`, which creates a `PRESTUS_env_4.6.0` environment with SimNIBS and the optional PlanTUS and [Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench) dependencies in one step:

```bash
bash examples/setup_PRESTUS_env.sh
```

Optional flags:

- `--no-plantus`: skip the PlanTUS Python packages (`nilearn`, `vtk`, `h5py`)
- `--no-workbench`: skip Connectome Workbench (e.g. when `wb_command` is available as an HPC module)

When it finishes, the script prints the `simnibs_bin_path` and `env_path` values to put in your config. Alternatively, follow the [official conda guide](https://simnibs.github.io/simnibs/build/html/installation/conda.html).

### Point PRESTUS to SimNIBS

Tell PRESTUS where the environment is, in the config or in MATLAB:

```matlab
parameters.startup.simnibs_bin_path = fullfile('/home', 'neuromod', 'USER', '.conda', 'envs', 'simnibs_v4.6.0', 'bin');
parameters.hpc.ld_library_path = "/opt/gcc/7.2.0/lib64";   % only if charm reports 'undefined symbol' errors
```

??? note "Manual SimNIBS installation on a Linux cluster"
    Some steps and paths differ between clusters.

    1. Log in and start an interactive session (adjust walltime and memory):

        QSUB: `qsub -I -X -N 'simnibs' -l walltime=01:00:00,mem=8g`
        [SLURM](https://slurm.schedmd.com/) (Donders default): `srun --pty --x11 -J simnibs -t 1:00:00 --mem=8G bash`

    2. Download the environment file for the SimNIBS version you want to your home directory:

        `wget https://github.com/simnibs/simnibs/releases/download/v4.6.0/environment_linux.yml -O ~/environment_linux_v4.6.0.yml`

        Other versions: replace `v4.6.0` with `v4.0.0`, `v4.1.0` or `v4.5.0`, or use `https://github.com/simnibs/simnibs/releases/latest/download/environment_linux.yml`.

    3. Load anaconda and create the environment (it lands in `~/.conda/envs/`):

        ```bash
        module load anaconda3/2020.07
        conda env create -f ~/environment_linux_v4.6.0.yml -n simnibs_v4.6.0
        source activate simnibs_v4.6.0
        ```

    4. Install SimNIBS into it:

        `pip install https://github.com/simnibs/simnibs/releases/download/v4.6.0/simnibs-4.6.0-cp311-cp311-linux_x86_64.whl`

        Other wheels: `v4.0.0` and `v4.1.0` are `cp39`, `v4.5.0` and `v4.6.0` are `cp311`; see the [releases](https://github.com/simnibs/simnibs/releases/latest).

    5. Optional: menu icons, file associations, the MATLAB library and the system path:

        ```bash
        mkdir $HOME/SimNIBS_v4.6.0
        postinstall_simnibs --setup-links -d $HOME/SimNIBS_v4.6.0
        ```

    To check the installation later: start an interactive session, `module load anaconda3/2020.07`, `source activate simnibs_v4.6.0`, `simnibs_gui`. The SimNIBS window should open. `conda deactivate` leaves the environment.

## Install Optional Components

PRESTUS itself is MATLAB, but several optional parts run in their own Python (or [R](https://www.r-project.org/)) environment. Keep them separate: their dependencies conflict, and most are only needed on the machine that runs that part.

| Environment | Needed for | Create with | Where PRESTUS finds it |
|---|---|---|---|
| SimNIBS (e.g. `simnibs_v4.6.0`, or `PRESTUS_env_4.6.0` incl. PlanTUS and Workbench) | Head segmentation (`charm`), PlanTUS placement | `examples/setup_PRESTUS_env.sh` or the SimNIBS `environment_linux.yml` (see [above](#install-simnibs)) | `simnibs_bin_path` |
| BabelBrain | `acoustic_backend` / `thermal_backend: babelbrain` | `external/BabelBrain/environment_<platform>.yml` ([BabelBrain](doc_babelbrain.md)) | `babelbrain.python_env` |
| jwave | `acoustic_backend: jwave` | `conda create` + `pip install` ([jwave](doc_jwave.md)) | `jwave.python_env`, `jwave.python_env_gpu` |
| `tus_entry` | TUS_entry placement (R) | `examples/environment_tus_entry.yml` ([Placement](doc_placement.md)) | `placement.tus_entry.env_path` |
| `prestus-nexus` | PRESTUS Nexus (config builder GUI) | `conda env create -f nexus/environment.yml` ([Nexus](doc_nexus.md#installation)) | `prestus_nexus(Python='<env>/bin/python')`; default `python3` |

The Nexus runs where you edit configs and needs a display; the others run where the pipeline runs (typically the cluster). After setting up or updating a BabelBrain or jwave environment, run `smoke_test_backends` (in `tests/`).

### k-Wave C++ Binaries

Only the `cpp_cpu` and `cpp_gpu` code types need these optional binaries: download the C++ binaries from the [k-Wave website](http://www.k-wave.org/download.php) and place them in `k-Wave/binaries/` ([k-Wave instructions](http://www.k-wave.org/documentation/kspaceFirstOrder3DC.php)). The [kspacefirstorder-unified](https://github.com/waltsims/kspacefirstorder-unified) repository is an alternative source for the same C++ code (OpenMP and CUDA variants, with macOS support). It also hosts the MLX port for Apple Silicon GPUs, which the `cpp_mlx` code type needs and which is not on the k-Wave website. Compiling and selecting binaries is described in [k-Wave](doc_kwave.md).

### PlanTUS Dependencies

The optional `plantus` placement mode needs extra Python packages in the SimNIBS environment (the setup script installs them unless `--no-plantus` was given):

```bash
source activate simnibs_v4.6.0
pip install nilearn vtk h5py pynput
python -c "import nilearn; import vtk; import h5py; import pynput; print('PlanTUS dependencies OK')"
```

PlanTUS also needs `wb_command` (Connectome Workbench). Without it PlanTUS only warns, but intermediate files are not created and the pipeline fails. Either install it into the SimNIBS environment, where PRESTUS finds it in the environment's `bin/` without further settings:

```bash
conda install -c conda-forge connectome-workbench
```

or, on a cluster where Workbench is a module that MATLAB's Python subprocess does not see, point PRESTUS at its binary folder:

```yaml
placement:
  plantus:
    connectome_wb_path: "/path/to/workbench/bin_linux64"
```

### MRI Defacing

Optional MRI defacing (`segmentation.deface: 1`) needs [pydeface](https://github.com/poldracklab/pydeface) in the SimNIBS environment:

```bash
source activate simnibs_v4.6.0
pip install pydeface
python -c "import pydeface; print('pydeface OK')"
```
