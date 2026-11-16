# Installation

PRESTUS runs in [MATLAB](https://www.mathworks.com/products/matlab.html) and simulates with [k-Wave](https://github.com/ucl-bug/k-wave) by default; head simulations also need [SimNIBS](https://simnibs.github.io/simnibs/build/html/index.html) to segment the MRI. MATLAB, PRESTUS and k-Wave are enough for a free-water simulation, and the Python environments further down serve optional features only. This page lists what each kind of simulation needs and how to install it, so you can install only what your study uses.

On this page:

- [Quick Install](#quick-install): three steps to a working free-water setup, and the Nexus.
- [Required Components](#required-components): which component each kind of simulation requires.
- [Install PRESTUS](#install-prestus): the clone with its submodules and the MATLAB path.
- [Install SimNIBS](#install-simnibs): the SimNIBS environment for head segmentation, and how PRESTUS finds it.
- [Install Optional Components](#install-optional-components): environments for other backends, placement tools, the Nexus, C++ binaries and defacing.

!!! note "At a glance"
    - **When you need it:** Before your [first simulation](doc_getting-started.md).
    - **Before you start:** MATLAB R2023b (other versions may work; R2023b is the tested default). Installing SimNIBS takes the most time, so start with the free-water simulation if you only want to try PRESTUS.

## Quick Install

1. Clone PRESTUS (or download it as a zip from GitHub):

    ```bash
    git clone https://github.com/Donders-Institute/PRESTUS.git
    ```

2. In MATLAB, go to the PRESTUS folder and run:

    ```matlab
    prestus_setup(Save=true)
    ```

    This puts PRESTUS on the MATLAB path for this and later sessions, fetches the required submodules (k-Wave, export_fig, FEX-minimize, xml2struct) with git, and prints what is installed and what is missing, with how to fix each item. Run `prestus_check_install` at any time to see the same report.

3. Optional, for the [PRESTUS Nexus](doc_nexus.md): install its Python environment once.

    | System | Run |
    |---|---|
    | macOS, Linux | `bash nexus/install.sh` in the PRESTUS folder |
    | Windows | double-click `nexus\install.bat` |

    The installer uses conda if you have it and otherwise offers to install [Miniforge](https://github.com/conda-forge/miniforge). Then `prestus_nexus` in MATLAB finds the environment by itself.

That is enough for a free-water simulation: continue with [Your first simulation](doc_getting-started.md). Head simulations also need [SimNIBS](#install-simnibs).

## Required Components

| Component | Needed for | Where it comes from |
|---|---|---|
| MATLAB R2023b | Everything | MathWorks. More recent versions can cause issues on HPCs with GPUs. |
| PRESTUS with its submodules (k-Wave 1.4.1, export_fig, FEX-minimize, xml2struct) | Everything | The recursive clone below |
| SimNIBS 4 | Head simulations (segmentation with `charm`) | [Install SimNIBS](#install-simnibs) |
| k-Wave C++ binaries | `cpp_cpu`, `cpp_gpu` and `cpp_mlx` code types | [k-Wave C++ binaries](#k-wave-c-binaries) |
| Python environments | [BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain) or [jwave](https://github.com/ucl-bug/jwave) backends, [PlanTUS](https://github.com/mlueckel/PlanTUS) and [TUS_entry](https://github.com/CyrilAtkinson/TUS_entry) placement, the Nexus | [Install Optional Components](#install-optional-components) |

## Install PRESTUS

Clone the repository and run `prestus_setup` in MATLAB, as in the [Quick Install](#quick-install). `prestus_setup` fetches only the required submodules; `--recurse-submodules` also fetches the optional ones (BabelBrain, PlanTUS):

```bash
git clone --recurse-submodules https://github.com/Donders-Institute/PRESTUS.git
# or the development branch:
git clone --recurse-submodules -b development https://github.com/Donders-Institute/PRESTUS.git
```

`prestus_setup` adds `functions/`, `config/` and `external/` to the MATLAB path from any folder (`run('/path/to/PRESTUS/prestus_setup')`). Without `Save=true` the path holds for the current session only. If MATLAB cannot save the path (no write access to `pathdef.m`, common on clusters), add the `run(...)` line to your `startup.m`.

??? note "Submodules and manual downloads"
    k-Wave 1.4 is supported, but version 1.4.1 (currently [GitHub only](https://github.com/ucl-bug/k-wave/releases/tag/v1.4.1)) adds GPU support for thermal simulations with `kWaveDiffusion`, so the submodule points there. If you updated an existing clone, `prestus_setup` fetches missing required submodules; to update all of them:

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

Head simulations need [SimNIBS 4](https://simnibs.github.io/simnibs/build/html/index.html#simnibs-4) for the `charm` segmentation.

On macOS and Windows, use the [SimNIBS installer](https://simnibs.github.io/simnibs/build/html/installation/simnibs_installer.html). `prestus_check_install` then reports the bin folder it installed `charm` into, as the value for `startup.simnibs_bin_path`.

On Linux and clusters, install SimNIBS in its own [conda](https://docs.conda.io/) environment; this also works on clusters that restrict user permissions and lets several SimNIBS versions coexist. The quickest route is the setup script in `examples/`, which creates a `PRESTUS_env_4.6.0` environment with SimNIBS and the optional PlanTUS and [Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench) dependencies in one step:

```bash
bash examples/setup_PRESTUS_env.sh
```

Optional flags:

- `--no-plantus`: skip the PlanTUS Python packages (`nilearn`, `vtk`, `h5py`)
- `--no-workbench`: skip Connectome Workbench (e.g. when `wb_command` is available as an HPC module)

When it finishes, the script prints the `simnibs_bin_path` value to put in your config; `hpc.ld_library_path` and the PlanTUS `env_path` follow it while they are empty. Alternatively, follow the [official conda guide](https://simnibs.github.io/simnibs/build/html/installation/conda.html).

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
| `prestus-nexus` | PRESTUS Nexus (config builder GUI) | `nexus/install.sh` or `nexus/install.bat` ([Nexus](doc_nexus.md#installation)) | found by `prestus_nexus`; override with `Python='<env>/bin/python'` |

The Nexus runs where you edit configs and needs a display; the others run where the pipeline runs (typically the cluster). After setting up or updating a BabelBrain or jwave environment, run `smoke_test_backends` (in `tests/`).

### k-Wave C++ Binaries

Only the `cpp_cpu` and `cpp_gpu` code types need these optional binaries: download the C++ binaries from the [k-Wave website](http://www.k-wave.org/download.php) and place them in `k-Wave/binaries/` ([k-Wave instructions](http://www.k-wave.org/documentation/kspaceFirstOrder3DC.php)). The [kspacefirstorder-unified](https://github.com/waltsims/kspacefirstorder-unified) repository is an alternative source for the same C++ code (OpenMP and CUDA variants, with macOS support). It also hosts the MLX port for Apple Silicon GPUs, which the `cpp_mlx` code type needs and which is not on the k-Wave website. Compiling and selecting binaries is described in [k-Wave](doc_kwave.md).

### PlanTUS Dependencies

The optional `plantus` placement mode needs [PlanTUS](https://github.com/mlueckel/PlanTUS) 2.0 or later (`git clone https://github.com/mlueckel/PlanTUS`, then `placement.plantus.script_path`) and extra Python packages in the SimNIBS environment (the setup script installs them unless `--no-plantus` was given):

```bash
source activate simnibs_v4.6.0
pip install nilearn vtk h5py trimesh pynput
python -c "import nilearn, vtk, h5py, trimesh, pynput; print('PlanTUS dependencies OK')"
```

Connectome Workbench is needed only to choose the position in Workbench (`wb_view`), and with PlanTUS 1.x, which also calls `wb_command`. Either install it into the SimNIBS environment, where PRESTUS finds it in the environment's `bin/` without further settings:

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
