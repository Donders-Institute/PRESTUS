# jwave Acoustic Backend

[jwave](https://github.com/ucl-bug/jwave) implements the [k-Wave](https://github.com/ucl-bug/k-wave) acoustic equations independently, in [JAX](https://github.com/jax-ml/jax), so a jwave run of the same setup and tissue properties cross-checks a k-Wave field. Its default Helmholtz solver computes the continuous-wave steady state directly, in one linear solve instead of thousands of time steps. jwave covers the acoustic stage only; the thermal simulation runs in k-Wave.

On this page:

- [Quick Start](#quick-start): selecting the backend, the Python environment, and the solver.
- [Python Environment Setup](#python-environment-setup): installing jwave and JAX for CPU or GPU.
- [Limitations](#limitations): what jwave does not support.

!!! warning "Experimental"
    The jwave backend is not yet validated against measurements or the [ITRUSST](https://itrusst.org/) benchmarks in PRESTUS. The conservative Helmholtz operator became the default in October 2026 ([Helmholtz operator](doc_jwave_configuration.md)). Use jwave to cross-check k-Wave, not to replace it.

!!! note "At a glance"
    - **When you need it:** To cross-check k-Wave with an independent implementation, or for a fast CW steady-state solve.
    - **Before you start:** A [Python](https://www.python.org/) 3.11 environment with `jax`, `jwave`, `h5py`, and `numpy` ([Python environment setup](#python-environment-setup)).

PRESTUS can run [jwave](https://github.com/ucl-bug/jwave) in place of k-Wave for the acoustic stage. Pre-processing (segmentation, grid, medium masks, placement) and post-processing (focal metrics, [NIfTI](https://nifti.nimh.nih.gov/) export, HTML report) run unchanged. jwave (University College London) implements the k-Wave pseudo-spectral equations in [JAX](https://github.com/google/jax), which gives JIT compilation, automatic differentiation, and NVIDIA GPU acceleration without compiled [MATLAB](https://www.mathworks.com/products/matlab.html)/C++ binaries. It offers a time-domain solver equivalent to `kspaceFirstOrder*` and a frequency-domain Helmholtz solver that computes the CW steady state directly.

| | PRESTUS (k-Wave) | jwave |
|---|---|---|
| **Solver** | Pseudo-spectral FDTD | Pseudo-spectral (time-domain) or Helmholtz GMRES (frequency-domain) |
| **Physics** | Scalar linear acoustic | Scalar linear acoustic |
| **Skull model** | Scalar ρ, c, α | Scalar ρ, c, α |
| **Differentiable** | No | Yes (JAX automatic differentiation) |
| **GPU backends** | [CUDA](https://developer.nvidia.com/cuda-toolkit) (via MATLAB `gpuArray` or C++ binary) | NVIDIA CUDA only |
| **Thermal simulation** | Yes (`kWaveDiffusion`) | Not supported |

jwave uses the same scalar acoustic model as k-Wave: no shear waves in bone and no Rayleigh pre-step. The `medium_properties` PRESTUS assigns apply to both equally, so results are directly comparable. See [k-Wave](doc_kwave.md) for the k-Wave backend options.

## Quick Start

Set three things: the backend, the Python environment, and the solver.

```yaml
acoustic_backend: jwave

jwave:
  python_env: '/home/user/.conda/envs/jwave_env/bin/python'   # CPU environment (see environment setup)
  solver: helmholtz          # 'helmholtz' (default) | 'pseudo_spectral' | 'finite_differences'
  python_env_gpu: ''         # GPU (CUDA) environment; leave empty for CPU-only
  use_gpu: ~                 # null = auto-detect via nvidia-smi | true = force GPU | false = force CPU
  dry_run: false             # write HDF5 inputs only, skip simulation (for debugging)

# PML thickness uses the same grid.pml_size field as k-Wave (default: 'auto'),
# but jwave resolves 'auto' to ~2.5 wavelengths of PML (clipped to [10, 40]
# grid points), not k-Wave's FFT-prime-factorization heuristic.
# grid:
#   pml_size: 20             # override if needed
```

Or in MATLAB:

```matlab
parameters.acoustic_backend = 'jwave';

parameters.jwave.solver          = 'helmholtz';
parameters.jwave.python_env      = '/home/user/.conda/envs/jwave_env/bin/python';
parameters.jwave.python_env_gpu  = '/home/user/.conda/envs/jwave_env_gpu/bin/python';
parameters.jwave.dry_run         = false;
% parameters.jwave.use_gpu = true;   % override auto-detection
% parameters.grid.pml_size = 20;     % override the wavelength-based 'auto' PML
```

`python_env` is the path of the CPU environment (its root, its `bin` directory, or the Python executable). `python_env_gpu` points to the GPU (CUDA) environment. Leaving both empty falls back to the system `python3`, which is unlikely to have jwave installed. You get a `sensor_data` struct with the same fields as from k-Wave, so all downstream PRESTUS modules work unchanged.

??? note "Which environment runs"
    | `use_gpu` | GPU detected (`nvidia-smi`) | `python_env_gpu` set | Environment used |
    |---|---|---|---|
    | `[]` (auto) | yes | yes | GPU |
    | `[]` (auto) | yes | no | CPU (warning) |
    | `[]` (auto) | no | – | CPU |
    | `true` | – | yes | GPU |
    | `true` | – | no | **error** |
    | `false` | – | – | CPU |

**Progress and live view.** The time-domain solvers report their steps to the [Nexus](doc_nexus_tabs.md#run-tab) progress bar. The Helmholtz solver reports nothing while it solves, unless [`simulation.live_view: 1`](doc_parameters_simulation.md#simulation-type-simulation) is set: it then solves in chunks of 10 GMRES restart cycles and reports the residual after each chunk. With live view, every solver saves the evolving pressure field (Helmholtz: the current iterate).

## Python Environment Setup

jwave needs a dedicated Python environment. A separate environment avoids version conflicts with other PRESTUS Python tools ([BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain), [SimNIBS](https://simnibs.github.io/simnibs/build/html/index.html)).

**1. Create the environment**

```bash
# conda
conda create -n jwave_env python=3.11
conda activate jwave_env

# venv (no conda)
python3.11 -m venv ~/.venvs/jwave_env
source ~/.venvs/jwave_env/bin/activate
```

**2. Install jwave and its dependencies.** `jax[cpu]` and `jax[cuda12]` install different builds of `jaxlib` and are mutually exclusive: installing one overwrites the other. Pick the variant that matches your hardware.

```bash
# CPU only (macOS, Windows, or CPU-only Linux)
pip install "jax[cpu]" jwave h5py numpy

# NVIDIA GPU (CUDA 12, Linux x86-64 only)
pip install "jax[cuda12]" jwave h5py numpy
```

JAX CUDA wheels are published for Linux x86-64 only; the `jaxlib-cuda` wheel is large (~2 GB). A CUDA install falls back to CPU when no GPU is present at runtime. On macOS or Windows, always use `jax[cpu]`.

**3. Verify the installation**

```bash
python -c "import jwave, jax, h5py; print('jwave:', jwave.__version__); print('jax:', jax.__version__)"
python -c "import jax; print(jax.devices())"
# Expected on GPU host: [CudaDevice(id=0)]
# Expected on CPU host: [CpuDevice(id=0)]
```

**4. Find the Python executable path** and copy it into `jwave.python_env`:

```bash
which python   # inside the activated environment
# e.g. /home/user/.conda/envs/jwave_env/bin/python
#      /home/user/.venvs/jwave_env/bin/python
```

??? note "Parallel CPU and GPU environments"
    If you work on both CPU and GPU machines (a laptop and an HPC cluster with shared home directories), maintain two environments and point `jwave.python_env` and `jwave.python_env_gpu` to them:

    ```bash
    # CPU environment (default, works everywhere)
    conda create -n jwave_env python=3.11
    conda activate jwave_env
    pip install "jax[cpu]" jwave h5py numpy

    # GPU environment (Linux + CUDA 12 hosts only)
    conda create -n jwave_env_gpu python=3.11
    conda activate jwave_env_gpu
    pip install "jax[cuda12]" jwave h5py numpy
    ```

    ```yaml
    jwave:
      python_env: '/home/user/.conda/envs/jwave_env/bin/python'       # CPU
      # python_env: '/home/user/.conda/envs/jwave_env_gpu/bin/python'  # GPU
    ```

## Limitations

- **Acoustic only.** jwave has no bio-heat solver. Thermal simulation (`run_heating_sims`) after a jwave acoustic run uses k-Wave (`thermal_backend: kwave`); BabelBrain thermal requires `acoustic_backend: babelbrain`.
- **No axisymmetric grids.** jwave requires a full 3-D domain.
- **GPU: NVIDIA CUDA only.** JAX does not officially support AMD [ROCm](https://rocm.docs.amd.com/) or Apple [Metal](https://developer.apple.com/metal/). CPU execution works on any platform.

The solvers and every `jwave` setting are described in [jwave Configuration](doc_jwave_configuration.md).
