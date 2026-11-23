# k-Wave

[k-Wave](https://github.com/ucl-bug/k-wave) is the default solver for the acoustic and the thermal simulation and is installed with PRESTUS. Its implementation (`simulation.code_type`) decides whether a run uses the CPU or the GPU, through [MATLAB](https://www.mathworks.com/products/matlab.html) or a C++ binary. The elastic skull model changes the physics instead: it adds shear waves in the skull, which fluid skull models omit and which lowered the intracranial pressure in ex vivo measurements; in PRESTUS it is a draft and not yet validated.

On this page:

- [Solver Implementations (`code_type`)](#solver-implementations-code_type): CPU or GPU, MATLAB or C++, and building the binaries.

!!! note "At a glance"
    - **When you need it:** To choose `simulation.code_type` for your hardware, or to model shear waves in the skull.
    - **Before you start:** Nothing; k-Wave is installed with PRESTUS. The C++ variants need the k-Wave binaries.

## Solver Implementations (`code_type`)

`simulation.code_type` selects the k-Wave implementation. This mainly affects the acoustic simulation; with either GPU variant, thermal simulations also use GPU acceleration (if k-Wave 1.4.1 is provided).

| `code_type` | Implementation | Requirements |
|---|---|---|
| `matlab_cpu` | `kspaceFirstOrder3D` / `kspaceFirstOrder2D` / `kspaceFirstOrderAS` in MATLAB | none |
| `matlab_gpu` | MATLAB with GPU acceleration via `DataCast` | [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html) |
| `cpp_cpu` | `kspaceFirstOrder3DC` wraps the C++ binary `kspaceFirstOrder-OMP` | The binary from [k-Wave downloads](http://www.k-wave.org/download.php) |
| `cpp_gpu` | `kspaceFirstOrder3DG` wraps the C++/[CUDA](https://developer.nvidia.com/cuda-toolkit) binary `kspaceFirstOrder-CUDA`, offloading computation to NVIDIA GPUs | A CUDA-capable Nvidia GPU and the binary from [k-Wave downloads](http://www.k-wave.org/download.php) |
| `cpp_mlx` | `kspaceFirstOrder3DG` wraps the C++/[MLX](https://github.com/ml-explore/mlx) binary `kspaceFirstOrder-MLX`, offloading computation to the Apple Silicon GPU | An Apple Silicon Mac (macOS 14 or later) and the binary with `mlx.metallib`, built as described below |

```yaml
simulation:
  code_type: cpp_gpu   # matlab_cpu | matlab_gpu (default) | cpp_cpu | cpp_gpu | cpp_mlx
```

**Recommended: the C++ codes.** Where the binaries are installed, use `cpp_gpu` on an NVIDIA GPU and `cpp_cpu` otherwise. They run several times faster than their MATLAB counterparts on the same hardware ([Run Time](doc_run-time.md)), keep peak memory low by streaming the source from disk ([Memory](doc_simulations-acoustic.md#memory-use)), and report their progress to the Nexus. In return they need a separate install and temporary disk space for their input and output files ([HPC](doc_hpc.md)). The configuration default is `matlab_gpu`, because it runs without the binaries.

??? note "Precompiled binaries for the Donders HPC"
    The [PRESTUS_bin repository](https://github.com/jkosciessa/PRESTUS_bin/tree/main/donders_kwave_bin) holds three Linux x86-64 binaries. Place the ones you need into `PRESTUS/external/k-wave/k-Wave/binaries/`.

    | Binary | What it is | Used by |
    |---|---|---|
    | [`kspaceFirstOrder-OMP`](https://github.com/jkosciessa/PRESTUS_bin/raw/refs/heads/main/donders_kwave_bin/kspaceFirstOrder-OMP) | k-Wave C++ 1.3 CPU solver ([OpenMP](https://www.openmp.org/), multi-threaded) from the k-Wave website | `cpp_cpu` |
    | [`kspaceFirstOrder-CUDA`](https://github.com/jkosciessa/PRESTUS_bin/raw/refs/heads/main/donders_kwave_bin/kspaceFirstOrder-CUDA) | k-Wave C++ 1.3 CUDA solver, compiled on the Donders HPC with CUDA 11.8 as described below. It loads cuFFT (CUDA 11) and HDF5 1.10 from the system, so it runs only where these libraries are installed. | `cpp_gpu` |
    | [`acousticFieldPropagator-OMP`](https://github.com/jkosciessa/PRESTUS_bin/raw/refs/heads/main/donders_kwave_bin/acousticFieldPropagator-OMP) | C++ backend of k-Wave's `acousticFieldPropagatorC`, which computes the steady-state field of a continuous-wave source in a homogeneous medium | not used by PRESTUS |

??? note "Compiling the CUDA binary for your GPU"
    To use C++ with GPU support you may need to recompile the binary for your GPU version (for example the A100 GPUs on the Donders HPC).

    - Download the [Linux C++/CUDA source files](http://www.k-wave.org/download.php) (version 1.3; the C++ code is versioned separately from the k-Wave 1.4.1 MATLAB toolbox). This includes the Makefile to compile the binary.
    - Follow the compilation instructions by [TU Delft](https://qiweb.tudelft.nl/sysman/kwave_hpc.html). Include your GPU type in the Makefile.
    - On the Donders HPC, `examples/hpc_compile_CUDA_a.sh` downloads and unpacks the source into `~/kwave-1.3-Donders/source`, and `examples/hpc_compile_CUDA_b.sh` loads `cuda/11.8` (the version sets the newest GPU architecture it can target) and compiles `kspaceFirstOrder-CUDA` into `~/kwave-1.3-Donders/binaries`.
    - For the Donders HPC, a precompiled binary is available in the PRESTUS_bin repository (see above).
    - Place the updated binary into `PRESTUS/external/k-wave/k-Wave/binaries/`

??? note "Building the MLX binary (macOS)"
    `cpp_mlx` is experimental and has not been validated on a PRESTUS head simulation. It runs the acoustic simulation on the Apple GPU; thermal simulations run on the CPU. Axisymmetric grids fall back to `matlab_cpu`, and `kwave_elastic` runs with `matlab_cpu`, because the port has no solver for either.

    The binary is the MLX port of the k-Wave CUDA solver by S. Pichardo ([kspacefirstorder-unified, pull request 31](https://github.com/waltsims/kspacefirstorder-unified/pull/31)). Building it needs [CMake](https://cmake.org/) and [Xcode](https://developer.apple.com/xcode/) with the [Metal](https://developer.apple.com/metal/) Toolchain (`xcodebuild -downloadComponent MetalToolchain`).

    ```bash
    git clone --branch MLX-k-Wave-binary https://github.com/ProteusMRIgHIFU/kspacefirstorder-unified.git
    cd kspacefirstorder-unified
    cmake -S repos/kspaceFirstOrder-mlx-darwin -B build/mlx-darwin \
      -DCMAKE_BUILD_TYPE=Release -DCMAKE_OSX_DEPLOYMENT_TARGET=14.0 -DKWAVE_BUNDLE_DEPS=ON
    cmake --build build/mlx-darwin --parallel
    ```

    Copy `kspaceFirstOrder-MLX` and `mlx.metallib` from `build/mlx-darwin` into `PRESTUS/external/k-wave/k-Wave/binaries/`. The two files must stay side by side, because the binary loads `mlx.metallib` from its own folder.

The elastic skull model (shear waves) is described in [Elastic Skull (Shear Waves)](doc_kwave_elastic.md).
