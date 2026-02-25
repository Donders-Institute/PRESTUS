# Backend code implementation

PRESTUS allows code deployment using different computing setups (`parameters.code_type`). This will primarily affect acoustic simulations. With either of the below GPU variants, thermal simulations will also use GPU acceleration (if kWave 1.4.1 is provided).

- `matlab_cpu`
    - `kspaceFirstOrder3D` | `kspaceFirstOrder2D` | `kspaceFirstOrderAS` 
- `matlab_gpu`
    - MATLAB with GPU acceleration via `DataCast='gpuArray-single'`.
    - Requires Parallel Computing Toolbox
- `cpp_cpu` (C++)
    - `kspaceFirstOrder3DC` warps the C++ binary `kspaceFirstOrder-OMP`. 
    - Requires the specific binary from [k-Wave downloads](http://www.k-wave.org/download.php). 
- `cpp_gpu`
    - `kspaceFirstOrder3DG` wraps the C++/CUDA binary `kspaceFirstOrder-CUDA`, offloading computation to NVIDIA GPUs.
    - Requires a CUDA-capable Nvidia GPU and the specific binary from [k-Wave downloads](http://www.k-wave.org/download.php).

To use C++ with GPU support you may need to recompile the binary for your GPU version (e.g., for the A100 GPUs on the Donders HPC). 

- Download the [LINUX source files for kWave 1.3](http://www.k-wave.org/download.php). This includes the Makefile to compile the binary.
- To compile the binaries, you can follow the instruction by [TU Delft](https://qiweb.tudelft.nl/sysman/kwave_hpc.html). Note that you need to include your GPU type in the Makefile.
- Example scripts for compiling binaries on the Donders HPC can be found in `PRESTUS/examples/hpc_compile...`
- Place the updated binary into `PRESTUS/toolboxes/k-Wave/binaries/`