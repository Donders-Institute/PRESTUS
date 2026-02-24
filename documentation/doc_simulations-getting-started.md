# How to design your own processing pipeline

### install packages (ensure the paths and subfolders are added in Matlab)
- SimNIBS 4
- k-Wave 1.4.1 (place in PRESTUS' toolbox folder)
    - Note that k-Wave 1.4 is supported, but version 1.4.1 (currently [GitHub exclusive](https://github.com/ucl-bug/k-wave/releases/tag/v1.4.1)) introduced GPU support for thermal simulations with kWaveDiffusion.

### start with the tutorial (documentation/PRESTUS_intro_tutorial.md)
- The acoustic profile is derived from the manufacturers measurements and adjusted to more closely match the actual characteristics of the transducer.
- The structural data is segemented with SimNIBS 4's charm.
- The segmented data is translated to a matrix of acoustic properties.
- k-wave performs acoustic and thermal simulations.
- Select outputs are plotted.

### create a study-specific config file
- Configuration metadata must be fully specified in (a set of) config files (in the .yaml format) placed in the 'configs' folder (see ```doc_config.md```). 
- The ```load_parameters.m``` function imports the configuration. 
- Some of the parameters found in the 'default_config' are mandatory for the pipeline to run while others can be left out based on the requirements of the analysis.
- The different parameters that can be used are found in the 'default_config.yaml'. The default config should generally not be altered. Instead, parameter changes should be specified via an application-specific config file or before performing the subject call in MATLAB).

### specify a transducer and target location
- These can be specified in your config `parameters.transducer.trans_pos` and `parameters.transducer.focus_pos`.
- See also `doc_neuronav.md`

### choose a simulation medium

- see `doc_medium.md`

### choose a simulation type

PRESTUS allows code deployment using different computing setups (`parameters.code_type`).

- `matlab_cpu`
    - `kspaceFirstOrder3D` | `kspaceFirstOrder2D` | `kspaceFirstOrderAS` 
- `matlab_gpu`
    - MATLAB with GPU acceleration via DataCast='gpuArray-single'.
    - Requires Parallel Computing Toolbox
- `cpp_cpu` (C++)
    - `kspaceFirstOrder3DC` warps the C++ binary `kspaceFirstOrder-OMP`. 
    - Requires the specific binary from [k-Wave downloads](http://www.k-wave.org/download.php). 
- `cpp_gpu`
    - `kspaceFirstOrder3DG` wraps the C++/CUDA binary `kspaceFirstOrder-CUDA`, offloading computation to NVIDIA GPUs.
    - Requires a CUDA-capable Nvidia GPU and the specific binary from [k-Wave downloads](http://www.k-wave.org/download.php). 


### run the single_subject_pipeline
- Ideally, the full pipeline only needs a subject id (a number consisting of no more than 3 digits).
- To run multiple subjects or setups in parallel, the single_subject_pipeline can be submitted using high performance computing jobs (see doc_hpc.md).
