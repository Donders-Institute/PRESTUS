# Simulation backends

PRESTUS runs the acoustic and the thermal simulation through interchangeable backends. Pre-processing (segmentation, medium mapping, grid, transducer placement) and post-processing (focal metrics, CEM43, NIfTI export, HTML report) are shared, so switching backend is a config change:

```yaml
acoustic_backend: kwave    # 'kwave' (default) | 'babelbrain' | 'jwave'
thermal_backend:  kwave    # 'kwave' (default) | 'babelbrain'
```

| | k-Wave (default) | [BabelBrain](doc_babelbrain.md) | [jwave](doc_jwave.md) |
|---|---|---|---|
| **Acoustic solver** | Pseudo-spectral time domain | Staggered-grid FDTD (BabelViscoFDTD) | Helmholtz (frequency domain) or pseudo-spectral, in JAX |
| **Thermal solver** | k-Wave `kWaveDiffusion` (Pennes) | BabelBrain Pennes BHTE | — (use k-Wave thermal) |
| **Tissue acoustics** | PRESTUS `medium_properties` | PRESTUS `medium_properties` (injected; fluid bone) — or BabelBrain's own HU mapping via `babelbrain.ct_mapping` | PRESTUS `medium_properties` |
| **Tissue thermal properties** | PRESTUS `medium_properties` | BabelBrain's built-in table (differs in absorption fraction, see below) | — |
| **GPU** | CUDA (C++ binary), MATLAB GPU arrays | CUDA, OpenCL, Metal | CUDA (JAX) |
| **Not supported** | — | Axisymmetric grids | Axisymmetric grids, `kWaveArray` sources |
| **Extra install** | k-Wave submodule (+ binaries for `cpp_*`) | BabelBrain submodule + conda env | Python env with jwave/JAX |

Valid combinations:

| `acoustic_backend` | `thermal_backend` | |
|---|---|---|
| `kwave` / `jwave` | `kwave` | Supported |
| `babelbrain` | `kwave` | Supported |
| `babelbrain` | `babelbrain` | Supported |
| `kwave` / `jwave` / external NIfTI | `babelbrain` | Rejected at config load (BabelBrain thermal needs BabelBrain's `DataForSim.h5`) |

### Are results comparable across backends?

- **Acoustics — designed to be comparable.** All three backends use the same drive phases and the same PRESTUS tissue properties. BabelBrain models bone as a fluid, like k-Wave and jwave (shear waves off), and its output is calibrated against the analytic free-water field (`babelbrain.water_amp_calibration`). What remains is solver discretisation and BabelBrain's bone-attenuation Q-correction. Compare against a k-Wave run before relying on absolute BabelBrain values. Setting `babelbrain.ct_mapping` switches BabelBrain to its own HU→property models, and the backends then no longer share tissue acoustics.
- **Thermal — same acoustic amplitude, different tissue table (known limitation).** BabelBrain thermal receives the in-brain peak ISPPA of the same scaled field that k-Wave thermal consumes, so the heat source amplitude matches. BabelBrain's built-in thermal table matches the PRESTUS defaults for conductivity, specific heat and perfusion, but uses a different **absorption fraction**: skin/brain 0.85 (PRESTUS 1.0), skull 0.15–0.16 (PRESTUS 0.28). Expect lower BabelBrain heating, particularly in skull. PRESTUS thermal `medium_properties` (and uncertainty variants) cannot currently be passed to BabelBrain thermal; use `thermal_backend: kwave` where they must apply. It also does not apply `timing.post_ptri_dur` and does not compute ISO CEM43. Details: [BabelBrain thermal](doc_babelbrain.md#thermal-simulation).

After setting up or updating a BabelBrain or jwave environment, run `smoke_test_backends` (in `tests/`) for a quick check of imports, dependencies and bridge scripts.

---

## k-Wave backend

PRESTUS allows k-Wave deployment using different computing setups (`parameters.simulation.code_type`). This will primarily affect acoustic simulations. With either of the GPU variants, thermal simulations will also use GPU acceleration (if k-Wave 1.4.1 is provided).

- `matlab_cpu`
    - `kspaceFirstOrder3D` | `kspaceFirstOrder2D` | `kspaceFirstOrderAS` 
- `matlab_gpu`
    - MATLAB with GPU acceleration via `DataCast`.
    - Requires Parallel Computing Toolbox
- `cpp_cpu` (C++)
    - `kspaceFirstOrder3DC` warps the C++ binary `kspaceFirstOrder-OMP`. 
    - Requires the specific binary from [k-Wave downloads](http://www.k-wave.org/download.php). 
- `cpp_gpu`
    - `kspaceFirstOrder3DG` wraps the C++/CUDA binary `kspaceFirstOrder-CUDA`, offloading computation to NVIDIA GPUs.
    - Requires a CUDA-capable Nvidia GPU and the specific binary from [k-Wave downloads](http://www.k-wave.org/download.php).

To use C++ with GPU support you may need to recompile the binary for your GPU version (e.g., for the A100 GPUs on the Donders HPC). 

- Download the [LINUX source files for k-Wave 1.3](http://www.k-wave.org/download.php). This includes the Makefile to compile the binary.
- To compile the binaries, you can follow the instruction by [TU Delft](https://qiweb.tudelft.nl/sysman/kwave_hpc.html). Note that you need to include your GPU type in the Makefile.
- Example scripts for compiling binaries on the Donders HPC can be found in `PRESTUS/examples/hpc_compile...`
- For the Donders HPC, you can download a precompiled binary [here](https://github.com/jkosciessa/PRESTUS_bin/raw/refs/heads/main/donders_kwave_bin/kspaceFirstOrder-CUDA).
- Place the updated binary into `PRESTUS/external/k-Wave/binaries/`