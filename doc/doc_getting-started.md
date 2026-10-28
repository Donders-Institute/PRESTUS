# Your first simulation

This page takes you from a fresh installation to a finished simulation and its outputs. It starts in free water (no MRI, no SimNIBS, a few minutes on a CPU) and then points to what changes for a head simulation.

Before you start, install PRESTUS and its MATLAB dependencies ([Installation](doc_installation.md)). SimNIBS is only needed from step 3 on. If terms such as *exit plane*, *ISPPA* or *free-water baseline* are new to you, keep the [glossary](doc_concepts.md#glossary) open alongside.

---

## 1. A free-water simulation

Free water needs no images: the transducer sits in a homogeneous water grid. The settings below are the ones the pipeline test `tests/test_integration_water.m` runs (a single-element 250 kHz bowl on a small 1 mm grid).

```matlab
% from the PRESTUS root folder
addpath(genpath('functions')); addpath(genpath('config')); addpath(genpath('external'));

out = fullfile(tempdir, 'prestus_first');  mkdir(out);

parameters = load_parameters();               % toolbox defaults only
parameters.subject_id = 1;
parameters.platform   = 'matlab';             % run in this MATLAB session
parameters.path.sim   = out;                  % outputs
parameters.path.anat  = out;                  % not used in water, but must be set
parameters.path.seg   = out;

parameters.simulation.medium    = 'water';
parameters.simulation.code_type = 'matlab_cpu';   % 'matlab_gpu' with a GPU

parameters.grid.resolution_mm = 1.0;          % coarse grid for speed
parameters.grid.default_dims  = [72, 72, 128];
parameters.modules.run_posthoc_water_sims = 0;

% single-element annular bowl, driven at elem_p_amp (no ISPPA scaling)
parameters.transducer.type = 'annular';
parameters.transducer.freq_hz = 250e3;
parameters.transducer.annular.elem_n          = 1;
parameters.transducer.annular.elem_id_mm      = 0;
parameters.transducer.annular.elem_od_mm      = 32.0;
parameters.transducer.annular.curv_radius_mm  = 63.2;
parameters.transducer.annular.dist_geom_ep_mm = 52.38;
parameters.transducer.annular.elem_p_amp      = 91590;
parameters.transducer.annular.elem_phase_deg  = 0;
parameters.transducer.target_isppa_wcm2       = [];

% transducer (bowl apex) and focus in grid voxels, beam along z
parameters.transducer.trans_pos = [36, 36, 5];
parameters.transducer.focus_pos = [36, 36, 64];

prestus_pipeline_start(parameters);
```

When it finishes, `out/sub-001/` holds the results:

- `sub-001_water.csv`: the metrics table (ISPPA, focal distance, focal volume, …);
- `sub-001_water_report.html`: the HTML report with the safety dashboard and figures;
- `nii/`, `img/`, `log/` and `cache/`: maps, figures, the run log and intermediate files.

[Reading your results](doc_results.md) explains what to check first; [Outputs](doc_outputs.md) lists every file and column.

## 2. Use a transducer from the library

Real studies rarely type in element geometry and drive. PRESTUS ships an equipment and calibration library: name the transducer and its driving system, set the focal depth and the free-water intensity you want, and PRESTUS fills in geometry, phases and drive. The library is resolved when the config is loaded, so put these keys in a YAML file rather than setting them in MATLAB afterwards:

```yaml
# config_first.yaml
transducer:
  serial: CTX_500_026            # config/equipment/CTX_500_026.yaml
  combo:
    ds_serial: SC_203_035        # config/calibration/CTX_500_026_SC_203_035.yaml
  focal_distance_ep: 50          # focal depth setting, mm from the exit plane
  target_isppa_wcm2: 30          # report results at 30 W/cm² free-water ISPPA
```

```matlab
parameters = load_parameters('config_first.yaml', pwd);
% ... then the same subject, path, simulation and grid settings as in step 1
%     (without the transducer geometry lines), and prestus_pipeline_start(parameters)
```

With a target ISPPA set, PRESTUS also runs a short free-water baseline and scales the field to that intensity ([ISPPA targeting](doc_simulations-acoustic.md#isppa-targeting)). The available devices, how the library works and how to calibrate your own are described in [Transducer Library](doc_transducer_library.md).

## 3. A head simulation

A head (`layered`) simulation adds an MRI, a SimNIBS segmentation and a target. Compared to step 2 you need:

| What | Where it is set | See |
|---|---|---|
| T1 (optionally T2/UTE) per subject | `path.anat`, `path.t1_pattern`, `path.t2_pattern` | [Configuring a study](doc_configuration.md) |
| SimNIBS installation and segmentation folder | `startup.simnibs_bin_path`, `path.seg` | [Installation](doc_installation.md#simnibs-installation) |
| Medium | `simulation.medium: layered` | [Medium Setup](doc_medium.md) |
| Transducer position and target | `transducer.trans_pos`/`focus_pos`, or a `placement.mode` | [Placement](doc_placement.md) |
| Sonication protocol (for heating) | `timing.*`, `modules.run_heating_sims: 1` | [Thermal](doc_simulations-thermal.md) |

On the first run PRESTUS runs SimNIBS `charm` (this takes an hour or more); later runs reuse the segmentation unless `io.overwrite_simnibs` is set. To segment first and simulate later, set `modules.segmentation_only: 1`. Head simulations are much larger than step 1: use a GPU (`simulation.code_type: matlab_gpu` or `cpp_gpu`, see [Backend](doc_backend.md)) or a cluster ([HPC](doc_hpc.md)). Thermal simulations are off by default (`modules.run_heating_sims: 0`); switch them on once the acoustic result looks right.

## 4. Next steps

- Keep your settings in a project folder and a study config: [Configuring a study](doc_configuration.md).
- Prefer forms to YAML: [PRESTUS Nexus](doc_nexus.md).
- Use subject-specific skull properties from a UTE scan: [pseudo-CT](doc_pseudoCT.md).
- Re-run only part of the pipeline: [Modules](doc_modules.md).
- Many subjects or conditions: [Group analysis](doc_group.md#study-setup).
- Sweeps over intensities, several transducers, sequences or tissue-property uncertainty: [Advanced modes](doc_advanced.md).

## Demos

- A simplified 2D benchmarking example: [PRESTUS_2D_demo](https://github.com/jkosciessa/PRESTUS_2D_demo) (runs on a local CPU).
- A 3D demo on SimNIBS' Ernie head as a [DataLad dataset](https://gin.g-node.org/PRESTUS/sim_ernie/). It was made with an earlier PRESTUS version, but shows one way to organise data.
- `examples/` in the toolbox holds runnable scripts (calibration, Localite, pseudo-CT, uncertainty); see `examples/README.md`.
