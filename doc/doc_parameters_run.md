# Run, Paths, and Outputs Parameters

These parameters say which subject runs where, which files it reads and writes, which pipeline stages run, and how cluster jobs are sized.

On this page:

- [Subject and Platform (`subject_id`, `platform`)](#subject-and-platform-subject_id-platform): subject identifier and execution platform.
- [Data Paths (`path`)](#data-paths-path): input, segmentation, simulation, and [Localite](https://www.localite.de/) folders.
- [Environment & Toolbox Paths (`startup`)](#environment-toolbox-paths-startup): [SimNIBS](https://simnibs.github.io/simnibs/build/html/index.html), [MATLAB](https://www.mathworks.com/products/matlab.html) and other external tool paths.
- [I/O Management (`io`)](#io-management-io): output affix, overwriting, cache checks, and which matrices are saved.
- [Pipeline Modules (`modules`)](#pipeline-modules-modules): which pipeline stages run.
- [Report (`report`)](#report-report): the HTML reports a run writes.
- [Study Setup (`study`)](#study-setup-study): study manifest and output folder name.
- [High-Performance Computing (`hpc`)](#high-performance-computing-hpc): GPU, partition, time and memory limits for cluster jobs.

## Subject and Platform (`subject_id`, `platform`)

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `subject_id` | Subject identifier used for file naming and folder management. | `0` | Placeholder; must be set before calling `prestus_pipeline_start(parameters)`. |
| `platform` | Execution platform. | `'auto'` | `auto` / `slurm` / `qsub` / `matlab`. `auto` detects [SLURM](https://slurm.schedmd.com/) or qsub and errors when neither is installed; use `matlab` to run in the current session. See [HPC](doc_hpc.md#where-a-run-executes). |

---
## Data Paths (`path`)

Give file and folder settings as absolute paths. A relative path in these settings, and in the other file and folder settings of a run, is taken from the PRESTUS folder; the `*_pattern` settings are relative to `path.anat`.

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `anat` | Absolute path to structural input data. | none | Mandatory |
| `sim` | Absolute path to simulation outputs. | none | Mandatory |
| `seg` | Absolute path to SimNIBS segmentations (m2m folders). | none | Mandatory |
| `placement` | Placement outputs: one folder per run, `<placement>/sub-XXX/<mode>/<target_name>/` (heuristic, [PlanTUS](https://github.com/mlueckel/PlanTUS), [TUS_entry](https://github.com/CyrilAtkinson/TUS_entry)). | `''` | Empty = `<sim>/placement`. Share one folder across setups to reuse placements. See [Where outputs live](doc_outputs.md#where-outputs-live). |
| `localite_raw` | Raw Localite exports, `<localite_raw>/sub-XXX/ses-XX/[<loc_folder>/]…`. | `''` | See [Localite I/O](doc_placement_neuronav_coordinates.md#localite-file-layout). |
| `localite_post` | Ingested Localite markers, `<localite_post>/sub-XXX/ses-XX/InstrumentMarker_sub-XXX_ses-XX.xml/.json`. | `''` | Written by `neuronav_ingest_markers` / `neuronav_gui`; read by `placement.mode: 'localite'` (automatic mode). |
| `localite_t1` | Folder with the planning T1s loaded into Localite (`[sub-XXX/]<localite_t1_pattern>`). | `''` | Empty = `localite_raw`. |
| `localite_t1_pattern` | File-name template of the planning T1; `%s` = `sub-XXX`. | `'%s_*T1*.nii*'` | Replaces `neuronav_t1_pattern` (still read). |
| `localite_export` | Output folder for files to import into Localite (heuristic placement), per subject `sub-XXX/`. | `''` | Resolved into `io.dir_localite_export`. |
| `transducer_library` | Your own transducer library for every study of this project: equipment entries (`*_transducer.yaml`, `*_driver.yaml`), calibration entries, and `profiles/`. | `''` | A folder or a list of folders. Set it in the project copy of `config_default.yaml`. Checked after `transducer.library_path` and before the `PRESTUS_TRANSDUCER_LIBRARY` environment variable and the shipped library. See [Your own transducers and calibrations](doc_transducer_library.md#your-own-transducers-and-calibrations). |
| `pct` | Pseudo-CT outputs, per subject `<pct>/sub-NNN/`. | `''` | Empty = `<parent of path.seg>/pct/sub-NNN/`, next to the SimNIBS folder (never inside it). Resolved into `io.dir_pct`. |
| `acoustic_input` | Start from a prior pressure map instead of simulating: a peak-pressure NIfTI [Pa] in this subject's T1w space, e.g. an earlier run's `nii/*_T1w*_pressure.nii.gz`. | `''` | Its JSON sidecar is checked against the current settings and used to rescale the map to `transducer.target_isppa_wcm2`; acoustic analysis and thermal simulation then run as usual. Replaces `io.external_acoustic_nifti` (still read). See [Starting from earlier stage outputs](doc_outputs_reuse.md#starting-from-earlier-stage-outputs). |
| `localite` | **Deprecated.** | `''` | Read as `localite_raw` and as `localite_export` when those are empty. |
| `t1_pattern` | T1 image path template relative to `path.anat`. | `'sub-%1$03d_T1w.nii*'` | Supports `%03d`-style subject ID substitution |
| `t2_pattern` | T2 image path template relative to `path.anat`. | `''` | Optional; empty skips the T2 (e.g. `'sub-%1$03d_T2w.nii*'`). |
| `ct_pattern` | Measured CT path template relative to `path.anat`, read by the skull mapping when `pct.enabled` is 1. | `''` | In Hounsfield units and co-registered to the T1w. Empty = the pseudo-CT generated from the UTE image in `t2_pattern`. See [Using a Measured CT](doc_pseudoCT.md#using-a-measured-ct). |
| `bone_mask_pattern` | Bone mask file for `pct.bone_source: mask`, relative to `path.anat` unless absolute. | `''` | Nonzero voxels are bone; co-registered to the T1w, any voxel grid. See [Skull from the CT](doc_pseudoCT.md#skull-from-the-ct). |

---
## Environment & Toolbox Paths (`startup`)

Segmentation and MNI conversion need `startup.simnibs_bin_path`. The other paths are needed for pseudo-CT generation or to add folders to the MATLAB path.

| **Parameter** | **Description** | **Comments** |
|---|---|---|
| `simnibs_bin_path` | Absolute path to the bin folder of a SimNIBS 4 install (holds `charm`). | Mandatory for segmentation and MNI conversion. Selects the SimNIBS version. |
| `matlab_bin_path` | Absolute path to MATLAB binary. | Required for pseudo-CT generation on HPC (e.g. `/opt/matlab/R2024a/bin/matlab`). Local execution falls back to `matlabroot`. |
| `fsl_bin_path` | Absolute path to [FSL](https://fsl.fmrib.ox.ac.uk/fsl/) bin directory. | Required for pseudo-CT generation (e.g. `/opt/fsl/bin`). Prepended to PATH at runtime. |
| `ants_bin_path` | Absolute path to [ANTs](https://github.com/ANTsX/ANTs) bin directory. | Required for pseudo-CT generation (e.g. `/opt/ants/bin`). Prepended to PATH at runtime. |
| `ants_simg_path` | ANTs [Singularity](https://apptainer.org/) image (`.simg`/`.sif`). | When set, ANTs tools run via `singularity exec`. |
| `modules_to_load` | HPC environment modules loaded before pseudo-CT generation. | e.g. `['matlab/R2023b', 'fsl/6.0.7', 'ants/2.5.3']` |
| `paths_to_add` | Paths to add with `addpath()`. | [cell] e.g. `{"path/to/x"}` |
| `subpaths_to_add` | Paths to add recursively with `addpath(genpath())`. | [cell]; relative to the MATLAB working directory |

---
## I/O Management (`io`)

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `output_affix` | Optional affix for output file names. | `''` | Differentiates outputs for the same subject/transducer (e.g. different intensities or targets) |
| `overwrite_files` | File overwrite behaviour. | `'always'` | `never` / `always` / `force` / `ask`. Does NOT apply to SimNIBS segmentations. In chained pipelines (async, multi-I<sub>SPPA</sub>, uncertainty), `always` re-runs only the output stages (CSV, NIfTI, HTML report) and reuses cached intermediates; `force` re-runs every stage. See [Overwrite behaviour](doc_outputs_reuse.md#overwrite-behaviour). |
| `cache_signature_check` | Response when a reused grid, medium, or acoustic cache, or a stage output a run starts from (placement run, pseudo-CT, pressure map), was built from different settings. | `'error'` | `error` / `warn` (reuse anyway) / `off`. Caches: applies when one is reused because of `overwrite_files: never` (or a declined overwrite). See [Cache reuse checks](doc_outputs_reuse.md#cache-reuse-checks) and [Starting from earlier stage outputs](doc_outputs_reuse.md#starting-from-earlier-stage-outputs). |
| `acoustic_cache_remap` | What a run does with a reused acoustic cache written on a different grid (e.g. the head model changed in between). | `0` | `0`: re-simulate, or stop when the run sets up no source; `1`: map the field onto this grid through world space (pressures resampled as p²) and record the mapping in the acoustic provenance. 3-D runs only. See [Cache reuse checks](doc_outputs_reuse.md#cache-reuse-checks). |
| `overwrite_simnibs` | Overwrite SimNIBS segmentation results? | `0` | `1 = yes`, `0 = no`. Re-runs charm (`--forcerun`), which deletes the m2m folder first. Without it, an existing segmentation is reused even when another SimNIBS version or other charm settings are selected (with a warning). |
| `save_matrices` | Save all intermediate matrices (`.mat` / cache files)? | `[]` | `1` = save all, `0` = save none, `[]` = use the per-step defaults. Applies to every per-step matrix flag left empty; a per-step flag set to `0`/`1` always wins. `0` also triggers cleanup of uncertainty intermediates. |
| `save_source_matrices` | Save `kwave_source.mat` (transducer source array)? | `[]` → `1` | `[]` inherits `save_matrices`, else default. Needed by multi-stage pipelines (e.g. uncertainty). |
| `save_acoustic_matrices` | Save acoustic results (`sensor_data`, `kgrid`, medium; can be very large)? | `[]` → `0` | `[]` inherits `save_matrices`, else default. Needed to reload acoustics (e.g. separate thermal runs, calibration debugging). |
| `save_thermal_matrices` | Save `heating_res.mat` (temperature time series and field maps)? | `[]` → `1` | `[]` inherits `save_matrices`, else default. |
| `save_grid_cache` / `save_medium_cache` | Save grid / medium caches? | unset → `0` | Inherit `save_matrices` when unset. |
| `save_property_maps` | Save per-property NIfTIs to nii/properties/ (layered medium only)? | `0` | Not affected by `save_matrices`. |
| `save_transducer_nii` | Save a model of the planned transducer as a T1-space NIfTI (`nii/sub-XXX_<medium>_T1w<affix>_transducer.nii.gz`)? | `1` | Layered medium, annular transducers. Voxel values: `1` = emitting face (smallest `elem_id_mm` to largest `elem_od_mm`), `0.5` = the thickness behind it (`transducer_nii_back_mm`), `0` = elsewhere, including the central opening. Built from the planned T1-voxel positions, on the T1 voxel grid, for placement checks next to the T1 and the acoustic maps. No MNI copy. Not affected by `save_matrices`. See [NIfTI Images](doc_outputs.md#nifti-images). |
| `transducer_nii_back_mm` | Thickness added behind the emitting face in the transducer NIfTI [mm]. | `3` | Display only, value `0.5`: the emitting face (value `1`) stays where PRESTUS places it. `0` gives the face alone. |
| `save_MNI` | Save NIfTI outputs in MNI space in addition to native T1w space (requires SimNIBS m2m folder)? | `1` | Not affected by `save_matrices`. |
| `save_heatingvideo` | Save a video of incremental heating? | `0` | `1 = yes`, `0 = no` |
| `save_Ipa` | Save the pulse-average intensity (Ipa) NIfTI map? | `1` | Suppressed automatically when `timing.pd` is absent. |
| `save_Ita` | Save the temporal-average intensity (Ita) NIfTI map? | `1` | Suppressed automatically when `timing.pd`/`timing.pri` are absent. |
| `save_p_complex` | Save the complex pressure field (`p_complex`) in the acoustic cache? | `0` | *Experimental* for superposition. Every backend. k-Wave only: fitted over the last quarter period, and `p_max_all`/`p_min_all` become ± the fitted amplitude ([Complex Pressure Field](doc_multitransducer.md#complex-pressure-field)). |
| `deface_plots` | Remove facial anatomy from all QC figures? | `1` | `1 = yes`, `0 = no`. Uses the SimNIBS eye label to anchor the cut plane ~30 mm posterior to the eye socket (ear-canal level). Requires the SimNIBS segmentation (`path.seg`). Has no effect on saved NIfTIs; see `segmentation.deface` for that. |
| `adopted_heatmap` | Path to an existing **temperature** heatmap NIfTI (`heating_end.nii.gz`) to use as the thermal starting point for a follow-up simulation. Set automatically by the sequential dispatcher; can be overridden manually. | none | Optional. Used for sequential multi-target runs. |
| `adopted_cem43` | Path to an existing CEM43 heatmap NIfTI to accumulate heating across runs. Set automatically by the sequential dispatcher. | none | Optional. Used for sequential multi-target runs. |
| `adopted_cem43_iso` | Path to an existing ISO-variant CEM43 heatmap NIfTI (`CEM43_iso_end.nii.gz`) to accumulate across runs. Set automatically alongside `adopted_cem43` by the sequential dispatcher. | none | Optional. Used for sequential multi-target runs. |

**Runtime-derived (set by `path_log_setup`):**

| **Field** | **Description** |
|---|---|
| `dir_output` | Resolved per-subject output directory. Set from `path.sim` (+ subject subfolder). |
| `dir_nii` | NIfTI output subfolder (`nii/`). Created automatically. |
| `dir_nii_T1w` | Alias for `dir_nii`; space encoded in filename. |
| `dir_nii_MNI` | Alias for `dir_nii`; space encoded in filename. |
| `dir_img` | Image/figure output subfolder (`img/`). Used for all acoustic, thermal, and preprocessing figures. Created automatically. |
| `dir_tabular` | Directory for per-subject CSV results tables (same as `dir_output`). |
| `dir_reports` | Directory for simulation reports (same as `dir_output`). |
| `dir_logs` | Log subfolder (`log/`). Created automatically. |
| `dir_cache` | Cache subfolder (`cache/`). Created automatically. |
| `dir_debug` | Debug subfolder (`debug/`). Created when `simulation.debug = 1`. |
| `dir_debug_preproc` | Debug subfolder for preprocessing (`debug/preproc/`). |
| `dir_debug_medium` | Debug subfolder for medium setup (`debug/medium/`). |
| `dir_debug_source` | Debug subfolder for source setup (`debug/source/`). |
| `dir_pct` | Resolved subject pCT directory (`get_pct_dir`): `<path.pct>/sub-NNN`, or `<parent of path.seg>/pct/sub-NNN`. |
| `dir_localite_export` | Per-subject Localite export folder (`<path.localite_export>/sub-XXX/`, fallback `path.localite`); `''` if unset. |
| `filename_table` | Path to the per-subject CSV results table. |
| `filename_kwave_source` | Path to cached k-Wave source `.mat` file. Set by `source_sensor_setup`. |

---
## Pipeline Modules (`modules`)

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `run_transducer_placement` | Resolve `trans_pos`/`focus_pos` from `placement.mode`. | `1` | Required when `placement.manual.coordinate_space` is `ras` or `mni` (positions are converted here). |
| `run_grid_setup` | Setup grid and run head processing. | `1` | Mandatory for simulations. |
| `run_medium_setup` | Map medium acoustic properties. | `1` | Mandatory for simulations. |
| `run_source_setup` | Set up acoustic source. | `1` | Mandatory for simulations. |
| `run_water_baseline` | Run a free-water simulation that records the I<sub>SPPA</sub> of the current drive, for scaling to `transducer.target_isppa_wcm2`. | `'auto'` | `auto`: when a target is set and the transducer has no calibration I<sub>SPPA</sub> or the target differs from it. `1`: always. `0`: only when the target differs from a calibration I<sub>SPPA</sub>. A target that cannot be met stops the run (`prestus:isppa:noBaseline`). See [I<sub>SPPA</sub> targeting](doc_simulations-acoustic.md#isppa-targeting). |
| `run_acoustic_sims` | Run acoustic simulations. | `1` | |
| `run_acoustic_analysis` | Run acoustic analysis. | `1` | |
| `run_heating_sims` | Run thermal simulations. | `0` | Enable after acoustic results are validated. |
| `run_thermal_analysis` | Run thermal analysis. | `1` | |
| `run_nifti_creation` | Export results as NIfTI files. | `1` | |
| `run_posthoc_water_sims` | Run free-water reference simulations after head simulations. | `1` | |
| `generate_report` | Generate self-contained HTML simulation report. | `1` | |
| `run_targeting_analysis` | Measure targeting accuracy from the focal-mask NIfTIs: simulated focus vs planned focus (native), vs MNI target / ROI (MNI). Writes `sub-XXX_<medium><affix>_targeting.csv`. | `1` | Layered or phantom media; needs `run_nifti_creation`. See [Targeting accuracy](doc_group.md#targeting-accuracy). |
| `generate_group_report` | Aggregate completed subjects into a group HTML report at the end of the run. | `1` | Best-effort; auto-discovers subjects. Or run `prestus_group_report_start` manually. Switched off for runs started by `prestus_study_start`, which builds one study report at the end. |
| `segmentation_only` | Stop after segmentation; skip grid setup and all simulations. | `0` | Stops after the segmentation stage for every medium; only layered media are segmented, so water and phantom runs compute nothing. |

All flags: `1 = yes`, `0 = no`.

---
## Report (`report`)

Whether the simulation report is written is set by `modules.generate_report`; `style` and `style_switcher` change only how the reports look.

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `style` | Visual skin of the HTML report: `instrument`, `classic` (flat cards), `compact`, `paper`. Changes presentation only, not content. | `'instrument'` | Try skins on an existing output folder with `preview_report_styles(<subject dir>)`. |
| `style_switcher` | Embed all skins plus a dropdown to switch between them in the report. | `0` | |
| `segmentation` | Write the [segmentation report](doc_head-model.md#segmentation-report) `sub-XXX_layered_segmentation_report<affix>.html`. | `1` | Layered media. Shows each step from the SimNIBS segmentation to the simulation medium. |
| `debug` | Write the [debug report](doc_outputs.md#debug-outputs) `sub-XXX_<medium>_debug_report<affix>.html`. | `1` | Only with `simulation.debug: 1`. Also written when the run fails. |

---
## Study Setup (`study`)

Multi-subject studies with per-run settings; see [Group analysis](doc_group_study-setup.md#study-setup).

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `manifest` | Study manifest (`.tsv`/`.csv`): one row per run (`subject_id`, optional `condition`, `include`, `config`) plus any config key as a column. | `''` | Relative paths resolve against the study config's folder. Run with `prestus_study_start(<config>)`; report with `prestus_group_report_start(<config>)`. |
| `name` | Name of the study outputs folder `<path.sim>/group/<name>/`. | `''` | `''` = manifest file name without extension. |

---
## High-Performance Computing (`hpc`)

See [k-Wave](doc_kwave.md) and [HPC](doc_hpc.md).

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `gpu` | Request a specific GPU. | `''` | e.g. `"nvidia_a100-sxm4-40gb:1"`; usually not needed. `scontrol show nodes \| egrep -o gres/gpu:.*=[0-9] \| egrep -o 'nvidia_.*=' \| sort \| uniq \| sed 's/=//'` lists the available types. See [GPU support](doc_hpc.md#gpu-support). |
| `partition` | Request a dedicated queue partition. | `''` | |
| `reservation` | Request a reserved queue. | `''` | |
| `name` | Cluster presets. | `'default'` | `'snellius'` selects the [Snellius](https://www.surf.nl/en/services/compute/snellius-the-national-supercomputer) partition presets (`snellius.*`, below). |
| `exclude` | Nodes to avoid, comma-separated (SLURM `sbatch --exclude`). | `''` | e.g. a node whose GPU is held by a stale process. |
| `cores` | CPU cores per job (SLURM `--cpus-per-task`). | `[]` | `[]` = cluster default (often 1). The k-Wave CPU builds (`cpp_cpu`, `matlab_cpu`) use this many threads. Set by the presets with `hpc.name: 'snellius'`. |
| `wait_for_job` | Block MATLAB until the submitted HPC job completes? | `false` | |
| `timelimit` | Job time limit. | `'04:00:00'` | |
| `memorylimit` | Memory limit [GB]. | `20` | |
| `posthoc_water_timelimit` | Time limit of the post-hoc free-water jobs (`modules.run_posthoc_water_sims`). | `'05:00:00'` | |
| `posthoc_water_memorylimit` | Memory limit of the post-hoc free-water jobs [GB]. | `40` | |
| `ld_library_path` | Prefixed to `LD_LIBRARY_PATH` for SimNIBS tools (charm, MNI conversion). | `''` | `''` = the libraries of the install in `startup.simnibs_bin_path` (`<simnibs package>/external/lib/linux`), so it follows the selected environment. Set only to override, e.g. `/opt/gcc/7.2.0/lib64` for `undefined symbol` errors. The pseudo-CT job uses only a set value. |
| `job_prefix` | Prefix string for HPC job names. | `'PS'` | Keep short; the job name is capped at 20 characters. Overridden to `'TP'` for transducer positioning jobs. |
| `max_wait_checks` | Maximum number of job status checks when `wait_for_job = true`. | `540` | At 1 check/20 s, `540` ≈ 3 hours. |
| `tmp_gb` | Node-local disk per SLURM job [GB] (`sbatch --tmp`). | `0` | `0` = off. With a value, the k-Wave C++ `.h5` files go to `scratch_dir` inside the job ([Node-local scratch](doc_hpc.md#node-local-scratch)). |
| `scratch_dir` | Scratch root used with `tmp_gb > 0`. | `'auto'` | `'auto'` = first existing of `$PRESTUS_SCRATCH`, `$SLURM_TMPDIR`, `/scratch/$USER/$SLURM_JOB_NAME`, `/scratch/$USER/$SLURM_JOB_ID` (else created under `/scratch/$USER`), then `$TMPDIR`. Or a path; `$VAR` is expanded and the folder created. |
| `depend_type` | Scheduler dependency type for `depend_job_id(s)`. | `'afterok'` | `'afterok'` or `'afterany'`. Internal (not in `config_default.yaml`); `prestus_study_start` uses `'afterany'` for the study-report job. |


### Snellius GPU Partitions (`snellius`)

Used when `hpc.name: 'snellius'`: `hpc.partition` selects one of these blocks, which sets the job's resources.

| **Parameter** | **Description** | **Default** | **Comments** |
|---|---|---|---|
| `gpu_a100.cores` | CPU cores requested on `gpu_a100`. | `18` | |
| `gpu_a100.memorylimit` | Memory [GB] on `gpu_a100`. | `120` | |
| `gpu_a100.timelimit` | Wall time on `gpu_a100`. | `'120:00:00'` | Partition maximum. |
| `gpu_a100.n_gpu` | GPUs requested on `gpu_a100`. | `1` | |
| `gpu_h100.cores` | CPU cores requested on `gpu_h100`. | `16` | |
| `gpu_h100.memorylimit` | Memory [GB] on `gpu_h100`. | `180` | |
| `gpu_h100.timelimit` | Wall time on `gpu_h100`. | `'120:00:00'` | Partition maximum. |
| `gpu_h100.n_gpu` | GPUs requested on `gpu_h100`. | `1` | |

---
