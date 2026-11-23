# Changelog

Notable changes to this project are documented here. 

---

## v0.7.0 (unreleased)

PRESTUS v0.7 adds the PRESTUS Nexus as the interface to configure, check, and run simulations, a shipped library of calibrated transducers with ISPPA targeting by default, and new placement modes (MNI, TUS_entry, PlanTUS 2.0, batch Localite ingestion). BabelBrain and jwave join k-Wave as acoustic backends. Reports gain a ready-to-use Methods section and an ITRUSST reporting overview, and uncertainty runs bound acoustic and thermal results with tissue-property variants aligned with the ITRUSST Acoustic Properties Group. Many fixes change results of v0.6.x runs: check the [rerun guide](#rerun-guide).

#### Upgrade notes

**Changed defaults**

| Key | v0.6.1 → v0.7 | What to do |
|---|---|---|
| `transducer.target_isppa_wcm2` | `[]` → `30` | Results are scaled to 30 W/cm² free-water ISPPA; set `[]` to report the drive level. A target that cannot be met stops the run |
| `modules.run_water_baseline` | `0` → `auto` | A free-water simulation is added when scaling needs one |
| `grid.mode` | `transducer_axis` → `focal_axis` | The grid follows the mean beam axis; multi-transducer grids change, so rerun them rather than reuse caches |
| `grid.default_dims` | `[144, 144, 400]` → `auto` | The free-water grid is sized to the transducer |
| `headmodel.csf_expansion` | `40` → `20` | Now the full dilation distance; halve values set in project and study configs |
| `io.save_source_matrices`, `save_acoustic_matrices`, `save_thermal_matrices` | → `[]` | Inherit `io.save_matrices`; default saving is unchanged |
| `thermal.pulse_model` (new) | v0.6.1 heated pulse by pulse → `train_averaged` | ⚠️ Heating results change: each pulse train heats as one segment, as k-Plan (peak rises within 0.8 %, thermal dose about 1 %). Set `pulse_resolved` to reproduce v0.6.1 heating |

**Changed behaviour**

- Coherent multi-transducer runs scale each transducer to its own `transducer(k).target_isppa_wcm2`; ISPPA sweeps with coherent runs stop with an error
- `io.overwrite_files: always` re-runs the output stages of chained pipelines (async, multi-ISPPA, uncertainty); new `force` also re-runs their intermediate stages
- `placement.reuse: auto` reuses a recorded heuristic, PlanTUS, or TUS_entry placement with matching settings, also under `io.overwrite_files: always`
- Caches and stage outputs (pseudo-CT, placement runs, T1w pressure maps) record the settings and segmentation they were built from; reuse with other settings stops the run (`io.cache_signature_check`)

**Configs and scripts**

- The toolbox `config_default.yaml` is always the base layer and a project copy is merged on top; refresh old project copies with `prestus_config_upgrade`
- Configs from v0.5.0 on load as they are: renamed keys are converted with one warning (`prestus:renamedConfigKey`), and keys without a counterpart stop the run and name the replacement (`prestus:legacyConfigKey`). Unknown study-config keys warn
- Renamed: `transducer.annular.elem_amp` and `transducer.matrix.elem_amp` → `elem_p_amp`; `placement.plantus.additional_offset_mm` → `placement.skin_gap_mm`; `headmodel.skull_wrap_radius` (grid voxels) → `headmodel.skull_wrap_radius_mm` (default 5 mm; an old value is converted at the run's resolution); CSV column `real_focal_distance_mm` → `real_focal_distance_bowl_mm` (next to `real_focal_distance_ep_mm`); update your own scripts that read it
- Removed: `thermal.sensor_xy_halfsize` and `thermal.record_t_at_every_step` (no effect; configs that set them load with a notice). A top-level `expected_focal_distance_mm` stops the run unless every transducer sets `transducer.focal_distance_ep` or `focal_distance_bowl` (then it is ignored with a warning, as in v0.6)
- Deprecated: `path.localite` (use `path.localite_raw` / `path.localite_export`), `placement.localite.markertype` (detected from the files), `pct.skull_mapping: carpino` (use `custom` with `custom_slope: -2194`, `custom_intercept: 2236`)
- Uncertainty: the liberal and conservative variants are replaced by Temp-H/L and Trans-H/L (affixes `_desc-temphigh`, `_desc-templow`, `_desc-transhigh`, `_desc-translow`); `options.liberal_config`, `options.conservative_config`, and their affixes stop with `uncertainty_pipeline:legacyVariant`
- Removed functions: `prestus_gui` and `gui_run_worker` (use `prestus_nexus`); `launch_web_alignment` and the browser alignment viewer (use the Nexus placement viewer); `create_group_MNI_plots`, `combine_plots_by_suffix`, and `changem_vectorized` (use `create_group_maps`); `plot_coronal_slices` and `plot_median_montage`; `calibrate_transducer`, `calibration_extract_correction`, `compute_analytical_solution`, `extract_analytical_profile`, `extract_simulated_profile`, `fit_velocity_to_intensity`, `generate_tran_ini_from_geometry`, `phase_optimization_annulus`, `plot_opt_sim_results`, `recompute_analytical_solution`, and `scale_real_intensity_profile` (see [doc_calibration.md](doc_calibration.md) for the current calibration workflow). `prestus_group_start` runs `prestus_group_report_start` and ignores its plot options
- SimNIBS 4 (charm) is required: SimNIBS 3 (`headreco`) segmentations and the affine transforms `12dof` and `6dof` of `mni2subject_coords_LDfix` / `subject2mni_coords_LDfix` are no longer supported (`prestus:transform:simnibs4Only`); segment such subjects again with charm. `convert_final_to_MNI_matlab` (headreco only) is removed; `convert_final_to_MNI_simnibs` writes MNI outputs
- Scripts from earlier versions keep running with warnings: the call form `prestus_pipeline_start(subject_id, parameters)`, the v0.5.0 functions `single_subject_pipeline`, `single_subject_pipeline_with_slurm`, and `single_subject_pipeline_with_qsub`, and old key names set on the parameters struct after `load_parameters` are converted. Keys of the development builds between v0.5.0 and v0.6.0 (`hpc_timelimit`, `hpc_memorylimit`, `hpc_wait_for_job`, `run_grid_setup` and the other stage flags, `use_pseudoCT`, `smooth_fwhm_mm`) are converted too
- Equipment: entries are grouped by device model (`config/equipment/model-<MODEL>/`, driving systems in `config/equipment/ds/`); flat folders are still read. `IS_PCD15473_01001_OPM` and `_01003_OPM` are removed: use `IS_PCD15473_01001` / `_01003` with `IGT_32_ch_comb_10_ch`

**Output locations**

- Heuristic, PlanTUS, and TUS_entry runs write to `<path.placement>/sub-XXX/<mode>/<target_name>/` (default `<path.sim>/placement`). Earlier PlanTUS results in `<path.seg>/PlanTUS/` and heuristic `tpars_*.csv` files are not found there: re-run them or copy them over
- Pseudo-CTs are written to `<path.pct>/sub-NNN`, or `<parent of path.seg>/pct/sub-NNN` when `path.pct` is empty, never into the SimNIBS folder; pseudo-CTs in the old location are still used, with a warning
- The subject's T1 in MNI space is written to `<parent of path.seg>/mni/sub-NNN/`

**Requirements**

- SimNIBS 4 or later for segmentation
- Fitting transducer library entries needs the Optimization Toolbox; loading them does not
- BabelBrain backend: BabelBrain 0.8.8 (`external/BabelBrain`); update its environment with `pip install babelviscofdtd==1.2.7 psutil`
- Submodules: export_fig v3.58 and PlanTUS v2.0.0; `prestus_setup` fetches missing submodules

**Experimental features**

Not validated in PRESTUS; the docs and `config_default.yaml` mark each one. Check results against a k-Wave fluid run or measurements before you report them.

- Backends: `kwave_elastic`, `jwave`, `babelbrain` (acoustic and thermal), and the `cpp_mlx` code type
- Transducer: dephased library conditions, the within-ring velocity profile, the matrix steering correction (on by default when the library holds one), and the `define_here` options
- Placement: PlanTUS, TUS_entry, the heuristic air criterion, and time reversal (including `placement_search` and `pose_workflow`)
- Head model: pseudo-CT mappings and the air-cavity check
- Runs and outputs: async multi-transducer runs, uncertainty variants (including `uncertainty_fov: partial`), superposition of saved complex fields (`io.save_p_complex`), and the group report and group MNI maps

#### Rerun guide

Results from v0.6.x change if your runs used any of the following. Fixes to features new in v0.7 are not listed.

| If your runs used | What changed | Rerun |
|---|---|---|
| Any head simulation | The simulation ended about 1.9× too early, before steady state (end time derived from the bone sound speed); MI used the peak positive instead of the peak rarefaction pressure; T1w maps were resampled linearly (focal peaks up to ~3 % low, now cubic) | Acoustic stage and NIfTI export |
| Layered head models | The rubber-wrap skull fill ignored brain and skin, and its radius was in grid voxels (now mm, so fills at resolutions other than 0.5 mm change); filled skull holes were labelled trabecular; the CSF-based crop dilated by half of `csf_expansion`; `smooth_img` pads by edge replication; caches were reused after re-segmentation; a reused grid cache (`io.preproc_affix`) replaced the run's own transducer settings | Grid stage onward |
| A T1 not stored RAS+ | An oblique T1 lost its rotation; configured voxel positions ignored the reorientation; the Localite planning image and heuristic Localite coordinates flipped left and right | Whole run; regenerate planning images and heuristic exports |
| Positions from RAS, MNI, or Localite coordinates | RAS → voxel and native → MNI conversions were one voxel off; MNI conversion near the edge of the SimNIBS field pulled points towards the origin; `analysis.mni_warp_method: native` was one voxel off; heuristic exports rounded down and Localite mm exports were one voxel high; voxel-derived focal distances rounded the voxel size to whole mm; `trans_pos_mm` / `focus_pos_mm` mapped 0 mm to voxel 0 | Placement onward (positions move by about one voxel); subtract one voxel from configs that compensated for the old `trans_pos_mm` mapping |
| Localite placement | The instrument matrix was transposed; with `tracker_to_bowl_mm` / `reference_distance_mm` (and always with PlanTUS) the focus was too shallow by `curv_radius_mm − dist_geom_ep_mm` | Placement onward |
| Heuristic or PlanTUS placement | All modes use one exit plane and one scalp gap (`placement.skin_gap_mm`, 5 mm): heuristic transducers sit about 2 mm (CTX-500) and PlanTUS transducers 5 mm further out. PlanTUS skull thickness is now a criterion, and PlanTUS 2.0 smooths the skin | Placement onward; `skin_gap_mm: 0` and `placement.plantus.skin_smoothing_iterations: 0` keep earlier PlanTUS positions |
| Axisymmetric grids | Rings of multi-ring transducers were misplaced; the kWaveArray source raised free-water focal intensity by 12–23 %; NIfTIs of off-centre transducers were cut at the grid centre | Acoustic stage; NIfTI export for off-centre runs |
| Matrix transducers | `define_here` grids used m instead of mm, swapped width and height, and had the wrong cap curvature; `use_kWaveArray: 0` did not switch to the kWaveArray source; `project_on_new_ROC` mixed mm and m; curved arrays without `dist_geom_ep_mm` put the exit plane at the apex; Clover leaves are `clover.az_sep_deg` apart (1- and 2-leaf layouts change) | Acoustic stage onward |
| 2D annular sources | The per-ring outer diameter was used instead of the full aperture | Acoustic stage |
| Several transducers | Coherent: transducers 2..N used transducer 1's drive, phases, geometry, and frequency. Async: heating applied all intensities during every pulse (up to N× too high) and now follows `timing.async_schedule` | Coherent: acoustic stage; async: thermal stage |
| Calibrations or library entries made with v0.6.x | The Rayleigh model integrated over the projected ring area (on-axis focal pressure 1–9 % low); the O'Neil model was evaluated 0.5 mm short; row-vector profiles could give degenerate fits; interpolation to unmeasured depths fixed the focus; calibration files recorded the default ISPPA | Refit; set `target_isppa_wcm2` in configs built from old calibration files |
| Skull mappings `aubry` (density, sound speed, attenuation) or `mueller` (attenuation) | The densest-bone HU was the largest skull HU inside the cropped simulation grid; it is now the 99.9th percentile of the skull over the whole head (`pct.hu_reference_percentile`), and HU above it map as the densest bone | Medium stage onward |
| Pseudo-CT | `treeby` was 27 HU low in the skull; the `kosciessa` fit used the wrong width crossings; the UTE soft-tissue peak was the bin edge; with `mapping_density: k-plan`, the thermal skull density landed on other tissue when layers were removed or reordered | Regenerate the pseudo-CT; thermal stage for the `k-plan` case |
| `pct.mapping_density: k-plan` with heating | The fixed thermal skull density (1850 kg/m³) also entered the heat source, so skull heating was too low where the mapped bone is less dense (0.74 of k-Plan on k-Plan's own pressure field); it now enters only the heat capacity and diffusion, and the heat source uses the mapped density and sound speed, as k-Plan | Thermal stage |
| `headmodel.smooth_properties: true` with heating | The heat source multiplied the separately smoothed attenuation and absorption-fraction maps, which overshoots at tissue boundaries (up to 8× in the last skin voxel before bone at 1 mm FWHM); it now uses the absorption map smoothed as one, so skin and skull heating fall | Thermal stage |
| Uncertainty runs | New variants (see Upgrade notes); variants shared one medium cache with `io.save_medium_cache: 1`; safety cards followed the conservative variant; with a pseudo-CT, the bound variants use a homogeneous skull | Rerun and regenerate reports |
| Reused or repeated outputs | A dotted `io.output_affix` truncated NIfTI names; `overwrite_files: never` / `ask` mixed the CSV and NIfTIs of two runs; `always` was ignored by chained pipelines; cached heating files lost `temp_0`; `io.save_p_complex` stored a column vector, and with k-Wave the field is now a fit over the last quarter period (`p_max_all` and ISPPA are ± the fitted amplitude); phantom NIfTIs had no affine | Rerun affected folders; re-save complex fields; re-export phantom NIfTIs |

#### New

**Interface**

- **PRESTUS Nexus** (`prestus_nexus`; see [doc_nexus.md](doc_nexus.md)): a standalone app to build, check, and run configs locally or on a cluster, also over SSH. It covers a QuickStart from a T1 and an MNI target, transducers and calibration, placement with a 3D viewer, medium and pseudo-CT settings, sonication sequences, run-time and memory prediction, study manifests, results with an ITRUSST reporting overview, and follow-up runs. `calibration_gui` and `neuronav_gui` remain
- Updating: *Help ▸ Check for updates…* in the Nexus and `prestus-nexus --update` (no window needed) move a PRESTUS clone to the latest release; see [doc_nexus.md](doc_nexus.md#updating-prestus)
- Installation: `prestus_setup` puts PRESTUS on the path and fetches submodules; `prestus_check_install` reports the installation; `nexus/install.sh` and `nexus/install.bat` set up the Nexus environment
- Runs check the external tools their stages need before starting (`prestus:missingTools`); `predict_outputs` lists what a run will write
- Run progress (stage, solver progress, time left) in the Nexus for local runs and cluster jobs; `simulation.live_view` plots the evolving pressure and temperature fields without the MATLAB desktop and saves them for the Nexus

**Transducers and calibration**

- Transducer library (`config/calibration/`; see [doc_transducer_library.md](doc_transducer_library.md)): 12 device models resolved from `transducer.serial`, with per-depth refinements fitted to hydrophone profiles. Your own library via `path.transducer_library` or `transducer.library_path`; `refine_transducer_library` fits your own profiles; drive limits warn
- Calibration: focal-depth convention (`calibration.focal_reference`), within-ring velocity profile, matrix steering correction (`transducer.matrix.use_steering_correction`), and dephased conditions (*experimental*)
- Transducers: flat annular elements; matrix element rotation, diameter, and phase override; Clover layouts with any number of leaves

**Placement**

- New modes `mni` (entry and focus in MNI mm) and `tus_entry` (TUS_entry toolbox); `placement.manual.coordinate_space` (`voxel`, `ras`, `mni`)
- PlanTUS 2.0 support, a distance window, and placement at a chosen vertex; heuristic placement of a chosen candidate
- A placement report for every mode; placement runs are recorded and reused (`path.placement`, `placement.reuse`)
- Localite: batch ingestion, target resolution and label maps, multi-coil and bilateral setups (`neuronav.*`, `path.localite_*`)
- Time reversal (`time_reversal.*`, *draft*): refocusing phases, aberration metrics, and a ranking of candidate poses

**Head model and pseudo-CT**

- Defacing (`segmentation.deface`, `io.deface_plots`); charm settings file and `segmentation.noneck`; one head model per subject across placements (`headmodel.reuse_across_placements`)
- Skull mappings `air-tissue`, `webb`, `pichardo`, `mcdannold`, `marsac-babelbrain`, and `custom`; HU → density tables (`pct.density_calibration`); a measured CT instead of a pseudo-CT (`path.ct_pattern`); `pct_predict_skull_properties`
- Skull from the CT: an HU threshold on the simulation grid instead of the segmentation's skull, as BabelBrain's CT path (`pct.bone_source: hu`), or a bone mask file (`pct.bone_source: mask`, `path.bone_mask_pattern`), both reconciled with the segmentation the same way
- CT interpolation (`pct.interpolation`, including BabelBrain's cubic spline), one regularisation scheme for HU and the mapped skull properties (`pct.limits`), and a whole-head HU reference for mappings that normalise by the densest bone (`pct.hu_reference_layers`, `pct.hu_reference_percentile`, `pct.hu_reference_resolution_mm`)
- Air-cavity check, on by default for layered media: reports and flags how much of the beam enters the sinuses or other internal air, from the SimNIBS charm air label (default), segmentation labels, or the CT (`headmodel.air_cavities`)

**Simulation**

- Backends: BabelBrain (acoustic and thermal), jwave (Helmholtz and time domain), and k-Wave elastic (*draft*); see [doc_backend.md](doc_backend.md)
- BabelBrain parity: build the PRESTUS grid at BabelBrain's own step (`babelbrain.match_grid`)
- BabelBrain on the PRESTUS medium: BabelViscoFDTD runs on the acoustic medium k-Wave sees, including smoothed properties and a pseudo-CT skull (`babelbrain.medium: prestus`)
- Head model: `pct.bone_open_mm` opens the CT-defined bone mask before the closing, and `pct.soft_hu_threshold` lets the CT separate soft tissue from water on the simulation grid, so a head model can follow density-threshold segmentation as in k-Plan
- Partial acoustic field of view along the beam (`grid.acoustic_fov_*`), a separate thermal grid, a free-water baseline on a box around the beam (`grid.water_baseline_grid`), and configurable time axis, PML, and band-limiting
- Heat deposition per tissue as `absorption_fraction` or `absorption_coeff`; pulse ramps (`timing.pulse_ramp_*`, applied to the pulse-average intensity, heating, and TIC); runs that start from an earlier T1w pressure map (`path.acoustic_input`)
- A reused acoustic cache from a different grid can be mapped onto the run's grid through world space instead of re-simulated (`io.acoustic_cache_remap`); a run without a source that cannot use the cache names the settings that differ from the run that wrote it
- Multi-transducer: transducers at different frequencies, async firing schedules (`timing.async_schedule`), firing groups (transducers that fire together within an async run, summed coherently, also on per-transducer grids), pulses that partly overlap at one frequency (interference while both are on, `timing.overlap_phase`), acoustic fields shared by repeated stances and kept for later runs, sonication sequences, and complex-field superposition (`superpose_complex_fields`)

**Studies, outputs, and reports**

- Study manifests (`study.manifest`, `prestus_study_start`), targeting accuracy against an MNI target or ROI (`analysis.*`), group MNI maps (`create_group_maps`), and a group report per condition (*experimental*)
- Uncertainty variants follow the ITRUSST Acoustic Properties Group bounds: Temp-H/L bound brain heating and Trans-H/L bound brain intensity. Uncertainty composes with multi-ISPPA and async runs; single bound values can be edited (`simulation.uncertainty_overrides`); the bound variants can run in a partial acoustic FOV, rescaled to the whole grid by a repeat of the default medium (`simulation.uncertainty_fov: partial`, *experimental*)
- Reports: four styles (`report.style`), a Methods section with copy buttons, calibration provenance, and segmentation, debug, and placement reports. New metrics: −3 dB focal widths, `MI` (free water), `Psptp`, and `riseT37_*`
- Transducer model NIfTI in T1 space for placement checks, written after medium setup (`io.save_transducer_nii`, layered media, annular transducers)
- Provenance: version, commit, and uncommitted-change state in parameters, logs, sidecars, and reports; JSON sidecars for stage outputs
- HPC: resource report with a CO₂e estimate, `hpc.exclude`, node-local disk for the k-Wave files (`hpc.tmp_gb`, `hpc.scratch_dir`), GPU requests that follow the solvers, and run statistics that refine run-time predictions

??? note "Full list of changes in v0.7"

    - k-Wave: a linear attenuation law (`fit_alpha_power: 0` with `alpha_power: 1` in every tissue) runs without dispersion on the MATLAB code types; the C++/CUDA code types, which do not accept it, stop at medium setup with an explanation.
    - ⚠️ **Changed:** heating heats each pulse train as one segment with the source weighted by PD/PRI, as k-Plan (new `thermal.pulse_model`, default `train_averaged`; `pulse_resolved` keeps the per-pulse segments; interleaved async schedules always resolve the pulses). With per-tissue masks no longer rebuilt at every segment, a 150 s protocol on a 0.5 mm head grid took 116 s instead of 88 min (400 s pulse-resolved); peak rises changed by at most 0.8 % and thermal dose by about 1 %. Set `thermal.pulse_model: pulse_resolved` to reproduce earlier heating results exactly.
- Placement: `placement.export_formats` writes the placed poses for Brainsight, 3D Slicer (ITK `.tfm`) and k-Plan (images without rotation, `.kps`, entry sheet).
- Nexus: *Import trajectory* (Placement) reads k-Plan plans (`.kdp`, `.kps`, `PLAN-*-RESULTS.h5`): transducer at the plan's bowl apex, focus at its target, focal depth from the selected transducer (needs `h5py`).
- Heating on the acoustic FOV (`grid.thermal_use_acoustic_fov: 1`) runs instead of stopping with an array-size error.
- Pseudo-CT: new `marsac-babelbrain` density mapping (`pct.mapping_density`), BabelBrain's CT-path density, for emulating BabelBrain with another backend.
    - Nexus: the page scrolls under plots; tabs share the window width with their full labels (in a window too narrow for them, the main tabs other than QuickStart show icons), and the label column of the cards follows the window width.
    - Nexus: Edit ▸ Undo / Redo cover the whole window in one history: config edits on every tab, placement viewer actions (whose arrows now use the same history instead of a separate one) and Study tab runs-table edits.
    - Installation: `prestus_setup` puts PRESTUS on the MATLAB path from any folder, fetches missing required submodules, and reports the installation; `prestus_check_install` prints that report on its own.
    - Nexus: `nexus/install.sh` (macOS, Linux) and `nexus/install.bat` (Windows) create or update the `prestus-nexus` environment, installing Miniforge if conda is missing; `prestus_nexus` finds and remembers the environment's Python instead of defaulting to `python3`.
    - Dependencies: export_fig v3.46 → v3.58 (bug fixes); the PlanTUS submodule moves to v2.0.0 with later fixes (5b69f66), and the PlanTUS setup also installs `pyyaml`, `PyQt5`, and `potpourri3d`.
    - `examples/setup_PRESTUS_env.sh` stops on macOS and Windows with a pointer to the SimNIBS installer, and downloads with `curl` when `wget` is missing.
    - Provenance: each run records its PRESTUS version, commit, commit date, and whether tracked files had uncommitted changes in `io.prestus_version` (parameter snapshots) and in the run log; sidecars and the simulation report show the same, instead of the `VERSION` tag alone.
    - Report methods text: an Output spaces paragraph states that MNI maps use the SimNIBS charm warp onto FSL's MNI152 template (`MNI152NLin6Asym`), citing Evans et al. (2012); the output docs name the template too.
    - Segmentation report: layered runs write `sub-XXX_layered_segmentation_report<affix>.html` with each step from the SimNIBS segmentation to the simulation medium (tissue mapping and smoothing, skull hole filling, eyes to water, CSF-based crop, grid crop): before, after, and changed-voxel slices through the focus, voxels moved between tissues, and tissue volumes (new `report.segmentation`, default on).
    - Debug report: with `simulation.debug`, `sub-XXX_<medium>_debug_report<affix>.html` lists every debug figure and volume of the run by stage, with caption and writing function, and opens with the error when the run fails (new `report.debug`). New `simulation.debug_stages`, `debug_figures`, and `debug_volumes` select what debug mode writes. 3D source and transducer figures are saved along the beam, from the side, and obliquely instead of as MATLAB `.fig` files, and the source mask is saved as a debug volume.

    - Transducers at different frequencies: each is simulated alone at its own frequency and medium, and the joint heat source uses the absorption at each frequency. Async coupling combines them as firing in turn; coherent coupling as firing together: complex fields summed per frequency, intensity and heat summed across frequencies, peak pressure from the summed complex fields.

    - Equipment listings (Nexus device lists, `load_equipment_config`) show each device model before its units.

    - Nexus: the *Sonication sequence* card (Sonication tab) sets up stances sonicated one after another, each with its own duration, pulse timing, and pause, and runs them as one joint async run or as a chain of sequential runs.
    - Async `custom` schedule: `transducer(i).timing.ptrd` sets each transducer's sonication duration.
    - ⚠️ Async `custom` schedule: the sonication ends at each transducer's start plus its PTRD, not at its last pulse train, so when PTRI is longer than PTD the end-of-run maps (`*_end`) are taken later. Re-run `custom` async runs whose end-of-run thermal results you report.

    - Nexus: the Transducer, Medium, Sonication, Solver, and Run tabs show the expected time of the run (each stage, the total, and GPU memory) for the current settings and GPU (`hpc.gpu`, or the local GPU for runs on this computer; `python -m prestus_nexus.workload`). `run_stats.jsonl` and telemetry record the simulated grid, time steps, and the duration of each stage, which refit the prediction (`python -m prestus_nexus.workload fit`). The new docs page *Run Time and Resources* has an interactive calculator of the same prediction.
- Run statistics record the modelled skull sound speed, density, and attenuation (all skull voxels, and on each beam path to the focus); telemetry sends their medians, with consent to the updated notice only.
    - Nexus calibration: the k-Wave drive check can be submitted as a cluster job (*Submit as a job*; `kwave_drive_submit`), here or over SSH, so its waiting runs on a compute node.

    - Nexus: *Transducer* tab, **Advanced**, **Equipment** sub-tab shows any library transducer and adds annular or matrix transducers and driving-system combinations (new driving systems included) to your equipment library, as `calibration_gui`'s **View geometry** and **New…** do. The Calibration page also calibrates any library transducer and driving system (**Any library transducer**), saving to a folder of your own library.
    - Nexus: the Results tab's *Report* list includes the placement report and group reports in `path.sim`; the Study tab opens the study's group report (*Open group report*); the Calibration page opens the library entry's HTML report (*Report*).
    - `predict_outputs` and the Nexus output preview list the placement report and the time-reversal tables.

    - Head model on fine grids: the rubber-wrap skull fill is evaluated only around voxels where brain meets skin instead of over the whole grid, its closing no longer needs the Image Processing Toolbox; the CSF-based crop dilates with one moving maximum per axis instead of a convolution whose cost grew with the cube volume (over an hour at 0.5 mm with `csf_expansion: 40`), with the same result.
    - ⚠️ Rubber-wrap skull fill: the radius is set in mm (`headmodel.skull_wrap_radius_mm`, default 5 mm, the old default of 10 voxels on the default 0.5 mm grid) instead of grid voxels, so the filled holes no longer depend on `grid.resolution_mm`; a config's `headmodel.skull_wrap_radius` is converted at its run's resolution. Rerun the grid stage of head runs at resolutions other than 0.5 mm that used the default radius.

    **Fixes**

    - Pseudo-CT: the rubber-wrapped skull mask (`pct_skullexpand`) is written with the header of `skull_mask.nii.gz` instead of without one.
    - Thermal grid: with `grid.thermal_resolution_mm` or `grid.thermal_fov_mm`, the pressure is resampled onto the thermal grid and the heating and CEM43 maps come back on the acoustic grid, instead of the thermal stage stopping on a size mismatch; both need a 3-D grid (`thermal_grid_setup:needs3D`).
    - Transducer library: per-depth refinements of a single-element transducer (e.g. BSonix) were read as one row of several elements and stopped the run (`load_transducer_from_library:correctionLength`).
    - Transducer library: a `transducer.serial` whose library entry is geometric (e.g. `CTX_500`) stopped the run with an unrecognized `elem_p_amp` field; it now runs with the nominal drive, scaled to `target_isppa_wcm2`.
    - Head model: voxels that no smoothed layer keeps are water in the medium masks and cached head models, as in the simulated medium; they, the head padding, and the fill of a rotated head model no longer widen the grid crop.
    - Simulation report: the debug section never found the debug folder; it is now a Diagnostics section linking the segmentation and debug reports.
    - ⚠️ **Changed:** `io.save_p_complex` with k-Wave (`kwave`, `kwave_elastic`; every `code_type`): the complex field is a least-squares fit over the last quarter period instead of an FFT of the last periods, which needed ~6× that recording in memory and ran out of 80 GB on 0.75 mm head grids. In these runs `p_max_all`/`p_min_all`, and the ISPPA and heating results that follow from them, are ± the fitted amplitude; rerun them to compare with runs without `save_p_complex`. Multi-frequency setups keep the steady-state window and stop above 8 GB of recording (`PRESTUS:pComplexTooLarge`).
    - `io.save_p_complex` with jwave (all solvers) and BabelBrain: writes the same `p_complex` file as k-Wave, so `superpose_complex_fields` combines fields from every backend; `p_complex_info.method` and `.phase_reference` state how each was computed.
    - Sequential chains: the follow-up's start state is the base run's T1w end map (heating, CEM43, CEM43 ISO), resampled into its own grid. The NIfTI stage is enabled for chained runs, and dispatch errors (`sequential_pipeline:missingThermalEnd`) when a map is missing, instead of starting from the water baseline.
    - `convert_final_to_MNI_simnibs`: `subject2mni` runs one volume at a time in the foreground. A failed conversion now raises `PRESTUS:MNIConversionFailed` and fails the job. The previous background path discarded output and reported success when volumes were missing.
    - Startup tool check (`check_tool_availability`): `prestus_pipeline_start` (before HPC submission) and `prestus_pipeline` (inside the job) check the external tools the enabled stages need: JVM, `sbatch`/`qsub`, the SimNIBS binaries (`charm`, `subject2mni`, `mni2subject_coords`) and `m2m/toMNI`, the k-Wave C++ binary for the selected `code_type`, and a GPU. A missing tool stops the run with `prestus:missingTools`, listing every missing tool. Another process holding more than half of a GPU's memory fails the check; more than 5% raises `prestus:gpuInUse`.
    - New `hpc.exclude` (SLURM): nodes passed to `sbatch --exclude`, e.g. a node whose GPU is held by a stale process.
    - HPC jobs request a GPU from the solvers that run and their compute backends (k-Wave `simulation.code_type`, `jwave.use_gpu`, `babelbrain.computing_backend`), not from `simulation.code_type` alone.
    - Nexus: Localite placement stays selectable while no pose location is set (needs setup) and is unavailable only when the set folders hold no poses for the subject; a placement field lists why each unavailable mode cannot run
    - ⚠️ **Changed:** async multi-transducer heating follows each transducer's firing schedule (new `timing.async_schedule`: `pulses` (default), `trains`, or `custom` with per-transducer `timing_offset_s` and `timing`) instead of applying all intensities during the pulses of the one `timing` block, which overestimated heating by up to N when the transducers shared that block's on-time; rerun the thermal stage of async runs
    - ⚠️ **Fixed:** MNI conversion through the SimNIBS coordinate fields (Localite MNI coordinates) interpolated with unmapped neighbours within one voxel of the field edge, pulling such points towards the origin by up to tens of mm; recompute MNI coordinates of affected recorded poses
    - MNI placement stops with `PRESTUS:placement:mniUnmapped` when SimNIBS cannot map a coordinate; the coordinate wrappers return NaN instead of 0 0 0 for such points
    - ⚠️ **Fixed:** time-reversal focal-shift predictions were clipped at ±10 mm; a peak on that boundary is now searched again within ±60 mm, and the summary CSV flags predictions on the search boundary (`shift_at_edge`). Rerun time-reversal runs whose predicted shift was 10–11 mm
    - ⚠️ **Fixed:** a partial acoustic FOV narrower than the transducer cropped the source, also with `acoustic_fov_diameter_mm: auto` (outer diameter without margin) and behind the bowl apex, where the box started; the diameter now includes 4 voxels per side, smaller values are enlarged with a warning, and the box starts 8 voxels behind the transducer. Rerun partial-FOV runs
    - ⚠️ **Fixed:** with `grid.mode: ras_plus` and a partial acoustic FOV, the FOV box of an oblique beam was aligned to the grid axes and could cut off part or all of the transducer; rerun ras_plus runs with a partial FOV
    - ⚠️ **Fixed:** the Rayleigh calibration model integrated over the projected ring area instead of the bowl surface, which lowered the on-axis focal pressure by 1–9 % (more for strongly curved bowls) and under-weighted outer rings; refit library entries and calibrations made with the default `rayleigh_phase_gain` preset or, in v0.6.1, with `calibration.forward_model: rayleigh`
    - ⚠️ **Fixed:** the rubber-wrap skull fill (`headmodel.skull_fill_method: rubberwrap`) ignored the brain and skin layers, so skull gaps where brain meets skin stayed open; rerun the grid stage of head runs
    - ⚠️ **Fixed:** with separate cortical and trabecular skull layers, filled skull holes were labelled trabecular instead of cortical; rerun the grid stage of runs with both layers
    - ⚠️ **Fixed:** the CSF-based crop dilated by about half of `headmodel.csf_expansion`; it now dilates by the full value along each axis, and the default is 20 mm (was 40), which keeps default head models (one 0.5 mm voxel layer differs on one side, ~20 mm from the CSF). Halve `csf_expansion` in project copies of `config_default.yaml` and study configs that set it
    - ⚠️ **Fixed:** with a pseudo-CT and `pct.mapping_density: k-plan`, the thermal skull density (1850 kg/m³) landed on non-skull tissue when layers were removed or reordered; rerun the thermal stage of such runs
    - ⚠️ **Fixed:** `export_babelbrain_inputs` wrote the pulse-train repetition duration (`ptrd`) as BabelBrain's `DurationUS`; it writes the thermal timing as `babelbrain_thermal_wrapper` uses it (`DurationUS` = `ptd`, plus `DurationOff`, `Repetitions`, `PauseBetweenGroupedSonications`). Re-export BabelBrain inputs used for standalone thermal runs
    - ⚠️ **Fixed:** jwave runs with a partial acoustic FOV padded the FOV by k-Wave's PML instead of the PML jwave applies; rerun them
    - ⚠️ **Fixed:** matrix transducers read from file with `project_on_new_ROC` mixed mm and m in the element depth; rerun these runs
    - Time reversal: the reverse run honours `grid.pml_size_range` and `grid.pml_alpha`, as the forward run does
    - HPC: `hpc_validate_parameters` returns the parameters, so batch jobs run with `simulation.interactive` off; qsub jobs write one log per job (`sub-XXX_qsub_output_<id>.log`); a second qsub wait in one MATLAB session no longer reports the job finished at once
    - `transducer_positioning_start`: `platform: matlab` runs locally without a scheduler, and waiting for the HPC job no longer errors after it ends
    - Pseudo-CT: the soft-tissue peak step runs with `matlab -batch`, so a MATLAB error ends the job instead of holding it to walltime; a missing `pseudoCT.nii.gz` error names the subject's own job
    - Placement: `neuronav_get_group_mean_mni` reads the `target` column that `neuronav_export_session_csv` writes; TUS_entry writes its target mask for int16 and uint8 T1s; heuristic placement stops with `PRESTUS:heuristic:noCandidate` when no candidate passes the expanded criterion
    - Head model: `smooth_img` box smoothing with a kernel below one voxel no longer fails, and anisotropic 2D box kernels keep the image size; `segmentation_run` no longer forces `charm --forcerun` when `path.seg` has no trailing separator
    - Uncertainty: the cleanup after the uncertainty report keeps the heating results (`heating_res`) and removes only the preprocessing caches
    - Thermal: the heating video works for runs with more than one transducer
    - Thermal: heating runs on 2D Cartesian grids; the thermal summary and the transducer mask of the thermal analysis no longer assume 3-D
    - Thermal: long protocols on large grids no longer run out of memory; the focal-plane time courses keep every k-th step within `thermal.max_focal_plane_gb` (default 2 GB; `HEAT.planetime`)
    - BabelBrain thermal: `babelbrain.computing_backend` names (`cuda`, `opencl`, `metal`) work as they do for the acoustic run, and `cpu` stops with the no-CPU-backend error
    - BabelBrain: `babelbrain.python_env` accepts the Python executable as well as the environment directory
    - Calibration: library entries store the version in `meta.prestus_version` and the commit in `meta.prestus_commit` (all writers); `extract_real_intensity_profile` applies its default `skip_front_peak_mm`
    - Plots and reports: `show_3d_head` with `align_to_focus: false`, `risk_color` for empty values, and the focus-box warning of `plot_t1_with_transducer` work; report flags are plain tinted boxes
    - `prestus_config_init(..., open_gui=true)` opens the config in the PRESTUS Nexus
    - ⚠️ **Changed:** T1w intensity, pressure, MI, and temperature maps are resampled with cubic instead of linear interpolation (`affine_resample_3d` now runs `'cubic'` instead of falling back to linear); linear underestimated focal peaks by up to ~3 % for narrow foci. Runs that start from a T1w pressure or temperature map read it with cubic too. Re-run the NIfTI export to update T1w maps
    - ⚠️ **Fixed:** a T1 with an oblique affine that is not stored RAS+ lost its rotation when PRESTUS reordered it to RAS+, which misplaced T1-space and MNI outputs and the conversion of RAS and Localite positions to voxels. Rerun runs whose T1 is both oblique and not stored RAS+
    - ⚠️ **Fixed:** `calibration_fit` without `FitOptions` fitted with `fit_transducer_model`'s amplitude-domain loss (`HuberDelta` 0.05) instead of the `rayleigh_phase_gain` preset (intensity, 0.8); the GUIs, the Nexus, and the library build were not affected. Refit entries from scripts that called `calibration_fit` without `FitOptions`
    - ⚠️ **Fixed:** the Localite planning image (`canonical_affine_transform`) and the heuristic placement's Localite coordinates flipped left and right for a T1 not stored RAS+; both are now reoriented to RAS+. Regenerate planning images and heuristic exports for such T1s
    - ⚠️ **Fixed:** heuristic placement exported transducer and target voxels rounded down while its table rounded to nearest; all exported positions use the nearest voxel, so positions can move by one voxel. Rerun heuristic placements
    - ⚠️ **Fixed:** the pseudo-CT soft-tissue peak was the lower edge of the modal UTE histogram bin instead of its centre; regenerate UTE pseudo-CTs
    - jwave: removed `jwave.precision`; jwave runs in single precision (float32/complex64), so `'double'` had no effect
    - Reports: banners, verdicts, safety cards, error sections, the licence note, and over-limit roster rows are plain tinted boxes instead of boxes with a coloured left bar (all skins)
    - The uncertainty medium configs cite the TUS Tissue Properties data set (Zenodo, https://doi.org/10.5281/zenodo.23811004)

    **Core / Config**

    - New warning `source_sensor_setup:staircasing` for head and phantom media with `grid.resolution_mm` above 0.75 mm: staircasing at tissue transitions may affect results
    - `load_parameters` accepts trailing `{key, value}` overrides, applied after all config layers
    - Unknown study-config keys warn with a suggestion (`prestus:unknownConfigKey`)
    - Configs from v0.5.0 on load as they are: renamed keys are converted to their current names with one warning (`prestus:renamedConfigKey`), and keys without an exact counterpart stop the run with the key to set instead (`prestus:legacyConfigKey`)
    - ⚠️ **Fixed:** a top-level `expected_focal_distance_mm` was read as a distance from the exit plane, while most earlier code read it from the bowl apex; it now stops the run. Set `transducer.focal_distance_ep` or `focal_distance_bowl` and rerun runs whose config set it
    - Grid and medium caches record the settings they depend on and refuse reuse when these differ (`io.cache_signature_check`: `error` | `warn` | `off`)
    - Stage outputs that a later run starts from (pseudo-CT, placement run, T1w pressure map) get a JSON sidecar that is checked on reuse
    - `path.acoustic_input` starts a run from an earlier run's T1w pressure map, rescaled to `transducer.target_isppa_wcm2`
    - `predict_outputs` / `print_predicted_outputs` show, before submission, which pipeline mode will run and which files it writes
    - ⚠️ **Fixed:** a run reusing another run's grid cache (`io.preproc_affix`) replaced its own transducer settings (drive, phases, ISPPA targets) with those of the run that built the grid; rerun runs on another run's grid whose transducer settings differ from it
    - ⚠️ **Fixed:** NIfTI outputs were truncated at the first dot in `io.output_affix` (e.g. `_th0.50mm` → `_th0`), so runs with dotted affixes overwrote each other's NIfTIs; rerun them
    - Fixed: `io.save_matrices: 1` had no effect with the default config
    - `io.save_kwave_medium` defaults to `0`: acoustic caches no longer hold the k-Wave medium (async per-transducer runs still store it)
    - Fixed: `get_alpha_coeff` read a non-existent medium field
    - Fixed: acoustic caches whose free-water ISPPA was stored as a gpuArray could not be read back
    - The parameter summary and HTML report state the free-water ISPPA a run models; a target that differs from the drive's calibration ISPPA, or the default target on an uncalibrated drive, warns (`prestus:isppa:targetVsCalibration`)

    **Transducer library & calibration**

    - `calibration_transducer` converts calibration settings from earlier versions and, for a name missing from the equipment library, says how to add it
    - Transducer library (`config/calibration/`): per-depth refinements of geometric steering, interpolated between measured depths; transducers resolved from their serial number via the equipment library. Fitting needs the Optimization Toolbox; loading does not
    - User calibration library: `refine_transducer_library` fits your own hydrophone profiles; set `transducer.library_path` to use it. An optional free-water k-Wave check corrects the drive per depth
    - Calibration at unmeasured depths from interpolated profiles (`interpolate_axial_profiles`); extrapolation is opt-in (`calibration.max_extrapolation_mm`)
    - `calibration.focal_reference` (`flhm_3db` default, `flhm_6db`, `peak`) sets which focal-depth convention a quoted depth refers to
    - Within-ring velocity profile (*ring profile*), optionally fitted with transverse or axial plane scans (`read_offaxis_scans`); off by default
    - Steering correction of matrix arrays from hydrophone focal peaks (`transducer.matrix.use_steering_correction`, default `true`; acts only when the library or config holds a correction)
    - Dephased conditions in library entries, selected with `transducer.condition` (*experimental*, not validated)
    - `calibration_gui`: add new transducers, predict profiles and calibrate at any desired depth, export profiles and settings
    - `import_babelbrain_weights` reads BabelBrain's native weights layout and reduces per-sub-element weights to rings
    - ⚠️ **Fixed:** `phase_optimization_annulus_full_curve` could accept a degenerate fit when the target profile was a row vector; re-run calibrations fitted to measured (row-vector) profiles
    - ⚠️ **Fixed:** interpolating a measured axial profile to an unmeasured focal depth put the focus at the same position for every desired depth; re-run calibrations at depths that were not measured
    - ⚠️ **Fixed:** the O'Neil forward model was evaluated 0.5 mm short of the requested depth; affects library entries fitted with the `oneil_phase_only` preset and `babelbrain.profile_calibration` (the default `rayleigh_phase_gain` preset is unchanged)
    - Fixed: `export_babelbrain_weights` wrote a file BabelBrain itself could not read
    - `min_foc`/`max_foc` are the hardware steering range of the device model; a unit without its own range uses its model's. A depth outside the measured range of a library entry warns that the entry extrapolates (`load_transducer_from_library:extrapolated`)
    - Imasonic PCD15287 / PCD15473: exit-plane distances confirmed from the manufacturer drawings; steering ranges 7.05–100.3 mm and 7.12–102.7 mm from the exit plane (were 8–110 and 15–110 mm)
    - Equipment entries can be marked `status: deprecated` or `excluded` (with `superseded_by`); they stay usable for earlier studies and warn when named, and excluded entries are not listed or built
    - `IS_PCD15473_01001_OPM` and `IS_PCD15473_01003_OPM` are removed; use `IS_PCD15473_01001` / `_01003` with `IGT_32_ch_comb_10_ch`, or add the OPM unit to your own transducer library
    - Equipment and shipped library entries are grouped by device model: `config/equipment/model-<MODEL>/tran-<SERIAL>_transducer.yaml`, driving systems in `config/equipment/ds/`, `config/calibration/model-<MODEL>/tran-<TRAN>_ds-<DS>_calib.yaml`. Units take geometry and steering range from their model entry. Flat folders (`<SERIAL>.yaml`, `<TRAN>_<DS>.yaml`) are still read, so your own equipment folders and `transducer.library_path` libraries work unchanged
    - Your own transducer libraries: equipment entries (`*_transducer.yaml`, `*_driver.yaml`), calibration entries, and `profiles/` in a folder of your own add transducers or replace shipped ones. PRESTUS checks `transducer.library_path` (study), new `path.transducer_library` (project), and the `PRESTUS_TRANSDUCER_LIBRARY` environment variable (machine), then the shipped library
    - `calibration_gui` saves new transducers into your own library and no longer writes into `config/equipment` (option `AllowShipped` for maintainers)
    - Calibration entries record the transducer geometry they were fitted for (`meta.geometry_fingerprint`); using one with a different geometry warns (`load_transducer_from_library:geometryMismatch`)
    - Drive limits: `transducer.max_elem_p_amp`, or a combo's `drive_limit` in the equipment entry, warns when the library drive for the target ISPPA exceeds what the driving system delivers (`resolve_transducer_from_serial:driveLimit`); estimates shipped for the IGT 128-channel driver
    - ⚠️ **Fixed:** calibration files (`*-F*mm-I*wpercm2.yaml`) recorded the default `target_isppa_wcm2` instead of their calibrated intensity, so a config built from one was rescaled to 30 W/cm². Set `target_isppa_wcm2` in configs built from older calibration files

    - A library transducer (`transducer.serial`) with an empty `target_isppa_wcm2` and no inline drive now stops with a message instead of failing during the simulation

    **Transducers & sources**

    - ⚠️ **Fixed:** matrix transducers with `grid.use_kWaveArray: 0` were not switched to the kWaveArray source as logged; rerun them
    - `elem_phase_deg` and `elem_phase_rad` that disagree warn (`prestus:phaseMismatch`); the source follows `elem_phase_deg`
    - Flat annular elements (`transducer.annular.curv_radius_mm: .inf`), kWaveArray source only
    - Matrix arrays: `matrix_shape.extract_from_file.rot_col` (per-element rotation), `matrix.elem_diameter_mm`, `matrix.elem_phase_override_deg`; Clover layouts with any number of leaves
    - ⚠️ **Changed:** Clover leaves are `clover.az_sep_deg` apart (default `360 / n_leaves`; was always 120°) and tilted by `clover.elevation_deg`; 1- and 2-leaf Clover layouts change
    - Large apertures and fine grids no longer run out of memory during source construction; GPU out-of-memory during source setup falls back to CPU with a warning
    - ⚠️ **Fixed:** matrix transducers (`define_here` grids) were not simulated as configured (element sizes in m instead of mm, width and height swapped, wrong cap curvature); re-run matrix simulations
    - ⚠️ **Fixed:** axisymmetric annular rings were placed in front of the bowl and aimed at the focus instead of the centre of curvature; rerun axisymmetric results of multi-ring transducers (single-element and 3D runs are unaffected)
    - ⚠️ **Fixed:** the axisymmetric kWaveArray source dropped part of its band-limited footprint, raising free-water focal intensity by 12–23 %; axisymmetric results at the drive level change, results scaled to a simulated free-water baseline much less (3D unaffected)
    - ⚠️ **Fixed:** 2D annular sources used the per-ring outer diameter vector instead of the full aperture; rerun 2D annular simulations
    - ⚠️ **Fixed:** in coherent multi-transducer runs, transducers 2..N were built from transducer 1's drive, phases, geometry, and frequency; re-run coherent multi-transducer acoustics (single-transducer, async, sequential, and BabelBrain runs are unaffected)
    - Fixed: flat matrix arrays (`is_curved: false`) failed; `extract_from_file.file_path` accepts absolute paths outside the MATLAB path
    - ⚠️ **Fixed:** without `focal_distance_ep` or `focal_distance_bowl`, the focal distance derived from voxel `trans_pos` / `focus_pos` rounded the T1 voxel size to whole mm (0.8 or 0.9 mm read as 1 mm, focal distance 11–25 % too long); rerun such runs on non-1 mm or anisotropic T1s
    - ⚠️ **Fixed:** curved matrix arrays without `dist_geom_ep_mm` put the exit plane at the bowl apex; it is now derived from the aperture as for annular arrays. Rerun curved matrix runs that left `dist_geom_ep_mm` empty
    - ⚠️ **Fixed:** `transducer.trans_pos_mm` / `focus_pos_mm` mapped 0 mm to the invalid voxel 0; they are now measured from the centre of voxel 1 (`round(mm / resolution_mm) + 1`). Subtract one voxel (`resolution_mm`) from configs that compensated for the old mapping

    **Placement**

    - `placement.mode: mni`: transducer entry point and focus given in MNI mm
    - `placement.mode: tus_entry`: scalp entry point from the TUS_entry R toolbox (Atkinson et al.), minimising coupling geometry error
    - `placement.manual.coordinate_space` (`voxel` default | `ras` | `mni`) for manual positions
    - Placement report `sub-XXX_placement_<mode>.html` for every mode (`placement.report`, default `true`); for PlanTUS, a criteria evaluation with a weight sweep (`placement.plantus.evaluate`, default `true`)
    - PlanTUS distance window can be set (`placement.plantus.min_distance_mm`, `max_distance_mm`, `optimal_distance_mm`)
    - PlanTUS headless and interactive modes write the same outputs, including `plantus_metrics.mat` (`placement.plantus.save_metrics`)
    - PlanTUS interactive mode: each click prints the vertex's criteria and rank; after closing Workbench, the transducer is shown at the chosen vertex and a dialog offers to use it or choose another; the skull-thickness map is loaded in the planning scene
    - `placement.plantus.interactive` selects interactive PlanTUS without `simulation.interactive` or the MATLAB desktop; without a display, PlanTUS runs headless
    - Fixed: the PlanTUS placement scene (`vtx<N>/scene.scene`) did not find the T1
    - ⚠️ **Fixed:** RAS → voxel conversion was one voxel off along each axis; placements derived from RAS/MNI coordinates (Localite, MNI, heuristic MNI target, PlanTUS, `coordinate_space: ras`) move by one voxel compared with earlier runs
    - ⚠️ **Fixed:** configured voxel positions ignored the reorientation of a non-RAS+ planning T1; rerun runs whose T1 is not stored RAS+
    - ⚠️ **Fixed:** Localite placement used the transposed instrument matrix (transducer at the image origin); rerun Localite placements in file mode and via `placement.localite.enabled` (the multi-target toolkit was not affected)
    - ⚠️ **Fixed:** with `tracker_to_bowl_mm` / `reference_distance_mm` set, and always in PlanTUS mode, the focus was too shallow by `curv_radius_mm − dist_geom_ep_mm`; PlanTUS placements computed before this fix have a focus that is too shallow by that amount
    - ⚠️ **Fixed:** TUS_entry placement ignored `focal_distance_ep` and put `focus_pos` `curv_radius_mm` below the scalp entry; rerun TUS_entry placements
    - ⚠️ **Changed:** all placement modes use one exit plane, `curv_radius_mm − dist_geom_ep_mm` in front of the bowl apex (as BabelBrain), and one gap to the scalp, `placement.skin_gap_mm` (default 5 mm; replaces `placement.mni.skin_gap_mm` and `placement.plantus.additional_offset_mm`). Heuristic and MNI placements measured their 5 mm gap from the bowl rim, so their transducers now sit by the rim-to-exit-plane distance further from the head (about 2 mm for a CTX-500); PlanTUS and TUS_entry placements had no gap and now sit 5 mm further out. Recorded placement runs are not reused; set `placement.skin_gap_mm: 0` to keep earlier PlanTUS and TUS_entry positions
    - ⚠️ **Changed:** PlanTUS skull thickness is now a placement criterion (`weight_skull_thickness` had no effect); with the default weights, placements can change
    - Fixed: PlanTUS interactive mode hung, did not find its own output, and crashed with a single calibration point; PlanTUS target masks failed on integer T1s
    - Fixed: automatic Localite placement (`placement.localite.file` empty) always errored; file mode now honours `position` / `target`
    - Fixed: `focal_distance_calculation` failed for Localite and PlanTUS placements
    - PlanTUS 2.0 support (Workbench-free). New `placement.plantus.skin_smoothing_iterations` (default 2, as in PlanTUS 2.0)
    - ⚠️ **Changed:** with PlanTUS 2.0, the skin surface is smoothed and the beam aims at the target centre; placements can differ from earlier PlanTUS runs. Set `skin_smoothing_iterations: 0` to keep the raw surface.
    - New `placement.plantus.vertex`: places a chosen skin vertex on the maps of an earlier run without recomputing them
    - New `placement.heuristic.candidate`: places a chosen candidate of the heuristic's table instead of applying its criteria
    - Fixed: heuristic placement without `transducer.focal_distance_ep` / `focal_distance_bowl` never finished; it now stops with a message
    - ⚠️ **Fixed:** the heuristic placement's Localite mm export (`loc_*_mm`) was one voxel high along each axis; re-export heuristic positions before entering them in Localite
    - Fixed: Localite placement with empty transducer geometry stopped with a dimension error instead of using `placement.localite.reference_distance_mm`
    - A transducer without `dist_geom_ep_mm` gets its exit plane at the bowl rim (outer diameter of the largest element, or `outer_diameter_mm` for matrix arrays), with a notice; the Nexus draws the same plane and no bowl for flat matrix arrays

    **Neuronavigation (Localite)**

    - Batch ingestion (`neuronav_ingest_batch`), session inventory, target-label maps (`neuronav.target_labels`) and target-name suggestions from planned targets
    - `neuronav_resolve_target` resolves targets over four confidence tiers, including a study-supplied lookup (`neuronav_read_target_lookup`)
    - Multi-coil ingestion via `neuronav.coil_map`
    - New `path.localite_*` and `placement.localite.*` keys in `config_default.yaml`; `placement.localite.session_folder_pattern` is now applied
    - `placement.localite.target` / `.position` take one entry per transducer for bilateral setups
    - ⚠️ **Fixed:** series QC statistics `rotation_dev_rad` and `position_dev_mm` were wrong (the latter in voxels); sidecar JSONs written before carry the old values
    - ⚠️ **Fixed:** `neuronav_convert_native_to_MNI` returned transducer MNI positions one MNI voxel (1 mm) off along each axis; recompute MNI transducer positions taken from it
    - Fixed: GUMMarkers files with `InstrumentMarker` elements failed to parse; T1 lookup in the conversion functions; `neuronav_convert_MNI_to_native` / `neuronav_convert_native_to_MNI` read an old path key

    **Head model & pseudo-CT**

    - Automated defacing via pydeface, optionally in the segmentation pipeline
    - ⚠️ **Fixed:** grid, medium, and acoustic caches, placement runs, and pseudo-CTs of a layered head were reused after the segmentation was re-made (charm re-run). They now record the segmentation's checksum and are checked on reuse; a defaced copy is made again. Re-run simulations of subjects re-segmented after their caches were written.
    - `segmentation.deface` and `io.deface_plots` are in `config_default.yaml`
    - `hpc.ld_library_path` and `placement.plantus.env_path` follow `startup.simnibs_bin_path` while empty: SimNIBS tools get the install's own libraries (`external/lib/linux`) prefixed to `LD_LIBRARY_PATH`, and PlanTUS runs in the install's Python; `env_path` is no longer required
    - Segmentation: SimNIBS 4 or later required (checked from the install); `segmentation.settings_file` (charm `--usesettings`) and `segmentation.noneck`; an existing segmentation made with another SimNIBS version or other charm settings is reused with a warning
    - Grid cropping includes the bowl of every transducer
    - `headmodel.reuse_across_placements` (default `false`) builds the smoothed head model once per subject and only rotates and crops it for each placement; in rotated grid modes the smoothing then precedes the rotation, so results differ slightly
    - New skull property mappings: `air-tissue` density, `webb` sound speed and attenuation, `pichardo` and `mcdannold` sound speed and attenuation
    - `pct.skull_mapping: custom` with `pct.custom_slope` / `pct.custom_intercept`
    - `pct.skull_mapping: carpino` is deprecated (`prestus:pct:deprecatedSkullMapping`); `custom` with slope -2194 and intercept 2236 gives the same line
    - `pct.density_calibration` selects the HU→density table of `mapping_density: k-plan` (`kplan_default`, `ucl_lowdose_80kvp`, `cph2025`, or a CSV); default `''` keeps results unchanged. PRESTUS warns when `treeby` is not paired with `ucl_lowdose_80kvp`
    - `pct_skullmapping_fit` exposes the `kosciessa` histogram fit and its settings
    - `pct_predict_skull_properties` predicts the skull properties a pseudo-CT run will model, without simulating
    - New `path.ct_pattern`: the skull mapping (`pct.enabled`) reads a measured CT, co-registered to the T1w, instead of generating a pseudo-CT from a UTE image
    - ⚠️ **Fixed:** `pct.skull_mapping: treeby` used the intercept 3247.9 HU instead of 3274.9 HU; treeby pseudo-CTs are 27 HU low in the skull, regenerate them
    - ⚠️ **Fixed:** the `kosciessa` fit took width crossings around the histogram maximum, not around the identified peaks; rerun `kosciessa` pCTs that warn or whose `pCT_skull_mapping.png` shows a near-zero trabecular mode
    - ⚠️ **Changed:** `smooth_img` pads by edge replication; smoothed masks and maps differ from earlier runs in the outermost voxels
    - Fixed: with `path.pct` empty the pseudo-CT was not found, and with `path.pct` set it was regenerated on every run; `pct_create_pseudoCT.sh` no longer renames files inside the SimNIBS folder
    - Fixed: defacing also covers `T2_reg.nii.gz`; skull hole filling works when bone has a single label; aligning to the natural focus without `focus_pos` failed

    **Acoustic / Thermal**

    - Each tissue sets heat deposition as `medium_properties.*.absorption_fraction` (default) or `absorption_coeff`; a config layer that sets one empties the other, and setting both in one layer errors
    - Logs, run summary, and HTML methods give each tissue's effective absorption coefficient and fraction; property maps include `absorption_np` and `absorption_fraction`
    - Partial acoustic FOV along the beam axis (`grid.acoustic_fov_length_mm`, `grid.acoustic_fov_diameter_mm`, `'auto'` from focal geometry) for deep targets; the thermal grid can match it (`thermal.match_acoustic_fov`)
    - Free-water grid `'auto'` sizing (`grid.water_margin_mm`, `grid.z_beyond_focus_mm`); matrix transducers placed by `trans_pos` / `focus_pos` keep `[144, 144, 400]`
    - `grid.water_baseline_grid: 'beam'` (default) runs the k-Wave free-water baseline on a box around the source and beam: the same free-field result at a fraction of the cost; `'full'` uses the whole grid
    - Time axis and steady state configurable: `grid.steady_state_transits`, `grid.steady_state_num_periods`, `grid.source_ppp`, `grid.sim_time_us`; PML and band-limiting options `grid.pml_size_range`, `grid.pml_alpha`, `grid.bli_tolerance`, `grid.bli_upsampling_rate`
    - Lower memory use for CW sources on `matlab_cpu` / `matlab_gpu`, also on axisymmetric grids
    - Axial-profile plot with tissue bands and the free-water O'Neil and k-Wave profiles; CSV gains `freewater_oneil_Wcm2` / `freewater_kwave_Wcm2`
    - Time reversal (`time_reversal.enabled`, *draft*): one reverse simulation per target gives refocusing phases, aberration metrics, and a ranking of candidate poses; `placement_search.jitter_mm` takes several distances, and `jitter_mode: scalp` slides the pose over the scalp and re-seats it; candidate poses keep `placement.skin_gap_mm` (heuristic, MNI, PlanTUS, TUS_entry) or the planned pose's gap
    - Pulse amplitude ramps in thermal simulations (`timing.pulse_ramp_dur`, `timing.pulse_ramp_shape`: `linear` | `tukey`)
    - New thermal metrics `riseT37_brain` / `riseT37_skull` / `riseT37_skin` (rise from 37 °C); new output columns `MI` (free water) and `Psptp`
    - ⚠️ **Fixed:** MI was computed from peak positive pressure instead of peak rarefaction pressure (IEC 62359); rerun acoustic analysis to update MI
    - ⚠️ **Fixed:** simulation end time was derived from the maximum (bone) sound speed, ending about 1.9× too early before steady state; rerun head acoustic simulations
    - ⚠️ **Changed:** a target ISPPA is never silently ignored: transducers without calibration ISPPA now get a free-water baseline and are scaled to the target (with the default 30 W/cm², one extra free-water simulation); set `target_isppa_wcm2: []` to report the drive level
    - ⚠️ **Changed:** coherent multi-transducer runs scale each transducer's drive to its own `target_isppa_wcm2` before the simulation; coherent multi-transducer runs with targets (including the default 30 W/cm²) give different fields, and older caches are refused when targets are set. ISPPA sweeps with coherent runs now error
    - ⚠️ **Changed:** `grid.default_dims: 'auto'` is the default; free-water results of annular transducers change only at discretisation level
    - ⚠️ **Fixed:** `io.save_p_complex` saved the complex field as a column vector, not on the grid, with its phase referred to the recording window; re-save fields you combine. New `superpose_complex_fields` combines per-transducer complex fields
    - ⚠️ **Fixed:** initial temperature (`temp_0`) was not loaded from cached heating files; rerun thermal analyses that reused cached heating files
    - Fixed: `fit_alpha_power: 0` failed in k-Wave; `convert_axisymmetric_to_3d` failed on label maps; three runtime crashes in `thermal_analysis`
    - Fixed: report and GUI texts for CEM43 now cite ISO 14708-3:2017 and Aubry et al. (2025) for the ISO form and state the k-Wave rule correctly; `thermal.cem43_iso` selects the form shown in the report (both are always computed)
    - ⚠️ **Fixed:** axisymmetric phantom runs that reloaded the acoustic cache crashed on a half-grid medium; the reload now expands grid, positions, and medium like a fresh simulation and honours `grid.axisymmetric_output`
    - ⚠️ **Fixed:** axisymmetric NIfTI outputs were cut at the grid centre instead of the beam axis when the transducer sat off centre. Re-export NIfTIs of off-centre runs
    - ⚠️ **Fixed:** phantom NIfTI outputs carried no affine and kept the simulation's voxel order; they now carry the phantom's affine in its voxel order (transposed or axial-first 2D phantoms, plane of a 3D phantom). Re-export phantom NIfTIs
    - ⚠️ **Fixed:** axisymmetric 2D output with the beam axis right of centre failed; beyond the simulated radius the plane repeats the edge value, as the 3D expansion does
    - Axisymmetric runs keep the bilateral transducer positions when reusing a grid cache, and the thermal stage uses the expanded positions on a separate thermal grid
    - Phantom voxels must be isotropic; an axisymmetric run on a 3D phantom cuts the plane at a 3D `transducer.trans_pos` (warns when the centre is ambiguous)
    - Focus positions with `focus_pos(1) > focus_pos(2)` are no longer turned into column vectors
    - Grid, medium, and acoustic caches of phantom runs record the phantom file; a reused acoustic cache is checked against the run's settings (`io.cache_signature_check`)
    - Removed `thermal.sensor_xy_halfsize` and `thermal.record_t_at_every_step`: neither had an effect on results, and configs that set them load with a notice

    **Multi-transducer**

    - Incoherent summation of independent-element arrays; exit-plane mask covers all transducers
    - Async mode (*experimental*): shared grid and medium stage before the per-transducer jobs, partial FOV supported, per-transducer ISPPA table in the report; HPC jobs chain via job dependencies, and `[job_id, jobs] = prestus_pipeline_start(...)` returns the job IDs per stage
    - ⚠️ **Fixed:** `io.overwrite_files: 'always'` was ignored by chained pipelines (async, multi-ISPPA, uncertainty), so stale outputs could block the joint job and its report; rerun chained runs written into folders with earlier outputs

    **Backends**

    - **BabelBrain** as acoustic and thermal backend (`acoustic_backend: babelbrain`, `thermal_backend: babelbrain`), with a CT path (`babelbrain.ct_mapping`), internal air regions (`babelbrain.air_regions`) and free-water FDTD by default (`babelbrain.use_rayleigh_for_water: false`). One transducer per run, beam along the grid z axis; BabelBrain thermal uses its own tissue thermal properties. Defaults: `babelbrain.ppw: 12`, `babelbrain.q_correction` 1.0 on the non-CT path (3.0 with CT), `babelbrain.tight_narrow_beam` off, `babelbrain.profile_calibration` off
    - BabelBrain 0.8.8 (`external/BabelBrain`). Update the BabelBrain environment: `pip install babelviscofdtd==1.2.7 psutil`
    - BabelBrain simulates the PRESTUS transducer's ring geometry, on which its calibration was fitted, unless `babelbrain.tx_system` names a BabelBrain preset (resolved from BabelBrain's registry, incl. `R16081`, `TheraPhas`; `REMOPD` still accepted)
    - BabelBrain `Single` and `BSonix` presets take their aperture and focal length from the PRESTUS transducer; with `thermal_backend: babelbrain`, multi-ISPPA thermal jobs run one after another, as they share one BabelBrain thermal file
    - **jwave** as acoustic backend (`acoustic_backend: jwave`), Helmholtz and time-domain solvers; `jwave.helmholtz_operator: conservative` (default) corrects transmission through skull, `jwave` restores the built-in operator. kWaveArray sources supported; no axisymmetric grids
    - **k-Wave elastic** (`acoustic_backend: kwave_elastic`, *draft*, not validated): shear waves in the skull with properties from the new `elastic:` section; MATLAB CPU/GPU, no axisymmetric grids
    - New `babelbrain.shear: prestus_elastic` runs BabelBrain with elastic bone, using the `kwave_elastic` skull shear speed and attenuation at the drive frequency (default `off`: fluid bone)
    - Backend settings are validated at config load; the HTML report has backend sections and backend-specific methods text

    - BabelBrain's narrow-beam domain uses `grid.z_beyond_focus_mm`, the same axial extent past the focus as the k-Wave and jwave FOV

    **Advanced modes (uncertainty, multi-ISPPA, sequential)**

    - Uncertainty composes with multi-ISPPA sweeps and async multi-transducer runs, with one report per sweep point; uncertainty is the outer mode on every platform
    - Uncertainty mode runs five variants: default, Temp-H / Temp-L (the ITRUSST tissue-property bounds that maximise / minimise brain heating) and Trans-H / Trans-L (the acoustic bounds that maximise / minimise brain intensity; acoustic only). See [doc_uncertainty.md](doc_uncertainty.md)
    - `simulation.uncertainty_overrides` edits single bound values of the variants; the Nexus *Values per variant* table edits them while unlocked and highlights edited values
    - ⚠️ The liberal / conservative variants are replaced: new medium bounds and output affixes (`_desc-temphigh`, `_desc-templow`, `_desc-transhigh`, `_desc-translow`), and `options.liberal_config`, `options.conservative_config`, and `options.affixes.liberal` / `.conservative` now stop with `uncertainty_pipeline:legacyVariant`. Re-run uncertainty runs
    - The uncertainty report ranges acoustic metrics over the Trans pair and thermal metrics over the Temp pair, and lists every variant's values
    - ⚠️ **Fixed:** uncertainty safety cards and range-table rows were coloured by the conservative variant; they now follow the highest value of the default and the metric's pair. Regenerate uncertainty reports
    - ⚠️ **Fixed:** uncertainty variants shared one medium cache file and overwrote each other with `io.save_medium_cache: 1`; rerun uncertainty runs that used it
    - Fixed: the uncertainty report is named after the default variant's affix, so runs with different `io.output_affix` no longer overwrite each other's report; HPC resume finds existing reports and variant CSVs
    - Fixed: the uncertainty CSV is saved without an extra `tabular/` folder; variant and report jobs no longer run their own free-water baseline
    - ⚠️ With `pct.enabled: 1`, the uncertainty bound variants use a homogeneous skull layer on the pseudo-CT run's head geometry; only the default maps the skull from the pseudo-CT. Re-run pseudo-CT uncertainty runs
    - The uncertainty report of a pseudo-CT run has a Pseudo-CT Skull vs. Bounds section: the default run's skull property histograms against the bound variants' skull values
    - Fixed: `multi_isppa_pipeline` warned about missing thermal output for every intensity (wrong CSV name)

    **Study setup & group analysis**

    - Study manifest (`study.manifest`, `study.name`): one row per run with subject, condition, and any config key as a column; `prestus_study_start` runs or submits every row and builds a study report
    - Targeting accuracy (`modules.run_targeting_analysis`, `analysis.target_mni_mm`, `analysis.target_roi_mni`): focal centroid and peak against the planned focus, an MNI target and an ROI
    - Group MNI maps (`create_group_maps`): per-condition focal-volume overlap, mean intensity, and heating
    - Group report per run (subject × condition), with condition comparison, targeting, group maps, and outliers sections
    - ⚠️ **Fixed:** `analysis.mni_warp_method: native` placed MNI outputs one subject voxel off along each axis; rerun MNI exports made with it

    **Outputs & reporting**

    - Metrics CSV and report: −3 dB focal widths through the focal peak, `focal_fwhm_axial_mm` along the transducer axis and `focal_fwhm_lateral_min_mm` / `focal_fwhm_lateral_max_mm` across it
    - Reports redesigned (`report.style`: `instrument` default, `classic`, `compact`, `paper`); `preview_report_styles` re-renders an output folder in every skin
    - Manuscript-ready Methods section with copy buttons and stated tissue properties
    - Report shows the calibration provenance of the transducer (library entry, k-Wave-verified, or analytical only)
    - ⚠️ **Changed:** focal-distance columns name their reference: `real_focal_distance_mm` → `real_focal_distance_bowl_mm` (bowl apex), next to `real_focal_distance_ep_mm` (exit plane). Reports and the Nexus read older CSVs; update your own scripts that read `real_focal_distance_mm`
    - Safety dashboards show `MI` and `MItc` side by side
    - Run statistics in `<path.sim>/run_stats.jsonl`, used for resource estimates in the Nexus and `python -m prestus_nexus.estimate`
    - Fixed: Methods references for Martin et al. (2024), Pichardo et al. (2017), and Aubry et al. (2025) had wrong details; regenerate reports to update them

    **HPC**

    - Post-job resource report via SLURM `sacct`, with peak memory (including BabelBrain) and a CO₂e estimate
    - `hpc.posthoc_water_timelimit` / `hpc.posthoc_water_memorylimit`; `hpc.depend_type` for dependent jobs
    - Fixed: every Snellius SLURM submission failed (wall-time check)
    - Fixed: qsub scripts wrote `walltime` as character codes and requested GPU nodes for CPU-only jobs; SLURM job IDs are read with `sbatch --parsable`
    - Fixed: `io.overwrite_simnibs` never passed `--forcerun` to charm; `simnibs_t1_to_mni` respects `overwrite_files`
    - ⚠️ **Fixed:** under `io.overwrite_files: never` or `ask`, a recomputed acoustic or thermal stage kept existing NIfTIs but rewrote the CSV, mixing outputs of two runs. Re-run output folders that were simulated more than once under `never`

    - Fixed: PBS jobs loaded MATLAB R2023b, which does not support the CUDA 11.2 of PBS nodes; they now load R2022b

    **Nexus & GUIs**

    - Label text in the Nexus can be selected and copied
    - Nexus Medium tab: a *SimNIBS report* card below the segmentation shows charm's report; the tissue property table gives each tissue's volume in the segmentation and greys tissues it does not model
    - Fixed: the Nexus suggested `transducer.target_isppa_wcm2` for the old `set_intensity_w_per_cm2`, which is the drive's calibration ISPPA (`transducer.tx_calibration.isppa_wcm2`)
    - **PRESTUS Nexus** (`nexus/`): standalone PyQt6 app to build, check, and run study configs, with tabs per workflow step, a transducer panel backed by the library, a placement viewer with beam checks, Calibration and Localite pages, a study manifest editor, a Results tab, and cluster submission (also over SSH). MATLAB entry point: `prestus_nexus`. See [doc_nexus.md](doc_nexus.md)
    - The Nexus Pseudo-CT sub-tab plots the in-skull property distributions against the band of the uncertainty bound variants' skull values, with the fraction of skull voxels outside it; **Predict skull properties** predicts them before a run
    - The Nexus Simulation tab is split into Medium, Sonication, and Solver tabs. The Medium tab has the medium choice above its Head model, Acoustic properties, and Pseudo-CT sub-tabs. Medium properties show as one tissue × property table for the tissues of the chosen medium, editable after unlocking in advanced mode, with a plot of the values against the uncertainty-variant bounds. The Pseudo-CT sub-tab shows the equations of the selected mappings. More fields appear only when their mode or switch is on (Localite folders, pseudo-CT mapping, thermal model, time reversal)
    - Nexus Transducer tab: the cross-section sits under the device panel; Geometry (collapsed for library devices) and Drive and steering are separate cards; calibration coverage moved to the Calibration sub-tab. Non-k-Wave solvers are marked experimental. The CEM43 form is a selector that shows its formula
    - The Nexus flags a grid coarser than the points-per-wavelength minimum (`grid.min_ppw`) before a run
    - `prestus_gui` is removed in favour of the PRESTUS Nexus. `neuronav_gui` and `calibration_gui` remain supported.
    - `neuronav_gui`: stand-alone Localite ingest GUI with series QC, batch window, and MNI view of batch targets
    - `calibration_gui` and the Nexus Calibration page expose every calibration fit option and save how an entry was fitted (`meta.fit`)
    - Nexus Calibration page: the library entry's fit quality, drive, and phases are shown before MATLAB starts; the Transducer tab has a calibration coverage card
    - Nexus: tissue-property card is read-only until unlocked; with uncertainty on, a card lists and plots the three variants' properties
    - Nexus: a Start tab sets up a first simulation in four steps (T1, transducer, MNI target, run); pages without their inputs say what is missing; the Issues list opens from the header (View ▸ Issues…, Ctrl+Shift+I)
    - Nexus: with `placement.mode: plantus`, the Placement card runs only the PlanTUS placement, choosing the position in Workbench or automatically, and opens the result in the Viewer
    - Nexus: Placement ▸ PlanTUS shows the PlanTUS maps on the skin; click a position to see its criteria and the transducer there, and place it
    - Nexus: PlanTUS and TUS_entry placement-only runs are submitted as cluster jobs when the Run tab submits to a cluster
    - Nexus: the placement method is chosen at the top of the Placement tab; new TUS entry sub-tab; the sub-tab of the method in use is labelled in bold
    - Nexus: Placement ▸ Heuristic plots and lists the heuristic's candidates; click one to see its criteria and place it
    - Nexus: the segmentation is a step of its own: status and *Run segmentation* (segmentation-only job) on the Medium tab, the Placement tab (when the method needs it) and QuickStart
    - Nexus: Results opens on an ITRUSST reporting overview (Copy as Markdown); the maps moved to a *Spatial* sub-tab
    - Nexus Results overview: axial and lateral −3 dB focal widths (from the CSV, else measured on the run's T1w intensity map) and ISPTA per tissue
    - Nexus: the 3D view has opacity sliders for scalp, brain, and mask, marks the target with a crosshair, and draws the PlanTUS/TUS_entry target mask
    - Nexus: the Segmentation card selects the SimNIBS environment from the installs found (conda environments, install folders, `/opt`; with SSH submission also on the cluster, including environment modules); settings linked to it follow a change of environment, shows its version and the existing segmentation (with its charm report), and has the re-run switch and a custom charm settings file
    - Nexus: Medium cards regrouped (Segmentation refinement, Smoothing, Grid refinement); *Tissue layers* is advanced and locked; section cards carry icons
    - Nexus: formulas are drawn as vector graphics at text size; figures follow the dark theme; the Start tab is now *QuickStart*
    - Nexus: Results ▸ *Follow-up…* derives a run from an existing one (new intensity, new timing / thermal only, new placement on the same head, a time-reversal pose, or a sequential run with heat carry-over), checks the parent's outputs, and records the parent in the new config; a *Poses* sub-tab lists time-reversal poses
    - Nexus: Run ▸ *Pipeline options* passes sequential chains, per-stage time and memory limits, and the summary report to `prestus_pipeline_start`, also for SSH submission
    - Nexus: Placement ▸ *Transducers* (two or more transducers) lists the transducers and shows the grid orientation (`grid.mode`) with figures and each beam's tilt; the viewer can draw the grid z axis
    - Nexus: every figure has an export button (SVG, PDF, PNG at 300/600 dpi, copy); result maps default to the magma colormap
    - Nexus: Help menu and ↗ buttons open the documentation page of each tab and card; Edit ▸ Undo / Redo; Ctrl+1…9 switch tabs; Tools ▸ *Usage statistics…* records the telemetry decision
    - Nexus: new checks for `path.acoustic_input`, multi-transducer runs with water or axisymmetric grids, heating with zero pulse timing, the thermal field of view, and elastic shear points per wavelength
    - Nexus: saves are atomic; an unexpected error shows a dialog instead of closing the app; cluster job status and MATLAB requests no longer freeze the window, and MATLAB requests time out

    **Telemetry**

    - `PRESTUS_TELEMETRY` (`off` / `on`) sets the telemetry decision for the current process
    - Fixed: `matlab -batch` and SLURM/PBS sessions no longer prompt for consent; they print the notice and record nothing

    **Documentation**

    - Navigation regrouped by workflow; new pages *Your First Simulation*, *Configuring a Study*, *How PRESTUS Works*, *Concepts*, *Glossary*, *k-Wave*, *Reading Your Results*, and *Troubleshooting*
    - Fixed: `config_tutorial.yaml` nested transducer positions under an unused key and referenced a missing element file; the tutorial matrix config gave the element size in m (5.317e-3) instead of mm

---

## v0.6.1
*(2026-05-22)*

Extends the calibration pipeline with a unified entry point, parametric model store, and HPC dispatch. Introduces async multi-transducer simulations, a sequential simulation dispatcher, PlanTUS-based automated placement, and RAS as the canonical coordinate space. Drops the Image Processing Toolbox requirement.

#### Core / Config
- MATLAB Image Processing Toolbox no longer required — morphological ops and resampling use base-MATLAB equivalents
- `transducer.target_isppa_wcm2` migrated from calibration struct to transducer struct
- Optional `transducer.name` field added to transducer config

#### Calibration
- Unified `fit_source_parameters` entry point with HPC dispatch, multi-intensity mode, and parametric model store
- Geometric steering mode with per-element hardware correction
- Multiple analytical forward models — O'Neil and Rayleigh
- Transducer library (draft)
- Free-water correction moved post-analytical to make optional
- Free-water validation can be deactivated or limited to first/last target ISPPA
- Pluggable analytical forward model interface
- Option to use manufacturer phases to initiate global search (`opt_use_initial_phases`)
- ⚠️ **Fixed:** amplitude scaling bug present since early PRESTUS versions — correction factor is the square root of intensity ratios, not the ratio itself

#### Placement
- PlanTUS placement mode wraps the external PlanTUS optimiser for automated transducer targeting — **experimental** (see [doc_placement_plantus.md](doc_placement_plantus.md))
- T1 overlay plot auto-generated after heuristic and PlanTUS placement
- Browser-based transducer alignment viewer (*experimental*)

#### Acoustic/Thermal
- Async multi-transducer pipeline — **experimental** (see [doc_multitransducer.md](doc_multitransducer.md))
- Sequential simulation dispatcher: runs multiple parameter configurations in series, chains thermal simulations, and produces a consolidated multi-run report (see [doc_advanced.md](doc_advanced.md))
- ⚠️ **Fixed:** sequential simulation produced incorrect outputs prior to this release due to mismatching grid->T1w transforms and back of the adopted heatmaps; the heatmaps were not interpolated with the baseline free-water value leading to interpolation artefacts
- Uncertainty pipeline produces uncertainty bands across sequential multi-run parameter sweeps
- Complex pressure field optionally saved (`save_p_complex`) (*experimental*)

#### Head pipeline
- RAS established as canonical entry-point coordinate space; RAS+ orientation enforced at NIfTI load time; `ras_plus` grid mode rescales without rotation
- `tissuemask_binary` exposes a `water` mask field
- `smooth_img` gains an `off` method (no smoothing); exposed in GUI

#### Phantom
- `trans_pos_mm` / `focus_pos_mm` config fields allow resolution-independent transducer and focus position specification
- ⚠️ **Fixed:** phantom NIfTI outputs used incorrect voxel dimensions and inconsistent world-space orientation across pipeline stages

#### NIfTI / Transform
- Batched `tformarray` with optional parallel workers for MNI transforms
- Explicit FOV mask used in MNI outputs instead of zero sentinel (unreliable for pCT Hounsfield values)
- ⚠️ **Fixed:** pCT path layout corrected; NaN padding used instead of zero fill

#### Telemetry
- System resources reported at run start: OS, CPU model/cores, RAM, GPU, CUDA version

#### Reporting
- Public usage statistics (see [doc_usage_statistics.md](doc_usage_statistics.md))
- Group HTML report — `prestus_group_report_start` / `generate_group_report` aggregate per-subject CSVs and PNGs (across all subjects sharing a `medium` + `output_affix`) into a self-contained HTML with a mean ± SD exposure dashboard, subject roster, split acoustic (Intensity, MI) and thermal (maxT, riseT, CEM43) box plots, per-subject cards, and a client-side subject filter. Subjects are auto-discovered via `discover_group_subjects`; generated iteratively after each subject finishes when `modules.generate_group_report = 1`, or on demand. See [doc_group.md](doc_group.md#group-html-report).

---

## v0.6.0 
*(2026-04-30)*

Introduces matrix transducer support, an uncertainty quantification pipeline, and a rudimentary GUI. The pipeline is further integrated with automated pCT generation, placement dispatch, multi-ISPPA targeting, and a BIDS-inspired output structure. Parameter names and output paths have changed throughout.

> ⚠️ **Migration note:** Several parameter names and all output paths have changed in this release. Existing configs and post-processing scripts will need updating.
>
> **Parameters:**
>
> - `parameters.submit_medium` renamed to `parameters.platform`
> - Segmentation smoothing now specified in FWHM mm (parameter name changed)
> - Transducer field suffixes renamed consistently — check transducer config fields
> - Further parameter renames throughout — consult `config_default.yaml` for current names
>
> **Output paths:**
>
> - `parameters.io` fields renamed: `*_dir` → `dir_*`, `*_filename` → `filename_*`
> - Preprocessing cache files follow a consistent `res_N_descriptor` naming scheme
> - `space-` prefix removed from all NIfTI output filenames
> - `output_affix` now appended at the end of all output filenames
> - HPC log directory renamed: `hpc_log/` → `log_hpc/`
> - Outputs reorganised into BIDS-inspired typed subdirectories (`acoustic/`, `thermal/`, `cache/`, `debug/`, etc.)
>
> **Directory and config naming:**
>
> - `toolboxes/` directory renamed to `external/`
> - Config files consistently prefixed with `config_`

#### Added

- [**feature**] Uncertainty quantification pipeline (see [doc_uncertainty.md](doc_uncertainty.md))
    - Three-variant simulation workflow (default / liberal / conservative medium properties) producing a unified HTML report with per-metric safety cards colour-coded against ITRUSST non-significant risk limits, acoustic and thermal range tables, side-by-side intensity and temperature maps, and an SVG temperature timeseries with uncertainty band
    - Supports MATLAB (sequential) and HPC (parallel, scheduler-dependent) execution with automatic resumption from partial runs
- [**feature**] Matrix transducer support (see [doc_transducer.md](doc_transducer.md), [PR #117](https://github.com/Donders-Institute/PRESTUS/pull/117))
    - Flexible element layout options: rectangular grid, Fibonacci, Fermat spiral, and Clover multi-array configuration
    - Element positions can be defined inline or loaded from file
- [**pCT**] pCT generation fully automated inside the pipeline; a T1w image, UTE image, and mapping algorithm are sufficient — pre-generation remains possible for debugging; pCT section included in HTML report when used
- [**calibration**] In-pipeline amplitude calibration to target ISPPA
- [**feature**] Multi-ISPPA targeting — run acoustics once, scale to N independent thermal targets via post-hoc pressure scaling (see [doc_multiintensity.md](doc_multiintensity.md))
- [**thermal**] Dual CEM43 output: both kWave and ISO formulations propagated through pipeline, CSV, NIfTIs, and report
- [**report**] ISO CEM43 shown in Exposure Dashboard and methods section when enabled
- [**feature**] Placement dispatch — unified entry point for manual, Localite, and heuristic transducer positioning (see [doc_placement.md](doc_placement.md))
- [**feature**] `prestus_config_init` — bootstraps a project-specific config directory and `config_default.yaml`; `load_parameters` now resolves the project config automatically so PRESTUS can be called from any working directory
- [**neuronav**] Unified Localite pipeline via `localite_matrix_to_positions`; `InstrumentMarker` support added; `tracker_to_bowl_mm` parameter
- [**feature**] PRESTUS GUI — tabbed parameter editor and simulation launcher (see `doc_gui.md`)
- [**feature**] `prestus_group_start` — new entry point for group-level MNI-space analysis
- [**feature**] `segmentation_only` pipeline mode
- [**report**] MItc (mechanical index transcranial) reported across all intracranial tissues
- [**thermal**] `simulation.log_thermal_updates` — suppress verbose per-step k-Wave progress output (default: off)
- [**refactor**] `simulation_nifti` split into per-pipeline-stage calls; bone mask separated from pseudoCT NIfTIs
- [**telemetry**] Opt-in anonymous usage telemetry (disabled by default; see [doc_telemetry.md](doc_telemetry.md))

#### Changed

- [**I/O**] Subject-specific output subfolders enforced
- [**output**] Intensity output naming cleaned up: redundant `max` labels removed; naming aligned with standardised pressure and intensity metric conventions
- [**grid**] `pml_size` defaults to `'auto'`; effective PML size resolved and tracked in log
- [**config**] `simulation.debug` defaults to `0`; property-map and MNI NIfTI export flags added to `config_default.yaml`
- [**validation**] Function inputs consistently validated via MATLAB `arguments` blocks; `help <function>` exposes header information inline
- [**report**] Safety Dashboard renamed to **Exposure Dashboard**; water-medium limits now populated; pressure units corrected; dynamic pressure display added; CEM and temperature removed from water-medium report

#### Fixed

- [**paths**] `load_parameters` resolves `config_default.yaml` via an absolute path, preventing stray output folders
- [**paths**] Cross-platform SimNIBS path compatibility
- [**transducer**] Annular array element format and default thickness corrected
- [**nifti**] Property map export: corrected datatype mismatch, scalar `alpha_power` handling, and `FillValue` assignment
- [**hpc**] Transducer index suffix omitted from positioning filenames for single-transducer runs
- [**plots**] T1 intensity overlay orientation corrected using NIfTI sform; deprecated `LineSmoothing` property removed
- [**logging**] Post-hoc water simulation writes to its own independent log file; PPW warning no longer emits a backtrace
- [**debug**] Debug directory creation and plots skipped when `simulation.debug = 0`

---

## v0.5.0
*(2026-02-27)*

Major restructuring of the pipeline into discrete processing stages. Adds parallel multi-transducer support, an HTML simulation report, C++ backend, skull rubber-wrap, and heuristic transducer placement. Parameter names have changed and existing configs will need updating.

> ⚠️ **Migration note** — configs from v0.4.x will need updating:
>
> - `layer_labels` renamed to `layers`
> - Transducer and focus position fields renamed to `trans_pos` / `focus_pos` consistently
> - `grid_step_m` and `default_grid_size` removed
> - Grid expansion parameter renamed
> - Attenuation parameters renamed — consult `config_default.yaml`
> - Layer indexing field migrated from `layers` to `medium` throughout
> - Pipeline refactored into discrete processing stages — consult the documentation for current function entry points

#### Added

- [**feature**] C++ backend CPU & GPU support by @jkosciessa
- [**feature**] Parallel multi-transducer support for layered simulations by @dreamstimlab (see [PR #100](https://github.com/Donders-Institute/PRESTUS/pull/100))
- [**feature**] HTML simulation report — self-contained per-subject summary of acoustic and thermal outputs by @sirmrmarty
- [**feature**] Heuristic transducer placement — automated skull-surface search with configurable positioning criteria and ear exclusion (see [PR #109](https://github.com/Donders-Institute/PRESTUS/pull/109))
- [**skull**] Rubber-wrap skull mask inflation by @meijer-s (see [PR #103](https://github.com/Donders-Institute/PRESTUS/pull/103))
- [**thermal**] Flexible and standardised protocol timing specification aligned with TUSCalculator; protocol plot output
- [**thermal**] Additional thermal timeseries outputs (max per tissue, heating at focus)
- [**feature**] Per-step pipeline benchmarking: wall time, peak RAM, and disk usage logged at each stage
- [**grid**] Minimum PPW validation at source setup (`grid.min_ppw`, default 6); warns with the maximum `resolution_mm` that would satisfy the requirement
- [**calibration**] `calibration.opt_phase_precession` — optional constraint on element phases during amplitude/phase optimisation
- [**neuronav**] Updated Localite GUMMarkers parsing
- [**doc**] GitHub Pages [documentation](https://donders-institute.github.io/PRESTUS/) — parameters, functions, installation, quick start, and processing steps

#### Changed

- [**refactor**] Pipeline and core processing steps refactored; functions organised by processing stage
- [**axisymmetric**] Axisymmetric acoustic simulations now default to 3D thermal output

#### Fixed

- [**pCT**] kPlan skull density and sound speed regularised to water minimum
- [**acoustic**] GPU array correctly moved to CPU before NIfTI export
- [**calibration**] Unoptimized parameters no longer passed to Oneil recomputation
- [**savemat**] k-Wave source matrix saving can now be correctly deactivated; water simulations force-save matrix outputs

#### Deprecated/Removed

- [**feature**] Explicit `skull`-only simulation mode removed; equivalent behaviour is achieved by removing layers from the simulation setup

---

## v0.4.2 
*(2025-12-04)*

Refactors the calibration and acoustic profiling pipeline and introduces multi-element axisymmetry support.

#### Added

- [**calibration**] Acoustic profiling refactor with documented standalone calibration pipeline (see [doc_calibration.md](doc_calibration.md))
- [**calibration**] Multi-element axisymmetry support (experimental)

#### Changed

- [**default**] kWaveArray used by default

---

## v0.4.1 
*(2025-08-26)*

Post-release additions to v0.4.0: heat trace outputs, Localite neuronavigation integration, and 3D thermal as the default for axisymmetric runs.

#### Added

- [**thermal**] Heat trace saving
- [**logging**] MATLAB log written to output directory
- [**neuronav**] Localite coordinate processing

#### Changed

- [**axisymmetric**] 3D thermal simulation is now the default for axisymmetric acoustic runs

---

## v0.4.0 
*(2025-06-04)*

This release features major updates. Please consult the documentation for current parameter choices. Some features remain unstable and new bugs may pop up. Let's remain vigilant and try to quickly introduce hotfixes in those cases.

- [**Feature**] Acoustic profiling (beta) @MaCuinea
- [**Feature**] pCT->acoustic skull medium incl. multiple mapping algorithms (beta, doc) @eleonoracarpino @jkosciessa
- [**Feature**] Run consecutive simulations @KTZ228
- [**Feature**] GPU updates (doc) @jkosciessa
- [**Feature**] CEM43 ISO definition, see doc @jkosciessa
- [**Feature**] Modeling of absorption (as fraction of attenuation) and perfusion @jkosciessa. Be advised that this renders the heating estimates more liberal than in the past. See the documentation. Also see some comments in the issue.
- [**Feature**] axisymmetry simulations with 2D input grids (doc) @jkosciessa
- [**Feature**] simulations with computational phantoms (doc) @jkosciessa
- [**Documentation**] Initial draft of documentation (e.g., parameters) (see links above and here) @jkosciessa
- Many smaller bug fixes and QOL improvements (e.g., ability to disable .mat saving, more debugging outputs)
- DOC: minor typo in README.md by @mekman in #58
- Run simnibs segmentation code also supports SLURM by @MaCuinea in #62
- Prepare acoustic_profiling for merge by @MaCuinea in #69
- Semi-automated acoustic profiling by @MaCuinea in #71
- Added note to default config regarding ld_library_path in slurm by @KTZ228 in #70
- Merge pull request #70 from KTZ228/development by @MaCuinea in #77
- Fixed tutorial by @MaCuinea in #73
- Feature: consecutive simulations by @KTZ228 in #79
- Documentation, pCT, HCP updates & more, see above by @jkosciessa in #59
- The last one summarizes more than half a year of updates to code and documentation.

---

## v0.3.0

- PseudoCT feature by @eleonoracarpino and @jkosciessa
- Clean-up repository by @jkosciessa
- Incorrect allocation of attenuation values bug fix by @jkosciessa
- Improve correct labeling of tissue by @jkosciessa
- Improvement of 3D plots by @jkosciessa
- Check simulation stability and adjust time step accordingly by @MaCuinea and @jkosciessa
- Enable tissue-specific temp_0 specification by @jkosciessa
- Fix heating plot bug by @jkosciessa

---

## v0.2.0

- DOC: add GNU license by @mekman in #17
- Documentation changes and bug fixes by @KTZ228 in #14
- Disabled masking of neural tissue in skin in smooth_and_crop_layered by @KTZ228 in #19
- Adjusted create_Group_MNI_plots by @KTZ228 in #22
- Documentation changes by @KTZ228 in #24
- Minor changes to comments for heating setup by @jkosciessa in #29
- Figure saving updates, some bug fixes by @jkosciessa in #8
- Kenneth master by @KTZ228 in #33
- Group plots fixes and improvements by @KTZ228 in #26
- FIX: revert to working version by @mekman in #35
- DOC: reworked installation instructions by @mekman in #36

---

## v0.1.0

Initial in-house development version by @achetverikov

- Dedicated segmentation path option, minor bug fixes, updated simnibs call by @jkosciessa in #2
- fix small errors by @jkosciessa in #3
- adjusted skull values by @eleonoracarpino in #4
- Bug fixes and documentation changes by @KTZ228 in #6
- DOC: simplify readme.md wording regarding requirements by @mekman in #7
