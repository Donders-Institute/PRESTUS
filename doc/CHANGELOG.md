# Changelog

Notable changes to this project are documented here. 

---

## v0.7-pre (development)

**Fixes**
- Sequential chains: the follow-up's start state is the base run's T1w end map (heating, CEM43, CEM43 ISO), resampled into its own grid. The NIfTI stage is enabled for chained runs, and dispatch errors (`sequential_pipeline:missingThermalEnd`) when a map is missing, instead of starting from the water baseline.
- `convert_final_to_MNI_simnibs`: `subject2mni` runs one volume at a time in the foreground. A failed conversion now raises `PRESTUS:MNIConversionFailed` and fails the job. The previous background path discarded output and reported success when volumes were missing.
- Startup tool check (`check_tool_availability`): `prestus_pipeline_start` (before HPC submission) and `prestus_pipeline` (inside the job) check the external tools the enabled stages need: JVM, `sbatch`/`qsub`, the SimNIBS binaries (`charm`, `subject2mni`, `mni2subject_coords`) and `m2m/toMNI`, the k-Wave C++ binary for the selected `code_type`, and a GPU. A missing tool stops the run with `prestus:missingTools`, listing every missing tool. Another process holding more than half of a GPU's memory fails the check; more than 5% raises `prestus:gpuInUse`.
- New `hpc.exclude` (SLURM): nodes passed to `sbatch --exclude`, e.g. a node whose GPU is held by a stale process.
- HPC jobs request a GPU from the solvers that run and their compute backends (k-Wave `simulation.code_type`, `jwave.use_gpu`, `babelbrain.computing_backend`), not from `simulation.code_type` alone.
- ⚠️ **Fixed:** the rubber-wrap skull fill (`headmodel.skull_fill_method: rubberwrap`) ignored the brain and skin layers, so skull gaps where brain meets skin stayed open; rerun the grid stage of head runs
- ⚠️ **Fixed:** with separate cortical and trabecular skull layers, filled skull holes were labelled trabecular instead of cortical; rerun the grid stage of runs with both layers

Adds BabelBrain and jwave as selectable acoustic backends, a transducer calibration library with per-depth refinement, partial-FOV acoustic simulations, study manifests with group analysis, and the PRESTUS Nexus as the new configuration and run interface. Placement gains MNI and TUS_entry modes, a placement report and reusable placement runs; pseudo-CT gains new skull mappings and CT calibrations. Several result-changing fixes affect MI, simulation duration, axisymmetric sources, coherent multi-transducer runs, RAS coordinate conversion and Localite placement; check the ⚠️ lines below.

> ⚠️ **Migration note:** Several defaults have changed, and some outputs have moved.
>
> **Defaults:**
> - `grid.mode`: `transducer_axis` → `focal_axis`, which aligns the grid with the mean beam axis of all transducers (opposed beams averaged as axes). Single-transducer runs are unchanged; multi-transducer runs get a different grid, so rerun them rather than reusing caches
> - `grid.default_dims` (free-water grid): `[144, 144, 400]` → `'auto'` (sized to the transducer)
> - `transducer.target_isppa_wcm2`: `[]` → `30`, and `modules.run_water_baseline`: `0` → `auto`. Results are scaled to 30 W/cm² free-water ISPPA unless you set `target_isppa_wcm2: []`; a target that cannot be met stops the run
> - In coherent multi-transducer runs, `transducer(k).target_isppa_wcm2` now sets the target of transducer k alone (only transducer 1's target was read before)
> - `io.save_source_matrices` / `save_acoustic_matrices` / `save_thermal_matrices` default to `[]` (inherit `io.save_matrices`); default saving is unchanged
> - `placement.reuse: auto`: a recorded heuristic, PlanTUS or TUS_entry placement with matching settings is reused, also under `io.overwrite_files: always`
> - `io.overwrite_files: always` now re-runs the output stages of chained pipelines (async, multi-ISPPA, uncertainty); new `force` also re-runs their intermediate stages
>
> **Config layering:** the toolbox `config_default.yaml` is always the base layer and a project copy is merged on top. Refresh an old project copy with `prestus_config_upgrade`. Unknown study-config keys now warn.
>
> **Output locations:**
> - Heuristic, PlanTUS and TUS_entry runs write to `<path.placement>/sub-XXX/<mode>/<target_name>/` (new `path.placement`, default `<path.sim>/placement`). Earlier PlanTUS results in `<path.seg>/PlanTUS/` and heuristic `tpars_*.csv` files are not found there: re-run them or copy them over
> - Pseudo-CTs are written to `<path.pct>/sub-NNN`, or `<parent of path.seg>/pct/sub-NNN` when `path.pct` is empty, never into the SimNIBS folder. Pseudo-CTs in the old `m2m_sub-NNN` location are still used, with a warning
>
> **Deprecated / removed:**
> - `path.localite` (use `path.localite_raw` / `path.localite_export`) and `placement.localite.markertype` (detected from the files) are deprecated
> - `pct.skull_mapping: carpino` is deprecated; use `kosciessa`, or `custom` with `custom_slope: -2194`, `custom_intercept: 2236`
> - `prestus_group_start` and `create_group_MNI_plots` are deprecated; `prestus_group_start` runs `prestus_group_report_start`
> - `launch_web_alignment` and the browser alignment viewer are removed; use the Nexus placement viewer

#### Core / Config
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

#### Transducer library & calibration
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
- `IS_PCD15473_01001_OPM` and `IS_PCD15473_01003_OPM` are excluded until the OPM setup is documented; use `IS_PCD15473_01001` / `_01003` with `IGT_32_ch_comb_10_ch`
- Equipment and shipped library entries are grouped by device model: `config/equipment/model-<MODEL>/tran-<SERIAL>_transducer.yaml`, driving systems in `config/equipment/ds/`, `config/calibration/model-<MODEL>/tran-<TRAN>_ds-<DS>_calib.yaml`. Units take geometry and steering range from their model entry. Flat folders (`<SERIAL>.yaml`, `<TRAN>_<DS>.yaml`) are still read, so your own equipment folders and `transducer.library_path` libraries work unchanged
- Your own transducer libraries: equipment entries (`*_transducer.yaml`, `*_driver.yaml`), calibration entries and `profiles/` in a folder of your own add transducers or replace shipped ones. PRESTUS checks `transducer.library_path` (study), new `path.transducer_library` (project) and the `PRESTUS_TRANSDUCER_LIBRARY` environment variable (machine), then the shipped library
- `calibration_gui` saves new transducers into your own library and no longer writes into `config/equipment` (option `AllowShipped` for maintainers)
- Calibration entries record the transducer geometry they were fitted for (`meta.geometry_fingerprint`); using one with a different geometry warns (`load_transducer_from_library:geometryMismatch`)
- Drive limits: `transducer.max_elem_p_amp`, or a combo's `drive_limit` in the equipment entry, warns when the library drive for the target ISPPA exceeds what the driving system delivers (`resolve_transducer_from_serial:driveLimit`); estimates shipped for the IGT 128-channel driver
- ⚠️ **Fixed:** calibration files (`*-F*mm-I*wpercm2.yaml`) recorded the default `target_isppa_wcm2` instead of their calibrated intensity, so a config built from one was rescaled to 30 W/cm². Set `target_isppa_wcm2` in configs built from older calibration files

- A library transducer (`transducer.serial`) with an empty `target_isppa_wcm2` and no inline drive now stops with a message instead of failing during the simulation
#### Transducers & sources
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
- ⚠️ **Fixed:** in coherent multi-transducer runs, transducers 2..N were built from transducer 1's drive, phases, geometry and frequency; re-run coherent multi-transducer acoustics (single-transducer, async, sequential and BabelBrain runs are unaffected)
- Fixed: flat matrix arrays (`is_curved: false`) failed; `extract_from_file.file_path` accepts absolute paths outside the MATLAB path
- ⚠️ **Fixed:** without `focal_distance_ep` or `focal_distance_bowl`, the focal distance derived from voxel `trans_pos` / `focus_pos` rounded the T1 voxel size to whole mm (0.8 or 0.9 mm read as 1 mm, focal distance 11–25 % too long); rerun such runs on non-1 mm or anisotropic T1s
- ⚠️ **Fixed:** curved matrix arrays without `dist_geom_ep_mm` put the exit plane at the bowl apex; it is now derived from the aperture as for annular arrays. Rerun curved matrix runs that left `dist_geom_ep_mm` empty
- ⚠️ **Fixed:** `transducer.trans_pos_mm` / `focus_pos_mm` mapped 0 mm to the invalid voxel 0; they are now measured from the centre of voxel 1 (`round(mm / resolution_mm) + 1`). Subtract one voxel (`resolution_mm`) from configs that compensated for the old mapping

#### Placement
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
- PlanTUS 2.0 support (Workbench-free); PlanTUS 1.x still works. New `placement.plantus.skin_smoothing_iterations` (default 2, as in PlanTUS 2.0)
- ⚠️ **Changed:** with PlanTUS 2.0, the skin surface is smoothed and the beam aims at the target centre; placements can differ from PlanTUS 1.x runs. Set `skin_smoothing_iterations: 0` to keep the raw surface.
- New `placement.plantus.vertex`: places a chosen skin vertex on the maps of an earlier run without recomputing them
- New `placement.heuristic.candidate`: places a chosen candidate of the heuristic's table instead of applying its criteria
- Fixed: heuristic placement without `transducer.focal_distance_ep` / `focal_distance_bowl` never finished; it now stops with a message
- ⚠️ **Fixed:** the heuristic placement's Localite mm export (`loc_*_mm`) was one voxel high along each axis; re-export heuristic positions before entering them in Localite
- Fixed: Localite placement with empty transducer geometry stopped with a dimension error instead of using `placement.localite.reference_distance_mm`
- A transducer without `dist_geom_ep_mm` gets its exit plane at the bowl rim (outer diameter of the largest element, or `outer_diameter_mm` for matrix arrays), with a notice; the Nexus draws the same plane and no bowl for flat matrix arrays

#### Neuronavigation (Localite)
- Batch ingestion (`neuronav_ingest_batch`), session inventory, target-label maps (`neuronav.target_labels`) and target-name suggestions from planned targets
- `neuronav_resolve_target` resolves targets over four confidence tiers, including a study-supplied lookup (`neuronav_read_target_lookup`)
- Multi-coil ingestion via `neuronav.coil_map`
- New `path.localite_*` and `placement.localite.*` keys in `config_default.yaml`; `placement.localite.session_folder_pattern` is now applied
- `placement.localite.target` / `.position` take one entry per transducer for bilateral setups
- ⚠️ **Fixed:** series QC statistics `rotation_dev_rad` and `position_dev_mm` were wrong (the latter in voxels); sidecar JSONs written before carry the old values
- ⚠️ **Fixed:** `neuronav_convert_native_to_MNI` returned transducer MNI positions one MNI voxel (1 mm) off along each axis; recompute MNI transducer positions taken from it
- Fixed: GUMMarkers files with `InstrumentMarker` elements failed to parse; T1 lookup in the conversion functions; `neuronav_convert_MNI_to_native` / `neuronav_convert_native_to_MNI` read an old path key

#### Head model & pseudo-CT
- Automated defacing via pydeface, optionally in the segmentation pipeline
- ⚠️ **Fixed:** grid, medium and acoustic caches, placement runs and pseudo-CTs of a layered head were reused after the segmentation was re-made (charm re-run). They now record the segmentation's checksum and are checked on reuse; a defaced copy is made again. Re-run simulations of subjects re-segmented after their caches were written.
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

#### Acoustic / Thermal
- Each tissue sets heat deposition as `medium_properties.*.absorption_fraction` (default) or `absorption_coeff`; a config layer that sets one empties the other, and setting both in one layer errors
- Logs, run summary and HTML methods give each tissue's effective absorption coefficient and fraction; property maps include `absorption_np` and `absorption_fraction`
- Partial acoustic FOV along the beam axis (`grid.acoustic_fov_length_mm`, `grid.acoustic_fov_diameter_mm`, `'auto'` from focal geometry) for deep targets; the thermal grid can match it (`thermal.match_acoustic_fov`)
- Free-water grid `'auto'` sizing (`grid.water_margin_mm`, `grid.z_beyond_focus_mm`); matrix transducers placed by `trans_pos` / `focus_pos` keep `[144, 144, 400]`
- `grid.water_baseline_grid: 'beam'` (default) runs the k-Wave free-water baseline on a box around the source and beam: the same free-field result at a fraction of the cost; `'full'` uses the whole grid
- Time axis and steady state configurable: `grid.steady_state_transits`, `grid.steady_state_num_periods`, `grid.source_ppp`, `grid.sim_time_us`; PML and band-limiting options `grid.pml_size_range`, `grid.pml_alpha`, `grid.bli_tolerance`, `grid.bli_upsampling_rate`
- Lower memory use for CW sources on `matlab_cpu` / `matlab_gpu`, also on axisymmetric grids
- Axial-profile plot with tissue bands and the free-water O'Neil and k-Wave profiles; CSV gains `freewater_oneil_Wcm2` / `freewater_kwave_Wcm2`
- Time reversal (`time_reversal.enabled`, *draft*): one reverse simulation per target gives refocusing phases, aberration metrics and a ranking of candidate poses; `placement_search.jitter_mm` takes several distances, and `jitter_mode: scalp` slides the pose over the scalp and re-seats it; candidate poses keep `placement.skin_gap_mm` (heuristic, MNI, PlanTUS, TUS_entry) or the planned pose's gap
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
- ⚠️ **Fixed:** axisymmetric phantom runs that reloaded the acoustic cache crashed on a half-grid medium; the reload now expands grid, positions and medium like a fresh simulation and honours `grid.axisymmetric_output`
- ⚠️ **Fixed:** axisymmetric NIfTI outputs were cut at the grid centre instead of the beam axis when the transducer sat off centre. Re-export NIfTIs of off-centre runs
- ⚠️ **Fixed:** phantom NIfTI outputs carried no affine and kept the simulation's voxel order; they now carry the phantom's affine in its voxel order (transposed or axial-first 2D phantoms, plane of a 3D phantom). Re-export phantom NIfTIs
- ⚠️ **Fixed:** axisymmetric 2D output with the beam axis right of centre failed; beyond the simulated radius the plane repeats the edge value, as the 3D expansion does
- Axisymmetric runs keep the bilateral transducer positions when reusing a grid cache, and the thermal stage uses the expanded positions on a separate thermal grid
- Phantom voxels must be isotropic; an axisymmetric run on a 3D phantom cuts the plane at a 3D `transducer.trans_pos` (warns when the centre is ambiguous)
- Focus positions with `focus_pos(1) > focus_pos(2)` are no longer turned into column vectors
- Grid, medium and acoustic caches of phantom runs record the phantom file; a reused acoustic cache is checked against the run's settings (`io.cache_signature_check`)
- Removed `thermal.sensor_xy_halfsize` and `thermal.record_t_at_every_step`: neither had an effect on results, and configs that set them load with a notice

#### Multi-transducer
- Incoherent summation of independent-element arrays; exit-plane mask covers all transducers
- Async mode (*experimental*): shared grid and medium stage before the per-transducer jobs, partial FOV supported, per-transducer ISPPA table in the report; HPC jobs chain via job dependencies, and `[job_id, jobs] = prestus_pipeline_start(...)` returns the job IDs per stage
- ⚠️ **Fixed:** `io.overwrite_files: 'always'` was ignored by chained pipelines (async, multi-ISPPA, uncertainty), so stale outputs could block the joint job and its report; rerun chained runs written into folders with earlier outputs

#### Backends
- **BabelBrain** as acoustic and thermal backend (`acoustic_backend: babelbrain`, `thermal_backend: babelbrain`), with a CT path (`babelbrain.ct_mapping`), internal air regions (`babelbrain.air_regions`) and free-water FDTD by default (`babelbrain.use_rayleigh_for_water: false`). One transducer per run, beam along the grid z axis; BabelBrain thermal uses its own tissue thermal properties. Defaults: `babelbrain.ppw: 12`, `babelbrain.q_correction` 1.0 on the non-CT path (3.0 with CT), `babelbrain.tight_narrow_beam` off, `babelbrain.profile_calibration` off
- BabelBrain 0.8.8 (`external/BabelBrain`). Update the BabelBrain environment: `pip install babelviscofdtd==1.2.7 psutil`
- BabelBrain simulates the PRESTUS transducer's ring geometry, on which its calibration was fitted, unless `babelbrain.tx_system` names a BabelBrain preset (resolved from BabelBrain's registry, incl. `R16081`, `TheraPhas`; `REMOPD` still accepted)
- BabelBrain `Single` and `BSonix` presets take their aperture and focal length from the PRESTUS transducer; with `thermal_backend: babelbrain`, multi-ISPPA thermal jobs run one after another, as they share one BabelBrain thermal file
- **jwave** as acoustic backend (`acoustic_backend: jwave`), Helmholtz and time-domain solvers; `jwave.helmholtz_operator: conservative` (default) corrects transmission through skull, `jwave` restores the built-in operator. kWaveArray sources supported; no axisymmetric grids
- **k-Wave elastic** (`acoustic_backend: kwave_elastic`, *draft*, not validated): shear waves in the skull with properties from the new `elastic:` section; MATLAB CPU/GPU, no axisymmetric grids
- Backend settings are validated at config load; the HTML report has backend sections and backend-specific methods text

- Fixed: `babelbrain.z_beyond_focus_mm` was ignored (40 mm was always used)
#### Advanced modes (uncertainty, multi-ISPPA, sequential)
- Uncertainty composes with multi-ISPPA sweeps and async multi-transducer runs, with one report per sweep point; uncertainty is the outer mode on every platform
- Uncertainty mode runs five variants: default, Temp-H / Temp-L (the ITRUSST tissue-property bounds that maximise / minimise brain heating) and Trans-H / Trans-L (the acoustic bounds that maximise / minimise brain intensity; acoustic only). See [doc_uncertainty.md](doc_uncertainty.md)
- `simulation.uncertainty_temp_high`, `_temp_low`, `_trans_high` and `_trans_low` choose the variants' medium files (default: the shipped `config/uncertainty/` files); `options.<variant>_config` takes precedence
- ⚠️ The liberal / conservative variants are replaced: new medium bounds and output affixes (`_desc-temphigh`, `_desc-templow`, `_desc-transhigh`, `_desc-translow`), and `options.liberal_config`, `options.conservative_config` and `options.affixes.liberal` / `.conservative` now stop with `uncertainty_pipeline:legacyVariant`. Re-run uncertainty runs
- The uncertainty report ranges acoustic metrics over the Trans pair and thermal metrics over the Temp pair, and lists every variant's values
- ⚠️ **Fixed:** uncertainty safety cards and range-table rows were coloured by the conservative variant; they now follow the highest value of the default and the metric's pair. Regenerate uncertainty reports
- ⚠️ **Fixed:** uncertainty variants shared one medium cache file and overwrote each other with `io.save_medium_cache: 1`; rerun uncertainty runs that used it
- Fixed: the uncertainty report is named after the default variant's affix, so runs with different `io.output_affix` no longer overwrite each other's report; HPC resume finds existing reports and variant CSVs
- Fixed: the uncertainty CSV is saved without an extra `tabular/` folder; variant and report jobs no longer run their own free-water baseline
- ⚠️ With `pct.enabled: 1`, the uncertainty bound variants use a homogeneous skull layer on the pseudo-CT run's head geometry; only the default maps the skull from the pseudo-CT. Re-run pseudo-CT uncertainty runs
- The uncertainty report of a pseudo-CT run has a Pseudo-CT Skull vs. Bounds section: the default run's skull property histograms against the bound variants' skull values
- Fixed: `multi_isppa_pipeline` warned about missing thermal output for every intensity (wrong CSV name)

#### Study setup & group analysis
- Study manifest (`study.manifest`, `study.name`): one row per run with subject, condition and any config key as a column; `prestus_study_start` runs or submits every row and builds a study report
- Targeting accuracy (`modules.run_targeting_analysis`, `analysis.target_mni_mm`, `analysis.target_roi_mni`): focal centroid and peak against the planned focus, an MNI target and an ROI
- Group MNI maps (`create_group_maps`): per-condition focal-volume overlap, mean intensity and heating
- Group report per run (subject × condition), with condition comparison, targeting, group maps and outliers sections
- ⚠️ **Fixed:** `analysis.mni_warp_method: native` placed MNI outputs one subject voxel off along each axis; rerun MNI exports made with it

#### Outputs & reporting
- Reports redesigned (`modules.report_style`: `instrument` default, `classic`, `compact`, `paper`); `preview_report_styles` re-renders an output folder in every skin
- Manuscript-ready Methods section with copy buttons and stated tissue properties
- Report shows the calibration provenance of the transducer (library entry, k-Wave-verified or analytical only)
- ⚠️ **Changed:** focal-distance columns name their reference: `real_focal_distance_mm` → `real_focal_distance_bowl_mm` (bowl apex), next to `real_focal_distance_ep_mm` (exit plane). Reports and the Nexus read older CSVs; update your own scripts that read `real_focal_distance_mm`
- Safety dashboards flag peak pressure against a 2 MPa limit and show `MI` and `MItc` side by side
- Run statistics in `<path.sim>/run_stats.jsonl`, used for resource estimates in the Nexus and `python -m prestus_nexus.estimate`
- Fixed: Methods references for Martin et al. (2024), Pichardo et al. (2017) and Aubry et al. (2025) had wrong details; regenerate reports to update them

#### HPC
- Post-job resource report via SLURM `sacct`, with peak memory (including BabelBrain) and a CO₂e estimate
- `hpc.posthoc_water_timelimit` / `hpc.posthoc_water_memorylimit`; `hpc.depend_type` for dependent jobs
- Fixed: every Snellius SLURM submission failed (wall-time check)
- Fixed: qsub scripts wrote `walltime` as character codes and requested GPU nodes for CPU-only jobs; SLURM job IDs are read with `sbatch --parsable`
- Fixed: `io.overwrite_simnibs` never passed `--forcerun` to charm; `simnibs_t1_to_mni` respects `overwrite_files`
- ⚠️ **Fixed:** under `io.overwrite_files: never` or `ask`, a recomputed acoustic or thermal stage kept existing NIfTIs but rewrote the CSV, mixing outputs of two runs. Re-run output folders that were simulated more than once under `never`

- Fixed: PBS jobs loaded MATLAB R2023b, which does not support the CUDA 11.2 of PBS nodes; they now load R2022b
#### Nexus & GUIs
- Fixed: the Nexus suggested `transducer.target_isppa_wcm2` for the old `set_intensity_w_per_cm2`, which is the drive's calibration ISPPA (`transducer.tx_calibration.isppa_wcm2`)
- **PRESTUS Nexus** (`nexus/`, *experimental*): standalone PyQt6 app to build, check and run study configs, with tabs per workflow step, a transducer panel backed by the library, a placement viewer with beam checks, Calibration and Localite pages, a study manifest editor, a Results tab and cluster submission (also over SSH). MATLAB entry point: `prestus_nexus`. See [doc_nexus.md](doc_nexus.md)
- The Nexus Pseudo-CT sub-tab plots the in-skull property distributions against the band of the uncertainty bound variants' skull values, with the fraction of skull voxels outside it; **Predict skull properties** predicts them before a run
- `prestus_gui` is deprecated in favour of the PRESTUS Nexus; it still works and warns on start. `neuronav_gui` and `calibration_gui` remain supported.
- `neuronav_gui`: stand-alone Localite ingest GUI with series QC, batch window and MNI view of batch targets
- `calibration_gui` and the Nexus Calibration page expose every calibration fit option and save how an entry was fitted (`meta.fit`)
- Nexus Calibration page: the library entry's fit quality, drive and phases are shown before MATLAB starts; the Transducer tab has a calibration coverage card
- Nexus: tissue-property card is read-only until unlocked; with uncertainty on, a card lists and plots the three variants' properties
- Nexus: a Start tab sets up a first simulation in four steps (T1, transducer, MNI target, run); pages without their inputs say what is missing; the Issues list opens from the header (View ▸ Issues…, Ctrl+Shift+I)
- Nexus: with `placement.mode: plantus`, the Placement card runs only the PlanTUS placement, choosing the position in Workbench or automatically, and opens the result in the Viewer
- Fixed (`prestus_gui`): NaN placeholders wrote `0` into `grid.source_ppw`; 2D phantom coordinates; thermal timing state on config load

#### Telemetry
- `PRESTUS_TELEMETRY` (`off` / `on`) sets the telemetry decision for the current process
- Fixed: `matlab -batch` and SLURM/PBS sessions no longer prompt for consent; they print the notice and record nothing

#### Documentation
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
> - `parameters.submit_medium` renamed to `parameters.platform`
> - Segmentation smoothing now specified in FWHM mm (parameter name changed)
> - Transducer field suffixes renamed consistently — check transducer config fields
> - Further parameter renames throughout — consult `config_default.yaml` for current names
>
> **Output paths:**
> - `parameters.io` fields renamed: `*_dir` → `dir_*`, `*_filename` → `filename_*`
> - Preprocessing cache files follow a consistent `res_N_descriptor` naming scheme
> - `space-` prefix removed from all NIfTI output filenames
> - `output_affix` now appended at the end of all output filenames
> - HPC log directory renamed: `hpc_log/` → `log_hpc/`
> - Outputs reorganised into BIDS-inspired typed subdirectories (`acoustic/`, `thermal/`, `cache/`, `debug/`, etc.)
>
> **Directory and config naming:**
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
- [**feature**] PRESTUS GUI — tabbed parameter editor and simulation launcher (see [doc_gui.md](doc_gui.md))
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
