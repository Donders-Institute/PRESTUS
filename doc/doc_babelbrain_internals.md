# BabelBrain backend: integration internals

Notes for developers and for debugging the BabelBrain bridge. For using the backend see [BabelBrain backend](doc_babelbrain.md).

> This section describes how the BabelBrain integration is wired inside PRESTUS. It is intended for developers working on the backend or extending GPU platform support.

## Call chain

```
prestus_pipeline.m
  └─ acoustic_stage.m            ← dispatches on parameters.acoustic_backend (run_acoustic_backend)
       └─ babelbrain_wrapper.m   ← MATLAB orchestrator
            ├─ babelbrain_transducer_map.m   (preset resolution)
            ├─ write_babelbrain_input.m       (NIfTI label mask)
            ├─ export_babelbrain_weights.m    (calibration HDF5)
            └─ system() → run_babelbrain_case.py   ← Python bridge
                 └─ BabelBrain RUN_SIM().RunCases()
```

`babelbrain_wrapper.m` never imports Python directly. It constructs a command string and calls `system()`, which invokes `run_babelbrain_case.py` as a subprocess. All BabelBrain Python dependencies are therefore isolated in the conda environment pointed to by `babelbrain.python_env`; the MATLAB process only needs to call `system()`.

`run_babelbrain_case.py` adds the BabelBrain source tree (`external/BabelBrain`) to `sys.path` at runtime, imports `TranscranialModeling.BabelIntegrationBASE` dynamically, and calls `RUN_SIM().RunCases()`. The output `*DataForSim.h5` is written to `<dir_cache>/<subject_id>/`; PRESTUS reads it back with `read_babelbrain_output.m`.

## GPU backend selection and platform support

`computing_backend` maps directly to BabelBrain's `ComputingBackend` integer:

| Value | Backend | Supported hardware | Typical environment |
|---|---|---|---|
| `1` | CUDA | NVIDIA GPUs | Linux HPC / workstation with CUDA toolkit |
| `2` | OpenCL | AMD GPUs, Intel GPUs, some NVIDIA | Linux, any OpenCL-capable GPU |
| `3` | Metal | Apple Silicon (M1/M2/M3/M4) | macOS on ARM |
| `4` | CPU | Any | Fallback; much slower |

**Mac / Apple Silicon.** Set `computing_backend: 3` (Metal). This is the recommended path on M-series Macs; `device_name` can be left at `'auto'`. k-Wave's CUDA binary is unavailable on Apple Silicon, so BabelBrain + Metal is the only GPU-accelerated acoustic solver available on that platform.

**Linux / NVIDIA.** Set `computing_backend: 1`. With `device_name: auto` the Python bridge auto-detects the first CUDA GPU via `nvidia-smi`, so you normally do not need to set `device_name` explicitly on SLURM clusters. Ensure the SLURM job requests a GPU (`--gres=gpu:1`).

**Linux / AMD or other OpenCL GPU.** Set `computing_backend: 2`. The bridge uses `pyopencl` to auto-detect the first available GPU device. Install `pyopencl` in the BabelBrain conda environment if not already present.

**macOS Intel.** Metal is available but performance is limited on Intel integrated graphics. Use `computing_backend: 3` with `environment_mac_x64-314.yml`; for CPU-only runs use `computing_backend: 4`.

**Windows.** PRESTUS has not been tested on Windows. BabelBrain ships an `environment_win-314.yml`; the MATLAB `system()` call uses Unix-style path separators which would need adjustment.

## Device name auto-detection

`run_babelbrain_case.py::resolve_device_name()` handles `device_name: auto` by detecting the appropriate device for the selected backend:

- **CUDA (backend 1):** calls `nvidia-smi --query-gpu=name` and uses the first result.
- **OpenCL (backend 2):** iterates `pyopencl.get_platforms()` and picks the first GPU device.
- **Metal (backend 3):** resolves to `'Apple'`; no detection needed.

Set `device_name` explicitly to target a specific GPU when multiple are present.

## Transducer placement, ZSteering, and the drive phases

BabelBrain builds its simulation domain with the transducer apex at `z_mm ≈ 0` and the **FocalSpot (target) voxel at `z_mm = ` preset `FocalLength`** (read from `Babel_{tx_system}/default.yaml`) — i.e. it assumes the transducer sits at its *natural* focal distance from the target. PRESTUS's geometry instead places the apex at `trans_pos` (bowl rear centre), i.e. `|focus_pos − trans_pos|` from the target. That distance equals `focal_distance_bowl` only when the target sits at the nominal focal depth; planned/pFOV placements generally differ (e.g. 71 mm vs 92 mm for R15473 at 85 mm EP with a 64 mm-deep target). On the way back, `resample_babelbrain_vol.m` anchors the `z_mm = FocalLength` plane onto PRESTUS's `focus_pos` (`prestus_dim3 = (z_mm − focal_length_mm)/res + focus_pos`), so any difference between the two distances **rigidly translates the pasted field by `|focus_pos − trans_pos| − FocalLength`** along the beam axis (−7.7 mm for R15473 with the target at 85 mm EP; −28.75 mm with the target at 64 mm EP). Electronic steering cannot compensate for that translation: a translated field of a differently-steered transducer is not the true field — the near-field lobe structure shifts with the translation while the drive physics differ.

The wrapper therefore corrects the geometry **mechanically**:

```
TxMechanicalAdjustmentZ = preset_FocalLength − |focus_pos − trans_pos| · res_mm   (+ any user standoff)
```

It is computed by `babelbrain_placement_correction.m` (shared with `export_babelbrain_inputs`). Only when `trans_pos` is unavailable does it fall back to `preset_FocalLength − focal_distance_bowl`, which assumes the target sits at the nominal focal depth. Using `focal_distance_bowl` for a target at a different depth places BabelBrain's transducer too far from (or too close to) the head, shifting the whole field — skull incidence, near-field lobes and focus — relative to the anatomy.

This shifts the Tx *inside BabelBrain's domain* (`BabelIntegrationANNULAR_ARRAY`, `CalculateRayleighFieldsForward`) so the apex sits at its true distance from the target; the focus-anchored resample is then also apex-aligned (distance-preserving). BabelBrain's internal Rayleigh phase reference aims at `center_z = ZDim[FocalSpot] + ZSteering + TxMechanicalAdjustmentZ`, so the mechanical adjustment shifts Tx and aim together — with `ZSteering = 0` the aim is always the **natural focus relative to the apex**.

On top of the corrected geometry, the drive is set in one of two modes:

- **Explicit element phases** (`tr.annular.elem_phase_rad` present — calibrated/analytical configs): `ZSteering = 0` and the weights file carries the **full drive** as a delta from the natural-focus reference (where BabelBrain's internal ring phases are flat, all rings being equidistant from the centre of curvature). BabelBrain then reproduces *exactly* the per-ring phases that k-Wave/jwave/the Rayleigh model use — solver-benchmark parity by construction. `focal_calibration` runs **locate-only** in this mode: it reports the predicted focus of the fixed drive but never replaces it.
- **No explicit phases** (geometric fallback): BabelBrain steers electronically, with the aim measured physically from the apex:

  ```
  ZSteering = aim_bowl − preset_FocalLength
  ```

  where `aim_bowl` is the calibrated steer `D*` when `focal_calibration: true` (default; `calibrate_focal_steering.m` root-finds the geometric steer whose achieved `flhm_center`/`peak` lands on `focal_distance_ep`, correcting the diffraction focal shift), or `focal_distance_bowl` (uncalibrated).

Pin `babelbrain.z_steering_m` to override the auto-derivation entirely (also bypasses calibration; the weight-delta reference follows the pinned aim `FocalLength + ZSteering`).

**Free-water cross-check.** With explicit phases, the BabelBrain axial profile must overlay the analytic Rayleigh profile of the same phases in *both* the near-field lobes and the focal lobe. A rigid axial offset of ≈ `FocalLength − |focus_pos − trans_pos|` (near-field lobes shifted by the same amount as the peak) indicates the mechanical placement correction is missing; a *shape* mismatch of the near-field lobes indicates the effective drive phases differ (e.g. weights delta referenced to the wrong aim).

**Bowl rim vs. the source plane.** BabelBrain injects the Rayleigh field on the first domain slice on the transducer side and requires every source element to lie behind it; if the bowl rim protrudes into the domain it snaps the transducer back in whole voxels (`Warning: Need to apply correction to reposition Tx for …` in the BabelBrain log, and a non-zero `source_residual_mm` in the `_prestus_anchor.json` sidecar). This happens whenever any part of the bowl lies inside the grid: the rim sits one bowl depth (5.4 mm for R15473) on the target side of the apex, so the bowl spans `d_apex … d_apex + bowl_depth` from the transducer-side face (`d_apex` = apex distance into the grid, negative if the apex is outside). Examples: a pFOV crop that starts at `trans_pos` (source pushed 5.5 mm too far back), or the bone-slab phantom with the apex 7.5 mm inside the grid (13.5 mm, which shifted the BabelBrain focus ~14 mm proximally). `write_babelbrain_input.m` therefore trims `ceil((d_apex + bowl_depth)/dx) + 1` transducer-side slices from the BabelBrain mask (it is water inside the bowl concavity; a warning is raised if it contains tissue). The trimmed slab is zero-filled when the field is pasted back.

**Realised bone absorption (`q_correction`).** On the non-CT path the bone attenuation BabelBrain actually realises depends on the Q-factor correction the runner applies to the bone slots, and scales ~1/Q. A flat bone-slab phantom (6 mm cortical, 78.5 Np/m at 300 kHz, PPW 6; realised attenuation from paired runs with bone `alpha_coeff` as specified and set to 0) gave realised/specified attenuation 0.91 at Q = 1 (k-Wave 0.91), 0.29 at 3, 0.14 at 6, 0.087 at 10 and 0.063 at 14. BabelBrain's CT-mode `QCorrection` default of 3.0 is its calibration for its own CT-derived bone attenuation; applied to PRESTUS-specified `alpha_coeff` it removed two thirds of the bone attenuation (sub-602 pFOV transcranial Isppa 10.6 vs 7.6 W/cm² for k-Wave). The non-CT runner therefore uses Q = 1 by default. `babelbrain.q_correction` > 0 overrides it (and BabelBrain's CT-path value, via `--q_correction` → `BABEL_PYTEST_QFACTOR`, read in `BabelFTD_Simulations_BASE.__init__`). MATLAB `setenv` cannot be used instead: `run_command` launches Python via Java's ProcessBuilder, whose environment is the JVM start-up snapshot.

## Output alignment: tissue-marker anchoring and the ZDomainStart fix

**The problem.** BabelBrain's own `MaterialMap` — its internal record of the tissue mask it actually simulated on — is byte-identical to the mask PRESTUS wrote (`write_babelbrain_input.m`); BabelBrain never alters tissue. So if the returned field and material map are pasted back onto the PRESTUS grid at the wrong position, *every* downstream consumer that assumes k-Wave/jwave's "one shared grid, no internal resampling" invariant — `acoustic_analysis.m`'s per-tissue metrics (`Isppa_brain`/`Isppa_skull`), `thermal_analysis.m`'s per-tissue thermal/CEM43 metrics, the axial-profile and intensity-overlay plots, and any exported NIfTI compared against `medium_masks` in an external viewer — is silently wrong at tissue boundaries, most visibly at the skull-brain interface where BabelBrain's elastic-skull solver produces a real, physically-expected reflection hotspot.

Two independent effects combined to create a multi-millimetre misalignment (observed: ~5.5 mm on a benchmark pFOV run):

1. **`ZDomainStart`.** `CalculateRayleighFieldsForward` places the transducer at `Tx_z += TxMechanicalAdjustmentZ + ZDomainStart`, where `ZDomainStart = CalculateDomainZReference()` is a skin-surface heuristic (`GlobalStartSkin − LineOfSightStartSkin`) meant for BabelBrain's own GUI workflow with oblique trajectories through a full-head domain. PRESTUS already computes the correct `TxMechanicalAdjustmentZ` itself (see above) and hands BabelBrain a beam-aligned mask, so this extra term silently moves the source off the requested position.
2. The (now-removed) resample anchor compensated for source mis-placement by shifting the paste-back to match the *source*, which correctly re-aligned the field to the transducer but dragged the (unchanged) tissue away from the original `medium_masks`.

**The fix.**

- `run_babelbrain_case.py` monkey-patches `SimulationConditionsBASE.CalculateDomainZReference` to return `0.0` (`_patch_zdomain_reference`, same policy as the existing skull-attenuation patch — no edits under `external/`), so the source is placed at exactly the requested `TxMechanicalAdjustmentZ`. Disable with the Python bridge's `--no_zdomain_fix` flag (not currently exposed as a `babelbrain.*` config option — same as `--no_skull_atten_fix`) to reproduce the un-patched behaviour.
- `read_babelbrain_output.m` / `read_babelbrain_thermal.m` anchor the resample on **`/TargetLocation`** — the label-5 focal-marker voxel PRESTUS itself wrote at `focus_pos` — rather than on the source position. This is the exact inverse of the forward mask transform, so the tissue map lands back on `medium_masks` by construction (verified on a benchmark run: overall label agreement 83% → 99%, bone 20% → 98%). With the source correctly placed (previous point), one anchor now aligns tissue *and* field simultaneously — the marker-anchored paste is source-consistent too.
- The passed-in `focal_length_mm` (preset `FocalLength`) is used only as a **fallback** when `/TargetLocation` is missing from the HDF5 (older files).
- A per-run sidecar (`<target>_prestus_anchor.json`) now records the **realized-minus-requested source position** as a diagnostic (not an anchor correction): with the `ZDomainStart` fix active, the only remaining residual should be BabelBrain's whole-voxel source snap-back (≤ 1 native voxel). `read_babelbrain_output.m` warns when the residual exceeds 1.5 native voxels, which points at a missing/disabled `ZDomainStart` patch rather than a coordinate bug.

This is fully general — no `medium_masks`/`segmentation` special-casing is needed anywhere downstream. `thermal_analysis.m` (and any other consumer that trusts `medium_masks`) automatically inherits correctly-aligned data because the BabelBrain acoustic step now returns output in the same PRESTUS-space convention k-Wave/jwave always provided.

## BabelBrain output format

The bridge writes `*DataForSim.h5` under `<dir_cache>/<subject_id>/`. `read_babelbrain_output.m` reads:

- `/p_amp` — peak positive pressure [Pa], shape `[Nx, Ny, Nz]`
- `/p_amp_refocus` — skull-corrected pressure (only when `do_refocusing: true`)

`p_max_all` in the returned `sensor_data` struct is set to `/p_amp` (or `/p_amp_refocus` when refocusing is active), matching the convention used by `acoustic_wrapper`.

## FOV/domain-coverage and focal-anchor debugging

`read_babelbrain_output.m` permutes and axis-flip-corrects the raw HDF5 field so that BabelBrain's
`[Nx, Ny, Nz]` axes line up with PRESTUS's `[lateral X, lateral Y, axial]` convention (dim 3 = beam
axis, matching `acoustic_grid_fov.m`'s `beam_dim = 3`). This axis handling is verified against the
shared tissue-layer ordering in benchmark axial profiles (e.g. `plot_axial_profile.m` output), which
lines up identically between the k-Wave and BabelBrain backends.

A BabelBrain axial profile with no reflections and exactly-zero pressure over a stretch nearest the
transducer/skull has two distinct possible causes that look identical in the final profile. **Inspect
the raw `*DataForSim.h5` directly before assuming either one** (`h5py`/`h5disp` on `/MaterialMap`,
`/p_amp`, `/TargetLocation`, `/z_vec`) — the BabelBrain log lines `ZShrink_L` / `ZShrink_R` (printed
by `BabelIntegrationBASE.py` during domain sizing) tell you immediately whether narrow-beam cropping
actually touched the proximal side for a given run; it is easy to assume it did when it didn't.

1. **Domain-coverage limit.** BabelBrain's own domain sizing (`--tight_narrow_beam` cropping when no
   PRESTUS-side FOV crop is active, or its baseline PML/padding-offset sizing when the grid is already
   FOV-cropped) determines how far its returned volume actually extends. `resample_babelbrain_vol.m`
   maps that volume onto the PRESTUS grid with **zero-fill for any voxel outside the BabelBrain domain**
   (`affine_resample_3d(..., method, 0)`), so a BabelBrain domain that genuinely stops short of the
   skull reads back as "no reflections" when the data there is simply absent rather than physically
   zero. Check the `ZShrink_L`/`ZShrink_R`/`XShrink_*`/`YShrink_*` log lines — nonzero values there
   only occur when narrow-beam mode is actually active (no PRESTUS-side FOV crop). When it is active,
   its ray-traced footprint estimate has been found to under-shoot the true beam width laterally near
   the skull (the aperture is wide near the transducer, narrowing only toward the focus) — cross-check
   against `BabelViscoInput.nii.gz`'s own voxel counts (does it actually contain skin/cortical/
   trabecular voxels?) vs. the exported `DataForSim.h5` `MaterialMap`'s voxel counts to see whether
   narrow-beam cropped away tissue that was genuinely present in the input. Also check the "BabelBrain
   native resolution" log line for a mismatch against `grid.resolution_mm` — BabelBrain reuses the
   input mask's voxel array shape directly as its own grid dimensions, so a resolution mismatch there
   silently rescales the whole domain independent of any shrink settings.

2. **Focal-anchor error.** The z-axis resample anchor is BabelBrain's own `/TargetLocation` — the
   index where it placed the label-5 focal marker `write_babelbrain_input.m` wrote into the input
   mask, i.e. exactly PRESTUS's `focus_pos` voxel (see
   [Output alignment](doc_babelbrain_internals.md#output-alignment-tissue-marker-anchoring-and-the-zdomainstart-fix)). The
   passed-in `focal_length_mm` (preset `FocalLength` / `TxMechanicalAdjustmentZ` / `ZSteering` chain
   — see [Transducer placement](#transducer-placement-zsteering-and-the-drive-phases)) is used only
   as a **fallback** when `/TargetLocation` is absent from the file (older BabelBrain HDF5s) — it is
   never used to override or cross-check the marker when present. `/TargetLocation` is computed by
   BabelBrain *before* it flips `p_amp`/`MaterialMap` along the depth axis, so reading it correctly
   means indexing the **raw, unflipped** `/z_vec`, not the flip-corrected coordinate array
   (`zdim`/`z_mm`) used for the data arrays — get this backwards and the anchor lands on the wrong
   voxel. When the fallback *is* in use and disagrees with what the marker would have given (only
   possible on files without `/TargetLocation`), the entire field is rigidly shifted along the beam
   axis when resampled — real, non-zero BabelBrain data (including skull reflections) landing at the
   wrong depth, some of it outside the PRESTUS grid (reads as the same "zero near the skull" symptom
   as case 1).

When `simulation.debug: 1` and a BabelBrain→PRESTUS resample occurs, `read_babelbrain_output.m` writes
diagnostics to `<dir_output>/debug/babelbrain/`:

- `sub-<id>_babelbrain_fov_debug<output_affix>.png` — a stacked figure:
  1. BabelBrain-native central-axis **pressure** profile (its own grid, before resampling), with the
     `/TargetLocation` tissue-marker anchor at `z=0` and, only when the passed fallback anchor
     differs from the marker actually used, a second line showing where that (unused) fallback would
     have placed it.
  2. The same axis, but for BabelBrain's own `/MaterialMap` (tissue index: 0=water, 1=skin,
     2=cortical, 3=trabecular, 4=brain) — shows directly where BabelBrain itself thinks the skull is,
     independent of any PRESTUS-side coordinate assumption.
  3. The resampled PRESTUS-grid pressure profile, with the region outside the BabelBrain domain
     shaded (zero-filled). The depth axis is beam-sign-corrected, so the shaded (truncated) side
     always corresponds to the actual proximal/distal side of the domain regardless of which way the
     beam runs along grid dimension 3.
- `sub-<id>_babelbrain_fov_debug<output_affix>.csv` — coverage + anchor numbers as a table:
  `bb_depth_rel_focus_min_mm` / `_max_mm`, `prestus_depth_rel_focus_min_mm` / `_max_mm`,
  `zero_filled_proximal_mm` / `_distal_mm`, and `anchor_correction_applied_mm` (how far the geometric
  estimate differs from `/TargetLocation`; normally ~0).
- `sub-<id>_babelbrain_material_map<output_affix>.nii.gz` — BabelBrain's `/MaterialMap`, resampled
  onto the PRESTUS grid with **nearest-neighbor** interpolation (never linear — this is integer label
  data) so it can be overlaid directly against `medium_masks` / the T1 in an external viewer.

Run-time console lines summarize both checks even without inspecting the files:

```
[read_babelbrain_output] FOV coverage check: BabelBrain domain covers depth-from-focus [...] mm; PRESTUS grid requires [...] mm.
[read_babelbrain_output] -> zero-filled (no BabelBrain data): X mm proximal (transducer side) / Y mm distal (far side) of the PRESTUS grid.
[read_babelbrain_output] Anchor check: focal_length_mm=... mm (estimate) vs /TargetLocation=... mm (BabelBrain ground truth, mismatch ... mm) — using /TargetLocation for the resample anchor.
```

## BabelBrain pipeline stages and what to check at each one

```
write_babelbrain_input.m        PRESTUS medium_masks → BabelBrain label NIfTI (0-5)
  Cortical/trabecular skull labels are derived from a single combined
  (skull+cortical+trabecular) smoothed mask in preproc_medium_mask.m, with the
  trabecular assignment constrained to that combined shell — trabecular cannot
  extend past the cortical envelope (e.g. touch skin directly). A central-line
  tissue-transition dump of the input NIfTI should never show skin→trabecular
  with no cortical between them; if it does, check the raw segmentation itself.
  The pre-flip (line ~132) must exactly cancel BabelBrain's own np.flip(axis=2)
  on load — verify both sides together if either changes.
  The NIfTI header (voxel -> RAS affine) is built by grid_to_ras_affine.m,
  composing planimg.inv_transf (head-alignment rotation), planimg.t1_header
  (T1's own affine), any acoustic-FOV-crop offset (planimg.fov_offset_ac,
  set in prestus_pipeline.m right after acoustic_fov_crop), and the dim-3
  flip above — in that order. This is what makes the written NIfTI overlay
  correctly onto the T1/medium_masks; a purely diagonal/origin-only affine
  (no rotation) will misalign whenever grid.mode rotates the head, which is
  the default.

run_babelbrain_case.py          Invokes RUN_SIM().RunCases(); patches material DB
  _inject_prestus_materials only runs when --medium_properties_json is passed AND
  --no_skull_atten_fix is NOT set — without the skull-attenuation fix, injected bone
  attenuation is never realised (QCorr stays 1.0), reproducing near-zero skull
  attenuation even though the Material table looks correctly populated.
  BabelBrain computes its own SpatialStep from PPW and reuses the input mask's voxel
  array shape directly as its own grid dimensions, without resampling to that
  SpatialStep — write_babelbrain_input.m resamples the mask to BabelBrain's own
  native resolution before writing so voxel spacing actually matches (see
  bb_native_res_mm there); check the log's "BabelBrain native resolution" line if the
  effective grid resolution is in question.
  _patch_zdomain_reference (default on; --no_zdomain_fix disables) forces
  CalculateDomainZReference to 0 so the source lands at exactly the requested
  TxMechanicalAdjustmentZ — see "Output alignment" above. Check the run log's
  "ZDomainStart neutralized: ... mm" line to see the magnitude BabelBrain would
  otherwise have added.

BabelIntegrationBASE.py         Domain sizing, Rayleigh forward-prop, FDTD, flip-before-save
  --tight_narrow_beam domain shrink (XShrink/YShrink/ZShrink_L/R): only active when
  no PRESTUS-side FOV crop is configured (babelbrain_wrapper.m disables it otherwise,
  since the mask is already cropped to a beam-centred box matching k-Wave/jwave).
  When active, check the run log's own `ZShrink_L <n>` / `ZShrink_R <n>` lines;
  ZShrink_L>0 means the proximal (transducer/skull) side was actually cropped. It is
  often 0 even with narrow-beam active, since only the far/beyond-focus side is
  targeted by z_beyond_focus_mm — don't assume cropping without checking the log.
  The lateral (XShrink/YShrink) crop has been found to under-estimate the beam's
  true width near the skull, cropping away skull tissue that is genuinely present
  in the input mask — this was the reason narrow-beam mode is disabled by default
  whenever a PRESTUS-side FOV crop already exists.
  /TargetLocation is BabelBrain's own record of where it placed the focal-spot
  marker PRESTUS wrote into the input mask — this IS the resample anchor (see
  "Output alignment" above), not merely a cross-check.
  p_amp/MaterialMap are np.flip(axis=2)'d before being added to DataForSim; z_vec is
  added AFTER that flip loop and is therefore NOT flipped — read_babelbrain_output.m's
  zdim = z_vec_bb(end:-1:1) compensates for exactly this asymmetry. Do not change one
  side without the other.

read_babelbrain_output.m / read_babelbrain_thermal.m   h5 read, flip-correction, resample
  Pressure/complex fields resample with 'linear' (resample_babelbrain_vol.m default);
  MaterialMap (an INTEGER label field) is resampled with method='mode' (majority-vote
  over the finer BabelBrain-native voxels landing in each PRESTUS cell) — 'nearest'
  would drop thin structures (e.g. sub-voxel cortical bone) that miss the single
  sample point.

babelbrain_wrapper.m             Free-water ISPPA baseline (separate --water_only run)
  The water-only baseline run is only regenerated if its Water_DataForSim.h5 doesn't
  already exist in the cache dir (is_water_medium branch) — a stale Water_DataForSim.h5
  left over from an earlier run (different transducer config, different geometry) will
  be silently reused for ISPPA calibration (unless io.overwrite_files is 'always'/'force',
  which deletes it first). If the ISPPA-normalized field looks wrong, check the
  Water_DataForSim.h5's mtime against the main DataForSim.h5's.
  The water-only run must receive the same --source_amp_pa as the transcranial run;
  without it BabelBrain drives at its ~50 kPa default and freefield_isppa_wcm2 comes out
  ~(156.6/50)^2 ≈ 10x too low for R15473 (the wrapper passes it since this fix — water
  baselines cached from before it are wrong and must be regenerated).
```
