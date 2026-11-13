# Time Reversal: Refocusing, Aberration Report and Pose Search

!!! warning "Draft"
    This feature is new and not yet validated against forward simulations. Treat its predictions as a ranking signal and check any pose or phase set you intend to use with a regular acoustic simulation.

!!! note "At a glance"
    - **What:** One extra reverse simulation per target that reports how well the transducer focuses through the skull, gives refocusing phases for arrays, and scores alternative transducer poses.
    - **When you need it:** When the focus through the skull is weaker or shifted and you want to correct phases (arrays) or choose a better pose (single element or annular array).
    - **Before you start:** A working acoustic simulation of the planned pose on a 3D Cartesian grid.

The skull delays and attenuates the parts of the beam differently, so the waves from different parts of the transducer no longer arrive at the target in phase: the focus weakens and shifts. With `time_reversal.enabled: true`, PRESTUS measures this with one extra simulation per target, using acoustic reciprocity.

## What to set

```yaml
time_reversal:
  enabled: true          # reverse run: aberration report and refocusing phases
  apply_phases: false    # true: use the time-reversal phases in the acoustic run
  placement_search:
    enabled: false       # true: score poses around the planned one (curved annular transducers, 3D grid)
  pose_workflow:
    enabled: false       # true: search, forward verification and final run in one call
```

Outputs per transducer: `<prefix>_timereversal_summary.csv` (metrics), `<prefix>_timereversal_elements.csv` (relative transfer, current and time-reversal phase, and the phase error with which each element arrives at the target) and, with a pose search, `<prefix>_timereversal_poses_T<i>.csv`.

## What each transducer gets

| Transducer | Phases it can use | Report |
|---|---|---|
| Matrix array | One per element: full refocusing | All metrics; `apply_phases` uses the per-element phases |
| Annular array (`elem_n > 1`) | One per ring: the axially symmetric part of the correction | All metrics; `apply_phases` uses the ring phases |
| Single element | None | Coherence and predicted focal shift; nothing to apply |

The metrics in `<prefix>_timereversal_summary.csv`:

- **coherence_elements**: |Σ contributions| / Σ |contributions| over the elements with the current phases. 1 means the elements arrive in phase; `refocus_gain_pressure` = 1 / coherence_elements is what refocusing gains at the target.
- **coherence_aperture**: the same over the whole surface, including the phase spread within each element, which no element phase can correct.
- **shift_*_mm**: predicted focal shift with the current drive (axial > 0: deeper than the target), and with the time-reversal drive (`shift_refocused_*`).

??? note "Details: how the reverse run works"
    A point source is placed at the target and run as a continuous wave on the grid, medium and time axis of the acoustic simulation (`reverse_point_source_field`). By reciprocity, the pressure this source produces at a point *x* equals the pressure that a source at *x* would produce at the target (up to one constant common to all points). Recorded on the transducer surface, the reverse field therefore gives, for every element, the amplitude and phase with which it reaches the target through the skull (`time_reversal_refocus`):

    ```
    G_e = sum over the element's grid points of (weight × reverse field)
    pressure at the target = sum_e  A_e · exp(i·phi_e) · G_e
    ```

    Driving each element with the phase `−angle(G_e)` (phase conjugation, i.e. time reversal) makes all elements arrive in phase: the largest target pressure the element amplitudes allow.

    **Predicted focal shift.** Near the target the forward field is approximated by `F(y) = Σ_j q_j · exp(−i·k·(|x_j − y| − |x_j − t|))`: the skull aberration seen from the target is taken to be the same for points within a few millimetres of it (isoplanatic patch), and only the path through the tissue around the target changes. The peak of |F| gives the shift (`predict_focal_shift`). The approximation degrades with the size of the shift and the strength of the aberration: large predicted shifts mean "this pose focuses badly", not an exact focus position.

## Applying the phases

With `apply_phases: true`, the time-reversal phases replace the element phases before the acoustic simulation, and the source is rebuilt. The free-water baseline has already run with the geometric phases, so ISPPA targeting (`target_isppa_wcm2`) keeps the element drive amplitudes of the geometric drive: the refocused run shows what the same transducer output delivers when refocused. Whether your hardware can apply arbitrary phases (and with which precision) is not checked.

## Pose search (aberration-corrected placement)

For transducers that cannot refocus, the pose is the main lever: where the transducer sits on the scalp and where it points. With `placement_search.enabled: true` (curved annular transducers, 3D grid), the same reverse run also scores poses around the planned one. Each pose costs an interpolation of the recorded field, not a simulation.

The result is `<prefix>_timereversal_poses_T<i>.csv`, sorted by `p_corrected_rel` (target pressure relative to the planned pose with geometric aim), with the apex, corrected aim point and target in grid voxels, T1 voxels and RAS mm. A pose is set up as a manual placement with `trans_pos` = apex and `focus_pos` = aim point (`placement.manual.coordinate_space: ras`, `transducer.focal_distance_ep: []`). The acoustic run keeps the planned pose. To verify the best poses with forward simulations and run the chosen one, use the [pose workflow](#pose-workflow), or set up a pose by hand with `time_reversal_pose_parameters(parameters, pose_row)`.

`time_reversal.only: true` stops the run after this stage (no water baseline, acoustic or thermal simulation), for a report or pose table on its own.

??? note "Details: candidate poses and scoring"
    **Candidates** (`reciprocity_pose_candidates`): beam directions tilted from the planned one up to `max_angle_deg`, in steps of `angle_step_deg`, all aimed at the target. For each direction the exit plane sits `placement.skin_gap_mm` off the scalp (heuristic, MNI, PlanTUS and TUS_entry placements) or at the planned pose's gap (manual and Localite), so the required focal depth changes with the direction. Where the scalp curves, the bowl is moved outward until its surface overlaps the scalp no more than the planned bowl does (`push_mm` in the pose table), so all poses are seated alike. A single element cannot steer: its apex is placed at the curvature radius from the target instead, and its aim correction is lateral only. With `jitter_mm` (one or more distances, e.g. `0.5` or `[0.5 1 2 5]`), the planned pose with its apex moved by each distance is scored too (rows of kind `jitter`), the local search often run with forward simulations. `jitter_mode: grid` moves the apex ± along each grid axis, as when editing voxel coordinates; `jitter_mode: scalp` slides the transducer over the scalp perpendicular to the beam in 8 directions and seats it against the scalp, which suits shifts beyond about 1 mm. Poses whose surface overlaps the head by more than the planned surface does plus 2 % of its points are dropped (`inside_fraction`); the planned pose is always kept. The reverse run records the field over the region these poses cover.

    **Scoring** (`reciprocity_score_poses`): for each pose, the surface is sampled (`annular_aperture_points`) and the target pressure is a sum over the reverse field there. Per pose:

    - target pressure and coherence with geometric aim and steering;
    - predicted focal shift;
    - **aberration-corrected aim**: aim and steer at *target − predicted shift* (tilting the bowl about its apex), repeated `aim_iterations` times, so that the predicted focus lands on the target; the target pressure after this correction ranks the poses;
    - the gain ring refocusing would add.

    **Grid.** The simulation grid ends close to the planned transducer, so tilted poses can lie outside it. The reverse run of a pose search therefore uses its own grid: the full grid (also with a [limited acoustic FOV](doc_limited_fov.md), whose crop the forward run keeps), padded with water on the sides the candidate poses need (the log reports the padding). The padding makes the reverse run larger, and so longer, than the forward run.

## Pose workflow

`time_reversal.pose_workflow.enabled: true` turns a run into the whole placement workflow, in one call of `prestus_pipeline_start` (`pose_workflow_pipeline`). It needs one curved annular transducer (ring array or single element) and a head model (`simulation.medium: layered`). The `placement_search` settings apply; `placement_search.enabled` is set by the workflow.

1. **Search** (`<affix>_posesearch`): the reverse run on the full grid scores the candidate poses; no forward simulation.
2. **Verify** (`<affix>_pose<id>_<geo|aim>`): forward runs, acoustics only, of the planned pose (the reference, aimed at the target) and the `n_top` best poses by predicted target pressure after aim correction, plus the pose ids in `poses`. `aim` sets the aim of each run: `corrected` (at the aberration-corrected aim point; the planned pose is then run with it too), `geometric` (at the target) or `both`. `verify_fov: partial` (default) runs them on an acoustic FOV around each beam. The poses run are listed in `<prefix><affix>_poseverify_selection.csv`.
3. **Results and final run**: per verification run, the mean intensity in a sphere of `target_radius_mm` around the target and the distance of the −3 dB focal centroid and peak from the target (axial and lateral), next to the predicted values, in `<prefix><affix>_poseverify_results.csv`. The runs are ranked by `rank_by` (`target_intensity` or `focus_error`), among those within `max_focus_error_mm` of the target if set; the log reports the rank agreement between predicted and simulated intensity. With `final_run: true`, the best pose then runs on the full grid with the configured modules (thermal, report): `<affix>_posebest`.

On MATLAB the stages run one after the other. On SLURM or qsub, each stage is submitted: a small CPU job after the search job submits the verification runs, and another after those (whether they succeeded or not) collects the results and submits the final run. This needs a cluster that lets jobs submit jobs. A stage whose outputs exist is skipped unless `io.overwrite_files: always`, so calling the workflow again after a failure continues where it stopped.

The verification and final runs place the transducer as a manual placement in RAS mm (apex = `trans_pos`, aim point = `focus_pos`), with the focal depth of the pose:

- a drive from the calibration library (`transducer.serial`) is resolved again for the pose's focal depth (calibrated phases and amplitudes);
- inline ring phases are replaced by geometric steering to the aim point in water, the drive the pose scores assume;
- a single element has no phases; its poses already place the geometric focus on the target or the corrected aim point.

The runs keep the config's intensity setting. With `target_isppa_wcm2`, each run is scaled to the same free-water ISPPA, so the simulated intensities compare poses at the same device setting; the predictions compare them at the same element drive. The two differ where the poses' focal depths differ much (the free-water focal intensity changes with steering depth).

Cost: one reverse run (as long as a full-grid acoustic run, more memory: `search_memorylimit`), `n_top + 1` (or twice that with `aim: both`) partial-FOV forward runs, and one full run.

## Cost

One k-Wave run per distinct target, about as long as the acoustic simulation (longer by the padding of a pose search), recording the pressure time series (`grid.steady_state_num_periods` periods) at the sensor points: the transducer surface (refocusing), plus the candidate region (pose search). Scoring takes seconds per pose.

Memory stays below that of the acoustic simulation: the forward source is released during the stage and rebuilt afterwards, the padded grid reuses the forward source's element weights (no second source is built), and the recording is reduced to one complex amplitude per point in blocks, without a double-precision copy (0.5 mm head runs of the forward simulation peak at 25–31 GB).

## Assumptions and limits

- Linear, fluid acoustics. With `acoustic_backend: kwave_elastic` the reverse run uses the fluid medium (warning). jwave forward runs use the same medium, so the k-Wave reverse run applies.
- The transducer surface lies in one medium (coupling water); k-Wave scales a pressure source by the local sound speed.
- Cartesian grids only (no axisymmetric runs).
- Not with `acoustic_backend: babelbrain`, which has its own refocusing (`babelbrain.do_refocusing`, phased-array presets only).
- The pose search ignores calibrated phase offsets of the library entry (it steers geometrically in water) and treats each pose as a rigid bowl. The verification runs use the library drive, so they include the offsets.
- The pose scores interpolate the recorded field (cubic) at the transducer surface; at ~7 points per wavelength the interpolation error is about 1 % of the target pressure, the size of the differences between ±0.5 mm poses.

??? note "References: related work and benchmarks"
    - Pinton (2026), *Whole-skull acoustic transparency from a single time-reversal solve for reciprocity-based transducer placement and aperture optimization in transcranial focused ultrasound*, [arXiv:2609.27852](https://arxiv.org/abs/2609.27852): the same principle (one time-reversal solve, then placement, aperture and orientation searched on the recorded map), with Fullwave 2 on a CT skull.
    - BabelBrain refocusing: phase correction for phased arrays from a back-propagated field (`babelbrain.do_refocusing`).
    - SCOUT (2024, [PMC11592240](https://pmc.ncbi.nlm.nih.gov/articles/PMC11592240/)): position of a single-element transducer optimised with repeated forward k-Wave simulations.
    - Park et al. (2019), *Comput Methods Programs Biomed* 179: placement by the lowest average reflection coefficient along straight beam lines (no wave simulation).
    - [PlanTUS](doc_placement_plantus.md): geometric placement criteria (skull angle, thickness, distance).
    - Benchmarks (predictions vs forward simulations): `code/report/12_time_reversal` of the PRESTUS v05 project.
