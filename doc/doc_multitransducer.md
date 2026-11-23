# Multi-Transducer Modelling

When several transducers sonicate the same head, the exposure depends on how their fields combine, which a run per transducer does not capture. PRESTUS simulates them in one layered setup: transducers that fire together (`coherent`) add as complex pressures and interfere, while for transducers that fire in separate windows (`async`) acoustic safety metrics take the per-voxel maximum of the fields and heating uses the sum of their intensities.

On this page:

- [How the Fields Combine](#how-the-fields-combine): an interactive comparison of the coupling modes.
- [Configuration](#configuration): listing transducers and choosing how they fire.
- [Grid Setup](#grid-setup): how the grid is aligned and cropped around several beams.
- [I<sub>SPPA</sub> Targets](#isppa-targets): per-transducer free-water targets and their baselines.
- [Complex Pressure Field](#complex-pressure-field): saving fields and superposing them afterwards at any relative phase.
- [Transducers at Different Frequencies](#transducers-at-different-frequencies): firing together or in turn at more than one frequency.
- [Partial Acoustic FOV](#partial-acoustic-fov): limited fields of view with several transducers.
- [Known Limitations](#known-limitations): what refers to the first transducer only, and unsupported setups.
- [Subpages](#subpages): async firing in separate windows.

!!! note "At a glance"
    - **When you need it:** Bilateral or dual-site stimulation, or any protocol with more than one transducer on the head.
    - **Before you start:** A layered (skull) simulation in a 3-D grid; one equipment entry and placement per transducer.

!!! note "Experimental support"
    See [this pull request](https://github.com/Donders-Institute/PRESTUS/pull/100) for examples.

## How the Fields Combine

One focused transducer as the reference, then two on opposite sides firing together at one frequency, together at different frequencies, or in turn. Switch between the setups, and for one frequency change the phase of transducer 2; the readout gives the peak pressure and time-averaged intensity at the focus relative to one transducer.

<div class="mt-fields"></div>

The illustration is schematic: a lossless, homogeneous 2-D medium, slowed down, with pulses shortened in the async setup. The combination rules are the ones PRESTUS applies; see [Transducers at different frequencies](#transducers-at-different-frequencies) and [Async Multi-Transducer Mode](doc_multitransducer_async.md).

## Configuration

List several transducers in a layered simulation:

```yaml
transducer:
  - name: right
    serial: IS_PCD15473_01001
    combo: {ds_serial: IGT_32_ch_comb_10_ch}
    # ... trans_pos, focus_pos, focal_distance_ep

  - name: left
    serial: IS_PCD15473_01001
    combo: {ds_serial: IGT_32_ch_comb_10_ch}
    # ...
```

`simulation.transducer_coupling` sets how the transducers fire:

| Firing | Setting |
|---|---|
| Simultaneously, same frequency, known phase | `coherent` (default): one [k-Wave](https://github.com/ucl-bug/k-wave) run with all sources active; the fields interfere. All transducers must share the fundamental frequency (a warning is issued otherwise). |
| In separate, non-overlapping windows | `async`: see [Async Multi-Transducer Mode](doc_multitransducer_async.md). |
| In sequence, or with partial overlap | Use [Sequential Simulations](doc_sequential.md): carry each heating run's temperature and CEM43 maps into the next via `adopted_heatmap` / `adopted_cem43`. This assumes the acoustic fields do not overlap in time. |

In the Nexus, *Multi-transducer coupling* (Sonication tab, shown once the config has more than one transducer) selects the row above, and *Target free-water I<sub>SPPA</sub>* is set per transducer on the Transducer tab. The card can show an analytic prediction of both forms: the at-target intensity when the beams add independently and its range over their relative phase when coherent, and free-water maps of both in the beams' plane. The maps use Gaussian focal spots from each transducer's aperture and focal distance, with no skull and no side lobes; they guide the choice and do not replace the simulation.

## Grid Setup

### Grid Alignment

The default `focal_axis` mode rotates the grid so that the mean beam axis of all transducers aligns with z. Beams are averaged as axes, not directions, so opposed (bilateral) beams do not cancel: z bisects their two axes. Each beam off z is drawn obliquely on the grid, which costs accuracy at coarse resolution (for async runs, see [`grid.async_layout: per_transducer`](doc_multitransducer_async.md#grid-geometry-mode)).

```yaml
grid:
  mode: 'focal_axis'   # default; omit to get the same behaviour
```

With three or more transducers whose axes are too spread to average (mean-axis norm below `grid.focal_axis_min_norm`, default 0.5), `focal_axis` falls back to transducer(1)'s axis and warns.

To align the grid to transducer(1)'s focal axis instead, set `transducer_axis`. `ras_plus` keeps the RAS+ orientation without rotation (isotropic rescaling only), so every beam is oblique; it does not work with `grid.axisymmetric`. See [Preprocessing: Orientation](doc_head-model.md#orientation) for the mode table. In the Nexus, *Placement ▸ Transducers* shows the grid z axis and each beam's tilt for the chosen mode.

### Grid Cropping

The grid is cropped to the bounding box that encloses the CSF-expanded head tissue **and the bowl geometry of every transducer**. See [Preprocessing: Cropping](doc_head-model.md#cropping).

## I<sub>SPPA</sub> Targets

`transducer(k).target_isppa_wcm2` is the free-water [I<sub>SPPA</sub>](doc_glossary.md#isppa) of **transducer k alone**, for every coupling mode. Each transducer has its own drive (`elem_p_amp`), so the target rescales that transducer's drive, not the combined field.

```yaml
simulation:
  transducer_coupling: coherent
transducer:
  - name: left
    serial: IS_PCD15473_01001
    target_isppa_wcm2: 20     # left transducer alone: 20 W/cm² in free water
  - name: right
    serial: IS_PCD15473_01003
    target_isppa_wcm2: 40     # right transducer alone: 40 W/cm²
```

| Coupling | How the targets are applied | Sweeps (vector targets) |
|---|---|---|
| ![](assets/icons/icon_coupling-coherent_glyph.svg){ .choice-icon } `coherent` | Before the simulation (`coherent_drive_calibration`): each transducer's `elem_p_amp` is multiplied by √(target / baseline), where baseline is the free-water I<sub>SPPA</sub> of that transducer alone at its configured drive. The combined field is simulated at these drives and not rescaled afterwards. | Not available. Errors (`prestus:isppa:coherentSweep`); run one config per point, or use `async`. |
| ![](assets/icons/icon_coupling-async_glyph.svg){ .choice-icon } `async` | After the per-transducer simulations: each field is scaled to its own target in `combine_async_intensity`, then combined. | Yes: one acoustic run per transducer, one combine + thermal pair per sweep point (see [Async Multi-Transducer Mode](doc_multitransducer_async.md)). |

**Coherent baseline.** The rule is the same as for one transducer ([I<sub>SPPA</sub> targeting](doc_simulations-acoustic.md#isppa-targeting)). If the target matches the transducer's calibration I<sub>SPPA</sub> (`tx_calibration.isppa_wcm2`) within 1 %, the calibration is the baseline. Otherwise PRESTUS simulates that transducer alone in water (others at zero amplitude, same grid). With N transducers that is up to N extra free-water simulations. A transducer with `target_isppa_wcm2: []` keeps its configured drive.

??? note "Provenance and cache checks"
    The result records targets, baselines, their source, and the drive scale in `acoustic_provenance.coherent_drive`. The combined free-water I<sub>SPPA</sub> (both firing) is reported but is not a target. A coherent cache without `coherent_drive` is refused when targets are set (`prestus:isppa:coherentWithoutDriveCalibration`): re-run its acoustic simulation. `path.acoustic_input` cannot be rescaled per transducer and errors when targets are set.

## Complex Pressure Field

```yaml
io:
  save_p_complex: true
```

This writes `sub-XXX_<medium>_p_complex<affix>.mat` with the complex amplitude at the fundamental on the simulation grid, at the simulated drive level, before I<sub>SPPA</sub> rescaling. Propagation in PRESTUS is linear, so fields simulated one at a time can be summed. `superpose_complex_fields` does this on one shared grid, back-projecting partial-FOV fields:

```matlab
out = superpose_complex_fields({'.../sub-001_layered_p_complex_trans-01.mat', ...
                                '.../sub-001_layered_p_complex_trans-02.mat'});
out.p_coherent    % |p1 + p2|          coherent (simultaneous) firing
out.p_incoherent  % sqrt(|p1|^2+|p2|^2) async heat source
out.p_max         % max(|p1|, |p2|)    async peak pressure
out = superpose_complex_fields(files, Weights=[1, exp(1i*pi/2)]);   % relative phase of transducer 2
```

One async run (`grid.async_layout: shared`) with `io.save_p_complex: true` therefore also gives the coherent field for any relative amplitude and phase, without another simulation. Fields on different grids (`per_transducer` mode) and fields at different frequencies are not combined. With k-Wave (`acoustic_backend: kwave` or `kwave_elastic`, every `code_type`), `p_complex` is fitted to the pressure recorded over the last quarter period, which keeps the recording small (about 3 GB on a 0.75 mm head grid). k-Wave takes `p_max_all` and `p_min_all` over the same window, so in runs with `save_p_complex` they are ± the fitted amplitude, and I<sub>SPPA</sub>, MI and heating follow from it. In a steady field the two agree within 0.1 % where the pressure exceeds half its peak and within 0.5 % above a tenth of it; in low-pressure regions the fitted amplitude can be a few percent lower than a peak over several periods. Compare runs with and without `save_p_complex` with this in mind. Multi-frequency setups record the whole steady-state window instead (the tones are only separated over several periods) and stop above 8 GB (`PRESTUS:pComplexTooLarge`); use a [partial FOV](doc_partial-fov.md).

The other backends write the same file without these limitations:

| Backend | How `p_complex` is computed | `p_max_all` | Phase reference |
|---|---|---|---|
| k-Wave (`kwave`, `kwave_elastic`) | fit over the last quarter period | ± fitted amplitude | t = 0 |
| jwave `pseudo_spectral`, `finite_differences` | running fit over the steady-state window, no time series stored | maximum over the window (unchanged) | t = 0 |
| jwave `helmholtz` | the solver's complex field | \|P\| (unchanged) | t = 0 |
| BabelBrain | BabelBrain's own field (one frequency bin over its last cycles) | BabelBrain's amplitude (unchanged) | BabelBrain's recording window |

`superpose_complex_fields` adds fields referred to t = 0 (k-Wave and jwave runs, in any mix), or BabelBrain fields with each other; it refuses to mix the two.

??? note "File contents"
    - `p_complex`: complex amplitude at the fundamental on the simulation grid (the acoustic FOV sub-grid with a [partial FOV](doc_partial-fov.md)), at the simulated drive level, before I<sub>SPPA</sub> rescaling. The phase refers to t = 0 of the simulation, when every source starts.
    - `p_complex_info`: `phase_reference` (`'simulation_t0'`), `freq_hz`, `record_t0_s`, `dt_s`, `method` (`'quarter_period'` or, multi-frequency, `'steady_state_window'`), `n_samples`, `p_max_from_fit` (whether `p_max_all`/`p_min_all` are ± the fitted amplitude), `drive_level` (`'simulated'`), `grid_dims`, and `fov_offset_ac` / `full_ac_dims` with a partial FOV.

    Files without `p_complex_info` are accepted only all together, with a warning (valid only for runs with the same time step and length).

## Transducers at Different Frequencies

When the transducers do not all share one frequency, each is simulated alone, at its own frequency and on a medium mapped for that frequency, and the results are combined. With `transducer_coupling: async` they fire in separate windows ([Async Multi-Transducer Mode](doc_multitransducer_async.md)). With `coherent` (the default) they fire together; this uses the same pipeline as async mode, on a shared grid (`grid.async_layout: shared`), with complex fields saved automatically. Runs where every transducer has the same frequency are unchanged: one joint simulation.

Firing together, the pressure is the sum of the complex fields, p(t) = Re(Σ_f P_f·e^{iω_f t}), where P_f is the coherent sum of the fields of the transducers at frequency f. Transducers at one frequency interfere with a fixed relative phase, so their fields add as complex numbers. Between different frequencies the relative phase drifts at the beat frequency |f1 − f2|, so the cross terms average out of the time-averaged quantities but not out of the peaks:

| Output | Firing together, more than one frequency |
|---|---|
| Pulse-average and time-average intensity (I<sub>SPPA</sub>, I<sub>SPTA</sub>) | Σ_f \|P_f\|²/(2ρc) |
| Heat source | Σ_f α(f)·\|P_f\|²/(ρc), all on together; α(f) from each tissue's power law |
| Peak pressure, MI | maximum and minimum of p(t) over the pulse or the common period 1/gcd(f), whichever is shorter; MI uses the lowest frequency |
| Temporal-peak intensity (Itp) | from the peak pressure |

When the common period needs more than 20000 time samples (frequencies close together or in a high-order ratio), the peak pressures are the upper bound ±Σ_f|P_f|; `acoustic_provenance.combination_peak` records which. Each pulse is assumed to start with the phases of t = 0. Nonlinear propagation (harmonics, sum and difference frequencies) is not modelled, as in every PRESTUS run.

??? developer "Implementation"
    `is_simultaneous_multifreq` routes these runs to `async_transducer_pipeline`. `combine_async_intensity` (mode `'simultaneous'`) sums the scaled complex fields per frequency and writes `p_intensity_equiv` (intensity sum, read by `acoustic_analysis` for Ipa and Ita), `p_max_async`/`p_min_async` (phase-sampled peaks), and `p_async_tx` with one field |P_f| per frequency (`p_async_freqs_hz`). `thermal_simulation` heats each field with its frequency's absorption, all on together (`thermal_schedule` mode `'simultaneous'`). Validation against one joint simulation with both tones: `code/report/08_multitransducer/multitransducer_multifreq_*` in PRESTUS_dev.

## Partial Acoustic FOV

A [partial acoustic FOV](doc_partial-fov.md) (`grid.acoustic_fov_diameter_mm`, `grid.acoustic_fov_length_mm`) works with several transducers:

- **Coherent**: the FOV box is the union of one box per transducer, each enclosing a cylinder of the FOV diameter along that transducer's beam axis. `'auto'` is resolved per transducer (its outer element diameter; its focal depth + `grid.z_beyond_focus_mm`). For opposite transducers the union spans most of the head, so the saving is mostly lateral.
- **Async**: each transducer's run gets its own FOV box. Fields are back-projected onto their full grid before combining, in both grid modes. Partial-FOV and full-grid caches can be mixed.
- **Complex superposition**: each field is zero outside its own box, so the sum is complete only where boxes overlap (`out.coverage` counts the fields per voxel). Make the boxes cover the region around both foci.

The Nexus placement viewer draws every transducer's FOV box.

## Known Limitations

- Simulations without a skull or layered setup (e.g. `medium: water`) only model the first transducer.
- `kWaveArray` mode supports multiple transducers only in non-axisymmetric mode (`grid.axisymmetric = 1` is single-transducer). Matrix-type transducers in `kWaveArray` mode assume a +Z beam axis; see [Transducer geometry](doc_transducer.md#beam-axis-constraint).
- **Different carrier frequencies**: see [Transducers at Different Frequencies](#transducers-at-different-frequencies).
- Exit-plane metrics (`real_focal_distance_ep_mm`, `real_focal_distance_bowl_mm`, `max_Isppa_after_exit_plane`) refer to the first transducer only.
- Focal-plane heating plots show only the first transducer's focal plane. Always inspect the 3D `maxT` and `CEM43` NIfTI volumes for a complete safety assessment.
- Post-hoc water-only simulations are not implemented for multiple transducers.
- The PRESTUS **[BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain) backend** (`acoustic_backend` / `thermal_backend: babelbrain`) supports one transducer per run. Multi-transducer runs need `kwave` or `jwave` with `thermal_backend: kwave`. See [BabelBrain: Geometry limits](doc_babelbrain.md#geometry-limits).

??? note "Comparison with BabelBrain"
    BabelBrain models one transducer per simulation. It treats several foci of one phased array as time-multiplexed, like PRESTUS async mode: the safety field is the per-voxel maximum of the pressure amplitudes (`CalculateFieldProcess.py`), and the thermal simulation switches between the spots' fields by duty-cycle share (`BHTEMultiplePressureFields` in `ThermalModeling/CalculateTemperatureEffects.py`). It does not add complex fields across foci.

## Subpages

- [Async Multi-Transducer Mode](doc_multitransducer_async.md): for transducers that fire in separate, non-overlapping windows; covers configuration, the physical model, grid modes, outputs, and constraints.

See also: [Transducer Calibration](doc_calibration.md), [Multi-I<sub>SPPA</sub> Mode](doc_multiintensity.md), [Sequential Simulations](doc_sequential.md), [Pipeline Modules](doc_modules.md), [Transducer Definition](doc_transducer.md)
