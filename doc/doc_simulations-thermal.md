# Thermal Simulations

Tissue absorbs part of the ultrasound and heats up, and temperature rise and thermal dose (CEM43) are the measures safety guidelines use for this effect. PRESTUS computes both from the acoustic field and the sonication protocol, so the predicted heating depends on the timing, step sizes and absorption settings described here. The report compares the results with published reference levels; PRESTUS does not make an independent safety recommendation.

On this page:

- [Protocol Timing](#protocol-timing): describing pulses, pulse trains, repetitions and the post-stimulation period.
- [Step, Temperature and Sensor Settings](#step-temperature-and-sensor-settings): step sizes, starting temperature and the temperature sensor.
- [Heat Deposition](#heat-deposition): how tissue absorption turns the acoustic field into heat.
- [Thermal Dose (CEM43) Calculation](#thermal-dose-cem43-calculation): the two CEM43 formulations and their outputs.

!!! note "At a glance"
    - **When you need it:** Whenever you simulate heating (`modules.run_heating_sims: 1`).
    - **Before you start:** A working acoustic simulation, and the timing of your protocol.

Thermal simulations use the output of the acoustic simulation; they cannot run without one. They run in [k-Wave](https://github.com/ucl-bug/k-wave) by default. With `acoustic_backend: babelbrain`, [BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain)'s own [Pennes](https://doi.org/10.1152/jappl.1948.1.2.93) solver can be used instead (`thermal_backend: babelbrain`). See [Simulation backends](doc_backend.md) for valid combinations and how the results compare, and [BabelBrain thermal](doc_babelbrain.md#thermal-simulation) for the timing mapping and tissue-property differences. **Known limitation:** BabelBrain thermal uses its own built-in tissue thermal properties, so the thermal fields in `medium_properties` apply only to the k-Wave thermal backend.

## Protocol Timing

Protocol timing follows the nomenclature of the [TUS Calculator](https://www.socsci.ru.nl/fusinitiative/tuscalculator/), which you can use to specify the protocol. All three levels of nesting (pulse, pulse train, pulse train repetition) must be set. For a single pulse train without repetitions, set `ptri` and `ptrd` equal to `ptd`.

Example: a 20 s, 5 Hz, 10 % duty-cycle protocol as a single pulse train repetition:

```yaml
timing:
  pd: 0.02               # pulse duration (PD) [s]
  pri: 0.2               # pulse repetition interval (PRI) [s]
  ptd: 20                # pulse train duration (PTD) [s]
  pt_timestep: 0.02      # modelling time step inside a PT, with equal_step_duration: 1 [s]
  ptri: 20               # pulse train repetition interval (PTRI) [s]
  ptrd: 20               # pulse train repetition duration (PTRD) [s]
  post_ptri_dur: 20      # steady-state duration following the PTRD [s]
  post_pt_timestep: 1    # modelling time step between PTs and after the PTRD [s]
  equal_step_duration: 0 # 0: one step per on and off period; 1: fixed step duration pt_timestep
```

Continuous pulse train protocols can be modelled either as a single pulse train repetition with multiple pulse trains (as above) or as multiple pulse train repetitions of a single pulse train.

A coarser protocol neglects the heating dynamics between pulses and models one longer pulse per train. Collapsing the off period lets the simulation use a coarser step size there:

```yaml
timing:
  pd: 2
  pri: 2
  ptd: 2
  pt_timestep: 0.02
  ptri: 20
  ptrd: 20
  post_ptri_dur: 0
  post_pt_timestep: 1
```

**Post-stimulation period.** Heat dissipation continues after the final pulse. Run the simulation until a steady state is reached by setting `timing.post_ptri_dur`; its time step is `timing.post_pt_timestep`, separate from `pt_timestep`. For a continuous pulse train, a steady-state phase can also be part of an extended PTR interval; `post_ptri_dur` also gives protocols with repetitions a steady-state phase. For long intervals between pulse trains, a coarser `post_pt_timestep` may be sufficient.

## Step, Temperature and Sensor Settings

| Parameter | Role |
|---|---|
| `timing.pt_timestep` | Step size for the on and off periods of each pulse within a pulse train. Used only with `equal_step_duration: 1`; PD and PRI − PD must be integer multiples of it. |
| `timing.equal_step_duration` | `1`: every on and off period is split into steps of `pt_timestep`. In the first example above (PD 20 ms, PRI 200 ms, PTD 20 s, `pt_timestep` 20 ms), each of the 100 pulses takes 1 step with stimulation and 9 without (10 × 100 = 1,000 steps). `0` (default): every on and off period is one simulation step (2 × 100 = 200 steps), which reduces the computational load. |
| `timing.post_ptri_dur` | Duration of the post-stimulation period after the PTRD. |
| `thermal.temp_0` | Starting temperature in degrees, per tissue. |
| `thermal.sensor_xy_halfsize` | Half-size (grid units) of the rectangular XY sensor window centred on the transducer position. Only voxels in this window enter the temperature sensor mask passed to `kWaveDiffusion`; the full Z extent is always included. A smaller value lowers memory use. |
| `thermal.record_t_at_every_step` | `0` (default): the sensor is discarded before the simulation and only the final temperature field is kept. `1`: the full sensor window is stored at every time step; needed for the complete temperature time series, but can cause `out of memory` errors for large grids or many time steps. |

## Heat Deposition

The acoustic field deposits heat at the rate Q = 2·α_abs·I = α_abs·p²/(ρc), with the absorption coefficient α_abs [Np/m] at the source frequency. Each tissue in `medium_properties` specifies α_abs in one of two ways (see [Medium properties](doc_medium.md)):

- `absorption_fraction`: a fraction of the attenuation coefficient. The defaults use this for every tissue, so a change of `alpha_coeff` also changes absorption.
- `absorption_coeff` [dB/(cm·MHz^y)]: a fixed coefficient, with the tissue's attenuation power law y (`alpha_power`). It does not follow changes of `alpha_coeff`.

A tissue sets one of the two; the other stays empty. When a config layer (a study [YAML](https://yaml.org/) file, a struct passed to `load_parameters`, or an override) sets one of them for a tissue, `load_parameters` empties the other one inherited from earlier layers, such as the default `absorption_fraction`. One run can combine both. For example, brain absorption as a fraction and skull absorption as a coefficient:

```yaml
medium_properties:
  brain:
    absorption_fraction: 0.77
  skull:
    absorption_coeff: 3.78   # clears the default skull absorption_fraction
```

A layer that sets both for one tissue stops with the error `medium_absorption_layer:bothSet`. If you change the struct in [MATLAB](https://www.mathworks.com/products/matlab.html) after loading, empty the other field yourself (`parameters.medium_properties.skull.absorption_fraction = [];`); otherwise the medium setup stops with `thermal_absorption_model:bothSet`. Empty and `NaN` values count as not set; a tissue with neither uses a fraction of 1.

A pseudo-CT skull with mapped attenuation has no matching absorption map. A skull `absorption_coeff` is then applied as the fraction `absorption_coeff` / `alpha_coeff` of the mapped attenuation, with a warning (`medium_setup:absorptionModelForced`).

The log, the run summary and the HTML methods text give each tissue's model with its effective coefficient and fraction (α_abs/α) at the configured attenuation. The medium log gives the mapped range for a pseudo-CT skull. With `io.save_property_maps`, both are written as `absorption_np` [Np/m at the source frequency] and `absorption_fraction` maps.

## Thermal Dose (CEM43) Calculation

PRESTUS computes **both** CEM43 formulations in every k-Wave thermal simulation. Both are written to the CSV output, exported as [NIfTI](https://nifti.nimh.nih.gov/) volumes and shown as separate figures. The HTML report (dashboard, safety verdict, methods text) uses the formulation selected by `thermal.cem43_iso`. With the BabelBrain thermal backend only BabelBrain's own dose is available (R = 0.25 below 43 °C, 0.5 at or above 43 °C): the ISO maps are NaN, and `load_parameters` sets `cem43_iso` to 0 with a warning, so the report uses BabelBrain's dose.

Both forms are based on the thermal isoeffective dose of [Sapareto & Dewey (1984)](https://doi.org/10.1016/0360-3016(84)90379-1), CEM43 = Σ Δt · R^(43 − T).

### k-Wave Form (Default)

CEM43 as implemented in k-Wave's `kWaveDiffusion`: R = 0.25 for 37 ≤ T < 43 °C, R = 0.5 for T ≥ 43 °C, and no dose below 37 °C.

```
cem43 = cem43 + ...
        dt ./ 60 .* ...
        (0.25 .* (T >= 37 & T < 43) + ...
         0.5  .* (T >= 43)).^(43 - T);
```

### ISO Form

Set `thermal.cem43_iso: 1` to use this form in the HTML report. It is computed by `cem43_iso_increment` and combines:

- **ISO 14708-3:2017**, Clause 17.1, Formula (1): R = 0.25 for T < 43 °C and R = 0.5 for T ≥ 43 °C, "valid for temperatures between 39 °C and 57 °C". Table 101 gives the tissue thresholds brain 2, bone 16 and skin 21 CEM43, which [ITRUSST](https://itrusst.com/) adopts ([Aubry et al., 2025](https://doi.org/10.1016/j.brs.2025.10.007)) and PRESTUS uses as reference levels (`get_risk_limits`).
- **ITRUSST biophysical safety consensus** (Aubry et al., 2025): dose is accumulated over the whole time the tissue is above 39 °C, including cooling; and because thermotolerance "is not maintained beyond 43 °C, nor does it return upon cooling", R = 0.5 is used for all time steps after 43 °C has been reached (citing Sapareto & Dewey, 1984). PRESTUS tracks this with T_max.
- **PRESTUS choice**: voxels that reach 57 °C (current T or T_max) are set to CEM43 = ∞. ISO 14708-3 only states that the formula is not valid above 57 °C; ∞ is a conservative flag for that range, not taken from either source.

```
cem43_iso = cem43_iso + ...
        dt ./ 60 .* ...
        (0    .* (T < 39  & T_max < 43) + ...
         0.25 .* (T >= 39 & T < 43 & T_max < 43) + ...
         0.5  .* (T >= 43 | T_max >= 43)).^(43 - T);
cem43_iso(T >= 57 | T_max >= 57) = Inf;
```

### Output Naming

| Suffix | Description |
|---|---|
| `CEM43` / `CEM43_end` | kWave form (max over protocol / at end of last pulse) |
| `CEM43_iso` / `CEM43_iso_end` | ISO form (max over protocol / at end of last pulse) |

Both suffixes apply to NIfTI filenames (`sub-XXX_<medium>_T1w<affix>_CEM43_iso.nii.gz`, etc.) and to CSV columns (`CEM43_brain`, `CEM43iso_brain`, etc.).

??? note "References"
    - Sapareto SA, Dewey WC (1984). Thermal dose determination in cancer therapy. *Int J Radiat Oncol Biol Phys* 10:787–800. https://doi.org/10.1016/0360-3016(84)90379-1
    - ISO 14708-3:2017. Implants for surgery. Active implantable medical devices. Part 3: Implantable neurostimulators. Clause 17.1, Formula (1), Table 101. https://www.iso.org/standard/60539.html
    - Aubry J-F, Attali D, et al. (2025). ITRUSST consensus on biophysical safety for transcranial ultrasound stimulation. *Brain Stimul* 18:1896–1905. https://doi.org/10.1016/j.brs.2025.10.007
    - Martin E, et al. (2024). ITRUSST consensus on standardised reporting for transcranial ultrasound stimulation. *Brain Stimul* 17:607–615 (Eq. 5: Sapareto & Dewey form only). https://doi.org/10.1016/j.brs.2024.04.013
