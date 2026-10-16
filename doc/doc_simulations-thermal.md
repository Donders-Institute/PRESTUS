# Thermal simulations

!!! note "At a glance"
    - **What:** How to describe the sonication protocol (pulses, trains, repetitions) and how heating and thermal dose are computed.
    - **When you need it:** Whenever you simulate heating (`modules.run_heating_sims: 1`).
    - **Before you start:** A working acoustic simulation, and the timing of your protocol.

Thermal simulations depend on the output of the acoustic simulations. As such, running thermal simulations without first running acoustic simulations is not an option.

Thermal simulations run in k-Wave by default. With `acoustic_backend: babelbrain`, BabelBrain's own Pennes solver can be used instead (`thermal_backend: babelbrain`). See [Simulation backends](doc_backend.md) for valid combinations and how the results compare, and [BabelBrain thermal](doc_babelbrain.md#thermal-simulation) for the timing mapping and tissue-property differences. **Known limitation:** BabelBrain thermal uses its own built-in tissue thermal properties, so the thermal fields in `medium_properties` apply only to the k-Wave thermal backend.

Protocol timing follows the nomenclature of the [TUS Calculator](https://www.socsci.ru.nl/fusinitiative/tuscalculator/). You can use this tool to specify the desired protocol. Note that all three levels of the nesting (pulse, pulse train, pulse train repetition) have to be specified. If you only want to simulate a pulse train without repetitions, set PTRI and PTRD to be equal to the PTD.

For protocols with long intervals between successive pulse trains, a coarser temporal sampling resolution may be sufficient. This can be specified with `post_pt_timestep`.

### Example: Timing Specification

Continuous pulse train protocols can be modelled either as (a) a single pulse train repetition with multiple pulse trains, or as (b) multiple pulse train repetitions of a single pulse train.

e.g., modeling a 20s, 5 Hz, 10% DC protocol as a single PTR (a):

- `pd` | 0.02 | Pulse Duration (PD) [seconds]
- `pri` | 0.2 | Pulse Repetition Interval (PRI) [seconds]
- `ptd` | 20 | Pulse Train Duration (PTD) [seconds]
- `pt_timestep` | 0.02 | Modeling time steps inside a PT [seconds]
- `ptri` | 20 | Pulse Train Repetition Interval (PTRI) [seconds]
- `ptrd` | 20 | Pulse Train Repetition Duration (PTRD) [seconds] 
- `post_ptri_dur` | 20 | Steady-state Duration following the PTRI [seconds]
- `post_pt_timestep` | 1 | Modeling time steps following PT & PTRI [seconds]
- `equal_step_duration` | 0 | [0] Variable step duration between pulse on & off [fixed step count = 1] [1] Variable step number between pulse on & off [fixed step duration: sim_time_steps]

You may also specify a protocol version that neglects the fine-grained heating dynamics between each pulse. 
In this scenario, you can model a single longer pulse train. Modeling a long PTR with a single Pulse - single PT has the added benefit of allowing a coarser estimation step size during the (collapsed) OFF period.

- `pd` | 2
- `pri` | 2
- `ptd` | 2
- `pt_timestep` | 0.02
- `ptri` | 20
- `ptrd` | 20
- `post_ptri_dur` | 0
- `post_pt_timestep` | 1

Heat dissipation continues beyond the final ultrasound pulse. To capture delayed heating, it is advised to run simulations until a steady-state has been reached. This can be specified via `timing.post_ptri_dur`. This acquires the corresponding duration following the end of the protocol. Time steps for the post-stimulation duration can be defined to be different from the time steps of the main simulation ```pt_timestep``` (set separately via `timing.post_pt_timestep`). Alternatively, for a continuous pulse train, a steady-state phase could be included as part of an extended PTR interval. The former option enables a steady-state phase also for protocols with PTRs.

### Additional parameters

These parameters influence the underlying mechanisms of how the thermal simulations are run.

- `pt_timestep` refer to the temporal step size for modeling pulses in each pulse train; this influences the length of the simulation.
- `equal_step_duration`
        - When `equal_step_duration` is enabled, the amount of simulated steps is calculated with the `pt_timestep` duration. For example, when using a `pt_timestep` duration of 10 ms, 6 steps of stimulation and 14 steps without stimulation are simulated 5 times for each `pt_timestep` (20\*5 = 100 steps).
        - When `equal_step_duration` is disabled, every on and off step will be simulated as one simulation step (2\*5 = 10 steps), reducing the computational load.
- `post_ptri_dur` refers to the duration of a post-stimulation period.
- `temp_0` refers to the starting temperature in degrees, and can be specified for each tissue separately.
- `sensor_xy_halfsize` defines the half-size (in grid units) of the rectangular XY sensor window centred on the transducer position. Only voxels within this window are included in the temperature sensor mask passed to `kWaveDiffusion`; voxels outside it are not recorded. The full Z extent of the grid is always included. Reducing this value lowers memory use during thermal simulation.
- `record_t_at_every_step` controls whether the full sensor window temperature is stored at every simulation time step. When disabled (default, `0`), the sensor is discarded before the simulation and only the final temperature field is retained, which is substantially more memory-efficient. Enable (`1`) only if you need the complete temperature time series, but note that this can cause `out of memory` errors for large grids or many time steps.

### Heat deposition

The acoustic field deposits heat at the rate Q = 2&alpha;<sub>abs</sub>I = &alpha;<sub>abs</sub>p<sup>2</sup>/(&rho;c), with the absorption coefficient &alpha;<sub>abs</sub> [Np/m] at the source frequency. `thermal.absorption_model` selects how &alpha;<sub>abs</sub> is specified per tissue (see [Medium properties](doc_medium.md)):

- `fraction` (default): `absorption_fraction` × the attenuation coefficient.
- `coefficient`: `absorption_coeff` [dB/(cm·MHz<sup>y</sup>)], with the tissue's attenuation power law y (`alpha_power`). Every tissue in `layers` needs a value; `config_default.yaml` provides them.

A run uses one model. A pseudo-CT skull with mapped attenuation has no matching absorption map, so `coefficient` falls back to `fraction` there, with a warning.

### Thermal dose (CEM43) calculation

PRESTUS always computes **both** CEM43 formulations in every k-Wave thermal simulation run. Both are written to the CSV output, exported as NIfTI volumes, and shown as separate figures. The HTML report (dashboard, safety verdict, methods text) uses the formulation selected by `parameters.thermal.cem43_iso`. With the BabelBrain thermal backend only BabelBrain's own dose is available (R = 0.25 below 43 °C, 0.5 at or above 43 °C): the ISO maps are NaN, and `load_parameters` sets `cem43_iso` to 0 with a warning, so the report uses BabelBrain's dose.

Both forms are based on the thermal isoeffective dose of Sapareto & Dewey (1984), CEM43 = Σ Δt · R^(43 − T).

#### kWave form (default)

CEM43 as implemented in k-Wave's `kWaveDiffusion`: R = 0.25 for 37 ≤ T < 43 °C, R = 0.5 for T ≥ 43 °C, and no dose below 37 °C.

```
cem43 = cem43 + ...
        dt ./ 60 .* ...
        (0.25 .* (T >= 37 & T < 43) + ...
         0.5  .* (T >= 43)).^(43 - T);
```

#### ISO form

Set `parameters.thermal.cem43_iso = 1` to use this form in the HTML report. It is computed by `cem43_iso_increment` and combines:

- **ISO 14708-3:2017**, Clause 17.1, Formula (1): R = 0.25 for T < 43 °C and R = 0.5 for T ≥ 43 °C, "valid for temperatures between 39 °C and 57 °C". Table 101 gives the tissue thresholds brain 2, bone 16 and skin 21 CEM43, which ITRUSST adopts (Aubry et al., 2025) and PRESTUS uses as reference levels (`get_risk_limits`).
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

#### References

- Sapareto SA, Dewey WC (1984). Thermal dose determination in cancer therapy. *Int J Radiat Oncol Biol Phys* 10:787–800. https://doi.org/10.1016/0360-3016(84)90379-1
- ISO 14708-3:2017. Implants for surgery — Active implantable medical devices — Part 3: Implantable neurostimulators. Clause 17.1, Formula (1), Table 101. https://www.iso.org/standard/60539.html
- Aubry J-F, Attali D, et al. (2025). ITRUSST consensus on biophysical safety for transcranial ultrasound stimulation. *Brain Stimul* 18:1896–1905. https://doi.org/10.1016/j.brs.2025.10.007
- Martin E, et al. (2024). ITRUSST consensus on standardised reporting for transcranial ultrasound stimulation. *Brain Stimul* 17:607–615 (Eq. 5: Sapareto & Dewey form only). https://doi.org/10.1016/j.brs.2024.04.013

#### Output naming

| Suffix | Description |
|---|---|
| `CEM43` / `CEM43_end` | kWave form (max over protocol / at end of last pulse) |
| `CEM43_iso` / `CEM43_iso_end` | ISO form (max over protocol / at end of last pulse) |

Both suffixes apply to NIfTI filenames (`sub-XXX_final_CEM43_iso_orig_coord.nii.gz`, etc.) and to CSV columns (`CEM43_brain`, `CEM43iso_brain`, etc.).