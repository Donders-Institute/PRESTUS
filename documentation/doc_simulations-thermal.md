# Thermal simulations

Thermal simulations depend on the output of the acoustic simulations. As such, running thermal simulations without first running acoustic simulations is not an option.

Protocol timing is oriented after the nomencalture in the [TUS Calculator](https://www.socsci.ru.nl/fusinitiative/tuscalculator/). You can use this tool to specify the desired protocol. Note that all three levels of the nesting (pulse, pule train, pulse train, pulse train repetition) have to be specified. If you only want to simulate a pulse train without repetitions, set PTRI and PTRD to be equal to the PTD.

For protocols with long intervals between successive pulse trains, a coarser temporal sampling resolution may be suffient. This can be specified with `post_pt_timestep`.

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

Heat dissipation continues beyond the final ultrasound pulse. To capture delayed heating, it is advised to run simulations until a steady-state has been reached. This can be specified via ```parameters.thermal.post_ptri_dur```. This acquires the corresponding duration following the end of the protocol. Time steps for the post-stimulation duration can be defined to be different from the time steps of the main simulation ```pt_timestep``` (by default they are equal, unless explicitly specified via ```parameters.thermal.post_pt_timestep```). Alternatively, for a continuous pulse train, a steady-state phase could be included as part of an extended PTR interval. The former option enables a steady-state phase also for protocols with PTRs.

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

### Thermal dose (CEM43) calculation

PRESTUS always computes **both** CEM43 formulations in every thermal simulation run. Both are written to the CSV output, exported as NIfTI volumes, and shown as separate figures. The HTML report displays the formulation selected by `parameters.thermal.cem43_iso`.

#### kWave form (default)

CEM43 as implemented in k-Wave's `kWaveDiffusion`:

```
cem43 = cem43 + ...
        dt ./ 60 .* ...
        (0.25 .* (T >= 37 & T < 43) + ...
         0.5  .* (T >= 43)).^(43 - T);
```

#### ISO form

Set `parameters.thermal.cem43_iso = 1` to display the ISO variant in the HTML report. The ISO form defines CEM43 as zero below 39 °C, infinite above 57 °C, and applies R = 0.5 to all voxels that have previously exceeded 43 °C (tracked via T_max), not just to those currently above 43 °C.

```
cem43_iso = cem43_iso + ...
        dt ./ 60 .* ...
        (0    .* (T < 39  & T_max < 43) + ...
         0.25 .* (T >= 39 & T < 43 & T_max < 43) + ...
         0.5  .* (T >= 43 | T_max >= 43) + ...
         Inf  .* (T >= 57)).^(43 - T);
```

#### Output naming

| Suffix | Description |
|---|---|
| `CEM43` / `CEM43_end` | kWave form (max over protocol / at end of last pulse) |
| `CEM43_iso` / `CEM43_iso_end` | ISO form (max over protocol / at end of last pulse) |

Both suffixes apply to NIfTI filenames (`sub-XXX_final_CEM43_iso_orig_coord.nii.gz`, etc.) and to CSV columns (`CEM43_brain`, `CEM43iso_brain`, etc.).