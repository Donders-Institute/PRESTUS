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
- `sensory_xy_halfsize` refers to the window along the transducer axis in which temperature changed are recorded.
- `record_t_at_every_step` can be enabled if there is a need to record temperature changes at every time step, but can lead to `out of memory` errors.

### Thermal dose (CEM43) calculation

There are two implementations of thermal dose calculations:

- By default, PRESTUS calculates CEM43 as implemented in k-Wave's kWaveDiffusion.

```
cem43 = cem43 + ...
        dt ./ 60 .* ...
        (0.25 .* (T >= 37 & T < 43 ) + ...
        0.5 .* (T >= 43)).^(43 - T);
```

- By setting ```paramaters.thermal.cem43_iso``` to 1, PRESTUS calculates CEM43 according to ISO standards. This variant defines CEM43 as zero below 39 degrees Celsius, sets it to infinite above 57 degrees, and applies a scaling of 0.5 to all voxels that have crossed 43 degrees in the past (not just following the most recent energy update).

```
 cem43_iso = cem43_iso + ...
        dt ./ 60 .* ...
        (0 .* (T < 39 & T_max < 43) + ...
        0.25 .* (T >= 39 & T < 43 & T_max < 43) + ...
        0.5 .* (T >= 43 | T_max >= 43) + ...
        Inf .* (T >= 57)).^(43 - T);
```