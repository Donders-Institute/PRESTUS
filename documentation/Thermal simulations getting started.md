# Tips on running the thermal simulations

- The thermal simulations themselves are based on the output of the acoustic simulations.
- So running the thermal simulations without any acoustic output to base it on is not an option.

Before feeding in parameters for the heating simulations it is important to note that these simulations don't model the individual pulses and their pulse-repetition frequencies (PRF) but concatenate the pulses over a trail into one block of stimulation and one block without.

# Examples

## Online study

Let's take a pulse-repetition frequency of 5Hz, a stimulation duration (`stim_duration` in a config file) of 1s, an inter-trial-interval of 15s (`iti`), 30 trials (`n_trials`) and a duty cycle of 30% (`duty_cycle`, 0.3).

When running the actual stimulation, these are the parameters of the stimulation sequence the transducer will produce: - A pulse-repetition frequency of 5Hz will result in 5 pulse repetition periods of 200ms (1s / 5Hz). - Within each pulse repetition period, the pulse parameters are determined by the duty cycle (0.3), the transducer will be 'on' for 60ms (0.3 \* 200ms) of stimulation, and 'off' for 140ms (0.7 \* 200ms) without stimulation. - Stimulation will occur 5 times (5 pulse repetition periods) within your stimulation duration of 1s. - After this, no stimulation will be administered for 14s. - This cycle will be repeated 30 times.

For the simulations, the number of recorded temperature readings is reduced to limit the computational load. This in an area where a compromise in temporal resolution can be made since the difference between having 10 or 100 cooling periods in a second (for example) is negligible for our application. Since the temperature is already recorded during every `on_off_repetition` cycle, the duty cycle will not be taken into account as well. So instead, all pulse repetition periods within a `stim_duration` are combined into a 300ms period of stimulation (60ms \* 5Hz) and 700ms without stimulation (140ms \* 5Hz) within each trial. This is why the pulse-repetition frequency does not have to be included in the thermal parameters.

- `duty_cycle` = 0.3
- `stim_duration` = 1
- `n_trials` = 30
- `iti` = 15s
- `continuous_stimulation` = 0

## Offline study

For this, you would use a continuous stimulation protocol that should be simulated as a single trial. Let's say you want to simulate an offline protocol with a duty cycle of 30% (`duty_cycle`, 0.3) and a total stimulation duration of 120s (`stim_duration`). You would use the following parameters:

- `duty_cycle` = 0.3
- `stim_duration` = 120
- `n_trials` = 1
- `continuous_stimulation` = 1

\*Again, the PRF is not taken into account since the cooling effects between these cycles are limited.

# Parameters explained

## Parameters that have to be changed for each experiment

- `duty_cycle` has to be a value between 0 and 1, and represents the percentage of time within a the trial during which stimulation is administered.
- `iti` is the inter-trial-interval given in seconds. Should always be higher than the `stim_duration` since `iti - stim_duration = post_stim_period`
- `n_trials` is representative of the amount of trials that your experiment will use.
        - Controls the amount of times the wrapper has to loop through the k-wave thermal simulation.
- `stim_duration` refers to the length of time in seconds in a trial during which stimulation is administered.
- `continuous_protocol` is a binary option that can be enabled to run continuous stimulation, mostly used in offline protocols. If this is enabled, the `iti` is not used anymore.

## Parameters that should be changed with caution

These parameters influence the underlying mechanisms of how the thermal simulations are run.

- `sim_time_steps` refer to the duration of each simulated step, so the length of the loop.

- `on_off_step_duration` is the duration of each on+off cycle [s]. For example, with a `stim_duration` of 120s and an `on_off_step_duration` of 0.01s, 1200 steps will be simulated. Can also be left blank if `duty_cycle` and `sim_time_steps` are specified. One issue that can result from this is that the on_off_step_duration is too small, resulting in the GPU running out of RAM. This can be resolved by manually setting a higher value.

- `equal_steps`
        - When `equal_steps` is enabled, the amount of simulated steps is calculated with the `sim_time_steps` duration. For example, when using a `sim_time_steps` duration of 10 ms, 6 steps of stimulation and 14 steps without stimulation are simulated 5 times for each `stim_duration` (20\*5 = 100 steps).
        - When `equal_steps` is disabled, every on and off step will be simulated as one simulation step (2\*5 = 10 steps), reducing the computational load.

- `post_stim_time_step_dur` refers to the duration of each step in the `post_stim_period`.

Not meeting the following assumptions can prevent the thermal simulations from running:
- `on_off_step_duration * duty_cycle / sim_time_steps` should produce an integer.
- `on_off_step_duration * (1 - duty_cycle) / sim_time_steps` should produce an integer.
- `cycle_duration / on_off_step_duration` should produce an integer.
- `cycle_duration` is the product of `stim_duration / on_off_step_duration`.

## Optional parameters

- `temp_0` refers to the starting temperature in degrees.
- `sensory_xy_halfsize` refers to the window along the transducer axis in which temperature changed are recorded.
- `record_t_at_every_step` can be enabled if there is a need to record temperature changes at every time step, but can lead to `out of memory` errors.
