# Tips on running the thermal simulations

Thermal simulations depend on the output of the acoustic simulations. As such, running thermal simulations without first running acoustic simulations is not an option.

# Examples

## Continuous pulse protocols

For a continuous pulse protocol, the default option is to model each pulse train individually. 

- `n_trials`: 400               | number of on-off pulses/trials in the k-wave thermal simulation
- `duty_cycle`: 0.2             | percentage of pulse duration within a PRI [0-1]
- `pri_duration`: 0.2           | duration of a pulse repetition interval [s]
- `sim_time_steps`: 0.04        | simulation time steps during the stimulation period [s]
- `post_stim_dur`: 10           | duration of post-stimulation diffusion [s]
- `equal_steps`: 0              | [0] Variable step duration between pulse on & off [fixed step count = 1] [1] Variable step number between pulse on & off [fixed step duration: sim_time_steps]

Note that you may also specify a version of your protocol that neglects the fine-grained heating dynamics between each pulse, as cooling effects during the off-cycle are limited. Here, you specify the protocol via a single virtual trial that assumes concatenated stim-on followed by concatenated stim-off. For example, you could model a protocol with a duty cycle of 20% (`duty_cycle`, 0.2) and a total stimulation duration of 80s, as follows (with a longer simulation time step):

- `n_trials`: 1
- `duty_cycle`: 0.2
- `pri_duration`: 80
- `sim_time_steps`: 0.1

Heat dissipation continues beyond the final ultrasound pulse. To capture delayed heating, it is advised to run simulations until a steady-state has been reached. This can be specified in two ways: (1) Specify ```parameters.thermal.post_stim_dur```. This acquires the corresponding number of trials at the end of the protocol (note: but not following each trial). (2) Specify trials beyond  stimulation, and set the stimulation to off by implementing a pause. The latter can be controlled via ```parameters.start_break_trials``` and ```parameters.stop_break_trials```. The latter specification is more flexible, as it can be used to emulate a pulse train sequence (or a longer break in the protocol).

## Trial-wise / nested pulse train protocols

The following is a possible trial-by-trial design.

- 30 pulse train intervals [**missing**]
- 15 s pulse train repetition interval [30 x 15 = 450 s protocol] [**missing**]
- 1 s pulse train duration
- 5 Hz pulse repetition frequency: 5 x 200 ms pulse repetition interval in 1 s pulse train
- 30 % duty cycle (`duty_cycle`: 0.3) | Within each 200 ms pulse repetition interval, the transducer will be 'on' for 60 ms (0.3 \* 200ms) of stimulation, and 'off' for 140 ms (0.7 \* 200ms) without stimulation. 

**Note**: PRESTUS currently does not support nested modeling of individual pulse train pulses. This protocol can be approximated by collapsing the timing acoss pulses and pulse trains, however.

- 30 (`n_trials`: 30) x 15 s pulse train repetition intervals (`pri_duration`: 15) = 450 s
- Duty cycle: 0.3 % (pulse DC) * 1 s (= 300 ms) each 15 s (DC = 0.3/15=0.02 = 2%) (`duty_cycle`: 0.02)
- 5 Hz pulse repetition frequency: 200 ms pulse repetition interval [not modeled individually]

# Parameters explained

## Parameters that should be changed with caution

These parameters influence the underlying mechanisms of how the thermal simulations are run.

- `sim_time_steps` refer to the duration of each simulated step, so the length of the loop.
- `equal_steps`
        - When `equal_steps` is enabled, the amount of simulated steps is calculated with the `sim_time_steps` duration. For example, when using a `sim_time_steps` duration of 10 ms, 6 steps of stimulation and 14 steps without stimulation are simulated 5 times for each `stim_duration` (20\*5 = 100 steps).
        - When `equal_steps` is disabled, every on and off step will be simulated as one simulation step (2\*5 = 10 steps), reducing the computational load.

Not meeting the following assumptions can prevent the thermal simulations from running:
- `pri_duration * duty_cycle / sim_time_steps` should produce an integer.
- `pri_duration * (1 - duty_cycle) / sim_time_steps` should produce an integer.
- `cycle_duration / pri_duration` should produce an integer.
- `cycle_duration` is the product of `stim_duration / pri_duration`.

## Optional parameters

- `post_stim_dur` refers to the duration of a post-stimulation period.
- `temp_0` refers to the starting temperature in degrees, and can be specified for each tissue separately.
- `sensory_xy_halfsize` refers to the window along the transducer axis in which temperature changed are recorded.
- `record_t_at_every_step` can be enabled if there is a need to record temperature changes at every time step, but can lead to `out of memory` errors.