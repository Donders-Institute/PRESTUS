# Tips on running the thermal simulations
- The thermal simulations themselves are based on the output of the acoustic simulations.
	- So running the thermal simulations without any acoustic output to base it on is not an option.

# Parameters explanation
## Parameters that have to be changed for each experiment
- duty_cycle has to be a value between 0 and 1, and represents the percentage of time of the trial during which stimulation is administered.
- iti is the inter-trial-interval given in seconds
- n_trials is representative of the aboumt of trials that your experiment will use.
	- Controls the amount of times the wrapper has to loop through the k-wave thermal simulation.
- stim_duration refers to the length of time in seconds in a trial during which stimulation is adminstered.

## Parameters that should be changed with caution
These parameters influence the underlying mechanisms of how the thermal simulations are run.

- sim_time_steps refer to the duration of each simulated step, so the length of the loop.
- on_off_step_duration is the duration of each duty cycle.
- equal_steps, when enabled, will concatenate different steps to reach a simulation duration similar to your duty cycle. If disabled, it will try to compute the step duration necessary for your given duty cycle.
- post_stim_time_step_dur refers to the amount of steps in between two trials.

Not meeting the following assumptions can prevent the thermal simulations from running:
- on_off_step_duration * duty_cycle / sim_time_steps should produce an integer.
- on_off_step_duration * (1 - duty_cycle) / sim_time_steps should produce an integer.
- cycle_duration / on_off_step_duration should produce an integer.

## Optional parameters
- temp_0 refers to the starting temperature in degrees.
- sensory_xy_halfsize refers to the window along the transducer axis in which temperature changed are recorded.
- record_t_at_every_step can be enabled if there is a need to record temperature changes at every timestep, but can be very memory intensive.
