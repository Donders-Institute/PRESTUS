function params_thermal = thermal_parameters(parameters, silent)
% THERMAL_PARAMETERS  Compute time-stepping parameters for a multi-level thermal sonication protocol
%
% Discretizes a hierarchical pulse protocol into integer time steps at two
% temporal resolutions: fine steps (pt_timestep) covering each pulse ON
% and OFF interval within a pulse train (PT), and coarse steps
% (post_pt_timestep) covering the inter-train OFF periods (PTRI-OFF) and
% the optional post-PTRD steady-state cool-off. All step counts are
% verified to be integers (error if not). The hierarchy is:
%   pulse (pd / pri) -> pulse train (ptd) -> PTRI reps (ptrd) -> post-PTRD
%
% Use as:
%   params_thermal = thermal_parameters(parameters)
%   params_thermal = thermal_parameters(parameters, silent)
%
% Input:
%   parameters - PRESTUS config; must contain timing sub-struct:
%                timing.pd [s], timing.pri [s], timing.ptd [s],
%                timing.pt_timestep [s], timing.ptri [s], timing.ptrd [s],
%                timing.post_pt_timestep [s], timing.post_ptri_dur [s],
%                timing.equal_step_duration
%   silent     - suppress console summary printout (optional, default: false)
%
% Output:
%   params_thermal - struct with fields:
%                    pri, pd, ptd, pt_dt [s], dc (duty cycle), prf [Hz],
%                    pt_on_steps_n, pt_on_steps_dur [s],
%                    pt_off_steps_n, pt_off_steps_dur [s],
%                    n_pulses_per_pt, ptri_off [s], ptri_off_steps_n,
%                    n_ptri_reps, post_ptri_steps_n, post_ptri_step_dur [s]
%
% See also: THERMAL_SIMULATION, THERMAL_PLOT_PROTOCOL

arguments
    parameters (1,1) struct
    silent     (1,1) logical = false
end

if ~silent
    disp('Implementing thermal protocol...');
end

thermal = parameters.thermal;
timing  = parameters.timing;

% === Pulse / Pulse Train (fine timestep: pt_timestep) ===

pri = timing.pri;
pd  = timing.pd;
ptd = timing.ptd;
pt_dt = timing.pt_timestep;
dc = pd / pri;  % Duty cycle (computed)

% encode the above in params_thermal
params_thermal.pri = pri;
params_thermal.pd = pd;
params_thermal.ptd = ptd;
params_thermal.pt_dt = pt_dt;
params_thermal.dc = dc;

if timing.equal_step_duration == 0
    % Different durations, fixed steps=1
    params_thermal.pt_on_steps_n   = 1;
    params_thermal.pt_on_steps_dur = pd;
    params_thermal.pt_off_steps_n  = 1;
    params_thermal.pt_off_steps_dur= pri - pd;
else
    % Equal duration, variable steps
    params_thermal.pt_on_steps_dur  = pt_dt;
    params_thermal.pt_off_steps_dur = pt_dt;
    
    params_thermal.pt_on_steps_n  = round_if_integer(pd  / pt_dt, 'ON steps must be integer');
    params_thermal.pt_off_steps_n = round_if_integer((pri-pd) / pt_dt, 'OFF steps must be integer');
end

% Pulses per PT
params_thermal.n_pulses_per_pt = round_if_integer(ptd / pri, 'Pulses per PT must be integer');

% === PTRI repetition ===

params_thermal.ptrd        = timing.ptrd;
params_thermal.ptri        = timing.ptri;
params_thermal.ptri_off    = timing.ptri - ptd;

params_thermal.n_ptri_reps = round_if_integer(params_thermal.ptrd / params_thermal.ptri, ...
                                             'PTRI reps in PTRD must be integer');

% === PTRI-off & post-PTRI (coarse timestep: post_pt_timestep) ===

params_thermal.post_pt_timestep = timing.post_pt_timestep;
post_dt = timing.post_pt_timestep;

% PTRI-OFF (between PTs)
if params_thermal.ptri_off <= 0
    params_thermal.ptri_off_steps_n  = 0;
    params_thermal.ptri_off_step_dur = 0;
else
    params_thermal.ptri_off_step_dur = post_dt;
    params_thermal.ptri_off_steps_n  = round_if_integer(params_thermal.ptri_off / post_dt, ...
                                                       'PTRI-OFF steps must be integer');
end

% Post PTRD steady-state
post_dur = timing.post_ptri_dur;
params_thermal.post_ptri_dur = post_dur;
if post_dur <= 0
    params_thermal.post_ptri_steps_n  = 0;
    params_thermal.post_ptri_step_dur = 0;
else
    params_thermal.post_ptri_step_dur = post_dt;
    params_thermal.post_ptri_steps_n  = round_if_integer(post_dur / post_dt, ...
                                                        'Post-PTRD steps must be integer');
end

% Store derived scalars for convenience (optional)
params_thermal.dc  = dc;
params_thermal.prf = 1 / pri;

%% Protocol Summary

if ~silent
    fprintf('\n========================================\n');
    fprintf('THERMAL PROTOCOL \n');
    fprintf('========================================\n');
    fprintf('Pulse:        PD=%.3fs, PRI=%.3fs (DC=%.1f%%, PRF=%.1fHz)\n', ...
        params_thermal.pd, params_thermal.pri, params_thermal.dc*100, params_thermal.prf);
    
    fprintf('Pulse Train:  PTD=%.2fs = %d pulses (%d ON + %d OFF steps @ %.2g s)\n', ...
        params_thermal.ptd, params_thermal.n_pulses_per_pt, ...
        params_thermal.pt_on_steps_n, params_thermal.pt_off_steps_n, params_thermal.pt_dt);
    
    fprintf('PTRI:         PTRI=%.1fs (%d reps/PTD, PTRI-OFF=%.1fs = %d coarse steps @ %.2fs)\n', ...
        params_thermal.ptri, params_thermal.n_ptri_reps, params_thermal.ptri_off, ...
        params_thermal.ptri_off_steps_n, params_thermal.post_pt_timestep);
    
    fprintf('PTRD:         PTRD=%.1fs s total sonication\n', params_thermal.ptrd);
    
    if params_thermal.post_ptri_steps_n > 0
        fprintf('Steady-state:     %.0fs post-PTRD (%d steps @ %.1fs)\n', ...
            params_thermal.post_ptri_dur, params_thermal.post_ptri_steps_n, params_thermal.post_ptri_step_dur);
    else
        fprintf('Steady-state:     None requested\n');
    end
    
    total_fine_steps = params_thermal.n_pulses_per_pt * ...
        (params_thermal.pt_on_steps_n + params_thermal.pt_off_steps_n) * params_thermal.n_ptri_reps;
    total_coarse_steps = params_thermal.ptri_off_steps_n * params_thermal.n_ptri_reps + params_thermal.post_ptri_steps_n;
    
    fprintf('Simulation steps:  %d fine + %d coarse = %d total\n', total_fine_steps, total_coarse_steps, total_fine_steps + total_coarse_steps);
    fprintf('========================================\n\n');
end
end