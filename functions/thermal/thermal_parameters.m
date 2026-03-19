function params_thermal = thermal_parameters(parameters, silent)
% THERMAL_PARAMETERS computes time stepping parameters for multi-level thermal ultrasound protocol.
%
% This function discretizes a hierarchical thermal sonication protocol into integer 
% time steps at two temporal resolutions: fine steps for pulse trains (pt_timestep) 
% and coarse steps for inter-train OFF periods (post_pt_timestep). The protocol 
% consists of pulses → pulse trains (PT) → pulse-train repetitions (PTRI within PTRD) 
% → optional post-PTRD steady-state.
%
% Use as:
%   params_thermal = thermal_parameters(cfg)
%
% Required input (in parameters.timing):
%   .pd              - pulse duration [s]
%   .pri             - pulse repetition interval [s]
%   .ptd             - pulse train duration [s]
%   .pt_timestep     - fine timestep for PT [s]
%   .ptri            - pulse train repetition interval [s]
%   .ptrd            - pulse train repetition duration [s]
%   .post_pt_timestep- coarse timestep for OFF periods [s]
%   .post_ptri_dur   - post-PTRD steady-state duration [s]
%   .equal_step_duration - 0: fixed 1-step ON/OFF; 1: equal dt steps
%
% Output structure params_thermal contains:
%   .pri, .pd, .ptd, .pt_dt              - input values
%   .dc                                  - duty cycle = pd/pri
%   .prf                                 - pulse repetition frequency = 1/pri [Hz]
%   .pt_on_steps_n, .pt_on_steps_dur     - pulse ON discretization within PT
%   .pt_off_steps_n, .pt_off_steps_dur   - pulse OFF discretization within PT  
%   .n_pulses_per_pt                     - integer pulses per pulse train
%   .ptri_off, .ptri_off_steps_n         - PTRI-OFF duration & coarse steps
%   .n_ptri_reps                         - integer PTRI repetitions in PTRD
%   .post_ptri_steps_n, .post_ptri_step_dur - post-PTRD cool-off discretization

arguments
    parameters struct
    silent (1,1) logical = false  % Default: show summary
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