function params_thermal = thermal_parameters(parameters)
% THERMAL_PARAMETERS Computes all thermal simulation step parameters.
%
% Input:
%   parameters - Struct with .thermal.pd, .pri, .ptd, .ptri, .ptrd, 
%                .pt_timestep, .post_pt_timestep, .post_ptri_dur, .equal_step_duration
%
% Output:
%   params_thermal - Struct with:
%     .on_steps_n, .on_steps_dur           % Pulse ON within PT
%     .off_steps_n, .off_steps_dur         % Pulse OFF within PT  
%     .n_pulses_per_pt                     % Pulses in one PT
%     .ptri_off_steps_n, .ptri_off_step_dur % PTRI-off period (coarse timestep)
%     .post_ptri_steps_n, .post_ptri_step_dur % Post whole-PTRD cooloff
%     .n_ptri_reps                         % PTRI reps in PTRD

arguments
    parameters struct
end

disp('Checking thermal parameters...');

thermal = parameters.thermal;

% === Pulse / Pulse Train (fine timestep: pt_timestep) ===

pri = thermal.pri;
pd  = thermal.pd;
ptd = thermal.ptd;
pt_dt = thermal.pt_timestep;
dc = pd / pri;  % Duty cycle (computed)

% encode the above in params_thermal
params_thermal.pri = pri;
params_thermal.pd = pd;
params_thermal.ptd = ptd;
params_thermal.pt_dt = pt_dt;
params_thermal.dc = dc;

if thermal.equal_step_duration == 0
    % Different durations, fixed steps=1
    params_thermal.pt_on_steps_n   = 1;
    params_thermal.pt_on_steps_dur = pd;
    params_thermal.pt_off_steps_n  = 1;
    params_thermal.pt_off_steps_dur= pri - pd;
    
    fprintf('PT cycle: 1 ON (%.3fs) + 1 OFF (%.3fs) [equal_step_duration=0]\n', ...
        params_thermal.pt_on_steps_dur, params_thermal.pt_off_steps_dur);
else
    % Equal duration, variable steps
    params_thermal.pt_on_steps_dur  = pt_dt;
    params_thermal.pt_off_steps_dur = pt_dt;
    
    params_thermal.pt_on_steps_n  = round_if_integer(pd  / pt_dt, 'ON steps must be integer');
    params_thermal.pt_off_steps_n = round_if_integer((pri-pd) / pt_dt, 'OFF steps must be integer');
    
    fprintf('PT cycle: %.0f ON + %.0f OFF steps @ %.3fs each [equal_step_duration=1]\n', ...
        params_thermal.pt_on_steps_n, params_thermal.pt_off_steps_n, pt_dt);
end

% Pulses per PT
params_thermal.n_pulses_per_pt = round_if_integer(ptd / pri, 'Pulses per PT must be integer');
fprintf('PTD=%.2fs contains %.0f pulses @ PRI=%.3fs\n', ptd, params_thermal.n_pulses_per_pt, pri);

% === PTRI repetition ===

params_thermal.ptrd        = thermal.ptrd;
params_thermal.ptri        = thermal.ptri;
params_thermal.ptri_off    = thermal.ptri - ptd;

params_thermal.n_ptri_reps = round_if_integer(params_thermal.ptrd / params_thermal.ptri, ...
                                             'PTRI reps in PTRD must be integer');

% === PTRI-off & post-PTRI (coarse timestep: post_pt_timestep) ===

post_dt = thermal.post_pt_timestep;

% PTRI-OFF (between PTs)
if params_thermal.ptri_off <= 0
    params_thermal.ptri_off_steps_n  = 0;
    params_thermal.ptri_off_step_dur = 0;
else
    params_thermal.ptri_off_step_dur = post_dt;
    params_thermal.ptri_off_steps_n  = round_if_integer(params_thermal.ptri_off / post_dt, ...
                                                       'PTRI-OFF steps must be integer');
    fprintf('PTRI-OFF=%.2fs → %.0f steps @ %.2fs\n', params_thermal.ptri_off, ...
        params_thermal.ptri_off_steps_n, post_dt);
end

% Post PTRD cooloff
post_dur = thermal.post_ptri_dur;
params_thermal.post_ptri_dur = post_dur;
if post_dur <= 0
    params_thermal.post_ptri_steps_n  = 0;
    params_thermal.post_ptri_step_dur = 0;
    disp('No post-PTRD modeling requested...');
else
    params_thermal.post_ptri_step_dur = post_dt;
    params_thermal.post_ptri_steps_n  = round_if_integer(post_dur / post_dt, ...
                                                        'Post-PTRD steps must be integer');
    fprintf('Post-PTRD=%.1fs → %.0f steps @ %.1fs\n', post_dur, ...
        params_thermal.post_ptri_steps_n, post_dt);
end

% Store derived scalars for convenience (optional)
params_thermal.dc  = dc;
params_thermal.prf = 1 / pri;

end