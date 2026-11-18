function eta = thermal_ramp_factors(params_thermal, sched)
% THERMAL_RAMP_FACTORS  Heat factor of a ramped pulse for each heat source
%
% With timing.pulse_ramp_dur > 0, every pulse rises and falls over that time.
% Heat scales with the squared amplitude, so each pulse of a source delivers
% COMPUTE_RAMP_FACTOR times the heat of a rectangular pulse, the same factor
% as Ipa (ACOUSTIC_ANALYSIS). The factor depends on the pulse duration, which
% is per source: the pd of the transducer whose windows the source follows
% (THERMAL_SCHEDULE SourceTx; transducer(i).timing.pd with async_schedule
% 'custom'). THERMAL_SIMULATION multiplies each source's heat by its factor.
%
% Use as:
%   eta = thermal_ramp_factors(params_thermal, sched)
%
% Input:
%   params_thermal - struct from THERMAL_PARAMETERS (ramp_dur [s], ramp_shape)
%   sched          - struct from THERMAL_SCHEDULE (source(i).pd [s])
%
% Output:
%   eta - [1 x n_sources] factors in (0, 1]; all 1 without a ramp
%
% See also: COMPUTE_RAMP_FACTOR, THERMAL_PARAMETERS, THERMAL_SCHEDULE, THERMAL_SIMULATION

n = numel(sched.source);
eta = ones(1, n);
if ~isfield(params_thermal, 'ramp_dur') || ~(params_thermal.ramp_dur > 0)
    return
end
for i = 1:n
    pd_i = sched.source(i).pd;
    if 2 * params_thermal.ramp_dur > pd_i * (1 + 1e-9)
        error('PRESTUS:thermal:ramp', ...
            'timing.pulse_ramp_dur (%g s) must satisfy 2*pulse_ramp_dur <= pd (%g s) of heat source %d.', ...
            params_thermal.ramp_dur, pd_i, i);
    end
    eta(i) = compute_ramp_factor(pd_i, params_thermal.ramp_dur, params_thermal.ramp_shape);
end
end
