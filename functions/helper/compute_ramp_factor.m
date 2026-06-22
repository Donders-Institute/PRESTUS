function f = compute_ramp_factor(pd, ramp_dur, shape)
% COMPUTE_RAMP_FACTOR  Pulse-average / temporal-peak intensity ratio for a ramped pulse
%
% Returns the factor f = ISPPA / ISPTP: the ratio of pulse-average intensity to
% temporal-peak intensity when a symmetric amplitude ramp is applied at the start
% and end of each pulse. Intensity scales as amplitude^2, so the integral over
% each ramp lobe is k_shape * ramp_dur rather than ramp_dur.
%
%   k_shape = 1/3  for 'linear'  (integral of (t/T)^2 = 1/3)
%   k_shape = 3/8  for 'tukey'   (integral of [0.5*(1-cos(pi*t/T))]^2 = 3/8)
%
% Use as:
%   f = compute_ramp_factor(pd, ramp_dur, shape)
%
% Input:
%   pd       - pulse duration [s]
%   ramp_dur - rise (= fall) ramp duration [s]; must satisfy 2*ramp_dur <= pd
%   shape    - 'linear' or 'tukey'
%
% Output:
%   f        - ISPPA/ISPTP factor in (0, 1]; = 1 when ramp_dur = 0

arguments
    pd       (1,1) double {mustBePositive}
    ramp_dur (1,1) double {mustBeNonnegative} = 0
    shape    (1,1) string {mustBeMember(shape, ["linear","tukey"])} = "linear"
end

if ramp_dur == 0
    f = 1;
    return
end

if 2 * ramp_dur > pd
    error('compute_ramp_factor:rampExceedsPulse', ...
        'Total ramp duration (%.4g s) exceeds pulse duration (%.4g s).', 2*ramp_dur, pd);
end

switch shape
    case "linear"
        k = 1/3;
    case "tukey"
        k = 3/8;
end

% f = [k*t_rise + t_flat + k*t_fall] / pd   (rise = fall = ramp_dur)
t_flat = pd - 2 * ramp_dur;
f = (k * ramp_dur + t_flat + k * ramp_dur) / pd;

end
