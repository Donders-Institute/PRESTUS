function bt = babelbrain_thermal_timing(timing)
% BABELBRAIN_THERMAL_TIMING  Convert PRESTUS timing to BabelBrain thermal parameters
%
% PRESTUS timing hierarchy (see thermal_parameters.m): pd/pri (single pulse
% within a train) -> ptd (ON duration of ONE pulse train) -> ptri (period of
% one pulse-train repetition, ON+OFF) -> ptrd (TOTAL span across ALL
% repetitions, = n_ptri_reps * ptri) -> post_ptri_dur (extra cooldown after
% the whole repeated sequence).
%
% BabelBrain's BHTE (CalculateTemperatureEffects.RunBHTECycles) runs
% NumberGroupedSonications x Repetitions cycles of DurationUS on / DurationOff
% off, and after every Repetitions cycles an extra PauseBetweenGroupedSonications
% off period. With NumberGroupedSonications = 1 that pause runs exactly once,
% after the last repetition's DurationOff -- i.e. at t = ptrd, which is where
% the native thermal path starts its post_ptri_dur cooldown.
%
% Use as:
%   bt = babelbrain_thermal_timing(parameters.timing)
%
% Input:
%   timing - PRESTUS timing struct; requires pd, pri, ptd; uses ptri,
%            n_ptri_reps and post_ptri_dur when present
%
% Output:
%   bt - struct with fields
%     .duty_cycle     pd / pri
%     .prf_hz         1 / pri                       [Hz]
%     .duration_s     ptd  (BabelBrain DurationUS)  [s]
%     .duration_off_s ptri - ptd, or 0 (DurationOff) [s]
%     .repetitions    n_ptri_reps, or 1 (Repetitions)
%     .pause_after_s  post_ptri_dur, or 0 (PauseBetweenGroupedSonications) [s]
%
% See also: BABELBRAIN_THERMAL_WRAPPER, THERMAL_PARAMETERS

arguments
    timing (1,1) struct
end

t = timing;
if ~isfield(t,'pd') || ~isfield(t,'pri') || t.pri <= 0
    error('babelbrain_thermal_timing:missingPulse', ...
        'parameters.timing.pd and .pri must be set (pri > 0).');
end
if ~isfield(t, 'ptd')
    error('babelbrain_thermal_timing:missingPtd', ...
        'parameters.timing.ptd must be set (single pulse-train ON duration).');
end

bt = struct();
bt.duty_cycle = t.pd / t.pri;
bt.prf_hz     = 1 / t.pri;
bt.duration_s = t.ptd;   % ON duration of a single pulse-train repetition (NOT ptrd, the total across all reps)

if bt.duty_cycle <= 0 || bt.duty_cycle > 1
    error('babelbrain_thermal_timing:invalidDutyCycle', ...
        ['duty_cycle=%.4f is outside (0, 1].\n' ...
        'Check parameters.timing.pd=%.4f s and .pri=%.4f s.'], bt.duty_cycle, t.pd, t.pri);
end
if bt.prf_hz <= 0 || ~isfinite(bt.prf_hz)
    error('babelbrain_thermal_timing:invalidPrf', ...
        ['PRF=%.2f Hz is invalid (<= 0 or non-finite).\n' ...
        'Check parameters.timing.pri=%.4f s.'], bt.prf_hz, t.pri);
end
if bt.duration_s <= 0
    error('babelbrain_thermal_timing:invalidDuration', ...
        'duration=%.4f s is <= 0. Check parameters.timing.ptd.', bt.duration_s);
end

% Cooling interval between repetitions (BabelBrain DurationOff). ptri is the
% period of one repetition (ON+OFF); the off-time is the remainder after ptd.
% If ptri is not set or is not longer than ptd, there is no inter-train gap.
if isfield(t, 'ptri') && t.ptri > t.ptd
    bt.duration_off_s = t.ptri - t.ptd;
else
    bt.duration_off_s = 0.0;
end

% Number of pulse-train repetitions: n_ptri_reps = round(ptrd / ptri)
% (thermal_parameters.m).
bt.repetitions = 1;
if isfield(t, 'n_ptri_reps') && t.n_ptri_reps > 1
    bt.repetitions = round(t.n_ptri_reps);
end

% Post-sequence cooldown (BabelBrain PauseBetweenGroupedSonications with a
% single group). post_ptri_dur may be NaN/empty when unset.
bt.pause_after_s = 0.0;
if isfield(t, 'post_ptri_dur') && ~isempty(t.post_ptri_dur) && ...
        isfinite(t.post_ptri_dur) && t.post_ptri_dur > 0
    bt.pause_after_s = t.post_ptri_dur;
end

end
