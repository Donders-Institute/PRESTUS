function ov = thermal_overlap_terms(sched, source_freq, source_tx, phase_model)
% THERMAL_OVERLAP_TERMS  Interference of heat sources whose pulses partly overlap
%
% Two heat sources at one frequency that are on at the same time interfere:
% while both are on, the heat source is |P_a + P_b e^(i theta)|^2 instead of
% |P_a|^2 + |P_b|^2, which adds the interference term
%
%   X(theta) = 2 Re(P_a conj(P_b) e^(-i theta)) = cos(theta) A + sin(theta) B,
%   A = 2 Re(P_a conj(P_b)),  B = 2 Im(P_a conj(P_b))
%
% A and B are fields (COMBINE_ASYNC_INTENSITY); this function says in which
% segments of the schedule each pair is on together, and with which phase
% theta. Sources of one firing group (same SourceTx) fire together by design
% and are summed when combined, so they never form a pair here; nor do sources
% at different frequencies, whose cross terms average out.
%
% theta, the phase of source b's carrier relative to the simulated fields:
%   'clock' - the channels share one clock: the carriers keep the phases they
%             were simulated with during every pulse; theta = 0.
%   'burst' - each pulse starts its carrier at its own onset (burst mode
%             triggered per pulse): theta = -2*pi*f*(t_on_b - t_on_a), with
%             t_on the onsets of the two pulses that overlap.
%
% Use as:
%   ov = thermal_overlap_terms(sched, source_freq, source_tx, phase_model)
%
% Input:
%   sched       - struct from THERMAL_SCHEDULE (t0, dur, active, source(i).pulse_windows,
%                 pulse_model)
%   source_freq - [1 x n_sources] frequency of each source [Hz]
%   source_tx   - [1 x n_sources] transducer whose windows each source follows
%   phase_model - 'clock' or 'burst'
%
% Output:
%   ov - struct with fields
%     pairs   [K x 2] source indices of each interfering pair (a < b)
%     active  [S x K] logical, both sources of the pair on in the segment
%     theta   [S x K] phase of source b relative to source a [rad] (NaN when not both on)
%     cos_w, sin_w [S x K] weights of the fields A and B (0 when not both on)
%     thetas  {1 x K} distinct phases of each pair (for peak pressures)
%
% See also: THERMAL_SCHEDULE, COMBINE_ASYNC_INTENSITY, THERMAL_SIMULATION

arguments
    sched        (1,1) struct
    source_freq  (1,:) double
    source_tx    (1,:) double
    phase_model  (1,:) char {mustBeMember(phase_model, {'clock', 'burst'})}
end

n = size(sched.active, 2);
S = size(sched.active, 1);
pairs = zeros(0, 2);
for a = 1:n
    for b = a+1:n
        same_f = isfinite(source_freq(a)) && isfinite(source_freq(b)) && ...
                 abs(source_freq(a) - source_freq(b)) <= 1e-6 * max(source_freq(a), source_freq(b));
        if same_f && source_tx(a) ~= source_tx(b) && any(sched.active(:, a) & sched.active(:, b))
            pairs(end+1, :) = [a, b]; %#ok<AGROW>
        end
    end
end
K = size(pairs, 1);
ov = struct('pairs', pairs, 'active', false(S, K), 'theta', NaN(S, K), ...
            'cos_w', zeros(S, K), 'sin_w', zeros(S, K), 'thetas', {cell(1, K)});
if K == 0
    return
end
if isfield(sched, 'pulse_model') && strcmp(sched.pulse_model, 'train_averaged')
    error('thermal_overlap_terms:trainAveraged', ...
        ['Sources %d and %d partly overlap at one frequency, which needs the pulses resolved: ' ...
         'set thermal.pulse_model to ''pulse_resolved''.'], pairs(1, 1), pairs(1, 2));
end

mid = sched.t0 + sched.dur / 2;
for kk = 1:K
    a = pairs(kk, 1); b = pairs(kk, 2);
    on = sched.active(:, a) & sched.active(:, b);
    ov.active(:, kk) = on;
    switch phase_model
        case 'clock'
            th = zeros(nnz(on), 1);
        case 'burst'
            ta = pulse_onset(sched.source(a).pulse_windows, mid(on));
            tb = pulse_onset(sched.source(b).pulse_windows, mid(on));
            th = wrap(-2 * pi * source_freq(a) * (tb - ta));
    end
    ov.theta(on, kk) = th;
    ov.cos_w(on, kk) = cos(th);
    ov.sin_w(on, kk) = sin(th);
    ov.thetas{kk} = uniquetol(th(:)', 1e-9);
end
end

% -------------------------------------------------------------------------
function t_on = pulse_onset(windows, t)
% Start of the pulse (row of windows = [start end]) that contains each t
    idx = discretize(t, [windows(:, 1); Inf]);
    t_on = windows(idx, 1);
end

function th = wrap(th)
% Phase in (-pi, pi]; exact multiples of 2*pi become 0
    th = mod(th + pi, 2 * pi) - pi;
    th(abs(th + pi) < 1e-12) = pi;
end
