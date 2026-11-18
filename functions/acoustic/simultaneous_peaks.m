function [pk_pos, pk_neg, info] = simultaneous_peaks(P, f, pd_s)
% SIMULTANEOUS_PEAKS  Peak pressures of tones at different frequencies firing together
%
% Peak positive and negative pressure of sum_i Re(P_i * exp(1i*2*pi*f_i*t))
% over t in [0, T]: T = the common period 1/gcd(f_i), or the pulse duration
% when shorter (each pulse starts at the phases of t = 0). Sampled at 64 points
% per period of the highest frequency (peak underestimated by at most 0.12 %).
% When that needs more than 20000 samples, returns the upper bound
% +/- sum_i |P_i| instead (info.method 'upper_bound').
%
% Use as:
%   [pk_pos, pk_neg, info] = simultaneous_peaks(P, f, pd_s)
%
% Input:
%   P    - cell array of complex amplitude arrays [Pa], all the same size
%   f    - frequencies [Hz], one per element of P
%   pd_s - pulse duration [s]; NaN or Inf = the full common period
%
% Output:
%   pk_pos, pk_neg - peak positive and negative pressure [Pa], size of P{1}
%   info           - .common_period_s, .evaluated_s, .n_samples, .frequencies_hz,
%                    .method ('phase_sampled' or 'upper_bound')
%
% See also: COMBINE_ASYNC_INTENSITY, STEADY_STATE_COMPLEX_AMPLITUDE
    g = round(f(1));
    for k = 2:numel(f)
        g = gcd(g, round(f(k)));
    end
    T_common = 1 / g;
    T = T_common;
    if isfinite(pd_s) && pd_s > 0 && pd_s < T_common
        T = pd_s;
    end
    n_t = ceil(T * max(f) * 64) + 1;
    info = struct('common_period_s', T_common, 'evaluated_s', T, 'n_samples', n_t, ...
                  'frequencies_hz', f, 'method', 'phase_sampled');
    if n_t > 20000
        info.method = 'upper_bound';
        pk_pos = zeros(size(P{1}));
        for k = 1:numel(P)
            pk_pos = pk_pos + abs(P{k});
        end
        pk_neg = -pk_pos;
        fprintf(['[simultaneous_peaks] Common period %.3g s needs %d samples: peak pressures are the ' ...
                 'upper bound sum of amplitudes.\n'], T_common, n_t);
        return
    end
    pk_pos = -inf(size(P{1}));
    pk_neg =  inf(size(P{1}));
    for it = 1:n_t
        t = (it - 1) * T / (n_t - 1);
        v = zeros(size(P{1}));
        for k = 1:numel(P)
            v = v + real(P{k} * exp(1i * 2 * pi * f(k) * t));
        end
        pk_pos = max(pk_pos, v);
        pk_neg = min(pk_neg, v);
    end
    fprintf(['[simultaneous_peaks] Peak pressure of the tones firing together: %d samples over ' ...
             '%.3g s (common period %.3g s).\n'], n_t, T, T_common);
end
