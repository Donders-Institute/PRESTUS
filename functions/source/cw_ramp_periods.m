function n = cw_ramp_periods(t_array, freq_hz)
% CW_RAMP_PERIODS  Start-ramp length for CREATECWSIGNALS that fits the simulated time
%
% k-Wave's createCWSignals tapers the first 4 periods of each source signal.
% A simulation shorter than that (e.g. a shortened run for timing) would index
% past the end of the signal. This returns 4 periods, or half the simulated
% periods when the run is shorter than 8 periods (0 when that ramp would span
% fewer than 2 time steps). Runs of 8 periods or more keep k-Wave's default,
% so their sources are unchanged.
%
% Use as:
%   cw = createCWSignals(t_array, f, amp, phase, cw_ramp_periods(t_array, f))
%
% Input:
%   t_array - kgrid.t_array [s]
%   freq_hz - source frequency [Hz]
%
% Output:
%   n - ramp length in periods (createCWSignals ramp_length)
%
% See also: SOURCE_CREATE, CREATECWSIGNALS

    n = 4;
    t_array = t_array(:);
    if numel(t_array) < 2
        n = 0;
        return
    end
    dt = t_array(2) - t_array(1);
    n_periods = (numel(t_array) * dt) * freq_hz;
    if n_periods >= 2 * n
        return
    end
    n = n_periods / 2;
    if round(n / freq_hz / dt) < 2
        n = 0;
    end
    fprintf('[cw_ramp_periods] %.2f simulated periods: start ramp shortened from 4 to %.2f periods.\n', ...
        n_periods, n);
end
