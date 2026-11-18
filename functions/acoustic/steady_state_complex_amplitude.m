function p_complex = steady_state_complex_amplitude(p_t, dt, freq_hz, record_start_index)
% STEADY_STATE_COMPLEX_AMPLITUDE  Complex pressure amplitude at given frequencies from a recorded window
%
% Least-squares fit of a*cos(w*t) + b*sin(w*t) to each point's recorded
% samples, one fit per frequency. Exact for a steady-state tone, also when
% a period is not a whole number of time steps or the window is shorter
% than a period. Content that is not the steady tone (residual transients)
% enters the estimate more strongly the shorter the window. Convention as
% extractAmpPhase: a signal A*cos(w*t + phi) gives A*exp(1i*phi), with phi
% referred to t = 0 of the simulation, so fields from separate runs (e.g.
% one per transducer, with different dt) can be added coherently.
% Processed in blocks of points, so the recording is never copied as a
% whole (memory: the recording + ~200 MB).
%
% Use as:
%   p_complex = steady_state_complex_amplitude(p_t, dt, freq_hz, record_start_index)
%
% Input:
%   p_t                - recorded pressure, [n_points x n_t] (k-Wave sensor_data.p)
%   dt                 - time step [s] (kgrid.dt)
%   freq_hz            - frequency or frequencies to extract [Hz]
%   record_start_index - time index of the first recorded sample (sensor.record_start_index;
%                        1 or [] when the whole simulation was recorded)
%
% Output:
%   p_complex - [n_points x numel(freq_hz)] complex double

    if nargin < 4 || isempty(record_start_index); record_start_index = 1; end
    n_pts = size(p_t, 1);
    n_t   = size(p_t, 2);

    nf = numel(freq_hz);
    k  = (0:n_t-1).';
    M  = zeros(n_t, 2 * nf);                                % least-squares projectors, [cos sin] per frequency
    for fi = 1:nf
        w = 2 * pi * freq_hz(fi) * dt;
        B = [cos(w * k), sin(w * k)];
        M(:, 2*fi-1:2*fi) = B / (B.' * B);
    end

    C = zeros(n_pts, 2 * nf);
    blk = max(1, floor(2e8 / (8 * n_t)));                   % ~200 MB of doubles per block
    for i0 = 1:blk:n_pts
        r = i0:min(i0 + blk - 1, n_pts);
        C(r, :) = double(p_t(r, :)) * M;
    end

    t0_s = (double(record_start_index) - 1) * dt;
    p_complex = complex(zeros(n_pts, nf));
    for fi = 1:nf
        % A*cos(w*t+phi) = A*cos(phi)*cos(w*t) - A*sin(phi)*sin(w*t)
        p_complex(:, fi) = (C(:, 2*fi-1) - 1i * C(:, 2*fi)) .* exp(-1i * 2 * pi * freq_hz(fi) * t0_s);
    end
end
