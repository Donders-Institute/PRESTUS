function [p_complex, info] = reverse_point_source_field(parameters, kgrid, kwave_medium, sensor, target_vox, sensor_mask, freq_hz)
% REVERSE_POINT_SOURCE_FIELD  Complex CW field of a point source at the target, at the sensor points
%
% One k-Wave run with a continuous-wave point source at the target voxel, on
% the grid, medium and time axis of the forward simulation. By acoustic
% reciprocity, the complex pressure this source produces at a point x equals
% (up to one complex constant) the pressure at the target produced by a unit
% source at x. Sampled over a transducer aperture, it therefore gives each
% element's transfer to the target through the skull (TIME_REVERSAL_REFOCUS),
% and sampled over a region it scores any transducer pose placed there
% (RECIPROCITY_SCORE_POSES), without one forward run per element or pose.
%
% The constant is the same for all points of one run, so phases and ratios
% between points or poses are exact; absolute values are not. Reciprocity
% assumes a linear, fluid medium, and that the source points (transducer
% surface) all lie in the same medium (coupling water): k-Wave scales an
% additive pressure source by the local sound speed.
%
% Phase convention as ACOUSTIC_WRAPPER (p_complex): extractAmpPhase phase
% (signal ~ cos(wt + phase)), referred to t = 0 of the simulation. An outgoing
% wave from the target is then ~ exp(-1i*k*r).
%
% Use as:
%   [p_complex, info] = reverse_point_source_field(parameters, kgrid, kwave_medium, ...
%                           sensor, target_vox, sensor_mask, freq_hz)
%
% Input:
%   parameters   - PRESTUS config (simulation.code_type, grid.pml_size, ...)
%   kgrid        - kWaveGrid of the forward run, time axis set
%   kwave_medium - k-Wave medium struct of the forward run
%   sensor       - forward sensor (record_start_index sets the steady-state window)
%   target_vox   - [1 x n_dims] target voxel (rounded to the nearest grid point)
%   sensor_mask  - logical grid, points at which the field is returned
%   freq_hz      - source frequency [Hz]
%
% Output:
%   p_complex - [n_points x 1] complex pressure at find(sensor_mask), in that order
%   info      - struct: target_vox, n_points, record_start_index, code_type
%
% See also: TIME_REVERSAL_REFOCUS, RECIPROCITY_SCORE_POSES, TIME_REVERSAL_STAGE,
%           ACOUSTIC_SIMULATION

arguments
    parameters   (1,1) struct
    kgrid        (1,1)
    kwave_medium (1,1) struct
    sensor       (1,1) struct
    target_vox   (1,:) {mustBeNumeric}
    sensor_mask  {mustBeNumericOrLogical}
    freq_hz      (1,1) {mustBeNumeric, mustBePositive}
end

    dims = parameters.grid.dims;
    if numel(dims) == 2 && isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1
        error('PRESTUS:timeReversal:axisym', ...
            'Time reversal needs a 2D or 3D Cartesian grid; axisymmetric grids are not supported.');
    end
    target_vox = round(target_vox);
    if any(target_vox < 1) || any(target_vox > dims)
        error('PRESTUS:timeReversal:targetOutsideGrid', ...
            'Target voxel [%s] lies outside the simulation grid [%s].', num2str(target_vox), num2str(dims));
    end

    % Point source at the target, unit amplitude, zero phase
    source = struct();
    source.p_mask = false(dims);
    idx = num2cell(target_vox);
    source.p_mask(idx{:}) = true;
    source.p = cast(createCWSignals(kgrid.t_array, freq_hz, 1, 0), parameters.simulation.precision);

    rev_sensor = struct();
    rev_sensor.mask   = logical(sensor_mask);
    rev_sensor.record = {'p'};
    rev_sensor.record_start_index = sensor.record_start_index;

    % PML as ACOUSTIC_WRAPPER (outside the grid, 'auto' resolved by k-Wave's rule)
    pml = parameters.grid.pml_size;
    if (ischar(pml) || isstring(pml)) && strcmp(pml, 'auto')
        pml = getOptimalPMLSize(kgrid);
    end
    input_args = struct('PMLInside', false, 'PMLSize', double(pml), 'PlotPML', false);

    % The reverse run is a fluid k-Wave run whatever the forward backend: the
    % jwave/BabelBrain forward runs use the same medium, an elastic forward run
    % is approximated by its fluid medium.
    rev_parameters = parameters;
    if isfield(rev_parameters, 'acoustic_backend') && ~strcmpi(rev_parameters.acoustic_backend, 'kwave')
        if strcmpi(rev_parameters.acoustic_backend, 'kwave_elastic')
            warning('PRESTUS:timeReversal:fluidReverse', ...
                'Time reversal uses a fluid k-Wave run; shear waves of the elastic forward model are not included.');
        end
        rev_parameters.acoustic_backend = 'kwave';
    end
    % Name the C++ scratch files apart from the forward run's
    if isfield(rev_parameters.io, 'output_affix')
        rev_parameters.io.output_affix = [rev_parameters.io.output_affix '_timereversal'];
    end

    fprintf('[time_reversal] reverse run: point source at voxel [%s], %d sensor points, %s\n', ...
        num2str(target_vox), nnz(rev_sensor.mask), rev_parameters.simulation.code_type);
    sensor_data = acoustic_simulation(kgrid, kwave_medium, source, rev_sensor, input_args, rev_parameters);
    p_t = sensor_data.p;                                  % [n_points x n_t], kept in its class (single)
    clear sensor_data

    % Complex amplitude at freq_hz by a single-frequency DFT over the last
    % whole number of periods of the recorded window (exact for a steady
    % state), in blocks of points, so the recording is never copied to double
    % as a whole. Same convention as extractAmpPhase: signal = a*cos(w*t + phi)
    % -> a*exp(1i*phi), referred to t = 0 of the simulation.
    n_t   = size(p_t, 2);
    ppp   = 1 / (freq_hz * kgrid.dt);                     % points per period
    n_use = max(1, round(floor(n_t / ppp) * ppp));
    if n_use > n_t, n_use = n_t; end
    k0    = n_t - n_use;                                  % first sample used (0-based in the window)
    w     = 2 * pi * freq_hz * kgrid.dt;
    E     = exp(-1i * w * (0:n_use-1)).' * (2 / n_use);   % [n_use x 1]
    t0_s  = (double(sensor.record_start_index) - 1 + k0) * kgrid.dt;
    p_complex = complex(zeros(size(p_t, 1), 1));
    blk = max(1, floor(2e8 / (8 * n_use)));               % ~200 MB of doubles per block
    for i0 = 1:blk:size(p_t, 1)
        r = i0:min(i0 + blk - 1, size(p_t, 1));
        p_complex(r) = double(p_t(r, k0+1:end)) * E;
    end
    clear p_t
    p_complex = p_complex .* exp(-1i * 2 * pi * freq_hz * t0_s);

    info = struct('target_vox', target_vox, 'n_points', numel(p_complex), ...
        'record_start_index', sensor.record_start_index, ...
        'code_type', rev_parameters.simulation.code_type);
end
