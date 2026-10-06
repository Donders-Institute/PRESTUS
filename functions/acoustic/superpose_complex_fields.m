function out = superpose_complex_fields(files_in, options)
% SUPERPOSE_COMPLEX_FIELDS  Coherent sum of complex pressure fields from separate runs.
%
% Acoustic propagation is linear, so the field of several transducers firing
% together equals the sum of their complex fields simulated one at a time
% (same medium, grid and frequency). This combines p_complex files written
% with io.save_p_complex (one per transducer, e.g. the per-transducer runs of
% async mode) into the coherent field, and returns the incoherent and
% per-voxel-maximum combinations of the same fields for comparison:
%
%   coherent   |sum_i w_i p_i|        transducers fire together (interference)
%   incoherent sqrt(sum_i |w_i p_i|^2) time-averaged heat source of
%                                     non-overlapping firing (async thermal)
%   max        max_i |w_i p_i|        peak pressure of non-overlapping firing
%                                     (async acoustic safety)
%
% Partial-FOV fields are back-projected onto their full grid first
% (p_complex_info.fov_offset_ac / full_ac_dims); voxels outside a field's
% FOV get zero from that field, and out.coverage counts how many fields
% cover each voxel. All fields must share one full grid (async_grid_mode
% 'shared'): complex fields are not resampled between grids, because
% interpolating the phase of a field sampled at a few points per wavelength
% is not accurate.
%
% Phases must refer to the same time origin. Files written by ACOUSTIC_WRAPPER
% since p_complex_info was added refer to t = 0 of the simulation
% (phase_reference 'simulation_t0'). Older files (no p_complex_info) refer to
% the start of their own recording window; they are accepted only when all
% files are of that kind, with a warning, which is valid only when every run
% had the same time step and simulation length (same grid, no partial FOV).
%
% Use as:
%   out = superpose_complex_fields({'..._p_complex_trans-01.mat', '..._p_complex_trans-02.mat'})
%   out = superpose_complex_fields(files, Weights=[1, exp(1i*pi/2)])   % relative phase sweep
%   out = superpose_complex_fields(files, Scale=[s1, s2])              % per-transducer pressure scale
%
% Input:
%   files_in         - cellstr of p_complex .mat files (variables p_complex
%                      and, optionally, p_complex_info)
%   options.Weights  - complex weight per file (default 1): relative drive
%                      amplitude and phase applied before summing
%   options.Scale    - real pressure scale per file (default 1), e.g.
%                      sqrt(target / baseline) from ISPPA_SCALE_FACTOR
%
% Output:
%   out.p_coherent   - |sum| on the full grid [Pa]
%   out.p_incoherent - sqrt(sum of squares) [Pa]
%   out.p_max        - per-voxel maximum [Pa]
%   out.p_sum        - complex sum [Pa]
%   out.coverage     - uint8, number of fields that cover each voxel
%   out.full_dims    - full grid dimensions
%   out.freq_hz      - common frequency [Hz] (NaN when unknown)
%   out.files        - files_in
%
% See also: ACOUSTIC_WRAPPER, COMBINE_ASYNC_INTENSITY, ISPPA_SCALE_FACTOR

arguments
    files_in
    options.Weights (1,:) double = []
    options.Scale   (1,:) double = []
end

files_in = cellstr(files_in);
N = numel(files_in);
if N < 1
    error('superpose_complex_fields:noFiles', 'At least one p_complex file is required.');
end
w = ones(1, N);
if ~isempty(options.Weights)
    assert(numel(options.Weights) == N, 'superpose_complex_fields:badWeights', ...
        'Weights needs one entry per file (%d).', N);
    w = options.Weights;
end
if ~isempty(options.Scale)
    assert(numel(options.Scale) == N, 'superpose_complex_fields:badScale', ...
        'Scale needs one entry per file (%d).', N);
    w = w .* options.Scale;
end

infos = cell(1, N);
for k = 1:N
    if ~isfile(files_in{k})
        error('superpose_complex_fields:missingFile', 'p_complex file not found: %s', files_in{k});
    end
    vars = who('-file', files_in{k});
    if ~ismember('p_complex', vars)
        error('superpose_complex_fields:noField', 'No p_complex variable in %s', files_in{k});
    end
    if ismember('p_complex_info', vars)
        s = load(files_in{k}, 'p_complex_info');
        infos{k} = s.p_complex_info;
    end
end

% Phase reference: all t0-referenced, or all legacy (window-referenced).
has_info = ~cellfun(@isempty, infos);
if any(has_info) && ~all(has_info)
    error('superpose_complex_fields:mixedPhaseReference', ...
        ['Some files refer phases to t = 0 of the simulation and some to their recording ' ...
         'window: re-run the acoustic stage of the older runs with io.save_p_complex.']);
end
if ~any(has_info)
    warning('superpose_complex_fields:legacyPhaseReference', ...
        ['p_complex files carry no p_complex_info: phases refer to each run''s recording ' ...
         'window. The sum is valid only if all runs had the same time step and length.']);
end

% Common frequency.
freq_hz = NaN;
if all(has_info)
    f = cellfun(@(s) s.freq_hz, infos);
    if any(abs(f - f(1)) > 1e-6 * f(1))
        error('superpose_complex_fields:frequencyMismatch', ...
            ['Fields at different frequencies (%s Hz) do not interfere on average; ' ...
             'combine them incoherently instead.'], num2str(f));
    end
    freq_hz = f(1);
end

full_dims = [];
for k = 1:N
    s = load(files_in{k}, 'p_complex');
    p = double(s.p_complex);
    clear s
    [offset, fd] = field_placement(infos{k}, size(p));
    if isempty(full_dims)
        full_dims = fd;
        p_sum  = complex(zeros(full_dims));
        p_sq   = zeros(full_dims);
        p_max  = zeros(full_dims);
        cover  = zeros(full_dims, 'uint8');
    elseif ~isequal(fd, full_dims)
        error('superpose_complex_fields:gridMismatch', ...
            ['File %d is on a %s grid, file 1 on %s. Complex fields are only added on one ' ...
             'shared grid (async_grid_mode ''shared'').'], k, mat2str(fd), mat2str(full_dims));
    end
    idx = arrayfun(@(d) offset(d):offset(d) + size(p, d) - 1, 1:numel(full_dims), ...
        'UniformOutput', false);
    pk = w(k) * p;
    p_sum(idx{:}) = p_sum(idx{:}) + pk;
    p_sq(idx{:})  = p_sq(idx{:}) + abs(pk).^2;
    p_max(idx{:}) = max(p_max(idx{:}), abs(pk));
    cover(idx{:}) = cover(idx{:}) + 1;
end

out = struct('p_coherent', abs(p_sum), 'p_incoherent', sqrt(p_sq), 'p_max', p_max, ...
    'p_sum', p_sum, 'coverage', cover, 'full_dims', full_dims, 'freq_hz', freq_hz, ...
    'files', {files_in});
end

function [offset, full_dims] = field_placement(info, sz)
% Where a (possibly FOV-cropped) field sits in its full grid.
offset    = ones(1, numel(sz));
full_dims = sz;
if isempty(info) || ~isfield(info, 'fov_offset_ac') || isempty(info.fov_offset_ac)
    return
end
offset    = double(info.fov_offset_ac(:)');
full_dims = double(info.full_ac_dims(:)');
if any(offset + sz - 1 > full_dims)
    error('superpose_complex_fields:fovExceedsGrid', ...
        'FOV field of size %s at offset %s exceeds the full grid %s.', ...
        mat2str(sz), mat2str(offset), mat2str(full_dims));
end
end
