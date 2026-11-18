function fg = firing_groups(parameters)
% FIRING_GROUPS  Which transducers fire together, and the heat sources they form
%
% Transducers that fire together form a group: their fields are added as
% complex amplitudes (they interfere) by COMBINE_ASYNC_INTENSITY. Groups fire
% in their own windows, like the transducers of an async run.
%
%   - Async coupling with timing.async_schedule 'custom': transducers with the
%     same timing_offset_s and the same effective pulse timing (the timing
%     block with transducer(i).timing applied: pd, pri, ptd, ptri, ptrd) form
%     one group. Transducers whose windows overlap without being identical
%     stay separate (THERMAL_SCHEDULE rejects that overlap at one frequency).
%   - Coherent transducers at different frequencies (IS_SIMULTANEOUS_MULTIFREQ):
%     one group of all transducers.
%   - Otherwise every transducer is its own group.
%
% Each group gives one heat source per frequency: members at one frequency are
% summed into one field, and fields at different frequencies heat with their
% own absorption (THERMAL_SIMULATION). Sources are ordered by group, then by
% frequency; a group is numbered by its first transducer.
%
% Use as:
%   fg = firing_groups(parameters)
%
% Output:
%   fg - struct with fields
%     group         [1 x N] group of each transducer
%     n_groups      number of groups
%     members       {1 x G} transducer indices of each group
%     coherent      [1 x G] logical, group has more than one transducer
%     source_tx     [1 x S] first transducer of each source's group: the
%                   source fires in that transducer's windows (THERMAL_SCHEDULE)
%     source_freq   [1 x S] frequency of each source [Hz]
%     source_group  [1 x S] group of each source
%
% See also: COMBINE_ASYNC_INTENSITY, THERMAL_SCHEDULE, ASYNC_TRANSDUCER_PIPELINE

arguments
    parameters (1,1) struct
end

N = 0;
if isfield(parameters, 'transducer')
    N = numel(parameters.transducer);
end
freqs = NaN(1, N);
for i = 1:N
    if isfield(parameters.transducer(i), 'freq_hz') && ~isempty(parameters.transducer(i).freq_hz)
        freqs(i) = double(parameters.transducer(i).freq_hz);
    end
end

if N <= 1
    group = ones(1, N);
elseif is_simultaneous_multifreq(parameters)
    group = ones(1, N);
elseif is_custom_async(parameters)
    key = zeros(N, 6);
    for i = 1:N
        key(i, :) = window_key(parameters, i);
    end
    group = zeros(1, N);
    g = 0;
    for i = 1:N
        if group(i) > 0, continue, end
        g = g + 1;
        group(i) = g;
        for j = i+1:N
            if group(j) == 0 && same_key(key(i, :), key(j, :))
                group(j) = g;
            end
        end
    end
else
    group = 1:N;
end

G = max([group, 0]);
members  = arrayfun(@(g) find(group == g), 1:G, 'UniformOutput', false);
coherent = cellfun(@numel, members) > 1;

source_tx = []; source_freq = []; source_group = [];
for g = 1:G
    m = members{g};
    f = unique_freqs(freqs(m));
    source_tx    = [source_tx,    repmat(m(1), 1, numel(f))]; %#ok<AGROW>
    source_freq  = [source_freq,  f];                          %#ok<AGROW>
    source_group = [source_group, repmat(g, 1, numel(f))];     %#ok<AGROW>
end

fg = struct('group', group, 'n_groups', G, 'members', {members}, 'coherent', coherent, ...
            'source_tx', source_tx, 'source_freq', source_freq, 'source_group', source_group);
end

% -------------------------------------------------------------------------
function tf = is_custom_async(parameters)
tf = isfield(parameters, 'simulation') && isfield(parameters.simulation, 'transducer_coupling') && ...
     strcmp(char(parameters.simulation.transducer_coupling), 'async') && ...
     isfield(parameters, 'timing') && isfield(parameters.timing, 'async_schedule') && ...
     strcmp(char(parameters.timing.async_schedule), 'custom');
end

function k = window_key(parameters, i)
% [offset, pd, pri, ptd, ptri, ptrd] of transducer i, NaN where unset
fields = {'pd', 'pri', 'ptd', 'ptri', 'ptrd'};
k = NaN(1, 6);
if isfield(parameters.transducer(i), 'timing_offset_s') && ~isempty(parameters.transducer(i).timing_offset_s)
    k(1) = double(parameters.transducer(i).timing_offset_s);
else
    k(1) = 0;
end
ov = struct();
if isfield(parameters.transducer(i), 'timing') && isstruct(parameters.transducer(i).timing)
    ov = parameters.transducer(i).timing;
end
for f = 1:numel(fields)
    name = fields{f};
    if isfield(ov, name) && ~isempty(ov.(name))
        k(f + 1) = double(ov.(name));
    elseif isfield(parameters.timing, name) && ~isempty(parameters.timing.(name))
        k(f + 1) = double(parameters.timing.(name));
    end
end
end

function tf = same_key(a, b)
% Equal within 1 us on the offset and a relative 1e-9 on the timing (NaN = NaN)
both_nan = isnan(a) & isnan(b);
d = abs(a - b);
tol = [1e-6, 1e-9 * max(1, abs(a(2:end)))];
tf = all(both_nan | d <= tol);
end

function f = unique_freqs(f)
% Sorted distinct frequencies; NaN (unknown) counts as one frequency
if all(isnan(f))
    f = NaN;
    return
end
f = f(~isnan(f));
f = reshape(uniquetol(sort(f), 1e-6, 'DataScale', max(f)), 1, []);
end
