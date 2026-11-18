function sched = thermal_schedule(parameters, n_sources, opts)
% THERMAL_SCHEDULE  Firing schedule of one or more heat sources as thermal time segments
%
% Builds the on-windows of every source (transducer) from the timing block and
% splits the protocol into segments between consecutive switching times. In
% each segment a fixed set of sources is on, so the heat source is constant:
% Q = sum of Q_i over the sources on. THERMAL_SIMULATION steps through the
% segments in order.
%
% One source reproduces the hierarchy of THERMAL_PARAMETERS exactly: one
% segment per pulse ON and per pulse OFF inside a pulse train (pt_timestep
% steps when timing.equal_step_duration is 1, else one step), and one segment
% per post_pt_timestep step between trains and after the PTRD.
%
% Several sources (async coupling) follow timing.async_schedule:
%   'pulses' (default) - every transducer uses the timing block; transducer i
%                        starts (i-1)*PRI/N later, so pulses alternate within
%                        each pulse repetition interval. Needs N*PD <= PRI.
%   'trains'           - every transducer uses the timing block; transducer i
%                        starts (i-1)*PTRI/N later, so pulse trains alternate.
%                        Needs N*PTD <= PTRI.
%   'simultaneous'     - every transducer uses the timing block and fires at the
%                        same time; set internally for coherent transducers at
%                        different frequencies (IS_SIMULTANEOUS_MULTIFREQ), whose
%                        time-averaged heat sources add.
%   'custom'           - transducer(i).timing_offset_s [s] (default 0) delays
%                        transducer i; transducer(i).timing may override pd,
%                        pri, ptd, ptri and ptrd of the timing block for that
%                        transducer, so transducers fired one after another
%                        (a sequence of stances) can each sonicate for their
%                        own duration. post_ptri_dur and the time steps stay
%                        global; the protocol ends post_ptri_dur after the
%                        last pulse train of any transducer (or timing.ptrd).
% Time between pulse trains of every source is stepped at post_pt_timestep;
% time inside any source's pulse train at pt_timestep (or one step per
% segment). Overlapping on-windows of transducers at the same frequency are an
% error: their fields interfere, which the incoherent sum does not model (fire
% them together as a group, FIRING_GROUPS, or use transducer_coupling
% 'coherent'). At different frequencies they warn.
%
% Firing groups: by default source i fires in transducer i's windows. With
% SourceTx, source i fires in transducer SourceTx(i)'s windows, so the fields
% of one group (one per frequency, FIRING_GROUPS) share that group's windows;
% their overlap is intended and neither warns nor errors.
%
% Use as:
%   sched = thermal_schedule(parameters)              % one source
%   sched = thermal_schedule(parameters, n_sources)
%   sched = thermal_schedule(parameters, n_sources, SourceTx=tx, SourceFreq=f)
%
% Input:
%   parameters - PRESTUS config; uses timing.*, timing.async_schedule,
%                transducer(i).timing_offset_s, transducer(i).timing,
%                transducer(i).freq_hz
%   n_sources  - number of heat sources (default 1)
%   SourceTx   - (optional) [1 x n_sources] transducer whose windows each
%                source follows (default 1:n_sources)
%   SourceFreq - (optional) [1 x n_sources] frequency of each source [Hz]
%                (default: the frequency of its transducer)
%
% Output:
%   sched - struct with fields
%     mode         'single' | 'pulses' | 'trains' | 'custom'
%     n_sources    number of sources
%     t0, dur      [S x 1] segment start and duration [s]
%     n_steps      [S x 1] kWaveDiffusion steps per segment
%     step_dur     [S x 1] duration of one step [s]
%     active       [S x n_sources] logical, source on during the segment
%     weight       [S x n_sources] heat-source weight: 1 when on ('pulse_resolved'),
%                  PD/PRI through each pulse train ('train_averaged')
%     pulse_model  'train_averaged' (default) | 'pulse_resolved'; interleaved async
%                  schedules always use 'pulse_resolved'
%     fine         [S x 1] logical, segment inside a pulse train
%     t_sonication end of the last pulse train (>= ptrd); with 'custom', also
%                  >= offset_s + ptrd of every source [s]
%     t_end        t_sonication + post_ptri_dur [s]
%     total_steps  sum(n_steps)
%     source(i)    pd, pri, ptd, ptri, ptrd, offset_s, n_pulses_per_pt,
%                  n_ptri_reps, duty_cycle (overall), on_time_s,
%                  pulse_windows [P x 2], train_windows [R x 2] [s]
%
% See also: THERMAL_PARAMETERS, THERMAL_SIMULATION, THERMAL_PLOT_PROTOCOL

arguments
    parameters (1,1) struct
    n_sources  (1,1) double {mustBeInteger, mustBePositive} = 1
    opts.SourceTx   (1,:) double = []
    opts.SourceFreq (1,:) double = []
end

timing = parameters.timing;
n = n_sources;
grouped = ~isempty(opts.SourceTx);
source_tx = 1:n;
if grouped
    source_tx = opts.SourceTx;
    if numel(source_tx) ~= n || any(source_tx < 1) || any(source_tx ~= round(source_tx))
        error('thermal_schedule:badSourceTx', ...
            'SourceTx must hold one transducer index per heat source (%d).', n);
    end
end
% Transducer slots of 'pulses' / 'trains': one per transducer whose windows a source follows
[~, ~, slot] = unique(source_tx, 'stable');
slot = slot(:)';
n_slots = max(slot);

n_tx = 0;
if isfield(parameters, 'transducer')
    n_tx = numel(parameters.transducer);
end
if n == 1 && ~(grouped && n_tx > 1)
    mode = 'single';
else
    mode = 'pulses';
    if isfield(timing, 'async_schedule') && ~isempty(timing.async_schedule)
        mode = char(timing.async_schedule);
    end
    if ~ismember(mode, {'pulses', 'trains', 'custom', 'simultaneous'})
        error('thermal_schedule:unknownMode', ...
            'timing.async_schedule must be ''pulses'', ''trains'' or ''custom'', got ''%s''.', mode);
    end
    if ~isfield(parameters, 'transducer') || n_tx < max(source_tx)
        error('thermal_schedule:tooFewTransducers', ...
            '%d heat sources but only %d transducers in parameters.transducer.', ...
            max(source_tx), n_tx);
    end
end
source_freq = opts.SourceFreq;
if isempty(source_freq)
    source_freq = NaN(1, n);
    if ~strcmp(mode, 'single')
        source_freq = arrayfun(@(j) tx_freq(parameters, j), source_tx);
    end
elseif numel(source_freq) ~= n
    error('thermal_schedule:badSourceFreq', 'SourceFreq must hold one frequency per heat source (%d).', n);
end

%% Per-source protocol and offset
OVERRIDABLE = {'pd', 'pri', 'ptd', 'ptri', 'ptrd'};
src = struct('pd', {}, 'pri', {}, 'ptd', {}, 'ptri', {}, 'ptrd', {}, 'offset_s', {}, ...
             'n_pulses_per_pt', {}, 'n_ptri_reps', {}, 'duty_cycle', {}, 'on_time_s', {}, ...
             'pulse_windows', {}, 'train_windows', {});
for i = 1:n
    j = source_tx(i);          % the transducer whose windows source i follows
    t = timing;
    ov = struct();
    if ~strcmp(mode, 'single') && isfield(parameters.transducer(j), 'timing') && ...
            isstruct(parameters.transducer(j).timing)
        ov = parameters.transducer(j).timing;
    end
    ov_set = fieldnames(ov);
    ov_set = ov_set(cellfun(@(f) ~isempty(ov.(f)), ov_set));
    if ~isempty(ov_set)
        if ~strcmp(mode, 'custom')
            error('thermal_schedule:overrideNeedsCustom', ...
                ['transducer(%d).timing sets %s, which only async_schedule ''custom'' reads. ' ...
                 'Set timing.async_schedule: custom, or remove the override.'], ...
                j, strjoin(ov_set, ', '));
        end
        bad = setdiff(ov_set, OVERRIDABLE);
        if ~isempty(bad)
            error('thermal_schedule:badOverride', ...
                'transducer(%d).timing may set only %s; got %s.', ...
                j, strjoin(OVERRIDABLE, ', '), strjoin(bad, ', '));
        end
        for f = ov_set'
            t.(f{1}) = ov.(f{1});
        end
    end
    p_i = parameters;
    p_i.timing = t;
    if ~isfield(p_i, 'thermal'), p_i.thermal = struct(); end
    pt = thermal_parameters(p_i, true);     % validates integer pulse / train counts

    switch mode
        case {'single', 'simultaneous'}
            offset = 0;
        case 'pulses'
            offset = (slot(i) - 1) * pt.pri / n_slots;
            if n_slots * pt.pd > pt.pri * (1 + 1e-9)
                error('thermal_schedule:pulsesOverlap', ...
                    ['async_schedule ''pulses'' fires the %d transducers in turn within each pulse ' ...
                     'repetition interval, which needs N x PD <= PRI (%d x %.4g ms > %.4g ms). ' ...
                     'Use ''trains'', a shorter PD, or ''custom''.'], ...
                    n_slots, n_slots, pt.pd * 1e3, pt.pri * 1e3);
            end
        case 'trains'
            offset = (slot(i) - 1) * pt.ptri / n_slots;
            if ~(pt.ptri > 0) || n_slots * pt.ptd > pt.ptri * (1 + 1e-9)
                error('thermal_schedule:trainsOverlap', ...
                    ['async_schedule ''trains'' fires the %d transducers in turn within each pulse ' ...
                     'train repetition interval, which needs N x PTD <= PTRI (%d x %.4g s > %.4g s). ' ...
                     'Use ''pulses'', a shorter PTD, or ''custom''.'], ...
                    n_slots, n_slots, pt.ptd, pt.ptri);
            end
        case 'custom'
            offset = 0;
            if isfield(parameters.transducer(j), 'timing_offset_s') && ...
                    ~isempty(parameters.transducer(j).timing_offset_s)
                offset = double(parameters.transducer(j).timing_offset_s);
            end
            if ~isscalar(offset) || ~isfinite(offset) || offset < 0
                error('thermal_schedule:badOffset', ...
                    'transducer(%d).timing_offset_s must be a finite scalar >= 0 [s].', j);
            end
    end

    np = pt.n_pulses_per_pt;
    nr = pt.n_ptri_reps;
    k  = (0:np-1)' * pt.pri;
    r  = (0:nr-1)  * pt.ptri;
    starts = offset + reshape(k + r, [], 1);            % pulses of train 1, then train 2, ...
    tr_starts = offset + r(:);

    s = struct();
    s.pd = pt.pd; s.pri = pt.pri; s.ptd = pt.ptd; s.ptri = pt.ptri; s.ptrd = pt.ptrd;
    s.offset_s        = offset;
    s.n_pulses_per_pt = np;
    s.n_ptri_reps     = nr;
    s.on_time_s       = pt.pd * np * nr;
    s.duty_cycle      = s.on_time_s / pt.ptrd;
    s.pulse_windows   = [starts, starts + pt.pd];
    s.train_windows   = [tr_starts, tr_starts + pt.ptd];
    src(i) = s; %#ok<AGROW>
end

%% Pulse model: each pulse resolved, or each pulse train as one window
% 'train_averaged' heats through the whole train with the source weighted by
% PD/PRI (same energy), as k-Plan's thermal solver does; the diffusion steps
% are then set by THERMAL_SIMULATION. Pulses (0.1-1 s cycles) are short
% against tissue thermal time constants, so peak and end temperatures agree.
pulse_model = 'train_averaged';
if isfield(parameters, 'thermal') && isfield(parameters.thermal, 'pulse_model') && ...
        ~isempty(parameters.thermal.pulse_model)
    pulse_model = lower(char(parameters.thermal.pulse_model));
end
if ~ismember(pulse_model, {'pulse_resolved', 'train_averaged'})
    error('thermal_schedule:pulseModel', ...
        'thermal.pulse_model must be ''train_averaged'' or ''pulse_resolved'', got ''%s''.', pulse_model);
end
if strcmp(pulse_model, 'train_averaged') && ~ismember(mode, {'single', 'simultaneous'})
    % interleaved transducers fire in turn within a train: resolve the pulses
    fprintf('[thermal_schedule] async schedule ''%s'': pulses resolved (thermal.pulse_model ''train_averaged'' needs sources firing together)\n', mode);
    pulse_model = 'pulse_resolved';
end
averaged = strcmp(pulse_model, 'train_averaged');
for i = 1:n
    if averaged
        src(i).heat_windows = src(i).train_windows;
        src(i).heat_weight  = src(i).pd / src(i).pri;
    else
        src(i).heat_windows = src(i).pulse_windows;
        src(i).heat_weight  = 1;
    end
end

%% Breakpoints and segments
all_pulse = vertcat(src.heat_windows);
all_train = vertcat(src.train_windows);
t_son = max([timing.ptrd; all_train(:, 2)]);
if strcmp(mode, 'custom')
    % Each transducer's protocol spans its own PTRD from its start, as a single
    % run spans timing.ptrd, so a sequence ends where its last run would.
    t_son = max([t_son; [src.offset_s]' + [src.ptrd]']);
end
post  = 0;
if isfield(timing, 'post_ptri_dur') && ~isempty(timing.post_ptri_dur)
    post = max(0, timing.post_ptri_dur);
end
t_end = t_son + post;

tol = min(1e-9 * max(t_end, 1), 1e-3 * min([src.pd]));
b = [0; t_son; t_end; all_pulse(:); all_train(:)];
b = b(b >= 0 & b <= t_end);
b = uniquetol(sort(b), tol, 'DataScale', 1);
t0  = b(1:end-1);
dur = diff(b);
mid = t0 + dur / 2;

active = false(numel(t0), n);
in_train = false(numel(t0), 1);
for i = 1:n
    active(:, i) = inside(mid, src(i).heat_windows);
    in_train = in_train | inside(mid, src(i).train_windows);
end
fine = in_train & mid < t_son;

%% Overlap between sources
% Sources of one group (same SourceTx) fire together by design; only sources of
% different transducers or groups count as overlapping.
n_groups_on = zeros(numel(t0), 1);
for row = 1:numel(t0)
    n_groups_on(row) = numel(unique(source_tx(active(row, :))));
end
multi = n_groups_on > 1;
if any(multi) && ~strcmp(mode, 'simultaneous')   % intended, at different frequencies
    pairs = unique(active(multi, :), 'rows');
    for row = 1:size(pairs, 1)
        on = find(pairs(row, :));
        f = source_freq(on);
        if numel(unique(f(~isnan(f)))) < numel(f(~isnan(f))) || all(isnan(f))
            error('thermal_schedule:coherentOverlap', ...
                ['Transducers %s are on at the same time at the same frequency. Their fields ' ...
                 'interfere, which async coupling (summed intensities) does not model. Give them ' ...
                 'the same offset and timing to fire them together as one group, change the ' ...
                 'offsets, or use simulation.transducer_coupling: coherent.'], ...
                mat2str(unique(source_tx(on))));
        end
    end
    warning('thermal_schedule:overlap', ...
        ['Transducers at different frequencies are on at the same time for %.4g s in total; ' ...
         'their intensities are summed (valid when the pulse is long against 1/|f1-f2|). ' ...
         'The peak pressure then reaches the sum of their amplitudes at the beat maxima, ' ...
         'more than the async p_max_all (per-voxel maximum) and MI report.'], ...
        sum(dur(multi)));
end

%% Time steps
eq_steps = isfield(timing, 'equal_step_duration') && ~isempty(timing.equal_step_duration) && ...
    timing.equal_step_duration;
n_steps = ones(numel(t0), 1);
if eq_steps
    n_steps(fine) = max(1, round(dur(fine) / timing.pt_timestep));
end
% Coarse segments (and train-averaged trains): one recorded segment per post_pt_timestep step
n_rec = ones(numel(t0), 1);
coarse = ~fine | averaged;
n_rec(coarse) = max(1, round(dur(coarse) / timing.post_pt_timestep));
idx = repelem((1:numel(t0))', n_rec);
first = cumsum([1; n_rec(1:end-1)]);
kk = (1:numel(idx))' - first(idx);                     % 0-based step within its segment
sub_dur = dur(idx) ./ n_rec(idx);

sched = struct();
sched.mode         = mode;
sched.n_sources    = n;
sched.t0           = t0(idx) + kk .* sub_dur;
sched.dur          = sub_dur;
sched.n_steps      = n_steps(idx);
sched.step_dur     = sched.dur ./ sched.n_steps;
sched.active       = active(idx, :);
sched.weight       = double(sched.active) .* [src.heat_weight];   % heat-source weight per source
sched.pulse_model  = pulse_model;
sched.fine         = fine(idx);
sched.t_sonication = t_son;
sched.t_end        = t_end;
sched.total_steps  = sum(sched.n_steps);
sched.source       = src;
end

function tf = inside(x, w)
% x within any of the sorted, non-overlapping windows w = [start end].
if isempty(w)
    tf = false(size(x));
    return;
end
n_started = discretize(x, [-Inf; w(:, 1); Inf]) - 1;
n_ended   = discretize(x, [-Inf; w(:, 2); Inf]) - 1;
tf = n_started > n_ended;
end

function f = tx_freq(parameters, j)
f = NaN;
if isfield(parameters.transducer(j), 'freq_hz') && ~isempty(parameters.transducer(j).freq_hz)
    f = double(parameters.transducer(j).freq_hz);
end
end
