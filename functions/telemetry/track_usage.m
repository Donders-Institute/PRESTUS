function track_usage(event, parameters, options, pipe_options)
% TRACK_USAGE  Send anonymous usage statistics to the PRESTUS telemetry endpoint.
%
% Sends only the levels the user shares (TELEMETRY_LEVELS: chosen via
% telemetry_setup or telemetry_set_consent, or PRESTUS_TELEMETRY for this
% process); nothing when none is shared. PRESTUS_TELEMETRY=off disables
% sending even after a recorded opt-in (see TELEMETRY_ENV_OVERRIDE).
% All network errors are silently swallowed – this function never interrupts
% the pipeline. No personal data is transmitted; see doc/doc_telemetry.md.
%
% Use as:
%   track_usage('run_end',   parameters, struct('run_id', rid, 'duration_s', t, 'status', 'success'))
%   track_usage('run_error', parameters, struct('run_id', rid, 'duration_s', t, 'error_id', me.identifier))
%
% Events:
%   run_end    - emitted on clean pipeline exit
%   run_error  - emitted from the pipeline catch block
%
% A per-run UUID is generated in prestus_pipeline via generate_run_id() and
% passed as options.run_id so individual runs are identifiable on the backend.
%
% pipe_options is the pipeline options struct (e.g. with sequential_configs)
% used for pipeline mode detection; it is optional and never transmitted.
%
% See also: TELEMETRY_SETUP, TELEMETRY_LEVELS, TELEMETRY_LEVEL_INFO, TELEMETRY_ENV_OVERRIDE, PRESTUS_VERSION

    arguments
        event        (1,:) char
        parameters   (1,1) struct
        options      (1,1) struct = struct()
        pipe_options (1,1) struct = struct()
    end

    [levels, notice] = telemetry_levels();
    if isempty(levels)
        return
    end

    try
        payload = build_payload(event, parameters, options, pipe_options, levels, notice);
        send_payload(payload);
    catch
        % silently drop – telemetry must never interrupt the pipeline
    end
end

% =========================================================================
%% Internal helpers
% =========================================================================

function payload = build_payload(event, parameters, options, pipe_options, levels, notice)
% Every record carries the envelope; each shared level adds its fields
% (FIELD_GROUPS). A field in several groups is sent when any of them is shared.

    has = @(level) any(strcmp(levels, level));
    p = struct();

    % --- envelope (all anonymous) ---
    p.event = event;
    % event format version; raise it when a field's meaning or generation
    % changes. 2: run_id and uuid drawn independently of MATLAB's rand.
    % 3: fields grouped by consent level; outcomes added.
    p.telemetry_schema = 3;
    p.consent_levels   = levels;
    p.consent_notice   = notice;
    p.run_type         = telemetry_run_type();
    p.uuid             = get_or_create_uuid();
    [p.prestus_hash, p.prestus_ver] = prestus_version();
    if isfield(options, 'run_id'),       p.run_id       = options.run_id;       end
    if isfield(options, 'status'),       p.status       = options.status;       end

    % --- usage and errors ---
    % development runs can be filtered out: run_type 'development'
    % (PRESTUS_RUN_TYPE), or a checkout off main or with local changes. The
    % provenance of this run (PATH_LOG_SETUP) is reused: no second git call.
    info = safe_get(parameters, {'io', 'prestus_version'}, struct());
    p.prestus_branch = branch_class(info);
    p.prestus_dirty  = safe_get(info, {'dirty'}, 'unknown');
    p.matlab_ver     = version('-release');
    p.platform       = computer('arch');
    if has('basic')
        try
            ci = getComputerInfo();
            p.kwave_ver = ci.kwave_version;
        catch
            p.kwave_ver = 'unknown';
        end
    end
    p.sim_platform     = safe_get(parameters, {'platform'}, 'unknown');
    p.code_type        = safe_get(parameters, {'simulation', 'code_type'}, 'unknown');
    p.acoustic_backend = safe_get(parameters, {'acoustic_backend'}, 'kwave');
    p.thermal_backend  = safe_get(parameters, {'thermal_backend'}, 'kwave');
    p.medium           = safe_get(parameters, {'simulation', 'medium'}, 'unknown');
    p.pct_enabled      = logical(safe_get(parameters, {'pct', 'enabled'}, false));
    p.pct_density      = safe_get(parameters, {'pct', 'mapping_density'}, '');
    p.pct_speed        = safe_get(parameters, {'pct', 'mapping_soundspeed'}, '');
    p.pct_atten        = safe_get(parameters, {'pct', 'mapping_attenuation'}, '');
    p.pipeline_mode    = detect_pipeline_mode(parameters, pipe_options);
    p.modules_enabled  = modules_on(parameters);
    p.acoustic_fov     = ~isempty(safe_get(parameters, {'grid', 'acoustic_fov_diameter_mm'}, [])) || ...
                         ~isempty(safe_get(parameters, {'grid', 'acoustic_fov_length_mm'}, []));
    p.headmodel_reuse  = logical(safe_get(parameters, {'headmodel', 'reuse_across_placements'}, false));
    % transducer: model and type only, no geometry
    p.transducer_type       = safe_get(parameters, {'transducer', 'type'}, 'unknown');
    p.transducer_model      = transducer_model(parameters);
    p.freq_hz               = safe_get(parameters, {'transducer', 'freq_hz'}, NaN);
    p.n_transducer_elements = get_elem_n(parameters);
    p.n_transducers         = numel(safe_get(parameters, {'transducer'}, struct([])));
    p.transducer_coupling   = char(safe_get(parameters, {'simulation', 'transducer_coupling'}, ''));
    if isfield(options, 'duration_s'),   p.duration_s   = options.duration_s;   end
    if isfield(options, 'error_id'),     p.error_id     = options.error_id;     end
    if isfield(options, 'failed_stage'), p.failed_stage = options.failed_stage; end

    % --- performance (Nexus run-time estimates): hardware, resources, grid ---
    if has('performance')
        p.cpu_model = get_cpu_model();
        p.ram_gb    = get_ram_gb();
        p.gpu_model = get_gpu_model();
    end
    p.precision             = safe_get(parameters, {'simulation', 'precision'}, 'unknown');
    p.hpc_timelimit_h       = timelimit_hours(safe_get(parameters, {'hpc', 'timelimit'}, ''));
    p.hpc_memory_gb         = safe_get(parameters, {'hpc', 'memorylimit'}, NaN);
    % partition class only: a partition name can identify the institution
    p.hpc_partition         = partition_class(safe_get(parameters, {'hpc', 'partition'}, ''));
    p.n_cpu_alloc           = str2double(getenv('SLURM_CPUS_PER_TASK'));   % NaN outside SLURM
    p.resolution_mm         = safe_get(parameters, {'grid', 'resolution_mm'}, NaN);
    p.fit_alpha_power       = safe_get(parameters, {'fit_alpha_power'}, 1);
    p.source_cfl            = safe_get(parameters, {'grid', 'source_cfl'}, NaN);
    p.source_ppw            = safe_get(parameters, {'grid', 'source_ppw'}, NaN);
    p.steady_state_transits = safe_get(parameters, {'grid', 'steady_state_transits'}, NaN);
    p.sim_time_us           = safe_get(parameters, {'grid', 'sim_time_us'}, NaN);
    rec = struct();
    if isfield(options, 'sim_stats') && isstruct(options.sim_stats)
        rec = options.sim_stats;
        p.sim_stats = sim_stats_fields(rec, levels);
    end

    % --- simulation outcomes ---
    if has('outcomes')
        p.outcomes = outcome_fields(parameters, rec);
    end

    % keep the envelope and the fields of the shared levels
    groups = field_groups();
    keep = groups.envelope;
    for k = 1:numel(levels)
        if isfield(groups, levels{k})
            keep = [keep, groups.(levels{k})]; %#ok<AGROW>
        end
    end
    payload = struct();
    fn = fieldnames(p);
    for k = 1:numel(fn)
        if any(strcmp(keep, fn{k}))
            payload.(fn{k}) = p.(fn{k});
        end
    end
end

% -------------------------------------------------------------------------
function g = field_groups()
% Top-level payload fields per consent level (TELEMETRY_LEVEL_INFO). The
% solver, medium and transducer fields sit in several levels: each level's
% numbers need them to be read.
    g.envelope = {'event', 'telemetry_schema', 'consent_levels', 'consent_notice', ...
                  'run_type', 'uuid', 'prestus_hash', 'prestus_ver', 'run_id', 'status'};
    shared = {'code_type', 'acoustic_backend', 'medium'};
    g.basic = [shared, {'prestus_branch', 'prestus_dirty', 'matlab_ver', 'platform', ...
               'kwave_ver', 'sim_platform', 'thermal_backend', 'pct_enabled', 'pct_density', ...
               'pct_speed', 'pct_atten', 'pipeline_mode', 'modules_enabled', 'acoustic_fov', ...
               'headmodel_reuse', 'transducer_type', 'transducer_model', 'freq_hz', ...
               'n_transducer_elements', 'n_transducers', 'transducer_coupling', ...
               'duration_s', 'error_id', 'failed_stage'}];
    g.performance = [shared, {'cpu_model', 'ram_gb', 'gpu_model', 'precision', ...
               'hpc_timelimit_h', 'hpc_memory_gb', 'hpc_partition', 'n_cpu_alloc', ...
               'resolution_mm', 'fit_alpha_power', 'source_cfl', 'source_ppw', ...
               'steady_state_transits', 'sim_time_us', 'sim_stats', 'duration_s'}];
    g.outcomes = [shared, {'thermal_backend', 'pct_enabled', 'transducer_type', ...
               'transducer_model', 'freq_hz', 'n_transducers', 'resolution_mm', ...
               'sim_stats', 'outcomes'}];
    g.skull = {'sim_stats'};      % notice-2 records: the skull medians only
end

% -------------------------------------------------------------------------
function o = outcome_fields(parameters, rec)
% Eight rounded summary statistics of the run (first transducer, first
% analysed field): brain intensity, its ratio to free water, brain MI,
% targeting error, peak brain and skull temperature rise, and the two
% protocol numbers that the temperature rise depends on (overall duty cycle
% and sonication duration). No maps, positions or target names; values a
% run did not compute are left out.
    sig = @(x) round(double(x), 2, 'significant');
    o = struct();
    map = {'outcome_isppa_brain_wcm2', 'isppa_brain_wcm2'; ...
           'outcome_mi_brain',         'mi_brain'; ...
           'outcome_focal_error_mm',   'focal_error_mm'};
    for k = 1:size(map, 1)
        v = safe_get(rec, map(k, 1), []);
        if isnumeric(v) && isscalar(v) && isfinite(v)
            o.(map{k, 2}) = sig(v);
        end
    end
    ff = safe_get(rec, {'outcome_freefield_isppa_wcm2'}, []);
    if isfield(o, 'isppa_brain_wcm2') && isnumeric(ff) && isscalar(ff) && isfinite(ff) && ff > 0
        o.brain_transmission = sig(rec.outcome_isppa_brain_wcm2 / ff);
    end
    for f = {'brain', 'skull'}
        v = safe_get(rec, {['outcome_riset_' f{1} '_c']}, []);
        if isnumeric(v) && isscalar(v) && isfinite(v)
            o.(['riset_' f{1} '_c']) = round(double(v), 1);
        end
    end
    [duty, dur] = protocol_summary(safe_get(parameters, {'timing'}, struct()));
    if isfinite(duty), o.duty_cycle = sig(duty); end
    if isfinite(dur),  o.sonication_s = sig(dur); end
end

% -------------------------------------------------------------------------
function [duty, dur] = protocol_summary(t)
% Overall duty cycle (pulse duty x train duty) and total sonication time [s]
% from timing.*; NaN when the protocol is not set.
    num = @(f) double(safe_get(t, {f}, 0));
    pd = num('pd'); pri = num('pri'); ptd = num('ptd'); ptri = num('ptri'); ptrd = num('ptrd');
    duty = NaN;
    dur = NaN;
    if ~(isscalar(pd) && pd > 0), return, end
    duty = 1;
    if isscalar(pri) && pri > 0, duty = min(pd / pri, 1); end
    if isscalar(ptd) && ptd > 0 && isscalar(ptri) && ptri > 0, duty = duty * min(ptd / ptri, 1); end
    if isscalar(ptrd) && ptrd > 0
        dur = ptrd;
    elseif isscalar(ptd) && ptd > 0
        dur = ptd;
    end
end

% -------------------------------------------------------------------------
function names = modules_on(parameters)
    mods = safe_get(parameters, {'modules'}, struct());
    fn = fieldnames(mods);
    names = {};
    for i = 1:numel(fn)
        if isequal(mods.(fn{i}), 1) || isequal(mods.(fn{i}), true)
            names{end+1} = fn{i}; %#ok<AGROW>
        end
    end
end

% -------------------------------------------------------------------------
function c = branch_class(info)
% 'main', 'development' or 'other': the branch class, never the branch name.
    b = char(safe_get(info, {'branch'}, ''));
    if any(strcmp(b, {'main', 'master'}))
        c = 'main';
    elseif strcmp(b, 'development')
        c = 'development';
    elseif isempty(b) || strcmp(b, 'unknown')
        c = 'unknown';
    else
        c = 'other';
    end
end

% -------------------------------------------------------------------------
function h = timelimit_hours(t)
% SLURM/PBS wall time 'D-HH:MM:SS', 'HH:MM:SS', 'MM:SS' or minutes -> hours (NaN if unset).
    h = NaN;
    if isnumeric(t) && isscalar(t)
        h = double(t) / 60;
        return
    end
    t = strtrim(char(t));
    if isempty(t)
        return
    end
    days = 0;
    parts = split(t, '-');
    if numel(parts) == 2
        days = str2double(parts{1});
        t = parts{2};
    end
    v = str2double(split(t, ':'))';
    if any(isnan(v)) || isnan(days)
        return
    end
    v = [zeros(1, 3 - numel(v)), v];           % pad to [h m s]
    h = days * 24 + v(1) + v(2) / 60 + v(3) / 3600;
end

% -------------------------------------------------------------------------
function uuid = get_or_create_uuid()
% The anonymous install ID in <prefdir_prestus>/uuid.txt, created on first use.
% Clients before 2026-11 drew it from MATLAB's rand at its default seed, so
% many installs share the same ID (COLLIDED_ID); such a file gets a new one.
    uuid_file = fullfile(prefdir_prestus(), 'uuid.txt');
    if isfile(uuid_file)
        fid  = fopen(uuid_file, 'r');
        uuid = strtrim(fgetl(fid));
        fclose(fid);
        if ischar(uuid) && ~isempty(uuid) && ~strcmp(uuid, collided_id())
            return
        end
    end
    uuid = generate_run_id();
    fid  = fopen(uuid_file, 'w');
    fprintf(fid, '%s\n', uuid);
    fclose(fid);
end

function id = collided_id()
% The first rand(1,16) draw of a fresh MATLAB session as a UUID: the ID that
% earlier clients gave every install and run created at the default seed.
    id = 'd0e720e9-a118-478c-b5f7-28f8f57ccc24';
end

% -------------------------------------------------------------------------
function val = safe_get(s, keys, default)
% Navigate nested struct with a cell array of field names; return default
% if any level is missing or not a struct.
    if nargin < 3, default = []; end
    val = s;
    for i = 1:numel(keys)
        if isstruct(val) && isfield(val, keys{i})
            val = val.(keys{i});
        else
            val = default;
            return
        end
    end
end

% -------------------------------------------------------------------------
function mode = detect_pipeline_mode(parameters, pipe_options)
    if nargin < 2, pipe_options = struct(); end
    if is_async_mode(parameters)
        mode = 'async';
    elseif is_multi_isppa_mode(parameters)
        mode = 'multi_isppa';
    elseif isfield(pipe_options, 'sequential_configs') || ...
           (isfield(pipe_options, 'is_sequential_run') && pipe_options.is_sequential_run)
        mode = 'sequential';
    else
        mode = 'single';
    end
end

% -------------------------------------------------------------------------
function n = get_elem_n(parameters)
    tr_type = safe_get(parameters, {'transducer', 'type'}, '');
    switch tr_type
        case 'annular'
            n = safe_get(parameters, {'transducer', 'annular', 'elem_n'}, NaN);
        case 'matrix'
            % matrix element count lives under the chosen arrangement sub-struct
            arrangements = {'fibonacci', 'fermat', 'file'};
            n = NaN;
            for k = 1:numel(arrangements)
                v = safe_get(parameters, {'transducer', 'matrix', 'elements', arrangements{k}, 'elem_n'}, NaN);
                if ~isnan(v), n = v; return; end
            end
            % rectangular stores row × col
            nr = safe_get(parameters, {'transducer', 'matrix', 'elements', 'rectangular', 'elem_n_row'}, NaN);
            nc = safe_get(parameters, {'transducer', 'matrix', 'elements', 'rectangular', 'elem_n_col'}, NaN);
            if ~isnan(nr) && ~isnan(nc), n = nr * nc; end
        otherwise
            n = NaN;
    end
end

% -------------------------------------------------------------------------
function m = get_cpu_model()
    try
        arch = computer('arch');
        if contains(arch, 'mac')
            [~, m] = run_command({'sysctl', '-n', 'machdep.cpu.brand_string'});
        elseif isunix
            % /proc/cpuinfo is a regular (virtual) file -- read it
            % directly instead of `grep ... | cut ...`, no process at all.
            cpuinfo = fileread('/proc/cpuinfo');
            line_match = regexp(cpuinfo, '^model name\s*:\s*(.*)$', 'once', 'tokens', 'lineanchors');
            if ~isempty(line_match)
                m = line_match{1};
            else
                m = 'unknown';
            end
        else
            [~, m] = run_command({'wmic', 'cpu', 'get', 'name', '/value'});
            m = regexprep(m, 'Name=', '');
        end
        m = strtrim(m);
    catch
        m = 'unknown';
    end
end

% -------------------------------------------------------------------------
function gb = get_ram_gb()
    try
        arch = computer('arch');
        if contains(arch, 'mac')
            [~, r] = run_command({'sysctl', '-n', 'hw.memsize'});
            gb = round(str2double(strtrim(r)) / 1024^3);
        elseif isunix
            % /proc/meminfo is a regular (virtual) file -- read it
            % directly instead of `grep MemTotal /proc/meminfo`.
            meminfo = fileread('/proc/meminfo');
            r = regexp(meminfo, '^MemTotal:.*$', 'once', 'match', 'lineanchors');
            kb = str2double(regexp(r, '\d+', 'match', 'once'));
            gb = round(kb / 1024^2);
        else
            [~, r] = run_command({'wmic', 'computersystem', 'get', 'TotalPhysicalMemory', '/value'});
            bytes = str2double(regexp(r, '\d+', 'match', 'once'));
            gb = round(bytes / 1024^3);
        end
    catch
        gb = NaN;
    end
end

% -------------------------------------------------------------------------
function name = get_gpu_model()
    try
        if ~license('test', 'Distrib_Computing_Toolbox')
            name = 'no_pct';
            return
        end
        g = gpuDevice();
        name = g.Name;
    catch
        name = 'none';
    end
end

% -------------------------------------------------------------------------
function c = partition_class(partition)
% '' (unset), 'gpu' (the name contains "gpu"), or 'other'; never the name.
    partition = char(string(partition));
    if isempty(partition)
        c = '';
    elseif contains(lower(partition), 'gpu')
        c = 'gpu';
    else
        c = 'other';
    end
end

% -------------------------------------------------------------------------
function m = transducer_model(parameters)
% Shipped equipment model of the first transducer (config/equipment/model-*),
% matched by its serial or name; 'custom' for anything else, '' when neither is
% set. Never the name or serial: a serial identifies one device, and a free-text
% name can name a lab or a person.
    m = '';
    t = safe_get(parameters, {'transducer'}, struct());
    if iscell(t) && ~isempty(t), t = t{1}; end
    if ~isstruct(t) || isempty(t), return, end
    t = t(1);
    squash = @(x) lower(regexprep(char(string(x)), '[^A-Za-z0-9]', ''));
    keys = {};
    for f = {'serial', 'name'}
        if isfield(t, f{1}) && (ischar(t.(f{1})) || isstring(t.(f{1}))) && strlength(string(t.(f{1}))) > 0
            keys{end+1} = squash(t.(f{1})); %#ok<AGROW>
        end
    end
    if isempty(keys), return, end
    m = 'custom';
    root  = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    files = dir(fullfile(root, 'config', 'equipment', 'model-*', 'tran-*_transducer.yaml'));
    for k = 1:numel(files)
        [~, folder] = fileparts(files(k).folder);
        model  = extractAfter(folder, 'model-');
        serial = regexprep(files(k).name, '^tran-(.*)_transducer\.yaml$', '$1');
        if any(strcmp(keys, squash(serial))) || any(strcmp(keys, squash(model)))
            m = model;
            return
        end
    end
end

% -------------------------------------------------------------------------
function out = sim_stats_fields(rec, levels)
% Allow-listed numeric fields of the run record (RUN_STATS); anything else is
% dropped. Grid size, time steps, stage times and memory with 'performance';
% the skull property medians (3 significant figures) with 'outcomes', or
% with a notice-2 consent ('skull').
    keep = {};
    if any(strcmp(levels, 'performance'))
        keep = {'acoustic_voxels', 'acoustic_dims', 'acoustic_pml', 'axisymmetric', ...
                'acoustic_nt', 'acoustic_s', 'thermal_voxels', 'thermal_steps', ...
                'thermal_s', 'peak_ram_gb', 'n_cpu', 'gpu_mem_gb', 'queue_wait_s'};
        fn = fieldnames(rec);
        stages = fn(~cellfun(@isempty, regexp(fn, '^stage_[a-z_]+_s$', 'once')));   % stage times
        keep = [keep, stages(:)'];
    end
    skull = {};
    if any(ismember(levels, {'outcomes', 'skull'}))
        skull = {'skull_c_median', 'skull_rho_median', 'skull_atten_median', ...
                 'skull_path_c_median', 'skull_path_rho_median', 'skull_path_atten_median'};
    end
    out = struct();
    for k = 1:numel(keep)
        if isfield(rec, keep{k}) && isnumeric(rec.(keep{k}))
            out.(keep{k}) = double(rec.(keep{k}));
        end
    end
    for k = 1:numel(skull)
        if isfield(rec, skull{k}) && isnumeric(rec.(skull{k}))
            out.(skull{k}) = round(double(rec.(skull{k})), 3, 'significant');
        end
    end
end

% -------------------------------------------------------------------------
function send_payload(payload)
    endpoint = 'https://lruelxwhkjibezkgipme.supabase.co/functions/v1/ingest-event';
    try
        opts = weboptions( ...
            'MediaType',     'application/json', ...
            'Timeout',        4, ...
            'RequestMethod', 'post', ...
            'HeaderFields',  {'Content-Type', 'application/json'});
        webwrite(endpoint, jsonencode(payload), opts);
    catch
        % Fire-and-forget: network errors are silently ignored.
    end
end
