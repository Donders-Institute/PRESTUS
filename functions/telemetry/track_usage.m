function track_usage(event, parameters, options, pipe_options)
% TRACK_USAGE  Send anonymous usage statistics to the PRESTUS telemetry endpoint.
%
% Does nothing unless the user has explicitly opted in via telemetry_setup
% or telemetry_set_consent, or PRESTUS_TELEMETRY=on/1/true is set for this
% process. PRESTUS_TELEMETRY=off/0/false disables sending even after a
% recorded opt-in (see TELEMETRY_ENV_OVERRIDE).
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
% See also: TELEMETRY_SETUP, TELEMETRY_SETUP_RESET, TELEMETRY_ENV_OVERRIDE, PRESTUS_VERSION

    arguments
        event        (1,:) char
        parameters   (1,1) struct
        options      (1,1) struct = struct()
        pipe_options (1,1) struct = struct()
    end

    [opted, notice] = telemetry_opted_in();
    if ~opted
        return
    end

    try
        payload = build_payload(event, parameters, options, pipe_options, notice);
        send_payload(payload);
    catch
        % silently drop – telemetry must never interrupt the pipeline
    end
end

% =========================================================================
%% Internal helpers
% =========================================================================

function [opted, notice] = telemetry_opted_in()
% PRESTUS_TELEMETRY=off always wins; a recorded decision comes next, so a
% recorded opt-out is never overridden; PRESTUS_TELEMETRY=on applies only
% when no decision is recorded. notice: the notice version the consent was
% given to (TELEMETRY_NOTICE_VERSION; 1 for records without it, the current
% version for PRESTUS_TELEMETRY=on).
    env    = telemetry_env_override();
    opted  = false;
    notice = 1;
    if strcmp(env, 'off'), return, end
    cfg_file = fullfile(prefdir_prestus(), 'telemetry.json');
    if isfile(cfg_file)
        try
            cfg   = jsondecode(fileread(cfg_file));
            opted = isfield(cfg, 'opt_in') && isequal(cfg.opt_in, true);
            if isfield(cfg, 'notice') && isnumeric(cfg.notice) && isscalar(cfg.notice)
                notice = double(cfg.notice);
            end
        catch
        end
        return
    end
    opted = strcmp(env, 'on');
    notice = telemetry_notice_version();
end

% -------------------------------------------------------------------------
function payload = build_payload(event, parameters, options, pipe_options, notice)

    % --- identity (all anonymous) ---
    payload.event         = event;
    payload.run_type      = telemetry_run_type();
    payload.timestamp_utc = posixtime(datetime('now', 'TimeZone', 'UTC'));
    payload.uuid          = get_or_create_uuid();
    [payload.prestus_hash, payload.prestus_ver] = prestus_version();
    % development runs can be filtered out: run_type 'development'
    % (PRESTUS_RUN_TYPE), or a checkout off main or with local changes. The
    % provenance of this run (PATH_LOG_SETUP) is reused: no second git call.
    info = safe_get(parameters, {'io', 'prestus_version'}, struct());
    payload.prestus_branch = branch_class(info);
    payload.prestus_dirty  = safe_get(info, {'dirty'}, 'unknown');
    payload.consent_notice = notice;
    payload.matlab_ver    = version('-release');
    payload.platform      = computer('arch');
    try
        info = getComputerInfo();
        payload.kwave_ver = info.kwave_version;
    catch
        payload.kwave_ver = 'unknown';
    end
    payload.os_version    = get_os_version();
    payload.cpu_model     = get_cpu_model();
    payload.n_cpu_cores   = feature('numcores');
    payload.ram_gb        = get_ram_gb();
    payload.gpu_model     = get_gpu_model();
    payload.toolboxes     = get_toolbox_names();

    % --- execution environment ---
    payload.sim_platform  = safe_get(parameters, {'platform'}, 'unknown');
    payload.code_type     = safe_get(parameters, {'simulation', 'code_type'}, 'unknown');
    payload.precision     = safe_get(parameters, {'simulation', 'precision'}, 'unknown');
    payload.hpc_name      = safe_get(parameters, {'hpc', 'name'}, 'unknown');
    payload.use_gpu       = ~isempty(safe_get(parameters, {'hpc', 'gpu'}, ''));
    % requested job resources (scheduler runs): headroom and the limits runs hit
    payload.hpc_timelimit_h = timelimit_hours(safe_get(parameters, {'hpc', 'timelimit'}, ''));
    payload.hpc_memory_gb   = safe_get(parameters, {'hpc', 'memorylimit'}, NaN);
    payload.hpc_gpu         = char(safe_get(parameters, {'hpc', 'gpu'}, ''));
    payload.hpc_partition   = char(safe_get(parameters, {'hpc', 'partition'}, ''));

    % --- simulation configuration ---
    payload.medium        = safe_get(parameters, {'simulation', 'medium'}, 'unknown');
    payload.pct_enabled   = logical(safe_get(parameters, {'pct', 'enabled'}, false));
    payload.pct_density   = safe_get(parameters, {'pct', 'mapping_density'}, '');
    payload.pct_speed     = safe_get(parameters, {'pct', 'mapping_soundspeed'}, '');
    payload.pct_atten     = safe_get(parameters, {'pct', 'mapping_attenuation'}, '');
    payload.layers        = get_layer_names(parameters);
    % solver and grid settings that change run time
    payload.thermal_backend       = safe_get(parameters, {'thermal_backend'}, 'kwave');
    payload.fit_alpha_power       = safe_get(parameters, {'fit_alpha_power'}, 1);
    payload.source_cfl            = safe_get(parameters, {'grid', 'source_cfl'}, NaN);
    payload.source_ppw            = safe_get(parameters, {'grid', 'source_ppw'}, NaN);
    payload.steady_state_transits = safe_get(parameters, {'grid', 'steady_state_transits'}, NaN);
    payload.sim_time_us           = safe_get(parameters, {'grid', 'sim_time_us'}, NaN);
    payload.acoustic_fov          = ~isempty(safe_get(parameters, {'grid', 'acoustic_fov_diameter_mm'}, [])) || ...
                                    ~isempty(safe_get(parameters, {'grid', 'acoustic_fov_length_mm'}, []));
    payload.headmodel_reuse       = logical(safe_get(parameters, {'headmodel', 'reuse_across_placements'}, false));

    % --- pipeline mode ---
    payload.pipeline_mode = detect_pipeline_mode(parameters, pipe_options);

    % --- transducer (model/type only, no geometry) ---
    payload.transducer_type       = safe_get(parameters, {'transducer', 'type'}, 'unknown');
    payload.freq_hz               = safe_get(parameters, {'transducer', 'freq_hz'}, NaN);
    payload.transducer_name       = safe_get(parameters, {'transducer', 'name'}, '');
    payload.n_transducer_elements = get_elem_n(parameters);
    payload.n_transducers         = numel(safe_get(parameters, {'transducer'}, struct([])));
    payload.transducer_coupling   = char(safe_get(parameters, {'simulation', 'transducer_coupling'}, ''));

    % --- modules enabled ---
    mods = safe_get(parameters, {'modules'}, struct());
    mod_fields = fieldnames(mods);
    modules_on = {};
    for i = 1:numel(mod_fields)
        if isequal(mods.(mod_fields{i}), 1) || isequal(mods.(mod_fields{i}), true)
            modules_on{end+1} = mod_fields{i}; %#ok<AGROW>
        end
    end
    payload.modules_enabled = modules_on;

    % --- simulation size and duration (Nexus resource prediction) ---
    % Allow-listed numeric fields of the run_stats record: grid size, time
    % steps, stage durations, peak RAM. Never the host, paths or timestamps.
    payload.resolution_mm    = safe_get(parameters, {'grid', 'resolution_mm'}, NaN);
    payload.acoustic_backend = safe_get(parameters, {'acoustic_backend'}, 'kwave');
    payload.n_cpu_alloc      = str2double(getenv('SLURM_CPUS_PER_TASK'));   % NaN outside SLURM
    if isfield(options, 'sim_stats') && isstruct(options.sim_stats)
        payload.sim_stats = sim_stats_fields(options.sim_stats, notice);
    end

    % --- outcome fields (populated on run_end / run_error) ---
    if isfield(options, 'run_id'),      payload.run_id      = options.run_id;       end
    if isfield(options, 'duration_s'),  payload.duration_s  = options.duration_s;  end
    if isfield(options, 'status'),      payload.status      = options.status;       end
    if isfield(options, 'error_id'),    payload.error_id    = options.error_id;     end
    if isfield(options, 'failed_stage'), payload.failed_stage = options.failed_stage; end
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
    uuid_file = fullfile(prefdir_prestus(), 'uuid.txt');
    if isfile(uuid_file)
        fid  = fopen(uuid_file, 'r');
        uuid = strtrim(fgetl(fid));
        fclose(fid);
        if ~isempty(uuid), return, end
    end
    % Generate a random UUID v4
    bytes = uint8(floor(rand(1,16) * 256));
    bytes(7) = bitor(bitand(bytes(7), uint8(15)), uint8(64));  % version 4
    bytes(9) = bitor(bitand(bytes(9), uint8(63)), uint8(128)); % variant
    uuid = sprintf('%02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x', bytes);
    fid  = fopen(uuid_file, 'w');
    fprintf(fid, '%s\n', uuid);
    fclose(fid);
end

% -------------------------------------------------------------------------
function names = get_layer_names(parameters)
    layers = safe_get(parameters, {'layers'}, struct());
    if isstruct(layers)
        names = fieldnames(layers);
    else
        names = {};
    end
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
function rid = generate_run_id()
    bytes = uint8(floor(rand(1,16) * 256));
    bytes(7) = bitor(bitand(bytes(7), uint8(15)), uint8(64));
    bytes(9) = bitor(bitand(bytes(9), uint8(63)), uint8(128));
    rid = sprintf('%02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x', bytes);
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
function v = get_os_version()
    try
        arch = computer('arch');
        if contains(arch, 'mac')
            [~, v] = run_command({'sw_vers', '-productVersion'});
        elseif isunix
            [~, v] = run_command({'uname', '-r'});
        else
            % 'ver' is a cmd.exe builtin, not a standalone executable --
            % it needs the one shell layer run_shell_script provides.
            [~, v] = run_shell_script('ver');
        end
        v = strtrim(v);
    catch
        v = 'unknown';
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
function names = get_toolbox_names()
    try
        info = ver;
        names = {info.Name};
        names = names(~strcmp(names, 'MATLAB'));
    catch
        names = {};
    end
end

% -------------------------------------------------------------------------
function out = sim_stats_fields(rec, notice)
% Numeric run_stats fields only (see RUN_STATS); anything else is dropped.
% Skull properties: the medians only, and only with consent to notice 2 or later.
    keep = {'acoustic_voxels', 'acoustic_dims', 'acoustic_pml', 'axisymmetric', ...
            'acoustic_nt', 'acoustic_s', 'thermal_voxels', 'thermal_steps', ...
            'thermal_s', 'peak_ram_gb', 'n_cpu', 'gpu_mem_gb', 'queue_wait_s'};
    if nargin > 1 && notice >= 2
        keep = [keep, {'skull_c_median', 'skull_rho_median', 'skull_atten_median', ...
                       'skull_path_c_median', 'skull_path_rho_median', 'skull_path_atten_median'}];
    end
    out = struct();
    if ~isstruct(rec)
        return
    end
    fn = fieldnames(rec);
    stages = fn(~cellfun(@isempty, regexp(fn, '^stage_[a-z_]+_s$', 'once')));   % stage times
    keep = [keep, stages(:)'];
    for k = 1:numel(keep)
        if isfield(rec, keep{k}) && isnumeric(rec.(keep{k}))
            out.(keep{k}) = double(rec.(keep{k}));
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
