function track_usage(event, parameters, options, pipe_options)
% TRACK_USAGE  Send anonymous usage statistics to the PRESTUS telemetry endpoint.
%
% Does nothing unless the user has explicitly opted in via telemetry_setup.
% All network errors are silently swallowed – this function never interrupts
% the pipeline. No personal data is transmitted; see docs/telemetry.md.
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
% See also: TELEMETRY_SETUP, TELEMETRY_SETUP_RESET, PRESTUS_VERSION

    arguments
        event        (1,:) char
        parameters   (1,1) struct
        options      (1,1) struct = struct()
        pipe_options (1,1) struct = struct()
    end

    if ~telemetry_opted_in()
        return
    end

    try
        payload = build_payload(event, parameters, options, pipe_options);
        send_payload(payload);
    catch
        % silently drop – telemetry must never interrupt the pipeline
    end
end

% =========================================================================
%% Internal helpers
% =========================================================================

function opted = telemetry_opted_in()
    cfg_file = fullfile(prefdir_prestus(), 'telemetry.json');
    opted = false;
    if ~isfile(cfg_file), return, end
    try
        cfg   = jsondecode(fileread(cfg_file));
        opted = isfield(cfg, 'opt_in') && isequal(cfg.opt_in, true);
    catch
    end
end

% -------------------------------------------------------------------------
function payload = build_payload(event, parameters, options, pipe_options)

    % --- identity (all anonymous) ---
    payload.event         = event;
    payload.timestamp_utc = posixtime(datetime('now', 'TimeZone', 'UTC'));
    payload.uuid          = get_or_create_uuid();
    [payload.prestus_hash, payload.prestus_ver] = prestus_version();
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

    % --- simulation configuration ---
    payload.medium        = safe_get(parameters, {'simulation', 'medium'}, 'unknown');
    payload.pct_enabled   = logical(safe_get(parameters, {'pct', 'enabled'}, false));
    payload.pct_density   = safe_get(parameters, {'pct', 'mapping_density'}, '');
    payload.pct_speed     = safe_get(parameters, {'pct', 'mapping_soundspeed'}, '');
    payload.pct_atten     = safe_get(parameters, {'pct', 'mapping_attenuation'}, '');
    payload.layers        = get_layer_names(parameters);

    % --- pipeline mode ---
    payload.pipeline_mode = detect_pipeline_mode(parameters, pipe_options);

    % --- transducer (model/type only, no geometry) ---
    payload.transducer_type       = safe_get(parameters, {'transducer', 'type'}, 'unknown');
    payload.freq_hz               = safe_get(parameters, {'transducer', 'freq_hz'}, NaN);
    payload.transducer_name       = safe_get(parameters, {'transducer', 'name'}, '');
    payload.n_transducer_elements = get_elem_n(parameters);

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

    % --- outcome fields (populated on run_end / run_error) ---
    if isfield(options, 'run_id'),      payload.run_id      = options.run_id;       end
    if isfield(options, 'duration_s'),  payload.duration_s  = options.duration_s;  end
    if isfield(options, 'status'),      payload.status      = options.status;       end
    if isfield(options, 'error_id'),    payload.error_id    = options.error_id;     end
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
    elseif isfield(pipe_options, 'sequential_configs')
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
            [~, v] = system('sw_vers -productVersion');
        elseif isunix
            [~, v] = system('uname -r');
        else
            [~, v] = system('ver');
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
            [~, m] = system('sysctl -n machdep.cpu.brand_string');
        elseif isunix
            [~, m] = system('grep -m1 "model name" /proc/cpuinfo | cut -d: -f2');
        else
            [~, m] = system('wmic cpu get name /value');
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
            [~, r] = system('sysctl -n hw.memsize');
            gb = round(str2double(strtrim(r)) / 1024^3);
        elseif isunix
            [~, r] = system('grep MemTotal /proc/meminfo');
            kb = str2double(regexp(r, '\d+', 'match', 'once'));
            gb = round(kb / 1024^2);
        else
            [~, r] = system('wmic computersystem get TotalPhysicalMemory /value');
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
