function nexus_bridge(bridge_dir, options)
% NEXUS_BRIDGE  Persistent MATLAB session that serves the PRESTUS Nexus
%
% The Nexus (Python/Qt) starts this once with matlab -batch and sends it
% requests as files, so interactive steps (Localite detect/write, transducer
% calibration predict/fit/save) run in MATLAB without a MATLAB start-up per
% click. Protocol, all inside bridge_dir:
%
%   req/<id>.json       request {"id", "fn", "args"}, written by the Nexus
%                       (to a temporary name, then renamed)
%   resp/<id>.json      response {"id", "ok", "result" | "error", "identifier",
%                       "report", "log"}; NaN/Inf are written as null
%   progress/<id>.txt   last progress message of a running request (NEXUS_PROGRESS)
%   ready.json          written once the session can take requests
%   heartbeat           touched by the Nexus; older than HeartbeatTimeout = exit
%   stop                the Nexus asks the session to end
%
% Only the functions in the allow-list can be called (NEXUS_CALIBRATION,
% NEXUS_NEURONAV, NEXUS_PING); arguments are decoded with jsondecode. The
% session ends after IdleTimeout seconds without requests, since it holds a
% MATLAB licence. Long work (batch ingestion, k-Wave) does not go through
% the bridge: the Nexus runs it as a script on its Run tab.
%
% Use as:
%   nexus_bridge(bridge_dir)
%   nexus_bridge(bridge_dir, IdleTimeout=1800, HeartbeatTimeout=60)
%
% Input:
%   bridge_dir - folder shared with the Nexus (created if missing)
%
% Name-value options:
%   IdleTimeout      - seconds without a request before the session ends (default 1800)
%   HeartbeatTimeout - seconds without a Nexus heartbeat before it ends (default 60)
%   Poll             - seconds between checks for requests (default 0.05)
%
% See also: NEXUS_CALIBRATION, NEXUS_NEURONAV, NEXUS_PROGRESS, PRESTUS_NEXUS

arguments
    bridge_dir (1,:) char
    options.IdleTimeout (1,1) double = 1800
    options.HeartbeatTimeout (1,1) double = 60
    options.Poll (1,1) double = 0.05
end

    allowed = {'nexus_calibration', 'nexus_neuronav', 'nexus_ping'};
    root = get_prestus_path();
    % genpath over a network drive takes seconds: add the trees only when missing, and of
    % external/ only what the façades use (YAML, Localite XML, the fexminimize optimiser)
    if isempty(which('load_parameters')); safe_addpath(fullfile(root, 'functions')); end
    if isempty(which('yaml.loadFile'))
        for ext = {'yaml', 'xml2struct', 'FEX-minimize'}
            safe_addpath(fullfile(root, 'external', ext{1}));
        end
    end

    for d = {'', 'req', 'resp', 'progress'}
        f = fullfile(bridge_dir, d{1});
        if ~isfolder(f); mkdir(f); end
    end
    write_json(fullfile(bridge_dir, 'ready.json'), struct('pid', feature('getpid'), ...
        'matlab', version, 'prestus_root', root));
    fprintf('[nexus_bridge] ready in %s\n', bridge_dir);

    last_request = tic;
    hb = fullfile(bridge_dir, 'heartbeat');
    while true
        if isfile(fullfile(bridge_dir, 'stop'))
            fprintf('[nexus_bridge] stop requested\n'); break;
        end
        if toc(last_request) > options.IdleTimeout
            fprintf('[nexus_bridge] idle for %.0f s, exiting\n', options.IdleTimeout); break;
        end
        if isfile(hb)
            info = dir(hb);
            if (now - info.datenum) * 86400 > options.HeartbeatTimeout %#ok<TNOW1>
                fprintf('[nexus_bridge] no heartbeat from the Nexus, exiting\n'); break;
            end
        end
        reqs = dir(fullfile(bridge_dir, 'req', '*.json'));
        if isempty(reqs)
            pause(options.Poll); continue;
        end
        [~, order] = sort({reqs.name});
        for r = reqs(order)'
            f = fullfile(r.folder, r.name);
            try
                txt = fileread(f);
            catch
                continue;               % being renamed; next round
            end
            delete(f);
            serve(bridge_dir, txt, allowed);
        end
        last_request = tic;
    end
    if isfile(hb); delete(hb); end
end

function serve(bridge_dir, txt, allowed)
    id = '';
    try
        req = jsondecode(txt);
        id = char(req.id);
        fn = char(req.fn);
        if ~any(strcmp(fn, allowed))
            error('nexus_bridge:notAllowed', '%s cannot be called through the bridge.', fn);
        end
        args = struct();
        if isfield(req, 'args') && isstruct(req.args); args = req.args; end
        prog = fullfile(bridge_dir, 'progress', [id '.txt']);
        setappdata(0, 'nexus_bridge_progress', prog);
        cleanup = onCleanup(@() clear_progress(prog));
        [log, out] = call(fn, args);
        resp = struct('id', id, 'ok', true, 'result', {out}, 'log', log);
    catch ME
        resp = struct('id', id, 'ok', false, 'error', ME.message, 'identifier', ME.identifier, ...
            'report', getReport(ME, 'extended', 'hyperlinks', 'off'), 'log', '');
    end
    if isempty(id); id = sprintf('invalid_%s', char(java.util.UUID.randomUUID())); resp.id = id; end
    try
        write_json(fullfile(bridge_dir, 'resp', [id '.json']), resp);
    catch ME
        write_json(fullfile(bridge_dir, 'resp', [id '.json']), struct('id', id, 'ok', false, ...
            'error', ['Result could not be encoded: ' ME.message], 'identifier', ME.identifier, ...
            'report', '', 'log', ''));
    end
end

function [log, out] = call(fn, args)
% Run fn(args), capturing what it prints
    out = [];
    log = evalc('out = feval(fn, args);');
end

function clear_progress(f)
    setappdata(0, 'nexus_bridge_progress', '');
    if isfile(f); delete(f); end
end

function write_json(f, s)
    tmp = [f '.tmp'];
    fid = fopen(tmp, 'w', 'n', 'UTF-8');
    assert(fid > 0, 'nexus_bridge:write', 'Cannot write %s', tmp);
    fwrite(fid, jsonencode(s, 'ConvertInfAndNaN', true), 'char');
    fclose(fid);
    movefile(tmp, f, 'f');
end
