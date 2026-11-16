function checks = prestus_check_install(options)
% PRESTUS_CHECK_INSTALL  Report what a PRESTUS installation has and what it lacks
%
% Checks the parts PRESTUS needs for every simulation (MATLAB, Java, the
% PRESTUS path, the required submodules) and the optional ones (SimNIBS for
% head simulations, a GPU, k-Wave C++ binaries, the Nexus Python), and
% prints one line per check with how to fix what is missing. Needs no
% config; CHECK_TOOL_AVAILABILITY checks a loaded config before a run.
%
% Use as:
%   prestus_check_install()
%   checks = prestus_check_install(Quiet=true)
%
% Name-value options:
%   Quiet - (1,1) logical, do not print the report (default false)
%   Nexus - (1,1) logical, look for a Python that starts the Nexus; this
%           starts Python once per candidate (default true)
%
% Output:
%   checks - struct array with fields
%            .name      what was checked
%            .required  true if every simulation needs it
%            .ok        true if present
%            .detail    what was found
%            .fix       how to install it ('' when ok)
%            .submodule external/ folder of a missing required submodule,
%                       '' otherwise (PRESTUS_SETUP fetches these)
%
% See also: PRESTUS_SETUP, CHECK_TOOL_AVAILABILITY, FIND_NEXUS_PYTHON

arguments
    options.Quiet (1,1) logical = false
    options.Nexus (1,1) logical = true
end

    root = get_prestus_path();
    checks = struct('name', {}, 'required', {}, 'ok', {}, 'detail', {}, ...
        'fix', {}, 'submodule', {});

    % ---- required ----
    rel = version('-release');
    checks(end + 1) = make('MATLAB', true, true, ['R' rel], '');
    if ~strcmp(rel, '2023b')
        checks(end).detail = sprintf('R%s (PRESTUS is tested on R2023b)', rel);
    end

    checks(end + 1) = make('Java (JVM)', true, usejava('jvm'), 'MATLAB''s Java', ...
        'Start MATLAB without -nojvm; the YAML loader and external commands need Java.');

    where = which('load_parameters');
    on_path = startsWith(where, fullfile(root, 'functions'));
    detail = root;
    if ~on_path && ~isempty(where)
        detail = sprintf('another copy comes first: %s', where);
    elseif ~on_path
        detail = 'not on the MATLAB path';
    end
    checks(end + 1) = make('PRESTUS on path', true, on_path, detail, ...
        sprintf('Run prestus_setup in %s.', root));

    subs = required_submodules();
    for k = 1:size(subs, 1)
        ok = isfile(fullfile(root, 'external', subs{k, 2}));
        chk = make(subs{k, 1}, true, ok, ['external/' subs{k, 3}], ...
            sprintf(['Run prestus_setup (fetches it with git), or download %s ' ...
                     'into external/%s.'], subs{k, 4}, subs{k, 3}));
        if ~ok
            chk.detail = sprintf('external/%s is empty', subs{k, 3});
            chk.submodule = ['external/' subs{k, 3}];
        end
        checks(end + 1) = chk; %#ok<AGROW>
    end

    % ---- optional ----
    bins = simnibs_bins();
    if isempty(bins)
        checks(end + 1) = make('SimNIBS 4 (head simulations)', false, false, ...
            'no charm found', ['Install SimNIBS 4 (doc_installation.md#install-simnibs), ' ...
            'then set startup.simnibs_bin_path to its bin folder.']);
    else
        checks(end + 1) = make('SimNIBS 4 (head simulations)', false, true, ...
            sprintf('charm in %s; set startup.simnibs_bin_path: ''%s''', ...
            strjoin(bins, ', '), bins{1}), '');
    end

    [n_gpu, gpu_detail] = gpu_info();
    checks(end + 1) = make('Nvidia GPU (matlab_gpu, cpp_gpu)', false, n_gpu > 0, ...
        gpu_detail, 'Simulations run on the CPU (code_type matlab_cpu) without one.');

    bin_dir = fullfile(root, 'external', 'k-wave', 'k-Wave', 'binaries');
    names = {'kspaceFirstOrder-OMP', 'kspaceFirstOrder-CUDA', 'kspaceFirstOrder-MLX'};
    found = names(cellfun(@(n) isfile(fullfile(bin_dir, n)) || ...
        isfile(fullfile(bin_dir, [n '.exe'])), names));
    detail = 'none';
    if ~isempty(found)
        detail = strjoin(found, ', ');
    end
    checks(end + 1) = make('k-Wave C++ binaries (cpp_* code types)', false, ...
        ~isempty(found), detail, 'Only for code_type cpp_cpu/cpp_gpu/cpp_mlx: see doc_kwave.md.');

    if options.Nexus
        python = find_nexus_python();
        detail = python;
        if isempty(python)
            detail = 'no Python with PyQt6 found';
        end
        if ispc
            fix = sprintf('Run %s once.', fullfile(root, 'nexus', 'install.bat'));
        else
            fix = sprintf('Run  bash %s  once.', fullfile(root, 'nexus', 'install.sh'));
        end
        checks(end + 1) = make('PRESTUS Nexus', false, ~isempty(python), detail, fix);
    end

    if ~options.Quiet
        print_report(checks, root);
    end
    if nargout == 0
        clear checks
    end
end

function subs = required_submodules()
% {name, file inside external/, folder in external/, source}
    subs = {
        'k-Wave',       fullfile('k-wave', 'k-Wave', 'kspaceFirstOrder3D.m'), 'k-wave', ...
                        'https://github.com/ucl-bug/k-wave (v1.4.1)'
        'export_fig',   fullfile('export_fig', 'export_fig.m'),     'export_fig', ...
                        'https://github.com/altmany/export_fig'
        'FEX-minimize', fullfile('FEX-minimize', 'minimize.m'),     'FEX-minimize', ...
                        'https://github.com/rodyo/FEX-minimize'
        'xml2struct',   fullfile('xml2struct', 'xml2struct.m'),     'xml2struct', ...
                        'https://github.com/joe-of-all-trades/xml2struct'
        };
end

function bins = simnibs_bins()
% Folders holding charm: PATH, conda environments and the SimNIBS
% installers' default folders (as candidate_bins in nexus/prestus_nexus/simnibs.py).
    home = char(java.lang.System.getProperty('user.home'));
    cands = strsplit(getenv('PATH'), pathsep);
    envs = conda_env_dirs();
    if ispc
        cands = [cands, cellfun(@(e) fullfile(e, 'Scripts'), envs, 'UniformOutput', false)];
        pats = {fullfile(home, 'SimNIBS*', 'bin'), ...
                fullfile(getenv('LOCALAPPDATA'), 'SimNIBS*', 'bin'), ...
                fullfile(getenv('PROGRAMFILES'), 'SimNIBS*', 'bin')};
        charm = {'charm.cmd', 'charm.exe', 'charm.bat', 'charm'};
    else
        cands = [cands, cellfun(@(e) fullfile(e, 'bin'), envs, 'UniformOutput', false)];
        pats = {fullfile(home, 'SimNIBS*', 'bin'), ...
                fullfile(home, 'Applications', 'SimNIBS*', 'bin'), ...
                '/Applications/SimNIBS*/bin', '/opt/SimNIBS*/bin', '/opt/simnibs*/bin', ...
                '/opt/simnibs/*/bin', '/opt/simnibs/*/*/bin', '/usr/local/SimNIBS*/bin'};
        charm = {'charm'};
    end
    for k = 1:numel(pats)
        hits = dir(fileparts(pats{k}));        % the SimNIBS* folders
        for j = 1:numel(hits)
            if hits(j).isdir && ~startsWith(hits(j).name, '.')
                cands{end + 1} = fullfile(hits(j).folder, hits(j).name, 'bin'); %#ok<AGROW>
            end
        end
    end
    has_charm = @(d) any(cellfun(@(c) isfile(fullfile(d, c)), charm));
    cands = cands(~cellfun(@isempty, cands));
    [~, keep] = unique(cands, 'stable');
    cands = cands(sort(keep));
    bins = cands(cellfun(has_charm, cands));
end

function [n, detail] = gpu_info()
    n = 0;
    if ~license('test', 'Distrib_Computing_Toolbox') || ~exist('gpuDeviceCount', 'file')
        detail = 'Parallel Computing Toolbox not available';
        return
    end
    try
        n = gpuDeviceCount('available');
    catch
        n = 0;
    end
    if n == 0
        detail = 'no supported GPU on this computer';
    else
        g = gpuDevice;
        detail = sprintf('%d available (%s)', n, g.Name);
    end
end

function c = make(name, required, ok, detail, fix)
    if ok
        fix = '';
    end
    c = struct('name', name, 'required', required, 'ok', logical(ok), ...
        'detail', detail, 'fix', fix, 'submodule', '');
end

function print_report(checks, root)
    fprintf('\nPRESTUS installation check (%s)\n\n', root);
    groups = {true, 'Required for every simulation'; false, 'Optional'};
    for g = 1:2
        fprintf('%s\n', groups{g, 2});
        sel = checks([checks.required] == groups{g, 1});
        for c = sel
            if c.ok
                mark = 'OK  ';
            elseif c.required
                mark = 'MISS';
            else
                mark = ' -  ';
            end
            fprintf('  [%s] %-40s %s\n', mark, c.name, c.detail);
            if ~isempty(c.fix)
                fprintf('         %s\n', c.fix);
            end
        end
        fprintf('\n');
    end
    missing = checks([checks.required] & ~[checks.ok]);
    if isempty(missing)
        fprintf('Everything a free-water simulation needs is installed.\n\n');
    else
        fprintf('%d required part(s) missing; see above.\n\n', numel(missing));
    end
end
