function prestus_setup(options)
% PRESTUS_SETUP  Put PRESTUS on the MATLAB path and check the installation
%
% Run once per MATLAB session (or once with Save=true). Adds functions/,
% config/ and external/ of this PRESTUS folder to the MATLAB path, from any
% working directory; fetches required submodules that are missing (k-Wave,
% export_fig, FEX-minimize, xml2struct) with git; and prints what is
% installed and what is missing (PRESTUS_CHECK_INSTALL). Optional
% submodules (BabelBrain, PlanTUS) are left alone.
%
% Use as:
%   prestus_setup                    % in the PRESTUS folder
%   run('/path/to/PRESTUS/prestus_setup')
%   prestus_setup(Save=true)         % also keep the path for later sessions
%
% Name-value options:
%   Save  - (1,1) logical, save the path with SAVEPATH (default false)
%   Fetch - (1,1) logical, fetch missing required submodules (default true)
%   Check - (1,1) logical, print the installation report (default true)
%
% See also: PRESTUS_CHECK_INSTALL, PRESTUS_NEXUS, LOAD_PARAMETERS

arguments
    options.Save (1,1) logical = false
    options.Fetch (1,1) logical = true
    options.Check (1,1) logical = true
end

    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, 'functions', 'helper'));
    add_paths(root);

    checks = prestus_check_install(Quiet=true, Nexus=false);
    missing = {checks(~cellfun(@isempty, {checks.submodule})).submodule};
    if options.Fetch && ~isempty(missing)
        fetch_submodules(root, missing);
        add_paths(root);
    end

    if options.Save
        save_path(root);
    end

    if options.Check
        prestus_check_install();
        if ~options.Save
            fprintf(['The path holds for this MATLAB session; prestus_setup(Save=true) ' ...
                     'keeps it.\nNext: doc/doc_first-simulation.md (First Simulation).\n']);
        end
    end
end

function save_path(root)
% The default pathdef.m lives in matlabroot, which needs admin rights on
% Windows and on shared installs; MATLAB also reads a pathdef.m in the
% userpath (startup) folder, so fall back to that.
    if savepath() == 0
        fprintf('MATLAB path saved; PRESTUS is on the path in new sessions.\n');
        return
    end
    up = userpath();
    if ~isempty(up) && isfolder(up) && savepath(fullfile(up, 'pathdef.m')) == 0
        fprintf('MATLAB path saved to %s; PRESTUS is on the path in new sessions.\n', ...
            fullfile(up, 'pathdef.m'));
        return
    end
    warning('prestus:savepath', ['Could not save the MATLAB path (no write access ' ...
        'to pathdef.m). Add  run(''%s'')  to your startup.m instead.'], ...
        fullfile(root, 'prestus_setup'));
end

function add_paths(root)
    addpath(root);                     % prestus_setup itself, from any folder
    for d = {'functions', 'config', 'external'}
        safe_addpath(fullfile(root, d{1}));
    end
end

function fetch_submodules(root, missing)
    if ~isfolder(fullfile(root, '.git')) && ~isfile(fullfile(root, '.git'))
        fprintf(['%s is not a git clone, so the missing submodules (%s) cannot be ' ...
                 'fetched; see the report below.\n'], root, strjoin(missing, ', '));
        return
    end
    fprintf('Fetching required submodules: %s\n', strjoin(missing, ', '));
    [git, env] = git_program();
    [status, out] = run_command([{git, '-C', root, 'submodule', 'update', ...
        '--init', '--'}, missing], 'Echo', true, 'Env', env);
    if status == 127
        fprintf('git was not found; install git or download the submodules by hand.\n');
    elseif status ~= 0
        warning('prestus:submodules', 'git submodule update failed (exit %d):\n%s', ...
            status, strtrim(out));
    end
end

function [git, env] = git_program()
% git, and its environment, for a call from MATLAB. Linux: MATLAB's LD_LIBRARY_PATH makes
% git's https helper load MATLAB's own libssh2/libcurl/OpenSSL, which fail. Windows: the
% git.exe found on PATH may be the bare one (Git\ucrt64\bin or mingw64\bin), whose shell
% scripts ('git submodule') cannot find git-sh-setup; Git\cmd\git.exe sets up their
% environment.
    git = 'git';
    env = struct();
    if isunix && ~ismac
        env.LD_LIBRARY_PATH = '';
    elseif ispc
        [status, exec_path] = run_command({'git', '--exec-path'});
        if status == 0
            git_root = fileparts(fileparts(fileparts(strtrim(exec_path))));  % ...\Git
            wrapper = fullfile(git_root, 'cmd', 'git.exe');
            if isfile(wrapper)
                git = wrapper;
            end
        end
        fprintf('Using %s\n', git);
    end
end
