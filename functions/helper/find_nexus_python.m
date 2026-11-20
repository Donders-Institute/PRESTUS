function [python, tried] = find_nexus_python()
% FIND_NEXUS_PYTHON  A Python that can start the PRESTUS Nexus
%
% Tries, in order, and returns the first Python that can import PyYAML and
% start a PyQt6 application:
%   1. $PRESTUS_NEXUS_PYTHON
%   2. the Python remembered from the last successful start
%      (getpref('prestus', 'nexus_python'), set by PRESTUS_NEXUS)
%   3. the conda environment 'prestus-nexus' (created by nexus/install.sh,
%      nexus/install.bat or nexus/environment.yml)
%   4. python3 and python on the system path
%
% Use as:
%   [python, tried] = find_nexus_python()
%
% Output:
%   python - char Python executable, '' if none of the candidates works
%   tried  - struct array with fields python, ok, detail: every candidate
%            tested, for error messages
%
% See also: PRESTUS_NEXUS, CONDA_ENV_DIRS, PRESTUS_CHECK_INSTALL

    candidates = {};
    candidates{end + 1} = getenv('PRESTUS_NEXUS_PYTHON');
    if ispref('prestus', 'nexus_python')
        candidates{end + 1} = getpref('prestus', 'nexus_python');
    end
    envs = conda_env_dirs();
    [~, names] = cellfun(@fileparts, envs, 'UniformOutput', false);
    hits = envs(strcmp(names, 'prestus-nexus'));
    for k = 1:numel(hits)                  % by index: a 0x1 selection must not loop
        candidates{end + 1} = env_python(hits{k}); %#ok<AGROW>
    end
    candidates = [candidates, {'python3', 'python'}];
    candidates = candidates(~cellfun(@isempty, candidates));
    [~, keep] = unique(candidates, 'stable');
    candidates = candidates(sort(keep));

    python = '';
    tried = struct('python', {}, 'ok', {}, 'detail', {});
    for k = 1:numel(candidates)
        c = candidates{k};
        if any(c == '/' | c == '\') && ~isfile(c)
            continue                       % a remembered path that has since gone
        end
        [ok, detail] = nexus_python_ok(c);
        tried(end + 1) = struct('python', c, 'ok', ok, 'detail', detail); %#ok<AGROW>
        if ok
            python = c;
            return
        end
    end
end

function py = env_python(env)
    if ispc
        py = fullfile(env, 'python.exe');
    else
        py = fullfile(env, 'bin', 'python');
    end
end

function [ok, detail] = nexus_python_ok(python)
% The Nexus runs in the background with its output discarded, so check that
% this Python can start a Qt app (PyQt6 installed, platform plugin loadable).
    [status, cmdout] = run_command({python, '-c', ...
        'import yaml; from PyQt6.QtWidgets import QApplication; QApplication([])'});
    ok = status == 0;
    detail = strtrim(cmdout);
    if status == 127
        detail = 'not found';
    end
end
