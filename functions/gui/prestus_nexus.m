function prestus_nexus(config_file, options)
% PRESTUS_NEXUS  Open the PRESTUS Nexus (standalone config builder)
%
% Starts the Python/Qt PRESTUS Nexus from nexus/ in the PRESTUS root
% without blocking MATLAB. The Nexus builds and checks a study config
% (override YAML for load_parameters), places transducers on the T1, and
% can validate or run the config in MATLAB. See nexus/README.md.
%
% Use as:
%   prestus_nexus()                          % untitled study (File > New study… / Open…)
%   prestus_nexus('config_study.yaml')       % open a study config
%   prestus_nexus(..., Python='/opt/miniconda3/envs/prestus-nexus/bin/python')
%
% Input:
%   config_file - (optional) char path of a study config YAML
%
% Name-value options:
%   Python - (1,:) char Python executable of the Nexus environment; default ''
%            finds it (FIND_NEXUS_PYTHON: $PRESTUS_NEXUS_PYTHON, the Python
%            of the last start, the conda env 'prestus-nexus', python3). The
%            Python that started the Nexus is remembered for the next call.
%
% See also: FIND_NEXUS_PYTHON, LOAD_PARAMETERS, PRESTUS_CONFIG_INIT

arguments
    config_file (1,:) char = ''
    options.Python (1,:) char = ''
end

    nexus_dir = fullfile(get_prestus_path(), 'nexus');
    if ~isfolder(fullfile(nexus_dir, 'prestus_nexus'))
        error('prestus_nexus:notFound', 'PRESTUS Nexus not found in %s.', nexus_dir);
    end

    if isempty(options.Python)
        [python, tried] = find_nexus_python();
        if isempty(python)
            if ispc
                install = fullfile(nexus_dir, 'install.bat');
            else
                install = ['bash ' fullfile(nexus_dir, 'install.sh')];
            end
            error('prestus_nexus:python', ...
                ['No Python that can start the PRESTUS Nexus was found.\n%s\n' ...
                 'Install the Nexus environment once:\n  %s\n' ...
                 'or pass a Python: prestus_nexus(Python=''<env python>'')'], ...
                tried_summary(tried), install);
        end
    else
        python = options.Python;
        [status, cmdout] = run_command({python, '-c', ...
            'import yaml; from PyQt6.QtWidgets import QApplication; QApplication([])'});
        if status ~= 0
            error('prestus_nexus:python', ...
                '%s cannot start the PRESTUS Nexus (exit %d):\n%s', ...
                python, status, strtrim(cmdout));
        end
    end
    setpref('prestus', 'nexus_python', python);

    argv = {python, '-m', 'prestus_nexus'};
    if ~isempty(config_file)
        argv{end + 1} = char(java.io.File(config_file).getAbsolutePath());
    end
    run_command(argv, 'Background', true, 'WorkingDir', nexus_dir, ...
        'Env', struct('PYTHONPATH', nexus_dir));
    fprintf('PRESTUS Nexus started (%s).\n', python);
end

function txt = tried_summary(tried)
    lines = arrayfun(@(t) sprintf('  %s: %s', t.python, ...
        first_line(t.detail)), tried, 'UniformOutput', false);
    txt = ['Tried:' newline strjoin(lines, newline)];
end

function s = first_line(s)
    s = strtrim(regexprep(s, '[\r\n].*$', ''));
    if isempty(s)
        s = 'failed';
    end
end
