function prestus_planner(config_file, options)
% PRESTUS_PLANNER  Open the PRESTUS Planner (standalone config builder)
%
% Starts the Python/Qt PRESTUS Planner from planner/ in the PRESTUS root
% without blocking MATLAB. The planner builds and checks a study config
% (override YAML for load_parameters), places transducers on the T1, and
% can validate or run the config in MATLAB. See planner/README.md.
%
% Use as:
%   prestus_planner()                          % untitled study (File > New study… / Open…)
%   prestus_planner('config_study.yaml')       % open a study config
%   prestus_planner(..., Python='/opt/miniconda3/envs/prestus-gui/bin/python')
%
% Input:
%   config_file - (optional) char path of a study config YAML
%
% Name-value options:
%   Python - (1,:) char Python executable with PyQt6, PyYAML (and numpy,
%            nibabel, matplotlib for the placement viewer); default 'python3'
%
% See also: PRESTUS_GUI, LOAD_PARAMETERS, PRESTUS_CONFIG_INIT

arguments
    config_file (1,:) char = ''
    options.Python (1,:) char = 'python3'
end

    planner_dir = fullfile(get_prestus_path(), 'planner');
    if ~isfolder(fullfile(planner_dir, 'prestus_planner'))
        error('prestus_planner:notFound', 'PRESTUS Planner not found in %s.', planner_dir);
    end

    % The planner runs in the background with its output discarded, so check
    % first that this Python can start a Qt app (PyQt6 installed, platform
    % plugin loadable) and report why not.
    [status, cmdout] = run_command({options.Python, '-c', ...
        'import yaml; from PyQt6.QtWidgets import QApplication; QApplication([])'});
    if status ~= 0
        error('prestus_planner:python', ...
            ['%s cannot start the PRESTUS Planner (exit %d):\n%s\n' ...
             'Pass a Python with PyQt6 and PyYAML, e.g. from the env in ' ...
             '%s:\n  conda env create -f environment.yml\n' ...
             '  prestus_planner(Python=''<conda>/envs/prestus-planner/bin/python'')'], ...
            options.Python, status, strtrim(cmdout), planner_dir);
    end

    argv = {options.Python, '-m', 'prestus_planner'};
    if ~isempty(config_file)
        argv{end + 1} = char(java.io.File(config_file).getAbsolutePath());
    end
    run_command(argv, 'Background', true, 'WorkingDir', planner_dir, ...
        'Env', struct('PYTHONPATH', planner_dir));
    fprintf('PRESTUS Planner started (%s).\n', options.Python);
end
