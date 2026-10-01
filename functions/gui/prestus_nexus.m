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
%   prestus_nexus(..., Python='/opt/miniconda3/envs/prestus-gui/bin/python')
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

    nexus_dir = fullfile(get_prestus_path(), 'nexus');
    if ~isfolder(fullfile(nexus_dir, 'prestus_nexus'))
        error('prestus_nexus:notFound', 'PRESTUS Nexus not found in %s.', nexus_dir);
    end

    % The Nexus runs in the background with its output discarded, so check
    % first that this Python can start a Qt app (PyQt6 installed, platform
    % plugin loadable) and report why not.
    [status, cmdout] = run_command({options.Python, '-c', ...
        'import yaml; from PyQt6.QtWidgets import QApplication; QApplication([])'});
    if status ~= 0
        error('prestus_nexus:python', ...
            ['%s cannot start the PRESTUS Nexus (exit %d):\n%s\n' ...
             'Pass a Python with PyQt6 and PyYAML, e.g. from the env in ' ...
             '%s:\n  conda env create -f environment.yml\n' ...
             '  prestus_nexus(Python=''<conda>/envs/prestus-nexus/bin/python'')'], ...
            options.Python, status, strtrim(cmdout), nexus_dir);
    end

    argv = {options.Python, '-m', 'prestus_nexus'};
    if ~isempty(config_file)
        argv{end + 1} = char(java.io.File(config_file).getAbsolutePath());
    end
    run_command(argv, 'Background', true, 'WorkingDir', nexus_dir, ...
        'Env', struct('PYTHONPATH', nexus_dir));
    fprintf('PRESTUS Nexus started (%s).\n', options.Python);
end
