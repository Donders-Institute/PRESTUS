function prestus_config_init(target_dir, options)
% PRESTUS_CONFIG_INIT  Bootstrap a project-specific PRESTUS configuration folder
%
% Copies config_default.yaml from the toolbox into a project config directory
% and optionally creates an empty project-specific override file.
% Optionally opens the PRESTUS GUI for interactive refinement.
%
% Use as:
%   prestus_config_init(target_dir)
%   prestus_config_init(target_dir, project_name='myStudy')
%   prestus_config_init(target_dir, project_name='myStudy', open_gui=true)
%
% Input:
%   target_dir   - directory where the config files will be written;
%                  created if it does not exist
%   project_name - (name-value, default '') short label for the override file,
%                  e.g. 'myStudy' produces config_myStudy.yaml;
%                  if empty, only config_default.yaml is copied
%   open_gui     - (name-value, default false) open the PRESTUS GUI after
%                  writing the files so parameters can be reviewed visually
%
% Output files written to target_dir:
%   config_default.yaml        - full copy of the toolbox default (base config)
%   config_<project_name>.yaml - empty project override (add only what differs)
%
% Typical subsequent usage:
%   parameters = load_parameters('config_myStudy.yaml', '/path/to/project/config');

    arguments
        target_dir             (1,:) char
        options.project_name   (1,:) char   = ''
        options.open_gui       (1,1) logical = false
    end

    %% Resolve toolbox config_default.yaml
    this_dir     = fileparts(mfilename('fullpath'));
    prestus_root = fileparts(fileparts(this_dir));
    src_default  = fullfile(prestus_root, 'config', 'config_default.yaml');

    assert(exist(src_default, 'file') == 2, ...
           'Cannot find toolbox config_default.yaml at: %s', src_default);

    %% Create target directory if needed
    if ~exist(target_dir, 'dir')
        mkdir(target_dir);
        fprintf('Created config directory: %s\n', target_dir);
    end

    %% Copy config_default.yaml
    dst_default = fullfile(target_dir, 'config_default.yaml');
    if exist(dst_default, 'file')
        warning('prestus_config_init:exists', ...
                'config_default.yaml already exists in %s — skipping copy.', target_dir);
    else
        copyfile(src_default, dst_default);
        fprintf('Copied config_default.yaml → %s\n', dst_default);
    end

    %% Optionally create project-specific override file
    project_config_name = '';
    if ~isempty(options.project_name)
        project_config_name = sprintf('config_%s.yaml', options.project_name);
        dst_project = fullfile(target_dir, project_config_name);
        if exist(dst_project, 'file')
            warning('prestus_config_init:exists', ...
                    '%s already exists — skipping creation.', dst_project);
        else
            fid = fopen(dst_project, 'w');
            assert(fid ~= -1, 'Cannot write to %s', dst_project);
            fprintf(fid, '# Project-specific PRESTUS configuration: %s\n', options.project_name);
            fprintf(fid, '# Override only the parameters that differ from config_default.yaml.\n');
            fprintf(fid, '# Load with:\n');
            fprintf(fid, '#   parameters = load_parameters(''%s'', ''%s'');\n', ...
                    project_config_name, target_dir);
            fprintf(fid, '\n');
            fclose(fid);
            fprintf('Created project config: %s\n', dst_project);
        end
    end

    %% Summary
    fprintf('\nThe default configuration has been copied from PRESTUS to %s.\n', target_dir);
    if ~isempty(project_config_name)
        fprintf('To load parameters in your scripts:\n');
        fprintf('  parameters = load_parameters(''%s'', ''%s'');\n\n', ...
                project_config_name, target_dir);
    else
        fprintf(['Create a project-specific configuration config_<project>.yaml in %s, ' ...
                 'then load with:\n'], target_dir);
        fprintf('  parameters = load_parameters(''config_<project>.yaml'', ''%s'');\n\n', target_dir);
    end

    %% Optionally open GUI
    if options.open_gui
        fprintf('Opening PRESTUS GUI for interactive refinement...\n');
        if ~isempty(project_config_name)
            gui_config_file = project_config_name;
        else
            gui_config_file = 'config_default.yaml';
        end
        parameters = load_parameters(gui_config_file, target_dir);
        parameters.config_path = fullfile(target_dir, gui_config_file);
        prestus_gui(parameters);
    end

end
