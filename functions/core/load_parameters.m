function parameters = load_parameters(varargin)
% LOAD_PARAMETERS  Load default PRESTUS config and merge with one or more extra YAML files or structs
%
% Reads config/config_default.yaml relative to the PRESTUS root (resolved from
% this file's location so CWD is irrelevant), then deep-merges any caller-supplied
% YAML files or parameter structs on top. After merging, runs validation:
% interactive-mode check, transducer parameter derivation, thermal timing check,
% output-affix sanitization, path existence checks, and MATLAB version guard.
%
% Use as:
%   parameters = load_parameters()
%   parameters = load_parameters(extra_config_file)
%   parameters = load_parameters(extra_config_file, extra_config_location)
%   parameters = load_parameters(file1, loc1, file2, loc2, ...)
%
% Input:
%   varargin - (optional) variable-length argument list:
%              * 0 args:  load default config only
%              * 1 arg:   path to an extra YAML file, or a parameter struct to merge
%              * 2 args:  (filename, folder) pair for one extra YAML file
%              * 2N args: N (filename, folder) pairs, merged left-to-right
%
% Output:
%   parameters - merged and validated PRESTUS simulation parameters
%
% See also: LOAD_TRANSDUCER_PARAMETERS, PRINT_PARAMETER_SUMMARY, PATH_LOG_SETUP

    %% Load default configuration file
    % Use an absolute path derived from this file's location so that
    % load_parameters works correctly regardless of the caller's CWD.
    this_dir     = fileparts(mfilename('fullpath'));          % functions/core/
    prestus_root = fileparts(fileparts(this_dir));            % PRESTUS root
    toolbox_default = fullfile(prestus_root, 'config', 'config_default.yaml');

    % Ensure yaml (and other external dependencies) are on the path.
    % This makes load_parameters callable before a full PRESTUS path setup.
    external_dir = fullfile(prestus_root, 'external');
    if exist(external_dir, 'dir') && ~exist('yaml.loadFile', 'file')
        addpath(genpath(external_dir));
    end

    % Prefer config_default.yaml from the project config directory over the
    % toolbox copy. Check (in order):
    %   1. Explicit location argument (filename, location) pair
    %   2. Directory of a single full-path config file argument
    %   3. Current working directory
    %   4. Fall back to toolbox default
    default_config_path = toolbox_default;
    using_toolbox_default = true;

    if nargin >= 2 && ischar(varargin{2})
        % (filename, location) form — look for config_default in that location
        candidate = fullfile(varargin{2}, 'config_default.yaml');
        if exist(candidate, 'file')
            default_config_path = candidate;
            using_toolbox_default = false;
        end
    elseif nargin >= 1 && ischar(varargin{1})
        % Single full-path argument — check the file's own directory
        candidate = fullfile(fileparts(varargin{1}), 'config_default.yaml');
        if exist(candidate, 'file')
            default_config_path = candidate;
            using_toolbox_default = false;
        end
    end

    if using_toolbox_default
        % Check current working directory as a last resort before toolbox fallback
        candidate = fullfile(pwd, 'config_default.yaml');
        if exist(candidate, 'file')
            default_config_path = candidate;
            using_toolbox_default = false;
        end
    end

    if using_toolbox_default
        warning('prestus:toolboxDefault', ...
            ['Loading config_default.yaml from the PRESTUS toolbox directory.\n' ...
             'For reproducible, project-specific setups, run prestus_config_init() once:\n' ...
             '  prestus_config_init(''/path/to/project/config'', project_name=''projectX'');\n' ...
             'Then load parameters with:\n' ...
             '  parameters = load_parameters(''config_projectX.yaml'', ''/path/to/project/config'');']);
    end

    parameters = yaml.loadFile(default_config_path, "ConvertToArray", true);

    %% Merge with additional configuration files or structures

    if nargin == 1
        extra_config_file = varargin{1};
        if isstruct(extra_config_file)
            extra_parameters = extra_config_file; % Use provided struct directly
        else
            extra_parameters = yaml.loadFile(fullfile(extra_config_file), "ConvertToArray", true);
        end
        parameters = MergeStruct(parameters, extra_parameters);
    elseif nargin == 2
        extra_config_file = varargin{1};
        extra_config_location = varargin{2};
        extra_parameters = yaml.loadFile(fullfile(extra_config_location, extra_config_file), "ConvertToArray", true);
        parameters = MergeStruct(parameters, extra_parameters);
        parameters = restore_struct_arrays(parameters, extra_parameters);
    elseif nargin > 2
        % Check that extra inputs come in pairs (filename, location)
        if mod(numel(varargin), 2) ~= 0
            error('Extra parameters must be provided as name/value pairs.');
        end
        % Loop through each pair and assign to the structure
        for i = 1:2:numel(varargin)
            extra_config_file = varargin{i};
            extra_config_location = varargin{i+1};
            extra_parameters = yaml.loadFile(fullfile(extra_config_location, extra_config_file), "ConvertToArray", true);
            parameters = MergeStruct(parameters, extra_parameters);
            parameters = restore_struct_arrays(parameters, extra_parameters);
        end
    end

    %% Backwards compatibility: migrate calibration.target_isppa_wcm2 → transducer(i).target_isppa_wcm2
    %
    % The global calibration.target_isppa_wcm2 field is deprecated. Each transducer
    % now carries its own target_isppa_wcm2. If the old field is present and a
    % per-transducer value has not already been set, copy the global value to every
    % transducer entry and remove the old field.
    if isfield(parameters, 'calibration') && ...
            isfield(parameters.calibration, 'target_isppa_wcm2') && ...
            ~isempty(parameters.calibration.target_isppa_wcm2) && ...
            any(isfinite(parameters.calibration.target_isppa_wcm2))
        warning('prestus:deprecatedField', ...
            ['calibration.target_isppa_wcm2 is deprecated. ' ...
             'Set transducer.target_isppa_wcm2 instead (per-transducer). ' ...
             'The value has been migrated automatically this run.']);
        if isfield(parameters, 'transducer')
            for ti = 1:numel(parameters.transducer)
                if ~isfield(parameters.transducer(ti), 'target_isppa_wcm2') || ...
                        isempty(parameters.transducer(ti).target_isppa_wcm2)
                    parameters.transducer(ti).target_isppa_wcm2 = ...
                        parameters.calibration.target_isppa_wcm2;
                end
            end
        end
        parameters.calibration = rmfield(parameters.calibration, 'target_isppa_wcm2');
    end

    %% Check interactive mode requirements

    assert(parameters.simulation.interactive == 0 || usejava('desktop'), ...
           'MATLAB should run in desktop mode if parameters.simulation.interactive is enabled in PRESTUS config');

    %% Transducer settings validation and derived calculations
    parameters = load_transducer_parameters(parameters);

    %% Validate thermal simulation settings
    % Request no timing overview here

    if parameters.modules.run_heating_sims
        thermal_parameters(parameters, true);
    end

    %% Output file settings validation and sanitization

    % Sanitize output file affix to ensure valid characters only
    % yaml.loadFile reads empty-string YAML values as [] (numeric empty); coerce to char.
    if ~ischar(parameters.io.output_affix)
        parameters.io.output_affix = '';
    end
    % Allow hyphens (needed for BIDS entity values like desc-liberal).
    sanitized_affix = regexprep(parameters.io.output_affix, '[^a-zA-Z0-9_-]', '_');
    if ~strcmp(sanitized_affix, parameters.io.output_affix)
        fprintf('The original `io.output_affix` was sanitized. "%s" will be used instead of "%s"\n', ...
                sanitized_affix, parameters.io.output_affix);
    end
    % Wrap user-supplied label in the BIDS desc entity (_desc-<label>).
    % Strip any leading underscore before wrapping so users can write either
    % 'myRun' or '_myRun' in the config.
    if ~isempty(sanitized_affix) && ~startsWith(sanitized_affix, '_desc-')
        label = regexprep(sanitized_affix, '^_+', '');
        sanitized_affix = ['_desc-' label];
    end
    parameters.io.output_affix = sanitized_affix;

    %% Validate paths for required libraries and binaries

    % Check LD_LIBRARY_PATH existence and warn user if missing
    if isfield(parameters, 'hpc') && isfield(parameters.hpc, 'ld_library_path') && ...
       ~strcmp(parameters.hpc.ld_library_path, "") && ~exist(parameters.hpc.ld_library_path, 'dir')
        assert(all(confirmation_dlg('The path `hpc.ld_library_path` has been specified but does not exist. Do you want to continue?', ...
                                    'Yes', 'No')), 'Exiting');
    end

    % Check segmentation software path existence and warn user if missing
    if isfield(parameters.startup, 'simnibs_bin_path') && ~strcmp(parameters.startup.simnibs_bin_path, "") && ...
       ~exist(fullfile(parameters.startup.simnibs_bin_path, 'charm'), 'file')
        assert(all(confirmation_dlg(sprintf('charm does not exist at %s. Do you want to continue?', ...
                                            parameters.startup.simnibs_bin_path), ...
                                    'Yes', 'No')), 'Exiting');
    elseif (~isfield(parameters.startup, 'simnibs_bin_path') || strcmp(parameters.startup.simnibs_bin_path, "")) && contains(parameters.simulation.medium, {'layered'})
        warn('prestus:noSimNIBS', 'No path to SimNIBS binaries provided. Segmentation and MNI-conversion may fail...');
    end

    %% Default segmentation path fallback

    if ~isfield(parameters.path, 'seg') || isempty(parameters.path.seg) || strcmp(parameters.path.seg, '')
        parameters.path.seg = parameters.path.anat;
    end

    %% Convert additional paths into cell arrays for processing

    % Convert `paths_to_add` into cell array format (split by semicolon)
    if isfield(parameters.startup, 'paths_to_add') && ~isempty(parameters.startup.paths_to_add) ...
            && ~iscell(parameters.startup.paths_to_add)
        parameters.startup.paths_to_add = cellstr(strsplit(parameters.startup.paths_to_add, ';'));
    end

    % Convert `subpaths_to_add` into cell array format (split by semicolon)
    if isfield(parameters.startup, 'subpaths_to_add') && ~isempty(parameters.startup.subpaths_to_add) ...
            && ~iscell(parameters.startup.subpaths_to_add)
        parameters.startup.subpaths_to_add = cellstr(strsplit(parameters.startup.subpaths_to_add, ';'));
    end

    %% MATLAB version check for compatibility

    % Warn user about outdated MATLAB versions (< R2022b)
    if verLessThan('matlab', '9.13')
       assert(all(confirmation_dlg('MATLAB appears to be outdated. Please update before continuing. Do you want to continue?', ...
                                   'Yes', 'No')), 'Exiting');
    end

    %% [DEBUG] Summary of requested parameters
    % With debugging off, the parameters are saved to the log at the
    % start of pipeline execution (see path_log_setup).

    if parameters.simulation.debug == 1
        fprintf('PRELIMINARY specified parameters. Note that this may not be the final specification...\n');
        print_parameter_summary(parameters)
    end

end

function merged = restore_struct_arrays(merged, extra)
% After MergeStruct, top-level struct-array fields (e.g. transducer) may have
% been collapsed to a scalar because getfield on a struct array returns only
% the first element.  Restore any field that the extra config defined as a
% larger array than what survived the merge.
%
% Strategy: element(1) of merged is already the correct merge of defaults +
% extra element(1).  For elements 2..N we merge that same default-merged
% base with each extra element so every element retains all default fields.
    if ~isstruct(extra); return; end
    fn = fieldnames(extra);
    for i = 1:numel(fn)
        f = fn{i};
        if isfield(merged, f) && isstruct(extra.(f)) && ...
                numel(extra.(f)) > numel(merged.(f))
            base = merged.(f);   % 1×1: correctly merged defaults + extra element(1)
            result = base;
            for k = 2:numel(extra.(f))
                result(k) = MergeStruct(base, extra.(f)(k));
            end
            merged.(f) = result;
        end
    end
end
