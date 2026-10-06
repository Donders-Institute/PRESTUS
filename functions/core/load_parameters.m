function parameters = load_parameters(varargin)
% LOAD_PARAMETERS  Load default PRESTUS config and merge with one or more extra YAML files or structs
%
% Reads config/config_default.yaml relative to the PRESTUS root (resolved from
% this file's location so CWD is irrelevant). When a project-specific
% config_default.yaml is found (see below), it is merged on top of the toolbox
% default, so keys added to the toolbox since the project copy was made still
% get their defaults (a warning lists them). Caller-supplied YAML files or
% parameter structs are then deep-merged on top; keys in those YAML files that
% PRESTUS does not know are reported (CONFIG_UNKNOWN_KEYS). After merging, runs
% validation:
% interactive-mode check, transducer parameter derivation, thermal timing check,
% output-affix sanitization, path existence checks, and MATLAB version guard.
%
% Use as:
%   parameters = load_parameters()
%   parameters = load_parameters(extra_config_file)
%   parameters = load_parameters(extra_config_file, extra_config_location)
%   parameters = load_parameters(file1, loc1, file2, loc2, ...)
%   parameters = load_parameters(..., {'subject_id', 3; 'io.output_affix', 'M1'})
%
% Input:
%   varargin - (optional) variable-length argument list:
%              * 0 args:  load default config only
%              * 1 arg:   path to an extra YAML file, or a parameter struct to merge
%              * 2 args:  (filename, folder) pair for one extra YAML file
%              * 2N args: N (filename, folder) pairs, merged left-to-right
%              * a trailing N×2 cell {'dotted.key', value; ...} of overrides,
%                applied after all merges and before validation/derivation
%                (keys may index struct arrays, e.g. 'transducer(2).serial';
%                see APPLY_CONFIG_OVERRIDES)
%
% Output:
%   parameters - merged and validated PRESTUS simulation parameters
%
% See also: LOAD_TRANSDUCER_PARAMETERS, PRINT_PARAMETER_SUMMARY, PATH_LOG_SETUP

    % Trailing cell = key/value overrides (applied after all merges)
    overrides = {};
    if ~isempty(varargin) && iscell(varargin{end})
        overrides = varargin{end};
        varargin(end) = [];
    end
    nargs = numel(varargin);

    %% Load default configuration file
    % Use an absolute path derived from this file's location so that
    % load_parameters works correctly regardless of the caller's CWD.
    this_dir     = fileparts(mfilename('fullpath'));          % functions/core/
    prestus_root = fileparts(fileparts(this_dir));            % PRESTUS root
    toolbox_default = fullfile(prestus_root, 'config', 'config_default.yaml');

    % Ensure yaml (and other external dependencies) are on the path.
    % This makes load_parameters callable before a full PRESTUS path setup.
    external_dir = fullfile(prestus_root, 'external');
    % (exist() is 0 for package functions such as yaml.loadFile, so use which)
    if exist(external_dir, 'dir') && isempty(which('yaml.loadFile'))
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

    if nargs >= 2 && ischar(varargin{2})
        % (filename, location) form — look for config_default in that location
        candidate = fullfile(varargin{2}, 'config_default.yaml');
        if exist(candidate, 'file')
            default_config_path = candidate;
            using_toolbox_default = false;
        end
    elseif nargs >= 1 && ischar(varargin{1})
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
            if ~same_file(candidate, toolbox_default)
                warning('prestus:cwdDefault', ...
                    ['Using config_default.yaml from the current folder (%s). Pass the ' ...
                     'config folder explicitly, load_parameters(file, folder), to make ' ...
                     'the project defaults independent of where MATLAB runs.'], candidate);
            end
        end
    end

    if ~using_toolbox_default && same_file(default_config_path, toolbox_default)
        using_toolbox_default = true;
    end

    if using_toolbox_default
        warning('prestus:toolboxDefault', ...
            ['Loading config_default.yaml from the PRESTUS toolbox directory.\n' ...
             'For reproducible, project-specific setups, run prestus_config_init() once:\n' ...
             '  prestus_config_init(''/path/to/project/config'', project_name=''projectX'');\n' ...
             'Then load parameters with:\n' ...
             '  parameters = load_parameters(''config_projectX.yaml'', ''/path/to/project/config'');']);
    end

    % Toolbox defaults are always the base layer; a project copy of
    % config_default.yaml (prestus_config_init) pins its values on top.
    parameters = yaml_load_validated(toolbox_default);
    if ~using_toolbox_default
        fprintf('Project defaults: %s\n', default_config_path);
        project_defaults = yaml_load_validated(default_config_path);
        report_stale_project_default(default_config_path, parameters, project_defaults);
        parameters = mergestruct(parameters, project_defaults);
        parameters = restore_struct_arrays(parameters, project_defaults);
    end
    known_keys  = parameters;          % for the unknown-key check below
    extra_keys  = config_known_extra_keys(prestus_root);
    calibration_cfg = fullfile(prestus_root, 'config', 'config_calibration.yaml');
    if isfile(calibration_cfg)
        % Calibration-workflow settings are loaded into parameters.calibration
        if ~isfield(known_keys, 'calibration') || ~isstruct(known_keys.calibration)
            known_keys.calibration = struct();
        end
        known_keys.calibration = mergestruct(known_keys.calibration, ...
            yaml_load_validated(calibration_cfg));
    end

    %% Merge with additional configuration files or structures

    if nargs == 1
        extra_config_file = varargin{1};
        if isstruct(extra_config_file)
            extra_parameters = extra_config_file; % Use provided struct directly
        else
            extra_parameters = yaml_load_validated(fullfile(extra_config_file));
            report_unknown_keys(extra_parameters, known_keys, extra_keys, extra_config_file);
        end
        parameters = mergestruct(parameters, extra_parameters);
        parameters = restore_struct_arrays(parameters, extra_parameters);
    elseif nargs == 2
        extra_config_file = varargin{1};
        extra_config_location = varargin{2};
        extra_parameters = yaml_load_validated(fullfile(extra_config_location, extra_config_file));
        report_unknown_keys(extra_parameters, known_keys, extra_keys, extra_config_file);
        parameters = mergestruct(parameters, extra_parameters);
        parameters = restore_struct_arrays(parameters, extra_parameters);
    elseif nargs > 2
        % Check that extra inputs come in pairs (filename, location)
        if mod(numel(varargin), 2) ~= 0
            error('Extra parameters must be provided as name/value pairs.');
        end
        % Loop through each pair and assign to the structure
        for i = 1:2:numel(varargin)
            extra_config_file = varargin{i};
            extra_config_location = varargin{i+1};
            extra_parameters = yaml_load_validated(fullfile(extra_config_location, extra_config_file));
            report_unknown_keys(extra_parameters, known_keys, extra_keys, extra_config_file);
            parameters = mergestruct(parameters, extra_parameters);
            parameters = restore_struct_arrays(parameters, extra_parameters);
        end
    end

    %% Key/value overrides (e.g. one row of a study manifest)
    if ~isempty(overrides)
        parameters = apply_config_overrides(parameters, overrides);
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
    % Users can write either 'myRun' or '_myRun' in the config.
    parameters.io.output_affix = sanitize_output_affix(sanitized_affix);

    %% Validate paths for required libraries and binaries

    % Check LD_LIBRARY_PATH existence and warn user if missing
    if isfield(parameters, 'hpc') && isfield(parameters.hpc, 'ld_library_path') && ...
       ~strcmp(parameters.hpc.ld_library_path, "") && ~exist(parameters.hpc.ld_library_path, 'dir')
        confirm_or_warn(parameters, 'prestus:missingLdLibraryPath', ...
            sprintf('The path `hpc.ld_library_path` (%s) has been specified but does not exist.', ...
                    parameters.hpc.ld_library_path));
    end

    % Check segmentation software path existence and warn user if missing
    if isfield(parameters.startup, 'simnibs_bin_path') && ~strcmp(parameters.startup.simnibs_bin_path, "") && ...
       ~exist(fullfile(parameters.startup.simnibs_bin_path, 'charm'), 'file')
        confirm_or_warn(parameters, 'prestus:missingCharm', ...
            sprintf('charm does not exist at %s (startup.simnibs_bin_path).', ...
                    parameters.startup.simnibs_bin_path));
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
    if isMATLABReleaseOlderThan('R2022b')
        confirm_or_warn(parameters, 'prestus:oldMatlab', ...
            sprintf('MATLAB %s is older than R2022b; please update.', version('-release')));
    end

    %% Validate acoustic_backend value
    if isfield(parameters, 'acoustic_backend') && ~isempty(parameters.acoustic_backend)
        valid_backends = {'kwave', 'babelbrain', 'jwave'};
        if ~any(strcmpi(parameters.acoustic_backend, valid_backends))
            error('prestus:unknownAcousticBackend', ...
                'Unknown acoustic_backend "%s". Valid options: %s.', ...
                parameters.acoustic_backend, strjoin(valid_backends, ', '));
        end
    end

    %% Validate thermal_backend value and backend combination
    % Fail here rather than after the (expensive) acoustic stage.
    if isfield(parameters, 'thermal_backend') && ~isempty(parameters.thermal_backend)
        valid_thermal = {'kwave', 'babelbrain'};
        if ~any(strcmpi(parameters.thermal_backend, valid_thermal))
            error('prestus:unknownThermalBackend', ...
                'Unknown thermal_backend "%s". Valid options: %s.', ...
                parameters.thermal_backend, strjoin(valid_thermal, ', '));
        end
        if strcmpi(parameters.thermal_backend, 'babelbrain')
            ac_is_bb = isfield(parameters, 'acoustic_backend') && ...
                strcmpi(parameters.acoustic_backend, 'babelbrain');
            has_ext_ac = isfield(parameters, 'io') && ...
                isfield(parameters.io, 'external_acoustic_nifti') && ...
                ~isempty(parameters.io.external_acoustic_nifti);
            if ~ac_is_bb || has_ext_ac
                error('prestus:thermalBackendRequiresBabelBrainAcoustics', ...
                    ['thermal_backend "babelbrain" requires acoustic_backend "babelbrain" ' ...
                     '(it runs on the *DataForSim.h5 BabelBrain writes) and cannot be combined ' ...
                     'with io.external_acoustic_nifti. Use thermal_backend "kwave" with ' ...
                     'k-Wave/jwave/external acoustics.']);
            end
        end
    end

    %% Normalise babelbrain.computing_backend names to BabelBrain's integer codes
    if isfield(parameters, 'babelbrain') && isfield(parameters.babelbrain, 'computing_backend') && ...
            (ischar(parameters.babelbrain.computing_backend) || isstring(parameters.babelbrain.computing_backend))
        cb_names = {'cuda', 'opencl', 'metal', 'cpu'};
        cb = lower(char(parameters.babelbrain.computing_backend));
        cb_idx = find(strcmp(cb, cb_names));
        if ~isempty(cb_idx)
            parameters.babelbrain.computing_backend = cb_idx;
        elseif ismember(str2double(cb), 1:4)
            parameters.babelbrain.computing_backend = str2double(cb);
        elseif ~strcmp(cb, 'auto')
            error('prestus:unknownBabelBrainComputingBackend', ...
                ['Unknown babelbrain.computing_backend "%s". Valid options: auto, cuda, opencl, ' ...
                 'metal, cpu (or 1-4).'], cb);
        end
    end

    %% jwave compatibility: kWaveArray source is not supported
    if isfield(parameters, 'acoustic_backend') && ...
            strcmpi(parameters.acoustic_backend, 'jwave') && ...
            isfield(parameters, 'grid') && isfield(parameters.grid, 'use_kWaveArray') && ...
            parameters.grid.use_kWaveArray
        warning('prestus:jwaveKWaveArrayUnsupported', ...
            ['jwave backend does not support kWaveArray sources. ' ...
             'Overriding grid.use_kWaveArray to 0.']);
        parameters.grid.use_kWaveArray = 0;
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
                result(k) = mergestruct(base, extra.(f)(k));
            end
            merged.(f) = result;
        end
    end

    % Verify that all elements of any restored struct array carry the same field
    % set. Heterogeneous arrays cause cryptic errors deep in the pipeline.
    fn_merged = fieldnames(merged);
    for i = 1:numel(fn_merged)
        f = fn_merged{i};
        if isstruct(merged.(f)) && numel(merged.(f)) > 1
            base_fields = sort(fieldnames(merged.(f)(1)));
            for k = 2:numel(merged.(f))
                elem_fields = sort(fieldnames(merged.(f)(k)));
                assert(isequal(base_fields, elem_fields), 'prestus:heterogeneousStructArray', ...
                    ['Field "%s" element %d has different ' ...
                     'fields from element 1 after config merge. ' ...
                     'Check that your override config defines the same sub-fields for every ' ...
                     'array element (or leave extra elements entirely to inherit defaults).'], f, k);
            end
        end
    end
end

function s = yaml_load_validated(yaml_path)
% Wrapper around yaml.loadFile that fails early with a clear message when
% the returned value is not a struct (malformed or empty YAML file).
    s = yaml.loadFile(yaml_path, "ConvertToArray", true);
    if ~isstruct(s)
        error('load_parameters:invalidConfig', ...
            'Config file did not parse to a struct (got %s):\n  %s', class(s), yaml_path);
    end
end

function confirm_or_warn(parameters, id, message)
% Ask whether to continue in an interactive desktop session; otherwise (batch,
% HPC, simulation.interactive = 0) warn and continue, as no dialog can be shown.
    interactive = isfield(parameters, 'simulation') && ...
        isfield(parameters.simulation, 'interactive') && ...
        parameters.simulation.interactive && usejava('desktop');
    if interactive
        assert(all(confirmation_dlg([message ' Do you want to continue?'], 'Yes', 'No')), ...
            id, 'Exiting: %s', message);
    else
        warning(id, '%s Continuing (non-interactive session).', message);
    end
end

function report_stale_project_default(project_path, toolbox, project)
% Warn when a project copy of config_default.yaml lacks keys that the toolbox
% default now has; those keys take their toolbox default values.
    missing = config_unknown_keys(toolbox, project);
    if isempty(missing)
        return
    end
    shown = strjoin(missing(1:min(end, 12)), ', ');
    if numel(missing) > 12
        shown = sprintf('%s, ... (%d more)', shown, numel(missing) - 12);
    end
    warning('prestus:staleProjectDefault', ...
        ['Project config_default.yaml (%s) lacks %d key(s) added to the toolbox default ' ...
         'since it was copied; their toolbox defaults are used: %s. Run ' ...
         'prestus_config_upgrade(''%s'') to refresh the copy.'], ...
        project_path, numel(missing), shown, fileparts(project_path));
end

function report_unknown_keys(user, known, extra_keys, source)
% Warn about study-config keys PRESTUS does not read (typos are otherwise
% silently ignored and the default value is used).
    [unknown, hints] = config_unknown_keys(user, known, extra_keys);
    if isempty(unknown)
        return
    end
    items = cell(size(unknown));
    for i = 1:numel(unknown)
        if isempty(hints{i})
            items{i} = unknown{i};
        else
            items{i} = sprintf('%s (did you mean ''%s''?)', unknown{i}, hints{i});
        end
    end
    warning('prestus:unknownConfigKey', ...
        'Config %s has key(s) PRESTUS does not use; they are ignored:\n  %s', ...
        char(source), strjoin(items, sprintf('\n  ')));
end

function tf = same_file(a, b)
% True when two paths name the same file (resolves relative paths and links).
    tf = false;
    if ~isfile(a) || ~isfile(b)
        return
    end
    ia = dir(a);
    ib = dir(b);
    tf = strcmp(fullfile(ia.folder, ia.name), fullfile(ib.folder, ib.name));
end
