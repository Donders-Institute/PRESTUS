function report_path = prestus_group_report_start(config, subject_list)
% PRESTUS_GROUP_REPORT_START  Manual entry point for the group HTML report
%
% Loads a PRESTUS configuration (YAML file path or pre-loaded struct),
% optionally accepts an explicit subject_list (otherwise auto-discovers via
% DISCOVER_GROUP_SUBJECTS), and calls GENERATE_GROUP_REPORT to produce a
% self-contained group HTML report at:
%
%     <path.sim>/group_<medium>_report<output_affix>.html
%
% This is the group-level analogue of PRESTUS_GROUP_START. It never runs any
% simulations — it only reads existing per-subject outputs (CSV + PNG) and
% writes one HTML file.
%
% Use as:
%   prestus_group_report_start(config)
%   prestus_group_report_start(config, subject_list)
%   report_path = prestus_group_report_start(...)
%
% Input:
%   config       - YAML config file path (char/string) OR a pre-loaded
%                  parameters struct. Must define at minimum:
%                    path.sim
%                    simulation.medium
%                  Optional:
%                    io.output_affix
%   subject_list - (optional) numeric array of subject IDs. If omitted or
%                  empty, the subject list is auto-discovered by scanning
%                  path.sim for sub-NNN folders that contain a matching
%                  per-subject CSV.
%
% Output:
%   report_path  - path to the generated HTML report, or '' if generation
%                  failed or no subjects were found.
%
% Examples:
%   % Auto-discover and report
%   prestus_group_report_start('config/config_study.yaml');
%
%   % Explicit subject list, using a preloaded parameters struct
%   parameters = load_parameters('config/config_study.yaml');
%   prestus_group_report_start(parameters, [1 2 3 5 7]);
%
% See also: GENERATE_GROUP_REPORT, DISCOVER_GROUP_SUBJECTS,
%           PRESTUS_GROUP_START, LOAD_PARAMETERS

arguments
    config
    subject_list (1,:) {mustBeNumeric} = []
end

    %% Load parameters
    if ischar(config) || isstring(config)
        fprintf('Loading config: %s\n', char(config));
        parameters = load_parameters(char(config));
    elseif isstruct(config)
        parameters = config;
    else
        error('prestus_group_report_start:badConfig', ...
            'config must be a YAML file path (char/string) or a parameters struct.');
    end

    %% Resolve subject list
    if isempty(subject_list)
        subject_list = discover_group_subjects(parameters);
    end

    %% Print summary
    fprintf('========================================\n');
    fprintf('PRESTUS GROUP HTML REPORT\n');
    fprintf('Subjects: %s (N=%d)\n', num2str(subject_list), numel(subject_list));
    fprintf('Medium:   %s\n', parameters.simulation.medium);
    if isfield(parameters, 'io') && isfield(parameters.io, 'output_affix') ...
            && ~isempty(parameters.io.output_affix)
        fprintf('Affix:    %s\n', parameters.io.output_affix);
    end
    fprintf('Output:   %s\n', parameters.path.sim);
    fprintf('========================================\n\n');

    if isempty(subject_list)
        fprintf('No subjects discovered — nothing to do.\n');
        report_path = '';
        return
    end

    %% Generate
    report_path = generate_group_report(parameters, subject_list);

    if ~isempty(report_path)
        fprintf('\nGroup report ready: %s\n', report_path);
    end
end
