function report_path = prestus_group_report_start(config, subject_list, options)
% PRESTUS_GROUP_REPORT_START  Manual entry point for the group HTML report (and group MNI maps)
%
% Loads a PRESTUS configuration (YAML file path or pre-loaded struct) and
% produces a self-contained group HTML report. It never runs any
% simulations — it only reads existing per-run outputs (CSV, targeting
% sidecars, PNG, MNI NIfTIs).
%
% Study mode — the config sets study.manifest (see PRESTUS_STUDY_START):
%   runs come from the manifest (subjects x conditions, each with its own
%   settings); group MNI maps per condition are built with CREATE_GROUP_MAPS;
%   everything goes to <path.sim>/group/<study.name>/:
%     group_<name>_report.html, group_<name>_runs.csv, maps/, plots/
%
% Subject mode — no manifest: all subjects sharing simulation.medium and
%   io.output_affix (auto-discovered, or subject_list) form one group:
%     <path.sim>/group_<medium>_report<output_affix>.html
%   with group maps in <path.sim>/group/<medium><output_affix>/.
%
% Use as:
%   prestus_group_report_start(config)
%   prestus_group_report_start(config, subject_list)
%   prestus_group_report_start(config, [], 'Maps', false)
%   report_path = prestus_group_report_start(...)
%
% Input:
%   config       - YAML config file path (char/string) OR a pre-loaded
%                  parameters struct. Must define at minimum:
%                    path.sim
%                    simulation.medium
%                  Optional:
%                    io.output_affix, study.manifest, study.name,
%                    analysis.target_roi_mni (outlined on the maps)
%   subject_list - (optional) numeric array of subject IDs (subject mode
%                  only). If omitted or empty, the subject list is
%                  auto-discovered by scanning path.sim for sub-NNN folders
%                  that contain a matching per-subject CSV.
%   options      - name-value pairs:
%     'Maps'     build the group MNI maps (default true; layered media)
%
% Output:
%   report_path  - path to the generated HTML report, or '' if generation
%                  failed or no subjects were found.
%
% Examples:
%   % Study: runs from the manifest named in the study config
%   prestus_group_report_start('config/config_study.yaml');
%
%   % Explicit subject list, using a preloaded parameters struct
%   parameters = load_parameters('config/config_study.yaml');
%   prestus_group_report_start(parameters, [1 2 3 5 7]);
%
% See also: GENERATE_GROUP_REPORT, CREATE_GROUP_MAPS, STUDY_SETUP,
%           PRESTUS_STUDY_START, DISCOVER_GROUP_SUBJECTS, LOAD_PARAMETERS

arguments
    config
    subject_list (1,:) {mustBeNumeric} = []
    options.Maps (1,1) logical = true
end

    %% Study mode: runs from the manifest
    if isempty(subject_list) && has_manifest(config)
        study = study_setup(config);
        p = study.parameters;
        runs = study.runs(study.runs.include, :);
        fprintf('========================================\n');
        fprintf('PRESTUS STUDY REPORT: %s\n', study.name);
        fprintf('Runs:     %d (of %d manifest rows)\n', height(runs), height(study.runs));
        fprintf('Output:   %s\n', study.out_dir);
        fprintf('========================================\n\n');
        maps = group_maps(p, runs, study.out_dir, options.Maps);
        report_path = generate_group_report(p, runs, struct('name', study.name, ...
            'out_dir', study.out_dir, 'maps', maps, 'manifest', study.manifest_path));
        if ~isempty(report_path)
            fprintf('\nStudy report ready: %s\n', report_path);
        end
        return
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
    runs = table(subject_list(:), strings(numel(subject_list), 1), true(numel(subject_list), 1), ...
        repmat(string(parameters.io.output_affix), numel(subject_list), 1), ...
        repmat(string(parameters.simulation.medium), numel(subject_list), 1), ...
        repmat(string(parameters.path.sim), numel(subject_list), 1), ...
        compose("sub-%03d", subject_list(:)), ...
        'VariableNames', {'subject_id', 'condition', 'include', 'affix', 'medium', 'sim_dir', 'label'});
    maps = group_maps(parameters, runs, fullfile(parameters.path.sim, 'group', ...
        [char(parameters.simulation.medium) char(parameters.io.output_affix)]), options.Maps);
    if isempty(maps)
        report_path = generate_group_report(parameters, subject_list);
    else
        report_path = generate_group_report(parameters, subject_list, struct('maps', maps));
    end

    if ~isempty(report_path)
        fprintf('\nGroup report ready: %s\n', report_path);
    end
end

% -------------------------------------------------------------------------
function tf = has_manifest(config)
% Does the config (path or struct) set study.manifest? Reads the YAML only.
    tf = false;
    if ischar(config) || isstring(config)
        if ~isfile(config); return; end
        raw = yaml.loadFile(char(config));
    else
        raw = config;
    end
    tf = isstruct(raw) && isfield(raw, 'study') && isstruct(raw.study) && ...
        isfield(raw.study, 'manifest') && ~isempty(raw.study.manifest) && ...
        strlength(string(raw.study.manifest)) > 0;
end

function maps = group_maps(p, runs, out_dir, enabled)
% Group MNI maps (best-effort; never blocks the report).
    maps = struct([]);
    if ~enabled || ~contains(char(p.simulation.medium), 'layered'); return; end
    roi = '';
    if isfield(p, 'analysis') && isfield(p.analysis, 'target_roi_mni')
        roi = char(string(p.analysis.target_roi_mni));
    end
    seg = '';
    if isfield(p.path, 'seg'); seg = char(string(p.path.seg)); end
    try
        maps = create_group_maps(runs, out_dir, 'SegDir', seg, 'RoiMni', roi);
    catch ME
        warning('prestus_group_report_start:maps', 'Group maps failed: %s', ME.message);
    end
end
