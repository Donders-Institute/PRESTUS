function subject_list = discover_group_subjects(parameters)
% DISCOVER_GROUP_SUBJECTS  Auto-scan path.sim for completed PRESTUS subjects
%
% Returns the unique sorted list of subject IDs that have a per-subject
% output table (CSV) matching the current simulation medium and output_affix.
% A subject is "found" iff its canonical per-subject CSV exists:
%
%     <path.sim>/sub-NNN/sub-NNN_<medium><output_affix>.csv
%
% This CSV is the same one produced by acoustic_analysis.m and (optionally)
% extended by thermal_analysis.m, so its presence is a reliable signal that
% the per-subject pipeline finished its analysis steps for that medium/affix.
%
% Deduplication: subject IDs are extracted from the sub-NNN folder name via
% regex, then passed through unique(), so each subject appears at most once
% even if e.g. sub-1 and sub-001 both exist (they'd collapse to the same int).
%
% Use as:
%   subject_list = discover_group_subjects(parameters)
%
% Input:
%   parameters - PRESTUS parameters struct. Must define:
%                  parameters.path.sim
%                  parameters.simulation.medium
%                Optional:
%                  parameters.io.output_affix  (defaults to '')
%
% Output:
%   subject_list - 1xN numeric row vector of unique sorted subject IDs.
%                  Empty if no matching subjects are found.
%
% Example:
%   parameters = load_parameters('config/config_study.yaml');
%   ids = discover_group_subjects(parameters);
%   fprintf('Found %d subjects: %s\n', numel(ids), num2str(ids));
%
% See also: GENERATE_GROUP_REPORT, PRESTUS_GROUP_REPORT_START

arguments
    parameters (1,1) struct
end

    subject_list = [];

    if ~isfield(parameters, 'path') || ~isfield(parameters.path, 'sim') ...
            || isempty(parameters.path.sim)
        warning('discover_group_subjects:noSimPath', ...
            'parameters.path.sim is not set; cannot scan for subjects.');
        return
    end
    if ~isfolder(parameters.path.sim)
        warning('discover_group_subjects:simPathMissing', ...
            'parameters.path.sim does not exist: %s', parameters.path.sim);
        return
    end

    medium = parameters.simulation.medium;
    if isfield(parameters, 'io') && isfield(parameters.io, 'output_affix')
        affix = parameters.io.output_affix;
    else
        affix = '';
    end

    % Enumerate sub-* directories
    entries = dir(fullfile(parameters.path.sim, 'sub-*'));
    if isempty(entries)
        return
    end

    ids = [];
    for k = 1:numel(entries)
        e = entries(k);
        if ~e.isdir, continue; end

        tok = regexp(e.name, '^sub-(\d+)$', 'tokens', 'once');
        if isempty(tok), continue; end
        id = str2double(tok{1});
        if ~isfinite(id), continue; end

        csv_path = fullfile(parameters.path.sim, sprintf('sub-%03d', id), ...
                            sprintf('sub-%03d_%s%s.csv', id, medium, affix));
        if isfile(csv_path)
            ids(end+1) = id; %#ok<AGROW>
        end
    end

    subject_list = unique(ids);  % unique() returns sorted ascending
    subject_list = subject_list(:).';  % force row vector
end
