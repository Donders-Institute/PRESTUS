function [sessions, labels] = neuronav_scan_sessions(parameters, subjects, session_ids)
% NEURONAV_SCAN_SESSIONS  Inventory raw Localite data and planned-target labels.
%
% Lists every sub-*/ses-* folder under path.localite_raw with the marker
% files found and the named planned targets (EntryTarget exports and
% GUMMarkers EntryTargetPairs) of each session. The label table is the
% basis of the label map: when neuronav.target_labels is empty, the labels
% are used as found; its .stem column proposes a grouping of numbered
% variants (e.g. 'L_caudateV2_1', 'L_caudateV2_2' -> 'L_caudateV2') that can be
% turned into a map with NEURONAV_LABEL_MAP_FROM_SCAN.
%
% The session table is also the selection input of NEURONAV_INGEST_BATCH.
%
% Use as:
%   [sessions, labels] = neuronav_scan_sessions(parameters)
%   [sessions, labels] = neuronav_scan_sessions(parameters, subjects)
%   [sessions, labels] = neuronav_scan_sessions(parameters, subjects, session_ids)
%
% Input:
%   subjects    - cell of 'sub-XXX' ([] / omitted = all sub-* folders)
%   session_ids - cell of 'ses-XX' ([] / omitted = all ses-* folders)
%
% Output:
%   sessions - table, one row per subject/session:
%                sub, ses, n_trigger (usable TriggerMarkers files, > 10 kB),
%                n_gum, n_entrytarget, labels (planned-target labels,
%                comma-separated), ingested (output XML exists)
%   labels   - table, one row per distinct label:
%                label, stem, target (after neuronav.target_labels),
%                n_subjects, n_sessions, sessions ('sub/ses, …')
%
% See also: NEURONAV_INGEST_BATCH, NEURONAV_MAP_LABEL, NEURONAV_LABEL_MAP_FROM_SCAN

    if nargin < 2, subjects = []; end
    if nargin < 3, session_ids = []; end

    paths = neuronav_resolve_paths(parameters);
    if isempty(paths.raw) || ~isfolder(paths.raw)
        error('neuronav_scan_sessions:paths', 'path.localite_raw is not set or not a folder.');
    end
    loc = 'loc';
    if isfield(parameters, 'placement') && isfield(parameters.placement, 'localite') && ...
            isfield(parameters.placement.localite, 'loc_folder')
        loc = char(parameters.placement.localite.loc_folder);
    end
    target_labels = [];
    if isfield(parameters, 'neuronav') && isfield(parameters.neuronav, 'target_labels')
        target_labels = parameters.neuronav.target_labels;
    end

    if isempty(subjects), subjects = list_dirs(paths.raw, 'sub-*'); end
    subjects = cellstr(subjects);

    rows = {};
    all_labels = {};  where = {};
    for i = 1:numel(subjects)
        sub = subjects{i};
        ses_list = list_dirs(fullfile(paths.raw, sub), 'ses-*');
        if ~isempty(session_ids)
            ses_list = intersect(ses_list, cellstr(session_ids), 'stable');
        end
        for j = 1:numel(ses_list)
            ses = ses_list{j};
            base = fullfile(paths.raw, sub, ses, loc);
            trig = dir(fullfile(base, '**', 'TMSTrigger', 'TriggerMarkers_Coil*.xml'));
            trig = trig([trig.bytes] > 10000);
            gum  = dir(fullfile(base, '**', 'GUMMarkers*.xml'));
            et   = dir(fullfile(base, '**', 'EntryTarget*.xml'));
            lab = {};
            if isfolder(base)
                evalc('c = neuronav_read_label_candidates(paths.raw, sub, ses, loc);');
                lab = {c.description};
            end
            ingested = ~isempty(paths.post) && isfile(fullfile(paths.post, sub, ses, ...
                sprintf('InstrumentMarker_%s_%s.xml', sub, ses)));
            rows(end+1, :) = {sub, ses, numel(trig), numel(gum), numel(et), ...
                strjoin(lab, ', '), ingested}; %#ok<AGROW>
            all_labels = [all_labels, lab]; %#ok<AGROW>
            where = [where, repmat({[sub '/' ses]}, 1, numel(lab))]; %#ok<AGROW>
        end
    end

    vars = {'sub', 'ses', 'n_trigger', 'n_gum', 'n_entrytarget', 'labels', 'ingested'};
    if isempty(rows)
        sessions = cell2table(cell(0, numel(vars)), 'VariableNames', vars);
    else
        sessions = cell2table(rows, 'VariableNames', vars);
    end

    [u, ~, idx] = unique(all_labels, 'stable');
    lrows = cell(numel(u), 6);
    for k = 1:numel(u)
        w = unique(where(idx == k), 'stable');
        subs = unique(cellfun(@(x) strtok(x, '/'), w, 'UniformOutput', false));
        lrows(k, :) = {u{k}, label_stem(u{k}), neuronav_map_label(u{k}, target_labels), ...
            numel(subs), numel(w), strjoin(w, ', ')};
    end
    labels = cell2table(lrows, 'VariableNames', ...
        {'label', 'stem', 'target', 'n_subjects', 'n_sessions', 'sessions'});
end


function s = label_stem(label)
% Drop a trailing separator + number ('L_MD_2' -> 'L_MD', 'target 3' -> 'target').
    s = regexprep(label, '[_\-\s]+\d+$', '');
    if isempty(s), s = label; end
end

function names = list_dirs(parent, pattern)
    names = {};
    if ~isfolder(parent), return; end
    d = dir(fullfile(parent, pattern));
    d = d([d.isdir]);
    names = sort({d.name});
end
