function target_map = neuronav_build_target_map(coil_map)
% NEURONAV_BUILD_TARGET_MAP  Expand a coil_map spec into a flat target_map.
%
% Converts the coil_map configuration (one entry per physical transducer) into
% the flat target_map struct array consumed by neuronav_ingest_markers and
% related functions.
%
% Each coil_map entry must have:
%   .coil     - hardware coil index (0-based, matches Localite trigger files)
%
% and either:
%   .target   - target label string for a single sequence (e.g. 'L_M1')
%
% or:
%   .sequences - struct array, one element per sequence on this coil, each with:
%                  .target  - target label for that sequence
%
% Use sequences when a coil records multiple position blocks in a session:
%   - Same target repeated (e.g. two stim blocks) → names get a _N suffix:
%       'L_M1', 'L_M1' → 'L_M1_1', 'L_M1_2'
%   - Different targets (coil repositioned) → names stay distinct:
%       'R_M1', 'R_V1' → 'R_M1', 'R_V1'
%   - Mixed: a target that appears once keeps its bare name; one that
%     appears multiple times gets _N on all of its occurrences.
%
% Optional fields forwarded to every target_map entry they belong to:
%   coil_map(i).localite_key          - passed to all sequences of this coil
%   coil_map(i).sequences(j).localite_key - overrides coil-level key for seq j
%   Each localite_key may be a single char/string, or a cell array of
%   candidate keys when more than one EntryTarget description convention is
%   in use (e.g. {'LEFT V', 'L_V'}) -- neuronav_match_entry_targets tries
%   each candidate in order and uses the first that matches.
%
% Output:
%   target_map - struct array with one entry per sequence, fields:
%                  .series_index  - 1-based position in the flat list
%                  .transducer_id - hardware coil index (.coil value)
%                  .target_name   - full label, e.g. 'L_M1' or 'L_M1_2'
%                  .color         - default plot colour (#rrggbb)
%                  .localite_key  - set from coil_map if present, else ''
%
% See also: NEURONAV_INGEST_MARKERS, NEURONAV_MATCH_ENTRY_TARGETS

colors = {'#ff0000','#00ff00','#0000ff','#ff00ff','#ffff00','#00ffff'};

% Every target_map entry is produced by copying this template and filling
% in its fields, so the field set (names, order, count) is identical across
% every entry regardless of which optional inputs were supplied.
template = struct( ...
    'series_index',  [], ...
    'transducer_id', [], ...
    'target_name',   '', ...
    'color',         '', ...
    'localite_key',  '');
target_map = repmat(template, 0, 0);

if isempty(coil_map)
    error('neuronav_build_target_map: coil_map is empty');
end
if ~isstruct(coil_map)
    error('neuronav_build_target_map: coil_map must be a struct array, got %s', class(coil_map));
end

for ci = 1:numel(coil_map)
    cm = coil_map(ci);

    if ~isfield(cm, 'coil') || isempty(cm.coil)
        error('neuronav_build_target_map: coil entry %d is missing required field .coil', ci);
    end

    % Normalise to a sequences list
    if isfield(cm, 'sequences') && ~isempty(cm.sequences)
        seqs = cm.sequences;
    elseif isfield(cm, 'target') && ~isempty(cm.target)
        seqs = struct('target', {cm.target});
    else
        error('neuronav_build_target_map: coil entry %d has neither .target nor .sequences', ci);
    end

    if ~isfield(seqs, 'target') || any(cellfun(@isempty, {seqs.target}))
        error('neuronav_build_target_map: coil entry %d has a sequence with an empty/missing .target', ci);
    end

    % Normalise every target name to plain char, so downstream comparisons
    % (strcmp, containers.Map keys) behave the same regardless of whether
    % the config supplied char or string.
    all_targets = cellfun(@char, {seqs.target}, 'UniformOutput', false);

    % Count how many times each target name appears across this coil's sequences.
    % Those appearing more than once will have _N appended to all occurrences.
    unique_targets = unique(all_targets);
    needs_suffix = containers.Map('KeyType','char','ValueType','logical');
    for ui = 1:numel(unique_targets)
        needs_suffix(unique_targets{ui}) = sum(strcmp(all_targets, unique_targets{ui})) > 1;
    end

    occurrence = containers.Map('KeyType','char','ValueType','double');

    for si = 1:numel(seqs)
        tname = all_targets{si};

        if needs_suffix(tname)
            if isKey(occurrence, tname)
                occurrence(tname) = occurrence(tname) + 1;
            else
                occurrence(tname) = 1;
            end
            full_name = sprintf('%s_%d', tname, occurrence(tname));
        else
            full_name = tname;
        end

        entry = template;
        entry.series_index  = numel(target_map) + 1;
        entry.transducer_id = cm.coil;
        entry.target_name   = full_name;
        entry.color         = colors{mod(numel(target_map), numel(colors)) + 1};

        % localite_key: sequence-level overrides coil-level. Falls back to
        % the template's '' default when neither is set. May be a single
        % char/string or a cell array of candidate keys (tried in order by
        % neuronav_match_entry_targets) -- pass cell arrays through as-is.
        if isfield(seqs(si), 'localite_key') && ~isempty(seqs(si).localite_key)
            entry.localite_key = normalize_localite_key(seqs(si).localite_key);
        elseif isfield(cm, 'localite_key') && ~isempty(cm.localite_key)
            entry.localite_key = normalize_localite_key(cm.localite_key);
        end

        target_map(end+1) = entry; %#ok<AGROW>
    end
end
end

% =========================================================================
%  LOCAL HELPER: normalize a localite_key value
% =========================================================================
function key = normalize_localite_key(key)
% Pass a cell array of candidate keys through unchanged; coerce anything
% else (char/string) to plain char, matching the rest of this file's
% target_name handling.
    if iscell(key)
        return;
    end
    key = char(key);
end
