function [localite_mm, localite_vox] = neuronav_match_entry_targets( ...
        entry_targets, target_map, t1_header)
% NEURONAV_MATCH_ENTRY_TARGETS  Match Localite EntryTarget entries to a target map.
%
% For each entry in target_map, finds the corresponding operator-confirmed brain
% target from the Localite EntryTarget XML.
%
% Matching strategy (in order of priority):
%   1. If target_map(i).localite_key is set, use it as a case-insensitive
%      substring search against EntryTarget descriptions.
%   2. Otherwise, use target_map(i).target_name as the search key. A trailing
%      '_<N>' suffix is stripped before matching because a coil may produce
%      multiple position sequences within one session (all aimed at the same
%      physical target). Each sequence gets its own target_map entry named
%      '<target>_1', '<target>_2', etc., but all of them correspond to the
%      same single EntryTarget entry.
%
% Use as:
%   [mm, vox] = neuronav_match_entry_targets(entry_targets, target_map)
%   [mm, vox] = neuronav_match_entry_targets(entry_targets, target_map, t1_header)
%
% Input:
%   entry_targets  - struct array returned by neuronav_read_entry_target,
%                    with fields .index, .description, .xyz [1x3 mm]
%   target_map     - struct array, one entry per coil sequence.
%                    Required field:
%                      .target_name  - label string used as the fallback search
%                                      key (case-insensitive substring match
%                                      against EntryTarget descriptions).
%                                      When a coil has multiple sequences for
%                                      the same target, name them '<target>_1',
%                                      '<target>_2', etc.; the trailing '_N'
%                                      is stripped automatically before matching.
%                    Optional field:
%                      .localite_key - explicit case-insensitive substring to
%                                      search for in descriptions. Takes
%                                      precedence over target_name matching.
%                                      Required when description conventions
%                                      differ from target_name (e.g. target
%                                      'L_M1' described as 'LEFT M1 Target').
%                                      May be a single char/string, or a cell
%                                      array of candidate keys when more than
%                                      one description convention is in use
%                                      (e.g. {'LEFT M', 'L_M'}); candidates
%                                      are tried in order, first match wins.
%   t1_header      - (optional) niftiinfo struct for the planning T1 image.
%                    When provided, matched mm positions are converted to voxel
%                    indices via ras_to_grid and returned in localite_vox.
%
% Output:
%   localite_mm    - [N x 3] EntryTarget xyz [mm] for each target_map entry.
%                    Rows without a match are all-NaN.
%   localite_vox   - [N x 3] Corresponding voxel coordinates (integer-rounded).
%                    Rows are all-NaN when unmatched or when t1_header is absent.
%
% See also: NEURONAV_READ_ENTRY_TARGET, RAS_TO_GRID

N = numel(target_map);
localite_mm  = nan(N, 3);
localite_vox = nan(N, 3);

if isempty(entry_targets)
    return;
end

for ti = 1:N
    tm = target_map(ti);

    if isfield(tm, 'localite_key') && ~isempty(tm.localite_key)
        raw_keys = tm.localite_key;
    else
        % Strip trailing numeric disambiguation suffix (e.g. 'L_M1_2' → 'L_M1')
        raw_keys = regexprep(strtrim(tm.target_name), '_\d+$', '');
    end
    if ischar(raw_keys) || isstring(raw_keys)
        keys = {char(raw_keys)};
    else
        keys = cellstr(raw_keys);
    end
    keys = upper(strtrim(keys));

    matched_description = '';
    matched_key = '';
    for ki = 1:numel(keys)
        for ei = 1:numel(entry_targets)
            if contains(upper(entry_targets(ei).description), keys{ki})
                localite_mm(ti,:) = entry_targets(ei).xyz;
                matched_description = entry_targets(ei).description;
                matched_key = keys{ki};
                break;
            end
        end
        if ~isempty(matched_description), break; end
    end

    if isfield(tm, 'transducer_id')
        coil_label = sprintf('coil %d', tm.transducer_id);
    else
        coil_label = 'coil ?';
    end
    keys_str = strjoin(keys, '" / "');
    if isempty(matched_description)
        fprintf('  [entry_target] %s (%s, search keys "%s"): no matching EntryTarget label\n', ...
            tm.target_name, coil_label, keys_str);
    else
        fprintf('  [entry_target] %s (%s, search keys "%s") -> matched label "%s" via key "%s"\n', ...
            tm.target_name, coil_label, keys_str, matched_description, matched_key);
    end
end

% Convert matched mm positions to voxel space
if nargin >= 3 && ~isempty(t1_header)
    for ti = 1:N
        if ~any(isnan(localite_mm(ti,:)))
            localite_vox(ti,:) = ras_to_grid(localite_mm(ti,:)', t1_header);
        end
    end
end
end
