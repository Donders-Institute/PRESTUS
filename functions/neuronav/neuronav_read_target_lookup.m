function target_lookup_vox = neuronav_read_target_lookup(lookup_csv, sub_id, target_map)
% NEURONAV_READ_TARGET_LOOKUP  Look up pre-specified target coordinates from
% a study-supplied lookup table, keyed by subject and target name.
%
% Some studies establish the intended stimulation target BEFORE the session
% (e.g. via individual anatomical segmentation and an optimization search
% over candidate transducer placements), rather than only reconstructing it
% after the fact from neuronavigation exports (EntryTarget/GUMMarkers, which
% depend on the operator faithfully re-confirming/re-exporting it, and are
% simply absent for some marker-export conventions). When available, this
% pre-specified value is the highest-confidence target coordinate source and
% should be preferred; EntryTarget/GUMMarkers-derived values remain useful
% as an independent cross-check (see NEURONAV_MATCH_ENTRY_TARGETS) and as
% the fallback when no lookup entry exists for a given subject/target.
%
% Lookup CSV format (header row required):
%   subject_id, target_name, space, x, y, z  [, source_file]
%     subject_id  - e.g. 'sub-003' (case-insensitive match)
%     target_name - matches target_map(i).target_name, e.g. 'L_AMY'
%                   (case-insensitive match; trailing '_<N>' sequence
%                   suffixes are NOT stripped -- provide one row per exact
%                   target_name if a coil has multiple sequences)
%     space       - coordinate space tag. Only 'native_vox' (T1 grid voxel
%                   space, matching target_native_vox in the exported
%                   session CSV) is currently supported; other values are
%                   skipped with a warning.
%     x, y, z     - the coordinate, in the space above
%     source_file - optional, for provenance/audit only; not used here
%
% A study builds this CSV with its own helper script (e.g. the amythal
% study's code/build_target_lookup.m, which parses per-subject/per-target
% planning files into this generic format) -- this function only reads the
% already-unified table, so it stays study-agnostic.
%
% Use as:
%   target_lookup_vox = neuronav_read_target_lookup(lookup_csv, sub_id, target_map)
%
% Input:
%   lookup_csv - path to the lookup CSV. Pass '' or a nonexistent path to
%                get an all-NaN result (i.e. this is a no-op when a study
%                hasn't set one up).
%   sub_id     - subject identifier string, e.g. 'sub-003'
%   target_map - struct array (from NEURONAV_BUILD_TARGET_MAP), one entry
%                per target; only .target_name is used.
%
% Output:
%   target_lookup_vox - [Nx3], N = numel(target_map). NaN rows for any
%                        target_map entry with no matching CSV row (or when
%                        lookup_csv is empty/missing).
%
% See also: NEURONAV_MATCH_ENTRY_TARGETS, NEURONAV_READ_ENTRY_TARGET,
%           NEURONAV_BUILD_TARGET_MAP

arguments
    lookup_csv (1,:) char
    sub_id     (1,:) char
    target_map struct
end

N = numel(target_map);
target_lookup_vox = nan(N, 3);

if isempty(lookup_csv) || ~isfile(lookup_csv)
    return;
end

T = readtable(lookup_csv, 'TextType', 'string', 'Delimiter', ',');
required_cols = {'subject_id', 'target_name', 'space', 'x', 'y', 'z'};
missing_cols = required_cols(~ismember(required_cols, T.Properties.VariableNames));
if ~isempty(missing_cols)
    warning('neuronav_read_target_lookup:badSchema', ...
        'Lookup CSV %s is missing column(s): %s -- returning all-NaN.', ...
        lookup_csv, strjoin(missing_cols, ', '));
    return;
end

T = T(strcmpi(T.subject_id, sub_id), :);

for ti = 1:N
    tname = target_map(ti).target_name;
    row = T(strcmpi(T.target_name, tname), :);
    if isempty(row)
        continue;
    end
    if height(row) > 1
        warning('neuronav_read_target_lookup:duplicateRow', ...
            '%s %s: %d matching rows in lookup CSV; using the first.', ...
            sub_id, tname, height(row));
        row = row(1, :);
    end
    if ~strcmpi(row.space(1), 'native_vox')
        warning('neuronav_read_target_lookup:unsupportedSpace', ...
            '%s %s: lookup space "%s" not supported (only native_vox); skipping.', ...
            sub_id, tname, row.space(1));
        continue;
    end
    target_lookup_vox(ti, :) = [row.x(1), row.y(1), row.z(1)];
end
end
