function [target_vox, target_mm, target_source] = neuronav_resolve_target( ...
        parameters, sub_id, session, sessions, target_map, t1_header, ...
        target_pos_projected, target_lookup_csv)
% NEURONAV_RESOLVE_TARGET  Resolve each target_map entry's brain-target
% coordinate from the best available source, in order of confidence.
%
% Multi-target Localite studies (see doc_placement_neuronav_multitarget.md)
% need one confirmed brain-target coordinate per named target
% (e.g. 'L_AMY', 'R_MD') per session. No single source is reliably present
% for every subject/session, so this tries four, in order, and records
% which one actually supplied each row:
%
%   1. target_lookup   - a study-supplied pre-hoc lookup table (the target
%                         as it was actually planned/known to the
%                         experimenter before the session -- see
%                         NEURONAV_READ_TARGET_LOOKUP). Highest confidence:
%                         not a post-hoc reconstruction at all. When this
%                         source and a lower tier both produce a value for
%                         the same target, the lookup value is used and the
%                         agreement/disagreement (in voxels) is logged as a
%                         cross-check rather than silently resolved.
%   2. entrytarget_current  - this session's own operator-confirmed
%                         EntryTarget export (NEURONAV_READ_ENTRY_TARGET +
%                         NEURONAV_MATCH_ENTRY_TARGETS).
%   3. entrytarget:<ses>    - another session's EntryTarget export, walked
%                         backward from the current session toward the
%                         first (some sessions never export their own --
%                         e.g. the operator only re-confirms targets on
%                         some visits -- even though that session's own
%                         entry-point triggers were recorded independently;
%                         trans_pos always comes from the current session's
%                         own triggers and is never affected by this).
%   4. projected       - transducer-axis-projected geometric estimate
%                         (target_pos_projected, computed by
%                         NEURONAV_CONVERT_TRIGGER_TO_VOXELS from this
%                         session's own tracked position). No operator
%                         confirmation exists anywhere for this
%                         subject/target (e.g. sessions that only ever used
%                         GUMMarkers, which has no separate target
%                         confirmation) -- always available, so target_vox
%                         is never NaN, but treat rows on this tier as
%                         estimates, not ground truth.
%
% Use as:
%   [target_vox, target_mm, target_source] = neuronav_resolve_target( ...
%       parameters, sub_id, session, sessions, target_map, t1_header, ...
%       target_pos_projected, target_lookup_csv)
%
% Input:
%   parameters      - PRESTUS config struct; uses parameters.path.localite_raw
%                      and parameters.placement.localite.loc_folder (default 'loc')
%   sub_id          - subject identifier, e.g. 'sub-003'
%   session         - current session identifier, e.g. 'ses-04'
%   sessions        - {1xM} cell array of ALL of this subject's session
%                      identifiers in chronological order, e.g.
%                      {'ses-02','ses-03','ses-04'} -- defines the walk-back
%                      order for tier 3
%   target_map      - struct array from NEURONAV_BUILD_TARGET_MAP
%   t1_header       - niftiinfo struct for the planning T1 (for EntryTarget
%                      mm -> vox conversion), or [] if unavailable
%   target_pos_projected - [Nx3] voxel-space projected target for every
%                      target_map entry (always available; tier 4)
%   target_lookup_csv - path to a lookup CSV (NEURONAV_READ_TARGET_LOOKUP
%                      format), or '' / a nonexistent path to skip tier 1
%
% Output:
%   target_vox    - [Nx3] resolved target, native T1 voxel space. Never
%                    NaN (see tier 4). This is what downstream simulation
%                    scripts should read as focus_pos.
%   target_mm     - [Nx3] operator-confirmed target in mm, from EntryTarget
%                    matching only (tiers 2-3). NaN wherever no operator
%                    confirmation exists -- deliberately NOT filled by
%                    tiers 1 or 4, so it stays a reliable "was this
%                    actually confirmed" signal independent of target_source.
%   target_source - {Nx1} cell array of strings, one of 'target_lookup',
%                    'entrytarget_current', 'entrytarget:<session>', or
%                    'projected'. Write this to your output table --
%                    checking it is the only reliable way to know whether
%                    target_vox is ground truth or an estimate.
%
% See also: NEURONAV_READ_TARGET_LOOKUP, NEURONAV_READ_ENTRY_TARGET,
%           NEURONAV_MATCH_ENTRY_TARGETS, NEURONAV_BUILD_TARGET_MAP,
%           NEURONAV_CONVERT_TRIGGER_TO_VOXELS, NEURONAV_EXPORT_SESSION_CSV

arguments
    parameters struct
    sub_id (1,:) char
    session (1,:) char
    sessions cell
    target_map struct
    t1_header
    target_pos_projected (:,3) double
    target_lookup_csv (1,:) char = ''
end

loc_folder_name = 'loc';
if isfield(parameters, 'placement') && isfield(parameters.placement, 'localite') && ...
        isfield(parameters.placement.localite, 'loc_folder')
    loc_folder_name = parameters.placement.localite.loc_folder;
end

% ---- Tiers 2-3: EntryTarget, current session then walking back ----------
entry_target_session = session;
entry_targets = neuronav_read_entry_target( ...
    parameters.path.localite_raw, sub_id, session, loc_folder_name);

if isempty(entry_targets)
    this_idx = find(strcmp(sessions, session), 1);
    for fb_idx = (this_idx - 1):-1:1
        fb_session = sessions{fb_idx};
        fprintf('  [entry_target] %s %s: no local EntryTarget export, trying %s ...\n', ...
            sub_id, session, fb_session);
        entry_targets = neuronav_read_entry_target( ...
            parameters.path.localite_raw, sub_id, fb_session, loc_folder_name);
        if ~isempty(entry_targets)
            fprintf('  [entry_target] %s %s: using target confirmation from %s\n', ...
                sub_id, session, fb_session);
            entry_target_session = fb_session;
            break;
        end
    end
end

[target_mm, target_vox] = neuronav_match_entry_targets(entry_targets, target_map, t1_header);

n_targets = size(target_vox, 1);
matched_mask = ~any(isnan(target_vox), 2);
target_source = repmat({''}, n_targets, 1);
if strcmp(entry_target_session, session)
    target_source(matched_mask) = {'entrytarget_current'};
else
    target_source(matched_mask) = {sprintf('entrytarget:%s', entry_target_session)};
end

% ---- Tier 4: axis-projected estimate, for anything still unmatched ------
no_match = ~matched_mask;
if any(no_match)
    fprintf('  [entry_target] %s %s: %d/%d target(s) unconfirmed -- using axis-projected position\n', ...
        sub_id, session, sum(no_match), numel(no_match));
    target_vox(no_match, :) = target_pos_projected(no_match, :);
    target_source(no_match) = {'projected'};
end
% target_mm is deliberately left NaN for no_match rows -- see docstring.

% ---- Tier 1: pre-hoc lookup, overrides everything above, cross-checked --
target_lookup_vox = neuronav_read_target_lookup(target_lookup_csv, sub_id, target_map);
has_lookup = ~any(isnan(target_lookup_vox), 2);
for ti = find(has_lookup)'
    tname = target_map(ti).target_name;
    if matched_mask(ti) || no_match(ti)
        delta = norm(target_lookup_vox(ti,:) - target_vox(ti,:));
        if delta > 2
            fprintf(['  [target_lookup] %s %s %s: MISMATCH vs %s by %.1f vox ' ...
                '(lookup=[%s] vs %s=[%s]) -- using lookup value\n'], ...
                sub_id, session, tname, target_source{ti}, delta, ...
                num2str(target_lookup_vox(ti,:)), target_source{ti}, num2str(target_vox(ti,:)));
        else
            fprintf('  [target_lookup] %s %s %s: confirms %s (Delta=%.1f vox)\n', ...
                sub_id, session, tname, target_source{ti}, delta);
        end
    end
    target_vox(ti, :) = target_lookup_vox(ti, :);
    target_source{ti} = 'target_lookup';
end
end
