function candidates = neuronav_read_label_candidates(raw_root, sub_id, session, loc_folder)
% NEURONAV_READ_LABEL_CANDIDATES  Collect named planned targets for a session.
%
% Gathers every named, set target point Localite recorded for a session, as
% candidates for naming ingested series (see NEURONAV_SUGGEST_TARGET_NAMES):
%   1. the latest EntryTarget export (NEURONAV_READ_ENTRY_TARGET), and
%   2. set EntryTargetPair targets in any GUMMarkers export of the session
%      (the pair's description names the target).
% Targets are often re-planned during a session, so each label is kept once,
% at its most recent position: EntryTarget first (operator-confirmed), then
% GUMMarkers files from newest to oldest (timestamp in the file name).
%
% Use as:
%   candidates = neuronav_read_label_candidates(raw_root, sub_id, session)
%   candidates = neuronav_read_label_candidates(raw_root, sub_id, session, loc_folder)
%
% Input:
%   raw_root   - root folder containing <sub_id>/<session>/<loc_folder>/...
%   sub_id     - subject identifier, e.g. 'sub-001'
%   session    - session string ('ses-02') or number (2)
%   loc_folder - (optional) subfolder name, default 'loc'; '' = none
%
% Output:
%   candidates - struct array with fields
%                  .description (char)  target label
%                  .xyz         [1×3]   position, Localite planning-image mm
%                  .source      (char)  'EntryTarget <file>' / 'GUMMarkers <file>'
%
% See also: NEURONAV_SUGGEST_TARGET_NAMES, NEURONAV_READ_ENTRY_TARGET

arguments
    raw_root   {mustBeTextScalar}
    sub_id     {mustBeTextScalar}
    session
    loc_folder = 'loc'
end

    if nargin < 4, loc_folder = 'loc'; end   % '' = markers directly under ses-XX/
    if isnumeric(session), session = sprintf('ses-%02d', double(session)); end
    session = char(session);

    candidates = struct('description', {}, 'xyz', {}, 'source', {});
    loc_root = fullfile(raw_root, sub_id, session, loc_folder);
    if ~isfolder(loc_root), return; end

    % --- EntryTarget export -------------------------------------------------
    et_files = dir(fullfile(loc_root, '**', 'EntryTarget*.xml'));
    if ~isempty(et_files)
        [~, order] = sort({et_files.name});   % same pick as neuronav_read_entry_target
        k = order(end);
        et = neuronav_read_entry_target(raw_root, sub_id, session, loc_folder);
        for i = 1:numel(et)
            candidates(end+1) = struct('description', et(i).description, ...
                'xyz', double(et(i).xyz), 'source', ['EntryTarget ' et_files(k).name]); %#ok<AGROW>
        end
    end

    % --- EntryTargetPairs in GUMMarkers exports -----------------------------
    gum_files = dir(fullfile(loc_root, '**', 'GUMMarkers*.xml'));
    [~, order] = sort({gum_files.name});
    gum_files = gum_files(fliplr(order));   % newest first
    for f = 1:numel(gum_files)
        fpath = fullfile(gum_files(f).folder, gum_files(f).name);
        try
            g = readstruct(fpath);
        catch
            continue;
        end
        if ~isfield(g, 'Element'), continue; end
        elements = g.Element;
        if ~iscell(elements), elements = num2cell(elements); end
        for i = 1:numel(elements)
            e = elements{i};
            if ~(isfield(e, 'EntryTargetPair') && isstruct(e.EntryTargetPair)), continue; end
            etp = e.EntryTargetPair;
            desc = '';
            if isfield(etp, 'descriptionAttribute'), desc = strtrim(char(etp.descriptionAttribute)); end
            if isempty(desc) || ~isfield(etp, 'Target'), continue; end
            m = etp.Target.Marker;
            if ~(isfield(m, 'setAttribute') && strcmpi(string(m.setAttribute), "true")), continue; end
            cv  = m.ColVec3D;
            xyz = double([cv.data0Attribute, cv.data1Attribute, cv.data2Attribute]);
            if all(xyz == 0), continue; end
            candidates(end+1) = struct('description', desc, 'xyz', xyz, ...
                'source', ['GUMMarkers ' gum_files(f).name]); %#ok<AGROW>
        end
    end

    % --- One entry per label: first seen = highest priority / newest -------
    if numel(candidates) > 1
        [~, first] = unique({candidates.description}, 'stable');
        candidates = candidates(sort(first));
    end
end
