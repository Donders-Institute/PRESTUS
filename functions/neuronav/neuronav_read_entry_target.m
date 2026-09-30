function targets = neuronav_read_entry_target(raw_root, sub_id, session, loc_folder)
% NEURONAV_READ_ENTRY_TARGET  Read planned target positions from Localite EntryTarget XML.
%
% Finds the final (latest-timestamp) EntryTarget XML in the raw Localite
% session folder and returns all Target entries that have been explicitly
% set (set="true", non-zero coordinates). Coordinates are in the Localite
% planning image coordinate space (zero-origin mm, same as TriggerMarker
% matrix translations).
%
% Use as:
%   targets = neuronav_read_entry_target(raw_root, sub_id, session)
%   targets = neuronav_read_entry_target(raw_root, sub_id, session, loc_folder)
%
% Input:
%   raw_root   - root folder containing <sub_id>/<session>/<loc_folder>/...
%   sub_id     - subject identifier, e.g. 'sub-001'
%   session    - session identifier, e.g. 'ses-02'
%   loc_folder - (optional) subfolder name, default 'loc'; '' = none
%
% Output:
%   targets    - struct array, one entry per set Target, with fields:
%                  .index       (double) Localite index attribute (0-based)
%                  .description (char)   Marker description string
%                  .xyz         [1×3]    (x,y,z) in Localite planning image mm
%
% Notes:
%   The coordinates are in the planning T1 image space with origin at voxel
%   (0,0,0) and axes aligned to the image grid (world = voxel × voxel_size).
%   They are NOT in physical scanner RAS space. Use the planning T1 header
%   and ras_to_grid to convert to voxel indices if needed.

    if nargin < 4, loc_folder = 'loc'; end   % '' = markers directly under ses-XX/

    targets = struct('index', {}, 'description', {}, 'xyz', {});

    loc_root = fullfile(raw_root, sub_id, session, loc_folder);
    % Localite writes these into an EntryTarget/ folder; exports are also
    % copied flat into the loc folder, so search any depth.
    et_files = dir(fullfile(loc_root, '**', 'EntryTarget*.xml'));

    if isempty(et_files)
        warning('neuronav_read_entry_target: no EntryTarget files for %s %s in\n  %s', ...
            sub_id, session, loc_root);
        return;
    end

    % Sort by filename (timestamp is encoded in the name) and take the last file,
    % which reflects the final confirmed state of the session.
    [~, sort_idx] = sort({et_files.name});
    latest = et_files(sort_idx(end));
    fpath  = fullfile(latest.folder, latest.name);

    try
        s = readstruct(fpath);
    catch e
        warning('neuronav_read_entry_target: failed to read %s:\n  %s', fpath, e.message);
        return;
    end

    if ~isfield(s, 'Target')
        warning('neuronav_read_entry_target: no Target elements in %s', fpath);
        return;
    end

    target_list = s.Target;   % struct array (1 per <Target> element)
    for i = 1:numel(target_list)
        t = target_list(i);
        if ~isfield(t, 'Marker'), continue; end
        m = t.Marker;

        % Only include explicitly set targets
        is_set = isfield(m, 'setAttribute') && strcmpi(char(m.setAttribute), 'true');
        if ~is_set, continue; end

        % Description (required for label matching)
        desc = '';
        if isfield(m, 'descriptionAttribute')
            desc = strtrim(char(m.descriptionAttribute));
        end
        if isempty(desc), continue; end

        % Coordinates
        xyz = [NaN, NaN, NaN];
        if isfield(m, 'ColVec3D')
            cv  = m.ColVec3D;
            xyz = [cv.data0Attribute, cv.data1Attribute, cv.data2Attribute];
        end
        if all(xyz == 0), continue; end   % zeros = not truly set

        idx = i - 1;   % default 0-based index
        if isfield(t, 'indexAttribute'), idx = t.indexAttribute; end

        targets(end+1).index       = idx;
        targets(end).description   = desc;
        targets(end).xyz           = xyz;
    end

    fprintf('  [entry_target] %s %s: %d set targets read from %s\n', ...
        sub_id, session, numel(targets), latest.name);
    if ~isempty(targets)
        fprintf('    labels: %s\n', strjoin({targets.description}, ', '));
    end
end
