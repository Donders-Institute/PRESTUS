function f = resolve_profile_file(rel, equip, serial, user_libraries, profile_root)
% RESOLVE_PROFILE_FILE  Hydrophone-profile CSV named by an equipment entry's char_data_path
%
% A relative char_data_path is looked up, in this order, in
%   1. profiles/ next to the equipment entry (and that folder itself), for an
%      entry from your own library
%   2. profiles/ of each of your own libraries (most specific first)
%   3. profile_root (default <PRESTUS>/local/profiles)
% An absolute path is used as is.
%
% Use as:
%   f = resolve_profile_file(rel, equip, serial)
%   f = resolve_profile_file(rel, equip, serial, user_libraries, profile_root)
%
% Output:
%   f - the existing file, or '' when none is found
%
% See also: TRANSDUCER_LIBRARY_FOLDERS, LOAD_EQUIPMENT_CONFIG

arguments
    rel {mustBeTextScalar}
    equip
    serial {mustBeTextScalar}
    user_libraries = {}
    profile_root = fullfile(get_prestus_path(), 'local', 'profiles')
end

    if nargin < 4 || isempty(user_libraries); user_libraries = transducer_library_folders(); end
    if nargin < 5 || isempty(profile_root); profile_root = fullfile(get_prestus_path(), 'local', 'profiles'); end
    rel = char(rel);
    f = '';
    if isempty(rel); return; end
    if java.io.File(rel).isAbsolute()
        if isfile(rel); f = rel; end
        return;
    end
    cands = {};
    if ~isempty(equip) && isfield(equip, 'files') && isfield(equip.files, serial) && ...
            isfield(equip, 'user') && isfield(equip.user, serial) && equip.user.(serial)
        d = fileparts(equip.files.(serial));
        cands = [cands, {fullfile(d, 'profiles', rel), fullfile(d, rel)}];
    end
    for k = 1:numel(user_libraries)
        cands{end+1} = fullfile(user_libraries{k}, 'profiles', rel); %#ok<AGROW>
    end
    cands{end+1} = fullfile(profile_root, rel);
    hit = find(cellfun(@isfile, cands), 1);
    if ~isempty(hit); f = cands{hit}; end
end
