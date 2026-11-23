function bone_file = pct_bone_mask_file(parameters, opts)
% PCT_BONE_MASK_FILE  This subject's bone mask file (pct.bone_source 'mask')
%
% path.bone_mask_pattern gives the file, as path.ct_pattern: sprintf
% specifiers (e.g. %03d for the subject) and shell-style globs, relative to
% path.anat unless it is an absolute path.
%
% Use as:
%   bone_file = pct_bone_mask_file(parameters)
%   bone_file = pct_bone_mask_file(parameters, MustExist=false)
%
% Input:
%   parameters - PRESTUS parameter struct (subject_id, path.anat, path.bone_mask_pattern)
%   MustExist  - false: return '' instead of erroring when the pattern is empty
%                or matches no file (default true)
%
% Output:
%   bone_file - full path of the bone mask. Errors when the pattern is empty or
%               matches no file.
%
% See also: PCT_LOAD_BONE_MASK, PCT_CT_FILE

arguments
    parameters (1,1) struct
    opts.MustExist (1,1) logical = true
end

    pattern = char(string(get_field(parameters, {'path', 'bone_mask_pattern'}, '')));
    bone_file = '';
    if isempty(pattern) && ~opts.MustExist
        return
    elseif isempty(pattern)
        error('PRESTUS:pct:noBoneMask', 'pct.bone_source ''mask'' needs path.bone_mask_pattern.');
    end
    pattern = sprintf(pattern, parameters.subject_id);
    if ~is_absolute(pattern)
        pattern = fullfile(char(parameters.path.anat), pattern);
    end
    d = dir(pattern);
    if isempty(d) && ~opts.MustExist
        return
    elseif isempty(d)
        error('PRESTUS:pct:boneMaskNotFound', 'No file matched path.bone_mask_pattern: %s', pattern);
    end
    if numel(d) > 1
        warning('PRESTUS:pct:boneMaskAmbiguous', '%d files match path.bone_mask_pattern (%s); using %s.', ...
            numel(d), pattern, d(1).name);
    end
    bone_file = fullfile(d(1).folder, d(1).name);
end

function tf = is_absolute(p)
    tf = startsWith(p, '/') || startsWith(p, '\') || ~isempty(regexp(p, '^[A-Za-z]:[\\/]', 'once'));
end
