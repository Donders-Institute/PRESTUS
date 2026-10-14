function pct_dir = get_pct_dir(parameters)
% GET_PCT_DIR  Per-subject pseudo-CT folder
%
% The pseudo-CT is a subject derivative (it depends only on the subject's
% images), so it sits next to the SimNIBS segmentation, never inside it
% (see doc_outputs.md, "Where outputs live"):
%
%   <path.pct>/sub-NNN                      when path.pct is set
%   <parent of path.seg>/pct/sub-NNN        otherwise
%
% The single resolver for where the pseudo-CT is written and expected
% (PATH_LOG_SETUP stores it as io.dir_pct). FIND_PSEUDOCT_DIR also accepts a
% pseudo-CT that older versions wrote into the SimNIBS m2m folder.
%
% Use as:
%   pct_dir = get_pct_dir(parameters)
%
% Input:
%   parameters - PRESTUS parameter struct (subject_id, path.pct or path.seg)
%
% Output:
%   pct_dir - char path (not created)
%
% See also: FIND_PSEUDOCT_DIR, PCT_CREATE_PCT_RUN, PATH_LOG_SETUP

arguments
    parameters (1,1) struct
end

    sub = sprintf('sub-%03d', parameters.subject_id);
    if isfield(parameters, 'path') && isfield(parameters.path, 'pct') && ...
            (ischar(parameters.path.pct) || isstring(parameters.path.pct)) && ...
            strlength(string(parameters.path.pct)) > 0
        pct_dir = fullfile(char(parameters.path.pct), sub);
    elseif isfield(parameters, 'path') && isfield(parameters.path, 'seg') && ...
            (ischar(parameters.path.seg) || isstring(parameters.path.seg)) && ...
            strlength(string(parameters.path.seg)) > 0
        % fileparts of '/a/seg/' would give '/a/seg': strip a trailing separator
        seg = regexprep(char(parameters.path.seg), '[\\/]+$', '');
        pct_dir = fullfile(fileparts(seg), 'pct', sub);
    else
        error('get_pct_dir:missingPath', 'Set path.pct or path.seg to locate the pseudo-CT.');
    end
end
