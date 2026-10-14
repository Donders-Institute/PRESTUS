function pct_dir = find_pseudoct_dir(parameters)
% FIND_PSEUDOCT_DIR  Folder holding this subject's existing pseudoCT.nii.gz
%
% GET_PCT_DIR when pseudoCT.nii.gz is there. Older versions wrote the
% pseudo-CT into the SimNIBS m2m_sub-NNN folder when path.pct was empty; such
% a file is still used, with a warning to move it, but nothing is written
% there any more. Returns '' when no pseudo-CT exists yet.
%
% Use as:
%   pct_dir = find_pseudoct_dir(parameters)
%
% Input:
%   parameters - PRESTUS parameter struct
%
% Output:
%   pct_dir - folder with pseudoCT.nii.gz (and tissues_mask.nii.gz), or ''
%
% See also: GET_PCT_DIR

arguments
    parameters (1,1) struct
end

    pct_dir = get_pct_dir(parameters);
    if isfile(fullfile(pct_dir, 'pseudoCT.nii.gz'))
        return
    end

    legacy = '';
    if isfield(parameters, 'path') && isfield(parameters.path, 'seg') && ~isempty(parameters.path.seg)
        legacy = fullfile(char(parameters.path.seg), sprintf('m2m_sub-%03d', parameters.subject_id));
    end
    if ~isempty(legacy) && isfile(fullfile(legacy, 'pseudoCT.nii.gz'))
        warning('PRESTUS:pct:legacyLocation', ...
            ['Using the pseudo-CT in the SimNIBS folder:\n  %s\n' ...
             'PRESTUS no longer writes there. Move pseudoCT.nii.gz, tissues_mask.nii.gz ' ...
             '(and pseudoCT.json, debug_pct/) to:\n  %s'], legacy, pct_dir);
        pct_dir = legacy;
        return
    end
    pct_dir = '';
end
