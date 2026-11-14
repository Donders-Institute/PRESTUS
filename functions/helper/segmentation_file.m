function f = segmentation_file(parameters)
% SEGMENTATION_FILE  The SimNIBS segmentation of this subject (final_tissues.nii.gz)
%
% Use as:
%   f = segmentation_file(parameters)
%
% Input:
%   parameters - PRESTUS parameter struct (subject_id, path.seg)
%
% Output:
%   f - <path.seg>/m2m_sub-XXX/final_tissues.nii.gz, '' when path.seg or
%       subject_id is not set (the file need not exist)
%
% See also: SEGMENTATION_FINGERPRINT, PREPROC_SEGMENTATION

    f = '';
    if ~isfield(parameters, 'path') || ~isfield(parameters.path, 'seg') || isempty(parameters.path.seg) ...
            || ~isfield(parameters, 'subject_id') || isempty(parameters.subject_id)
        return
    end
    f = fullfile(char(parameters.path.seg), sprintf('m2m_sub-%03d', parameters.subject_id), ...
        'final_tissues.nii.gz');
end
