function [mni_file, m2m_file] = get_t1_mni_file(seg_dir, subject_id)
% GET_T1_MNI_FILE  Path of a subject's T1-in-MNI image, outside the SimNIBS m2m folder
%
% The T1-in-MNI image is a subject derivative, so like the pseudo-CT
% (GET_PCT_DIR) it sits next to the SimNIBS segmentation, never inside it:
%
%   <parent of seg_dir>/mni/sub-NNN/T1_to_MNI_post-hoc.nii.gz
%
% The second output is the location inside the m2m folder used by earlier
% outputs; read it only, never write there.
%
% Use as:
%   mni_file = get_t1_mni_file(parameters.path.seg, parameters.subject_id)
%   [mni_file, m2m_file] = get_t1_mni_file(seg_dir, subject_id)
%
% Input:
%   seg_dir    - folder with the m2m_sub-NNN folders (path.seg)
%   subject_id - numeric subject ID
%
% Output:
%   mni_file - char path where the T1-in-MNI image is written (not created)
%   m2m_file - char path of the same image inside m2m_sub-NNN/toMNI (read only)
%
% See also: SIMNIBS_T1_TO_MNI, GET_PCT_DIR, CREATE_GROUP_MAPS

arguments
    seg_dir {mustBeTextScalar, mustBeNonzeroLengthText}
    subject_id (1,1) double
end

    sub  = sprintf('sub-%03d', subject_id);
    name = 'T1_to_MNI_post-hoc.nii.gz';
    % fileparts of '/a/seg/' would give '/a/seg': strip a trailing separator
    seg = regexprep(char(seg_dir), '[\\/]+$', '');
    mni_file = fullfile(fileparts(seg), 'mni', sub, name);
    m2m_file = fullfile(seg, ['m2m_' sub], 'toMNI', name);
end
