function f = study_run_files(subject_id, medium, affix, sim_dir)
% STUDY_RUN_FILES  Output file names of one run (subject + medium + affix)
%
% Single place that mirrors the per-subject naming of PATH_LOG_SETUP,
% NIFTI_ACOUSTIC / NIFTI_THERMAL, TARGETING_ANALYSIS and the per-subject
% report, for the group tools.
%
% Use as:
%   f = study_run_files(subject_id, medium, affix, sim_dir)
%
% Input:
%   subject_id - numeric subject ID
%   medium     - simulation.medium
%   affix      - sanitised io.output_affix ('' or '_desc-...')
%   sim_dir    - path.sim
%
% Output:
%   f - struct: sub_dir, csv, targeting_csv, report, img_dir, nii_dir, and
%       function handles mni(type) / t1w(type) giving the NIfTI path of a
%       data type (e.g. 'focal_mask_3dB', 'Ipa', 'heatrise')
%
% See also: STUDY_RUN_TABLE, CREATE_GROUP_MAPS, GENERATE_GROUP_REPORT

    medium = char(medium); affix = char(affix);
    sub = sprintf('sub-%03d', subject_id);
    f.sub_dir       = fullfile(char(sim_dir), sub);
    f.csv           = fullfile(f.sub_dir, sprintf('%s_%s%s.csv', sub, medium, affix));
    f.targeting_csv = fullfile(f.sub_dir, sprintf('%s_%s%s_targeting.csv', sub, medium, affix));
    f.report        = fullfile(f.sub_dir, sprintf('%s_%s_report%s.html', sub, medium, affix));
    f.img_dir       = fullfile(f.sub_dir, 'img');
    f.nii_dir       = fullfile(f.sub_dir, 'nii');
    nii_dir = f.nii_dir;
    f.mni = @(type) fullfile(nii_dir, sprintf('%s_%s_MNI%s_%s.nii.gz', sub, medium, affix, type));
    f.t1w = @(type) fullfile(nii_dir, sprintf('%s_%s_T1w%s_%s.nii.gz', sub, medium, affix, type));
end
