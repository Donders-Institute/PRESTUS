function nifti_acoustic(parameters, planimg, results_acoustic, acoustic_Ipa, acoustic_MI, acoustic_pressure, highlighted_pos)
% NIFTI_ACOUSTIC  Export acoustic-stage data to NIfTI in T1 and MNI space
%
% Writes intensity (Isppa), mechanical index, and pressure maps back-transformed
% to T1 space. For the intensity map, also generates an overlay plot over the T1
% image. Optionally converts all outputs to MNI space.
%
% Use as:
%   nifti_acoustic(parameters, planimg, results_acoustic, ...
%                  acoustic_Ipa, acoustic_MI, acoustic_pressure, highlighted_pos)
%
% Input:
%   parameters       - PRESTUS config struct
%   planimg          - struct with fields: t1_image_orig, inv_transf, t1_header, transf
%   results_acoustic - struct with fields: Isppa, Isppa_brain, …
%   acoustic_Ipa     - (:,:,:) single, Isppa map [W/cm²]
%   acoustic_MI      - (:,:,:) single, mechanical index map [-]
%   acoustic_pressure- (:,:,:) single, temporal peak pressure map [Pa]
%   highlighted_pos  - [1x3] grid-space peak intensity position [voxels]
%
% See also: NIFTI_MEDIUM, NIFTI_THERMAL, CONVERT_FINAL_TO_MNI_SIMNIBS

arguments
    parameters        (1,1) struct
    planimg           (1,1) struct
    results_acoustic  (1,1) struct
    acoustic_Ipa
    acoustic_MI
    acoustic_pressure
    highlighted_pos   (1,:) {mustBeNumeric}
end

    if isfield(parameters.modules, 'run_nifti_creation') && ...
            parameters.modules.run_nifti_creation == 0
        disp('No nifti creation requested...')
        return
    end

    if ~contains(parameters.simulation.medium, {'layered'; 'phantom'})
        return
    end

    is_layered = strcmp(parameters.simulation.medium, 'layered');
    m2m_folder = fullfile(parameters.path.seg, sprintf('m2m_sub-%03d', parameters.subject_id));

    data_types = ["intensity", "MI", "pressure"];

    for data_type = data_types
        orig_file    = fullfile(parameters.io.dir_nii_T1w, sprintf('sub-%03d_%s_T1w%s_%s', ...
            parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix, data_type));
        orig_file_gz = strcat(orig_file, '.nii.gz');
        mni_file     = fullfile(parameters.io.dir_nii_MNI, sprintf('sub-%03d_%s_MNI%s_%s.nii.gz', ...
            parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix, data_type));

        switch data_type
            case "intensity", data = single(acoustic_Ipa);
            case "MI",        data = single(acoustic_MI);
            case "pressure",  data = single(acoustic_pressure);
        end

        nifti_to_t1w(data, orig_file, parameters, planimg, 'IsLayered', is_layered);
        nifti_to_mni(orig_file_gz, mni_file, parameters, is_layered, m2m_folder);

        if strcmp(data_type, "intensity") && is_layered && isfile(orig_file_gz)
            plot_intensity_t1_overlay(niftiread(orig_file_gz), planimg, parameters, results_acoustic, highlighted_pos);
        end

        clear data;
    end
end
