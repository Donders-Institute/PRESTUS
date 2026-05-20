function nifti_acoustic(parameters, planimg, results_acoustic, acoustic_Ipa, acoustic_MI, acoustic_pressure, highlighted_pos)
% NIFTI_ACOUSTIC  Export acoustic-stage data to NIfTI in T1 and MNI space
%
% Writes intensity (Isppa), mechanical index, and pressure maps back-transformed
% to T1 space. For the intensity map, also generates an overlay plot over the T1
% image. Optionally converts all outputs to MNI space.
%
% The T1w back-transform (tformarray) is batched across all three volumes, and
% the MNI transform (subject2mni) is run in parallel.
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
    n = numel(data_types);

    % Build file paths up front
    orig_files    = cell(1, n);
    orig_files_gz = cell(1, n);
    mni_files     = cell(1, n);
    for k = 1:n
        dt = data_types(k);
        orig_files{k}    = fullfile(parameters.io.dir_nii_T1w, ...
            sprintf('sub-%03d_%s_T1w%s_%s', parameters.subject_id, ...
                    parameters.simulation.medium, parameters.io.output_affix, dt));
        orig_files_gz{k} = strcat(orig_files{k}, '.nii.gz');
        mni_files{k}     = fullfile(parameters.io.dir_nii_MNI, ...
            sprintf('sub-%03d_%s_MNI%s_%s.nii.gz', parameters.subject_id, ...
                    parameters.simulation.medium, parameters.io.output_affix, dt));
    end

    % -------------------------------------------------------------------------
    % T1w back-transform — single batched tformarray call across all 3 volumes
    % -------------------------------------------------------------------------
    if is_layered
        orig_hdr          = planimg.t1_header;
        orig_hdr.Datatype = 'single';

        stack = cat(4, single(acoustic_Ipa), single(acoustic_MI), single(acoustic_pressure));
        backtransf = tformarray(stack, planimg.inv_transf, ...
            makeresampler('cubic', 'fill'), [1 2 3], [1 2 3], ...
            size(planimg.t1_image_orig), [], 0);

        for k = 1:n
            orig_file_gz = orig_files_gz{k};
            if ~confirm_overwriting(orig_file_gz, parameters); continue; end
            niftiwrite(backtransf(:,:,:,k), orig_files{k}, orig_hdr, 'Compressed', true);
        end
    else
        data_list = {single(acoustic_Ipa), single(acoustic_MI), single(acoustic_pressure)};
        for k = 1:n
            if ~confirm_overwriting(orig_files_gz{k}, parameters); continue; end
            nifti_to_t1w(data_list{k}, orig_files{k}, parameters, planimg, 'IsLayered', false);
        end
    end

    % Overlay plot on T1 for intensity (needs the file to exist)
    if is_layered && isfile(orig_files_gz{1})
        plot_intensity_t1_overlay(niftiread(orig_files_gz{1}), planimg, parameters, ...
            results_acoustic, highlighted_pos);
    end

    % -------------------------------------------------------------------------
    % MNI transform — batch all three volumes into a single parallel run
    % -------------------------------------------------------------------------
    if ~is_layered || ~should_save_output(parameters.io, 'save_MNI')
        return
    end

    mni_in  = {};
    mni_out = {};
    for k = 1:n
        if ~confirm_overwriting(mni_files{k}, parameters); continue; end
        if ~isfile(orig_files_gz{k}); continue; end
        mni_in{end+1}  = orig_files_gz{k}; %#ok<AGROW>
        mni_out{end+1} = mni_files{k};     %#ok<AGROW>
    end

    if ~isempty(mni_in)
        convert_final_to_MNI_simnibs(mni_in, m2m_folder, mni_out, parameters);
    end
end
