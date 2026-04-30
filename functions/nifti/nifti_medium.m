function nifti_medium(parameters, planimg, medium_masks, kwave_medium, pseudoCT)
% NIFTI_MEDIUM  Export medium-stage data to NIfTI
%
% Pipeline outputs (dir_nii_T1w + optional MNI):
%   medium_masks  - voxel-wise tissue label map (always)
%   pseudoCT      - Hounsfield-unit skull map (only when parameters.pct.enabled == 1)
%
% Debug outputs (dir_debug_medium):
%   sound_speed, density, alpha_coeff, alpha_power, thermal_conductivity,
%   specific_heat, perfusion_coeff, absorption_fraction
%   All back-transformed to T1 space; skipped if file already exists.
%
% Use as:
%   nifti_medium(parameters, planimg, medium_masks, kwave_medium)
%   nifti_medium(parameters, planimg, medium_masks, kwave_medium, pseudoCT)
%
% Input:
%   parameters   - PRESTUS config struct
%   planimg      - struct with fields: t1_image_orig, inv_transf, t1_header, transf
%   medium_masks - (:,:,:) uint8, simulation-grid tissue labels
%   kwave_medium - struct with medium property fields
%   pseudoCT     - (:,:,:) numeric, simulation-grid Hounsfield values (optional)
%
% See also: NIFTI_ACOUSTIC, NIFTI_THERMAL, CONVERT_FINAL_TO_MNI_SIMNIBS

arguments
    parameters   (1,1) struct
    planimg      (1,1) struct
    medium_masks
    kwave_medium (1,1) struct
    pseudoCT               = []
end

    if isfield(parameters.modules, 'run_nifti_creation') && ...
            parameters.modules.run_nifti_creation == 0
        disp('No nifti creation requested...')
        return
    end

    if ~contains(parameters.simulation.medium, {'layered'; 'phantom'})
        return
    end

    is_layered = ~strcmp(parameters.simulation.medium, 'phantom');
    output_mni = ~isfield(parameters, 'analysis') || ...
                 ~isfield(parameters.analysis, 'output_mni') || ...
                 parameters.analysis.output_mni ~= 0;
    m2m_folder = fullfile(parameters.path.seg, sprintf('m2m_sub-%03d', parameters.subject_id));

    % ---- pipeline outputs ------------------------------------------------

    masks_file    = fullfile(parameters.io.dir_nii_T1w, sprintf('sub-%03d_%s_T1w%s_medium_masks', ...
        parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));
    masks_file_gz = strcat(masks_file, '.nii.gz');
    mni_masks_file = fullfile(parameters.io.dir_nii_MNI, sprintf('sub-%03d_%s_MNI%s_medium_masks.nii.gz', ...
        parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));

    nifti_write_volume(uint8(medium_masks), masks_file, parameters, planimg, ...
        'IsLayered', is_layered, 'Datatype', 'uint8', 'BitsPerPixel', 8, 'Resampler', 'nearest');
    nifti_to_mni(masks_file_gz, mni_masks_file, parameters, is_layered, output_mni, m2m_folder);

    % pseudoCT (pipeline output: goes to dir_nii_T1w)
    if is_layered && isfield(parameters, 'pct') && isfield(parameters.pct, 'enabled') && ...
            parameters.pct.enabled == 1 && ~isempty(pseudoCT)
        pct_file    = fullfile(parameters.io.dir_nii_T1w, sprintf('sub-%03d_%s_T1w%s_pseudoCT', ...
            parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));
        pct_file_gz = strcat(pct_file, '.nii.gz');
        pct_mni_file = fullfile(parameters.io.dir_nii_MNI, sprintf('sub-%03d_%s_MNI%s_pseudoCT.nii.gz', ...
            parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));

        nifti_write_volume(single(pseudoCT), pct_file, parameters, planimg, ...
            'IsLayered', is_layered, 'Resampler', 'nearest');
        nifti_to_mni(pct_file_gz, pct_mni_file, parameters, is_layered, output_mni, m2m_folder);
    end

    % ---- debug outputs ---------------------------------------------------

    if is_layered
        debug_properties = {'sound_speed', 'density', 'alpha_coeff', 'alpha_power', ...
                            'thermal_conductivity', 'specific_heat', 'perfusion_coeff', 'absorption_fraction'};
        try
            out_dir = char(parameters.io.dir_debug_medium);
            if ~isfolder(out_dir); mkdir(out_dir); end

            orig_hdr = planimg.t1_header;
            orig_hdr.Datatype = 'single';

            for prop = debug_properties
                prop = char(prop); %#ok<FXSET>
                if ~isfield(kwave_medium, prop)
                    continue
                end
                file_name = fullfile(out_dir, [prop '_t1']);
                if isfile([file_name '.nii.gz']) || isfile([file_name '.nii'])
                    continue
                end
                data_t1 = single(tformarray(kwave_medium.(prop), planimg.inv_transf, ...
                    makeresampler('nearest', 'fill'), [1 2 3], [1 2 3], ...
                    orig_hdr.ImageSize, [], 0));
                niftiwrite(data_t1, file_name, orig_hdr);
                nii_file = [file_name '.nii'];
                if isfile(nii_file)
                    gzip(nii_file);
                    delete(nii_file);
                end
            end
        catch ME
            prev = warning('off', 'backtrace');
            warn('nifti_medium:debugWrite', ...
                'Could not write medium property debug NIfTIs: %s', ME.message);
            warning(prev);
        end
    end

    % T1 in MNI space (charm does not produce one; written once, existence-guarded)
    if is_layered && output_mni
        path_to_input_img  = fullfile(m2m_folder, 'T1.nii.gz');
        path_to_output_img = fullfile(m2m_folder, 'toMNI', 'T1_to_MNI_post-hoc.nii.gz');
        if ~exist(path_to_output_img, 'file')
            convert_final_to_MNI_simnibs(path_to_input_img, m2m_folder, path_to_output_img, parameters);
        end
    end
end
