function nifti_thermal(parameters, planimg, results_heating, kwave_medium)
% NIFTI_THERMAL  Export thermal-stage data to NIfTI in T1 and MNI space
%
% Writes heating maps, heat-rise maps, and CEM43 maps back-transformed to T1
% space, with edge-artefact correction applied to the temperature map.
% Optionally converts all outputs to MNI space.
%
% Use as:
%   nifti_thermal(parameters, planimg, results_heating, kwave_medium)
%
% Input:
%   parameters      - PRESTUS config struct
%   planimg         - struct with fields: t1_image_orig, inv_transf, t1_header, transf
%   results_heating - struct with fields: maxT, heating_endT, CEM43, CEM43_end,
%                     CEM43_iso, CEM43_iso_end
%   kwave_medium    - struct with field: temp_0
%
% See also: NIFTI_MEDIUM, NIFTI_ACOUSTIC, CONVERT_FINAL_TO_MNI_SIMNIBS

arguments
    parameters      (1,1) struct
    planimg         (1,1) struct
    results_heating (1,1) struct
    kwave_medium    (1,1) struct
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

    data_types = ["heating", "heating_end", "heatrise", "heatrise_end", ...
                  "CEM43", "CEM43_end", "CEM43_iso", "CEM43_iso_end"];

    for data_type = data_types
        orig_file    = fullfile(parameters.io.dir_nii_T1w, sprintf('sub-%03d_%s_T1w%s_%s', ...
            parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix, data_type));
        orig_file_gz = strcat(orig_file, '.nii.gz');
        mni_file     = fullfile(parameters.io.dir_nii_MNI, sprintf('sub-%03d_%s_MNI%s_%s.nii.gz', ...
            parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix, data_type));

        switch data_type
            case "heating",      data = single(results_heating.maxT);
            case "heating_end",  data = single(results_heating.heating_endT);
            case "heatrise",     data = single(results_heating.maxT - kwave_medium.temp_0);
            case "heatrise_end", data = single(results_heating.heating_endT - kwave_medium.temp_0);
            case "CEM43",        data = single(results_heating.CEM43);
            case "CEM43_end",    data = single(results_heating.CEM43_end);
            case "CEM43_iso",    data = single(results_heating.CEM43_iso);
            case "CEM43_iso_end",data = single(results_heating.CEM43_iso_end);
        end

        if confirm_overwriting(orig_file_gz, parameters)
            if is_layered
                orig_hdr = planimg.t1_header;
                orig_hdr.Datatype = 'single';
                data_backtransf = tformarray(data, planimg.inv_transf, ...
                    makeresampler('cubic', 'fill'), [1 2 3], [1 2 3], ...
                    size(planimg.t1_image_orig), [], 0);

                if strcmp(data_type, "heating")
                    % Remove back-transform edge artefacts: replace the outer
                    % shell of the heated region with the median boundary value.
                    heatmap_mask = data_backtransf > 0;
                    heatmap_mask_eroded = imerode(heatmap_mask, strel('sphere', 3));
                    heatmap_mask_eroded_edge = imdilate(heatmap_mask_eroded, strel('sphere', 1)) & heatmap_mask & ~heatmap_mask_eroded;
                    heatmap_band_values = data_backtransf(heatmap_mask_eroded_edge);
                    if ~isempty(heatmap_band_values)
                        heatmap_band_value = median(heatmap_band_values);
                    else
                        heatmap_band_value = parameters.thermal.temp_0.water;
                    end
                    data_backtransf(~heatmap_mask_eroded) = heatmap_band_value;
                    data_backtransf(1:2,:,:)       = parameters.thermal.temp_0.water;
                    data_backtransf(end-1:end,:,:) = parameters.thermal.temp_0.water;
                    data_backtransf(:,1:2,:)       = parameters.thermal.temp_0.water;
                    data_backtransf(:,end-1:end,:) = parameters.thermal.temp_0.water;
                    data_backtransf(:,:,1:2)       = parameters.thermal.temp_0.water;
                    data_backtransf(:,:,end-1:end) = parameters.thermal.temp_0.water;

                elseif any(strcmp(data_type, ["CEM43", "CEM43_end", "CEM43_iso", "CEM43_iso_end"]))
                    data_backtransf(data_backtransf <= 0) = 0.0000001;
                end

                niftiwrite(data_backtransf, orig_file, orig_hdr, 'Compressed', true);
            else
                niftiwrite(data, orig_file, 'Compressed', true);
            end
        end

        nifti_to_mni(orig_file_gz, mni_file, parameters, is_layered, output_mni, m2m_folder);

        clear data;
    end
end
