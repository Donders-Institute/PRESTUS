function simulation_nifti(stage, parameters, planimg, varargin)
% SIMULATION_NIFTI  Export k-Wave simulation results to NIfTI in native and MNI space
%
% Called once per pipeline stage; each stage writes only the data available
% at that point. Back-transforms to T1 space via tformarray and optionally
% converts to MNI via SimNIBS subject2mni. Skips silently when
% parameters.modules.run_nifti_creation == 0.
%
% Use as:
%   simulation_nifti('medium',   parameters, planimg, medium_masks, kwave_medium)
%   simulation_nifti('acoustic', parameters, planimg, results_acoustic, ...
%                    acoustic_Ipa, acoustic_MI, acoustic_pressure, highlighted_pos)
%   simulation_nifti('thermal',  parameters, planimg, results_heating, kwave_medium)
%
% Input:
%   stage           - 'medium' | 'acoustic' | 'thermal'
%   parameters      - (1,1) struct, PRESTUS config
%   planimg         - (1,1) struct with fields: t1_image_orig, inv_transf, t1_header, transf
%
%   Stage 'medium':
%     medium_masks  - (:,:,:) uint8, voxel-wise tissue label map
%     kwave_medium  - (1,1) struct with fields: sound_speed, density, alpha_coeff
%
%   Stage 'acoustic':
%     results_acoustic - (1,1) struct (Isppa_brain, Isppa, …)
%     acoustic_Ipa     - (:,:,:) single, Isppa map [W/cm²]
%     acoustic_MI      - (:,:,:) single, mechanical index map [-]
%     acoustic_pressure- (:,:,:) single, temporal peak pressure map [Pa]
%     highlighted_pos  - [1x3] grid-space position of peak intensity [voxels]
%
%   Stage 'thermal':
%     results_heating - (1,1) struct (maxT, heating_endT, CEM43, …)
%     kwave_medium    - (1,1) struct with field: temp_0
%
% See also: ACOUSTIC_ANALYSIS, CONVERT_FINAL_TO_MNI_SIMNIBS, PLOT_OVERLAY

    if isfield(parameters.modules, 'run_nifti_creation') && ...
            parameters.modules.run_nifti_creation == 0
        disp('No nifti creation requested...')
        return
    end

    if ~contains(parameters.simulation.medium, {'layered'; 'phantom'})
        return
    end

    switch stage
        case 'medium'
            [medium_masks, kwave_medium] = varargin{:};
            data_types = ["medium_masks", "sound_speed", "density", "alpha_coeff"];

        case 'acoustic'
            [results_acoustic, acoustic_Ipa, acoustic_MI, acoustic_pressure, highlighted_pos] = varargin{:};
            data_types = ["intensity", "MI", "pressure"];

        case 'thermal'
            [results_heating, kwave_medium] = varargin{:};
            data_types = ["heating", "heating_end", "heatrise", "heatrise_end", ...
                          "CEM43", "CEM43_end", "CEM43_iso", "CEM43_iso_end"];

        otherwise
            error('simulation_nifti: unknown stage "%s"', stage);
    end

    for data_type = data_types
        orig_file = fullfile(parameters.io.dir_nii_T1w, sprintf('sub-%03d_%s_T1w%s_%s',...
            parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix, data_type));
        mni_file  = fullfile(parameters.io.dir_nii_MNI, sprintf('sub-%03d_%s_MNI%s_%s.nii.gz',...
            parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix, data_type));

        switch data_type
            case "medium_masks",   data = medium_masks;
            case "sound_speed",    data = single(kwave_medium.sound_speed);
            case "density",        data = single(kwave_medium.density);
            case "alpha_coeff",    data = single(kwave_medium.alpha_coeff);
            case "intensity",      data = single(acoustic_Ipa);
            case "MI",             data = single(acoustic_MI);
            case "pressure",       data = single(acoustic_pressure);
            case "heating",        data = single(results_heating.maxT);
            case "heating_end",    data = single(results_heating.heating_endT);
            case "heatrise",       data = single(results_heating.maxT - kwave_medium.temp_0);
            case "heatrise_end",   data = single(results_heating.heating_endT - kwave_medium.temp_0);
            case "CEM43",          data = single(results_heating.CEM43);
            case "CEM43_end",      data = single(results_heating.CEM43_end);
            case "CEM43_iso",      data = single(results_heating.CEM43_iso);
            case "CEM43_iso_end",  data = single(results_heating.CEM43_iso_end);
        end
        orig_file_with_ext = strcat(orig_file, '.nii.gz');

        if confirm_overwriting(orig_file_with_ext, parameters)
            if ~strcmp(parameters.simulation.medium, 'phantom')
                orig_hdr = planimg.t1_header;

                if strcmp(data_type, "medium_masks")
                    data_backtransf = tformarray(uint8(data), planimg.inv_transf, ...
                                            makeresampler('nearest', 'fill'), [1 2 3], [1 2 3], size(planimg.t1_image_orig), [], 0);
                    orig_hdr.Datatype = 'uint8';
                    orig_hdr.BitsPerPixel = 8;
                elseif any(strcmp(data_type, ["sound_speed", "density", "alpha_coeff"]))
                    % Nearest-neighbour preserves piecewise-constant tissue boundaries;
                    % fill with the water background value so edges are physically correct.
                    switch data_type
                        case "sound_speed",  fill_val = parameters.medium_properties.water.sound_speed;
                        case "density",      fill_val = parameters.medium_properties.water.density;
                        case "alpha_coeff",  fill_val = parameters.medium_properties.water.alpha_coeff;
                    end
                    data_backtransf = tformarray(data, planimg.inv_transf, ...
                                            makeresampler('nearest', 'fill'), [1 2 3], [1 2 3], size(planimg.t1_image_orig), [], fill_val);
                    orig_hdr.Datatype = 'single';
                else
                    data_backtransf = tformarray(data, planimg.inv_transf, ...
                                            makeresampler('cubic', 'fill'), [1 2 3], [1 2 3], size(planimg.t1_image_orig), [], 0);
                    orig_hdr.Datatype = 'single';

                    if strcmp(data_type, "heating")
                        % Remove back-transform edge artifacts
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
                end
                niftiwrite(data_backtransf, orig_file, orig_hdr, 'Compressed', true);
            else
                niftiwrite(data, orig_file, 'Compressed', true);
            end
        else
            data_backtransf = niftiread(orig_file_with_ext);
        end

        if strcmp(data_type, "intensity") && ~strcmp(parameters.simulation.medium, 'phantom')

            max_plots = min(2, numel(parameters.transducer));
            if numel(parameters.transducer) > max_plots
                warn('More than two transducers: intensity-over-T1 plots will be created only for the first 2 transducers');
            end

            if isfield(results_acoustic, 'Isppa_brain') && ~isempty(results_acoustic.Isppa_brain) && ~isnan(results_acoustic.Isppa_brain)
                max_val = results_acoustic.Isppa_brain;
            else
                max_val = results_acoustic.Isppa;
            end

            for ti = 1:max_plots
                tpos_sim = parameters.transducer(ti).trans_pos;
                fpos_sim = parameters.transducer(ti).focus_pos;

                backtransf_coordinates = round(tformfwd(...
                    [tpos_sim; fpos_sim; highlighted_pos], ...
                    planimg.inv_transf));

                [~, source_labels] = transducer_setup(...
                    parameters.transducer(1), ...
                    backtransf_coordinates(1,:), ...
                    backtransf_coordinates(2,:), ...
                    size(planimg.t1_image_orig), ...
                    planimg.t1_header.PixelDimensions(1), parameters);

                if ~isempty(planimg.t1_header) && isfield(planimg.t1_header, 'Transform')
                    [do_t, flip_r, flip_c] = nifti_slice_orientation(planimg.t1_header, 'y');
                else
                    [do_t, flip_r, flip_c] = deal(false, false, false);
                end

                [~,~,~,~,~,~,~,h]=plot_overlay(...
                    data_backtransf, ...
                    planimg.t1_image_orig, ...
                    source_labels, ...
                    parameters, ...
                    {'y', backtransf_coordinates(2,2)}, ...
                    backtransf_coordinates(1,:), ...
                    backtransf_coordinates(2,:), ...
                    backtransf_coordinates(3,:), ...
                    'show_rectangles', 0, ...
                    'grid_step', planimg.t1_header.PixelDimensions(1), ...
                    'overlay_color_range', [0 max_val], ...
                    'overlay_threshold_low', 0, ...
                    'overlay_threshold_high', max_val, ...
                    'rotation', 0, ...
                    'do_transpose', do_t, ...
                    'flip_rows', flip_r, ...
                    'flip_cols', flip_c);

                trans_suffix = '';
                if max_plots > 1; trans_suffix = sprintf('_T%02d', ti); end
                output_plot_filename = fullfile(parameters.io.dir_img, ...
                    sprintf('sub-%03d_%s%s_intensity_t1%s.png', ...
                    parameters.subject_id, parameters.simulation.medium, ...
                    parameters.io.output_affix, trans_suffix));
                saveas(h, output_plot_filename, 'png')
                close(h);
            end
        end

        m2m_folder = fullfile(parameters.path.seg, sprintf('m2m_sub-%03d', parameters.subject_id));

        output_mni = ~isfield(parameters, 'analysis') || ...
                     ~isfield(parameters.analysis, 'output_mni') || ...
                     parameters.analysis.output_mni ~= 0;
        if ~output_mni || ~confirm_overwriting(mni_file, parameters) || strcmp(parameters.simulation.medium, 'phantom')
            continue
        end
        convert_final_to_MNI_simnibs(orig_file_with_ext, m2m_folder, mni_file, parameters, 'interpolation_order', 0);

        clear data;
    end

    % Since charm does not transform the T1 into MNI space, one is manually created here.
    % Only needs to run once; guard with existence check.
    output_mni = ~isfield(parameters, 'analysis') || ...
                 ~isfield(parameters.analysis, 'output_mni') || ...
                 parameters.analysis.output_mni ~= 0;
    if output_mni && ~strcmp(parameters.simulation.medium, 'phantom')
        m2m_folder = fullfile(parameters.path.seg, sprintf('m2m_sub-%03d', parameters.subject_id));
        path_to_input_img  = fullfile(m2m_folder, 'T1.nii.gz');
        path_to_output_img = fullfile(m2m_folder, 'toMNI', 'T1_to_MNI_post-hoc.nii.gz');
        if ~exist(path_to_output_img, 'file')
            convert_final_to_MNI_simnibs(path_to_input_img, m2m_folder, path_to_output_img, parameters);
        end
    end
end
