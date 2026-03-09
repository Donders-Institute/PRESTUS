function kwave_medium = medium_setup(parameters, medium_masks, planimg, pseudoCT)

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %                         Setup k-wave medium                       %
    %                                                                   %
    % This function sets up the medium for k-wave simulations to take   %
    % place in. The way this mask is set up is by starting with an      %
    % empty grid and adjust the acoustic and thermal parameters in the  %
    % shape of the provided mask. This can be repeated to incorporate   %
    % different kinds of tissue (see 'layered').                        %
    %                                                                   %
    % Some notes:                                                       %
    % - The alpha value is converted at the end of this script.         %
    % - pseudoCT is an optional input when using pseudoCTs              %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    % Set residual unmapped medium voxels to water
    medium_masks(medium_masks==0) = find(strcmp(fieldnames(parameters.medium), 'water'));
    
    % Loads the medium settings from the config file
    medium = parameters.medium;
    
    % Create empty matrices for medium properties
    empty_grid = NaN(parameters.grid_dims);
    sound_speed = empty_grid; 
    density = empty_grid;
    alpha_coeff = empty_grid;
    alpha_power = empty_grid;
    thermal_conductivity = empty_grid;
    specific_heat = empty_grid;
    perfusion = empty_grid;
    absorption_fraction = ones(parameters.grid_dims); % default: convert 100% of attenuation into absorption
    temp_0 = empty_grid;

    % Get layer and medium labels
    layer_labels = fieldnames(parameters.layers);
    medium_labels = fieldnames(parameters.medium);

    % Iterate through each layer & assign medium ID to create medium mask
    for label_i = 1:length(layer_labels)
        label_name = layer_labels{label_i};
        medium_i = find(strcmp(medium_labels, label_name));
        if parameters.use_pseudoCT == 1 && strcmp(label_name, 'skull')
            skull_idx = find(ismember(medium_masks,medium_i));

            % set skull thermal conductivity
            thermal_conductivity(skull_idx) = medium.(label_name).thermal_conductivity;

            % set skull heat capacity     
            specific_heat(skull_idx) = medium.(label_name).specific_heat_capacity;

            % set skull perfusion
            if isfield(medium.(label_name), 'perfusion')
                perfusion(skull_idx) = medium.(label_name).perfusion;
            end

            % set skull absorption fraction
            if isfield(medium.(label_name), 'absorption_fraction')
                absorption_fraction(skull_idx) = medium.(label_name).absorption_fraction;
            end

            % set skull starting temperature
            if isfield(parameters.thermal.temp_0, label_name)
                temp_0(skull_idx) = parameters.thermal.temp_0.(label_name);
            end

            % map skull bone density with the desired algorithm
            if isfield(parameters, "pct_mapping_density")
                pct_mapping_density = parameters.pct_mapping_density;
            else
                pct_mapping_density = 'none';
            end
            [density] = medium_pct_density(parameters, medium, density, pseudoCT, skull_idx, pct_mapping_density);

            % map skull bone sound speed with the desired algorithm
            if isfield(parameters, "pct_mapping_soundspeed")
                pct_mapping_soundspeed = parameters.pct_mapping_soundspeed;
            else
                pct_mapping_soundspeed = 'none';
            end
            [sound_speed] = medium_pct_soundspeed(parameters, medium, sound_speed, pseudoCT, skull_idx, pct_mapping_soundspeed);
            
            % map skull bone attenuation with the desired algorithm
            if isfield(parameters, "pct_mapping_attenuation")
                pct_mapping_attenuation = parameters.pct_mapping_attenuation;
            else
                pct_mapping_attenuation = 'none';
            end
            [alpha_coeff, alpha_power] = medium_pct_attenuation(parameters, medium, alpha_coeff, alpha_power, pseudoCT, skull_idx, pct_mapping_attenuation);
            
            % [DEBUG] save pCT mapping overview
            if parameters.debug == 1
                h = figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.25, 1]);
                subplot(4,1,1); 
                    hold on; histogram(pseudoCT(skull_idx)); xlabel("pseudo-HU")
                    title(['pseudoCT tissue property ranges: ', parameters.pseudoCT_variant]);
                subplot(4,1,2); 
                    hold on; histogram(density(skull_idx)); xlabel("Density [kg/m3]")
                    % add lines for the fixed parameters
                    xline(medium.skull_trabecular.density, 'r', 'LineWidth', 2);
                    xline(medium.skull_cortical.density, 'r', 'LineWidth', 2);
                subplot(4,1,3); 
                    hold on; histogram(sound_speed(skull_idx)); xlabel("Sound speed [m/s]")
                    xline(medium.skull_trabecular.sound_speed, 'r', 'LineWidth', 2);
                    xline(medium.skull_cortical.sound_speed, 'r', 'LineWidth', 2);
                subplot(4,1,4); 
                    hold on; histogram(alpha_coeff(skull_idx)); xlabel("Attenuation [dB/(cm.MHzy)]")
                    xline(medium.skull_trabecular.alpha_coeff, 'r', 'LineWidth', 2);
                    xline(medium.skull_cortical.alpha_coeff, 'r', 'LineWidth', 2);
                output_plot = fullfile(parameters.debug_dir, ...
                    sprintf('pCT_histograms_%s.png',parameters.pseudoCT_variant));
                exportgraphics(h, output_plot, 'Resolution', 150);
                close(h);
            end

            clear skull_idx

        else
            thermal_conductivity(medium_masks==medium_i) = medium.(label_name).thermal_conductivity;
            specific_heat(medium_masks==medium_i) = medium.(label_name).specific_heat_capacity;
            sound_speed(medium_masks==medium_i) = medium.(label_name).sound_speed; 
            density(medium_masks==medium_i) = medium.(label_name).density;
            alpha_coeff(medium_masks==medium_i) = medium.(label_name).alpha_coeff; 
            alpha_power(medium_masks==medium_i) = medium.(label_name).alpha_power;
            if isfield(medium.(label_name), 'perfusion')
                perfusion(medium_masks==medium_i) = medium.(label_name).perfusion;
            end
            if isfield(medium.(label_name), 'absorption_fraction')
                absorption_fraction(medium_masks==medium_i) = medium.(label_name).absorption_fraction;
            end
            if isfield(parameters.thermal.temp_0, label_name)
                temp_0(medium_masks==medium_i) = parameters.thermal.temp_0.(label_name);
            end
        end
    end

    % convert perfusion rate [mL/min/kg] into perfusion coefficient [1/s]
    perfusion_coeff = (perfusion ./ 60) .* density * 1e-6;

    %% account for k-Wave's actual attenuation behaviour
    % limits discrepancies for high attenuation estimates (see https://doi.org/10.1121/1.4894790).
    % 'alpha_coeff' is rescaled to match the specified alpha_coeff and alpa_power_true for the center frequency

    if ~isfield(parameters, 'fit_alpha_power')
        parameters.fit_alpha_power = 1;
    end

    if parameters.fit_alpha_power == 1
        alpha_power_fixed = 2;
    
        if parameters.debug == 1
            plot_fit = true;
        else
            plot_fit = false;
        end
    
        if numel(unique([parameters.transducer.source_freq_hz])) > 1
            warning('Multiple source frequencies not yet supported. Using %i Hz from first transducer.', ...
                    parameters.transducer(1).source_freq_hz)
        end
        
        alpha_coeff_fixed = fitPowerLawParamsMulti(...
            alpha_coeff, ...
            alpha_power, ...
            sound_speed, ...
            parameters.transducer(1).source_freq_hz, ...
            alpha_power_fixed, ...
            plot_fit);
    
        % DEBUG mode: save plot of fitted attenuation values
        if parameters.debug == 1
            fig_path = fullfile(parameters.debug_dir, ...
            ['attenuation_fit', char(parameters.results_filename_affix), '.png']);
            saveas(gcf, fig_path);
            close(gcf);
        end
    else
        warning('Attenuation power-law is not remapped (as requested).')
        alpha_coeff_fixed = alpha_coeff;
        alpha_power_fixed = alpha_power;
    end

    %% smooth medium masks

    if isfield(parameters, 'smooth_properties') && parameters.smooth_properties == true
        disp("Smoothing acoustic proprty maps ...");

        tmp_density = density; % keep unsmoothed image for figure

        sound_speed = smooth_img(sound_speed, parameters.smooth_window, 0, parameters.smooth_method);
        density = smooth_img(density, parameters.smooth_window, 0, parameters.smooth_method);
        alpha_coeff_fixed = smooth_img(alpha_coeff_fixed, parameters.smooth_window, 0, parameters.smooth_method);
        thermal_conductivity = smooth_img(thermal_conductivity, parameters.smooth_window, 0, parameters.smooth_method);
        specific_heat = smooth_img(specific_heat, parameters.smooth_window, 0, parameters.smooth_method);
        perfusion_coeff = smooth_img(perfusion_coeff, parameters.smooth_window, 0, parameters.smooth_method);
        absorption_fraction = smooth_img(absorption_fraction, parameters.smooth_window, 0, parameters.smooth_method);
        temp_0 = smooth_img(temp_0, parameters.smooth_window, 0, parameters.smooth_method);
    
        % [DEBUG] Plot unsmoothed and smoothed density
        if parameters.debug == 1
            h = figure('Position', [100 100 800 400]);
            if numel(size(density))==3
                density_pre = squeeze(tmp_density(:,round(size(tmp_density,2)/2),:));
                density_post = squeeze(density(:,round(size(density,2)/2),:));
            else
                density_pre = squeeze(tmp_density);
                density_post = squeeze(density);
            end
            subplot(1,2,1); imagesc(density_pre); title('Original density')
            subplot(1,2,2); imagesc(density_post); title('Smoothed density')
            output_plot_filename = fullfile(parameters.debug_dir, ...
                sprintf('sub-%03d_%s_density_smoothing_changes%s.png', ...
                parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
            saveas(h, output_plot_filename, 'png')
            close(h);
            clear density_pre density_post;
        end; clear tmp_density;
    else
        disp('No smoothing applied to acoustic property maps ...')
    end


    %% specify the medium as a kWave-compatible structure 
    % Note: absorption_fraction and temp_0 need to be removed later)
    kwave_medium = struct('sound_speed', sound_speed, ...
                          'density', density, ...
                          'alpha_coeff', alpha_coeff_fixed,...
                          'alpha_power', alpha_power_fixed, ...
                          'thermal_conductivity', thermal_conductivity,...
                          'specific_heat', specific_heat,...
                          'perfusion_coeff', perfusion_coeff,...
                          'absorption_fraction', absorption_fraction,...
                          'temp_0', temp_0);
    
    %% [debug] save acoustic property images
    if parameters.debug == 1
        % save raw medium matrices as niftis
        try
            filename_density = fullfile(parameters.debug_dir, sprintf('matrix_density'));
            niftiwrite(density, filename_density, 'Compressed',true); pause(0.1);
            filename_sound_speed = fullfile(parameters.debug_dir, sprintf('matrix_sound_speed'));
            niftiwrite(sound_speed, filename_sound_speed, 'Compressed',true); pause(0.1);
            filename_alpha_coeff = fullfile(parameters.debug_dir, sprintf('matrix_alpha_coeff'));
            niftiwrite(alpha_coeff, filename_alpha_coeff, 'Compressed',true); pause(0.1);
            filename_alpha_power = fullfile(parameters.debug_dir, sprintf('matrix_alpha_power'));
            niftiwrite(alpha_power, filename_alpha_power, 'Compressed',true); pause(0.1);
            filename_alpha_coeff_fixed = fullfile(parameters.debug_dir, sprintf('matrix_alpha_coeff_fixed'));
            niftiwrite(alpha_coeff_fixed, filename_alpha_coeff_fixed, 'Compressed',true); pause(0.1);
            filename_perfusion = fullfile(parameters.debug_dir, sprintf('matrix_perfusion'));
            niftiwrite(perfusion_coeff, filename_perfusion, 'Compressed',true); pause(0.1);
            filename_absorption = fullfile(parameters.debug_dir, sprintf('matrix_absorption'));
            niftiwrite(absorption_fraction, filename_absorption, 'Compressed',true); pause(0.1);
        catch
            warning("Error with saving debug images: medium mapping. May result from concurrent write attempts...")
        end

        % save images of assigned medium properties with proper headers
        if contains(parameters.simulation_medium, {'layered'}) && (exist('planimg') & ~isempty(planimg))
            medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'sound_speed')
            medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'density')
            medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'alpha_coeff')
            medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'alpha_power')
            medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'thermal_conductivity')
            medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'specific_heat')
            medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'perfusion_coeff')
            medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'absorption_fraction')
        end

        % save a pCT (if used)
        if parameters.use_pseudoCT == 1
            filename_pct = fullfile(parameters.debug_dir, sprintf('pct%s', parameters.results_filename_affix));
            niftiwrite(pseudoCT, filename_pct, 'Compressed',true);
        end
    end

end