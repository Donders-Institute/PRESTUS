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
    thermal_conductivity = empty_grid;       % [W/(m.K)]
    specific_heat = empty_grid;            % [J/(kg.K)]
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
            thermal_conductivity(medium_masks==medium_i) = medium.(label_name).thermal_conductivity;                 % [W/(m.K)]
            specific_heat(medium_masks==medium_i) = medium.(label_name).specific_heat_capacity;                      % [J/(kg.K)]
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
                temp_0(medium_masks==medium_i) = parameters.thermal.temp_0.(label_name);                             % [degC]
            end
        end
    end

    % account for k-Wave's actual attenuation behaviour
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

    % convert perfusion rate [mL/min/kg] into perfusion coefficient [1/s]
    perfusion_coeff = (perfusion ./ 60) .* density * 1e-6; % [1/s]

    % specify the medium as a kWave-compatible structure 
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
    
    % save images for debugging
    if parameters.debug == 1
        try
            filename_density = fullfile(parameters.debug_dir, sprintf('matrix_density'));
            niftiwrite(density, filename_density, 'Compressed',true);
            pause(0.1);
            filename_sound_speed = fullfile(parameters.debug_dir, sprintf('matrix_sound_speed'));
            niftiwrite(sound_speed, filename_sound_speed, 'Compressed',true);
            pause(0.1);
            filename_alpha_coeff = fullfile(parameters.debug_dir, sprintf('matrix_alpha_coeff'));
            niftiwrite(alpha_coeff, filename_alpha_coeff, 'Compressed',true);
            pause(0.1);
            filename_alpha_power = fullfile(parameters.debug_dir, sprintf('matrix_alpha_power'));
            niftiwrite(alpha_power, filename_alpha_power, 'Compressed',true);
            pause(0.1);
            filename_alpha_coeff_fixed = fullfile(parameters.debug_dir, sprintf('matrix_alpha_coeff_fixed'));
            niftiwrite(alpha_coeff_fixed, filename_alpha_coeff_fixed, 'Compressed',true);
            pause(0.1);
            filename_perfusion = fullfile(parameters.debug_dir, sprintf('matrix_perfusion'));
            niftiwrite(perfusion_coeff, filename_perfusion, 'Compressed',true);
            pause(0.1);
            filename_absorption = fullfile(parameters.debug_dir, sprintf('matrix_absorption'));
            niftiwrite(absorption_fraction, filename_absorption, 'Compressed',true);
        catch
            warning("Error with saving debug images: medium mapping. May result from concurrent write attempts...")
        end
    end

    % save images of assigned medium properties
    if contains(parameters.simulation_medium, {'layered'}) && parameters.debug == 1 && (exist('planimg') & ~isempty(planimg))
        medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'sound_speed')
        medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'density')
        medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'alpha_coeff')
        medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'alpha_power')
        medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'thermal_conductivity')
        medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'specific_heat')
        medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'perfusion_coeff')
        medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'absorption_fraction')
    end

    % save a pCT if used
    if parameters.use_pseudoCT == 1 && parameters.debug == 1
        filename_pct = fullfile(parameters.debug_dir, sprintf('pct%s', parameters.results_filename_affix));
        niftiwrite(pseudoCT, filename_pct, 'Compressed',true);
    end

end