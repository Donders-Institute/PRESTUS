close all; clear;

% Prepare running simulations
func_path = fileparts(mfilename('fullpath')); % get the path of the current script
main_folder = fileparts(func_path); % get the main folder path

cd(main_folder) % change directory to the main folder

% Add paths for necessary functions and toolboxes
addpath('functions')
addpath(genpath('toolboxes')) 
addpath('/home/common/matlab/fieldtrip/qsub') % uncomment if you are using Donders HPC

% Check for GPU availability or the presence of the HPC cluster
%if gpuDeviceCount==0 && ~exist('/home/common/matlab/fieldtrip/qsub','dir')
%    error('This script assumes that you have a GPU available for computations or that you''re using the Donders HPC cluster.')
%end

% Load configuration file
parameters = yaml.loadFile('configs/config.yaml');

% Get the available combinations from the configuration file
available_combos = fieldnames(parameters.combos);
disp(available_combos)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

combinations = {'CTX_250_001_SC_105_010'};
focal_depths = {[24, 45], [34, 39]}; % [mm]
desired_intensities = {[20, 40], [21, 41]}; % [W/cm2], every [] will be performed for same index combination and for each foci in same index focal depth set []

% Set location for output
sim_param.data_path = fullfile(main_folder, 'acoustic_profiling_data');

%% Determine virtual parameters for chosen input parameters

% Set simulation parameters
sim_param.simulation_medium = 'water';
sim_param.overwrite_files = 'always';
sim_param.interactive = 0;

for i = 1:size(combinations)
    combo_name = combinations{i};
    combo = parameters.combos.(combo_name);

    ds_serial = combo.ds_serial;
    tran_serial = combo.tran_serial;

    % Extract transducer information
    tran = parameters.trans.(tran_serial);
    sim_param.transducer = tran.prestus.transducer;
    
    % Load default parameters including additional chosen parameters like
    % transducer
    sim_param = load_parameters(sim_param);

    % Extract driving system information
    ds = parameters.ds.(ds_serial);
    fprintf('Equipment: transducer %s and driving system %s \n', [tran.name, ds.name])

    % Extract characterization data
    charac_data = readmatrix(combo.char_data_path);
    available_foci = round(charac_data(1, 2:end), 2);
    dist = charac_data(2:end, 1);

    % Perform acoustic profiling for every focal depth listed specifically for each equipment combination 
    for j = 1:size(focal_depths{i}, 2)
        focus = round(focal_depths{i}(j), 2);
        fprintf('Focus: %.2f \n', focus)

        % Check if focal depth is within the range
        if focus < tran.min_foc || focus > tran.max_foc
            warning('Chosen focus of %.2f is not within the range of %.2f and %.2f for transducer %s. Focus will be skipped.', focus, tran.min_foc, tran.max_foc, tran.name)
            continue
        end

        % Check if exact focus is available as axial profile
        col_index = find(available_foci == focus);
        
        % If exact focus is not found, find closest foci
        if isempty(col_index)
            [~, closestIndex1] = min(abs(available_foci - focus));

            % Shift index so values can be found in charac_data
            closestIndex1 = closestIndex1 + 1;
            
            % Find second measurement to perform an interpolation
            closestIndex2 = closestIndex1 + 1; 

            % Check if second measurement exists
            if closestIndex2 > size(available_foci, 2)
                closestIndex2 = closestIndex1;
                closestIndex1 = closestIndex1 - 1;

                warning('Focus is higher than last available axial profile. Therefore, interpolation between two profiles is not possible. Axial profiles of foci %.2f and %.2f are used instead. \n', available_foci(closestIndex1), available_foci(closestIndex2))
            end

            % Interpolate to get profile at desired focus
            focus_1 = round(charac_data(1, closestIndex1), 2);
            focus_2 = round(charac_data(1, closestIndex2), 2);
            profile_1 = charac_data(2:end, closestIndex1)';
            profile_2 = charac_data(2:end, closestIndex2)';
    
            profile_focus = interp1([focus_1, focus_2], [profile_1; profile_2], focus);

            figure;
            plot(dist, profile_1, 'DisplayName', ...
                ['Measurement 1, focus at ' num2str(focus_1)]);
            hold on;
            plot(dist, profile_2, 'DisplayName', ...
                ['Measurement 2, focus at ' num2str(focus_2)]);
            hold on;
            plot(dist, profile_focus, 'DisplayName', ...
                ['Interpolated, focus at ' num2str(focus)]);
            legend
            xlabel('Distance wrt exit plane [mm]');
            ylabel('Intensity [W/cm^2]');
            title('Interpolation input and results')
        else
            profile_focus = charac_data(2:end, col_index);

            figure;
            plot(dist, profile_focus);
            xlabel('Distance wrt exit plane [mm]');
            ylabel('Intensity [W/cm^2]');
            title(['Axial profile at focus ' num2str(focus)])
        end
        
        max_intens = max(profile_focus);

        for k = 1:size(desired_intensities{i}, 2)
            % Scale profile according to desired intensity
            desired_intensity = desired_intensities{i}(k);
            fprintf('Current maximum intensity: %.2f \n Desired maximum intensity: %.2f \n Adjust profile to match desired maximum intensity. \n ', [max_intens, desired_intensity])

            adjustment_factor_intensity = max_intens / desired_intensity;
            adjusted_profile_focus = profile_focus ./ adjustment_factor_intensity;

            % All parameters are calculated to perform simulations
            sim_id = i*j*k;
            disp('Run initial simulation...')

            single_subject_pipeline_with_qsub(sim_id, parameters);

            % Load initial results
            outputs_folder = sprintf('%s/sim_outputs/sub-%03d', parameters.data_path, sim_id);
            load(sprintf('%s/sub-%03d_water_results%s.mat', outputs_folder, sim_id, parameters.results_filename_affix),'sensor_data','parameters');
            
            % Get maximum pressure
            p_max = gather(sensor_data.p_max_all); % transform from GPU array to normal array
            
            % Plot 2D intensity map
            imagesc((1:size(p_max, 1)) * parameters.grid_step_mm, ...
                (1:size(p_max, 3)) * parameters.grid_step_mm, ...
                squeeze(p_max(:, parameters.transducer.pos_grid(2), :))')
            axis image;
            colormap(getColorMap);
            xlabel('Lateral Position [mm]');
            ylabel('Axial Position [mm]');
            axis image;
            cb = colorbar;
            title('Pressure for the focal plane')

            % Simulated pressure along the focal axis
            pred_axial_pressure = squeeze(p_max(parameters.transducer.pos_grid(1), parameters.transducer.pos_grid(2),:)); % get the values at the focal axis
            
            [p_axial_oneil, simulated_grid_adj_factor, velocity, axial_position] = compute_oneil_solution(parameters, pred_axial_pressure);
            
            % Optimize phases and source amplitude to match real profile
            [opt_phases, opt_velocity] = perform_global_search(parameters, dist, adjusted_profile_focus, velocity);

            % Recalculate analytical solution based on optimized phases and
            % velocity
            p_axial_oneil_opt = recalculate_analytical_sol(parameters, p_axial_oneil, opt_phases, opt_velocity, dist, adjusted_profile_focus, axial_position);
            
            % Calculate optimized source amplitude
            opt_source_amp = round(opt_velocity / velocity * parameters.transducer.source_amp / simulated_grid_adj_factor);
            sprintf('the optimised source_amp = %i', opt_source_amp(1))

            % Redo simulation with optimized phases and source amplitude
            opt_param = sim_param;
            opt_param.transducer.source_amp = opt_source_amp;
            opt_param.transducer.source_phase_rad = [0 opt_phases];
            opt_param.transducer.source_phase_deg = [0 opt_phases]/pi*180;
            opt_param.results_filename_affix = '_optimized';

            single_subject_pipeline_with_qsub(sim_id, opt_param);

            plot_opt_sim_results(opt_param, outputs_folder, sim_id, axial_position, dist, adjusted_profile_focus, desired_intensity, p_axial_oneil_opt)

            save_optimized_values(parameters, char_data_path, focus, desired_intensity, opt_phases, opt_source_amp);
        end

    end


end

%% Function definitions

function [p_axial_oneil, simulated_grid_adj_factor, velocity, axial_position] = compute_oneil_solution(parameters, pred_axial_pressure)
    % Compute O'Neil solution and plot it along with comparisons
    %
    % Arguments:
    % - parameters: Structure containing simulation and transducer parameters.
    % - pred_axial_pressure: Predicted pressure along the beam axis [Pa].
    %
    % Returns:
    % - p_axial_oneil: Computed O'Neil solution for pressure along the beam axis [Pa].
    % - simulated_grid_adj_factor: Adjustment factor to align simulated pressure with analytical solution.
    % - velocity: Particle velocity [m/s].
    % - axial_position: Axial position vector [mm].

    % Define transducer parameters
    velocity = parameters.transducer.source_amp(1) / (parameters.medium.water.density*parameters.medium.water.sound_speed); % [m/s]
    
    % Define position vectors
    % TODO: should the resolution be retrieved from somewhere else?
    axial_position   = (1:parameters.default_grid_dims(3))*0.5; % [mm]
    
    % Evaluate pressure analytically
    % focusedAnnulusONeil provides an analytic solution for the pressure at the
    % focal (beam) axis
    [p_axial_oneil] = focusedAnnulusONeil(parameters.transducer.curv_radius_mm / 1e3, ...
        [parameters.transducer.Elements_ID_mm; parameters.transducer.Elements_OD_mm] / 1e3, repmat(velocity, 1, parameters.transducer.n_elements), ...
        parameters.transducer.source_phase_rad, parameters.transducer.source_freq_hz, parameters.medium.water.sound_speed, ...
        parameters.medium.water.density, (axial_position - 0.5) * 1e-3);
    
    % Convert pressure to intensity
    i_axial_oneil = p_axial_oneil .^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4;
    pred_axial_intensity = pred_axial_pressure.^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4;

    % Plot focal axis pressure
    figure('Position', [10 10 900 500]);
    plot(axial_position, i_axial_oneil);
    xlabel('Distance wrt exit plane [mm]');
    ylabel('Intensity [W/cm^2]');
    hold on
    plot(axial_position-(parameters.transducer.pos_grid(3)-1)*0.5, pred_axial_intensity,'--');
    plot(dist, adjusted_profile_focus)
    hold off
    xline(parameters.expected_focal_distance_mm, '--');
    legend('O''neil''s analytical solution - Initial', 'Simulated - Initial', 'Desired profile')
    title('Pressure along the beam axis')
    
    % What is distance to the maximum pressure?
    fprintf('Estimated distance to the point of maximum pressure: %.2f mm\n', axial_position(p_axial_oneil == max(p_axial_oneil)))
    
    % Compute the approximate adjustment from simulated (on a grid) to analytic solution
    simulated_grid_adj_factor = max(pred_axial_pressure) / max(p_axial_oneil);
end

function [opt_phases, opt_velocity] = perform_global_search(parameters, dist, adjusted_profile_focus, velocity)
    % Perform a global search to optimize phases and velocity
    %
    % Arguments:
    % - parameters: Structure containing simulation and transducer parameters.
    % - dist: Distance vector for desired profile focus [mm].
    % - adjusted_profile_focus: Adjusted desired intensity profile [W/cm^2].
    % - velocity: Initial velocity estimate [m/s].
    %
    % Returns:
    % - opt_phases: Optimized phases for each transducer element [rad].
    % - opt_velocity: Optimized particle velocity [m/s].

    % Initialize global search
    gs = GlobalSearch;
    %opt_velocity = desired_pressure/max_pressure*velocity;
    
    % Define optimization objective function
    %optimize_phases = @(phases) phase_optimization_annulus(phases, parameters, velocity, axial_position, parameters.expected_focal_distance_mm);
    optimize_phases = @(phases_and_velocity) phase_optimization_annulus_full_curve(phases_and_velocity(1:(parameters.transducer.n_elements-1)),...
        parameters, phases_and_velocity(parameters.transducer.n_elements),...
        dist, adjusted_profile_focus);
    
    % Set random seed for consistency
    rng(195,'twister');

    % Create optimization problem
    problem = createOptimProblem('fmincon','x0', [randi(360, [1 parameters.transducer.n_elements-1])/180*pi velocity],...
        'objective',optimize_phases,'lb',zeros(1,parameters.transducer.n_elements),'ub',[2*pi*ones(1,parameters.transducer.n_elements-1) 0.2],...
        'options', optimoptions('fmincon','OptimalityTolerance', 1e-8)); 
    
    % Run global search
    [opt_phases_and_velocity, min_err] = run(gs,problem);

    % Plot optimization results
    phase_optimization_annulus_full_curve(opt_phases_and_velocity(1:(parameters.transducer.n_elements-1)), parameters,...
        opt_phases_and_velocity(parameters.transducer.n_elements),...
        dist, adjusted_profile_focus, 1);
    
    fprintf('Optimal phases: %s deg.; velocity: %.2f; optimization error: %.2f \n', mat2str(round(opt_phases_and_velocity(1:(parameters.transducer.n_elements-1))/pi*180)),...
        opt_phases_and_velocity((parameters.transducer.n_elements)), min_err);

    opt_phases = opt_phases_and_velocity(1:(parameters.transducer.n_elements-1));
    opt_velocity = opt_phases_and_velocity(parameters.transducer.n_elements);
end

function p_axial_oneil_opt = recalculate_analytical_sol(parameters, p_axial_oneil, opt_phases, opt_velocity, dist, adjusted_profile_focus, axial_position)
    % Recalculate analytical solution based on optimized phases and velocity
    %
    % Arguments:
    % - parameters: Structure containing simulation and transducer parameters.
    % - p_axial_oneil: Initial O'Neil solution for pressure [Pa].
    % - opt_phases: Optimized phases for each transducer element [rad].
    % - opt_velocity: Optimized particle velocity [m/s].
    % - dist: Distance vector for desired profile focus [mm].
    % - adjusted_profile_focus: Adjusted desired intensity profile [W/cm^2].
    % - axial_position: Axial position vector [mm].
    %
    % Returns:
    % - p_axial_oneil_opt: Optimized O'Neil solution for pressure along the beam axis [Pa].

    [p_axial_oneil_opt] = focusedAnnulusONeil(parameters.transducer.curv_radius_mm / 1e3, ...
        [parameters.transducer.Elements_ID_mm; parameters.transducer.Elements_OD_mm] / 1e3, repmat(opt_velocity, 1, parameters.transducer.n_elements), ...
        [0 opt_phases], parameters.transducer.source_freq_hz, parameters.medium.water.sound_speed, ...
        parameters.medium.water.density, (axial_position - 0.5) * 1e-3);
    
    % Convert pressure to intensity
    i_axial_oneil = p_axial_oneil .^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4;
    i_axial_oneil_opt = p_axial_oneil_opt .^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4;

    figure('Position', [10 10 900 500]);
    plot(axial_position, i_axial_oneil);
    xlabel('Distance wrt exit plane [mm]');
    ylabel('Intensity [W/cm^2]');
    hold on
    plot(axial_position, i_axial_oneil_opt);
    plot(dist, adjusted_profile_focus)
    hold off
    xline(parameters.expected_focal_distance_mm, '--');
    yline(30, '--');
    legend('Original analytical profile', sprintf('Optimized analytical profile to match the desired profile'),'Desired profile')
    title('Pressure along the beam axis')
    
    fprintf('Estimated distance to the point of maximum pressure: %.2f mm\n',axial_position(p_axial_oneil_opt==max(p_axial_oneil_opt)))
    
    fprintf('Estimated distance to the center of half-maximum range: %.2f mm\n', get_flhm_center_position(axial_position, p_axial_oneil_opt))
end

function plot_opt_sim_results(opt_param, outputs_folder, sim_id, axial_position, dist, adjusted_profile_focus, desired_intensity, p_axial_oneil_opt)
    % Plot optimized simulation results and compare with desired profiles
    %
    % Arguments:
    % - opt_param: Structure containing optimized parameters.
    % - outputs_folder: Directory containing output simulation results.
    % - sim_id: Simulation ID for loading specific results.
    % - axial_position: Axial position vector [mm].
    % - dist: Distance vector for desired profile focus [mm].
    % - adjusted_profile_focus: Adjusted desired intensity profile [W/cm^2].
    % - desired_intensity: Target intensity for optimization [W/cm^2].
    % - p_axial_oneil_opt: Optimized O'Neil solution for pressure [Pa].  
    
    % Load optimized simulation results    
    opt_res = load(sprintf('%s/sub-%03d_water_results%s.mat', outputs_folder, sim_id, opt_param.results_filename_affix),'sensor_data','parameters');

    % Get maximum pressure
    p_max = gather(opt_res.sensor_data.p_max_all);

    % Plot 2D intensity map
    imagesc((1:size(p_max,1)) * opt_res.parameters.grid_step_mm, ...
        (1:size(p_max,3)) * opt_res.parameters.grid_step_mm , ...
        squeeze(p_max(:, opt_res.parameters.transducer.pos_grid(2), :))')
    axis image;
    colormap(getColorMap);
    xlabel('Lateral Position [mm]');
    ylabel('Axial Position [mm]');
    axis image;
    colorbar;
    title('Pressure for the focal plane')

    % Simulated pressure along the focal axis
    pred_axial_pressure_opt = squeeze(p_max(opt_res.parameters.transducer.pos_grid(1), opt_res.parameters.transducer.pos_grid(2),:));
  
    % Compare optimized profile
    figure('Position', [10 10 900 500]);
    hold on
    plot(axial_position, p_axial_oneil.^2/(2*opt_param.medium.water.sound_speed*opt_param.medium.water.density) .* 1e-4);
    xlabel('Distance wrt exit plane [mm]');
    ylabel('Intensity [W/cm^2]');
    plot(axial_position, p_axial_oneil_opt .^2/(2*opt_param.medium.water.sound_speed*opt_param.medium.water.density) .* 1e-4);
    
    sim_res_axial_position = axial_position-(opt_res.parameters.transducer.pos_grid(3)-1)*0.5; % axial position for the simulated results, relative to transducer position
    plot(sim_res_axial_position, ...
        pred_axial_pressure_opt .^2/(2*opt_param.medium.water.sound_speed*opt_param.medium.water.density) .* 1e-4);
    plot(dist, adjusted_profile_focus)
    hold off
    xline(opt_res.parameters.expected_focal_distance_mm, '--');
    yline(desired_intensity, '--');
    legend('Original simulation', sprintf('Optimized for %2.f mm distance, analytical', opt_res.parameters.expected_focal_distance_mm), ...
        sprintf('Optimized for %2.f mm distance, simulated', opt_res.parameters.expected_focal_distance_mm),'Desired profile','Location', 'best')
    title('Desired vs optimized profiles')

    fprintf('Estimated distance to the point of maximum pressure: %.2f mm\n', sim_res_axial_position(pred_axial_pressure_opt==max(pred_axial_pressure_opt)))
end

function save_optimized_values(parameters, char_data_path, focus, desired_intensity, opt_phases, opt_source_amp)
    % Save optimized phases and amplitude values to an Excel file
    %
    % Arguments:
    % - parameters: Structure containing simulation and transducer parameters.
    % - char_data_path: Path to the Excel file for saving optimized values.
    % - focus: Target focal distance [mm].
    % - desired_intensity: Target intensity for optimization [W/cm^2].
    % - opt_phases: Optimized phases for each transducer element [rad].
    % - opt_source_amp: Optimized source amplitude [Pa].

    disp('Save optimized values in Excel file...')

    % Construct the file path
    prestus_dir = parameters.gen.prestus_virt_path;
    [~, filename, ext] = fileparts(char_data_path);
    equipment_name = erase(filename, parameters.gen.axial_prof_name);
    prestus_path = fullfile(prestus_dir, strcat(parameters.gen.prestus_virt_name, equipment_name, ext));
        
    fprintf('Excel file can be found here: %s \n', prestus_path)
    
    % Check if file exists
    if isfile(prestus_path)
        % Read existing data 
        virtual_data = readcell(prestus_path);
        prestus_foci = round(cell2mat(virtual_data(1, 2:end)), 2);
        prestus_int = cell2mat(virtual_data(2:end, 1));

        % Find or add focus
        col_index_foc = find(prestus_foci == focus);
        if isempty(col_index_foc)
            col_index_foc = size(virtual_data, 2) + 1;
            virtual_data{1, col_index_foc} = focus;
        end

        % Find or add desired intensity
        row_index_int = find(prestus_int == desired_intensity);
        
        if isempty(row_index_int)
            row_index_int = size(virtual_data, 1) + 1;
            virtual_data{row_index_int, 1} = desired_intensity;
        else
            row_index_int = row_index_int + 1; % adjust for header row
        end

        % Save optimized values, old values will be overwritten
        virtual_data{row_index_int, col_index_foc} = mat2str([opt_phases, opt_source_amp]);
        
        % Handle missing values
        mask = cellfun(@(x) any(isa(x,'missing')),virtual_data);
        virtual_data(mask) = {[]};

        % Sort rows by desired intensity
        first_col = [virtual_data{2:end, 1}];
        [~, sortIdxCol] = sort(first_col);
        virtual_data_sort = [virtual_data(1, :); virtual_data(1 + sortIdxCol, :)];

        % Sort columns by focus
        first_row = [virtual_data_sort{1, 2:end}];
        [~, sortIdxRow] = sort(first_row);
        virtual_data_sort = [virtual_data_sort(:, 1), virtual_data_sort(:, 1 + sortIdxRow)];

        % Write sorted data back to file
        writecell(virtual_data_sort, prestus_path)
    else
        % Create a new file
        virtual_data = {'Desired intensity [W/cm^2]', focus; desired_intensity, mat2str([opt_phases, opt_source_amp])};

        writecell(virtual_data, prestus_path);
    end
end