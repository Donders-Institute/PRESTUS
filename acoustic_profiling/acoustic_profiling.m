close all; clear;
format long; % to display the intensity values correclty

SOUND_SPEED_WATER = 1500; % [m/s]
DENSITY_WATER = 997; % [kg/m3]

% Prepare running simulations
func_path = fileparts(mfilename('fullpath')); % get the path of the current script
main_folder = fileparts(func_path); % get the main folder path

cd(main_folder) % change directory to the main folder

% Add paths for necessary functions and toolboxes
addpath('functions')
addpath('acoustic_profiling')
addpath(genpath('toolboxes'))

% Load configuration file
equip_param = yaml.loadFile('acoustic_profiling/config.yaml', 'ConvertToArray', true);

% Get the available combinations from the configuration file
available_combos = fieldnames(equip_param.combos);
disp(available_combos)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Path containing 'Axial profiles' and 'PRESTUS virtual parameter' folders 
% of /home/common/matlab/PRESTUS.
data_path = '/home/common/matlab/PRESTUS/';
submit_medium = 'matlab'; % run scripts via 'matlab' (debugging) or via a job using 'slurm' (recommended) or 'qsub'

combinations = {'IS_PCD15287_01001_IGT_128_ch_comb_10_ch'};
%{'CTX_250_001_SC_203_035', 'CTX_250_009_SC_203_035', 'CTX_250_014_SC_203_035', 'CTX_250_026_SC_203_035', 'CTX_500_006_SC_203_035', 'CTX_500_026_SC_203_035',
    %'CTX_250_001_SC_105_010', 'CTX_250_009_SC_105_010', 'CTX_250_014_SC_105_010', 'CTX_250_026_SC_105_010', 'CTX_500_006_SC_105_010', 'CTX_500_026_SC_105_010'};
focal_depths_wrt_exit_plane = {30}; %{[], [], [], [], [], []}; % [mm], every [] will be performed for same index combination
desired_intensities = {[30]}; %{[30, 60], [30, 60], [30, 60], [30, 60], [30, 60], [30, 60]}; % [W/cm2], every [] will be performed for same index combination and for each foci in same index focal depth set []

% Set location for output
sim_param.output_location = '/home/neuromod/marcorn/Documents/acoustic_profiling/output/';
sim_param.data_path = sim_param.output_location;

%% Determine virtual parameters for chosen input parameters

% Set simulation parameters
sim_param.simulation_medium = 'water';
if strcmp(submit_medium, 'matlab') == true
    sim_param.code_type = 'matlab_cpu';
    sim_param.using_donders_hpc = 0;
else
    sim_param.overwrite_files = 'always';
    sim_param.interactive = 0;
end

for i = 1:size(combinations, 2)
    combo_name = combinations{i};
    combo = equip_param.combos.(combo_name);
    
    ds_serial = combo.ds_serial;
    tran_serial = combo.tran_serial;
    
    % Extract transducer information
    tran = equip_param.trans.(tran_serial);
    sim_param.transducer = tran.prestus.transducer;
    sim_param.transducer.source_phase_deg = zeros(1, tran.prestus.transducer.n_elements); % initial phases
    sim_param.transducer.source_amp = repmat([200000], 1, tran.prestus.transducer.n_elements); % initial amplitude

    % Extract phase information
    tran_manufact = equip_param.trans.(tran_serial).manufact;
    combo.phase_table = fullfile(data_path, combo.phase_table);
    if isequal(tran_manufact, "Sonic Concepts")
        phase_table = readtable(combo.phase_table);
        distance = phase_table.Distance;
        phases = phase_table(:, 2:end-1); % remove TPOSetpoint column
    elseif isequal(tran_manufact, "Imasonic")
        % Read transducer ini
        phase_table = readIniFile(combo.phase_table);
    end

    % Combine directory, folder and file name to one path
    combo.char_data_path = fullfile(data_path, combo.char_data_path);

    % Load default parameters including additional chosen parameters like
    % transducer
    sim_param = load_parameters(sim_param);

    % Extract driving system information
    ds = equip_param.ds.(ds_serial);
    fprintf('Equipment: transducer %s and driving system %s \n', [tran.name, ds.name])

    % Extract characterization data
    charac_data = readmatrix(combo.char_data_path);

    % Extract available foci
    available_foci_wrt_exit_plane = round(charac_data(2, 2:end), 2);
        
    % Extract distance wrt exit plane
    dist_from_exit_plane = charac_data(3:end, 1);

    % Calculate distance wrt mid bowl of transducer instead of exit plane
    % Simulation results are wrt mid bowl of transducer
    dist_tran_exit_plane = sim_param.transducer.curv_radius_mm - sim_param.transducer.dist_to_plane_mm;
    dist_from_tran = dist_from_exit_plane + dist_tran_exit_plane;

    % Extract intensity data 
    intens_data = charac_data(3:end, 2:end);

    % Construct the file path
    % Combine directory, folder and file name to one path
    prestus_dir = fullfile(data_path, equip_param.gen.prestus_virt_folder);

    [~, filename, ext] = fileparts(combo.char_data_path);
    equipment_name = erase(filename, equip_param.gen.axial_prof_name);
    prestus_virtual_path = fullfile(prestus_dir, strcat(equip_param.gen.prestus_virt_name, equipment_name, '.csv'));

    % If no focal depths are chosen, perform the simulations for all
    % available focal depths
    if size(focal_depths_wrt_exit_plane{i}, 2) == 0
        focal_depths_wrt_exit_plane{i} = available_foci_wrt_exit_plane;
    end

    % Perform acoustic profiling for every focal depth listed specifically for each equipment combination 
    for j = 1:size(focal_depths_wrt_exit_plane{i}, 2)
        focus_wrt_exit_plane = round(focal_depths_wrt_exit_plane{i}(j), 2);
        fprintf('Focus: %.2f \n', focus_wrt_exit_plane)

        % Extract phase information per focal depth
        if isequal(tran_manufact, "Sonic Concepts")
            % Extract row from table
            foc_index = find(distance == focus_wrt_exit_plane);
            init_phases = table2array(phases(foc_index, :));

        elseif isequal(tran_manufact, "Imasonic")
            % Calculate phases belonging to focal depth  
            init_phases = computePhases(SOUND_SPEED_WATER, tran, focus_wrt_exit_plane, phase_table);
        end

        % Convert real number of elements to virtual number of elements
        if tran.n_elem ~= tran.prestus.transducer.n_elements
            original_pos = linspace(1, tran.n_elem, tran.n_elem);

            % Unwrap the phases
            unwrapped_phases = unwrap(init_phases * pi/180) * 180/pi;  % Convert to radians, unwrap, then back to degrees
            
            virtual_pos = linspace(1, tran.n_elem, tran.prestus.transducer.n_elements);

            virtual_phases = interp1(original_pos, unwrapped_phases, virtual_pos, 'linear'); 

            wrapped_virtual_phases = mod(virtual_phases, 360);

            init_phases = wrapped_virtual_phases;
        end
        
        sim_param.transducer.source_phase_deg = init_phases; % initial phases
        
        % TODO check if expected_focal_distance_mm is wrt exit plane
        sim_param.expected_focal_distance_mm = focus_wrt_exit_plane;

        % Check if focal depth is within the range
        if focus_wrt_exit_plane < tran.min_foc || focus_wrt_exit_plane > tran.max_foc
            warning('Chosen focus of %.2f is not within the range of %.2f and %.2f for transducer %s. Focus will be skipped.', focus_wrt_exit_plane, tran.min_foc, tran.max_foc, tran.name)
            continue
        end

        % Check if exact focus is available as axial profile
        col_index = find(available_foci_wrt_exit_plane == focus_wrt_exit_plane);
        
        % If exact focus is not found, find closest foci
        if isempty(col_index)
            [~, closestIndex] = min(abs(available_foci_wrt_exit_plane - focus_wrt_exit_plane));

            closest_foci_wrt_exit_plane = available_foci_wrt_exit_plane(closestIndex);

            if closest_foci_wrt_exit_plane > focus_wrt_exit_plane
                closestIndex2 = closestIndex;

                % Find second measurement to perform an interpolation
                closestIndex1 = closestIndex2 - 1;
            else
                closestIndex1 = closestIndex;
                
                % Find second measurement to perform an interpolation
                closestIndex2 = closestIndex1 + 1; 
            end

            % Check if second measurement exists
            if closestIndex2 > size(available_foci_wrt_exit_plane, 2)
                error('Focus is higher than last available axial profile. Therefore, interpolation between two profiles is not possible.')
            end

            % Interpolate to get profile at desired focus
            focus_wrt_exit_plane_1 = round(available_foci_wrt_exit_plane(closestIndex1), 2);
            focus_wrt_exit_plane_2 = round(available_foci_wrt_exit_plane(closestIndex2), 2);
            profile_1 = intens_data(:, closestIndex1)';
            profile_2 = intens_data(:, closestIndex2)';
    
            profile_focus = interp1([focus_wrt_exit_plane_1, focus_wrt_exit_plane_2], [profile_1; profile_2], focus_wrt_exit_plane);

            figure;
            plot(dist_from_tran, profile_1, 'DisplayName', ...
                ['Measurement 1, focus at ' num2str(focus_wrt_exit_plane_1)]);
            hold on;
            plot(dist_from_tran, profile_2, 'DisplayName', ...
                ['Measurement 2, focus at ' num2str(focus_wrt_exit_plane_2)]);
            hold on;
            plot(dist_from_tran, profile_focus, 'DisplayName', ...
                ['Interpolated, focus at ' num2str(focus_wrt_exit_plane)]);
            legend
            xlabel('Distance wrt exit plane [mm]');
            ylabel('Intensity [W/cm^2]');
            title('Interpolation input and results')
        else
            profile_focus = intens_data(:, col_index)';

            figure;
            plot(dist_from_tran, profile_focus);
            xlabel('Distance wrt exit plane [mm]');
            ylabel('Intensity [W/cm^2]');
            title(['Axial profile at focus wrt exit plane: ' num2str(focus_wrt_exit_plane) ' [mm]'])
        end
        
        fig_path = fullfile(prestus_dir, strcat('Interpolation_at_F_', num2str(focus_wrt_exit_plane), '_', equipment_name, '.png'));
        saveas(gcf, fig_path);

        % TODO bring 10 mm to front-end
        % Ignore first 10 mm of signal to prevent catching max peak in 
        % near field peak and search for the maximum from that point on.
        [minValue, closestIndex] = min(abs(dist_from_tran-10));

        max_intens = max(profile_focus(closestIndex:end));

        for k = 1:size(desired_intensities{i}, 2)
            % Scale profile according to desired intensity
            desired_intensity = desired_intensities{i}(k);
            fprintf('Current maximum intensity: %.2f \n Desired maximum intensity: %.2f \n Adjust profile to match desired maximum intensity. \n ', [max_intens, desired_intensity])
        
            p_pa = sqrt(2*desired_intensity*10^4*DENSITY_WATER*SOUND_SPEED_WATER);
            sim_param.transducer.source_amp = repmat(p_pa, 1, tran.prestus.transducer.n_elements); % initial amplitude

            adjustment_factor_intensity = max_intens / desired_intensity;
            adjusted_profile_focus = profile_focus' ./ adjustment_factor_intensity;

            % All parameters are calculated to perform simulations
            sim_id = i*j*k;
            disp('Run initial simulation...')
            
            if strcmp(submit_medium, 'qsub') == true
                single_subject_pipeline_with_qsub(sim_id, sim_param, true);
            elseif strcmp(submit_medium, 'slurm') == true
                single_subject_pipeline_with_slurm(sim_id, sim_param, true);
            elseif strcmp(submit_medium, 'matlab') == true
                single_subject_pipeline(sim_id, sim_param);
            else
                error('Submit medium does not correspond to available options.')
            end

            % Load initial results
            outputs_folder = sprintf('%s/sub-%03d', sim_param.data_path, sim_id);
            initial_res = load(sprintf('%s/sub-%03d_water_results%s.mat', outputs_folder, sim_id, sim_param.results_filename_affix),'sensor_data','parameters');
            
            % Get maximum pressure
            p_max = gather(initial_res.sensor_data.p_max_all); % transform from GPU array to normal array
            
            % Plot 2D intensity map
            figure;
            imagesc((1:size(p_max, 1)) * initial_res.parameters.grid_step_mm, ...
                (1:size(p_max, 3)) * initial_res.parameters.grid_step_mm, ...
                squeeze(p_max(:, initial_res.parameters.transducer.pos_grid(2), :))')
            axis image;
            colormap(getColorMap);
            xlabel('Lateral Position [mm]');
            ylabel('Axial Position [mm]');
            axis image;
            cb = colorbar;
            title('Pressure for the focal plane')
            
            fig_path = fullfile(prestus_dir, strcat('Intensity_map_2D_at_F_', num2str(focus_wrt_exit_plane), '_at_I_', num2str(desired_intensity), '_', equipment_name, '.png'));
            saveas(gcf, fig_path);

            % Simulated pressure along the focal axis
            pred_axial_pressure = squeeze(p_max(initial_res.parameters.transducer.pos_grid(1), initial_res.parameters.transducer.pos_grid(2),:)); % get the values at the focal axis
            
            [p_axial_oneil, simulated_grid_adj_factor, velocity, axial_position] = compute_oneil_solution(initial_res.parameters, pred_axial_pressure, dist_from_tran, adjusted_profile_focus, focus_wrt_exit_plane, desired_intensity, prestus_dir, equipment_name);
            
            % Optimize phases and source amplitude to match real profile
            [opt_phases, opt_velocity, min_err] = perform_global_search(initial_res.parameters, dist_from_tran, adjusted_profile_focus, velocity);

            % Recalculate analytical solution based on optimized phases and
            % velocity
            p_axial_oneil_opt = recalculate_analytical_sol(initial_res.parameters, p_axial_oneil, opt_phases, opt_velocity, dist_from_tran, adjusted_profile_focus, axial_position, focus_wrt_exit_plane, desired_intensity, prestus_dir, equipment_name);
            
            % Calculate optimized source amplitude
            opt_source_amp = round(opt_velocity / velocity * initial_res.parameters.transducer.source_amp / simulated_grid_adj_factor);
            sprintf('the optimised source_amp = %i', opt_source_amp(1))

            % Redo simulation with optimized phases and source amplitude
            opt_param = sim_param;
            opt_param.transducer.source_amp = opt_source_amp;
            opt_param.transducer.source_phase_rad = [0 opt_phases];
            opt_param.transducer.source_phase_deg = [0 opt_phases]/pi*180;
            opt_param.results_filename_affix = '_optimized';

            if strcmp(submit_medium, 'qsub') == true
                single_subject_pipeline_with_qsub(sim_id, opt_param, true);
            elseif strcmp(submit_medium, 'slurm') == true
                single_subject_pipeline_with_slurm(sim_id, opt_param, true);
            elseif strcmp(submit_medium, 'matlab') == true
                single_subject_pipeline(sim_id, opt_param);
            else
                error('Submit medium does not correspond to available options.')
            end

            plot_opt_sim_results(opt_param, outputs_folder, sim_id, axial_position, dist_from_tran, adjusted_profile_focus, p_axial_oneil_opt, p_axial_oneil, focus_wrt_exit_plane, desired_intensity, prestus_dir, equipment_name, min_err)

            save_optimized_values(prestus_virtual_path, focus_wrt_exit_plane, desired_intensity, opt_param, ds_serial, tran_serial);
        end

    end


end

%% Function definitions

function iniData = readIniFile(filename)
    % READINIFILE Reads an INI file into a MATLAB struct.
    % 
    % INPUT:
    %   filename - Path to the INI file.
    % OUTPUT:
    %   iniData  - Struct containing the parsed data.

    % Open the file
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open the file: %s', filename);
    end

    iniData = struct();
    currentSection = '';
    
    % Read the file line by line
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        
        % Skip empty lines or comments
        if isempty(line) || startsWith(line, ';') || startsWith(line, '#')
            continue;
        end
        
        % Check if the line is a section header
        if startsWith(line, '[') && endsWith(line, ']')
            currentSection = line(2:end-1); % Remove [ and ]
            iniData.(currentSection) = struct(); % Create a new section
            continue;
        end
        
        % Parse key-value pairs
        if contains(line, '=')
            tokens = strsplit(line, '=');
            key = strtrim(tokens{1});
            value = strtrim(tokens{2});
            
            % Convert numeric values if possible
            numericValue = str2double(value);
            if ~isnan(numericValue)
                value = numericValue;
            end
            
            % Assign the key-value pair to the current section
            % Convert key to a valid field name
            validKey = matlab.lang.makeValidName(key);
            if isempty(currentSection)
                iniData.(validKey) = value;
            else
                iniData.(currentSection).(validKey) = value;
            end
        end
    end
    
    % Close the file
    fclose(fid);
end

function phases = computePhases(SOUND_SPEED_WATER, tran, focus_wrt_exit_plane, tran_ini_data)
    % COMPUTEPHASES Computes the phases necessary to aim at the specified point.
    %   Modifies the 'pulse' object by updating its phases to aim at the target.
    %
    % Inputs:
    %   pulse            - Struct or object with frequency information.
    %   point_mm         - [x, y, z] Cartesian coordinates (in mm) of the target.
    %   set_focus_mm     - Focal depth [mm] relative to the exit plane.
    %   dephasing_degree - Array of dephasing degrees, or empty for no dephasing.
    %
    % Outputs:
    %   success          - Boolean indicating if phase computation succeeded.

    wavelen = SOUND_SPEED_WATER / tran.prestus.transducer.source_freq_hz; % Wavelength for single frequency [m]
    
    % Initialize phases
    phases = zeros(1, tran.n_elem);
    
    % Convert focus wrt exit plane to aim wrt natural focus
    dist_to_mid_bowl = tran.prestus.transducer.curv_radius_mm - tran.prestus.transducer.dist_to_plane_mm;
    focus_wrt_mid_bowl = focus_wrt_exit_plane + dist_to_mid_bowl;
    aim_wrt_natural_focus = tran.prestus.transducer.curv_radius_mm - focus_wrt_mid_bowl; 

    point_mm = [0, 0, aim_wrt_natural_focus];

    % Target point in meters
    x = point_mm(1) / 1000.0;
    y = point_mm(2) / 1000.0;
    z = point_mm(3) / 1000.0;
    
    % Compute phases for each channel
    for i = 1:tran.n_elem
        field_name = sprintf('x%d', i); 

        str_elem = tran_ini_data.elements.(field_name); % Get element position [x, y, z]

        elem_cell = strsplit(str_elem, '|');
        
        % Convert to numeric array
        elem = str2double(elem_cell) / 1000; % convert m to mm

        % Distance from element to target
        dist = sqrt((elem(1) - x)^2 + (elem(2) - y)^2 + (elem(3) - z)^2);
        
        % Fractional part of distance in wavelengths
        rem = mod(dist / wavelen, 1);
        
        % Compute phase in degrees
        phases(i) = rem * 360.0;
    end

end

function [p_axial_oneil, simulated_grid_adj_factor, velocity, axial_position] = compute_oneil_solution(parameters, pred_axial_pressure, dist_exit_plane, adjusted_profile_focus, focus_wrt_exit_plane, desired_intensity, prestus_dir, equipment_name)
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
    plot(dist_exit_plane, adjusted_profile_focus)
    hold off
    xline(parameters.expected_focal_distance_mm, '--');
    legend('O''neil''s analytical solution - Initial', 'Simulated - Initial', 'Desired profile')
    title('Pressure along the beam axis')

    fig_path = fullfile(prestus_dir, strcat('Initial_simulation_at_F_', num2str(focus_wrt_exit_plane), '_at_I_', num2str(desired_intensity), '_', equipment_name, '.png'));
    saveas(gcf, fig_path);
    
    % What is distance to the maximum pressure?
    fprintf('Estimated distance to the point of maximum pressure: %.2f mm\n', axial_position(p_axial_oneil == max(p_axial_oneil)))
    
    % Compute the approximate adjustment from simulated (on a grid) to analytic solution
    simulated_grid_adj_factor = max(pred_axial_pressure(:)) / max(p_axial_oneil(:));
end

function [opt_phases, opt_velocity, min_err] = perform_global_search(parameters, dist_exit_plane, adjusted_profile_focus, velocity)
    % Perform a global search to optimize phases and velocity
    %
    % Arguments:
    % - parameters: Structure containing simulation and transducer parameters.
    % - dist_exit_plane: Distance vector for desired profile focus [mm].
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
        dist_exit_plane, adjusted_profile_focus, [], [], ones(length(dist_exit_plane), 1));
    
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
        dist_exit_plane, adjusted_profile_focus, 1, [], ones(length(dist_exit_plane), 1));
    
    fprintf('Optimal phases: %s deg.; velocity: %.2f; optimization error: %.2f \n', mat2str(round(opt_phases_and_velocity(1:(parameters.transducer.n_elements-1))/pi*180)),...
        opt_phases_and_velocity((parameters.transducer.n_elements)), min_err);

    opt_phases = opt_phases_and_velocity(1:(parameters.transducer.n_elements-1));
    opt_velocity = opt_phases_and_velocity(parameters.transducer.n_elements);
end

function p_axial_oneil_opt = recalculate_analytical_sol(parameters, p_axial_oneil, opt_phases, opt_velocity, dist_from_tran, adjusted_profile_focus, axial_position, focus_wrt_exit_plane, desired_intensity, prestus_dir, equipment_name)
    % Recalculate analytical solution based on optimized phases and velocity
    %
    % Arguments:
    % - parameters: Structure containing simulation and transducer parameters.
    % - p_axial_oneil: Initial O'Neil solution for pressure [Pa].
    % - opt_phases: Optimized phases for each transducer element [rad].
    % - opt_velocity: Optimized particle velocity [m/s].
    % - dist_exit_plane: Distance vector for desired profile focus [mm].
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
    
    focus_wrt_mid_bowl = focus_wrt_exit_plane + (parameters.transducer.curv_radius_mm - parameters.transducer.dist_to_plane_mm);

    figure('Position', [10 10 900 500]);
    plot(axial_position, i_axial_oneil);
    xlabel('Distance wrt mid bowl of transducer [mm]');
    ylabel('Intensity [W/cm^2]');
    hold on
    plot(axial_position, i_axial_oneil_opt);
    plot(dist_from_tran, adjusted_profile_focus)
    hold off
    xline(focus_wrt_mid_bowl, '--');
    yline(desired_intensity, '--');
    legend('Original analytical profile', sprintf('Optimized analytical profile to match the desired profile'),'Desired profile')
    title('Pressure along the beam axis')
    
    fig_path = fullfile(prestus_dir, strcat('Recalculated_oneil_at_F_', num2str(focus_wrt_exit_plane), '_at_I_', num2str(desired_intensity), '_', equipment_name, '.png'));
    saveas(gcf, fig_path);

    fprintf('Estimated distance to the point of maximum pressure: %.2f mm\n',axial_position(p_axial_oneil_opt==max(p_axial_oneil_opt)))
    
    fprintf('Estimated distance to the center of half-maximum range: %.2f mm\n', get_flhm_center_position(axial_position, p_axial_oneil_opt))
end

function plot_opt_sim_results(opt_param, outputs_folder, sim_id, axial_position, dist_exit_plane, adjusted_profile_focus, p_axial_oneil_opt, p_axial_oneil, focus_wrt_exit_plane, desired_intensity, prestus_dir, equipment_name, min_err)
    % Plot optimized simulation results and compare with desired profiles
    %
    % Arguments:
    % - opt_param: Structure containing optimized parameters.
    % - outputs_folder: Directory containing output simulation results.
    % - sim_id: Simulation ID for loading specific results.
    % - axial_position: Axial position vector [mm].
    % - dist_exit_plane: Distance vector for desired profile focus [mm].
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
    
    fig_path = fullfile(prestus_dir, strcat('Opt_intensity_map_2D_at_F_', num2str(focus_wrt_exit_plane), '_at_I_', num2str(desired_intensity), '_', equipment_name, '.png'));
    saveas(gcf, fig_path);

    % Simulated pressure along the focal axis
    pred_axial_pressure_opt = squeeze(p_max(opt_res.parameters.transducer.pos_grid(1), opt_res.parameters.transducer.pos_grid(2),:));
  
    focus_wrt_mid_bowl = focus_wrt_exit_plane + (parameters.transducer.curv_radius_mm - parameters.transducer.dist_to_plane_mm);

    % Compare optimized profile
    figure('Position', [10 10 900 500]);
    hold on
    plot(axial_position, p_axial_oneil.^2/(2*opt_param.medium.water.sound_speed*opt_param.medium.water.density) .* 1e-4);
    xlabel('Distance wrt mid bowl of transducer [mm]');
    ylabel('Intensity [W/cm^2]');
    plot(axial_position, p_axial_oneil_opt .^2/(2*opt_param.medium.water.sound_speed*opt_param.medium.water.density) .* 1e-4);

    sim_res_axial_position = axial_position-(opt_res.parameters.transducer.pos_grid(3)-1)*0.5; % axial position for the simulated results, relative to transducer position
    plot(sim_res_axial_position, ...
        pred_axial_pressure_opt .^2/(2*opt_param.medium.water.sound_speed*opt_param.medium.water.density) .* 1e-4);
    plot(dist_exit_plane, adjusted_profile_focus)
    hold off
    xline(focus_wrt_mid_bowl, '--');
    yline(desired_intensity, '--');
    legend('Original simulation', sprintf('Optimized for %2.f mm distance, analytical', opt_res.parameters.expected_focal_distance_mm), ...
        sprintf('Optimized for %2.f mm distance, simulated', opt_res.parameters.expected_focal_distance_mm),'Desired profile','Location', 'best')
    title(strcat('Desired vs optimized profiles - optimization error', {' '}, num2str(min_err)))
    
    fig_path = fullfile(prestus_dir, strcat('Opt_simulation_at_F_', num2str(focus_wrt_exit_plane), '_at_I_', num2str(desired_intensity), '_', equipment_name, '.png'));
    saveas(gcf, fig_path);

    fprintf('Estimated distance to the point of maximum pressure: %.2f mm\n', sim_res_axial_position(pred_axial_pressure_opt==max(pred_axial_pressure_opt)))
end

function save_optimized_values(prestus_path, focus_wrt_exit_plane, desired_intensity, opt_param, ds_serial, tran_serial)
    % Save optimized phases and amplitude values to a CSV file
    %
    % Arguments:
    % - prestus_path: Path to the PRESTUS file for saving optimized values.
    % - focus: Target focal distance [mm].
    % - desired_intensity: Target intensity for optimization [W/cm^2].

    disp('Save optimized values in CSV file...')
        
    fprintf('CSV file can be found here: %s \n', prestus_path)

    % Extract phases
    opt_phases = opt_param.transducer.source_phase_deg;
    opt_amp = opt_param.transducer.source_phase_deg;
    
    % Round down phases
    opt_phases = round(opt_phases, 2);
    source_amp = double(opt_amp(1));
    
    % Check if file exists
    if isfile(prestus_path)
        % Read existing data 
        virtual_data = readcell(prestus_path);
        prestus_foci_wrt_exit_plane = round(cell2mat(virtual_data(1, 2:end)), 2);
        prestus_int = cell2mat(virtual_data(2:end, 1));

        % Find or add focus
        col_index_foc = find(prestus_foci_wrt_exit_plane == focus_wrt_exit_plane);
        if isempty(col_index_foc)
            col_index_foc = size(virtual_data, 2) + 1;
            virtual_data{1, col_index_foc} = focus_wrt_exit_plane;
        else
            col_index_foc = col_index_foc + 1; % adjust for header row
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
        virtual_data{row_index_int, col_index_foc} = mat2str([opt_phases, source_amp]);
        
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
        writecell(virtual_data_sort, prestus_path, 'FileType', 'text');
    else
        % Create a new file
        virtual_data = {'Desired intensity [W/cm^2] - Focus wrt exit plane [mm]', focus_wrt_exit_plane; desired_intensity, mat2str([opt_phases, source_amp])};

        writecell(virtual_data, prestus_path, 'FileType', 'text');
    end

    % Save transducer object with optimized values in yaml file for easy
    % integration in PRESTUS config file
    yaml_file = strcat(tran_serial, '-', ds_serial, '-F', num2str(focus_wrt_exit_plane), 'mm-I', num2str(intensity), 'wpercm2.yaml');
    yaml_path = fullfile(opt_param.output_location, yaml_file);

    opt_param.transducer.set_focus_wrt_exit_plane_mm = focus_wrt_exit_plane;
    opt_param.transducer.set_intensity_w_per_cm2 = intensity;

    % Wrap it in a parent structure
    data = struct('transducer', opt_param.transducer);

    yaml.dumpFile(yaml_path, data);
    fprintf('Transducer including optimized values is saved in a yaml file for easy integration into PRESTUS config file. It can be found here: %s \n', yaml_path);
end