%% PREPROCESSING - Main MATLAB script for acoustic profiling simulations
% This script configures and executes acoustic profiling simulations 
% using equipment configurations and user-defined input parameters. 
% By performing focal adjustments, intensity scaling, and interpolating
% acoustic profiles, it determines the real axial profile. This profile
% is then used to optimize the virtual transducer's amplitude and phases
% to replicate the behavior of the actual transducer measured in a water tank.
% Due to model constraints, additional virtual elements might be required 
% to mimic the actual transducer behavior accurately. The script
% calc_virtual_elem.m in acoustic_profiling/functions can be used for this purpose.
% 
% This preprocessing step generates optimized amplitude and phase values, 
% which serve as inputs for subsequent acoustic and/or heating simulations.

% Clear workspace and close all figures
close all; clear;
format long; % Ensures accurate display of intensity values

%% Initialization
% Constants
SOUND_SPEED_WATER = 1500; % Speed of sound in water [m/s]
DENSITY_WATER = 997; % Density of water [kg/m^3]

% Ignore first part of the measured axial profile to not find the maximum
% in the near-field peak
skip_front_peak_mm = 10; % [mm]

% Set up paths
% Determine the current and main folder paths
func_path = fileparts(mfilename('fullpath'));
main_folder = fileparts(func_path);
cd(main_folder); % Change directory to the main folder

% Add necessary paths for functions and toolboxes
addpath('functions');
addpath(genpath('acoustic_profiling'));
addpath(genpath('configs'));
addpath(genpath('toolboxes'));

% Load configuration
% Load equipment parameters from the YAML configuration file
equip_param = yaml.loadFile('equipment_config.yaml', 'ConvertToArray', true);

% Display available equipment combinations
available_combos = fieldnames(equip_param.combos);
disp('Available Equipment Combinations:');
disp(available_combos);

%% User-defined input parameters
% Define data path and simulation medium
% Path containing 'Axial profiles' and 'PRESTUS virtual parameter' folders 
% of /home/common/matlab/PRESTUS.
data_path = '\\ru.nl\WrkGrp\FUS_Researchers\';
submit_medium = 'matlab'; % Options: 'matlab' (debugging), 'slurm' (recommended), 'qsub'

% Specify equipment combinations and simulation parameters
combinations = {'IS_PCD15287_01001_IGT_128_ch_comb_10_ch', 'IS_PCD15287_01002_IGT_128_ch_comb_10_ch'};

% Every [] will be performed for same index equipment combination. Example:
% all foci within the first [] will be executed for the first equipment 
% combination. If array is empty ([]), focal depths of available 
% characterization data will be used.
focal_depths_wrt_exit_plane = {[27], [30, 50]}; % Focal depths in mm

% Every [] will be performed for same index equipment combination and for 
% each foci in same index focal depth set []. Example: for the first
% equipment combination, all foci within the first [] will be executed per
% intensity within the first desired intensity set [].
desired_intensities = {[30, 60], [60]}; % Desired intensities in W/cm^2

% Output location for simulation results
sim_param.output_location = 'C:\Temp\PRESTUS\';
sim_param.data_path = sim_param.output_location;

% Save data in general PRESTUS output folder
save_in_general_folder = false;

%% Set up simulation parameters
sim_param.simulation_medium = 'water';
if strcmp(submit_medium, 'matlab') == true
    sim_param.code_type = 'matlab_cpu';
    sim_param.using_donders_hpc = 0;
else
    sim_param.overwrite_files = 'always';
    sim_param.interactive = 0;
end

%% Iterate through equipment combinations
for i = 1:length(combinations)
    combo_name = combinations{i};
    combo = equip_param.combos.(combo_name);
    
    % Extract equipment details
    ds_serial = combo.ds_serial;
    ds = equip_param.ds.(ds_serial);
    tran_serial = combo.tran_serial;
    tran = equip_param.trans.(tran_serial);
    
    % Configure transducer parameters
    sim_param.transducer = tran.prestus.transducer;
    sim_param.transducer.source_phase_deg = zeros(1, tran.prestus.transducer.n_elements); % initial phases
    sim_param.transducer.source_amp = repmat(200000, 1, tran.prestus.transducer.n_elements); % initial amplitude

    % Construct output directory
    prestus_dir = fullfile(data_path, equip_param.gen.prestus_virt_folder);
    [~, filename, ext] = fileparts(combo.char_data_path);
    equipment_name = erase(filename, equip_param.gen.axial_prof_name);
    prestus_virtual_path = fullfile(prestus_dir, strcat(equip_param.gen.prestus_virt_name, equipment_name, '.csv'));

    % Load phase data
    combo.phase_table = fullfile(data_path, combo.phase_table);
    if isequal(tran.manufact, "Sonic Concepts")
        phase_table = readtable(combo.phase_table);
    elseif isequal(tran.manufact, "Imasonic")
        phase_table = read_ini_file(combo.phase_table);
    else
        error('Unsupported transducer manufacturer: %s', tran.manufact);
    end

    % Load characterization data
    combo.char_data_path = fullfile(data_path, combo.char_data_path);
    charac_data = readmatrix(combo.char_data_path);
    available_foci_wrt_exit_plane = round(charac_data(2, 2:end), 2);
    dist_from_exit_plane = charac_data(3:end, 1);
    intens_data = charac_data(3:end, 2:end);

    % Convert exit plane reference to mid-bowl reference
    % Simulation results are w.r.t. mid-bowl of transducer
    dist_tran_exit_plane = sim_param.transducer.curv_radius_mm - sim_param.transducer.dist_to_plane_mm;
    dist_from_tran = dist_from_exit_plane + dist_tran_exit_plane;

    % Ensure focal depths are specified. If no focal depths, perform the 
    % simulations for all available focal depths
    if isempty(focal_depths_wrt_exit_plane{i})
        focal_depths_wrt_exit_plane{i} = available_foci_wrt_exit_plane;
    end

    % Load default parameters including additional chosen parameters like
    % transducer
    sim_param = load_parameters(sim_param);

    fprintf('Equipment: transducer %s and driving system %s \n', [tran.name, ds.name])

    %% Iterate through focal depths
    for j = 1:length(focal_depths_wrt_exit_plane{i})
        focus_wrt_exit_plane = round(focal_depths_wrt_exit_plane{i}(j), 2);
        fprintf('Focus: %.2f \n', focus_wrt_exit_plane)

        % Verify focal range
        if focus_wrt_exit_plane < tran.min_foc || focus_wrt_exit_plane > tran.max_foc
            warning('Focus %.2f mm is outside the range [%.2f, %.2f] mm. Skipping.\n', ...
                focus_wrt_exit_plane, tran.min_foc, tran.max_foc);
            continue;
        end

        % TODO double-check if expected_focal_distance_mm is wrt exit plane
        sim_param.expected_focal_distance_mm = focus_wrt_exit_plane;

        sim_param = set_real_phases(phase_table, tran, focus_wrt_exit_plane, sim_param, SOUND_SPEED_WATER);

        % Interpolate or select axial profile
        [profile_focus, max_intens] = extract_real_intensity_profile(available_foci_wrt_exit_plane, focus_wrt_exit_plane, intens_data, prestus_dir, equipment_name, dist_from_tran, skip_front_peak_mm, sim_param.output_location, save_in_general_folder);

        for k = 1:length(desired_intensities{i})
            desired_intensity = desired_intensities{i}(k);

            adjusted_profile_focus = scale_real_intensity_profile(sim_param, desired_intensity, DENSITY_WATER, SOUND_SPEED_WATER, profile_focus);

            %% Run initial simulation
            sim_id = i*j*k;
            disp('Run initial simulation...')
            
            % Run the simulation based on the submission method
            switch submit_medium
                case 'qsub'
                    single_subject_pipeline_with_qsub(sim_id, sim_param, true);
                case 'slurm'
                    single_subject_pipeline_with_slurm(sim_id, sim_param, true);
                case 'matlab'
                    single_subject_pipeline(sim_id, sim_param);
                otherwise
                    error('Submit medium does not correspond to available options.');
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
            
            % Save the intensity map
            if save_in_general_folder
                fig_path = fullfile(prestus_dir, strcat('Intensity_map_2D_at_F_', num2str(focus_wrt_exit_plane), '_at_I_', num2str(desired_intensity), '_', equipment_name, '.png'));
            else
                fig_path = fullfile(sim_param.output_location, strcat('Intensity_map_2D_at_F_', num2str(focus_wrt_exit_plane), '_at_I_', num2str(desired_intensity), '_', equipment_name, '.png'));     
            end

            saveas(gcf, fig_path);

            % Simulated pressure along the focal axis
            pred_axial_pressure = squeeze(p_max(initial_res.parameters.transducer.pos_grid(1), initial_res.parameters.transducer.pos_grid(2),:)); % get the values at the focal axis
            
            % Compute O'Neil solution and related parameters
            [p_axial_oneil, simulated_grid_adj_factor, velocity, axial_position] = compute_oneil_solution(initial_res.parameters, pred_axial_pressure, dist_from_tran, adjusted_profile_focus, focus_wrt_exit_plane, desired_intensity, prestus_dir, equipment_name, save_in_general_folder);
            
            % Optimize phases and source amplitude to match real profile
            [opt_phases, opt_velocity, min_err] = perform_global_search(initial_res.parameters, dist_from_tran, adjusted_profile_focus, velocity);

            % Recalculate analytical solution based on optimized phases and
            % velocity
            p_axial_oneil_opt = recalculate_analytical_sol(initial_res.parameters, p_axial_oneil, opt_phases, opt_velocity, dist_from_tran, adjusted_profile_focus, axial_position, focus_wrt_exit_plane, desired_intensity, prestus_dir, equipment_name, save_in_general_folder);
            
            % Calculate optimized source amplitude
            opt_source_amp = round(opt_velocity / velocity * initial_res.parameters.transducer.source_amp / simulated_grid_adj_factor);
            fprintf('The optimized source_amp = %i\n', opt_source_amp(1));

            % Redo simulation with optimized phases and source amplitude
            opt_param = sim_param;
            opt_param.transducer.source_amp = opt_source_amp;
            opt_param.transducer.source_phase_rad = [0 opt_phases];
            opt_param.transducer.source_phase_deg = [0 opt_phases]/pi*180;
            opt_param.results_filename_affix = '_optimized';

            % Run the simulation again with optimized parameters
            switch submit_medium
                case 'qsub'
                    single_subject_pipeline_with_qsub(sim_id, opt_param, true);
                case 'slurm'
                    single_subject_pipeline_with_slurm(sim_id, opt_param, true);
                case 'matlab'
                    single_subject_pipeline(sim_id, opt_param);
                otherwise
                    error('Submit medium does not correspond to available options.');
            end

            % Plot optimized simulation results
            plot_opt_sim_results(opt_param, outputs_folder, sim_id, axial_position, dist_from_tran, adjusted_profile_focus, p_axial_oneil_opt, p_axial_oneil, focus_wrt_exit_plane, desired_intensity, prestus_dir, equipment_name, min_err, save_in_general_folder)

            % Save optimized values
            save_optimized_values(prestus_virtual_path, focus_wrt_exit_plane, desired_intensity, opt_param, ds_serial, tran_serial, save_in_general_folder);
        end
    end
end