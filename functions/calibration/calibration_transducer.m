function [opt_source_amp, opt_source_phase_deg, opt_source_phase_rad] = calibration_transducer(...
    profile_empirical,...
    equipment_name, ...
    desired_intensity, ...
    desired_focal_distance_ep, ...
    parameters, ...
    sim_id)

% CALIBRATION_TRANSDUCER  Calibrate transducer phases and amplitude to match a target intensity profile
%
% Scales an empirical intensity profile to the desired intensity, runs an
% initial water simulation, optimises element phases and amplitude to match
% the target, reruns with optimised parameters, and saves the results.
%
% desired_intensity may be a scalar or a row vector of target intensities.
% When a vector is provided the initial simulation and phase optimisation run
% once (using the first intensity), and the intensity-dependent steps
% (velocity fit, validation rerun, save) loop over each entry — avoiding
% redundant expensive simulations.
%
% Use as:
%   [opt_source_amp, opt_source_phase_deg, opt_source_phase_rad] = ...
%       calibration_transducer(profile_empirical, equipment_name, ...
%                              desired_intensity, desired_focal_distance_ep, parameters, sim_id)
%
% Input:
%   profile_empirical         - struct with axial_intensity and axial_distance_bowl [mm]
%   equipment_name            - name/identifier for the transducer or equipment
%   desired_intensity         - target focal intensity/intensities [W/cm²] (scalar or row vector)
%   desired_focal_distance_ep - desired focal distance from transducer exit plane [mm]
%   parameters                - PRESTUS config with calibration settings
%   sim_id                    - numeric subject or simulation ID
%
% Output:
%   opt_source_amp        - optimised source amplitude(s), one per desired intensity
%   opt_source_phase_deg  - optimised element phases [°]  (shared across intensities)
%   opt_source_phase_rad  - optimised element phases [rad] (shared across intensities)
%
% See also: SCALE_REAL_INTENSITY_PROFILE, COMPUTE_PHASES, PERFORM_GLOBAL_SEARCH,
%           SAVE_OPTIMIZED_VALUES

    desired_intensities = desired_intensity(:)';   % ensure row vector
    N_intensities       = numel(desired_intensities);

    %% Attach more information to parameters

    parameters.calibration.equipment_name         = equipment_name;
    parameters.calibration.desired_focal_distance_ep = desired_focal_distance_ep;

    % Use the first intensity as the reference for the initial simulation
    % and phase optimisation. Since I ∝ v², optimal phases are the same for
    % every intensity level — only velocity (and hence amplitude) changes.
    parameters.calibration.desired_intensity = desired_intensities(1);

    % Scale the requested profile to the reference intensity
    [profile_target_ref, ~] = scale_real_intensity_profile(...
        parameters, ...
        profile_empirical, ...
        desired_intensities(1));

    %% Initial water simulation

    disp('Run initial simulation...')

    % Run all simulations in calibration folder (default)
    if parameters.calibration.save_in_calibration_folder
        parameters.path.anat = parameters.calibration.path_output;
        parameters.path.seg  = parameters.calibration.path_output;
        parameters.path.sim  = parameters.calibration.path_output;
    end
    disp(['Saving free-water calibration in ', parameters.path.sim]);

    parameters.io.outputs_folder = fullfile(parameters.path.sim, sprintf('sub-%03d', sim_id));

    % Copy calibration settings to relevant entries in simulation config
    sim_param = parameters;
    % Force water medium
    sim_param.simulation.medium = 'water';
    % Force save result matrices. save_acoustic_matrices is set explicitly
    % because config_default.yaml sets it to 0, and per-field flags take
    % precedence over the global save_matrices fallback in should_save_output.
    % Calibration requires the acoustic sensor_data .mat between runs.
    sim_param.io.save_matrices = 1;
    sim_param.io.save_acoustic_matrices = 1;
    % Never run nifti creation for calibration 
    sim_param.modules.run_nifti_creation = 0;
    % Never run thermal simulation for calibration 
    sim_param.modules.run_heating_sims = 0;
    % Overwrite transducer kwavearray modeling (if specified)
    if isfield(parameters.calibration, 'force_kwavearray') && ...
            parameters.calibration.force_kwavearray == 1
        sim_param.grid.use_kWaveArray = 1; % force to run with kwavearray setup
    end
    % Convert from default 3D to 2D axisymmetric simulation (if requested)
    if isfield(parameters.calibration, 'axisymmetric2D') && ...
            parameters.calibration.axisymmetric2D == 1
        parameters.grid.axisymmetric = 1;
        if numel(parameters.grid.default_dims)==3
            parameters.grid.default_dims(2) = [];
        end
    end
    % Force deactivate interactive mode
    sim_param.simulation.interactive = 0;

    % Run the simulation based on the submission method
    sim_param.subject_id = sim_id;
    sim_param.hpc.wait_for_job = true;
    prestus_pipeline_start(sim_param);

    %% Load initial results

    initial_res = load(fullfile(sim_param.io.outputs_folder, 'cache', ...
        sprintf('sub-%03d_water_results%s.mat', sim_id, sim_param.io.output_affix)));

    initial_params = initial_res.acoustic_info.parameters;
    initial_params.calibration.prefix = 'Initial_';

    %% Extract simulated intensity along the focal axis

    [profile_sim] = extract_simulated_profile(initial_res, initial_params);

    clear initial_res;

    %% Phase determination and O'Neil solution (run once — phases are intensity-independent)
    %
    % Two modes, selected by parameters.calibration.elem_phase_correction_deg:
    %
    %   Geometric correction (non-empty):
    %     Phases = geometric phases at this depth (already in initial_params)
    %              + pre-calibrated per-element hardware correction.
    %     Skips perform_global_search; only velocity is fit per intensity.
    %     Enable by setting parameters.calibration.elem_correction_file in
    %     calibration_standalone.m (loaded by calibration_setup).
    %
    %   Global search (default, empty correction):
    %     Existing behaviour: perform_global_search optimises phases and
    %     velocity jointly against the empirical profile.
    %     To also save the correction for future use at other depths, set
    %     parameters.calibration.save_elem_correction = true.

    use_geo_correction = isfield(parameters.calibration, 'elem_phase_correction_deg') && ...
        ~isempty(parameters.calibration.elem_phase_correction_deg);

    if use_geo_correction
        % --- Geometric steering + pre-calibrated hardware correction ----------
        correction_deg = parameters.calibration.elem_phase_correction_deg(:)';

        if numel(correction_deg) ~= parameters.transducer.annular.elem_n
            error(['calibration_transducer: elem_phase_correction_deg has %d elements ' ...
                'but transducer has %d elements.'], ...
                numel(correction_deg), parameters.transducer.annular.elem_n);
        end

        geo_phases_deg = initial_params.transducer.annular.elem_phase_deg;
        combined_deg   = mod(geo_phases_deg + correction_deg, 360);
        combined_rad   = combined_deg * pi / 180;

        % Update a parameter copy so O'Neil is evaluated with the corrected phases
        params_corrected = initial_params;
        params_corrected.transducer.annular.elem_phase_deg = combined_deg;
        params_corrected.transducer.annular.elem_phase_rad = combined_rad;

        [profile_oneil, simulated_analytical_scaling] = compute_oneil_solution(...
            params_corrected, profile_sim, profile_target_ref);

        opt_source_phase_rad = combined_rad;
        opt_source_phase_deg = combined_deg;
        opt_velocity_ref     = profile_sim.velocity;
        min_err              = NaN;

        fprintf('Geometric correction mode: geo + correction phases [deg]: %s\n', ...
            mat2str(round(combined_deg)));

    else
        % --- Existing global search path --------------------------------------
        [profile_oneil, simulated_analytical_scaling] = compute_oneil_solution(...
            initial_params, profile_sim, profile_target_ref);

        [opt_source_phase_rad, opt_velocity_ref, min_err, ~] = perform_global_search(...
            initial_params, profile_target_ref, profile_sim.velocity);

        opt_source_phase_deg = opt_source_phase_rad / pi * 180;
    end

    %% Save per-element correction (only during a reference-depth global-search run)
    if ~use_geo_correction && ...
            isfield(parameters.calibration, 'save_elem_correction') && ...
            parameters.calibration.save_elem_correction
        save_elem_correction(parameters.calibration.equipment_yaml_path, opt_source_phase_deg, ...
            initial_params.transducer.annular.elem_phase_deg, desired_focal_distance_ep);
    end

    %% Intensity calibration loop

    opt_source_amp = zeros(1, N_intensities);

    for k = 1:N_intensities
        desired_intensity_k = desired_intensities(k);

        % Scale empirical profile to target intensity
        [profile_target_k, ~] = scale_real_intensity_profile(...
            parameters, profile_empirical, desired_intensity_k);

        % Scale velocity to match ISPPA
        [opt_velocity_k, ~, ~] = fit_velocity_to_intensity(...
            initial_params, profile_oneil, opt_source_phase_rad, opt_velocity_ref, ...
            desired_intensity_k);

        % Scale optimised source amplitude to match ISPPA
        velocity_scaling = opt_velocity_k / profile_sim.velocity;
        opt_source_amp_k = round(velocity_scaling * ...
            (initial_params.transducer.annular.elem_amp / ...
            sqrt(simulated_analytical_scaling)));
        opt_source_amp(k) = opt_source_amp_k(1);

        fprintf('The optimized elem_amp = %i (intensity %g W/cm²)\n', ...
            opt_source_amp_k(1), desired_intensity_k);

        % Update parameters for this intensity
        params_k = initial_params;
        params_k.calibration.desired_intensity = desired_intensity_k;

        % Recalculate analytical solution with optimized phases and velocity
        profile_oneil_opt_k = recompute_oneil_solution(...
            params_k, ...
            profile_oneil, ...
            profile_target_k, ...
            opt_source_phase_rad, ...
            opt_velocity_k);

        % Use a per-intensity output affix when calibrating multiple intensities
        if N_intensities > 1
            intensity_affix = sprintf('_optimized_I%g', desired_intensity_k);
        else
            intensity_affix = '_optimized';
        end

        opt_param = sim_param;
        opt_param.transducer.annular.elem_amp       = opt_source_amp_k;
        opt_param.transducer.annular.elem_phase_rad = opt_source_phase_rad;
        opt_param.transducer.annular.elem_phase_deg = opt_source_phase_deg;
        opt_param.io.output_affix                   = intensity_affix;
        opt_param.calibration.desired_intensity     = desired_intensity_k;

        % Choose whether to run a free-water simulation to validate
        % intensity or rescale profile analytically
        if isfield(opt_param.calibration, 'opt_amp_validation')
            opt_amp_validation = opt_param.calibration.opt_amp_validation;
        else
            opt_amp_validation = 'always';
        end

        if strcmp(opt_amp_validation, 'always') || ...
                (strcmp(opt_amp_validation, 'initial') && k == 1) || ...
                (strcmp(opt_amp_validation, 'final') && k == N_intensities)

            % Run free-water validation simulation

            opt_param.subject_id       = sim_id;
            opt_param.hpc.wait_for_job = true;
            prestus_pipeline_start(opt_param);

            opt_res = load(fullfile(opt_param.io.outputs_folder, 'cache', ...
                sprintf('sub-%03d_water_results%s.mat', sim_id, opt_param.io.output_affix)));

            opt_params_ref              = opt_res.acoustic_info.parameters;
            opt_params_ref.calibration.prefix = 'Opt_';
            profile_sim_opt_k           = extract_simulated_profile(opt_res, opt_params_ref);

            clear opt_res;

        else

            % Create dummy simulation fields (copy structure from profile_sim if
            % profile_sim_opt_k does not yet exist, e.g. when validation is skipped
            % on the very first iteration)
            if ~exist('profile_sim_opt_k', 'var')
                profile_sim_opt_k = profile_sim;
            end

            % Intensity fields
            intensity_fields = {'Isppa', 'Ispta', 'Ita'};
            for f = intensity_fields
                fn = f{1};
                if isfield(profile_sim_opt_k, fn)
                    profile_sim_opt_k.(fn) = NaN;
                end
            end

            % Pressure fields
            pressure_fields = {'p_max', 'p_min', 'p_rms'};
            for f = pressure_fields
                fn = f{1};
                if isfield(profile_sim_opt_k, fn)
                    profile_sim_opt_k.(fn) = NaN;
                end
            end
        end

        % Plot optimized simulation results
        plot_opt_sim_results(...
            opt_param, ...
            profile_target_k, ...
            profile_oneil, ...
            profile_oneil_opt_k, ...
            profile_sim, ...
            profile_sim_opt_k, ...
            min_err)

        % Save optimized values (full CSV + YAML)
        save_optimized_values(opt_param);

    end % intensity loop

    %% Clean up cache folder (optional)
    % The free-water simulations write large .mat files into the cache folder.
    % We forced save_acoustic_matrices=1 so the pipeline saves them between
    % runs, but the whole folder can be removed once calibration is complete.
    % Set parameters.io.save_acoustic_matrices = 1 to keep the cache for
    % debugging (e.g. to inspect intermediate simulation results).
    if ~isfield(parameters.io, 'save_acoustic_matrices') || ...
            parameters.io.save_acoustic_matrices == 0 || ...
            ~isfield(parameters.io, 'save_matrices') || ...
            parameters.io.save_matrices == 0
        cache_folder = fullfile(sim_param.io.outputs_folder, 'cache');
        if isfolder(cache_folder)
            rmdir(cache_folder, 's');
            disp(['Calibration cache removed: ', cache_folder]);
        end
    end

end
