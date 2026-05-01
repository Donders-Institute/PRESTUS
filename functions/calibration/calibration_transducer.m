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

    %% Optimization (run once — phases are intensity-independent)

    % Compute analytical O'Neil solution and scaling factor to simulated intensity
    [profile_oneil, simulated_analytical_scaling] = ...
        compute_oneil_solution(...
        initial_params, ...
        profile_sim, ...
        profile_target_ref);

    % Optimize phases [rad] and source amplitude to match real profile shape
    [opt_phases, opt_velocity_ref, min_err, opt_params_raw] = ...
        perform_global_search(...
        initial_params, ...
        profile_target_ref, ...
        profile_sim.velocity);

    % Shared phase outputs (same for all intensities)
    opt_source_phase_rad = opt_phases;
    opt_source_phase_deg = opt_phases / pi * 180;

    %% Determine amplitude calibration strategy
    % 'scale'   (default) — skip validation reruns; save compact precession
    %                       model. Requires opt_phase_precession = 'linear'
    %                       or 'monotonic'.
    % 'compute'           — run validation simulation per intensity (original
    %                       high-fidelity path).
    if isfield(parameters.calibration, 'amp_calibration')
        amp_calibration = parameters.calibration.amp_calibration;
    else
        amp_calibration = 'scale';
    end

    if strcmp(amp_calibration, 'scale')
        precession_mode = '';
        if isfield(parameters.calibration, 'opt_phase_precession')
            precession_mode = parameters.calibration.opt_phase_precession;
        end
        if ~ismember(precession_mode, {'linear', 'monotonic'})
            error(['calibration_transducer: amp_calibration=''scale'' requires ' ...
                'calibration.opt_phase_precession to be ''linear'' or ''monotonic''.']);
        end
    end

    %% Per-intensity loop: velocity fit, (one) validation rerun, save

    opt_source_amp = zeros(1, N_intensities);

    % profile_sim_opt_ref: cached from the first validation rerun.
    % For k > 1, the simulated field scales linearly with source amplitude
    % (k-Wave operates in the linear regime), so I ∝ v² means
    %   profile_sim_opt_k.Isppa = profile_sim_opt_ref.Isppa * (I_k / I_ref)
    % No additional k-Wave runs are needed.
    profile_sim_opt_ref = [];
    intensity_ref       = desired_intensities(1);

    for k = 1:N_intensities
        desired_intensity_k = desired_intensities(k);

        % Scale empirical profile to this intensity
        [profile_target_k, ~] = scale_real_intensity_profile(...
            parameters, profile_empirical, desired_intensity_k);

        % Fit velocity to match desired peak intensity exactly
        if ~isfield(parameters.calibration, 'fit_velocity_to_intensity') || ...
                parameters.calibration.fit_velocity_to_intensity
            [opt_velocity_k, ~, ~] = fit_velocity_to_intensity(...
                initial_params, profile_oneil, opt_phases, opt_velocity_ref, ...
                desired_intensity_k, simulated_analytical_scaling);
        else
            opt_velocity_k = opt_velocity_ref;
        end

        % Calculate optimised source amplitude for this intensity
        opt_source_amp_k = round(opt_velocity_k / profile_sim.velocity * ...
            initial_params.transducer.annular.elem_amp / ...
            simulated_analytical_scaling);
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
            opt_phases, ...
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

        if strcmp(amp_calibration, 'scale')
            %% Parametric mode: no validation rerun
            % Save the compact precession model and amplitude; skip simulation.
            save_parametric_model(opt_param, opt_params_raw, ...
                parameters.calibration.opt_phase_precession);

        else
            %% Full mode: run validation simulation (first intensity only)

            if k == 1
                opt_param.subject_id       = sim_id;
                opt_param.hpc.wait_for_job = true;
                prestus_pipeline_start(opt_param);

                opt_res = load(fullfile(opt_param.io.outputs_folder, 'cache', ...
                    sprintf('sub-%03d_water_results%s.mat', sim_id, opt_param.io.output_affix)));

                opt_params_ref              = opt_res.acoustic_info.parameters;
                opt_params_ref.calibration.prefix = 'Opt_';
                profile_sim_opt_ref         = extract_simulated_profile(opt_res, opt_params_ref);

                profile_sim_opt_k = profile_sim_opt_ref;

            else
                %% Derive scaled profile analytically — no simulation needed
                % Intensity scales as (I_k / I_ref); all profile fields that carry
                % intensity units are scaled accordingly.  Fields with pressure units
                % scale as sqrt(I_k / I_ref).

                scale_I  = desired_intensity_k / intensity_ref;
                scale_p  = sqrt(scale_I);

                profile_sim_opt_k = profile_sim_opt_ref;

                % Intensity fields
                intensity_fields = {'Isppa', 'Ispta', 'Ita'};
                for f = intensity_fields
                    fn = f{1};
                    if isfield(profile_sim_opt_k, fn)
                        profile_sim_opt_k.(fn) = profile_sim_opt_ref.(fn) * scale_I;
                    end
                end

                % Pressure fields
                pressure_fields = {'p_max', 'p_min', 'p_rms'};
                for f = pressure_fields
                    fn = f{1};
                    if isfield(profile_sim_opt_k, fn)
                        profile_sim_opt_k.(fn) = profile_sim_opt_ref.(fn) * scale_p;
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

        end % amp_calibration branch

    end % intensity loop

    %% Clean up cache folder (optional)
    % The free-water simulations write large .mat files into the cache folder.
    % We forced save_acoustic_matrices=1 so the pipeline saves them between
    % runs, but the whole folder can be removed once calibration is complete.
    % Set parameters.io.save_acoustic_matrices = 1 to keep the cache for
    % debugging (e.g. to inspect intermediate simulation results).
    if ~isfield(parameters.io, 'save_acoustic_matrices') || ...
            parameters.io.save_acoustic_matrices == 0
        cache_folder = fullfile(sim_param.io.outputs_folder, 'cache');
        if isfolder(cache_folder)
            rmdir(cache_folder, 's');
            disp(['Calibration cache removed: ', cache_folder]);
        end
    end

end
