function [opt_source_amp, opt_source_phase_deg, opt_source_phase_rad] = calibration_transducer(...
    profile_empirical,...
    equipment_name, ...
    desired_intensity, ...
    desired_focal_distance_ep, ...
    parameters, ...
    sim_id)

% CALIBRATION_TRANSDUCER  Calibrate transducer phases and amplitude to match a target intensity profile
%
% Determines optimal element phases and source amplitude to reproduce a
% measured empirical intensity profile at one or more target intensities.
%
% Pipeline overview:
%   1. Phase optimisation runs purely analytically using the configured
%      forward model (O'Neil or Rayleigh; see calibration.forward_model).
%   2. An optional free-water k-Wave simulation (calibration.run_free_water_sim,
%      default true) then corrects the analytical-to-simulation scaling factor.
%      Set to false for a fully simulation-free calibration — appropriate when
%      the analytical model alone is considered sufficient (e.g. single-element
%      transducers with no inter-element phase steering).
%   3. For each target intensity the velocity and source amplitude are scaled
%      analytically (I ∝ v²), followed by an optional per-intensity validation
%      simulation (calibration.opt_amp_validation).
%
% Because optimal phases are independent of intensity level, phase optimisation
% runs once and the per-intensity loop only adjusts velocity and amplitude.
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

    % =========================================================================
    % STAGE 1 — Prepare target profile
    %
    % Attach run-specific metadata to the parameter struct and scale the
    % empirical hydrophone measurement to the first (reference) target intensity.
    % This scaled profile is the shape the transducer must reproduce; all
    % subsequent optimisation steps work to match it.
    % The first intensity is used as the reference because optimal phases are
    % independent of level (I ∝ v²) — only velocity changes across intensities.
    % =========================================================================

    parameters.calibration.equipment_name            = equipment_name;
    parameters.calibration.desired_focal_distance_ep = desired_focal_distance_ep;
    parameters.calibration.desired_intensity         = desired_intensities(1);

    [profile_target_ref, ~] = scale_real_intensity_profile(...
        parameters, ...
        profile_empirical, ...
        desired_intensities(1));

    % =========================================================================
    % STAGE 2 — Analytical bootstrap
    %
    % Compute an axial intensity profile from the configured analytical forward
    % model (O'Neil or Rayleigh; see calibration.forward_model) using the
    % nominal transducer parameters from the config.  This gives us a velocity
    % estimate and an intensity shape that phase optimisation can work with —
    % no k-Wave simulation is needed at this stage.
    %
    % If run_free_water_sim = true (default, Stage 4), this profile is later
    % replaced with the actual k-Wave result.
    % =========================================================================

    bootstrap_params = parameters;
    bootstrap_params.calibration.prefix = 'Initial_';
    if ~isfield(bootstrap_params.transducer.annular, 'elem_phase_rad')
        bootstrap_params.transducer.annular.elem_phase_rad = ...
            bootstrap_params.transducer.annular.elem_phase_deg * pi / 180;
    end

    profile_sim = extract_analytical_profile(bootstrap_params);

    % =========================================================================
    % STAGE 3 — Phase determination
    %
    % Find the per-element phases that steer the focus to the desired depth.
    % Two modes are available, selected by calibration.elem_phase_correction_deg:
    %
    %   Geometric correction (field is non-empty):
    %     Applies pre-measured per-element hardware correction offsets on top of
    %     the depth-dependent geometric steering phases already in the config.
    %     No search is run; the offsets encode systematic hardware delays
    %     (e.g. cable length mismatches) measured at a reference depth.
    %
    %   Global search (default, field absent or empty):
    %     perform_global_search jointly optimises phases and velocity to
    %     minimise the difference between the analytical forward model and the
    %     scaled empirical profile.  Only the analytical velocity and target
    %     profile are needed — no simulation data.
    %     Set calibration.save_elem_correction = true to persist the resulting
    %     per-element offsets for reuse at other focal depths.
    % =========================================================================

    % Configure output folder before any simulation or file I/O
    if parameters.calibration.save_in_calibration_folder
        parameters.path.anat = parameters.calibration.path_output;
        parameters.path.seg  = parameters.calibration.path_output;
        parameters.path.sim  = parameters.calibration.path_output;
    end
    disp(['Saving free-water calibration in ', parameters.path.sim]);
    parameters.io.outputs_folder = fullfile(parameters.path.sim, sprintf('sub-%03d', sim_id));

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

        geo_phases_deg       = bootstrap_params.transducer.annular.elem_phase_deg;
        combined_deg         = mod(geo_phases_deg + correction_deg, 360);
        combined_rad         = combined_deg * pi / 180;
        opt_source_phase_rad = combined_rad;
        opt_source_phase_deg = combined_deg;
        opt_velocity_ref     = profile_sim.velocity;
        min_err              = NaN;

        fprintf('Geometric correction mode: geo + correction phases [deg]: %s\n', ...
            mat2str(round(combined_deg)));

    else
        % --- Global search ----------------------------------------------------
        [opt_source_phase_rad, opt_velocity_ref, min_err, ~] = perform_global_search(...
            bootstrap_params, profile_target_ref, profile_sim.velocity);

        opt_source_phase_deg = opt_source_phase_rad / pi * 180;
    end

    % =========================================================================
    % STAGE 4 — Free-water correction simulation (optional)
    %
    % Runs a k-Wave simulation in free water using the nominal (pre-optimised)
    % source parameters.  The simulated peak intensity is compared to the
    % analytical prediction to compute simulated_analytical_scaling — the ratio
    % by which the analytical model under- or over-estimates the simulation.
    % This factor corrects the source amplitude in Stage 5.
    %
    % Set calibration.run_free_water_sim = false to skip this stage.  The
    % scaling factor then defaults to 1 (pure analytical amplitude prediction),
    % which is appropriate for single-element transducers where the analytical
    % model closely matches the simulation.
    % =========================================================================

    sim_param = parameters;
    sim_param.simulation.medium = 'water';
    % save_acoustic_matrices must be set explicitly: the per-field flag takes
    % precedence over the global save_matrices fallback and the default is 0.
    % The correction and validation simulations both need the .mat result files.
    sim_param.io.save_matrices        = 1;
    sim_param.io.save_acoustic_matrices = 1;
    sim_param.modules.run_nifti_creation = 0;
    sim_param.modules.run_heating_sims   = 0;
    if isfield(parameters.calibration, 'force_kwavearray') && ...
            parameters.calibration.force_kwavearray == 1
        sim_param.grid.use_kWaveArray = 1;
    end
    if isfield(parameters.calibration, 'axisymmetric2D') && ...
            parameters.calibration.axisymmetric2D == 1
        sim_param.grid.axisymmetric = 1;
        if numel(sim_param.grid.default_dims) == 3
            sim_param.grid.default_dims(2) = [];
        end
    end
    sim_param.simulation.interactive = 0;
    sim_param.subject_id             = sim_id;
    sim_param.hpc.wait_for_job       = true;

    run_free_water_sim = ~isfield(parameters.calibration, 'run_free_water_sim') || ...
        parameters.calibration.run_free_water_sim;

    if run_free_water_sim

        disp('Run free-water correction simulation...')
        prestus_pipeline_start(sim_param);

        initial_res = load(fullfile(sim_param.io.outputs_folder, 'cache', ...
            sprintf('sub-%03d_water_results%s.mat', sim_id, sim_param.io.output_affix)));

        % resolved_params replaces bootstrap_params: the simulation resolves
        % geometry-dependent fields (e.g. steering phases) that were not
        % fully determined from the config alone.
        resolved_params = initial_res.acoustic_info.parameters;
        resolved_params.calibration.prefix = 'Initial_';
        profile_sim = extract_simulated_profile(initial_res, resolved_params);
        clear initial_res;

        % In geometric correction mode the simulation also computes the true
        % steering phases for the requested focal depth; update combined phases.
        if use_geo_correction
            geo_phases_deg       = resolved_params.transducer.annular.elem_phase_deg;
            combined_deg         = mod(geo_phases_deg + correction_deg, 360);
            combined_rad         = combined_deg * pi / 180;
            opt_source_phase_rad = combined_rad;
            opt_source_phase_deg = combined_deg;
        end

    else
        % No simulation run; resolved_params is identical to bootstrap_params.
        resolved_params = bootstrap_params;
    end

    % =========================================================================
    % STAGE 5 — Compute analytical reference profile with optimised phases
    %
    % Evaluates the forward model once with the final optimised phases to get
    % a clean analytical intensity profile (profile_analytical) and the
    % simulated_analytical_scaling factor.  Both are used in the per-intensity
    % amplitude calculation below.
    % =========================================================================

    params_for_analytical = resolved_params;
    params_for_analytical.transducer.annular.elem_phase_deg = opt_source_phase_deg;
    params_for_analytical.transducer.annular.elem_phase_rad = opt_source_phase_rad;
    [profile_analytical, simulated_analytical_scaling] = compute_analytical_solution(...
        params_for_analytical, profile_sim, profile_target_ref);

    % Optionally save the per-element phase offsets found by global search so
    % they can be reused as hardware corrections at other focal depths.
    if ~use_geo_correction && ...
            isfield(parameters.calibration, 'save_elem_correction') && ...
            parameters.calibration.save_elem_correction
        save_elem_correction(parameters.calibration.equipment_yaml_path, opt_source_phase_deg, ...
            resolved_params.transducer.annular.elem_phase_deg, desired_focal_distance_ep);
    end

    % =========================================================================
    % STAGE 6 — Per-intensity amplitude calibration
    %
    % Phases are shared across all intensities (determined once above).
    % For each target intensity:
    %   a) Scale the empirical profile to the target level.
    %   b) Correct the particle velocity analytically so the forward-model peak
    %      matches the target (I ∝ v², so v_new = v_ref * sqrt(I_target / I_ref)).
    %   c) Derive the source amplitude from the corrected velocity, normalised
    %      by simulated_analytical_scaling to account for any model–sim offset.
    %   d) Optionally run a free-water validation simulation to confirm the
    %      achieved intensity (calibration.opt_amp_validation: 'always' |
    %      'initial' | 'final' | 'none').
    %   e) Plot and save results.
    % =========================================================================

    if isfield(parameters.calibration, 'opt_amp_validation')
        opt_amp_validation = parameters.calibration.opt_amp_validation;
    else
        opt_amp_validation = 'always';
    end

    opt_source_amp = zeros(1, N_intensities);

    for k = 1:N_intensities
        desired_intensity_k = desired_intensities(k);

        % Scale empirical profile to this target intensity
        [profile_target_k, ~] = scale_real_intensity_profile(...
            parameters, profile_empirical, desired_intensity_k);

        % Correct velocity so the analytical peak matches the target intensity
        [opt_velocity_k, ~, ~] = fit_velocity_to_intensity(...
            resolved_params, profile_analytical, opt_source_phase_rad, opt_velocity_ref, ...
            desired_intensity_k);

        % Derive source amplitude: scale from reference velocity, then apply
        % the simulation correction factor (sqrt because amplitude ∝ v ∝ sqrt(I))
        velocity_scaling = opt_velocity_k / profile_sim.velocity;
        opt_source_amp_k = round(velocity_scaling * ...
            (resolved_params.transducer.annular.elem_amp / ...
            sqrt(simulated_analytical_scaling)));
        opt_source_amp(k) = opt_source_amp_k(1);

        fprintf('The optimized elem_amp = %i (intensity %g W/cm²)\n', ...
            opt_source_amp_k(1), desired_intensity_k);

        % Build per-intensity parameter struct for plotting and saving
        params_k = resolved_params;
        params_k.calibration.desired_intensity = desired_intensity_k;

        % Evaluate forward model with the optimised phases and corrected velocity
        profile_analytical_opt_k = recompute_analytical_solution(...
            params_k, ...
            profile_analytical, ...
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

        % --- Optional validation simulation -----------------------------------
        % Runs a free-water k-Wave simulation with the optimised parameters to
        % verify that the achieved intensity matches the target.  Controlled by
        % opt_amp_validation: 'always' runs for every intensity; 'initial' and
        % 'final' run only for the first or last; 'none' skips all.
        if strcmp(opt_amp_validation, 'always') || ...
                (strcmp(opt_amp_validation, 'initial') && k == 1) || ...
                (strcmp(opt_amp_validation, 'final') && k == N_intensities)

            prestus_pipeline_start(opt_param);

            opt_res = load(fullfile(opt_param.io.outputs_folder, 'cache', ...
                sprintf('sub-%03d_water_results%s.mat', sim_id, opt_param.io.output_affix)));

            opt_params_ref                    = opt_res.acoustic_info.parameters;
            opt_params_ref.calibration.prefix = 'Opt_';
            profile_sim_opt_k                 = extract_simulated_profile(opt_res, opt_params_ref);

            clear opt_res;

        else
            % No validation sim; carry forward the previous profile_sim_opt_k
            % (seeded from profile_sim on the first iteration) with NaN intensity
            % and pressure fields so downstream plotting treats them as absent.
            if ~exist('profile_sim_opt_k', 'var')
                profile_sim_opt_k = profile_sim;
            end
            for fn = {'Isppa','Ispta','Ita','p_max','p_min','p_rms'}
                if isfield(profile_sim_opt_k, fn{1})
                    profile_sim_opt_k.(fn{1}) = NaN;
                end
            end
        end

        plot_opt_sim_results(...
            opt_param, ...
            profile_target_k, ...
            profile_analytical, ...
            profile_analytical_opt_k, ...
            profile_sim, ...
            profile_sim_opt_k, ...
            min_err)

        save_optimized_values(opt_param);

    end % intensity loop

    % =========================================================================
    % STAGE 7 — Cache cleanup
    %
    % k-Wave writes large .mat result files to the cache folder during
    % calibration.  Remove them once all simulations are complete unless the
    % caller explicitly requested to keep them for debugging by setting
    % parameters.io.save_acoustic_matrices = 1 or save_matrices = 1.
    % =========================================================================
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
