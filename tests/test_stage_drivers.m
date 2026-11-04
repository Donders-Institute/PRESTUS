classdef test_stage_drivers < matlab.unittest.TestCase
% TEST_STAGE_DRIVERS  Stage drivers on tiny synthetic inputs (no subject data, no GPU).
%
% Covers THERMAL_SIMULATION (k-Wave diffusion on a 16^3 homogeneous brain grid
% with a synthetic Gaussian pressure field: heating at the focus, no rise for a
% zero source, temperature rise linear in intensity, analytic early-time rise),
% THERMAL_STAGE (module/acoustics gating, run -> cache -> reload, missing cache
% with overwrite 'never'), CACHE_ACOUSTIC (save/load round trip, ISPPA rescaling
% on load, stale-grid and corrupt-cache errors, kwave_medium flag) and
% ACOUSTIC_STAGE's cache paths (no simulation requested, cache
% hit, baseline taken from the live run, stale cache without a source, source
% label size mismatch).
%
% Not covered here: an actual acoustic simulation (ACOUSTIC_WRAPPER, k-Wave),
% which test_integration_water runs at the 'water' level, and the jwave /
% BabelBrain backends (standalone suites, need Python and/or a GPU).
%
%   results = runtests('tests/test_stage_drivers.m');

    properties
        Dir
    end

    methods (TestClassSetup)
        function add_paths(tc)
            repo = fileparts(fileparts(mfilename('fullpath')));
            addpath(genpath(fullfile(repo, 'functions')));
            addpath(genpath(fullfile(repo, 'external')));
            addpath(fullfile(repo, 'config'));
            tc.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture( ...
                {'prestus:toolboxDefault', 'prestus:noSimNIBS', 'prestus:missingCharm'}));
        end
    end

    methods (TestMethodSetup)
        function make_dir(tc)
            tc.Dir = tempname;
            mkdir(tc.Dir);
            tc.addTeardown(@() rmdir(tc.Dir, 's'));
            tc.addTeardown(@() close(findall(groot, 'Type', 'figure')));
        end
    end

    % ------------------------------------------------------------------ %
    %% thermal_simulation
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'thermal_simulation'})

        function test_thermal_heats_at_focus(tc)
            p = stage_params(tc.Dir);
            [kgrid, med, masks] = brain_medium(p);
            [~, tss, res] = run_thermal(p, gaussian_p(p, 1e6), kgrid, med, masks);
            rise = res.maxT - 37;
            [mx, idx] = max(rise(:));
            [i, j, k] = ind2sub(size(rise), idx);
            tc.verifyGreaterThan(mx, 1e-3);
            tc.verifyEqual([i j k], p.transducer.focus_pos);
            tc.verifyGreaterThanOrEqual(min(rise(:)), -1e-6);       % nothing cools below baseline
            % Protocol: 2 pulses (ON+OFF) then 2 post-PTRD steps -> 1.4 s
            tc.verifyEqual([tss.time], sort([tss.time]));
            tc.verifyEqual(tss(end).time, 1.4, 'AbsTol', 1e-9);
            tc.verifyEqual(size(res.focal_planeT, 3), numel(tss));
        end

        function test_thermal_early_rise_matches_analytic(tc)
            % Early heating at the focus ~ Q * t_on / (rho * c_p): diffusion length
            % over 0.4 s (~0.3 mm) << Gaussian width (2 mm), perfusion ~0.4 %.
            p = stage_params(tc.Dir);
            [kgrid, med, masks] = brain_medium(p);
            p0 = 1e6;
            [~, ~, res] = run_thermal(p, gaussian_p(p, p0), kgrid, med, masks);
            b = p.medium_properties.brain;
            alpha_np = b.alpha_coeff * (p.transducer.freq_hz / 1e6)^b.alpha_power ...
                * 100 / (20 * log10(exp(1)));                              % [Np/m]
            Q = alpha_np * p0^2 / (b.density * b.sound_speed);             % [W/m^3]
            t_on = 2 * p.timing.pd;
            expected = Q * t_on / (b.density * b.specific_heat_capacity);
            f = p.transducer.focus_pos;
            tc.verifyEqual(res.maxT(f(1), f(2), f(3)) - 37, expected, 'RelTol', 0.1);
        end

        function test_thermal_zero_source_no_rise(tc)
            p = stage_params(tc.Dir);
            [kgrid, med, masks] = brain_medium(p);
            [~, ~, res] = run_thermal(p, gaussian_p(p, 0), kgrid, med, masks);
            tc.verifyLessThan(max(abs(res.maxT(:) - 37)), 1e-6);
            tc.verifyEqual(max(res.CEM43_iso(:)), 0);
        end

        function test_thermal_rise_linear_in_intensity(tc)
            % Pennes bio-heat is linear in Q ~ p^2: doubling pressure -> 4x rise
            p = stage_params(tc.Dir);
            [kgrid, med, masks] = brain_medium(p);
            [~, ~, r1] = run_thermal(p, gaussian_p(p, 1e6), kgrid, med, masks);
            [~, ~, r2] = run_thermal(p, gaussian_p(p, 2e6), kgrid, med, masks);
            d1 = r1.maxT - 37;  d2 = r2.maxT - 37;
            f = p.transducer.focus_pos;
            tc.verifyEqual(d2(f(1), f(2), f(3)) / d1(f(1), f(2), f(3)), 4, 'RelTol', 1e-3);
        end

        function test_thermal_absorption_coeff_matches_fraction(tc)
            % absorption_coeff = a * alpha_coeff (same power law) heats like absorption_fraction = a
            p = stage_params(tc.Dir);
            [kgrid, med, masks] = brain_medium(p);
            p.medium_properties.brain.absorption_fraction = 0.5;
            [~, ~, r_frac] = run_thermal(p, gaussian_p(p, 1e6), kgrid, med, masks);
            p.medium_properties.brain.absorption_fraction = [];
            p.medium_properties.brain.absorption_coeff = 0.5 * p.medium_properties.brain.alpha_coeff;
            [~, ~, r_coeff] = run_thermal(p, gaussian_p(p, 1e6), kgrid, med, masks);
            tc.verifyGreaterThan(max(r_frac.maxT(:)) - 37, 1e-3);
            tc.verifyEqual(r_coeff.maxT, r_frac.maxT, 'AbsTol', 1e-9);
        end

    end

    % ------------------------------------------------------------------ %
    %% thermal_stage
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'thermal_stage'})

        function test_thermal_stage_skips_when_disabled(tc)
            p = stage_params(tc.Dir);
            p.modules.run_heating_sims = 0;
            out = call_thermal_stage(p, [], [], struct(), struct());
            tc.verifyEqual(out.parameters.state.heating_available, 0);
            tc.verifyEqual(out.results_heating, struct());
            tc.verifyEmpty(dir(fullfile(tc.Dir, '*heating_res*')));
        end

        function test_thermal_stage_skips_without_acoustics(tc)
            p = stage_params(tc.Dir);
            p.modules.run_heating_sims = 1;
            p.state.acoustics_available = 0;
            out = call_thermal_stage(p, [], [], struct(), struct());
            tc.verifyEqual(out.parameters.state.heating_available, 0);
            tc.verifyEmpty(dir(fullfile(tc.Dir, '*heating_res*')));
        end

        function test_thermal_stage_runs_then_reloads_cache(tc)
            p = stage_params(tc.Dir);
            p.modules.run_heating_sims = 1;
            p.state.acoustics_available = 1;
            p.io.save_thermal_matrices = 1;
            [kgrid, med, masks, plus] = brain_medium(p);
            sd = struct('p_max_all', gaussian_p(p, 1e6));
            out1 = call_thermal_stage(p, sd, kgrid, med, plus, masks);
            tc.verifyEqual(out1.parameters.state.heating_available, 1);
            cache = fullfile(tc.Dir, 'sub-001_water_heating_res.mat');
            tc.verifyTrue(isfile(cache));
            tc.verifyGreaterThan(max(out1.results_heating.maxT(:)), 37);

            % Second run: overwrite 'never' -> load, even with a zero field passed in
            p.io.overwrite_files = 'never';
            sd0 = struct('p_max_all', zeros(p.grid.dims));
            out2 = call_thermal_stage(p, sd0, kgrid, med, plus, masks);
            tc.verifyEqual(out2.parameters.state.heating_available, 1);
            tc.verifyEqual(out2.results_heating.maxT, out1.results_heating.maxT);
            tc.verifyEqual([out2.time_status_seq.time], [out1.time_status_seq.time]);
        end

        function test_thermal_stage_never_without_cache_simulates(tc)
            % overwrite_files 'never' only protects an existing cache: with
            % none, CONFIRM_OVERWRITING allows the run. (The
            % 'thermalNoAcoustics' warning branch is reachable only when an
            % interactive user declines the confirmation dialog.)
            p = stage_params(tc.Dir);
            p.modules.run_heating_sims = 1;
            p.state.acoustics_available = 1;
            p.io.overwrite_files = 'never';
            [kgrid, med, masks, plus] = brain_medium(p);
            sd = struct('p_max_all', gaussian_p(p, 1e6));
            out = call_thermal_stage(p, sd, kgrid, med, plus, masks);
            tc.verifyEqual(out.parameters.state.heating_available, 1);
        end

    end

    % ------------------------------------------------------------------ %
    %% cache_acoustic
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'cache_acoustic'})

        function test_cache_round_trip_without_rescale(tc)
            p = stage_params(tc.Dir);
            f = fullfile(tc.Dir, 'c.mat');
            sd = struct('p_max_all', rand(p.grid.dims) * 1e5, 'p_min_all', -rand(p.grid.dims));
            prov = struct('freefield_isppa_wcm2', 30, 'freefield_source', 'kwave_simulation');
            cache_acoustic('save', f, p, sd, struct('dx', 1e-3), prov, struct('backend', 'kwave'));
            [sd2, kg2, prov2, info2] = cache_acoustic('load', f, p);
            tc.verifyEqual(sd2.p_max_all, sd.p_max_all);
            tc.verifyEqual(sd2.p_min_all, sd.p_min_all);
            tc.verifyFalse(isfield(sd2, 'isppa_scale_p'));
            tc.verifyEqual(kg2.dx, 1e-3);
            tc.verifyEqual(prov2.freefield_isppa_wcm2, 30);
            tc.verifyEqual(info2.backend, 'kwave');
        end

        function test_cache_load_rescales_to_target(tc)
            p = stage_params(tc.Dir);                     % target 30 W/cm^2
            f = fullfile(tc.Dir, 'c.mat');
            sd = struct('p_max_all', ones(p.grid.dims), 'p_min_all', -ones(p.grid.dims));
            prov = struct('freefield_isppa_wcm2', 15, 'freefield_source', 'kwave_simulation');
            cache_acoustic('save', f, p, sd, [], prov);
            sd2 = cache_acoustic('load', f, p);
            tc.verifyEqual(sd2.p_max_all, sqrt(2) * ones(p.grid.dims), 'RelTol', 1e-12);
            tc.verifyEqual(sd2.p_min_all, -sqrt(2) * ones(p.grid.dims), 'RelTol', 1e-12);
            tc.verifyEqual(sd2.isppa_scale_p, sqrt(2), 'RelTol', 1e-12);
        end

        function test_cache_stale_grid_errors_and_keeps_file(tc)
            p = stage_params(tc.Dir);
            f = fullfile(tc.Dir, 'c.mat');
            cache_acoustic('save', f, p, struct('p_max_all', ones(p.grid.dims)), [], ...
                struct('freefield_isppa_wcm2', 30));
            q = p;
            q.grid.dims = [16 16 18];
            tc.verifyError(@() cache_acoustic('load', f, q), 'cache_acoustic:staleGrid');
            tc.verifyTrue(isfile(f));
        end

        function test_cache_missing_variable_is_corrupt(tc)
            p = stage_params(tc.Dir);
            f = fullfile(tc.Dir, 'c.mat');
            sensor_data = struct('p_max_all', 1); %#ok<NASGU>
            save(f, 'sensor_data');
            tc.verifyError(@() cache_acoustic('load', f, p), 'cache_acoustic:corruptCache');
        end

        function test_cache_kwave_medium_follows_flag(tc)
            p = stage_params(tc.Dir);
            med = struct('density', ones(2));
            p.io.save_kwave_medium = 0;
            f1 = fullfile(tc.Dir, 'slim.mat');
            cache_acoustic('save', f1, p, struct('p_max_all', 1), [], struct(), struct(), med);
            tc.verifyFalse(ismember('kwave_medium', {whos('-file', f1).name}));
            p.io.save_kwave_medium = 1;
            f2 = fullfile(tc.Dir, 'full.mat');
            cache_acoustic('save', f2, p, struct('p_max_all', 1), [], struct(), struct(), med);
            tc.verifyTrue(ismember('kwave_medium', {whos('-file', f2).name}));
        end

        function test_cache_kwave_medium_not_saved_by_default(tc)
            % KNOWN DEVIATION (2026-10): cache_acoustic.m documents a slim cache
            % that includes kwave_medium only "when explicitly requested"
            % (async runs). io.save_kwave_medium is not in config_default.yaml,
            % and SHOULD_SAVE_OUTPUT returns true for unlisted fields, so every
            % acoustic cache also stores the full medium. Fails until the
            % default or the documentation changes.
            p = stage_params(tc.Dir);
            if isfield(p.io, 'save_kwave_medium'); p.io = rmfield(p.io, 'save_kwave_medium'); end
            f = fullfile(tc.Dir, 'default.mat');
            cache_acoustic('save', f, p, struct('p_max_all', 1), [], struct(), struct(), ...
                struct('density', ones(2)));
            tc.verifyFalse(ismember('kwave_medium', {whos('-file', f).name}), ...
                ['Known deviation: kwave_medium is saved when io.save_kwave_medium is unset ' ...
                 '(should_save_output default true; cache_acoustic.m:84).']);
        end

        function test_cache_rejects_unknown_mode(tc)
            p = stage_params(tc.Dir);
            tc.verifyError(@() cache_acoustic('delete', fullfile(tc.Dir, 'c.mat'), p), ?MException);
        end

    end

    % ------------------------------------------------------------------ %
    %% acoustic_stage (cache paths; no simulation)
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'acoustic_stage'})

        function test_acoustic_stage_nothing_requested(tc)
            p = stage_params(tc.Dir);
            p.modules.run_acoustic_sims = 0;
            [sd, p2] = call_acoustic_stage(p, [], struct());
            tc.verifyEmpty(sd);
            tc.verifyEqual(p2.state.acoustics_available, 0);
            tc.verifyEqual(p2.modules.run_acoustic_analysis, 0);
        end

        function test_acoustic_stage_loads_cache_and_scales(tc)
            p = stage_params(tc.Dir);
            p.io.overwrite_files = 'never';
            write_results_cache(p, ones(p.grid.dims), ...
                struct('freefield_isppa_wcm2', 15, 'freefield_source', 'kwave_simulation'));
            [sd, p2, ~, prov] = call_acoustic_stage(p, [], struct());
            tc.verifyEqual(p2.state.acoustics_available, 1);
            tc.verifyEqual(sd.p_max_all, sqrt(2) * ones(p.grid.dims), 'RelTol', 1e-12);
            tc.verifyEqual(prov.input_signature, cache_signature(p, 'acoustic'));
        end

        function test_acoustic_stage_baseline_from_live_run(tc)
            % Cache has no readable baseline: the current run's water baseline is used
            p = stage_params(tc.Dir);
            p.io.overwrite_files = 'never';
            write_results_cache(p, ones(p.grid.dims), struct('freefield_isppa_wcm2', NaN));
            live = struct('freefield_isppa_wcm2', 7.5, 'freefield_source', 'kwave_simulation');
            % KNOWN BUG (2026-10): CACHE_ACOUSTIC applies ISPPA scaling on load
            % (cache_acoustic.m:193) before ACOUSTIC_STAGE substitutes the live
            % baseline (acoustic_stage.m:117), so with a target set the load
            % errors 'prestus:isppa:noBaseline' and the fallback never runs.
            try
                [sd, ~, ~, prov] = call_acoustic_stage(p, [], live);
            catch ME
                tc.verifyFail(sprintf(['Known bug: cached field is scaled before the live ' ...
                    'baseline is filled in (%s: %s)'], ME.identifier, ME.message));
                return
            end
            tc.verifyEqual(prov.freefield_isppa_wcm2, 7.5);
            tc.verifyEqual(sd.p_max_all, 2 * ones(p.grid.dims), 'RelTol', 1e-12);   % sqrt(30/7.5)
        end

        function test_acoustic_stage_stale_cache_without_source(tc)
            p = stage_params(tc.Dir);
            p.io.overwrite_files = 'never';
            write_results_cache(p, ones(p.grid.dims), struct('freefield_isppa_wcm2', 30));
            q = p;
            q.grid.dims = [16 16 18];
            tc.verifyError(@() call_acoustic_stage(q, [], struct()), ...
                'acoustic_stage:staleCacheNoSource');
        end

        function test_acoustic_stage_axisym_cache_expands_inputs(tc)
            % Reloading an axisymmetric phantom field (cached in 3D) expands the
            % half-plane grid, positions and medium like a fresh simulation
            p = stage_params(tc.Dir);
            p.io.overwrite_files = 'never';
            p.simulation.medium = 'phantom';
            p.grid.axisymmetric = 1;
            p.grid.dims = [12 4];
            p.grid.axisym_bilateral_dims = [12 8];
            p.grid.axisym_transposed = true;
            p.transducer.trans_pos = [1 1];
            p.transducer.focus_pos = [9 1];
            p.transducer.trans_pos_bilateral = [1 5];
            p.transducer.focus_pos_bilateral = [9 5];
            prov = struct('freefield_isppa_wcm2', 30, 'input_signature', cache_signature(p, 'acoustic'));
            q = p; q.grid.dims = [8 8 12];
            f = fullfile(p.io.dir_cache, sprintf('sub-%03d_%s_results%s.mat', ...
                p.subject_id, p.simulation.medium, p.io.output_affix));
            cache_acoustic('save', f, q, struct('p_max_all', ones(8, 8, 12)), [], prov, struct());

            masks = 2 * ones(12, 4);
            med = struct('sound_speed', 1500 * ones(12, 4), 'density', 1000);
            [~, sd, p2, kgrid, med2, masks2, seg2, src2] = evalc_out(@() acoustic_stage(p, [], med, [], [], ...
                masks, masks, [], struct(), struct('transf', eye(4)), ...
                p.transducer.trans_pos, p.transducer.focus_pos, []));
            tc.verifyEqual(size(sd.p_max_all), [8 8 12]);
            tc.verifyEqual(p2.grid.dims, [8 8 12]);
            tc.verifyEqual(p2.transducer.trans_pos, [5 5 1]);
            tc.verifyEqual(p2.transducer.focus_pos, [5 5 9]);
            tc.verifyEqual(size(masks2), [8 8 12]);
            tc.verifyEqual(size(seg2), [8 8 12]);
            tc.verifyEqual(size(med2.sound_speed), [8 8 12]);
            tc.verifyEqual(med2.density, 1000);
            tc.verifyEqual([kgrid.Nx kgrid.Ny kgrid.Nz], [8 8 12]);
            tc.verifyEmpty(src2);                           % no source set up: left as is
        end

        function test_acoustic_stage_cache_signature_checked(tc)
            % A field cached with other settings is not reused silently
            p = stage_params(tc.Dir);
            p.io.overwrite_files = 'never';
            q = p; q.grid.resolution_mm = 2;
            write_results_cache(p, ones(p.grid.dims), ...
                struct('freefield_isppa_wcm2', 30, 'input_signature', cache_signature(q, 'acoustic')));
            tc.verifyError(@() call_acoustic_stage(p, [], struct()), 'prestus:cache:signatureMismatch');
            p.io.cache_signature_check = 'off';
            sd = call_acoustic_stage(p, [], struct());
            tc.verifyEqual(sd.p_max_all, ones(p.grid.dims));
        end

        function test_acoustic_stage_clears_mismatched_source_labels(tc)
            p = stage_params(tc.Dir);
            p.io.overwrite_files = 'never';
            write_results_cache(p, ones(p.grid.dims), struct('freefield_isppa_wcm2', 30));
            labels = ones(8, 8, 8);
            tc.verifyWarning(@() call_acoustic_stage(p, labels, struct()), ...
                'prestus_pipeline:sourceLabelsMismatch');
            [~, ~, labels2] = call_acoustic_stage(p, labels, struct());
            tc.verifyEqual(labels2, zeros(p.grid.dims));
        end

    end

end

% =========================================================================
%% Local helpers
% =========================================================================

function p = stage_params(root)
% Defaults with a 16^3, 1 mm, 500 kHz single-transducer setup on CPU.
    p = load_parameters(struct('simulation', struct('debug', 0, 'interactive', 0)));
    p.subject_id = 1;
    p.platform = 'matlab';
    p.simulation.medium = 'water';
    p.simulation.code_type = 'matlab_cpu';
    p.simulation.precision = 'double';
    p.grid.dims = [16 16 16];
    p.grid.resolution_mm = 1;
    p.transducer.freq_hz = 500e3;
    p.transducer.trans_pos = [8 8 1];
    p.transducer.focus_pos = [8 8 10];
    p.transducer.target_isppa_wcm2 = 30;
    p.timing = struct('pd', 0.1, 'pri', 0.2, 'ptd', 0.4, 'pt_timestep', 0.1, ...
        'ptri', 0.4, 'ptrd', 0.4, 'post_ptri_dur', 1, 'post_pt_timestep', 0.5, ...
        'equal_step_duration', 0);
    p.io.dir_cache = root;
    p.io.dir_img = root;
    p.io.dir_output = root;
    p.io.output_affix = '';
    p.io.overwrite_files = 'always';
    p.state = struct('acoustics_available', 1);
end

function [kgrid, kwave_medium, masks, medium_plus] = brain_medium(p)
% Homogeneous brain on the stage grid, built by MEDIUM_SETUP as in the pipeline.
    labels = fieldnames(p.medium_properties);
    masks = find(strcmp(labels, 'brain')) * ones(p.grid.dims);
    [kwave_medium, medium_plus] = medium_setup(p, masks, struct());
    d = p.grid.resolution_mm * 1e-3;
    kgrid = kWaveGrid(p.grid.dims(1), d, p.grid.dims(2), d, p.grid.dims(3), d);
end

function pr = gaussian_p(p, p0)
% Peak pressure p0 [Pa] at the focus, sigma = 2 voxels.
    [x, y, z] = ndgrid(1:p.grid.dims(1), 1:p.grid.dims(2), 1:p.grid.dims(3));
    f = p.transducer.focus_pos;
    pr = p0 * exp(-((x - f(1)).^2 + (y - f(2)).^2 + (z - f(3)).^2) / (2 * 2^2));
end

function [obj, tss, res] = run_thermal(p, pr, kgrid, kwave_medium, masks)
    [~, plus] = medium_setup(p, masks, struct());
    kwave_medium.temp_0 = plus.temp_0;
    if plus.absorption_from_coeff   % as thermal_stage
        kwave_medium.absorption_np = plus.absorption_np;
    else
        kwave_medium.absorption_fraction = plus.absorption_fraction;
    end
    sensor = struct('mask', ones(p.grid.dims));
    source = struct('Q', zeros(p.grid.dims));
    [~, obj, tss, res] = evalc_out(@() thermal_simulation(p, struct('p_max_all', pr), kgrid, ...
        kwave_medium, sensor, source, eye(4), masks));
end

function out = call_thermal_stage(p, sd, kgrid, kwave_medium, medium_plus, masks)
    if nargin < 6; masks = ones(p.grid.dims); end
    if isempty(fieldnames(kwave_medium)); kwave_medium = struct('density', 1); end
    planimg = struct('transf', eye(4));
    sensor = struct('mask', ones(p.grid.dims));
    source = struct('Q', zeros(p.grid.dims));
    fov = struct('active', false);
    [~, out.results_heating, out.time_status_seq, ~, ~, ~, ~, ~, out.parameters] = evalc_out( ...
        @() thermal_stage(p, sd, kgrid, kwave_medium, medium_plus, masks, masks, planimg, ...
        source, sensor, p.transducer.trans_pos, p.transducer.focus_pos, struct(), fov));
end

function write_results_cache(p, pmax, prov)
    f = fullfile(p.io.dir_cache, sprintf('sub-%03d_%s_results%s.mat', ...
        p.subject_id, p.simulation.medium, p.io.output_affix));
    cache_acoustic('save', f, p, struct('p_max_all', pmax), [], prov, struct());
end

function [sd, p2, labels, prov] = call_acoustic_stage(p, labels, live_prov)
    [~, sd, p2, ~, ~, ~, ~, ~, labels, prov] = evalc_out(@() acoustic_stage(p, [], struct(), [], [], ...
        ones(p.grid.dims), ones(p.grid.dims), labels, live_prov, struct('transf', eye(4)), ...
        p.transducer.trans_pos, p.transducer.focus_pos, []));
end

function varargout = evalc_out(fn)
% Run fn with its console output captured; varargout{1} is the text.
    n = nargout - 1;
    outs = cell(1, n);
    varargout{1} = evalc('[outs{1:n}] = fn();');
    varargout(2:n+1) = outs;
end
