classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_elastic_backend < matlab.unittest.TestCase
% TEST_ELASTIC_BACKEND  Elastic (shear-wave) skull backend: acoustic_backend 'kwave_elastic'.
%
% Covers MEDIUM_ELASTIC_PROPERTIES (Pichardo et al. 2017 shear mapping, density
% threshold, Poisson-ratio floor, Kelvin-Voigt absorption conversion), the
% matched BabelBrain skull properties (BABELBRAIN_SKULL_ELASTIC), config
% validation, and — with a small 2D free-water run — that the pressure→stress
% source conversion in ACOUSTIC_SIMULATION_ELASTIC reproduces the fluid
% k-Wave focal pressure.
%
%   results = runtests('tests/test_elastic_backend.m');


    methods (Static, Access = private)
        function p = params()
            p = struct();
            p.medium_properties = struct( ...
                'water', struct('sound_speed', 1500, 'density', 994), ...
                'brain', struct('sound_speed', 1546, 'density', 1046), ...
                'skull', struct('sound_speed', 2800, 'density', 1850), ...
                'skull_trabecular', struct('sound_speed', 2300, 'density', 1700), ...
                'skull_cortical', struct('sound_speed', 2800, 'density', 1850));
            p.transducer = struct('freq_hz', 500e3);
            p.pct = struct('enabled', 0);
            p.grid = struct('axisymmetric', 0, 'resolution_mm', 0.4);
            p.simulation = struct('code_type', 'matlab_cpu', 'precision', 'single', 'debug', 0);
        end

        function [medium, masks] = skull_medium()
            % labels: 1 water | 3 skull | 4 trabecular | 5 cortical | 2 brain
            masks  = [1 3 4 5 2];
            medium = struct('sound_speed', [1500 2800 2300 2800 1546], ...
                'density', [994 1850 1700 1850 1046], ...
                'alpha_coeff', [0.0022 13.3 13.3 13.3 0.59], 'alpha_power', 1);
        end
    end

    methods (Test, TestTags = {'elastic'})

        function test_pichardo_shear_speed_at_source_frequency(tc)
            % published fits, interpolated to 500 kHz, applied to layer densities
            [medium, masks] = tc.skull_medium();
            em = medium_elastic_properties(tc.params(), medium, masks);
            t = (500e3 - 270e3) / (836e3 - 270e3);
            a = 0.458 + (0.386 - 0.458) * t; b = 545.52 + (815.51 - 545.52) * t;
            tc.verifyEqual(em.sound_speed_shear, [0, a*[1850 1700 1850] + b, 0], 'AbsTol', 1e-9);
        end

        function test_pichardo_fits_reproduced_at_reference_frequencies(tc)
            [medium, masks] = tc.skull_medium();
            p = tc.params(); p.grid.resolution_mm = 0.1;
            p.transducer.freq_hz = 270e3;
            em = medium_elastic_properties(p, medium, masks);
            tc.verifyEqual(em.sound_speed_shear(4), 0.458*1850 + 545.52, 'AbsTol', 1e-9);
            p.transducer.freq_hz = 836e3;
            em = medium_elastic_properties(p, medium, masks);
            tc.verifyEqual(em.sound_speed_shear(4), 0.386*1850 + 815.51, 'AbsTol', 1e-9);
        end

        function test_shear_speed_extrapolation_warns(tc)
            [medium, masks] = tc.skull_medium();
            p = tc.params(); p.transducer.freq_hz = 250e3;
            tc.verifyWarning(@() medium_elastic_properties(p, medium, masks), ...
                'medium_elastic_properties:shearSpeedExtrapolated');
        end

        function test_density_threshold(tc)
            [medium, masks] = tc.skull_medium();
            medium.density(3) = 1150;   % trabecular voxel below shear_density_min
            em = medium_elastic_properties(tc.params(), medium, masks);
            tc.verifyEqual(em.sound_speed_shear(3), 0);
        end

        function test_poisson_floor(tc)
            [medium, masks] = tc.skull_medium();
            medium.sound_speed(3) = 1800;   % cS 1300 would give nu < 0.1
            em = medium_elastic_properties(tc.params(), medium, masks);
            cs = em.sound_speed_shear(3); cl = 1800;
            nu = (cl^2 - 2*cs^2) / (2*(cl^2 - cs^2));
            tc.verifyEqual(nu, 0.1, 'AbsTol', 1e-9);
            t = 230 / 566;
            tc.verifyEqual(em.sound_speed_shear(4), (0.458 - 0.072*t)*1850 + 545.52 + 269.99*t, 'AbsTol', 1e-9);   % untouched
        end

        function test_mapping_off_and_no_mask(tc)
            [medium, masks] = tc.skull_medium();
            p = tc.params(); p.elastic = struct('shear_mapping', 'off');
            em = medium_elastic_properties(p, medium, masks);
            tc.verifyEqual(nnz(em.sound_speed_shear), 0);
            tc.verifyEqual(nnz(em.alpha_coeff_shear), 0);
            em = medium_elastic_properties(tc.params(), medium, []);
            tc.verifyEqual(nnz(em.sound_speed_shear), 0);
        end

        function test_attenuation_exact_at_source_frequency(tc)
            % every coefficient a*f0^y must equal a2*f0^2 (Kelvin-Voigt, y = 2)
            [medium, masks] = tc.skull_medium();
            p = tc.params(); f = 0.5;
            em = medium_elastic_properties(p, medium, masks);
            % shear: 43.63 dB/(cm MHz^1.662); skull longitudinal: 27.31 dB/(cm MHz^1.722)
            tc.verifyEqual(em.alpha_coeff_shear(2:4) * f^2, repmat(43.63 * f^1.662, 1, 3), 'RelTol', 1e-12);
            tc.verifyEqual(em.alpha_coeff_compression(2:4) * f^2, repmat(27.31 * f^1.722, 1, 3), 'RelTol', 1e-12);
            % outside the skull: fluid value
            tc.verifyEqual(em.alpha_coeff_compression(5) * f^2, 0.59 * f, 'RelTol', 1e-12);
            tc.verifyEqual(em.alpha_coeff_shear([1 5]), [0 0]);
        end

        function test_long_attenuation_medium_keeps_fluid_values(tc)
            [medium, masks] = tc.skull_medium();
            p = tc.params(); p.elastic = struct('long_attenuation', 'medium');
            em = medium_elastic_properties(p, medium, masks);
            tc.verifyEqual(em.alpha_coeff_compression(2:4) * 0.5^2, repmat(13.3 * 0.5, 1, 3), 'RelTol', 1e-12);
        end

        function test_attenuation_reproduces_pichardo_values(tc)
            % defaults must give the published Np/m at 270 and 836 kHz
            [medium, masks] = tc.skull_medium();
            p = tc.params(); p.grid.resolution_mm = 0.1;
            to_npm = @(a2, f) a2 * f^2 * 100 / (20*log10(exp(1)));
            for c = {{270e3, 57, 33}, {836e3, 373, 231}}
                [f0, a_s, a_l] = c{1}{:};
                p.transducer.freq_hz = f0;
                em = medium_elastic_properties(p, medium, masks);
                tc.verifyEqual(to_npm(em.alpha_coeff_shear(2), f0/1e6), a_s, 'RelTol', 2e-3);
                tc.verifyEqual(to_npm(em.alpha_coeff_compression(2), f0/1e6), a_l, 'RelTol', 2e-3);
            end
        end

        function test_babelbrain_skull_elastic_matches_kwave_elastic(tc)
            % babelbrain.shear 'prestus_elastic': per skull layer, the kwave_elastic
            % shear speed, and both attenuations in Np/m at f0
            [medium, masks] = tc.skull_medium();
            p = tc.params();
            em = medium_elastic_properties(p, medium, masks);
            sk = babelbrain_skull_elastic(p);
            tc.verifyEqual(sort(fieldnames(sk))', {'skull', 'skull_cortical', 'skull_trabecular'});
            np = 100 / (20*log10(exp(1)));
            % voxels 2/3/4 of skull_medium are skull / trabecular / cortical
            cases = {'skull', 2; 'skull_trabecular', 3; 'skull_cortical', 4};
            for k = 1:size(cases, 1)
                s = sk.(cases{k, 1}); v = cases{k, 2};
                tc.verifyEqual(s.shear_sound_speed, em.sound_speed_shear(v), 'AbsTol', 1e-9);
                tc.verifyEqual(s.density, medium.density(v));
                tc.verifyEqual(s.sound_speed, medium.sound_speed(v));
                tc.verifyEqual(s.shear_att_np_m, 43.63 * 0.5^1.662 * np, 'RelTol', 1e-12);
                tc.verifyEqual(s.long_att_np_m,  27.31 * 0.5^1.722 * np, 'RelTol', 1e-12);
            end
        end

        function test_babelbrain_skull_elastic_reproduces_pichardo_values(tc)
            % Np/m conversion: the defaults give Pichardo et al. 2017's values
            p = tc.params(); p.transducer.freq_hz = 270e3;
            sk = babelbrain_skull_elastic(p);
            tc.verifyEqual(sk.skull_cortical.shear_att_np_m, 57, 'RelTol', 2e-3);
            tc.verifyEqual(sk.skull_cortical.long_att_np_m,  33, 'RelTol', 2e-3);
            tc.verifyEqual(sk.skull_cortical.shear_sound_speed, 0.458*1850 + 545.52, 'AbsTol', 1e-9);
        end

        function test_babelbrain_skull_elastic_long_attenuation_medium(tc)
            p = tc.params(); p.elastic = struct('long_attenuation', 'medium');
            p.medium_properties.skull_cortical.alpha_coeff = 13.3;
            p.medium_properties.skull_cortical.alpha_power = 1;
            sk = babelbrain_skull_elastic(p);
            tc.verifyEqual(sk.skull_cortical.long_att_np_m, db_cm_mhz_to_np_m(13.3, 500e3, 1), 'RelTol', 1e-12);
        end

        function test_babelbrain_skull_elastic_needs_skull(tc)
            p = tc.params();
            p.medium_properties = rmfield(p.medium_properties, {'skull', 'skull_cortical', 'skull_trabecular'});
            tc.verifyError(@() babelbrain_skull_elastic(p), 'babelbrain_skull_elastic:noSkull');
        end

        function test_shear_ppw_warning(tc)
            [medium, masks] = tc.skull_medium();
            p = tc.params(); p.grid.resolution_mm = 0.5;   % ~1384 m/s at 500 kHz -> 5.5 PPW
            tc.verifyWarning(@() medium_elastic_properties(p, medium, masks), ...
                'medium_elastic_properties:shearPPW');
        end

        function test_elastic_matches_fluid_in_water(tc)
            % Free-water 2D run: the stress source (sxx = syy = -p) must give the
            % same pressure field as the fluid pressure source.
            Nx = 96; Ny = 96; dx = 0.5e-3;
            kgrid = kWaveGrid(Nx, dx, Ny, dx);
            c0 = 1500; rho0 = 1000; f0 = 500e3;
            kgrid.makeTime(c0, 0.1, 50e-6);
            p_mask = false(Nx, Ny);
            p_mask(12, 36:60) = true;             % line source
            sig = 1e5 * sin(2*pi*f0*kgrid.t_array);
            source = struct('p_mask', p_mask, 'p', sig);
            sensor = struct('mask', true(Nx, Ny), 'record', {{'p_max_all'}});
            args = {'PMLInside', false, 'PMLSize', 20, 'PlotSim', false, 'DataCast', 'single'};

            fluid = struct('sound_speed', c0, 'density', rho0);
            sd_f = kspaceFirstOrder2D(kgrid, fluid, source, sensor, args{:});

            p = tc.params();
            p.acoustic_backend = 'kwave_elastic';
            sd_e = acoustic_simulation_elastic(kgrid, fluid, source, sensor, ...
                struct('PMLInside', false, 'PMLSize', 20, 'PlotSim', false), p);

            % Compare beyond the near field (pstd vs k-space dispersion differs slightly)
            roi_f = double(gather(sd_f.p_max_all(40:end, 30:66)));
            roi_e = double(gather(sd_e.p_max_all(40:end, 30:66)));
            tc.verifyEqual(max(roi_e(:)), max(roi_f(:)), 'RelTol', 0.05);
            tc.verifyLessThan(norm(roi_e(:) - roi_f(:)) / norm(roi_f(:)), 0.05);
        end

        function test_config_accepts_backend_and_overrides_cpp(tc)
            o = struct('acoustic_backend', 'kwave_elastic');
            o.simulation = struct('debug', 0, 'code_type', 'cpp_gpu');
            o.grid = struct('axisymmetric', 0);
            p = tc.verifyWarning(@() load_parameters(o), 'prestus:elasticNoCpp');
            tc.verifyEqual(char(p.simulation.code_type), 'matlab_gpu');
            tc.verifyEqual(char(p.elastic.shear_mapping), 'pichardo');
        end

        function test_config_overrides_cpp_mlx_to_matlab_cpu(tc)
            % The Apple GPU backend has no gpuArray: the elastic override must land on matlab_cpu
            o = struct('acoustic_backend', 'kwave_elastic');
            o.simulation = struct('debug', 0, 'code_type', 'cpp_mlx');
            o.grid = struct('axisymmetric', 0);
            p = tc.verifyWarning(@() load_parameters(o), 'prestus:elasticNoCpp');
            tc.verifyEqual(char(p.simulation.code_type), 'matlab_cpu');
        end

        function test_config_rejects_axisymmetric(tc)
            o = struct('acoustic_backend', 'kwave_elastic');
            o.simulation = struct('debug', 0);
            o.grid = struct('axisymmetric', 1);
            tc.verifyError(@() load_parameters(o), 'prestus:elasticAxisymmetric');
        end
    end
end
