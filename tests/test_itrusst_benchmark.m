classdef test_itrusst_benchmark < matlab.unittest.TestCase
% TEST_ITRUSST_BENCHMARK  ITRUSST intercomparison benchmark 1 (water, bowl source).
%
% Runs benchmark PH1-BM1-SC1 of Aubry et al. (2022), "Benchmark problems for
% transcranial ultrasound simulation: Intercomparison of compressional wave
% models", JASA 152(2):1003, through PRESTUS_PIPELINE and compares the
% pressure amplitude in the central plane with the published k-Wave reference
% (60 points per wavelength, axisymmetric):
%
%   source  focused bowl, ROC 64 mm, aperture 64 mm, 500 kHz, 60 kPa
%   medium  lossless water, c = 1500 m/s, rho = 1000 kg/m^3
%   plane   70 mm (lateral) x 120 mm (axial) from the bowl apex, 0.5 mm grid
%
% Metrics follow the intercomparison (computeDifferenceMetrics): differences
% in focal pressure and focal position, and the L-infinity / L2 norms of the
% field difference relative to the reference maximum, beyond the exit plane.
%
% Data: download INTERCOMPARISON-PH1-DATA-V1.0 (doi:10.5281/zenodo.6020543)
% and set PRESTUS_ITRUSST_DATA to that folder (it holds KWAVE/). Tests skip
% when it is not set.
%
% Tags:
%   'itrusst_axisym'  PRESTUS on a 2D axisymmetric grid (minutes on CPU)
%   'itrusst_3d'      PRESTUS on a 3D grid (GPU/HPC; PRESTUS_PHYSICS_CODE_TYPE)
%
% Before the 2026-10 fix of the axisymmetric source mask, PRESTUS exceeded the
% reference focal pressure by 6 % (0.5 mm) to 11 % (0.25 mm). Tolerances: see
% Tol below and doc_testing.md.
%
%   run_all_tests('itrusst')
%   runtests('tests/test_itrusst_benchmark.m', 'Tag', 'itrusst_axisym')

    properties (Constant)
        RocMm      = 64
        ApertureMm = 64
        FreqHz     = 500e3
        AmpPa      = 60e3
        ResMm      = 0.5
        % Tolerances per geometry. Axisymmetric: first passing run (2026-10-02,
        % R2023b CPU: peak +0.2 %, shift 0 mm, Linf 2.1 %, L2 1.5 %) with margin.
        % 3D: provisional until its first run on the reference GPU node.
        %   peak_rel  focal pressure;  focus_mm  focal position shift
        %   linf      max |p - p_ref| / max(p_ref);  l2  ||p - p_ref|| / ||p_ref||
        Tol = struct( ...
            'axisym', struct('peak_rel', 0.02, 'focus_mm', 0.5, 'linf', 0.04, 'l2', 0.03), ...
            'threed', struct('peak_rel', 0.05, 'focus_mm', 1.0, 'linf', 0.10, 'l2', 0.06))
    end

    properties
        Reference   % [axial x lateral] p_amp [Pa]
    end

    methods (TestClassSetup)
        function load_reference(tc)
            here = fileparts(mfilename('fullpath'));
            repo = fileparts(here);
            addpath(genpath(fullfile(repo, 'functions')));
            addpath(genpath(fullfile(repo, 'external')));
            addpath(fullfile(here, 'fixtures'));

            data = getenv('PRESTUS_ITRUSST_DATA');
            tc.assumeNotEmpty(data, 'PRESTUS_ITRUSST_DATA not set — skipping ITRUSST benchmark.');
            ref_file = fullfile(data, 'KWAVE', 'PH1-BM1-SC1_KWAVE.mat');
            tc.assumeTrue(isfile(ref_file), sprintf('Reference not found: %s', ref_file));
            ref = load(ref_file, 'p_amp');
            p = double(ref.p_amp);
            if size(p, 1) < size(p, 2)
                p = p.';                      % -> [axial, lateral]
            end
            tc.Reference = p;
        end
    end

    methods (Test, TestTags = {'itrusst_axisym'})
        function test_bm1_axisym(tc)
            tc.verify_against_reference(tc.run_bm1('axisym'), tc.Tol.axisym);
        end
    end

    methods (Test, TestTags = {'itrusst_3d'})
        function test_bm1_3d(tc)
            tc.verify_against_reference(tc.run_bm1('3d'), tc.Tol.threed);
        end
    end

    methods (Access = private)

        function plane = run_bm1(tc, geometry)
            % PRESTUS field in the reference plane, [axial x lateral]
            [n_ax, n_lat] = size(tc.Reference);
            half = (n_lat - 1) / 2;

            sim_dir = tempname;
            mkdir(sim_dir);
            tc.addTeardown(@() rmdir(sim_dir, 's'));
            prev = warning('off', 'all');
            restore = onCleanup(@() warning(prev));

            p = make_minimal_parameters(sim_dir);
            p.modules.run_nifti_creation     = 0;
            p.modules.store_acoustic_results = 1;
            p.modules.run_water_baseline     = 0;
            p.transducer.target_isppa_wcm2   = NaN;
            p.transducer.freq_hz                 = tc.FreqHz;
            p.transducer.annular.elem_n          = 1;
            p.transducer.annular.elem_id_mm      = 0;
            p.transducer.annular.elem_od_mm      = tc.ApertureMm;
            p.transducer.annular.curv_radius_mm  = tc.RocMm;
            p.transducer.annular.dist_geom_ep_mm = [];
            p.transducer.annular.elem_p_amp      = tc.AmpPa;
            p.transducer.annular.elem_phase_deg  = 0;
            p.transducer.annular.elem_phase_rad  = 0;
            % Lossless water as in the benchmark
            p.medium_properties.water.sound_speed = 1500;
            p.medium_properties.water.density     = 1000;
            p.medium_properties.water.alpha_coeff = 0;
            p.grid.resolution_mm = tc.ResMm;
            code_type = getenv('PRESTUS_PHYSICS_CODE_TYPE');
            if ~isempty(code_type); p.simulation.code_type = code_type; end

            margin = 8;          % voxels around the plane; keeps the bowl off the grid edge
            focus_vox = round(tc.RocMm / tc.ResMm);
            switch geometry
                case 'axisym'
                    c = half + 1 + margin;
                    p.grid.axisymmetric  = 1;
                    p.grid.default_dims  = [2 * (half + margin) + 1, n_ax + 2 * margin];
                    p.transducer.trans_pos = [c, 1 + margin];
                    p.transducer.focus_pos = [c, 1 + margin + focus_vox];
                case '3d'
                    c = half + 1 + margin;
                    p.grid.axisymmetric  = 0;
                    p.grid.default_dims  = [2 * (half + margin) + 1, 2 * (half + margin) + 1, n_ax + 2 * margin];
                    p.transducer.trans_pos = [c, c, 1 + margin];
                    p.transducer.focus_pos = [c, c, 1 + margin + focus_vox];
            end

            evalc('q = prestus_pipeline(p);');
            pm = double(gather(q.benchmark_sensor_data.p_max_all));
            tp = q.transducer(1).trans_pos;
            if ismatrix(pm)
                sub = pm(tp(1) - half : tp(1) + half, tp(2) : tp(2) + n_ax - 1);     % [lateral, axial]
            else
                sub = squeeze(pm(tp(1) - half : tp(1) + half, tp(2), tp(3) : tp(3) + n_ax - 1));
            end
            plane = sub.';
        end

        function verify_against_reference(tc, plane, tol)
            ref = tc.Reference;
            tc.assertSize(plane, size(ref));
            % Exclude the region behind the exit plane (inside the bowl)
            ep_depth_mm = tc.RocMm - sqrt(tc.RocMm^2 - (tc.ApertureMm / 2)^2);
            ep = ceil(ep_depth_mm / tc.ResMm) + 1;
            ref_f = ref(ep:end, :);
            sim_f = plane(ep:end, :);

            [ref_max, i_ref] = max(ref_f(:));
            [sim_max, i_sim] = max(sim_f(:));
            [ra, rl] = ind2sub(size(ref_f), i_ref);
            [sa, sl] = ind2sub(size(sim_f), i_sim);
            shift_mm = hypot(sa - ra, sl - rl) * tc.ResMm;
            linf = max(abs(sim_f(:) - ref_f(:))) / ref_max;
            l2   = norm(sim_f(:) - ref_f(:)) / norm(ref_f(:));

            fprintf(['ITRUSST BM1: focal p %.3f MPa (ref %.3f), focus shift %.2f mm, ' ...
                     'Linf %.1f %%, L2 %.1f %%\n'], sim_max / 1e6, ref_max / 1e6, shift_mm, ...
                     100 * linf, 100 * l2);
            tc.verifyEqual(sim_max, ref_max, 'RelTol', tol.peak_rel, 'Focal pressure');
            tc.verifyLessThanOrEqual(shift_mm, tol.focus_mm, 'Focal position shift [mm]');
            tc.verifyLessThanOrEqual(linf, tol.linf, 'L-infinity difference');
            tc.verifyLessThanOrEqual(l2, tol.l2, 'L2 difference');
        end

    end

end
