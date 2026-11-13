classdef (SharedTestFixtures = {prestus_paths}) test_physics_water < matlab.unittest.TestCase
% TEST_PHYSICS_WATER  Free-water acoustic benchmarks of the full PRESTUS pipeline.
%
% A single-element focused bowl in lossless water (the ITRUSST BM1 source: ROC
% 64 mm, aperture 64 mm, 500 kHz, 60 kPa; strongly focused, so O'Neil's
% solution is accurate) is simulated through PRESTUS_PIPELINE and compared with
% the analytic on-axis solution of O'Neil (k-Wave focusedBowlONeil):
%   - focal position (on-axis pressure maximum)       within FocusTolMm
%   - focal peak pressure                             within PeakRelTol
%   - -3 dB (intensity) axial focal length            within LengthRelTol
% a multi-ring annular array (in phase) against the O'Neil solution summed over
% its rings (axisymmetric; regression for rings placed at the apex plane instead
% of on the bowl, fixed 2026-10), plus two invariants of the linear solver and
% the ISPPA targeting:
%   - doubling the drive doubles the pressure          (axisymmetric)
%   - a target ISPPA on an uncalibrated transducer is met in water, and with
%     modules.run_water_baseline = 0 the run stops instead (axisymmetric)
%
% Tags:
%   'physics_axisym'  2D axisymmetric grid (kspaceFirstOrderAS, 0.5 mm), CPU, minutes
%   'physics_3d'      3D grid (0.5 mm), O'Neil comparison only; hours on CPU,
%                     run on a GPU / HPC node (set PRESTUS_PHYSICS_CODE_TYPE, e.g.
%                     'matlab_gpu' or 'cpp_gpu')
%
% Regression for the axisymmetric source mask that dropped the negative BLI
% side lobes and raised the focal pressure by 6-11 % (fixed 2026-10).
% Tolerances: see Tol below and doc_testing.md (fix them to the first passing
% run on the reference machine and keep them fixed thereafter).
%
%   run_all_tests('physics')                                   % axisymmetric
%   runtests('tests/test_physics_water.m', 'Tag', 'physics_3d')

    properties (Constant)
        RocMm        = 64
        ApertureMm   = 64
        FreqHz       = 500e3
        AmpPa        = 60e3
        % Ring array for the multi-element case: CTX-500 ring diameters [mm]
        RingsIdMm    = [0, 31.6988, 44.2688, 53.6688]
        RingsOdMm    = [31.14, 43.71, 53.11, 60.83]
        % Tolerances per geometry. Axisymmetric: from the first passing run
        % (2026-10-02, R2023b CPU: focus 0.06 mm, peak +0.2 %) with margin.
        % 3D: provisional until its first run on the reference GPU node.
        % Rings: the summed O'Neil peak is flat (±0.5 mm costs 0.08 % pressure),
        % so the peak position is fixed to within a grid step only (first run
        % 2026-10-08, R2023b CPU: focus 0.50 mm, peak +0.10 %, length +0.04 %).
        Tol = struct( ...
            'axisym', struct('focus_mm', 0.5, 'peak_rel', 0.02, 'length_rel', 0.05), ...
            'axisym_rings', struct('focus_mm', 1.0, 'peak_rel', 0.02, 'length_rel', 0.05), ...
            'threed', struct('focus_mm', 1.0, 'peak_rel', 0.05, 'length_rel', 0.10))
    end


    %% Axisymmetric
    methods (Test, TestTags = {'physics_axisym'})

        function test_axisym_matches_oneil(tc)
            out = tc.run_case('axisym', struct());
            tc.verify_against_oneil(out, tc.Tol.axisym);
        end

        function test_axisym_annular_rings_match_oneil(tc)
            rings = [tc.RingsIdMm; tc.RingsOdMm];
            out = tc.run_case('axisym', struct('rings', rings));
            tc.verify_against_oneil(out, tc.Tol.axisym_rings, rings);
        end

        function test_axisym_pressure_is_linear_in_drive(tc)
            a = tc.run_case('axisym', struct('amp', tc.AmpPa));
            b = tc.run_case('axisym', struct('amp', 2 * tc.AmpPa));
            ratio = max(b.p_axis) / max(a.p_axis);
            tc.verifyEqual(ratio, 2, 'RelTol', 0.01, ...
                'Peak pressure must scale linearly with the source amplitude.');
        end

        function test_axisym_target_isppa_met_without_calibration(tc)
            % Uncalibrated custom transducer + target: the 'auto' baseline runs
            % and the field is scaled to the target (regression for the silent
            % raw-drive-level result).
            target = 20;
            out = tc.run_case('axisym', struct('target', target, 'baseline', 'auto'));
            isppa = max(out.p_field(:))^2 / (2 * out.rho * out.c) * 1e-4;
            tc.verifyEqual(isppa, target, 'RelTol', 0.02);
        end

        function test_axisym_unmet_target_stops_run(tc)
            tc.verifyError(@() tc.run_case('axisym', struct('target', 20, 'baseline', 0)), ...
                'prestus:isppa:noBaseline');
        end

    end

    %% 3D
    methods (Test, TestTags = {'physics_3d'})

        function test_3d_matches_oneil(tc)
            out = tc.run_case('3d', struct());
            tc.verify_against_oneil(out, tc.Tol.threed);
        end

    end

    methods (Access = private)

        function out = run_case(tc, geometry, opts)
            if ~isfield(opts, 'amp');      opts.amp = tc.AmpPa; end
            if ~isfield(opts, 'target');   opts.target = NaN; end
            if ~isfield(opts, 'baseline'); opts.baseline = 'auto'; end
            if ~isfield(opts, 'rings');    opts.rings = [0; tc.ApertureMm]; end   % [ID; OD] [mm] per element
            n_el = size(opts.rings, 2);

            sim_dir = tempname;
            mkdir(sim_dir);
            tc.addTeardown(@() rmdir(sim_dir, 's'));
            prev = warning('off', 'all');
            restore = onCleanup(@() warning(prev));

            p = make_minimal_parameters(sim_dir);
            p.modules.run_nifti_creation     = 0;
            p.modules.store_acoustic_results = 1;
            p.modules.run_water_baseline     = opts.baseline;
            p.transducer.annular.curv_radius_mm = tc.RocMm;
            p.transducer.annular.elem_n         = n_el;
            p.transducer.annular.elem_id_mm     = opts.rings(1, :);
            p.transducer.annular.elem_od_mm     = opts.rings(2, :);
            p.transducer.annular.elem_phase_deg = zeros(1, n_el);
            p.transducer.annular.elem_phase_rad = zeros(1, n_el);
            p.transducer.annular.dist_geom_ep_mm = [];   % derived from ROC and aperture
            p.transducer.freq_hz                = tc.FreqHz;
            p.transducer.annular.elem_p_amp     = repmat(opts.amp, 1, n_el);
            p.transducer.target_isppa_wcm2      = opts.target;
            p.medium_properties.water.sound_speed = 1500;   % lossless water, as O'Neil
            p.medium_properties.water.density     = 1000;
            p.medium_properties.water.alpha_coeff = 0;
            code_type = getenv('PRESTUS_PHYSICS_CODE_TYPE');
            if ~isempty(code_type); p.simulation.code_type = code_type; end

            switch geometry
                case 'axisym'
                    res = 0.5;
                    dims = [157, 257];            % [2 x radial, axial] ~ 78 x 128 mm
                    c = (dims(1) + 1) / 2;
                    p.grid.resolution_mm = res;
                    p.grid.axisymmetric  = 1;
                    p.grid.default_dims  = dims;
                    p.transducer.trans_pos = [c, 9];
                    p.transducer.focus_pos = [c, 9 + round(tc.RocMm / res)];
                case '3d'
                    res = 0.5;
                    dims = [157, 157, 257];
                    c = (dims(1) + 1) / 2;
                    p.grid.resolution_mm = res;
                    p.grid.axisymmetric  = 0;
                    p.grid.default_dims  = dims;
                    p.transducer.trans_pos = [c, c, 9];
                    p.transducer.focus_pos = [c, c, 9 + round(tc.RocMm / res)];
            end

            evalc('q = prestus_pipeline(p);');
            pm = double(gather(q.benchmark_sensor_data.p_max_all));
            tp = q.transducer(1).trans_pos;
            if ismatrix(pm)
                line = pm(tp(1), tp(2):end);
            else
                line = squeeze(pm(tp(1), tp(2), tp(3):end));
            end
            out = struct();
            out.res_mm  = res;
            out.z_mm    = (0:numel(line) - 1)' * res;   % from the bowl apex
            out.p_axis  = line(:);
            out.p_field = pm;
            out.rho     = q.medium_properties.water.density;
            out.c       = q.medium_properties.water.sound_speed;
        end

        function verify_against_oneil(tc, out, tol, rings)
            u0 = tc.AmpPa / (out.rho * out.c);              % source particle velocity
            z_fine = (1:0.05:out.z_mm(end))';               % mm
            if nargin < 4
                p_an = focusedBowlONeil(tc.RocMm * 1e-3, tc.ApertureMm * 1e-3, u0, ...
                    tc.FreqHz, out.c, out.rho, z_fine * 1e-3, []);
            else
                p_an = oneil_rings_axial(tc.RocMm * 1e-3, rings * 1e-3, u0, ...
                    tc.FreqHz, out.c, out.rho, z_fine * 1e-3);
            end
            p_an = abs(p_an(:));

            % Compare beyond the near field (z > ROC / 2)
            far = out.z_mm > tc.RocMm / 2;
            [z_sim, p_sim] = refine_peak(out.z_mm(far), out.p_axis(far));
            far_an = z_fine > tc.RocMm / 2;
            [p_an_max, i_an] = max(p_an(far_an));
            z_an_far = z_fine(far_an);
            z_an = z_an_far(i_an);

            L_sim = focal_length_3db(out.z_mm, out.p_axis, tc.RocMm / 2);
            L_an  = focal_length_3db(z_fine, p_an, tc.RocMm / 2);
            fprintf(['O''Neil check: focus %.2f vs %.2f mm, peak %.0f vs %.0f Pa (%+.2f %%), ' ...
                     '-3 dB length %.2f vs %.2f mm (%+.2f %%)\n'], z_sim, z_an, p_sim, p_an_max, ...
                     100 * (p_sim / p_an_max - 1), L_sim, L_an, 100 * (L_sim / L_an - 1));

            tc.verifyEqual(z_sim, z_an, 'AbsTol', tol.focus_mm, ...
                sprintf('Focal position: simulated %.2f mm, O''Neil %.2f mm from the apex.', z_sim, z_an));
            tc.verifyEqual(p_sim, p_an_max, 'RelTol', tol.peak_rel, ...
                sprintf('Focal pressure: simulated %.0f Pa, O''Neil %.0f Pa.', p_sim, p_an_max));

            tc.verifyEqual(L_sim, L_an, 'RelTol', tol.length_rel, ...
                sprintf('-3 dB focal length: simulated %.1f mm, O''Neil %.1f mm.', L_sim, L_an));
        end

    end

end

function p = oneil_rings_axial(R, rings, u0, f, c, rho, z)
% On-axis O'Neil pressure of in-phase rings on one spherical bowl: each ring
% [ID; OD] [m] is the cap of its outer minus the cap of its inner aperture,
% P = sum(exp(-ik r_e(Ri)) - exp(-ik r_e(Ro))) / (1 - z/R), with
% r_e(a) = sqrt((z - h(a))^2 + a^2) and h(a) = R - sqrt(R^2 - a^2)
% (focusedBowlONeil for a single cap, ID = 0).
    k  = 2 * pi * f / c;
    h  = @(a) R - sqrt(R^2 - a.^2);
    re = @(a) sqrt((z - h(a)).^2 + a.^2);
    P = zeros(size(z));
    P_R = 0;                                        % limit at the centre of curvature
    for i = 1:size(rings, 2)
        ri = rings(1, i) / 2; ro = rings(2, i) / 2;
        P = P + exp(-1i * k * re(ri)) - exp(-1i * k * re(ro));
        P_R = P_R + k * (h(ro) - h(ri));
    end
    P = P ./ (1 - z / R);
    P(abs(z - R) < eps) = P_R;
    p = rho * c * u0 * abs(P);
end

function [z_peak, p_peak] = refine_peak(z, p)
% Sub-sample peak by a parabola through the maximum and its neighbours.
    [p_peak, i] = max(p);
    z_peak = z(i);
    if i > 1 && i < numel(p)
        y = p(i-1:i+1);
        denom = y(1) - 2 * y(2) + y(3);
        if denom < 0
            d = 0.5 * (y(1) - y(3)) / denom;
            z_peak = z(i) + d * (z(2) - z(1));
            p_peak = y(2) - 0.25 * (y(1) - y(3)) * d;
        end
    end
end

function L = focal_length_3db(z, p, z_min)
% Axial extent around the far-field peak where intensity >= half its maximum
% (pressure >= max / sqrt(2)), with linear interpolation of the crossings.
    sel = z > z_min;
    z = z(sel); p = p(sel);
    [pk, i] = max(p);
    thr = pk / sqrt(2);
    lo = i;
    while lo > 1 && p(lo - 1) >= thr; lo = lo - 1; end
    hi = i;
    while hi < numel(p) && p(hi + 1) >= thr; hi = hi + 1; end
    z_lo = z(lo);
    if lo > 1
        z_lo = interp1(p([lo - 1, lo]), z([lo - 1, lo]), thr);
    end
    z_hi = z(hi);
    if hi < numel(p)
        z_hi = interp1(p([hi + 1, hi]), z([hi + 1, hi]), thr);
    end
    L = z_hi - z_lo;
end
