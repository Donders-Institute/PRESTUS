classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_time_reversal < matlab.unittest.TestCase
% TEST_TIME_REVERSAL  Unit tests for the time-reversal / reciprocity functions
%
% Synthetic fields only (no k-Wave run): free-space Green's functions over the
% surface of a bowl, with and without a phase screen in front of it. Checks
%   ANNULAR_APERTURE_POINTS  ring areas
%   PREDICT_FOCAL_SHIFT      no shift in water; the shift a phase screen
%                            causes, against a brute-force forward sum
%   TIME_REVERSAL_REFOCUS    conjugate phases undo the screen; in water the
%                            ring phases equal geometric steering
%   TIME_REVERSAL_POSE_PARAMETERS  placement, focal depth and phases of a
%                            pose (geometric, library, single element)
%
%   Run with:   results = runtests('tests/test_time_reversal.m');

    properties
        Geom   % annular geometry
        Apex   % [m]
        U      % beam direction
        K      % wavenumber [rad/m]
        Pts; El; Area
    end

    methods (TestMethodSetup)
        function make_bowl(tc)
            tc.Geom = struct('curv_radius_mm', 64, ...
                'elem_id_mm', [0 20 32 41], 'elem_od_mm', [19 31 40 48]);
            tc.Apex = [0 0 0];
            tc.U    = [0 0 1];
            tc.K    = 2*pi*500e3 / 1500;
            [tc.Pts, tc.El, tc.Area] = annular_aperture_points(tc.Geom, tc.Apex, tc.U, 0.75e-3);
        end
    end

    methods (Test, TestTags = {'time_reversal'})

        function test_ring_areas(tc)
            R = tc.Geom.curv_radius_mm * 1e-3;
            for e = 1:numel(tc.Geom.elem_od_mm)
                th0 = asin(tc.Geom.elem_id_mm(e) / 2e3 / R);
                th1 = asin(tc.Geom.elem_od_mm(e) / 2e3 / R);
                expected = 2*pi*R^2 * (cos(th0) - cos(th1));
                tc.verifyEqual(sum(tc.Area(tc.El == e)), expected, 'RelTol', 1e-2);
            end
            % surface on the sphere around the geometric focus
            r = sqrt(sum((tc.Pts - [0 0 R]).^2, 2));
            tc.verifyEqual(r, repmat(R, size(r)), 'AbsTol', 1e-9);
        end

        function test_water_no_shift(tc)
            t = [0 0 0.060];                                   % steered 4 mm short of the natural focus
            [q, ~] = forward_terms(tc, t, zeros(size(tc.El)));
            pred = predict_focal_shift(tc.Pts, q, t, tc.U, tc.K);
            tc.verifyLessThan(pred.shift_mm, 0.3);
            tc.verifyEqual(abs(sum(q)) / sum(abs(q)), 1, 'AbsTol', 0.02);
        end

        function test_phase_screen_shift_matches_forward(tc)
            t = [0 0 0.064];
            psi = 100 * tc.Pts(:, 1);                          % linear phase screen [rad], across x
            [q, c] = forward_terms(tc, t, psi);                % expected shift ~ 100 * F / k = 3 mm
            pred = predict_focal_shift(tc.Pts, q, t, tc.U, tc.K);
            % brute force: forward field of the same drive through the same screen
            ax = -9e-3:0.1e-3:9e-3;
            Y = t + [ax(:), zeros(numel(ax), 1), zeros(numel(ax), 1)];
            F = zeros(numel(ax), 1);
            for i = 1:numel(ax)
                r = sqrt(sum((tc.Pts - Y(i, :)).^2, 2));
                F(i) = sum(c .* exp(1i * psi) .* exp(-1i * tc.K * r) ./ r);
            end
            [~, im] = max(abs(F));
            tc.verifyGreaterThan(abs(ax(im)), 1e-3, 'screen should move the focus by > 1 mm');
            tc.verifyEqual(pred.shift_m(1), ax(im), 'AbsTol', 0.3e-3);
            tc.verifyEqual(pred.shift_lateral_mm, abs(ax(im)) * 1e3, 'AbsTol', 0.3);
        end

        function test_refocus_undoes_ring_screen(tc)
            t = [0 0 0.064];
            n_el = numel(tc.Geom.elem_od_mm);
            ring_err = [0; 1.2; -0.8; 2.0];                    % per-ring phase error [rad]
            [~, ~, P] = forward_terms(tc, t, ring_err(tc.El));
            W = sparse(1:numel(tc.El), tc.El, tc.Area, numel(tc.El), n_el);
            C = ones(n_el, 1) .* exp(1i * geometric_phases(tc, t));
            res = time_reversal_refocus(P, W, C, tc.Pts, t, tc.U, tc.K, 'annular');
            tc.verifyEqual(res.capability, 'rings');
            tc.verifyLessThan(res.coherence_elements, 0.9);
            tc.verifyEqual(res.refocus_gain, 1 / res.coherence_elements, 'RelTol', 1e-10);
            % refocused drive arrives in phase
            C_tr = abs(C) .* exp(1i * res.phase_tr_rad);
            arrive = C_tr .* res.G;
            tc.verifyEqual(abs(sum(arrive)) / sum(abs(arrive)), 1, 'AbsTol', 1e-10);
            % the correction is the screen (relative to ring 1)
            d = angle(exp(1i * (res.phase_tr_rad - angle(C))));
            % (up to the difference between the mean path and the phase of the
            % ring sum, small for these ring widths)
            tc.verifyEqual(angle(exp(1i * (d - d(1)))), angle(exp(-1i * (ring_err - ring_err(1)))), 'AbsTol', 0.05);
        end

        function test_water_ring_phases_are_steering(tc)
            t = [0 0 0.074];                                   % steered 10 mm beyond the natural focus
            n_el = numel(tc.Geom.elem_od_mm);
            [~, ~, P] = forward_terms(tc, t, zeros(size(tc.El)));
            W = sparse(1:numel(tc.El), tc.El, tc.Area, numel(tc.El), n_el);
            C = ones(n_el, 1);                                 % unsteered
            res = time_reversal_refocus(P, W, C, tc.Pts, t, tc.U, tc.K, 'annular', 'PredictShift', false);
            geo = geometric_phases(tc, t);
            d = angle(exp(1i * (res.phase_tr_rad - geo)));
            tc.verifyLessThan(max(abs(angle(exp(1i * (d - d(1)))))), deg2rad(10));
        end

        function test_single_element_capability(tc)
            t = [0 0 0.064];
            [~, ~, P] = forward_terms(tc, t, zeros(size(tc.El)));
            W = sparse(1:numel(tc.El), 1, tc.Area, numel(tc.El), 1);
            res = time_reversal_refocus(P, W, 1, tc.Pts, t, tc.U, tc.K, 'annular', 'PredictShift', false);
            tc.verifyEqual(res.capability, 'none');
            tc.verifyEqual(res.refocus_gain, 1, 'AbsTol', 1e-12);
        end

        function test_pose_parameters_geometric(tc)
            [p, pose] = pose_config(tc);
            [q, info] = time_reversal_pose_parameters(p, pose, 'Aim', 'corrected', 'Affix', '_pose007_aim');
            apex = [10 20 30]; aim = [12 21 95];
            tc.verifyEqual(q.transducer.trans_pos, apex);
            tc.verifyEqual(q.transducer.focus_pos, aim);
            tc.verifyEqual(q.transducer.focal_distance_ep, norm(aim - apex) - (64 - 58), 'AbsTol', 1e-12);
            tc.verifyEmpty(q.transducer.focal_distance_bowl);
            tc.verifyEqual(info.phase_source, 'geometric');
            expected = annular_geometric_phases(tc.Geom, apex, aim, 500e3, 1500);
            tc.verifyEqual(q.transducer.annular.elem_phase_deg, expected, 'AbsTol', 1e-9);
            tc.verifyEqual(q.transducer.annular.elem_phase_rad, deg2rad(expected), 'AbsTol', 1e-12);
            tc.verifyEqual(q.placement.mode, 'manual');
            tc.verifyEqual(q.placement.manual.coordinate_space, 'ras');
            tc.verifyEqual(q.io.output_affix, '_base_pose007_aim');
            tc.verifyFalse(q.time_reversal.enabled);
            tc.verifyFalse(q.time_reversal.pose_workflow.enabled);
            % geometric aim: steered at the target
            q = time_reversal_pose_parameters(p, pose, 'Aim', 'geometric');
            tc.verifyEqual(q.transducer.focus_pos, [10 20 94]);
        end

        function test_pose_parameters_library_and_single(tc)
            [p, pose] = pose_config(tc);
            p.state.library_drive = {struct('phases', true, 'amp', true)};
            [q, info] = time_reversal_pose_parameters(p, pose);
            tc.verifyEqual(info.phase_source, 'library');
            tc.verifyEmpty(q.transducer.annular.elem_phase_deg);   % resolved again for the new depth
            tc.verifyEmpty(q.transducer.annular.elem_p_amp);
            tc.verifyFalse(isfield(q.state, 'library_drive'));
            % single element: phases untouched
            [p, pose] = pose_config(tc);
            p.transducer.annular.elem_id_mm = 0;
            p.transducer.annular.elem_od_mm = 48;
            p.transducer.annular.elem_phase_deg = 0;
            [q, info] = time_reversal_pose_parameters(p, pose);
            tc.verifyEqual(info.phase_source, 'none');
            tc.verifyEqual(q.transducer.annular.elem_phase_deg, 0);
        end

        function test_pose_parameters_rejects_matrix(tc)
            [p, pose] = pose_config(tc);
            p.transducer.type = 'matrix';
            tc.verifyError(@() time_reversal_pose_parameters(p, pose), 'PRESTUS:posePlacement:transducerType');
        end
    end
end

function [p, pose] = pose_config(tc)
% Minimal config of one curved annular transducer and one pose-table row (RAS mm)
    tr = struct('type', 'annular', 'freq_hz', 500e3, 'trans_pos', [0 0 0], 'focus_pos', [0 0 60], ...
        'focal_distance_ep', 54, 'focal_distance_bowl', 60);
    tr.annular = tc.Geom;
    tr.annular.dist_geom_ep_mm = 58;
    tr.annular.elem_phase_deg  = [0 10 20 30];
    tr.annular.elem_phase_rad  = deg2rad([0 10 20 30]);
    tr.annular.elem_p_amp      = [1 1 1 1] * 1e5;
    p = struct('transducer', tr, 'subject_id', 1);
    p.medium_properties.water.sound_speed = 1500;
    p.placement.mode = 'plantus';
    p.modules.run_transducer_placement = 1;
    p.io.output_affix = '_base';
    p.time_reversal = struct('enabled', true, 'only', false, ...
        'placement_search', struct('enabled', true), 'pose_workflow', struct('enabled', true));
    pose = struct('pose', 7, 'apex_ras_x', 10, 'apex_ras_y', 20, 'apex_ras_z', 30, ...
        'aim_ras_x', 12, 'aim_ras_y', 21, 'aim_ras_z', 95, ...
        'target_ras_x', 10, 'target_ras_y', 20, 'target_ras_z', 94);
end

function phi = geometric_phases(tc, t)
% Ring phases that steer to t in water (area-weighted mean path)
    r = sqrt(sum((tc.Pts - t).^2, 2));
    n_el = max(tc.El);
    phi = tc.K * accumarray(tc.El, tc.Area .* r, [n_el 1]) ./ accumarray(tc.El, tc.Area, [n_el 1]);
end

function [q, c, P] = forward_terms(tc, t, psi)
% Reverse field of a point source at t through a phase screen psi at the
% surface (reciprocal: the same screen on the way out), the geometric drive
% c (area-weighted) and the contributions q = c .* P
    r = sqrt(sum((tc.Pts - t).^2, 2));
    P = exp(1i * psi) .* exp(-1i * tc.K * r) ./ r;
    phi = geometric_phases(tc, t);
    c = tc.Area .* exp(1i * phi(tc.El));
    q = c .* P;
end
