classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_fill_annular_drive < matlab.unittest.TestCase
% TEST_FILL_ANNULAR_DRIVE  A custom annular transducer is driven by its focal depth
% setting and target ISPPA; given phases and amplitudes (calibration) are kept.
% Pure logic: no k-Wave, no data.
%
%   results = runtests('tests/test_fill_annular_drive.m');

    methods (Test)
        function test_phases_from_the_depth_setting(tc)
            tr = custom_annular();
            tr.focal_distance_ep = 40;
            out = fill_annular_drive(tr, 1500, 1);
            tran = struct('transducer', struct('freq_hz', tr.freq_hz, 'annular', out.annular));
            tc.verifyEqual(out.annular.elem_phase_deg, compute_phases(1500, tran, 40), 'AbsTol', 1e-9);
            tc.verifyEqual(out.annular.elem_phase_rad, deg2rad(out.annular.elem_phase_deg), 'AbsTol', 1e-12);
        end

        function test_bowl_distance_gives_the_same_phases(tc)
            tr = custom_annular();
            tr.annular.dist_geom_ep_mm = 56.9;
            tr.focal_distance_ep = 40;
            a = fill_annular_drive(tr, 1500, 1);
            tr = rmfield(tr, 'focal_distance_ep');
            tr.focal_distance_bowl = 40 + (62.9 - 56.9);
            b = fill_annular_drive(tr, 1500, 1);
            tc.verifyEqual(b.annular.elem_phase_deg, a.annular.elem_phase_deg, 'AbsTol', 1e-9);
        end

        function test_given_phases_and_amplitude_are_kept(tc)
            tr = custom_annular();
            tr.focal_distance_ep = 40;
            tr.annular.elem_phase_deg = [0 123];
            tr.annular.elem_p_amp = [2e5 3e5];
            out = fill_annular_drive(tr, 1500, 1);
            tc.verifyEqual(out.annular.elem_phase_deg, [0 123]);
            tc.verifyEqual(out.annular.elem_p_amp, [2e5 3e5]);
        end

        function test_nominal_amplitude_only_with_a_target_isppa(tc)
            tr = custom_annular();
            tr.target_isppa_wcm2 = 20;
            out = fill_annular_drive(tr, 1500, 1);
            tc.verifySize(out.annular.elem_p_amp, [1 2]);
            tr.target_isppa_wcm2 = [];
            out = fill_annular_drive(tr, 1500, 1);
            tc.verifyEmpty(out.annular.elem_p_amp);           % the amplitude sets the scale
        end

        function test_without_depth_nothing_is_steered(tc)
            out = fill_annular_drive(custom_annular(), 1500, 1);
            tc.verifyEmpty(out.annular.elem_phase_deg);
        end

        function test_flat_rings_are_left_alone(tc)
            tr = custom_annular();
            tr.annular.curv_radius_mm = Inf;
            tr.focal_distance_ep = 40;
            out = fill_annular_drive(tr, 1500, 1);
            tc.verifyEmpty(out.annular.elem_phase_deg);
        end
    end
end

function tr = custom_annular()
    tr = struct('type', 'annular', 'freq_hz', 500e3, 'target_isppa_wcm2', [], ...
                'focal_distance_ep', [], 'focal_distance_bowl', []);
    tr.annular = struct('elem_n', 2, 'elem_id_mm', [0 30], 'elem_od_mm', [29 42], ...
                        'curv_radius_mm', 62.9, 'elem_phase_deg', [], 'elem_p_amp', []);
end
