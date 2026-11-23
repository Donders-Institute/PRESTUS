classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_babelbrain_match_grid < matlab.unittest.TestCase
% TEST_BABELBRAIN_MATCH_GRID  babelbrain.match_grid and BabelBrain's spatial step.
%
% BabelBrain simulates at SmallestSOS / freq / PPW, SmallestSOS floored by its
% Pichardo shear speed (1102.515 m/s). With babelbrain.match_grid the PRESTUS
% grid takes that step, so write_babelbrain_input does not resample the medium.
%
%   Run with:   results = runtests('tests/test_babelbrain_match_grid.m');

    methods (Test, TestTags = {'unit', 'quick', 'nogpu'})

        function test_step_formula(tc)
            % Pichardo floor binds for the usual media; a slower medium takes over
            p = tc.base();
            for f = [250e3 500e3 750e3]
                for ppw = [6 9 12]
                    tc.verifyEqual(babelbrain_grid_step_mm(p, f, ppw), 1102.515 / f / ppw * 1e3, ...
                        'RelTol', 1e-12, sprintf('%g Hz, %d PPW', f, ppw));
                end
            end
            tc.verifyEqual(babelbrain_grid_step_mm(p, 500e3, 6), 0.367505, 'AbsTol', 1e-6);
            p.medium_properties.air_like = struct('sound_speed', 900);
            [s, sos] = babelbrain_grid_step_mm(p, 500e3, 6);
            tc.verifyEqual(sos, 900);
            tc.verifyEqual(s, 900 / 500e3 / 6 * 1e3, 'RelTol', 1e-12);
            [s0, sos0] = babelbrain_grid_step_mm(tc.base());
            tc.verifyTrue(isnan(s0));
            tc.verifyEqual(sos0, 1102.515, 'AbsTol', 1e-9);
        end

        function test_match_grid_sets_resolution(tc)
            % every frequency and PPW: the grid takes BabelBrain's step, and
            % write_babelbrain_input's PPW fallback lands on the same PPW (no resampling)
            for f = [250e3 500e3 750e3]
                for ppw = [6 9 12]
                    p = tc.base(); p.transducer.freq_hz = f; p.babelbrain.ppw = ppw;
                    p.babelbrain.match_grid = true;
                    q = babelbrain_match_grid(p);
                    tc.verifyEqual(q.grid.resolution_mm, 1102.515 / f / ppw * 1e3, 'RelTol', 1e-12);
                    tc.verifyEqual(round(1102.515 / (f * q.grid.resolution_mm * 1e-3)), ppw);
                    % idempotent
                    tc.verifyEqual(babelbrain_match_grid(q).grid.resolution_mm, q.grid.resolution_mm);
                end
            end
        end

        function test_off_leaves_grid(tc)
            % negative control: default (false) and absent key keep grid.resolution_mm
            p = tc.base(); p.babelbrain.match_grid = false;
            tc.verifyEqual(babelbrain_match_grid(p).grid.resolution_mm, 0.5);
            p.babelbrain = rmfield(p.babelbrain, 'match_grid');
            tc.verifyEqual(babelbrain_match_grid(p).grid.resolution_mm, 0.5);
            p.acoustic_backend = 'kwave';   % off: no backend check
            tc.verifyEqual(babelbrain_match_grid(p).grid.resolution_mm, 0.5);
        end

        function test_needs_babelbrain_backend(tc)
            for b = {'kwave', 'jwave', ''}
                p = tc.base(); p.babelbrain.match_grid = true; p.acoustic_backend = b{1};
                tc.verifyError(@() babelbrain_match_grid(p), 'babelbrain_match_grid:backend');
            end
        end
    end

    methods (Static)
        function p = base()
            p = struct();
            p.acoustic_backend = 'babelbrain';
            p.grid.resolution_mm = 0.5;
            p.babelbrain.ppw = 6;
            p.babelbrain.match_grid = false;
            p.transducer = struct('freq_hz', 500e3);
            p.medium_properties = struct( ...
                'water', struct('sound_speed', 1500), ...
                'skull', struct('sound_speed', 2800), ...
                'brain', struct('sound_speed', 1546));
        end
    end
end
