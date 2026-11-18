classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_acoustic_store_key < matlab.unittest.TestCase
% TEST_ACOUSTIC_STORE_KEY  Which transducer entries of an async run share an
% acoustic field (simulated once, and kept for later runs), and STRUCT_MD5.
% Pure logic: no k-Wave, no data.
%
%   results = runtests('tests/test_acoustic_store_key.m');

    methods (Static, Access = private)
        function p = params()
            % A sequence: T1, then T2 elsewhere, then T1's stance again (T3)
            p = struct();
            p.subject_id = 1;
            p.simulation = struct('medium', 'water', 'code_type', 'matlab_gpu', 'precision', 'single', ...
                                  'transducer_coupling', 'async');
            p.grid = struct('resolution_mm', 0.5, 'async_layout', 'shared');
            p.layers = struct('skull', 4);
            p.medium_properties = struct('skull', struct('sound_speed', 2800));
            p.placement = struct('mode', 'manual', 'manual', struct('coordinate_space', 'mni'));
            p.modules = struct('run_transducer_placement', 1);
            p.timing = struct('pd', 0.02, 'async_schedule', 'custom');
            p.transducer = struct( ...
                'serial',            {'CTX_500_026', 'CTX_500_026', 'CTX_500_026'}, ...
                'trans_pos',         {[-52 -8 58], [52 -8 58], [-52 -8 58]}, ...
                'focus_pos',         {[-29 -6 7], [29 -6 7], [-29 -6 7]}, ...
                'freq_hz',           {5e5, 5e5, 5e5}, ...
                'target_isppa_wcm2', {30, 30, 20}, ...
                'timing_offset_s',   {0, 50, 100}, ...
                'timing',            {struct(), struct(), struct('ptrd', 40)}, ...
                'name',              {'stance 1', 'stance 2', 'stance 3'});
            p.path = struct('sim', '/a');
            p.io = struct('overwrite_files', 'never');
        end

        function [key, share] = keys(p)
            n = numel(p.transducer);
            key = cell(1, n); share = false(1, n);
            for ti = 1:n
                pa = p;                                % as make_acoustic_params
                pa.transducer = p.transducer(ti);
                pa.transducer.target_isppa_wcm2 = NaN;
                [key{ti}, share(ti)] = acoustic_store_key(p, ti, pa);
            end
        end
    end

    methods (Test)

        function test_returning_stance_shares_the_field(tc)
            % Same stance and settings: same key, whatever the timing, target and name
            [key, share] = tc.keys(tc.params());
            tc.verifyEqual(key{1}, key{3});
            tc.verifyNotEqual(key{1}, key{2});
            tc.verifyEqual(share, [true true true]);
            tc.verifyEqual(strlength(string(key{1})), 16);
        end

        function test_changed_setting_changes_the_key(tc)
            % Negative controls: anything the field depends on separates the entries
            base = tc.params();
            changes = {
                @(p) setfield(p, 'transducer', {3}, 'focus_pos', [-29 -6 9])             % stance
                @(p) setfield(p, 'transducer', {3}, 'freq_hz', 2.5e5)                    % frequency
                @(p) setfield(p, 'medium_properties', 'skull', 'sound_speed', 3000)      % medium (uncertainty variant)
                @(p) setfield(p, 'grid', 'resolution_mm', 1)                             % grid
                @(p) setfield(p, 'acoustic_backend', 'jwave')                            % backend
                };
            [k0, ~] = tc.keys(base);
            for i = 1:numel(changes)
                [k, ~] = tc.keys(changes{i}(base));
                tc.verifyNotEqual(k{3}, k0{3}, func2str(changes{i}));
            end
        end

        function test_shared_grid_depends_on_every_transducer(tc)
            % The shared grid is sized for all transducers: moving T2 changes T1's field
            p = tc.params();
            [k0, ~] = tc.keys(p);
            q = p; q.transducer(2).trans_pos = [60 -8 50];
            [k1, ~] = tc.keys(q);
            tc.verifyNotEqual(k1{1}, k0{1});
            % ... but not on its own grid
            p.grid.async_layout = 'per_transducer';
            q.grid.async_layout = 'per_transducer';
            [k0, ~] = tc.keys(p);
            [k1, ~] = tc.keys(q);
            tc.verifyEqual(k1{1}, k0{1});
            tc.verifyEqual(k1{1}, k1{3});
        end

        function test_resolved_placement_is_never_shared(tc)
            % Positions resolved per index at run time: every entry its own field
            for mode = {'mni', 'heuristic', 'localite', 'plantus', 'tus_entry'}
                p = tc.params();
                p.placement.mode = mode{1};
                [key, share] = tc.keys(p);
                tc.verifyNotEqual(key{1}, key{3}, mode{1});
                tc.verifyFalse(any(share), mode{1});
            end
            % without a placement stage the entries are explicit again
            p.modules.run_transducer_placement = 0;
            [key, share] = tc.keys(p);
            tc.verifyEqual(key{1}, key{3});
            tc.verifyTrue(all(share));
        end

        function test_layered_without_segmentation_is_not_kept(tc)
            % The segmentation's content is part of the key: unknown before charm ran
            p = tc.params();
            p.simulation.medium = 'layered';
            [key, share] = tc.keys(p);
            tc.verifyEqual(key{1}, key{3});        % one run: still simulated once
            tc.verifyFalse(any(share));            % but not kept for later runs
        end

        function test_struct_md5_ignores_field_order(tc)
            a = struct('x', 1, 'y', struct('b', [1 2], 'a', 'z'));
            b = struct('y', struct('a', 'z', 'b', [1 2]), 'x', 1);
            tc.verifyEqual(struct_md5(a), struct_md5(b));
            b.x = 2;
            tc.verifyNotEqual(struct_md5(a), struct_md5(b));
            tc.verifyEqual(strlength(string(struct_md5(a))), 32);
        end

    end

end
