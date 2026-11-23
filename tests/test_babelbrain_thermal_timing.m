classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_babelbrain_thermal_timing < matlab.unittest.TestCase
% TEST_BABELBRAIN_THERMAL_TIMING  Unit tests for the PRESTUS -> BabelBrain thermal timing map.
%
% BabelBrain runs Repetitions x (DurationUS on, DurationOff off) and then, with
% a single sonication group, one PauseBetweenGroupedSonications off period.
% The mapping must reproduce the native PRESTUS protocol, including the
% post_ptri_dur cooldown after ptrd.
%
%   Run with:   results = runtests('tests/test_babelbrain_thermal_timing.m');

    methods (Test, TestTags = {'unit', 'quick', 'nogpu'})

        function test_two_train_protocol_with_cooldown(tc)
            t = tc.base_timing();
            bt = babelbrain_thermal_timing(t);
            tc.verifyEqual(bt.duty_cycle,     0.3,  'AbsTol', 1e-12);
            tc.verifyEqual(bt.prf_hz,         5,    'AbsTol', 1e-12);
            tc.verifyEqual(bt.duration_s,     10);
            tc.verifyEqual(bt.duration_off_s, 20);
            tc.verifyEqual(bt.repetitions,    2);
            tc.verifyEqual(bt.pause_after_s,  60);
            % BabelBrain total simulated time equals the native protocol span.
            total = bt.repetitions * (bt.duration_s + bt.duration_off_s) + bt.pause_after_s;
            tc.verifyEqual(total, t.ptrd + t.post_ptri_dur, 'AbsTol', 1e-12);
        end

        function test_repetitions_from_config_timing(tc)
            % parameters.timing as configured has ptrd and ptri but no
            % n_ptri_reps (only THERMAL_PARAMETERS derives it): the
            % repetitions come from ptrd / ptri, for every protocol shape
            shapes = {struct('ptd', 10, 'ptri', 30, 'ptrd', 60), ...   % 2 trains
                      struct('ptd', 40, 'ptri', 50, 'ptrd', 150), ...  % 3 trains
                      struct('ptd', 10, 'ptri', 10, 'ptrd', 10)};      % one continuous train
            for i = 1:numel(shapes)
                t = shapes{i}; t.pd = 0.06; t.pri = 0.2; t.post_ptri_dur = 0;
                bt = babelbrain_thermal_timing(t);
                tc.verifyEqual(bt.repetitions, round(t.ptrd / t.ptri), sprintf('shape %d', i));
                % the derived field, when present, gives the same answer
                t.n_ptri_reps = t.ptrd / t.ptri;
                tc.verifyEqual(babelbrain_thermal_timing(t).repetitions, bt.repetitions, sprintf('shape %d', i));
            end
        end

        function test_no_cooldown_when_unset_or_nonpositive(tc)
            t = rmfield(tc.base_timing(), 'post_ptri_dur');
            tc.verifyEqual(babelbrain_thermal_timing(t).pause_after_s, 0);
            for v = {0, NaN, [], -5}
                t.post_ptri_dur = v{1};
                tc.verifyEqual(babelbrain_thermal_timing(t).pause_after_s, 0);
            end
        end

        function test_single_continuous_train(tc)
            t = tc.base_timing();
            t.ptri = t.ptd;  t.ptrd = t.ptd;  t.n_ptri_reps = 1;
            bt = babelbrain_thermal_timing(t);
            tc.verifyEqual(bt.duration_off_s, 0);
            tc.verifyEqual(bt.repetitions,    1);
            tc.verifyEqual(bt.pause_after_s,  60);
        end

        function test_invalid_timing_errors(tc)
            t = rmfield(tc.base_timing(), 'ptd');
            tc.verifyError(@() babelbrain_thermal_timing(t), 'babelbrain_thermal_timing:missingPtd');
            t = tc.base_timing();  t.pd = 2 * t.pri;
            tc.verifyError(@() babelbrain_thermal_timing(t), 'babelbrain_thermal_timing:invalidDutyCycle');
        end

        function test_exported_yaml_matches_timing_map(tc)
            % export_babelbrain_inputs writes the thermal_inputs block from
            % babelbrain_thermal_timing; one protocol with trains and
            % cooldown, one single continuous train.
            prestus_root = fileparts(fileparts(mfilename('fullpath')));
            tc.assumeTrue(isfolder(fullfile(prestus_root, 'external', 'BabelBrain', 'TranscranialModeling')), ...
                'BabelBrain submodule not present.');
            t1 = tc.base_timing();
            t2 = t1;  t2.ptri = t2.ptd;  t2.ptrd = t2.ptd;  t2.n_ptri_reps = 1;  t2.post_ptri_dur = NaN;
            for t = {t1, t2}
                tmp = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
                p = tc.export_params(tmp, t{1});
                masks = ones(24, 24, 40);  masks(:, :, 20:end) = 2;
                out = evalc_free(@() export_babelbrain_inputs(p, masks, [12 12 30], ...
                    struct('origin_ras_mm', [0 0 0])));
                txt = fileread(out.config_yaml);
                bt = babelbrain_thermal_timing(t{1});
                tc.verifyEqual(yaml_num(txt, 'DutyCycle'),   bt.duty_cycle,     'RelTol', 1e-5);
                tc.verifyEqual(yaml_num(txt, 'PRF'),         bt.prf_hz,         'RelTol', 1e-5);
                tc.verifyEqual(yaml_num(txt, 'DurationUS'),  bt.duration_s,     'RelTol', 1e-5);
                tc.verifyEqual(yaml_num(txt, 'DurationOff'), bt.duration_off_s, 'AbsTol', 1e-9);
                tc.verifyEqual(yaml_num(txt, 'Repetitions'), bt.repetitions);
                tc.verifyEqual(yaml_num(txt, 'PauseBetweenGroupedSonications'), bt.pause_after_s, 'AbsTol', 1e-9);
            end
            % Negative control: DurationUS is the train ON time, not the total span.
            tc.verifyNotEqual(t1.ptd, t1.ptrd);
        end

    end

    methods (Static)
        function p = export_params(out_dir, timing)
            p = struct();
            p.subject_id = 1;
            p.io.output_affix = '_tt';
            p.io.dir_babelbrain = out_dir;
            p.grid.resolution_mm = 1;
            p.medium_properties.water.sound_speed = 1500;
            p.medium_properties.brain.sound_speed = 1540;
            p.babelbrain.computing_backend = 4;
            p.babelbrain.device_name = 'CPU';
            p.transducer(1).freq_hz = 500e3;
            p.transducer(1).trans_pos = [12 12 2];
            p.transducer(1).annular = struct('elem_id_mm', 0, 'elem_od_mm', 64, ...
                'curv_radius_mm', 63.2, 'elem_n', 1);
            p.timing = timing;
        end

        function t = base_timing()
            t = struct('pd', 0.06, 'pri', 0.2, 'ptd', 10, 'ptri', 30, 'ptrd', 60, ...
                       'n_ptri_reps', 2, 'post_ptri_dur', 60);
        end
    end
end

function out = evalc_free(fn)
% Run fn with its console output suppressed.
out = [];
evalc('out = fn();');
end

function v = yaml_num(txt, key)
tok = regexp(txt, ['\n\s*' key ':\s*([-+0-9.eEnaNIinf]+)'], 'tokens', 'once');
assert(~isempty(tok), 'Key %s not found in exported YAML.', key);
v = str2double(tok{1});
end
