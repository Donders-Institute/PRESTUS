classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_predict_outputs < matlab.unittest.TestCase
% TEST_PREDICT_OUTPUTS  Single-mode artifacts of PREDICT_OUTPUTS: the placement
% report and the time-reversal tables (no simulation).
%
% The Nexus mirrors these cases in nexus/tests/test_predict_outputs.py
% (test_single_mode_placement_report, test_single_mode_time_reversal_tables);
% if one side changes, the other should too.
%
%   Run with:   results = runtests('tests/test_predict_outputs.m');
%
% See also: PREDICT_OUTPUTS, PREPROC_TRANSDUCER_PLACEMENT, TIME_REVERSAL_STAGE

    properties (TestParameter)
        mode   = {'manual', 'localite', 'heuristic', 'mni', 'plantus', 'tus_entry', 'MNI'}
        stage  = {0, 1}
        report = {0, 1, []}
        medium = {'layered', 'water'}
        backend = {'kwave', 'jwave', 'babelbrain'}
    end

    methods (Test)

        function placementReport(tc, mode, stage, report)
            % Written with the placement stage on and placement.report on or
            % empty; the mode is lower-cased
            p = base_parameters();
            p.placement.mode = mode;
            p.modules.run_transducer_placement = stage;
            p.placement.report = report;
            names = {predict_outputs(p).artifacts.name};
            expected = sprintf('sub-003_placement_%s_desc-cond.html', lower(mode));
            on = stage == 1 && (isempty(report) || report == 1);
            tc.verifyEqual(any(strcmp(names, expected)), on);
            tc.verifyEqual(nnz(contains(names, '_placement_')), double(on));
        end

        function timeReversalTables(tc, medium, backend)
            % time_reversal_stage writes two tables; BabelBrain skips the stage
            for enabled = [false true]
                p = base_parameters();
                p.simulation.medium = medium;
                p.time_reversal.enabled = enabled;
                p.acoustic_backend = backend;
                names = {predict_outputs(p).artifacts.name};
                on = enabled && ~strcmp(backend, 'babelbrain');
                for t = {'summary', 'elements'}
                    f = sprintf('sub-003_%s_desc-cond_timereversal_%s.csv', medium, t{1});
                    tc.verifyEqual(any(strcmp(names, f)), on, f);
                end
            end
        end

        function absentFieldsUseDefaults(tc)
            % Without placement / time_reversal sections: report on, no tables
            p = base_parameters();
            p = rmfield(p, {'placement', 'time_reversal'});
            p.modules = rmfield(p.modules, 'run_transducer_placement');
            names = {predict_outputs(p).artifacts.name};
            tc.verifyTrue(any(strcmp(names, 'sub-003_placement_manual_desc-cond.html')));
            tc.verifyFalse(any(contains(names, 'timereversal')));
        end

    end
end

function p = base_parameters()
    p = struct();
    p.subject_id = 3;
    p.path.sim = tempdir;
    p.platform = 'matlab';
    p.simulation = struct('medium', 'layered', 'uncertainty', 0, 'debug', 0, ...
        'transducer_coupling', 'coherent');
    p.transducer = struct('target_isppa_wcm2', 30);
    p.io.output_affix = '_desc-cond';
    p.modules = struct('run_transducer_placement', 1, 'run_targeting_analysis', 1);
    p.placement = struct('mode', 'manual', 'report', true);
    p.time_reversal = struct('enabled', false);
    p.acoustic_backend = 'kwave';
end
