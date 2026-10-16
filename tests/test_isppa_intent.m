classdef test_isppa_intent < matlab.unittest.TestCase
% TEST_ISPPA_INTENT  The run states which free-water ISPPA it models.
%
% Covers ISPPA_INTENT (target vs. the calibration ISPPA of the drive) and the
% record LOAD_PARAMETERS keeps of a target left at the toolbox default.
% Regression for a 60 W/cm² calibrated drive that a config without
% target_isppa_wcm2 rescaled to the 30 W/cm² default without notice.
% Pure logic: no k-Wave, no data.
%
%   results = runtests('tests/test_isppa_intent.m');

    methods (TestClassSetup)
        function add_paths(~)
            repo = fileparts(fileparts(mfilename('fullpath')));
            addpath(genpath(fullfile(repo, 'functions')));
        end
    end

    methods (Static, Access = private)
        function p = params(target, cal, default)
            p = struct();
            p.transducer = struct('target_isppa_wcm2', target);
            if ~isempty(cal)
                p.transducer.tx_calibration = struct('isppa_wcm2', cal, 'isppa_source', 'library');
            end
            if nargin > 2
                p.state.target_isppa_default_wcm2 = default;
            end
        end
    end

    methods (Test, TestTags = {'isppa_policy'})

        function test_no_target_reports_drive_level(tc)
            [msg, bad] = isppa_intent(tc.params([], 60));
            tc.verifyFalse(bad);
            tc.verifySubstring(msg, 'no scaling');
            tc.verifySubstring(msg, 'calibrated for 60');
        end

        function test_target_matching_calibration(tc)
            [msg, bad] = isppa_intent(tc.params(60, 60));
            tc.verifyFalse(bad);
            tc.verifySubstring(msg, 'calibrated intensity');
        end

        function test_default_target_differs_from_calibration(tc)
            [msg, bad] = isppa_intent(tc.params(30, 60, 30));
            tc.verifyTrue(bad);
            tc.verifySubstring(msg, 'toolbox default');
            tc.verifySubstring(msg, 'calibrated for 60');
            tc.verifySubstring(msg, 'target_isppa_wcm2: 60');
        end

        function test_explicit_target_differs_from_calibration(tc)
            [msg, bad] = isppa_intent(tc.params(20, 60));
            tc.verifyTrue(bad);
            tc.verifyTrue(~contains(msg, 'toolbox default'));
        end

        function test_default_target_without_calibration_is_flagged(tc)
            [msg, bad] = isppa_intent(tc.params(30, [], 30));
            tc.verifyTrue(bad);
            tc.verifySubstring(msg, 'free-water baseline');
        end

        function test_explicit_target_without_calibration_not_flagged(tc)
            [~, bad] = isppa_intent(tc.params(30, []));
            tc.verifyFalse(bad);
        end

        function test_script_changed_default_is_not_default(tc)
            % a value set after load_parameters no longer counts as the default
            [msg, ~] = isppa_intent(tc.params(60, [], 30));
            tc.verifyTrue(~contains(msg, 'toolbox default'));
        end

        function test_load_parameters_records_default_only_when_unset(tc)
            w = warning('off', 'all');
            restore = onCleanup(@() warning(w));
            p = load_parameters(struct('simulation', struct('interactive', 0)));
            tc.verifyTrue(isfield(p, 'state') && isfield(p.state, 'target_isppa_default_wcm2'));
            p = load_parameters(struct('transducer', struct('target_isppa_wcm2', 60)));
            tc.verifyFalse(isfield(p, 'state') && isfield(p.state, 'target_isppa_default_wcm2'));
            p = load_parameters(struct('simulation', struct('interactive', 0)), ...
                {'transducer.target_isppa_wcm2', 45});
            tc.verifyFalse(isfield(p, 'state') && isfield(p.state, 'target_isppa_default_wcm2'));
        end
    end
end
