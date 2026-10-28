classdef test_config_tutorial < matlab.unittest.TestCase
% TEST_CONFIG_TUTORIAL  config_tutorial.yaml loads through LOAD_PARAMETERS with its positions (audit B4).
%
% The tutorial config once nested trans_pos/focus_pos under
% transducer.position, a path nothing reads, so the positions were silently
% lost. test_audit_fixes checks the raw YAML; this test checks the merged
% parameters that SIMPLE_MAIN and the pipeline actually receive.
%
%   results = runtests('tests/test_config_tutorial.m');

    properties
        P   % parameters loaded from config_tutorial.yaml
    end

    methods (TestClassSetup)
        function load_tutorial(tc)
            repo = fileparts(fileparts(mfilename('fullpath')));
            addpath(genpath(fullfile(repo, 'functions')));
            addpath(genpath(fullfile(repo, 'external')));
            addpath(fullfile(repo, 'config'));
            addpath(fullfile(repo, 'config', 'equipment'));
            % The tutorial's anat/SimNIBS paths are placeholders: missing-path
            % warnings are expected here.
            tc.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture( ...
                {'prestus:toolboxDefault', 'prestus:noSimNIBS', 'prestus:missingCharm'}));
            ws = warning('off', 'all');
            restore = onCleanup(@() warning(ws));
            evalc(['tc.P = load_parameters(''config_tutorial.yaml'', fullfile(repo, ''config''), ' ...
                '{''simulation.interactive'', 0; ''simulation.debug'', 0});']);
            clear restore
        end
    end

    methods (Test, TestTags = {'config_tutorial'})

        function test_positions_reach_parameters(tc)
            tr = tc.P.transducer;
            tc.verifyEqual(double(tr.trans_pos(:)'), [128 139 15]);
            tc.verifyEqual(double(tr.focus_pos(:)'), [128 142 75]);
            tc.verifyEqual(double(tr.focal_distance_ep), 60);
        end

        function test_no_dead_position_block(tc)
            tc.verifyFalse(isfield(tc.P.transducer, 'position'));
        end

        function test_transducer_settings_merged(tc)
            tr = tc.P.transducer;
            tc.verifyEqual(char(tr.type), 'matrix');
            tc.verifyEqual(double(tr.freq_hz), 250e3);
            tc.verifyEqual(char(tr.matrix.matrix_shape.type), 'define_here');
            tc.verifyEqual(double(tr.matrix.matrix_shape.define_here.grid_shape.fibonacci.elem_n), 52);
        end

        function test_single_transducer(tc)
            tc.verifyNumElements(tc.P.transducer, 1);
        end

        function test_timing_and_grid_merged(tc)
            t = tc.P.timing;
            tc.verifyEqual([t.pd t.pri t.ptd t.ptri t.ptrd], [0.3 0.6 0.6 1.6 14.4], 'AbsTol', 1e-12);
            tc.verifyEqual(double(tc.P.grid.default_dims(:)'), [144 144 256]);
        end

        function test_thermal_protocol_is_valid(tc)
            % The tutorial timing must pass THERMAL_PARAMETERS' integer checks
            pt = thermal_parameters(tc.P, true);
            tc.verifyEqual(pt.n_pulses_per_pt, 1);
            tc.verifyEqual(pt.n_ptri_reps, 9);
        end

    end

end
