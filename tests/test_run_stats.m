classdef test_run_stats < matlab.unittest.TestCase
% TEST_RUN_STATS  Unit tests for functions/helper/run_stats
%
% run_stats collects what a run simulated and appends one JSON line to
% <path.sim>/run_stats.jsonl; nexus/prestus_nexus/estimate.py reads it.
%
%   Run with:   results = runtests('tests/test_run_stats.m');

    properties
        TmpDir
        Parameters
    end

    methods (TestMethodSetup)
        function setup(tc)
            tc.TmpDir = tempname;
            mkdir(tc.TmpDir);
            tc.Parameters.path.sim = tc.TmpDir;
            tc.Parameters.simulation.code_type = 'matlab_cpu';
            tc.Parameters.simulation.medium = 'water';
            tc.Parameters.grid.resolution_mm = 0.5;
            run_stats('reset');
        end
    end

    methods (TestMethodTeardown)
        function teardown(tc)
            run_stats('reset');
            if isfolder(tc.TmpDir), rmdir(tc.TmpDir, 's'); end
        end
    end

    methods (Test)
        function testNothingSimulatedWritesNothing(tc)
            run_stats('write', tc.Parameters, 5);
            tc.verifyFalse(isfile(tc.statsFile()));
        end

        function testRecordAppended(tc)
            run_stats('add', 'acoustic_s', 10);
            run_stats('add', 'acoustic_s', 2.5);                    % summed
            run_stats('set', 'acoustic_voxels', 1e6, 'acoustic_nt', 300);
            run_stats('add', 'thermal_s', 4, 'thermal_steps', 40);
            run_stats('write', tc.Parameters, 20);
            rec = jsondecode(fileread(tc.statsFile()));
            tc.verifyEqual(rec.acoustic_s, 12.5);
            tc.verifyEqual(rec.acoustic_voxels, 1e6);
            tc.verifyEqual(rec.thermal_steps, 40);
            tc.verifyEqual(rec.total_s, 20);
            tc.verifyEqual(rec.medium, 'water');
            tc.verifyEqual(rec.acoustic_backend, 'kwave');
            tc.verifyEqual(rec.resolution_mm, 0.5);
        end

        function testWriteClears(tc)
            % A nested run (post-hoc water, sequential follow-up) must not
            % inherit the outer run's figures.
            run_stats('add', 'acoustic_s', 10);
            run_stats('write', tc.Parameters, 20);
            tc.verifyEmpty(fieldnames(run_stats('get')));
            run_stats('add', 'acoustic_s', 3);
            run_stats('write', tc.Parameters, 4);
            lines = strsplit(strtrim(fileread(tc.statsFile())), newline);
            tc.verifyNumElements(lines, 2);
            tc.verifyEqual(jsondecode(lines{2}).acoustic_s, 3);
        end

        function testMissingSimDirIsSilent(tc)
            p = tc.Parameters;
            p.path.sim = fullfile(tc.TmpDir, 'does_not_exist');
            run_stats('add', 'acoustic_s', 1);
            run_stats('write', p, 1);   % must not error
            tc.verifyFalse(isfile(fullfile(p.path.sim, 'run_stats.jsonl')));
        end
    end

    methods
        function f = statsFile(tc)
            f = fullfile(tc.TmpDir, 'run_stats.jsonl');
        end
    end
end
