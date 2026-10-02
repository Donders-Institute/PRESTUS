classdef test_demo_golden < matlab.unittest.TestCase
% TEST_DEMO_GOLDEN  Compare demo-subject key metrics with a stored reference.
%
% The demo tests (test_integration_demo, tags demo_acoustic / demo_thermal)
% check that outputs are written; this test checks that their numbers have not
% drifted. It reads the result table of the manual-placement demo run
% (sim_outputs/sub-009/.../sub-009_layered*.csv) and compares key metrics with
% golden/sub-009_manual_metrics.json in PRESTUS_DEMO_DATA:
%
%   Isppa_brain, Psptp_brain, MI_brain        relative tolerance (RelTol)
%   Ix/Iy/Iz_brain_vox (focus position)       within one voxel
%   maxT_brain, CEM43iso_brain (if thermal)   absolute / relative tolerance
%
% Without a reference, the current metrics are written to
% golden/sub-009_manual_metrics.candidate.json and the test is marked
% incomplete: review the candidate (e.g. against a trusted release), then
% rename it to sub-009_manual_metrics.json. A reference is never accepted
% automatically. Update it only for intended changes, and say so in the
% changelog.
%
%   run_all_tests('demo_acoustic'); run_all_tests('demo_thermal');
%   run_all_tests('demo_golden')

    properties (Constant)
        SUBJECT_ID = 9
        RelTol     = struct('Isppa_brain', 0.05, 'Psptp_brain', 0.03, 'MI_brain', 0.03, ...
                            'CEM43iso_brain', 0.20)
        AbsTol     = struct('Ix_brain_vox', 1, 'Iy_brain_vox', 1, 'Iz_brain_vox', 1, ...
                            'maxT_brain', 0.2)
    end

    properties
        DataRoot
    end

    methods (TestClassSetup)
        function setup(tc)
            here = fileparts(mfilename('fullpath'));
            addpath(genpath(fullfile(fileparts(here), 'functions')));
            tc.DataRoot = getenv('PRESTUS_DEMO_DATA');
            tc.assumeNotEmpty(tc.DataRoot, 'PRESTUS_DEMO_DATA not set — skipping golden metrics.');
        end
    end

    methods (Test, TestTags = {'demo_golden'})

        function test_manual_metrics_match_reference(tc)
            current = tc.current_metrics();
            golden_dir  = fullfile(tc.DataRoot, 'golden');
            golden_file = fullfile(golden_dir, sprintf('sub-%03d_manual_metrics.json', tc.SUBJECT_ID));
            if ~isfile(golden_file)
                if ~isfolder(golden_dir); mkdir(golden_dir); end
                cand = strrep(golden_file, '.json', '.candidate.json');
                write_json(cand, current);
                tc.assumeFail(sprintf(['No reference metrics yet. Wrote %s — review it, then ' ...
                    'rename it to %s.'], cand, golden_file));
            end
            golden = jsondecode(fileread(golden_file));

            keys = fieldnames(golden);
            for i = 1:numel(keys)
                k = keys{i};
                if ~isfield(current, k) || ~isfinite(current.(k))
                    tc.verifyFail(sprintf('Metric %s is missing from the current results.', k));
                    continue
                end
                msg = sprintf('%s: current %.6g, reference %.6g', k, current.(k), golden.(k));
                if isfield(tc.RelTol, k)
                    tc.verifyEqual(current.(k), golden.(k), 'RelTol', tc.RelTol.(k), msg);
                elseif isfield(tc.AbsTol, k)
                    tc.verifyEqual(current.(k), golden.(k), 'AbsTol', tc.AbsTol.(k), msg);
                else
                    tc.verifyEqual(current.(k), golden.(k), 'RelTol', 1e-6, msg);
                end
            end
        end

    end

    methods (Access = private)

        function m = current_metrics(tc)
            sub_dir = fullfile(tc.DataRoot, 'sim_outputs', sprintf('sub-%03d', tc.SUBJECT_ID));
            files = dir(fullfile(sub_dir, '**', sprintf('sub-%03d_layered*.csv', tc.SUBJECT_ID)));
            files = files(~contains({files.name}, {'targeting', 'timeseries', 'water', 'baseline'}));
            tc.assumeNotEmpty(files, sprintf(['No demo result table in %s; run ' ...
                'run_all_tests(''demo_acoustic'') (and ''demo_thermal'') first.'], sub_dir));
            [~, newest] = max([files.datenum]);
            T = readtable(fullfile(files(newest).folder, files(newest).name));

            m = struct();
            names = [fieldnames(tc.RelTol); fieldnames(tc.AbsTol)];
            for i = 1:numel(names)
                k = names{i};
                if ismember(k, T.Properties.VariableNames)
                    v = T.(k);
                    if iscell(v); v = str2double(v); end
                    if ~isempty(v) && isfinite(v(1))
                        m.(k) = double(v(1));
                    end
                end
            end
            tc.assumeTrue(isfield(m, 'Isppa_brain'), ...
                'The demo result table has no Isppa_brain column (acoustic analysis not run).');
        end

    end

end

function write_json(path, s)
    fid = fopen(path, 'w');
    fprintf(fid, '%s\n', jsonencode(s, 'PrettyPrint', true));
    fclose(fid);
end
