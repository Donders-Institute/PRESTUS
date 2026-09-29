classdef test_parse_memory_tags < matlab.unittest.TestCase
% TEST_PARSE_MEMORY_TAGS  Unit tests for functions/helper/parse_memory_tags
%
%   Run with:   results = runtests('tests/test_parse_memory_tags.m');

    properties (Constant)
        MatlabLine = ['[memory] after BabelBrain exit: MATLAB rss=2.9 GB (peak 3.1 GB) | ' ...
                      'cgroup current=3.0 GB, peak=9.4 GB, limit=64.0 GB'];
        BbLine     = ['[memory_bb] t=40s: python rss=6.0 GB (peak 7.2 GB) | ' ...
                      'cgroup current=9.0 GB, peak=9.8 GB'];
    end

    methods (Test, TestTags = {'parse_memory_tags'})

        function test_empty_text_returns_nan(tc)
            mem = parse_memory_tags('');
            tc.verifyTrue(all(isnan([mem.bb_peak_gb, mem.matlab_peak_gb, ...
                                     mem.cgroup_peak_gb, mem.cgroup_limit_gb])));
            tc.verifyEqual([mem.n_bb_lines, mem.n_matlab_lines], [0 0]);
        end

        function test_tags_are_kept_apart(tc)
            mem = parse_memory_tags(strjoin({tc.MatlabLine, tc.BbLine}, newline));
            tc.verifyEqual(mem.bb_peak_gb, 7.2, 'AbsTol', 1e-9);
            tc.verifyEqual(mem.matlab_peak_gb, 3.1, 'AbsTol', 1e-9);
            tc.verifyEqual(mem.cgroup_limit_gb, 64, 'AbsTol', 1e-9);
            tc.verifyEqual([mem.n_bb_lines, mem.n_matlab_lines], [1 1]);
        end

        function test_cgroup_peak_is_max_over_both_tags(tc)
            mem = parse_memory_tags(strjoin({tc.MatlabLine, tc.BbLine}, newline));
            tc.verifyEqual(mem.cgroup_peak_gb, 9.8, 'AbsTol', 1e-9);
        end

        function test_duplicate_sources_are_idempotent(tc)
            % Pipeline log + SLURM log carry the same lines
            once  = parse_memory_tags(strjoin({tc.MatlabLine, tc.BbLine}, newline));
            twice = parse_memory_tags(strjoin({tc.MatlabLine, tc.BbLine, ...
                                               tc.MatlabLine, tc.BbLine}, newline));
            tc.verifyEqual(twice.bb_peak_gb,     once.bb_peak_gb);
            tc.verifyEqual(twice.matlab_peak_gb, once.matlab_peak_gb);
            tc.verifyEqual(twice.cgroup_peak_gb, once.cgroup_peak_gb);
        end

        function test_nan_values_are_skipped(tc)
            % MATLAB prints NaN, Python prints nan when a figure is unavailable
            txt = strjoin({ ...
                ['[memory] x: MATLAB rss=NaN GB (peak NaN GB) | ' ...
                 'cgroup current=NaN GB, peak=NaN GB, limit=NaN GB'], ...
                ['[memory_bb] exit: python rss=1.0 GB (peak 2.0 GB) | ' ...
                 'cgroup current=nan GB, peak=nan GB']}, newline);
            mem = parse_memory_tags(txt);
            tc.verifyEqual(mem.bb_peak_gb, 2.0, 'AbsTol', 1e-9);
            tc.verifyTrue(isnan(mem.matlab_peak_gb));
            tc.verifyTrue(isnan(mem.cgroup_peak_gb));
            tc.verifyEqual(mem.n_matlab_lines, 1);
        end

        function test_summary_and_unrelated_lines_ignored(tc)
            txt = strjoin({ ...
                ['[memory_summary] BabelBrain: BabelBrain peak=99.0 GB | MATLAB peak=99.0 GB | ' ...
                 'cgroup peak=99.0 GB, limit=64.0 GB'], ...
                '[babelbrain_wrapper] BabelBrain finished in 12.0 s.', ...
                tc.BbLine}, newline);
            mem = parse_memory_tags(txt);
            tc.verifyEqual(mem.bb_peak_gb, 7.2, 'AbsTol', 1e-9);
            tc.verifyEqual(mem.cgroup_peak_gb, 9.8, 'AbsTol', 1e-9);
            tc.verifyEqual(mem.n_matlab_lines, 0);
        end

        function test_string_input(tc)
            mem = parse_memory_tags(string(tc.BbLine));
            tc.verifyEqual(mem.bb_peak_gb, 7.2, 'AbsTol', 1e-9);
        end

    end

end
