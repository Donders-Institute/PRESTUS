classdef prestus_test_resources < matlab.unittest.plugins.TestRunnerPlugin
% PRESTUS_TEST_RESOURCES  Test-runner plugin that records the resources of
% every test and every test class: wall time, MATLAB CPU time, and memory.
%
%   Rows has one row per test (kind 'test') and per class (kind 'class'; its
%   time includes class setup and teardown, so it exceeds the sum of its
%   tests). run_all_tests adds the outcome and writes the table to CSV.
%
%   Memory is the resident set size of the MATLAB process, read from
%   /proc/self/status (Linux only; NaN elsewhere):
%     rss_start_mb   resident memory when the test starts
%     peak_mb        highest resident memory during the test
%     peak_delta_mb  peak_mb - rss_start_mb: what the test itself adds
%   The peak is reset before each test by writing 5 to /proc/self/clear_refs.
%   Where the kernel refuses that reset, peak_mb is the process peak so far
%   and peak_reset is false in the table.
%
%   cpu_s is process CPU time (all MATLAB threads); cpu_s / wall_s is the
%   number of cores the test kept busy.
%
% See also: RUN_ALL_TESTS

    properties (SetAccess = private)
        Rows = struct('kind', {}, 'name', {}, 'wall_s', {}, 'cpu_s', {}, ...
                      'rss_start_mb', {}, 'peak_mb', {}, 'peak_delta_mb', {}, 'peak_reset', {})
    end

    properties (Access = private)
        ClassPeak = NaN
    end

    methods (Access = protected)

        function runTestClass(plugin, pluginData)
            reset_ok = reset_peak();
            [rss0, ~] = read_rss();
            plugin.ClassPeak = rss0;
            t0 = tic; c0 = cputime;
            runTestClass@matlab.unittest.plugins.TestRunnerPlugin(plugin, pluginData);
            [~, hwm] = read_rss();
            plugin.ClassPeak = max(plugin.ClassPeak, hwm);
            plugin.add_row('class', pluginData.Name, toc(t0), cputime - c0, rss0, plugin.ClassPeak, reset_ok);
        end

        function runTest(plugin, pluginData)
            % The class peak so far, before this test's reset clears it
            [~, hwm] = read_rss();
            plugin.ClassPeak = max(plugin.ClassPeak, hwm);
            reset_ok = reset_peak();
            [rss0, ~] = read_rss();
            t0 = tic; c0 = cputime;
            runTest@matlab.unittest.plugins.TestRunnerPlugin(plugin, pluginData);
            wall = toc(t0); cpu = cputime - c0;
            [~, peak] = read_rss();
            plugin.ClassPeak = max(plugin.ClassPeak, peak);
            plugin.add_row('test', pluginData.Name, wall, cpu, rss0, peak, reset_ok);
        end

    end

    methods (Access = private)
        function add_row(plugin, kind, name, wall, cpu, rss0, peak, reset_ok)
            plugin.Rows(end+1) = struct('kind', kind, 'name', name, 'wall_s', wall, ...
                'cpu_s', cpu, 'rss_start_mb', rss0, 'peak_mb', peak, ...
                'peak_delta_mb', peak - rss0, 'peak_reset', reset_ok);
        end
    end
end

function ok = reset_peak()
% Linux >= 4.0: writing 5 to clear_refs resets VmHWM to the current RSS
ok = false;
fid = fopen('/proc/self/clear_refs', 'w');
if fid < 0; return; end
ok = fprintf(fid, '5') == 1;
ok = fclose(fid) == 0 && ok;
end

function [rss_mb, hwm_mb] = read_rss()
% Current (VmRSS) and peak (VmHWM) resident memory of this process, in MB
rss_mb = NaN; hwm_mb = NaN;
txt = fileread_or_empty('/proc/self/status');
if isempty(txt); return; end
rss_mb = field_mb(txt, 'VmRSS');
hwm_mb = field_mb(txt, 'VmHWM');
end

function v = field_mb(txt, key)
tok = regexp(txt, [key ':\s*(\d+)\s*kB'], 'tokens', 'once');
if isempty(tok); v = NaN; else; v = str2double(tok{1}) / 1024; end
end

function txt = fileread_or_empty(f)
try
    txt = fileread(f);
catch
    txt = '';
end
end
