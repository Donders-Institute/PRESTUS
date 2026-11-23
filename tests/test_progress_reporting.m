classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_progress_reporting < matlab.unittest.TestCase
% TEST_PROGRESS_REPORTING  Unit tests for PRESTUS_PROGRESS and LIVE_VIEW (read by the Nexus)
%
%   Run with:   results = runtests('tests/test_progress_reporting.m');

    methods (Test)

        function test_progress_line_format(tc)
            out = evalc('prestus_progress(''acoustic'', 0, 1, ''solver cpp_gpu'')');
            tok = regexp(strtrim(out), '^\[progress\] acoustic 0/1 @([\d.]+) solver cpp_gpu$', 'tokens', 'once');
            tc.assertNotEmpty(tok, out);
            tc.verifyLessThan(abs(str2double(tok{1}) - posixtime(datetime('now', 'TimeZone', 'UTC'))), 60);
        end

        function test_progress_rate_limited_but_start_and_end_printed(tc)
            stage = sprintf('t%d', randi(1e6));
            out = evalc(sprintf(['prestus_progress(''%1$s'', 0, 1000); ' ...
                'for k = 1:999, prestus_progress(''%1$s'', k, 1000); end; ' ...
                'prestus_progress(''%1$s'', 1000, 1000);'], stage));
            lines = splitlines(strtrim(out));
            tc.verifyLessThan(numel(lines), 10, 'a fast loop prints only a few lines');
            tc.verifyTrue(startsWith(lines{1}, ['[progress] ' stage ' 0/1000']));
            tc.verifyTrue(startsWith(lines{end}, ['[progress] ' stage ' 1000/1000']));
        end

        function test_progress_detail_handle_only_evaluated_when_printed(tc)
            stage = sprintf('h%d', randi(1e6));
            evalc(sprintf('prestus_progress(''%s'', 1, 100)', stage));
            calls = containers.Map({'n'}, {0});          % handle: counts the evaluations
            detail = @() bump(calls);
            evalc('prestus_progress(stage, 1, 100, detail)');   % same k: rate-limited, skipped
            tc.verifyEqual(calls('n'), 0);
            out = evalc('prestus_progress(stage, 100, 100, detail)');
            tc.verifyEqual(calls('n'), 1);
            tc.verifySubstring(out, ' x');
        end

        function test_live_view_off_by_default_and_frame_written_when_on(tc)
            p = struct('simulation', struct('live_view', 0), 'subject_id', 3, ...
                'io', struct('dir_output', tempname, 'output_affix', '_lv'));
            tc.addTeardown(@() rmdir_if(p.io.dir_output));
            tc.verifyFalse(live_view('enabled', p));
            fig = figure('Visible', 'off');
            tc.addTeardown(@() close(fig));
            plot(1:3);
            evalc('live_view(''frame'', p, fig, ''thermal'', true)');
            file = fullfile(p.io.dir_output, 'live', 'thermal_sub-003_lv.png');
            tc.verifyFalse(isfile(file), 'live_view 0 writes nothing');
            p.simulation.live_view = 1;
            out = evalc('live_view(''frame'', p, fig, ''thermal'', true)');
            tc.verifyTrue(isfile(file));
            tc.verifySubstring(out, ['[live] thermal ' file]);
            out = evalc('live_view(''frame'', p, fig, ''thermal'', true)');
            tc.verifyEmpty(strtrim(out), 'the path is announced once');
            tc.verifyEmpty(dir(fullfile(p.io.dir_output, 'live', '*.part.png')));
        end

        function test_live_view_frame_of_a_figure_without_number(tc)
            % k-Wave's plot is captured as the newest figure; one without an integer
            % handle (empty Number) must not stop the capture
            p = struct('simulation', struct('live_view', 1), 'subject_id', 5, ...
                'io', struct('dir_output', tempname, 'output_affix', ''));
            tc.addTeardown(@() rmdir_if(p.io.dir_output));
            fig = figure('Visible', 'off', 'IntegerHandle', 'off');
            tc.addTeardown(@() close(fig));
            plot(1:3);
            t = live_view('watch', p, 'acoustic');
            out = evalc('live_view(''stop'', p, t)');
            png = fullfile(p.io.dir_output, 'live', 'acoustic_sub-005.png');
            tc.verifyTrue(isfile(png), out);
        end

        function test_live_view_data_frame(tc)
            p = struct('simulation', struct('live_view', 1), 'subject_id', 4, ...
                'io', struct('dir_output', tempname, 'output_affix', '_d'));
            tc.addTeardown(@() rmdir_if(p.io.dir_output));
            plane = single(reshape(1:12, 3, 4)) + 36;
            meta = struct('quantity', 'temperature', 'units', '°C', 'label', 'Temperature', ...
                'spacing_mm', [0.5 1], 'markers', struct('trans', [1 0]));
            evaluated = containers.Map({'n'}, {0});
            out = evalc('live_view(''data'', p, ''thermal'', @() bump_plane(evaluated, plane), meta, true)');
            js = fullfile(p.io.dir_output, 'live', 'thermal_sub-004_d.json');
            tc.verifySubstring(out, ['[live] thermal ' js]);
            tc.verifyEqual(live_view('file', p, 'thermal'), js, 'file: the data frame path');
            m = jsondecode(fileread(js));
            tc.verifyEqual(m.shape(:)', [3 4]);
            tc.verifyEqual(m.order, 'F');
            tc.verifyEqual(m.units, '°C');
            fid = fopen(fullfile(p.io.dir_output, 'live', m.data), 'r', 'ieee-le');
            back = reshape(fread(fid, Inf, 'single=>single'), m.shape(:)');
            fclose(fid);
            tc.verifyEqual(back, plane);
            evalc('live_view(''data'', p, ''thermal'', @() bump_plane(evaluated, plane), meta)');
            tc.verifyEqual(evaluated('n'), 1, 'within 5 s: not due, the plane is not computed');
        end

    end
end

function plane = bump_plane(calls, plane)
    calls('n') = calls('n') + 1; %#ok<NASGU> (handle object)
end

function s = bump(calls)
    calls('n') = calls('n') + 1; %#ok<NASGU> (handle object)
    s = 'x';
end

function rmdir_if(d)
    if isfolder(d); rmdir(d, 's'); end
end
