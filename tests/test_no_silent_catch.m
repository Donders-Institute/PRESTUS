classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_no_silent_catch < matlab.unittest.TestCase
% TEST_NO_SILENT_CATCH  No catch block in functions/ discards an error without a reason.
%
% An empty catch hides failures that change results (a cache judged valid, a
% file not read). Prefer an explicit check (ISFILE, ISFIELD, WHO('-file')) so
% no try is needed; when catching is well motivated, either handle the error
% (warn with its message, rethrow, fall back visibly) or state in a comment
% inside the catch why discarding it is safe, e.g.
%
%   catch
%       % Telemetry must never interrupt the pipeline.
%   end
%
%   results = runtests('tests/test_no_silent_catch.m');

    methods (Test)

        function test_empty_catch_blocks_state_a_reason(tc)
            root = fileparts(fileparts(mfilename('fullpath')));
            files = [dir(fullfile(root, 'functions', '**', '*.m')); dir(fullfile(root, 'prestus_pipeline.m'))];
            bad = {};
            for k = 1:numel(files)
                f = fullfile(files(k).folder, files(k).name);
                bad = [bad, silent_catches(f, root)]; %#ok<AGROW>
            end
            tc.verifyEmpty(bad, sprintf(['Empty catch blocks without a comment stating why ' ...
                'the error is safe to discard:\n  %s'], strjoin(bad, '\n  ')));
        end

        function test_detector_finds_a_silent_catch(tc)
            % Negative control: the scan must flag an empty catch and accept
            % one with a reason.
            f = [tempname '.m'];
            cleanup = onCleanup(@() delete(f));
            writelines(["function x = f()", "x = 1;", "try", "    x = 2;", "catch", "end", ...
                        "try", "    x = 3;", "catch", "    % Safe: optional.", "end", "end"], f);
            hits = silent_catches(f, fileparts(f));
            tc.verifyNumElements(hits, 1);
            tc.verifySubstring(hits{1}, ':5');
        end

    end

end

function hits = silent_catches(file, root)
% 'relative/path.m:line' for each catch followed (blank lines aside) by its end
    lines = strtrim(readlines(file));
    hits = {};
    for i = 1:numel(lines)
        % a comment on the catch line itself also states a reason
        if isempty(regexp(lines(i), '^catch(\s+\w+)?\s*$', 'once'))
            continue
        end
        j = i + 1;
        while j <= numel(lines) && lines(j) == ""
            j = j + 1;
        end
        % \> is MATLAB's end-of-word anchor (\b would be a backspace character)
        if j <= numel(lines) && ~isempty(regexp(lines(j), '^end\>', 'once'))
            rel = erase(file, [root filesep]);
            hits{end+1} = sprintf('%s:%d', rel, i); %#ok<AGROW>
        end
    end
end
