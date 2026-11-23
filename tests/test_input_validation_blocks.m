classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_input_validation_blocks < matlab.unittest.TestCase
% TEST_INPUT_VALIDATION_BLOCKS  Functions in hpc/, neuronav/, calibration/ and gui/ validate their inputs.
%
% The primary function of every .m file in these folders opens with an
% arguments block (directly after its help comment), so a wrong input fails
% at the call with a clear message. A file that cannot have one (a script, a
% function without inputs, varargin/inputParser parsing, inputs whose type
% varies too much to validate, or a large GUI app) goes in ALLOWLIST with a
% one-line reason. Add a block to a new function rather than an allowlist
% entry; remove the entry when a listed file gains a block.
%
%   results = runtests('tests/test_input_validation_blocks.m');

    properties (Constant)
        Folders = {'hpc', 'neuronav', 'calibration', 'gui'}

        % file (relative to functions/) -> reason it has no arguments block
        Allowlist = {
            'hpc/hpc_detect_system.m',                 'no inputs'
            'calibration/library_provenance.m',        'no inputs'
            'calibration/build_transducer_library.m',  'script (build_geometric_entry is a local function)'
            'calibration/calc_virtual_elem.m',         'script'
            'neuronav/neuronav_ingest_batch.m',        'varargin parsed with inputParser'
            'neuronav/neuronav_summary_to_mni.m',      'varargin parsed with inputParser'
            'neuronav/neuronav_convert_native_to_MNI.m', 'detects the optional side_order with exist(), which a default would defeat'
            'neuronav/neuronav_build_target_map.m',    'coil_map may be [] or a struct array'
            'neuronav/neuronav_check_target_map.m',    'target_map may be []'
            'neuronav/neuronav_label_map_from_scan.m', 'labels is a table; mode may be []'
            'neuronav/neuronav_map_label.m',           'label of any type; target_labels may be [], a struct or a cell'
            'neuronav/neuronav_match_entry_targets.m', 'all inputs may be [] or struct arrays'
            'gui/calibration_gui.m',                   'large GUI app'
            'gui/neuronav_gui.m',                      'large GUI app'
            }
    end

    methods (Test)

        function test_every_function_has_an_arguments_block(tc)
            root = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'functions');
            missing = {};
            for d = tc.Folders
                files = dir(fullfile(root, d{1}, '*.m'));
                for k = 1:numel(files)
                    rel = [d{1} '/' files(k).name];
                    if any(strcmp(rel, tc.Allowlist(:, 1)))
                        continue
                    end
                    if ~has_arguments_block(fullfile(files(k).folder, files(k).name))
                        missing{end + 1} = rel; %#ok<AGROW>
                    end
                end
            end
            tc.verifyEmpty(missing, sprintf(['Functions without an arguments block after the help ' ...
                'comment (add one, or an ALLOWLIST entry with the reason):\n  %s'], strjoin(missing, '\n  ')));
        end

        function test_allowlist_is_current(tc)
            % An entry for a file that is gone or now has a block is stale.
            root = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'functions');
            stale = {};
            for k = 1:size(tc.Allowlist, 1)
                f = fullfile(root, tc.Allowlist{k, 1});
                if ~isfile(f) || has_arguments_block(f)
                    stale{end + 1} = tc.Allowlist{k, 1}; %#ok<AGROW>
                end
            end
            tc.verifyEmpty(stale, sprintf('Stale ALLOWLIST entries (remove them):\n  %s', strjoin(stale, '\n  ')));
        end

        function test_detector_tells_blocks_apart(tc)
            % Negative control: the scan accepts a block after the help
            % comment (also after a continued function line) and rejects a
            % function without one, a block after other code, and a script.
            d = tempname;
            mkdir(d);
            cleanup = onCleanup(@() rmdir(d, 's'));
            cases = {
                ["function y = f(x)", "% F  Help", "%", "% More help", "", "arguments", "    x (1,1) double", "end", "", "y = x;", "end"], true
                ["function [a, b] = f(x, ...", "    z)", "% F  Help", "", "arguments", "    x", "    z = 1", "end", "a = x; b = z;", "end"], true
                ["function y = f(x)", "% F  Help", "", "y = x;", "end"], false
                ["function y = f(x)", "% F  Help", "", "y = x;", "arguments", "    x", "end", "end"], false
                ["% a script", "x = 1;", "arguments = 2;"], false
                };
            for k = 1:size(cases, 1)
                f = fullfile(d, sprintf('f%d.m', k));
                writelines(cases{k, 1}, f);
                tc.verifyEqual(has_arguments_block(f), cases{k, 2}, sprintf('case %d', k));
            end
        end

    end

end

function tf = has_arguments_block(file)
% True when the file's first function opens with an arguments block: after
% the (possibly continued) function line and the help comment, the first
% code line starts with 'arguments'. False for scripts.
    lines = strtrim(readlines(file));
    n = numel(lines);
    is_skip = @(s) s == "" || startsWith(s, "%");
    k = 1;
    while k <= n && is_skip(lines(k))
        k = k + 1;
    end
    tf = false;
    if k > n || isempty(regexp(lines(k), '^function\>', 'once'))
        return
    end
    while k <= n && ~isempty(regexp(lines(k), '\.\.\.', 'once'))   % continued signature
        k = k + 1;
    end
    k = k + 1;
    while k <= n && is_skip(lines(k))
        k = k + 1;
    end
    tf = k <= n && ~isempty(regexp(lines(k), '^arguments\>', 'once'));
end
