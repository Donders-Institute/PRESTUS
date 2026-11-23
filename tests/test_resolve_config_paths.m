classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_resolve_config_paths < matlab.unittest.TestCase
% TEST_RESOLVE_CONFIG_PATHS  Relative file and folder settings become absolute.
%
% PATH_LOG_SETUP takes relative settings from the PRESTUS root through
% RESOLVE_CONFIG_PATHS instead of changing MATLAB's working folder. Checked:
% relative text, string, and cell values are resolved; absolute, empty, ~,
% and $ values, path templates, and unlisted keys are left as they are;
% struct arrays (several transducers) are resolved per element; PATH_LOG_SETUP
% no longer changes the working folder.
%
%   results = runtests('tests/test_resolve_config_paths.m');

    properties (Constant)
        Base = fullfile(tempdir, 'prestus_root')
    end

    methods (Test)

        function test_relative_text_values(tc)
            p.path = struct('sim', 'out/sims', 'seg', "seg", 'transducer_library', {{'lib_a', "lib_b"}});
            [r, lines] = resolve_config_paths(p, tc.Base);
            tc.verifyEqual(r.path.sim, fullfile(tc.Base, 'out/sims'));
            tc.verifyEqual(r.path.seg, string(fullfile(tc.Base, 'seg')));
            tc.verifyEqual(r.path.transducer_library, {fullfile(tc.Base, 'lib_a'), string(fullfile(tc.Base, 'lib_b'))});
            tc.verifyNumElements(lines, 3);
        end

        function test_left_unchanged(tc)
            % Negative control: absolute, empty, home, variable, template, and
            % unlisted values pass through, and nothing is reported.
            p.path = struct('sim', '/data/sims', 'anat', '', 'seg', '~/seg', 'pct', '$SCRATCH/pct', ...
                't1_pattern', 'sub-%1$03d_T1w.nii*', 'placement', 'C:\study\placement');
            p.jwave = struct('python_env', 'prestus-jwave');
            p.transducer = struct('library_path', []);
            [r, lines] = resolve_config_paths(p, tc.Base);
            tc.verifyEqual(r, p);
            tc.verifyEmpty(lines);
        end

        function test_struct_array(tc)
            ef1 = struct('file_path', 'config/elems_a.csv');
            ef2 = struct('file_path', '/abs/elems_b.csv');
            mk = @(ef) struct('matrix', struct('matrix_shape', struct('extract_from_file', ef)));
            p.transducer = [mk(ef1), mk(ef2)];
            r = resolve_config_paths(p, tc.Base);
            tc.verifyEqual(r.transducer(1).matrix.matrix_shape.extract_from_file.file_path, ...
                fullfile(tc.Base, 'config/elems_a.csv'));
            tc.verifyEqual(r.transducer(2).matrix.matrix_shape.extract_from_file.file_path, '/abs/elems_b.csv');
        end

        function test_path_log_setup_keeps_working_folder(tc)
            % path_log_setup runs only with k-Wave and a full config, so its
            % source is checked: no change of MATLAB's working folder.
            src = fileread(which('path_log_setup'));
            code = regexprep(src, '%[^\n]*', '');     % drop comments
            tc.verifyEmpty(regexp(code, '\<cd\s*\(', 'once'));
            tc.verifyNotEmpty(regexp(code, 'resolve_config_paths\(', 'once'));
        end

    end

end
