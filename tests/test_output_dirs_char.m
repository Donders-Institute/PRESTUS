classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_output_dirs_char < matlab.unittest.TestCase
% TEST_OUTPUT_DIRS_CHAR  Output folders are char vectors when path.sim is a string.
%
%   YAML values load as MATLAB strings. Writers build file names as
%   [io.dir_xxx 'name.nii'], which on a string gives a 1x2 string array that
%   niftiwrite rejects. GET_OUTPUT_DIR and PATH_LOG_SETUP return char folders
%   for string and char path.sim alike.
%
%   results = runtests('tests/test_output_dirs_char.m');

    methods (Test)

        function get_output_dir_returns_char(tc)
            for sim = {'/data/sim', "/data/sim"}
                p = struct('path', struct('sim', sim{1}), 'subject_id', 7);
                d = get_output_dir(p);
                tc.verifyClass(d, 'char');
                tc.verifyEqual(d, fullfile('/data/sim', 'sub-007'));
            end
        end

        function path_log_setup_dirs_are_char(tc)
            root = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            here = pwd;
            tc.addTeardown(@() cd(here));
            tc.addTeardown(@() diary('off'));
            prestus_root = fileparts(fileparts(mfilename('fullpath')));
            p = load_parameters();
            tc.assertClass(p.path.sim, 'string', 'YAML paths no longer load as strings; test premise changed');
            p.subject_id = 1;
            p.simulation.code_type = 'matlab_cpu';      % gpu code types query the GPU here
            p.path.sim = string(fullfile(root, 'sim'));
            p.path.seg = string(fullfile(root, 'seg'));
            evalc('p = path_log_setup(p, prestus_root);');
            diary('off');
            f = fieldnames(p.io);
            dirs = f(startsWith(f, 'dir_'));
            tc.assertNotEmpty(dirs);
            for i = 1:numel(dirs)
                v = p.io.(dirs{i});
                if isempty(v); continue; end
                tc.verifyClass(v, 'char', sprintf('io.%s', dirs{i}));
            end
            % the failing pattern from nifti_acoustic: one file name, not a string array
            name = [fullfile(p.io.dir_nii_T1w, 'x') '.nii'];
            tc.verifyClass(name, 'char');
        end

    end

end
