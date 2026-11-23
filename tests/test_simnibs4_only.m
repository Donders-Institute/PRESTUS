classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_simnibs4_only < matlab.unittest.TestCase
% PRESTUS uses charm segmentations (SimNIBS 4 and later), which hold only the
% nonlinear MNI warp: the coordinate wrappers refuse the affine transforms of
% SimNIBS 3 (headreco), and no code asks for them or reads headreco's files.

    properties (TestParameter)
        wrapper = {'mni2subject_coords_LDfix', 'subject2mni_coords_LDfix'}
        affine  = {'12dof', '6dof'}
    end

    methods (Test)
        function test_wrapper_refuses_affine_transforms(tc, wrapper, affine)
            f = str2func(wrapper);
            tc.verifyError(@() f([0 0 0], tempdir, struct(), affine), 'prestus:transform:simnibs4Only');
        end

        function test_wrapper_accepts_nonl(tc, wrapper)
            % negative control: 'nonl' passes the check (and fails later, on the
            % missing SimNIBS install, with another error)
            f = str2func(wrapper);
            id = '';
            try
                f([0 0 0], fullfile(tempdir, 'no_m2m_here'), struct(), 'nonl');
            catch ME
                id = ME.identifier;
            end
            tc.verifyNotEqual(id, 'prestus:transform:simnibs4Only');
        end

        function test_no_headreco_transforms_in_code(tc)
            % Every MATLAB and Nexus source file: no affine transform type and no
            % file that only headreco writes
            root = fileparts(fileparts(mfilename('fullpath')));
            files = [dir(fullfile(root, 'functions', '**', '*.m')); ...
                     dir(fullfile(root, 'nexus', 'prestus_nexus', '**', '*.py'))];
            pattern = '[''"](12dof|6dof)[''"]|MNI2conform_12DOF|T1fs_nu_12DOF|T1fs_conform';
            hits = {};
            for k = 1:numel(files)
                txt = fileread(fullfile(files(k).folder, files(k).name));
                if ~isempty(regexp(txt, pattern, 'once'))
                    hits{end+1} = fullfile(files(k).folder, files(k).name); %#ok<AGROW>
                end
            end
            tc.verifyGreaterThan(numel(files), 100, 'source scan found too few files');
            tc.verifyEmpty(hits, sprintf('SimNIBS 3 (headreco) transforms or files in: %s', strjoin(hits, ', ')));
        end

        function test_scan_pattern_finds_a_headreco_call(tc)
            % negative control for the scan: the pattern matches the old call
            pattern = '[''"](12dof|6dof)[''"]|MNI2conform_12DOF|T1fs_nu_12DOF|T1fs_conform';
            tc.verifyNotEmpty(regexp('mni2subject_coords_LDfix(p, m2m, parameters, ''12dof'')', pattern, 'once'));
            tc.verifyNotEmpty(regexp('readmatrix(fullfile(m, ''toMNI'', ''MNI2conform_12DOF.txt''))', pattern, 'once'));
        end
    end
end
