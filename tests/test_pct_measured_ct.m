classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_pct_measured_ct < matlab.unittest.TestCase
% TEST_PCT_MEASURED_CT  A measured CT (path.ct_pattern) as the skull mapping's CT.
%
% PCT_CT_FILE resolves the pattern, PCT_TISSUES_MASK labels the skull as
% pct_create_pseudoCT.sh does, and PCT_LOAD_CT resamples the CT onto the
% segmentation grid through the world affines.
%
% Run with:   results = runtests('tests/test_pct_measured_ct.m');

    properties
        TmpDir
    end

    methods (TestMethodSetup)
        function setup(tc)
            tc.TmpDir = tempname;
            mkdir(tc.TmpDir);
            tc.addTeardown(@() rmdir(tc.TmpDir, 's'));
        end
    end

    methods (Test)
        function test_ct_file_empty_pattern(tc)
            p = base_params(tc.TmpDir);
            tc.verifyEqual(pct_ct_file(p), '');
        end

        function test_ct_file_missing_errors(tc)
            p = base_params(tc.TmpDir);
            p.path.ct_pattern = 'sub-%1$03d_CT.nii.gz';
            tc.verifyError(@() pct_ct_file(p), 'PRESTUS:pct:ctNotFound');
        end

        function test_tissues_mask(tc)
            seg = zeros(9, 9, 9, 'uint8');
            seg(3:7, 3:7, 3:7) = 5;          % scalp
            seg(5, 5, 5) = 7;                % compact bone inside the head
            seg(3, 5, 5) = 8;                % spongy bone at the head's edge
            m = pct_tissues_mask(seg);
            tc.verifyClass(m, 'uint8');
            tc.verifyEqual(m(5, 5, 5), uint8(4));
            tc.verifyEqual(m(4, 4, 4), uint8(4));     % one-voxel dilation of the inner bone
            tc.verifyEqual(m(1, 5, 5), uint8(4));     % skull-air interface (2,5,5), dilated once more
            tc.verifyEqual(m(7, 7, 7), uint8(5));     % scalp away from bone keeps its label
            tc.verifyEqual(m(9, 9, 9), uint8(0));
        end

        function test_load_ct_resamples_to_segmentation_grid(tc)
            p = base_params(tc.TmpDir);
            m2m = fullfile(p.path.seg, 'm2m_sub-009');
            mkdir(m2m);
            seg = zeros(10, 10, 10, 'uint8');
            seg(4:6, 4:6, 4:6) = 7;
            T_seg = eye(4); T_seg(4, 1:3) = [-5 -5 -5];          % 1 mm, row-vector convention
            write_nii(fullfile(m2m, 'final_tissues.nii.gz'), seg, T_seg);
            % CT at 0.5 mm covering the same world box; HU = 1000 inside the bone cube
            ct = -1000 * ones(20, 20, 20, 'single');
            ct(7:12, 7:12, 7:12) = 1000;   % world -2.25..0.25 mm: covers seg voxels 4:6 (world -2..0 mm)
            T_ct = diag([0.5 0.5 0.5 1]); T_ct(4, 1:3) = [-5.25 -5.25 -5.25];
            mkdir(p.path.anat);
            write_nii(fullfile(p.path.anat, 'sub-009_CT.nii.gz'), ct, T_ct);
            p.path.ct_pattern = 'sub-%1$03d_CT.nii.gz';

            [img, hdr, mask, mhdr, f] = pct_load_ct(p);
            tc.verifyEqual(size(img), [10 10 10]);
            tc.verifyEqual(hdr.Transform.T, mhdr.Transform.T, 'AbsTol', 1e-6);
            tc.verifySubstring(f, 'sub-009_CT.nii.gz');
            tc.verifyEqual(double(img(5, 5, 5)), 1000, 'AbsTol', 1);
            tc.verifyEqual(double(img(1, 1, 1)), -1000, 'AbsTol', 1);
            tc.verifyEqual(mask(5, 5, 5), uint8(4));
        end
    end
end

function p = base_params(root)
    p = struct('subject_id', 9, 'path', struct('seg', fullfile(root, 'simnibs'), ...
        'anat', fullfile(root, 'anat'), 'pct', ''), 'pct', struct('enabled', 1));
end

function write_nii(file, img, T)
    [folder, name] = fileparts(file);
    [~, name] = fileparts(name);                    % strip .nii of .nii.gz
    niftiwrite(img, fullfile(folder, name), 'Compressed', true);
    info = niftiinfo(fullfile(folder, [name '.nii.gz']));
    info.Transform = affine3d(T);
    info.TransformName = 'Sform';
    niftiwrite(img, fullfile(folder, name), info, 'Compressed', true);
end
