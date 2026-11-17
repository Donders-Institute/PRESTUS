classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_position_transducer_tus_entry < matlab.unittest.TestCase
% TEST_POSITION_TRANSDUCER_TUS_ENTRY  Target mask written by POSITION_TRANSDUCER_TUS_ENTRY.
%
% The spherical target mask is single data under the T1 header. Whatever the
% T1 datatype and intensity scaling, the mask is written as single with unit
% scaling. TUS_entry itself is not run: the stand-in Rscript cannot execute,
% so the call errors after the masks are written.
%
%   results = runtests('tests/test_position_transducer_tus_entry.m');

    methods (Test, TestTags = {'tus_entry'})

        function test_target_mask_ignores_t1_datatype(tc)
            % int16 and uint8 T1s (with scaling) failed in niftiwrite;
            % a single T1 is the negative control.
            for dtype = {'int16', 'uint8', 'single'}
                root = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
                [p, hdr] = make_case(root, dtype{1});
                try %#ok<TRYNC>
                    evalc('position_transducer_tus_entry(p, hdr)');
                end
                mask = fullfile(p.placement.tus_entry.work_dir, 'sub-001_target_target_mask.nii.gz');
                tc.assertTrue(isfile(mask), sprintf('no target mask for a %s T1', dtype{1}));
                info = niftiinfo(mask);
                tc.verifyEqual(info.Datatype, 'single', dtype{1});
                tc.verifyEqual(double(info.MultiplicativeScaling), 1, dtype{1});
                tc.verifyEqual(double(info.AdditiveOffset), 0, dtype{1});
                vol = niftiread(mask);
                tc.verifyEqual(vol(9, 9, 9), single(1), dtype{1});      % centre
                tc.verifyEqual(vol(1, 1, 1), single(0), dtype{1});      % outside r = 5 mm
                tc.verifyEqual(unique(vol(:))', single([0 1]), dtype{1});
            end
        end

    end
end

% =========================================================================

function [p, hdr] = make_case(root, dtype)
    n = 18;
    t1_path = fullfile(root, 'sub-001_T1w.nii');
    write_volume(t1_path, cast(100 * ones(n, n, n), dtype), 2);
    hdr = niftiinfo(t1_path);
    head_path = fullfile(root, 'head.nii');
    write_volume(head_path, ones(n, n, n, 'uint8'), 1);

    script = fullfile(root, 'TUS_entry.R');
    fclose(fopen(script, 'w'));
    env = fullfile(root, 'env');
    mkdir(fullfile(env, 'bin'));
    fclose(fopen(fullfile(env, 'bin', 'Rscript'), 'w'));   % not executable

    p = struct();
    p.subject_id = 1;
    p.path.anat = root;
    p.path.t1_pattern = 'sub-%03d_T1w.nii';
    p.path.seg = root;
    p.placement.tus_entry = struct('script_path', script, 'env_path', env, ...
        'scalp_mask', head_path, 'work_dir', fullfile(root, 'run'));
    p.transducer = struct('focus_pos', [9 9 9]);
end

function write_volume(path, vol, slope)
    % 1 mm isotropic NIfTI of class(vol) with the given scl_slope
    tmp = [tempname '.nii'];
    niftiwrite(vol, tmp);
    info = niftiinfo(tmp);
    delete(tmp);
    info.PixelDimensions = [1 1 1];
    info.MultiplicativeScaling = slope;
    info.AdditiveOffset = 0;
    niftiwrite(vol, path, info);
end
