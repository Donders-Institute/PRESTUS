classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_segmentation_run < matlab.unittest.TestCase
% TEST_SEGMENTATION_RUN  charm call built by SEGMENTATION_RUN (platform 'matlab', stub charm).
%
%   A stub charm records its arguments. charm gets --forcerun only when the
%   finished segmentation (m2m_sub-NNN/final_tissues.nii.gz) is missing,
%   with or without a trailing separator on path.seg.
%
%   results = runtests('tests/test_segmentation_run.m');

    methods (Test, TestTags = {'segmentation'})

        function test_forcerun_only_without_final_tissues(tc)
            tc.assumeFalse(ispc, 'stub charm is a shell script');
            root = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            bin = fullfile(root, 'bin');
            mkdir(bin);
            charm = fullfile(bin, 'charm');
            fid = fopen(charm, 'w');
            fprintf(fid, '#!/bin/sh\nprintf ''%%s\\n'' "$@" > "%s"\n', fullfile(bin, 'args.txt'));
            fclose(fid);
            fileattrib(charm, '+x');

            seg = fullfile(root, 'seg');
            m2m = fullfile(seg, 'm2m_sub-001');
            mkdir(m2m);
            t1 = fullfile(root, 'T1.nii.gz');
            fclose(fopen(t1, 'w'));
            p = struct('path', struct('seg', seg), 'platform', 'matlab', ...
                'startup', struct('simnibs_bin_path', bin), 'io', struct());

            final = fullfile(m2m, 'final_tissues.nii.gz');
            fclose(fopen(final, 'w'));
            for seg_arg = {seg, [seg filesep]}
                q = p;
                q.path.seg = seg_arg{1};
                args = run_stub(bin, seg_arg{1}, t1, q);
                tc.verifyFalse(ismember('--forcerun', args), ...
                    sprintf('--forcerun with final_tissues present (path.seg = %s)', seg_arg{1}));
                tc.verifyTrue(ismember(t1, args));
            end

            % Negative control: segmentation incomplete
            delete(final);
            args = run_stub(bin, seg, t1, p);
            tc.verifyTrue(ismember('--forcerun', args));
        end

    end

end

function args = run_stub(bin, seg, t1, p)
    out = fullfile(bin, 'args.txt');
    if isfile(out); delete(out); end
    w = warning('off', 'all');
    restore = onCleanup(@() warning(w));
    evalc('segmentation_run(seg, 1, t1, '''', p)');
    args = strtrim(splitlines(string(fileread(out))));
    args = cellstr(args(args ~= ""));
end
