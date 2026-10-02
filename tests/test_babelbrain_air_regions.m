classdef test_babelbrain_air_regions < matlab.unittest.TestCase
% TEST_BABELBRAIN_AIR_REGIONS  Air mask written for BabelBrain's pseudoCT path (issue #102).
%
% Runs _write_air_regions from run_babelbrain_case.py on a synthetic CT:
% a soft-tissue head with an internal air pocket, surrounded by air, with
% a wall that splits the exterior air into two pieces (as cropping can).
% Needs Python 3 with numpy, scipy, nibabel, pyyaml (BABELBRAIN_PYTHON or
% python3 on the PATH); skipped otherwise. No BabelBrain or GPU required.
%
%   results = runtests('tests/test_babelbrain_air_regions.m');

    properties
        PyBin
        Tmp
        Runner
    end

    methods (TestClassSetup)
        function setup(tc)
            repo_root = fileparts(fileparts(mfilename('fullpath')));
            tc.Runner = fullfile(repo_root, 'functions', 'babelbrain');
            tc.PyBin = getenv('BABELBRAIN_PYTHON');
            if isempty(tc.PyBin); tc.PyBin = 'python3'; end
            [st, ~] = system(sprintf('"%s" -c "import numpy, scipy, nibabel, yaml"', tc.PyBin));
            tc.assumeEqual(st, 0, 'Python with numpy/scipy/nibabel/pyyaml not available.');
        end
    end

    methods (TestMethodSetup)
        function make_tmp(tc)
            tc.Tmp = tempname; mkdir(fullfile(tc.Tmp, 'sub'));
            tc.addTeardown(@() rmdir(tc.Tmp, 's'));
        end
    end

    methods (Test)
        function test_internal_air_kept_exterior_removed(tc)
            out = tc.run_py(false);
            tc.verifyEqual(out.outside_pocket, 0, 'exterior air must not enter the mask');
            tc.verifyGreaterThan(out.in_pocket, 0.6 * out.pocket_vox, 'internal pocket should be kept');
            tc.verifyTrue(out.file_exists);
        end

        function test_keep_boundary_retains_split_exterior(tc)
            % BabelBrain's own recipe (largest component only) keeps the second
            % exterior piece -- the reason PRESTUS drops boundary components
            out = tc.run_py(true);
            tc.verifyGreaterThan(out.outside_pocket, 1000);
        end
    end

    methods (Access = private)
        function out = run_py(tc, keep_boundary)
            script = fullfile(tc.Tmp, 'air_test.py');
            res    = fullfile(tc.Tmp, 'res.txt');
            kb = 'False'; if keep_boundary; kb = 'True'; end
            fid = fopen(script, 'w');
            fprintf(fid, '%s\n', ...
                'import sys, os, numpy as np, nibabel as nib', ...
                sprintf('sys.path.insert(0, r"%s")', tc.Runner), ...
                'import run_babelbrain_case as r', ...
                sprintf('base = r"%s"', tc.Tmp), ...
                'n = 60; x, y, z = np.mgrid[:n, :n, :n]', ...
                'ct = np.full((n, n, n), -1000.0)', ...
                'ct[(x-30)**2 + (y-30)**2 + (z-34)**2 < 24**2] = 42', ...
                'pocket = (x-30)**2 + (y-22)**2 + (z-34)**2 < 4**2', ...
                'ct[pocket] = -950', ...
                'ct[:, 28:33, :] = 42', ...
                'ct[:4, :4, :4] = np.nan', ...
                'nib.Nifti1Image(ct, np.eye(4)).to_filename(os.path.join(base, "sub", "T_250kHz_6PPW_CT.nii.gz"))', ...
                sprintf('path, kept = r._write_air_regions(base, "sub", "T", 250e3, 6, keep_boundary=%s)', kb), ...
                'm = nib.load(path).get_fdata() > 0', ...
                sprintf('open(r"%s", "w").write("%%d %%d %%d %%d" %% (int(m[pocket].sum()), int(m[~pocket].sum()), int(pocket.sum()), int(os.path.isfile(path))))', res));
            fclose(fid);
            [st, msg] = system(sprintf('"%s" "%s"', tc.PyBin, script));
            tc.assertEqual(st, 0, msg);
            v = sscanf(fileread(res), '%d');
            out = struct('in_pocket', v(1), 'outside_pocket', v(2), 'pocket_vox', v(3), 'file_exists', v(4) == 1);
        end
    end
end
