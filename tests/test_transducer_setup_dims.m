classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_transducer_setup_dims < matlab.unittest.TestCase
% TEST_TRANSDUCER_SETUP_DIMS  Annular transducer masks on 2D and 3D grids.
%
% TRANSDUCER_SETUP builds the annular mask with makeBowl on a 3D grid and
% makeArc on a 2D grid (as SOURCE_CREATE does for the source). Every
% combination of grid dimensionality and central hole is compared with the
% k-Wave primitive it should use; the hole must remove voxels (negative
% control). No k-Wave run, no subject data.
%
%   results = runtests('tests/test_transducer_setup_dims.m');

    properties (TestParameter)
        n_dims = struct('two_d', 2, 'three_d', 3);
        id_mm  = struct('no_hole', 0, 'hole', 4);
    end

    methods (Test)

        function test_mask_matches_primitive(tc, n_dims, id_mm)
            [tr, grid_dims, trans_pos, focus_pos, res] = setup_case(n_dims, id_mm);

            [mask, labels] = transducer_setup(tr, trans_pos, focus_pos, grid_dims, res, struct());

            tc.verifySize(mask, grid_dims);
            tc.verifyGreaterThan(nnz(mask), 0);
            tc.verifyEqual(double(mask > 0), double(labels > 0));
            tc.verifyEqual(double(mask > 0), double(expected_mask(n_dims, tr, grid_dims, trans_pos, focus_pos, res) > 0));
        end

        function test_hole_removes_voxels(tc, n_dims)
            [tr0, grid_dims, trans_pos, focus_pos, res] = setup_case(n_dims, 0);
            tr1 = setup_case(n_dims, 4);
            full  = transducer_setup(tr0, trans_pos, focus_pos, grid_dims, res, struct());
            holed = transducer_setup(tr1, trans_pos, focus_pos, grid_dims, res, struct());
            tc.verifyLessThan(nnz(holed), nnz(full));
            tc.verifyTrue(all(holed(:) <= full(:)));
        end

    end
end

function [tr, grid_dims, trans_pos, focus_pos, res] = setup_case(n_dims, id_mm)
res = 0.5;                                            % mm per voxel
tr = struct('type', 'annular', 'annular', struct( ...
    'elem_n', 1, 'elem_od_mm', 16, 'elem_id_mm', id_mm, 'curv_radius_mm', 15));
if n_dims == 2
    grid_dims = [48 48];
    trans_pos = [24 4];
    focus_pos = [24 40];
else
    grid_dims = [40 40 48];
    trans_pos = [20 20 4];
    focus_pos = [20 20 40];
end
end

function m = expected_mask(n_dims, tr, grid_dims, trans_pos, focus_pos, res)
od = 2*floor(tr.annular.elem_od_mm / res / 2) + 1;
id = 2*floor(tr.annular.elem_id_mm / res / 2) + 1;
r  = round(tr.annular.curv_radius_mm / res);
if n_dims == 3, f = @makeBowl; else, f = @makeArc; end
m = f(grid_dims, trans_pos, r, od, focus_pos);
if tr.annular.elem_id_mm > 0
    m = m - f(grid_dims, trans_pos, r, id, focus_pos);
end
end
