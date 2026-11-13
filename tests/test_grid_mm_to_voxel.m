classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_grid_mm_to_voxel < matlab.unittest.TestCase
% TEST_GRID_MM_TO_VOXEL  transducer.*_pos_mm is measured from the centre of voxel 1.
%
% Regression for trans_pos_mm = 0 mapping to the invalid voxel index 0.
% Pure logic: no k-Wave, no data.
%
%   results = runtests('tests/test_grid_mm_to_voxel.m');


    methods (Test)
        function test_zero_mm_is_voxel_one(tc)
            tc.verifyEqual(grid_mm_to_voxel([0 0], 0.5), [1 1]);
        end

        function test_voxel_spacing(tc)
            % ITRUSST layout: axis at 34.5 mm, focus 64 mm deep
            tc.verifyEqual(grid_mm_to_voxel([34.5 64], 0.5), [70 129]);
        end

        function test_column_input_gives_row(tc)
            tc.verifySize(grid_mm_to_voxel([1; 2; 3], 1), [1 3]);
        end
    end
end
