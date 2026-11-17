classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_extract_element_pos < matlab.unittest.TestCase
% TEST_EXTRACT_ELEMENT_POS  Matrix element positions read from file (EXTRACT_ELEMENT_POS).
%
% Elements of a curved array lie on a sphere of radius curv_radius_mm. With
% project_on_new_ROC, the rim stays in place and the rim-to-apex depth
% becomes the sagitta of the new radius at the rim radius (all in mm).
%
%   results = runtests('tests/test_extract_element_pos.m');

    methods (Test, TestTags = {'extract_element_pos'})

        function test_projection_depth_is_new_sagitta(tc)
            root = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            roc = 60;
            for new_roc = [45 80 120]
                for r_rim = [10 30]
                    [tr, n] = make_transducer(root, roc, r_rim);
                    tr.matrix.matrix_shape.extract_from_file.project_on_new_ROC = true;
                    tr.matrix.matrix_shape.extract_from_file.ROC_projection.new_ROC_mm = new_roc;
                    [~, pos] = evalc('extract_element_pos(struct(), tr, [0; 0; 0])');
                    tc.verifySize(pos, [3 n]);
                    z_mm = 1000 * pos(3, :);
                    expected = new_roc - sqrt(new_roc^2 - r_rim^2);
                    tc.verifyEqual(max(z_mm) - min(z_mm), expected, 'AbsTol', 1e-9, ...
                        sprintf('ROC %g -> %g mm, rim radius %g mm', roc, new_roc, r_rim));
                    % The rim stays where it was
                    tc.verifyEqual(max(z_mm), roc - sqrt(roc^2 - r_rim^2), 'AbsTol', 1e-9);
                end
            end
        end

        function test_no_projection_keeps_sphere(tc)
            % Negative control: without projection the depth is the
            % sagitta of curv_radius_mm, apex at z = 0.
            root = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            [tr, ~] = make_transducer(root, 60, 30);
            pos = extract_element_pos(struct(), tr, [0.01; 0.02; 0.03]);
            z_mm = 1000 * (pos(3, :) - 0.03);
            tc.verifyEqual(min(z_mm), 0, 'AbsTol', 1e-9);
            tc.verifyEqual(max(z_mm), 60 - sqrt(60^2 - 30^2), 'AbsTol', 1e-9);
            tc.verifyEqual(pos(1, 1), 0.01, 'AbsTol', 1e-12);    % apex element at the transducer position
        end

    end
end

% =========================================================================

function [tr, n] = make_transducer(root, roc, r_rim)
    % Elements on a sphere (origin at its centre, z towards the array):
    % the apex and three rings out to r_rim [mm].
    xyz = [0 0 roc];
    for r = r_rim * [1/3 2/3 1]
        a = (0:5)' * pi / 3;
        xyz = [xyz; r * cos(a), r * sin(a), repmat(sqrt(roc^2 - r^2), 6, 1)]; %#ok<AGROW>
    end
    n = size(xyz, 1);
    f = fullfile(root, sprintf('elements_%g_%g.csv', roc, r_rim));
    writetable(array2table(xyz, 'VariableNames', {'x', 'y', 'z'}), f);
    tr = struct();
    tr.matrix.curv_radius_mm = roc;
    tr.matrix.is_curved = true;
    tr.matrix.matrix_shape.extract_from_file = struct('file_path', f, 'start_row', 1, ...
        'elem_n', n, 'start_col', 1, 'select_random_subset', false, ...
        'project_on_new_ROC', false);
end
