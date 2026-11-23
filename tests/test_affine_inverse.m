classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_affine_inverse < matlab.unittest.TestCase
% TEST_AFFINE_INVERSE  4x4 affine inversion in both conventions.
%
% AFFINE_INVERSE is the only 4x4 affine inversion in PRESTUS. Checked on
% axis-aligned, oblique, and flipped affines in row-vector and column-vector
% form: M * Minv is the identity, the fixed row or column is exactly
% [0 0 0 1], and a maketform struct is read. Negative controls: a projective
% matrix, a singular linear part, and a non-4x4 input are refused. A source
% scan keeps inv() out of functions/.
%
%   results = runtests('tests/test_affine_inverse.m');

    properties (TestParameter)
        affine = struct( ...
            'identity', eye(4), ...
            'scaled_shifted', [0.5 0 0 0; 0 0.75 0 0; 0 0 1.25 0; -90 -120 -70 1], ...
            'oblique', oblique_affine(), ...
            'flipped_x', [-1 0 0 0; 0 1 0 0; 0 0 1 0; 90 -126 -72 1]);
    end

    methods (Test)

        function test_row_vector(tc, affine)
            Minv = affine_inverse(affine);
            tc.verifyEqual(affine * Minv, eye(4), 'AbsTol', 1e-12);
            tc.verifyEqual(Minv(:, 4), [0; 0; 0; 1]);          % exact, no residue
            pts = [1 2 3; -40 15 70];
            tc.verifyEqual(affine_apply_row(affine_apply_row(pts, affine), Minv), pts, 'AbsTol', 1e-9);
        end

        function test_column_vector(tc, affine)
            C = affine';
            Cinv = affine_inverse(C);
            tc.verifyEqual(C * Cinv, eye(4), 'AbsTol', 1e-12);
            tc.verifyEqual(Cinv(4, :), [0 0 0 1]);
            tc.verifyEqual(Cinv, affine_inverse(affine)', 'AbsTol', 1e-12);
        end

        function test_maketform_struct(tc, affine)
            s = struct('tdata', struct('T', affine));
            tc.verifyEqual(affine_inverse(s), affine_inverse(affine));
        end

        function test_matches_inv(tc, affine)
            tc.verifyEqual(affine_inverse(affine), inv(affine), 'AbsTol', 1e-10);
        end

        function test_rejects_projective(tc)
            M = eye(4); M(1, 4) = 0.1; M(4, 1) = 2;
            tc.verifyError(@() affine_inverse(M), 'PRESTUS:transform:badAffine');
        end

        function test_rejects_singular(tc)
            M = eye(4); M(3, 3) = 0;
            tc.verifyError(@() affine_inverse(M), 'PRESTUS:transform:singularAffine');
        end

        function test_no_inv_in_functions(tc)
            % Every 4x4 inversion goes through affine_inverse: no inv() call
            % in functions/ outside comments.
            root = fileparts(fileparts(which('affine_inverse')));
            files = dir(fullfile(root, '**', '*.m'));
            hits = {};
            for i = 1:numel(files)
                f = fullfile(files(i).folder, files(i).name);
                if has_inv_call(fileread(f)); hits{end+1} = f; end %#ok<AGROW>
            end
            tc.verifyEmpty(hits, 'Use affine_inverse instead of inv().');
            % negative control: the scan sees a call and ignores a comment
            tc.verifyTrue(has_inv_call(sprintf('x = 1;\nMinv = inv(M);\n')));
            tc.verifyFalse(has_inv_call(sprintf('%% inv(M) in a comment\nMinv = affine_inverse(M);\n')));
        end

        function test_rejects_wrong_size(tc)
            tc.verifyError(@() affine_inverse(eye(3)), 'PRESTUS:transform:badAffine');
            tc.verifyError(@() affine_inverse([eye(3) zeros(3, 1); 0 0 0 NaN]), 'PRESTUS:transform:badAffine');
        end

    end

end

function tf = has_inv_call(src)
    code = regexprep(src, '%[^\n]*', '');
    tf = ~isempty(regexp(code, '(?<![\w.])inv\s*\(', 'once'));
end

function out = affine_apply_row(pts, M)
    h = [pts, ones(size(pts, 1), 1)] * M;
    out = h(:, 1:3);
end

function T = oblique_affine()
% 10 and 20 degree rotations about z and x, 0.8 mm voxels, an offset
    a = deg2rad(10); b = deg2rad(20);
    Rz = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
    Rx = [1 0 0; 0 cos(b) -sin(b); 0 sin(b) cos(b)];
    M = Rz * Rx * diag([0.8 0.8 0.8]);
    T = [M' zeros(3, 1); -80 -110 -60 1];
end
