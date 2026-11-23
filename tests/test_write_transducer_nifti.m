classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_write_transducer_nifti < matlab.unittest.TestCase
% TEST_WRITE_TRANSDUCER_NIFTI  Transducer model NIfTI on the voxel grid of a T1.
%
% WRITE_TRANSDUCER_NIFTI builds an annular bowl with makeBowl and writes it
% on the T1 grid: 1 = emitting face, 0.5 = display thickness behind it,
% 0 = elsewhere, including the central opening. Each test runs on an
% isotropic and an anisotropic synthetic T1. Checks: the output keeps the
% T1 grid; face voxels lie within the aperture; the opening is empty (with a
% solid bowl as negative control); back_mm = 0 leaves the face unchanged;
% a bowl partly outside the field of view is clipped and reported. No
% k-Wave run, no subject data.
%
%   results = runtests('tests/test_write_transducer_nifti.m');

    properties (TestParameter)
        pixdim = struct('iso', [1 1 1], 'aniso', [1 1 2]);
    end

    methods (Test)

        function test_output_keeps_t1_grid(tc, pixdim)
            [hdr, params, dir_out] = setup_case(tc, pixdim, 8);
            out = write_and_read(params, hdr, dir_out, 3);
            info = niftiinfo(fullfile(dir_out, 'transducer.nii.gz'));
            tc.verifyEqual(double(info.ImageSize(1:3)), double(hdr.ImageSize(1:3)));
            tc.verifyEqual(double(info.PixelDimensions(1:3)), double(hdr.PixelDimensions(1:3)));
            tc.verifyEqual(info.Transform.T, hdr.Transform.T, 'AbsTol', 1e-5);
            tc.verifyClass(out, 'single');
            tc.verifyTrue(all(ismember(unique(out(:)), single([0 0.5 1]))));
            tc.verifyGreaterThan(nnz(out == 1), 0);
            tc.verifyGreaterThan(nnz(out == 0.5), 0);
        end

        function test_face_within_aperture(tc, pixdim)
            [hdr, params, dir_out] = setup_case(tc, pixdim, 8);
            out = write_and_read(params, hdr, dir_out, 3);
            tr  = params.transducer;
            rho = axis_distance_mm(find(out == 1), size(out), tr, pixdim);
            tol = norm(pixdim) / 2 + min(pixdim);     % voxel mapping + bowl voxelisation
            tc.verifyLessThanOrEqual(max(rho), tr.annular.elem_od_mm / 2 + tol);
            tc.verifyGreaterThanOrEqual(min(rho), tr.annular.elem_id_mm / 2 - tol);
        end

        function test_opening_is_empty(tc, pixdim)
            % Bowl apex (trans_pos) lies in the central opening; a solid bowl
            % (no opening) fills it: negative control
            [hdr, params, dir_out] = setup_case(tc, pixdim, 8);
            holed = write_and_read(params, hdr, dir_out, 0);
            [~, params0] = setup_case(tc, pixdim, 0);
            solid = write_and_read(params0, hdr, dir_out, 0);
            p = num2cell(params.transducer.trans_pos);
            tc.verifyEqual(holed(p{:}), single(0));
            tc.verifyEqual(solid(p{:}), single(1));
            tc.verifyLessThan(nnz(holed), nnz(solid));
        end

        function test_back_thickness_keeps_face(tc, pixdim)
            [hdr, params, dir_out] = setup_case(tc, pixdim, 8);
            face_only = write_and_read(params, hdr, dir_out, 0);
            thick     = write_and_read(params, hdr, dir_out, 3);
            tc.verifyTrue(all(ismember(unique(face_only(:)), single([0 1]))));
            tc.verifyEqual(find(thick == 1), find(face_only == 1));
            tc.verifyGreaterThan(nnz(thick), nnz(face_only));
        end

        function test_partly_outside_fov_is_clipped(tc, pixdim)
            % Bowl centred near the x = 1 edge: about half lies outside the T1
            [hdr, params, dir_out] = setup_case(tc, pixdim, 8);
            params.transducer.trans_pos(1) = 3;
            params.transducer.focus_pos(1) = 3;
            log = evalc('out = write_and_read(params, hdr, dir_out, 3);');
            tc.verifySubstring(log, 'outside the T1 field of view');
            tc.verifySize(out, double(hdr.ImageSize(1:3)));
            tc.verifyGreaterThan(nnz(out == 1), 0);
        end

    end
end

function [hdr, params, dir_out] = setup_case(tc, pixdim, id_mm)
% Synthetic T1 (about 64 mm per side), beam along +z
dir_out = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
sz      = round(64 ./ pixdim);
t1_file = fullfile(dir_out, 't1.nii');
niftiwrite(zeros(sz, 'int16'), t1_file);
info = niftiinfo(t1_file);
info.PixelDimensions = pixdim;
info.Transform.T = diag([pixdim 1]);
info.Transform.T(4, 1:3) = -(sz / 2) .* pixdim;
niftiwrite(zeros(sz, 'int16'), t1_file, info);
hdr = niftiinfo(t1_file);

tr = struct('type', 'annular', 'annular', struct( ...
    'elem_n', 1, 'elem_od_mm', 30, 'elem_id_mm', id_mm, 'curv_radius_mm', 30));
tr.trans_pos = [round(sz(1) / 2), round(sz(2) / 2), round(6 / pixdim(3))];
tr.focus_pos = [tr.trans_pos(1:2), round(50 / pixdim(3))];
params = struct('transducer', tr);
end

function out = write_and_read(params, hdr, dir_out, back_mm)
out_file = fullfile(dir_out, 'transducer.nii.gz');
write_transducer_nifti(params, hdr, out_file, back_mm);
out = niftiread(out_file);
end

function rho = axis_distance_mm(idx, sz, tr, pixdim)
% Distance [mm] of voxel centres from the beam axis through trans_pos
[ii, jj, kk] = ind2sub(sz, idx);
X = ([ii jj kk] - tr.trans_pos) .* pixdim;
u = (tr.focus_pos - tr.trans_pos) .* pixdim;
u = u / norm(u);
rho = vecnorm(X - (X * u') * u, 2, 2);
end
