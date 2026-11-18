classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_placement_export < matlab.unittest.TestCase
% TEST_PLACEMENT_EXPORT  Placed poses written for Brainsight, 3D Slicer, and k-Plan.
%
% Each exported file is read back with the conventions of the tool's reader:
% Brainsight and Slicer as BabelBrain (ConvMatTransform: location = target,
% transducer along +z; itk_to_BSight for Slicer), k-Plan as k-Plan (world =
% 0-based voxel x voxel size of the exported images, apex in row 4 [m], beam
% along -z). The fixture T1 is oblique with anisotropic voxels. The Nexus
% reads the same files in nexus/tests/test_placement.py.
%
%   results = runtests('tests/test_placement_export.m');

    properties
        Dir
        T1file
        A           % 0-based voxel (column) -> RAS
        Apex = [-70 -10 20];
        Focus = [-10 -5 15];
    end

    methods (TestMethodSetup)
        function make_t1(tc)
            tc.Dir = tempname;
            mkdir(tc.Dir);
            tc.addTeardown(@() rmdir(tc.Dir, 's'));
            img = single(reshape(1:30*32*28, [30 32 28]));
            tc.T1file = fullfile(tc.Dir, 't1.nii');
            niftiwrite(img, tc.T1file);
            info = niftiinfo(tc.T1file);
            pixdim = [0.9 1.0 1.2];
            a = deg2rad(5);
            R = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
            T = eye(4);
            T(1:3, 1:3) = diag(pixdim) * R';            % row-vector convention
            T(4, 1:3) = [-90 -20 -10];
            info.PixelDimensions = pixdim;
            info.Transform = affine3d(T);
            info.TransformName = 'Sform';
            niftiwrite(img, tc.T1file, info);
            tc.A = T';
        end
    end

    methods (Test)

        function test_slicer_reads_back_as_babelbrain(tc)
            f = export_slicer_transform(tc.Apex, tc.Focus, fullfile(tc.Dir, 'p.tfm'));
            txt = fileread(f);
            tok = regexp(txt, 'Parameters:([^\n]*)', 'tokens', 'once');
            p = sscanf(tok{1}, '%f')';
            itk = eye(4);
            itk(1:3, 1:3) = reshape(p(1:9), 3, 3)';
            itk(1:3, 4) = p(10:12)';
            L = diag([-1 -1 1 1]);
            B = inv(L * itk * L);
            B(1:3, 1:3) = -B(1:3, 1:3);                  % itk_to_BSight
            tc.verifyEqual(B(1:3, 4)', tc.Focus, 'AbsTol', 1e-6);
            tc.verifyEqual(B(1:3, 3)', unit(tc.Apex - tc.Focus), 'AbsTol', 1e-6);
        end

        function test_brainsight_reads_back_as_babelbrain(tc)
            f = export_brainsight_trajectory(tc.Apex, tc.Focus, 'p', fullfile(tc.Dir, 'p.txt'));
            rows = splitlines(strtrim(fileread(f)));
            row = rows{find(~startsWith(rows, '#'), 1)};
            v = str2double(strsplit(row, '\t'));
            tc.verifyEqual(v(2:4), tc.Focus, 'AbsTol', 1e-4);
            tc.verifyEqual(v(11:13), unit(tc.Apex - tc.Focus), 'AbsTol', 1e-5);   % m2 column
        end

        function test_brainsight_label_follows_sform_code(tc)
            f = export_brainsight_trajectory(tc.Apex, tc.Focus, 'p', fullfile(tc.Dir, 'q.txt'), t1_file=tc.T1file);
            codes = {'Scanner', 'Aligned', 'Talairach', 'MNI'};
            tc.verifySubstring(fileread(f), ['# Coordinate system: NIfTI:S:' ...
                codes{niftiinfo(tc.T1file).raw.sform_code}]);
        end

        function test_kplan_pose_in_kplan_world(tc)
            ct = fullfile(tc.Dir, 'ct.nii');
            copyfile(tc.T1file, ct);
            files = export_kplan_plan(tc.T1file, tc.Apex, tc.Focus, fullfile(tc.Dir, 'kp'), 'sub-001', ...
                ct_file=ct, kplan_offset_mm=10.82, transducer_name='TX');
            zooms = sqrt(sum(tc.A(1:3, 1:3).^2, 1));
            world = @(r) ((tc.A \ [r(:); 1])' * [eye(3); 0 0 0]) .* zooms;
            T = reshape(double(h5read(files.kps, '/1/position_transform')), 4, 4)';
            tc.verifyEqual(T(4, 1:3) * 1e3, world(tc.Apex), 'AbsTol', 1e-3);
            tc.verifyEqual(T(3, 1:3), unit(world(tc.Apex) - world(tc.Focus)), 'AbsTol', 1e-5);
            tc.verifyEqual(det(T(1:3, 1:3)), 1, 'AbsTol', 1e-5);                       % right-handed axes
            for f = {files.t1, files.ct}
                info = niftiinfo(f{1});
                tc.verifyEqual(info.Transform.T(1:3, 1:3), diag(zooms), 'AbsTol', 1e-5);  % rotation removed
                tc.verifyEqual(info.Transform.T(4, 1:3), [0 0 0], 'AbsTol', 1e-6);         % voxel 0 at world 0
                tc.verifyEqual(niftiread(f{1}), niftiread(tc.T1file));                     % same voxels
            end
            entry = fileread(files.entry);
            d = norm(world(tc.Focus) - world(tc.Apex));
            tc.verifySubstring(entry, sprintf('%.2f mm from k-Plan''s exit plane', d - 10.82));
            tc.verifySubstring(entry, 'sub-001_transducer.kps');
        end

        function test_kplan_rejects_non_ras_and_off_grid_ct(tc)
            img = niftiread(tc.T1file);
            info = niftiinfo(tc.T1file);
            T = info.Transform.T;  T(1, :) = -T(1, :);                 % flip the first voxel axis
            info.Transform = affine3d(T);
            flipped = fullfile(tc.Dir, 'flipped.nii');
            niftiwrite(img, flipped, info);
            tc.verifyError(@() export_kplan_plan(flipped, tc.Apex, tc.Focus, fullfile(tc.Dir, 'k1'), 's'), ...
                'PRESTUS:kplan:notRAS');
            small = fullfile(tc.Dir, 'small.nii');
            niftiwrite(img(1:10, :, :), small, setfield(niftiinfo(tc.T1file), 'ImageSize', [10 32 28])); %#ok<SFLD>
            files = tc.verifyWarning(@() export_kplan_plan(tc.T1file, tc.Apex, tc.Focus, ...
                fullfile(tc.Dir, 'k2'), 's', ct_file=small), 'PRESTUS:kplan:ctGrid');
            tc.verifyEmpty(files.ct);
        end

        function test_export_placement_writes_each_format(tc)
            p = struct();
            p.subject_id = 1;
            p.path.anat = tc.Dir;
            p.path.t1_pattern = 't1.nii';
            p.io.dir_output = fullfile(tc.Dir, 'out');
            p.io.output_affix = '_x';
            p.placement.export_formats = {'brainsight', 'slicer', 'kplan', 'bogus'};
            p.transducer = struct('trans_pos', {[3 16 14], [5 16 14]}, 'focus_pos', {[20 17 14], [20 18 14]});
            files = tc.verifyWarning(@() export_placement(p), 'prestus:placement:exportFormat');
            tc.verifyNumElements(files, 6);                           % 3 formats x 2 transducers
            tc.verifyTrue(isfile(fullfile(p.io.dir_output, 'sub-001_placement_x_tx2_slicer.tfm')));
            tc.verifyTrue(isfolder(fullfile(p.io.dir_output, 'sub-001_placement_x_tx1_kplan')));
            p.placement.export_formats = {};
            tc.verifyEmpty(export_placement(p));
        end
    end
end

function u = unit(v)
u = v / norm(v);
end
