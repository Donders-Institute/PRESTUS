classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_show_3d_head < matlab.unittest.TestCase
% TEST_SHOW_3D_HEAD  show_3d_head renders every orientation branch on a synthetic head.
%
% Draws into an invisible figure: align_to_focus true, false (natural focus
% along +z), and absent; and a head without brain/skull labels, where isocaps
% returns no cap surface. Needs gpuDeviceCount (Parallel Computing Toolbox)
% and a graphics-capable session; skipped otherwise.
%
%   results = runtests('tests/test_show_3d_head.m');

    properties
        Fig
    end

    methods (TestMethodSetup)
        function open_hidden_figure(tc)
            tc.assumeTrue(exist('gpuDeviceCount', 'file') > 0, ...
                'show_3d_head calls gpuDeviceCount (Parallel Computing Toolbox).');
            tc.Fig = figure('Visible', 'off');
            tc.addTeardown(@() close(tc.Fig));
        end
    end

    methods (Test, TestTags = {'plot'})

        function test_align_to_focus_branches(tc)
            % align_to_focus = false used an undefined trans_pos_grid and errored
            for align = {true, false, []}
                [seg, targ, trans, p, mesh] = tc.scene(true);
                if isempty(align{1})
                    p.transducer = rmfield(p.transducer, 'align_to_focus');
                else
                    p.transducer.align_to_focus = align{1};
                end
                clf(tc.Fig);
                tc.verifyWarningFree(@() show_3d_head(seg, targ, trans, p, 1, mesh, ...
                    [0 0 0], [0 0], 0, false));
                tc.verifyNotEmpty(findobj(tc.Fig, 'Type', 'patch'));
            end
        end

        function test_head_without_cap_layers(tc)
            % skin only: no cap surface; rendering continues
            [seg, targ, trans, p, mesh] = tc.scene(false);
            clf(tc.Fig);
            show_3d_head(seg, targ, trans, p, 1, mesh, [0 0 0], [0 0], 0, false);
            tc.verifyNotEmpty(findobj(tc.Fig, 'Type', 'patch'));
        end
    end

    methods (Access = private)
        function [seg, targ, trans, p, mesh] = scene(~, with_brain)
            dims = [48 48 48];
            seg = make_synthetic_segmentation(dims);
            if ~with_brain
                seg(seg > 0) = 5;                           % skin label only
            end
            targ  = [24 24 24];
            trans = [24 24 2];
            p.grid.resolution_mm = 1;
            p.transducer = struct('type', 'annular', 'align_to_focus', true, ...
                'annular', struct('elem_od_mm', [10 20], 'curv_radius_mm', 30));
            [x, y, z] = ndgrid(1:dims(1), 1:dims(2), 1:dims(3));
            mesh = [x(:) y(:) z(:)];
        end
    end
end
