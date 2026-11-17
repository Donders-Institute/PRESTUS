classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_acoustic_fov_crop < matlab.unittest.TestCase
% TEST_ACOUSTIC_FOV_CROP  ACOUSTIC_FOV_CROP pads the jwave FOV by the PML that
% JWAVE_WRAPPER applies (RESOLVE_JWAVE_PML_SIZE on the cropped medium) and
% leaves the k-Wave FOV unpadded.
%
%   Run with:   results = runtests('tests/test_acoustic_fov_crop.m');

    properties (TestParameter)
        pml_size = {'auto', 12}
    end

    methods (Test)

        function test_jwave_padding_equals_jwave_pml(tc, pml_size)
            [p, med, masks, tp, fp] = tc.setup_grid('jwave', pml_size);
            [p_out, med_out] = acoustic_fov_crop(p, med, masks, [], tp, fp);
            pad_vox = (p_out.grid.acoustic_fov_diameter_mm - p.grid.acoustic_fov_diameter_mm) ...
                / 2 / p.grid.resolution_mm;
            pad_len = (p_out.grid.acoustic_fov_length_mm - p.grid.acoustic_fov_length_mm) ...
                / 2 / p.grid.resolution_mm;
            % What jwave_wrapper resolves for the cropped grid it receives
            jw = max(resolve_jwave_pml_size(p_out.grid.pml_size, p_out.transducer(1).freq_hz, ...
                max(med_out.sound_speed(:)), p_out.grid.resolution_mm / 1e3));
            tc.verifyEqual(pad_vox, jw, 'AbsTol', 1e-9);
            tc.verifyEqual(pad_len, jw, 'AbsTol', 1e-9);
        end

        function test_kwave_fov_not_padded(tc, pml_size)
            [p, med, masks, tp, fp] = tc.setup_grid('kwave', pml_size);
            p_out = acoustic_fov_crop(p, med, masks, [], tp, fp);
            tc.verifyEqual(p_out.grid.acoustic_fov_diameter_mm, p.grid.acoustic_fov_diameter_mm);
            tc.verifyEqual(p_out.grid.acoustic_fov_length_mm,   p.grid.acoustic_fov_length_mm);
        end

    end

    methods (Static)
        function [p, med, masks, tp, fp] = setup_grid(backend, pml_size)
            dims = [80 80 120];
            tp = [40 40 5];
            fp = [40 40 60];
            p = struct();
            p.acoustic_backend = backend;
            p.grid.dims = dims;
            p.grid.resolution_mm = 1;
            p.grid.pml_size = pml_size;
            p.grid.acoustic_fov_diameter_mm = 40;
            p.grid.acoustic_fov_length_mm = 80;
            p.transducer(1).freq_hz = 500e3;
            p.transducer(1).annular.elem_od_mm = 30;
            p.transducer(1).trans_pos = tp;
            p.transducer(1).focus_pos = fp;
            % Fastest voxels (skull-like) lie inside the FOV, near the focus
            med.sound_speed = 1500 * ones(dims);
            med.sound_speed(35:45, 35:45, 30:35) = 3000;
            med.density = 1000 * ones(dims);
            masks = ones(dims, 'uint8');
        end
    end
end
