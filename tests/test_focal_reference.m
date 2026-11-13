classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_focal_reference < matlab.unittest.TestCase
% TEST_FOCAL_REFERENCE  Unit tests for the calibration.focal_reference convention
% (FOCAL_REFERENCE_TO_FRAC, COMPUTE_FLHM_CENTER, COMPUTE_FOCAL_METRICS).
%
%   Run with:   results = runtests('tests/test_focal_reference.m');

    % ------------------------------------------------------------------ %
    %% focal_reference_to_frac
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'focal_reference'})

        function test_frac_3db(tc)
            tc.verifyEqual(focal_reference_to_frac('flhm_3db'), 0.5);
        end

        function test_frac_6db(tc)
            tc.verifyEqual(focal_reference_to_frac('flhm_6db'), 0.25);
        end

        function test_frac_peak_errors(tc)
            tc.verifyError(@() focal_reference_to_frac('peak'), ...
                'focal_reference_to_frac:peakHasNoFrac');
        end

        function test_frac_unknown_errors(tc)
            tc.verifyError(@() focal_reference_to_frac('fwhm'), ...
                'focal_reference_to_frac:unknownReference');
        end

    end

    % ------------------------------------------------------------------ %
    %% compute_flhm_center
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'focal_reference'})

        function test_flhm_center_default_is_3db(tc)
            [z, I] = gaussian_profile(50, 5);
            [c_def, lo_def, hi_def] = compute_flhm_center(I, z);
            [c_3db, lo_3db, hi_3db] = compute_flhm_center(I, z, 'flhm_3db');
            tc.verifyEqual([c_def lo_def hi_def], [c_3db lo_3db hi_3db]);
        end

        function test_flhm_span_matches_threshold(tc)
            sigma = 5;
            [z, I] = gaussian_profile(50, sigma);
            [c3, lo3, hi3] = compute_flhm_center(I, z, 'flhm_3db');
            [c6, lo6, hi6] = compute_flhm_center(I, z, 'flhm_6db');
            % Symmetric profile: both centres at the peak
            tc.verifyEqual(c3, 50, 'AbsTol', 0.1);
            tc.verifyEqual(c6, 50, 'AbsTol', 0.1);
            % Analytic full widths: 2*sigma*sqrt(2*ln(1/frac))
            tc.verifyEqual(hi3 - lo3, 2*sigma*sqrt(2*log(2)), 'AbsTol', 0.25);
            tc.verifyEqual(hi6 - lo6, 2*sigma*sqrt(2*log(4)), 'AbsTol', 0.25);
        end

        function test_flhm_centre_shifts_with_threshold_on_skewed_profile(tc)
            % Distally skewed lobe (as for a real focus): the lower threshold
            % reaches further into the long tail, so the centre moves distal.
            z = (0:0.1:100)';
            I = exp(-(z - 40).^2 / (2*4^2)) .* (z <= 40) + ...
                exp(-(z - 40).^2 / (2*10^2)) .* (z > 40);
            c3 = compute_flhm_center(I, z, 'flhm_3db');
            c6 = compute_flhm_center(I, z, 'flhm_6db');
            tc.verifyGreaterThan(c6, c3);
        end

        function test_flhm_center_peak_returns_peak_position(tc)
            z = (0:0.1:100)';
            I = exp(-(z - 40).^2 / (2*4^2)) .* (z <= 40) + ...
                exp(-(z - 40).^2 / (2*10^2)) .* (z > 40);
            [c, lo, hi] = compute_flhm_center(I, z, 'peak');
            [~, lo3, hi3] = compute_flhm_center(I, z, 'flhm_3db');
            tc.verifyEqual(c, 40, 'AbsTol', 1e-9);
            tc.verifyEqual([lo hi], [lo3 hi3]);   % edges fall back to -3 dB
        end

    end

    % ------------------------------------------------------------------ %
    %% compute_focal_metrics
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'focal_reference'})

        function test_focal_metrics_6db_span_wider(tc)
            [Ipa, mask] = gaussian_volume();
            fm3 = compute_focal_metrics(Ipa, max(Ipa(:)), mask, [11 11 1], 1, 'flhm_3db');
            fm6 = compute_focal_metrics(Ipa, max(Ipa(:)), mask, [11 11 1], 1, 'flhm_6db');
            tc.verifyGreaterThan(fm6.flhm_span_mm, fm3.flhm_span_mm);
            % ITRUSST volumes are independent of focal_reference
            tc.verifyEqual(fm6.vol_3dB_mm3, fm3.vol_3dB_mm3);
            tc.verifyEqual(fm6.vol_6dB_mm3, fm3.vol_6dB_mm3);
        end

        function test_focal_metrics_accepts_peak(tc)
            % 'peak' is a documented calibration.focal_reference value and is
            % passed through unchanged by acoustic_analysis, so it must not
            % error; the FLHM centre should then be the peak itself.
            [Ipa, mask] = gaussian_volume();
            fm  = compute_focal_metrics(Ipa, max(Ipa(:)), mask, [11 11 1], 1, 'peak');
            fm3 = compute_focal_metrics(Ipa, max(Ipa(:)), mask, [11 11 1], 1, 'flhm_3db');
            tc.verifyEqual(fm.flhm_center_vox, fm.peak_pos_vox);
            tc.verifyEqual(fm.flhm_span_mm, fm3.flhm_span_mm);   % span at -3 dB
        end

    end


    % plot_axial_profile: beam direction
    methods (Test, TestTags = {'focal_reference'})

        function test_minus_z_beam_gives_positive_depths(tc)
            tmp = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            % Transducer on the high-z face, beam travelling towards low z.
            [p, Ipa, P, masks] = axial_profile_inputs(tmp);
            exitplane_pos = [5 5 36];
            focus_pos     = [5 5 16];
            plot_axial_profile(Ipa, P, masks, p, exitplane_pos, focus_pos, focus_pos);

            T = readtable(axial_csv(p));
            d = T.depth_from_exit_plane_mm;
            tc.verifyEqual(d(16), 20 * p.grid.resolution_mm, 'AbsTol', 1e-9);   % focus 20 vox deep
            tc.verifyEqual(d(36), 0, 'AbsTol', 1e-9);                          % exit plane
            tc.verifyLessThan(d(40), 0);                                       % behind transducer
        end

        function test_plus_z_beam_unchanged(tc)
            tmp = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            [p, Ipa, P, masks] = axial_profile_inputs(tmp);
            exitplane_pos = [5 5 5];
            focus_pos     = [5 5 25];
            plot_axial_profile(Ipa, P, masks, p, exitplane_pos, focus_pos, focus_pos);

            T = readtable(axial_csv(p));
            d = T.depth_from_exit_plane_mm;
            tc.verifyEqual(d, ((1:40)' - 5) * p.grid.resolution_mm, 'AbsTol', 1e-9);
        end

    end
end

% ====================================================================== %
%% Local helpers
% ====================================================================== %

function [z, I] = gaussian_profile(mu, sigma)
    z = (0:0.05:100)';
    I = exp(-(z - mu).^2 / (2*sigma^2));
end

function [Ipa, mask] = gaussian_volume()
    [x, y, z] = ndgrid(1:21, 1:21, 1:60);
    Ipa  = 10 * exp(-((x-11).^2 + (y-11).^2) / (2*2^2) - (z-30).^2 / (2*6^2));
    mask = true(size(Ipa));
end

function [p, Ipa, P, masks] = axial_profile_inputs(out_dir)
    p = struct();
    p.grid.resolution_mm = 0.5;
    p.medium_properties = struct('water', struct(), 'brain', struct());
    p.io.dir_img       = out_dir;
    p.io.output_affix  = '_test';
    p.subject_id       = 1;
    p.simulation.medium = 'water';
    Ipa   = rand(9, 9, 40);
    P     = rand(9, 9, 40) * 1e5;
    masks = ones(9, 9, 40);
end

function f = axial_csv(p)
    f = fullfile(p.io.dir_img, sprintf('sub-%03d_%s_axial_profile%s.csv', ...
        p.subject_id, p.simulation.medium, p.io.output_affix));
end
