classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_extract_real_intensity_profile < matlab.unittest.TestCase
% TEST_EXTRACT_REAL_INTENSITY_PROFILE  EXTRACT_REAL_INTENSITY_PROFILE uses
% calibration.skip_front_peak_mm when set and its 10 mm default otherwise, so
% the near-field peak is excluded from max_intens either way.
%
%   Run with:   results = runtests('tests/test_extract_real_intensity_profile.m');

    methods (Test)

        function test_skip_distance_default_and_explicit(tc)
            out = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            z = (0:0.5:80)';
            I = 20 * exp(-((z - 50) / 6).^2) + 40 * exp(-((z - 3) / 1).^2);   % near-field peak at 3 mm
            p.calibration = struct('path_output_profiles', out);
            [prof, max_i] = extract_real_intensity_profile(p, 50, 50, I, 'TEST', z);
            tc.verifyEqual(prof(:), I, 'AbsTol', 1e-9);
            tc.verifyEqual(max_i, 20, 'RelTol', 1e-6);              % default 10 mm skip
            p.calibration.skip_front_peak_mm = 0;
            [~, max_i] = extract_real_intensity_profile(p, 50, 50, I, 'TEST', z);
            tc.verifyEqual(max_i, max(I), 'RelTol', 1e-6);          % no skip: near-field peak
        end

    end
end
