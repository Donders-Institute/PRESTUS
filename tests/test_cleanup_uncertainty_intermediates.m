classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_cleanup_uncertainty_intermediates < matlab.unittest.TestCase
% TEST_CLEANUP_UNCERTAINTY_INTERMEDIATES  Cleanup removes preprocessing caches and keeps heating results.
%
%   Heating results (cache/*_heating_res*.mat) hold thermal results that
%   later analyses read; CLEANUP_UNCERTAINTY_INTERMEDIATES must keep them,
%   whatever io.save_thermal_matrices says.
%
%   results = runtests('tests/test_cleanup_uncertainty_intermediates.m');

    methods (Test, TestTags = {'uncertainty'})

        function test_keeps_heating_results_removes_preproc_caches(tc)
            root = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            for save_thermal = {0, 1, []}
                p = struct('subject_id', 1, ...
                    'simulation', struct('medium', 'layered'), ...
                    'path', struct('sim', root), ...
                    'io', struct('save_matrices', 0, 'save_thermal_matrices', save_thermal{1}));
                cache_dir = fullfile(root, 'sub-001', 'cache');
                if ~isfolder(cache_dir); mkdir(cache_dir); end

                % Heating results: thermal_stage name (sub-NNN_<medium>_heating_res<affix>.mat)
                % for several variant affixes, plus the affix-before-tag form
                kept = {'sub-001_layered_heating_res.mat', ...
                        'sub-001_layered_heating_res_default.mat', ...
                        'sub-001_layered_heating_res_Temp.mat', ...
                        'sub-001_layered_default_heating_res.mat'};
                removed = {'sub-001_layered_rotated_scaled.mat', ...
                           'sub-001_layered_cropped_smoothed.mat'};
                cellfun(@(f) touch(fullfile(cache_dir, f)), [kept, removed]);

                evalc('cleanup_uncertainty_intermediates(p)');

                for f = kept
                    tc.verifyTrue(isfile(fullfile(cache_dir, f{1})), ...
                        sprintf('heating result deleted: %s', f{1}));
                end
                for f = removed   % positive control: cleanup does run
                    tc.verifyFalse(isfile(fullfile(cache_dir, f{1})), ...
                        sprintf('intermediate kept: %s', f{1}));
                end
            end
        end

        function test_missing_files_are_not_an_error(tc)
            root = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            p = struct('subject_id', 2, 'simulation', struct('medium', 'water'), ...
                'path', struct('sim', root), 'io', struct('save_matrices', 0));
            tc.verifyWarningFree(@() cleanup_uncertainty_intermediates(p));
        end

    end

end

function touch(f)
    fid = fopen(f, 'w');
    fclose(fid);
end
