classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_cache_signature < matlab.unittest.TestCase
% TEST_CACHE_SIGNATURE  Reused grid/medium caches must match the run's settings.
%
% Covers CACHE_SIGNATURE and CHECK_CACHE_SIGNATURE: which settings are part of
% a stage's signature, the error/warn/off modes of io.cache_signature_check,
% and legacy caches without a signature. Pure logic: no k-Wave, no data.
%
%   results = runtests('tests/test_cache_signature.m');


    methods (Static, Access = private)
        function p = params()
            p = struct();
            p.subject_id = 1;
            p.simulation = struct('medium', 'layered', 'debug', 0);
            p.grid = struct('resolution_mm', 0.5, 'default_dims', [128 128 256]);
            p.layers = struct('skull', 4);
            p.medium_properties = struct('skull', struct('sound_speed', 2800));
            p.transducer = struct('trans_pos', [10 20 30], 'focus_pos', [10 20 90], ...
                'freq_hz', 500e3, 'target_isppa_wcm2', 30);
            p.path = struct('sim', '/a');
            p.io = struct('overwrite_files', 'never', 'cache_signature_check', 'error');
        end
    end

    methods (Test, TestTags = {'cache_signature'})

        function test_identical_settings_match(tc)
            p = tc.params();
            ok = check_cache_signature('x.mat', cache_signature(p, 'grid'), ...
                cache_signature(p, 'grid'), p);
            tc.verifyTrue(ok);
        end

        function test_focal_axis_min_norm_legacy_cache(tc)
            % caches from before grid.focal_axis_min_norm was a default match
            % its built-in 0.5, but not another value
            p = tc.params();
            legacy = cache_signature(p, 'grid');
            p.grid.focal_axis_min_norm = 0.5;
            tc.verifyTrue(check_cache_signature('x.mat', legacy, cache_signature(p, 'grid'), p));
            p.grid.focal_axis_min_norm = 0.3;
            tc.verifyError(@() check_cache_signature('x.mat', legacy, cache_signature(p, 'grid'), p), ...
                'prestus:cache:signatureMismatch');
        end

        function test_changed_transducer_position_is_detected(tc)
            p = tc.params(); q = p;
            q.transducer.trans_pos = [10 20 31];
            tc.verifyError(@() check_cache_signature('x.mat', cache_signature(p, 'grid'), ...
                cache_signature(q, 'grid'), q), 'prestus:cache:signatureMismatch');
        end

        function test_changed_resolution_is_detected(tc)
            p = tc.params(); q = p;
            q.grid.resolution_mm = 1;
            tc.verifyError(@() check_cache_signature('x.mat', cache_signature(p, 'grid'), ...
                cache_signature(q, 'grid'), q), 'prestus:cache:signatureMismatch');
        end

        function test_drive_level_and_paths_do_not_invalidate(tc)
            % Target intensity scales the field, paths may move: same grid
            p = tc.params(); q = p;
            q.transducer.target_isppa_wcm2 = 60;
            q.path.sim = '/elsewhere';
            q.io.overwrite_files = 'always';
            tc.verifyTrue(check_cache_signature('x.mat', cache_signature(p, 'grid'), ...
                cache_signature(q, 'grid'), q));
        end

        function test_acoustic_fov_does_not_invalidate_grid_or_medium(tc)
            % The partial FOV is cropped after the grid and medium stages
            p = tc.params(); q = p;
            q.grid.acoustic_fov_diameter_mm = 'auto';
            q.grid.acoustic_fov_length_mm   = 'auto';
            for kind = {'grid', 'medium'}
                tc.verifyTrue(check_cache_signature('x.mat', cache_signature(p, kind{1}), ...
                    cache_signature(q, kind{1}), q));
            end
            % a cache written while the signature still recorded them
            old = cache_signature(p, 'grid');
            old.grid.acoustic_fov_length_mm = 'auto';
            tc.verifyTrue(check_cache_signature('x.mat', old, cache_signature(p, 'grid'), p));
            % a reused pressure map does depend on it
            tc.verifyError(@() check_cache_signature('x.mat', cache_signature(p, 'acoustic'), ...
                cache_signature(q, 'acoustic'), q), 'prestus:cache:signatureMismatch');
        end

        function test_medium_properties_only_in_medium_signature(tc)
            p = tc.params(); q = p;
            q.medium_properties.skull.sound_speed = 3000;
            tc.verifyTrue(check_cache_signature('x.mat', cache_signature(p, 'grid'), ...
                cache_signature(q, 'grid'), q));
            tc.verifyError(@() check_cache_signature('x.mat', cache_signature(p, 'medium'), ...
                cache_signature(q, 'medium'), q), 'prestus:cache:signatureMismatch');
        end

        function test_warn_mode_reuses_with_warning(tc)
            p = tc.params(); q = p;
            q.grid.resolution_mm = 1;
            q.io.cache_signature_check = 'warn';
            tc.verifyWarning(@() check_cache_signature('x.mat', cache_signature(p, 'grid'), ...
                cache_signature(q, 'grid'), q), 'prestus:cache:signatureMismatch');
        end

        function test_off_mode_skips(tc)
            p = tc.params(); q = p;
            q.grid.resolution_mm = 1;
            q.io.cache_signature_check = 'off';
            tc.verifyTrue(check_cache_signature('x.mat', cache_signature(p, 'grid'), ...
                cache_signature(q, 'grid'), q));
        end

        function test_unquoted_yaml_off_means_off(tc)
            % YAML 1.1 reads an unquoted 'off' as false
            p = tc.params(); q = p;
            q.grid.resolution_mm = 1;
            q.io.cache_signature_check = false;
            tc.verifyTrue(check_cache_signature('x.mat', cache_signature(p, 'grid'), ...
                cache_signature(q, 'grid'), q));
        end

        function test_legacy_cache_warns_unverified(tc)
            p = tc.params();
            tc.verifyWarning(@() check_cache_signature('x.mat', [], cache_signature(p, 'grid'), p), ...
                'prestus:cache:unverified');
        end

        function test_mismatch_message_names_keys(tc)
            p = tc.params(); q = p;
            q.grid.resolution_mm = 1;
            try
                check_cache_signature('x.mat', cache_signature(p, 'grid'), cache_signature(q, 'grid'), q);
                tc.verifyFail('expected an error');
            catch ME
                tc.verifySubstring(ME.message, 'grid.resolution_mm');
            end
        end

        function test_bad_mode_errors(tc)
            p = tc.params();
            p.io.cache_signature_check = 'maybe';
            tc.verifyError(@() check_cache_signature('x.mat', cache_signature(p, 'grid'), ...
                cache_signature(p, 'grid'), p), 'prestus:cache:badCheckMode');
        end

        function test_regenerated_phantom_is_detected(tc)
            % A phantom is read from final_tissues.nii.gz: a new file
            % invalidates grid and medium caches built from the old one
            root = tempname; mkdir(root);
            tc.addTeardown(@() rmdir(root, 's'));
            m2m = fullfile(root, 'm2m_sub-001'); mkdir(m2m);
            f = fullfile(m2m, 'final_tissues.nii.gz');
            fid = fopen(f, 'w'); fwrite(fid, uint8(1:10)); fclose(fid);
            p = tc.params();
            p.simulation.medium = 'phantom';
            p.path.seg = root;
            old = struct('grid', cache_signature(p, 'grid'), 'medium', cache_signature(p, 'medium'));
            fid = fopen(f, 'w'); fwrite(fid, uint8(1:20)); fclose(fid);
            for kind = {'grid', 'medium'}
                tc.verifyError(@() check_cache_signature('x.mat', old.(kind{1}), ...
                    cache_signature(p, kind{1}), p), 'prestus:cache:signatureMismatch');
            end
            % a cache written before the phantom file was recorded is not compared on it
            legacy = rmfield(cache_signature(p, 'grid'), 'phantom_file');
            tc.verifyTrue(check_cache_signature('x.mat', legacy, cache_signature(p, 'grid'), p));
            % layered runs do not record it
            p.simulation.medium = 'layered';
            tc.verifyFalse(isfield(cache_signature(p, 'grid'), 'phantom_file'));
        end

        function test_resegmented_head_is_detected(tc)
            % A layered head records the segmentation's content: a charm re-run
            % that rewrites final_tissues (same size here) invalidates every cache
            root = tempname; mkdir(root);
            tc.addTeardown(@() rmdir(root, 's'));
            m2m = fullfile(root, 'm2m_sub-001'); mkdir(m2m);
            f = fullfile(m2m, 'final_tissues.nii.gz');
            fid = fopen(f, 'w'); fwrite(fid, uint8(1:10)); fclose(fid);
            p = tc.params();
            p.path.seg = root;
            kinds = {'grid', 'medium', 'acoustic', 'headmodel'};
            old = struct();
            for k = kinds, old.(k{1}) = cache_signature(p, k{1}); end
            tc.verifyTrue(isfield(old.grid, 'segmentation'));
            % the same content written again (e.g. a copy) still matches
            fid = fopen(f, 'w'); fwrite(fid, uint8(1:10)); fclose(fid);
            java.io.File(f).setLastModified(java.lang.System.currentTimeMillis() + 5000);
            for k = kinds
                tc.verifyTrue(check_cache_signature('x.mat', old.(k{1}), cache_signature(p, k{1}), p));
            end
            fid = fopen(f, 'w'); fwrite(fid, uint8(11:20)); fclose(fid);
            java.io.File(f).setLastModified(java.lang.System.currentTimeMillis() + 10000);
            for k = kinds
                tc.verifyError(@() check_cache_signature('x.mat', old.(k{1}), ...
                    cache_signature(p, k{1}), p), 'prestus:cache:signatureMismatch');
            end
        end

        function test_legacy_cache_compared_by_segmentation_date(tc)
            % Caches from before the fingerprint: a segmentation written after
            % the cache counts as a change; the head model's recorded date is compared
            root = tempname; mkdir(root);
            tc.addTeardown(@() rmdir(root, 's'));
            m2m = fullfile(root, 'm2m_sub-001'); mkdir(m2m);
            f = fullfile(m2m, 'final_tissues.nii.gz');
            fid = fopen(f, 'w'); fwrite(fid, uint8(1:10)); fclose(fid);
            cache = fullfile(root, 'cache.mat');
            fid = fopen(cache, 'w'); fclose(fid);
            p = tc.params();
            p.path.seg = root;
            now_ms = java.lang.System.currentTimeMillis();

            legacy = rmfield(cache_signature(p, 'grid'), 'segmentation');
            java.io.File(f).setLastModified(now_ms - 60000);        % segmentation older than the cache
            java.io.File(cache).setLastModified(now_ms);
            tc.verifyTrue(check_cache_signature(cache, legacy, cache_signature(p, 'grid'), p));
            java.io.File(f).setLastModified(now_ms + 60000);        % charm ran after the cache
            [ok, diffs] = check_cache_signature(cache, legacy, cache_signature(p, 'grid'), ...
                setfield(p, 'io', struct('cache_signature_check', 'warn')));  %#ok<SFLD>
            tc.verifyFalse(ok);
            tc.verifyEqual(diffs, {'segmentation'});

            hm = rmfield(cache_signature(p, 'headmodel'), 'segmentation');
            d = dir(f);
            hm.segmentation_date = d.datenum;
            tc.verifyTrue(check_cache_signature(cache, hm, cache_signature(p, 'headmodel'), p));
            hm.segmentation_date = d.datenum - 1;
            tc.verifyError(@() check_cache_signature(cache, hm, cache_signature(p, 'headmodel'), p), ...
                'prestus:cache:signatureMismatch');
        end

        function test_unknown_kind_errors(tc)
            tc.verifyError(@() cache_signature(tc.params(), 'thermal'), 'cache_signature:badKind');
        end

    end

end
