classdef test_cache_signature < matlab.unittest.TestCase
% TEST_CACHE_SIGNATURE  Reused grid/medium caches must match the run's settings.
%
% Covers CACHE_SIGNATURE and CHECK_CACHE_SIGNATURE: which settings are part of
% a stage's signature, the error/warn/off modes of io.cache_signature_check,
% and legacy caches without a signature. Pure logic: no k-Wave, no data.
%
%   results = runtests('tests/test_cache_signature.m');

    methods (TestClassSetup)
        function add_paths(~)
            repo = fileparts(fileparts(mfilename('fullpath')));
            addpath(genpath(fullfile(repo, 'functions')));
        end
    end

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

        function test_unknown_kind_errors(tc)
            tc.verifyError(@() cache_signature(tc.params(), 'acoustic'), 'cache_signature:badKind');
        end

    end

end
