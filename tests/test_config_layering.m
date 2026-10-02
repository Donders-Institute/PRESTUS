classdef test_config_layering < matlab.unittest.TestCase
% TEST_CONFIG_LAYERING  Toolbox → project → study config layering and key checks.
%
% Covers LOAD_PARAMETERS' three-layer merge (toolbox default, project copy of
% config_default.yaml, study config), the stale-project-default and
% unknown-key warnings (CONFIG_UNKNOWN_KEYS, CONFIG_KNOWN_EXTRA_KEYS), the
% current-folder default warning, non-interactive handling of missing paths,
% PRESTUS_CONFIG_UPGRADE and GET_PRESTUS_PATH. Writes only to temporary folders.
%
%   results = runtests('tests/test_config_layering.m');

    properties
        Repo
        Dir
    end

    methods (TestClassSetup)
        function add_paths(tc)
            tc.Repo = fileparts(fileparts(mfilename('fullpath')));
            addpath(genpath(fullfile(tc.Repo, 'functions')));
            addpath(genpath(fullfile(tc.Repo, 'external')));
        end
    end

    methods (TestMethodSetup)
        function make_project(tc)
            tc.Dir = tempname;
            mkdir(tc.Dir);
            tc.addTeardown(@() rmdir(tc.Dir, 's'));
            % A deliberately old, partial project copy of config_default.yaml
            write_text(fullfile(tc.Dir, 'config_default.yaml'), { ...
                'simulation:', '  medium: water', '  debug: 0', ...
                'grid:', '  resolution_mm: 0.75'});
            write_text(fullfile(tc.Dir, 'config_study.yaml'), { ...
                'subject_id: 7'});
        end
    end

    methods (Test, TestTags = {'config_layering'})

        function test_project_values_override_toolbox(tc)
            p = tc.load_quiet('config_study.yaml');
            tc.verifyEqual(p.grid.resolution_mm, 0.75);
            tc.verifyEqual(p.subject_id, 7);
        end

        function test_keys_missing_from_project_copy_get_toolbox_defaults(tc)
            % The project copy predates io.cache_signature_check and most other keys
            p = tc.load_quiet('config_study.yaml');
            tc.verifyTrue(isfield(p.io, 'cache_signature_check'));
            tc.verifyTrue(isfield(p.modules, 'run_water_baseline'));
        end

        function test_stale_project_copy_warns(tc)
            tc.verifyWarning(@() load_parameters('config_study.yaml', tc.Dir), ...
                'prestus:staleProjectDefault');
        end

        function test_unknown_study_key_warns_with_suggestion(tc)
            write_text(fullfile(tc.Dir, 'config_typo.yaml'), {'grid:', '  resolutoin_mm: 1'});
            % The suggestion text itself is checked in test_unknown_keys_helper
            tc.verifyWarning(@() load_parameters('config_typo.yaml', tc.Dir), ...
                'prestus:unknownConfigKey');
        end

        function test_unknown_keys_helper(tc)
            ref  = struct('grid', struct('resolution_mm', 1, 'pml_size', 10), ...
                          'transducer', struct('freq_hz', 1, 'annular', struct('elem_n', 1)));
            user = struct('grid', struct('resolutoin_mm', 1), ...
                          'transducer', struct('freq_hz', {1, 2}, 'serial', {'a', 'b'}), ...
                          'acoustic_backend', 'kwave');
            [u, hints] = config_unknown_keys(user, ref, {'acoustic_backend', 'transducer.serial'});
            tc.verifyEqual(u, {'grid.resolutoin_mm'});
            tc.verifyEqual(hints, {'resolution_mm'});
        end

        function test_free_form_values_not_inspected(tc)
            % A struct where the default holds a scalar/empty value is a map
            ref  = struct('neuronav', struct('coil_map', []));
            user = struct('neuronav', struct('coil_map', struct('Coil0', 'T1')));
            tc.verifyEmpty(config_unknown_keys(user, ref));
        end

        function test_extra_keys_read_from_overlay(tc)
            keys = config_known_extra_keys(tc.Repo);
            tc.verifyTrue(all(ismember({'acoustic_backend', 'transducer.serial', ...
                'calibration.target_isppa_wcm2'}, keys)));
        end

        function test_shipped_example_configs_use_known_keys(tc)
            % config_tutorial.yaml must not contain keys PRESTUS ignores
            known = yaml.loadFile(fullfile(tc.Repo, 'config', 'config_default.yaml'), ...
                'ConvertToArray', true);
            extra = config_known_extra_keys(tc.Repo);
            cfg   = yaml.loadFile(fullfile(tc.Repo, 'config', 'config_tutorial.yaml'), ...
                'ConvertToArray', true);
            tc.verifyEmpty(config_unknown_keys(cfg, known, extra));
        end

        function test_cwd_default_warns(tc)
            tc.applyFixture(matlab.unittest.fixtures.CurrentFolderFixture(tc.Dir));
            tc.verifyWarning(@() load_parameters(struct('simulation', struct('debug', 0))), ...
                'prestus:cwdDefault');
        end

        function test_missing_charm_warns_instead_of_dialog(tc)
            % Non-interactive: no questdlg, just a warning
            s = struct('simulation', struct('debug', 0, 'interactive', 0), ...
                       'startup', struct('simnibs_bin_path', fullfile(tc.Dir, 'no_simnibs')));
            tc.verifyWarning(@() load_parameters(s), 'prestus:missingCharm');
        end

        function test_upgrade_dry_run_reports_without_changes(tc)
            before = fileread(fullfile(tc.Dir, 'config_default.yaml'));
            r = evalc_report(@() prestus_config_upgrade(tc.Dir, DryRun=true));
            tc.verifyNotEmpty(r.added);
            tc.verifyTrue(any(strcmp(r.changed.key, 'grid.resolution_mm')));
            tc.verifyEqual(fileread(fullfile(tc.Dir, 'config_default.yaml')), before);
            tc.verifyEqual(r.backup, '');
        end

        function test_upgrade_replaces_copy_and_keeps_backup(tc)
            before = fileread(fullfile(tc.Dir, 'config_default.yaml'));
            r = evalc_report(@() prestus_config_upgrade(tc.Dir));
            tc.verifyTrue(isfile(r.backup));
            tc.verifyEqual(fileread(r.backup), before);
            tc.verifyEqual(fileread(fullfile(tc.Dir, 'config_default.yaml')), ...
                fileread(fullfile(tc.Repo, 'config', 'config_default.yaml')));
            % Now up to date: loading no longer warns about a stale copy
            % (promote that warning to an error so it would fail the call)
            prev = warning('error', 'prestus:staleProjectDefault');
            tc.addTeardown(@() warning(prev));
            p = tc.load_quiet_except('config_study.yaml', 'prestus:staleProjectDefault');
            tc.verifyEqual(p.subject_id, 7);
        end

        function test_get_prestus_path_finds_root(tc)
            root = get_prestus_path();
            tc.verifyTrue(isfile(fullfile(root, 'config', 'config_default.yaml')));
            tc.verifyEqual(root, tc.Repo);
        end

    end

    methods (Access = private)
        function p = load_quiet(tc, file)
            s = warning('off', 'all');
            restore = onCleanup(@() warning(s));
            p = load_parameters(file, tc.Dir);
        end

        function p = load_quiet_except(tc, file, id)
            % Silence all warnings except id (which keeps its current state)
            state = warning('query', id);
            s = warning('off', 'all');
            warning(state.state, id);
            restore = onCleanup(@() warning(s));
            p = load_parameters(file, tc.Dir);
        end
    end

end

function write_text(path, lines)
    fid = fopen(path, 'w');
    fprintf(fid, '%s\n', lines{:});
    fclose(fid);
end

function r = evalc_report(fn)
    r = [];
    evalc('r = fn();');
end
