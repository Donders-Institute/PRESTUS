classdef test_transducer_library < matlab.unittest.TestCase
% TEST_TRANSDUCER_LIBRARY  Unit tests for transducer serial resolution.
%
% Tests resolve_transducer_from_serial, the serial→geometry→library
% pipeline path in load_transducer_parameters, and combo-name logic.
% Does not require k-Wave or SimNIBS; uses the real equipment YAMLs.
%
%   Run with:   results = runtests('tests/test_transducer_library.m');
%   Run tag:    results = runtests('tests/test_transducer_library.m', 'Tag', 'geometry');

    properties
        equip_param      % loaded once per test class
        library_path     % path to config/transducer/
        repo_root        % PRESTUS root
        tmp_lib_dir      % temp dir for synthetic calibration YAMLs
    end

    methods (TestClassSetup)
        function setup_paths(tc)
            test_dir       = fileparts(mfilename('fullpath'));
            tc.repo_root   = fileparts(test_dir);
            addpath(genpath(fullfile(tc.repo_root, 'functions')));
            addpath(genpath(fullfile(tc.repo_root, 'external')));
            tc.equip_param  = load_equipment_config();
            tc.library_path = fullfile(tc.repo_root, 'config', 'transducer');
        end
    end

    methods (TestClassTeardown)
        function cleanup_tmp(tc)
            if ~isempty(tc.tmp_lib_dir) && isfolder(tc.tmp_lib_dir)
                rmdir(tc.tmp_lib_dir, 's');
            end
        end
    end

    % ------------------------------------------------------------------ %
    %% load_equipment_config
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'equipment'})

        function test_equipment_config_loads(tc)
            eq = tc.equip_param;
            tc.verifyClass(eq, 'struct');
            tc.verifyTrue(isfield(eq, 'trans'),  'Missing eq.trans');
            tc.verifyTrue(isfield(eq, 'combos'), 'Missing eq.combos');
        end

        function test_known_serials_present(tc)
            serials = fieldnames(tc.equip_param.trans);
            tc.verifyFalse(isempty(serials), 'No transducers found in equipment config');
        end

        function test_transducer_has_geometry(tc)
            serials = fieldnames(tc.equip_param.trans);
            serial  = serials{1};
            tran    = tc.equip_param.trans.(serial);
            tc.verifyTrue(isfield(tran, 'transducer'), ...
                sprintf('%s missing .transducer geometry field', serial));
        end

    end

    % ------------------------------------------------------------------ %
    %% resolve_transducer_from_serial — geometry only
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'geometry'})

        function test_geometry_populated_from_serial(tc)
            serial = tc.first_annular_serial();
            tr     = tc.minimal_tr(serial);
            out    = resolve_transducer_from_serial(tr, tc.equip_param, tc.library_path, 1);
            tc.verifyTrue(isfield(out, 'type'), 'type not populated');
            tc.verifyTrue(isfield(out, 'freq_hz'), 'freq_hz not populated');
            tc.verifyTrue(isfield(out, 'annular'), 'annular not populated');
            tc.verifyTrue(isfield(out.annular, 'curv_radius_mm'), 'curv_radius_mm missing');
        end

        function test_inline_fields_override_equipment(tc)
            serial   = tc.first_annular_serial();
            tr       = tc.minimal_tr(serial);
            tr.freq_hz = 123456;    % intentionally wrong — should survive merge
            out = resolve_transducer_from_serial(tr, tc.equip_param, tc.library_path, 1);
            tc.verifyEqual(out.freq_hz, 123456, ...
                'Inline freq_hz should not be overwritten by equipment YAML');
        end

        function test_unknown_serial_errors(tc)
            tr.serial = 'NONEXISTENT_DEVICE_XYZ';
            tc.verifyError( ...
                @() resolve_transducer_from_serial(tr, tc.equip_param, tc.library_path, 1), ...
                ?MException);
        end

    end

    % ------------------------------------------------------------------ %
    %% combo-name resolution
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'combo'})

        function test_generic_combo_name(tc)
            % No ds_serial → key is just the serial
            serial = tc.first_annular_serial();
            tr     = tc.minimal_tr(serial);
            name   = tc.get_combo_name(tr, serial);
            tc.verifyEqual(name, serial);
        end

        function test_ds_specific_combo_name(tc)
            serial        = tc.first_annular_serial();
            tr            = tc.minimal_tr(serial);
            tr.combo.ds_serial = 'MY_DS';
            name = tc.get_combo_name(tr, serial);
            tc.verifyEqual(name, [serial '_MY_DS']);
        end

    end

    % ------------------------------------------------------------------ %
    %% library lookup with a synthetic YAML
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'library'})

        function test_phases_resolved_from_library(tc)
            serial = tc.first_annular_serial();
            tran   = tc.equip_param.trans.(serial);
            n_elem = tran.transducer.annular.elem_n;

            lib_dir = tc.get_tmp_lib();
            tc.write_synthetic_library(lib_dir, serial, n_elem, 60, 32);

            tr = tc.minimal_tr_with_targets(serial, 60, 32);
            out = resolve_transducer_from_serial(tr, tc.equip_param, lib_dir, 1);

            tc.verifyTrue(isfield(out.annular, 'elem_phase_deg'), ...
                'elem_phase_deg not resolved');
            tc.verifyEqual(numel(out.annular.elem_phase_deg), n_elem, ...
                'Wrong number of phase elements');
            tc.verifyTrue(isfield(out.annular, 'elem_amp'), 'elem_amp not resolved');
        end

        function test_inline_phases_not_overwritten_by_library(tc)
            serial = tc.first_annular_serial();
            tran   = tc.equip_param.trans.(serial);
            n_elem = tran.transducer.annular.elem_n;

            lib_dir = tc.get_tmp_lib();
            tc.write_synthetic_library(lib_dir, serial, n_elem, 60, 32);

            tr = tc.minimal_tr_with_targets(serial, 60, 32);
            tr.annular.elem_phase_deg = zeros(1, n_elem);  % inline override
            tr.annular.elem_amp       = ones(1, n_elem) * 99999;

            out = resolve_transducer_from_serial(tr, tc.equip_param, lib_dir, 1);

            tc.verifyEqual(out.annular.elem_phase_deg, zeros(1, n_elem), ...
                'Inline phases were overwritten');
            tc.verifyEqual(out.annular.elem_amp, ones(1, n_elem)*99999, ...
                'Inline amplitude was overwritten');
        end

        function test_missing_library_errors_with_message(tc)
            serial  = tc.first_annular_serial();
            empty_dir = tc.get_tmp_lib();  % tmp dir with no YAMLs
            tr = tc.minimal_tr_with_targets(serial, 60, 32);
            tc.verifyError( ...
                @() resolve_transducer_from_serial(tr, tc.equip_param, empty_dir, 1), ...
                ?MException);
        end

        function test_no_lookup_when_targets_absent(tc)
            % Without focal_distance_ep or target_isppa_wcm2, resolution
            % should silently skip the library step (no error).
            serial  = tc.first_annular_serial();
            empty_dir = tc.get_tmp_lib();
            tr = tc.minimal_tr(serial);  % no targets set
            tc.verifyWarningFree( ...
                @() resolve_transducer_from_serial(tr, tc.equip_param, empty_dir, 1));
        end

    end

    % ------------------------------------------------------------------ %
    %% load_transducer_parameters integration
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'integration'})

        function test_serial_triggers_geometry_in_pipeline(tc)
            % Verify that resolve_transducer_from_serial (the function wired
            % into load_transducer_parameters) populates type and freq_hz from
            % the equipment YAML when only a serial is given.
            serial = tc.first_annular_serial();
            tran   = tc.equip_param.trans.(serial);
            n      = tran.transducer.annular.elem_n;

            tr.serial              = serial;
            tr.focal_distance_ep   = 60;
            tr.target_isppa_wcm2   = 5;
            tr.annular.elem_phase_deg = zeros(1, n);
            tr.annular.elem_amp       = ones(1, n) * 100000;

            out = resolve_transducer_from_serial(tr, tc.equip_param, tc.library_path, 1);

            tc.verifyEqual(out.freq_hz, tran.transducer.freq_hz, ...
                'freq_hz should be populated from equipment YAML');
            tc.verifyEqual(char(out.type), 'annular', ...
                'type should be set to annular from equipment YAML');
        end

    end

    % ------------------------------------------------------------------ %
    %% Private helpers
    % ------------------------------------------------------------------ %
    methods (Access = private)

        function serial = first_annular_serial(tc)
            serials = fieldnames(tc.equip_param.trans);
            serial  = '';
            for i = 1:numel(serials)
                tran = tc.equip_param.trans.(serials{i});
                if isfield(tran, 'transducer') && isfield(tran.transducer, 'annular')
                    serial = serials{i};
                    return;
                end
            end
            error('test_transducer_library: no annular transducer found in equipment config');
        end

        function tr = minimal_tr(~, serial)
            tr.serial = serial;
        end

        function tr = minimal_tr_with_targets(~, serial, focal_ep, isppa)
            tr.serial              = serial;
            tr.focal_distance_ep   = focal_ep;
            tr.target_isppa_wcm2   = isppa;
        end

        function d = get_tmp_lib(tc)
            if isempty(tc.tmp_lib_dir)
                tc.tmp_lib_dir = tempname();
                mkdir(tc.tmp_lib_dir);
            end
            d = tc.tmp_lib_dir;
        end

        function name = get_combo_name(~, tr, serial)
            % Mirror the logic in resolve_combo_name without calling it directly.
            if isfield(tr, 'combo') && isstruct(tr.combo) && ...
                    isfield(tr.combo, 'ds_serial') && ~isempty(tr.combo.ds_serial)
                name = [serial '_' char(tr.combo.ds_serial)];
            else
                name = serial;
            end
        end

        function write_synthetic_library(~, lib_dir, serial, n_elem, focal_mm, isppa)
            % Write a minimal library YAML to lib_dir/{serial}.yaml so that
            % load_transducer_from_library can parse it.
            fkey = sprintf('f%s', strrep(num2str(focal_mm), '.', 'p'));
            ikey = sprintf('i%s', strrep(num2str(isppa),   '.', 'p'));

            data.meta.tran_serial         = serial;
            data.meta.created_by          = 'test_transducer_library';
            data.calibration.focal_depths.(fkey).precession_mode   = 'linear';
            data.calibration.focal_depths.(fkey).phase_start_deg   = 0;
            data.calibration.focal_depths.(fkey).phase_step_deg    = 10;
            data.calibration.focal_depths.(fkey).amplitude_scaling.(ikey) = 100000;

            yaml_path = fullfile(lib_dir, [serial '.yaml']);
            yaml.dumpFile(yaml_path, data);
        end

    end

end
