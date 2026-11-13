classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_library_conditions < matlab.unittest.TestCase
% TEST_LIBRARY_CONDITIONS  Unit tests for experimental dephased library conditions (issue #124).
%
% Synthetic only: no dephased hydrophone measurements exist yet, so the
% "measured" profile is generated from the model itself.
%
%   results = runtests('tests/test_library_conditions.m');

    properties
        equip
        tran
        lib
        tmp
    end

    methods (TestClassSetup)
        function setup_paths(tc)
            repo_root = fileparts(fileparts(mfilename('fullpath')));
            tc.equip = load_equipment_config();
            tc.tran  = tc.equip.trans.CTX_250_001;
            tc.lib   = yaml.loadFile(fullfile(repo_root, 'config', 'calibration', ...
                'CTX_250_001_SC_203_035.yaml'), 'ConvertToArray', true);
        end
    end

    methods (TestMethodSetup)
        function make_tmp(tc)
            tc.tmp = tempname; mkdir(tc.tmp);
            tc.addTeardown(@() rmdir(tc.tmp, 's'));
        end
    end

    methods (Test)

        function test_shift_phases_add_to_calibrated_setpoint(tc)
            sp = 50.02;
            cond = struct('label', 'geom_180', 'phase_shift_deg', [0 180 0 180], 'setpoint_ep_mm', sp);
            lib = library_add_conditions(tc.lib, tc.tran, cond);
            m = library_model_at_depth(tc.lib, tc.tran, sp, NaN);
            c = lib.conditions.geom_180;
            tc.verifyEqual(c.status, 'experimental');
            tc.verifyFalse(c.validated);
            tc.verifyEqual(mod(c.elem_phase_deg - m.elem_phase_deg - [0 180 0 180] + 180, 360) - 180, ...
                zeros(1, 4), 'AbsTol', 1e-3);
        end

        function test_absolute_without_correction(tc)
            cond = struct('label', 'abs_cond', 'elem_phase_deg', [0 120 240 0], 'setpoint_ep_mm', 50.02);
            lib = library_add_conditions(tc.lib, tc.tran, cond, struct('InheritCorrection', false));
            c = lib.conditions.abs_cond;
            tc.verifyEqual(c.elem_phase_deg, [0 120 240 0], 'AbsTol', 1e-9);
            tc.verifyEqual(c.elem_amp_weight, ones(1, 4));
        end

        function test_drive_reaches_reference_isppa(tc)
            cond = struct('label', 'geom_180', 'phase_shift_deg', [0 180 0 180], 'setpoint_ep_mm', 50.02);
            lib = library_add_conditions(tc.lib, tc.tran, cond, struct('IsppaWcm2', 20));
            c = lib.conditions.geom_180;
            ann = tc.tran.transducer.annular;
            z = (10:0.25:1.5 * ann.curv_radius_mm)';
            H = axial_transfer_matrix(tc.tran, z + ann.curv_radius_mm - ann.dist_geom_ep_mm, 1500, 994);
            w = c.elem_p_amp / (994 * 1500) * c.elem_amp_weight(:) .* exp(1i * c.elem_phase_deg(:) * pi/180);
            I = abs(H * w).^2 / (2 * 994 * 1500) * 1e-4;
            tc.verifyEqual(max(I), 20, 'RelTol', 2e-3);
        end

        function test_ingest_scores_measured_profile(tc)
            % "Measured" = model prediction -> R² ~ 1, focal error ~ 0
            [csv, ~] = tc.write_dephased_profile();
            lib_file = fullfile(tc.tmp, 'entry.yaml');
            yaml.dumpFile(lib_file, tc.lib);
            [lib, report] = tc.verifyWarning(@() library_ingest_conditions(lib_file, csv, tc.equip), ...
                'library_ingest_conditions:experimental');
            tc.verifyGreaterThan(report.r2, 0.999);
            tc.verifyLessThan(abs(report.focal_err_mm), 0.3);
            reread = yaml.loadFile(lib_file, 'ConvertToArray', true);
            tc.verifyTrue(isfield(reread.conditions, 'geom_180'));
            tc.verifyEqual(lib.conditions.geom_180.setpoint_ep_mm, 50.02, 'AbsTol', 1e-9);
        end

        function test_loader_selects_condition_and_warns(tc)
            cond = struct('label', 'geom_180', 'phase_shift_deg', [0 180 0 180], 'setpoint_ep_mm', 50.02);
            lib = library_add_conditions(tc.lib, tc.tran, cond, struct('IsppaWcm2', 30));
            yaml.dumpFile(fullfile(tc.tmp, 'CTX_250_001_SC_203_035.yaml'), lib);
            p = tc.verifyWarning(@() load_transducer_from_library('CTX_250_001_SC_203_035', [], 15, ...
                tc.equip, tc.tmp, 'geom_180'), 'library_condition:experimental');
            tc.verifyEqual(p.transducer.focal_distance_ep, 50.02, 'AbsTol', 1e-9);
            tc.verifyEqual(p.transducer.tx_calibration.phase_model, 'condition');
            tc.verifyEqual(p.transducer.tx_calibration.condition, 'geom_180');
            tc.verifyEqual(p.transducer.annular.elem_phase_deg, mod(lib.conditions.geom_180.elem_phase_deg, 360), 'AbsTol', 1e-9);
            expected = round(round(lib.conditions.geom_180.elem_p_amp * sqrt(15/30)) * lib.conditions.geom_180.elem_amp_weight);
            tc.verifyEqual(p.transducer.annular.elem_p_amp, expected);
        end

        function test_validated_condition_does_not_warn(tc)
            cond = struct('label', 'geom_180', 'phase_shift_deg', [0 180 0 180], 'setpoint_ep_mm', 50.02);
            lib = library_add_conditions(tc.lib, tc.tran, cond);
            lib.conditions.geom_180.validated = true;
            tc.verifyWarningFree(@() library_condition(lib, tc.tran, 'geom_180', 30));
        end

        function test_resolver_uses_condition(tc)
            cond = struct('label', 'geom_180', 'phase_shift_deg', [0 180 0 180], 'setpoint_ep_mm', 50.02);
            lib = library_add_conditions(tc.lib, tc.tran, cond);
            yaml.dumpFile(fullfile(tc.tmp, 'CTX_250_001_SC_203_035.yaml'), lib);
            tr = struct('serial', 'CTX_250_001', 'combo', struct('ds_serial', 'SC_203_035'), ...
                        'target_isppa_wcm2', 30, 'condition', 'geom_180', 'library_path', tc.tmp);
            tr = tc.verifyWarning(@() resolve_transducer_from_serial(tr, tc.equip, tc.tmp, 1), ...
                'library_condition:experimental');
            tc.verifyEqual(tr.focal_distance_ep, 50.02, 'AbsTol', 1e-9);
            tc.verifyEqual(tr.annular.elem_phase_deg, mod(lib.conditions.geom_180.elem_phase_deg, 360), 'AbsTol', 1e-9);
        end

        function test_errors(tc)
            tc.verifyError(@() library_condition(tc.lib, tc.tran, 'nope', 30), ...
                'library_condition:unknownCondition');
            bad = struct('label', 'x', 'phase_shift_deg', [0 180], 'setpoint_ep_mm', 50);
            tc.verifyError(@() library_add_conditions(tc.lib, tc.tran, bad), 'library_add_conditions:phaseLength');
            both = struct('label', 'x', 'phase_shift_deg', [0 0 0 0], 'elem_phase_deg', [0 0 0 0], 'setpoint_ep_mm', 50);
            tc.verifyError(@() library_add_conditions(tc.lib, tc.tran, both), 'library_add_conditions:phaseSpec');
            nosp = struct('label', 'x', 'phase_shift_deg', [0 180 0 180]);
            tc.verifyError(@() library_add_conditions(tc.lib, tc.tran, nosp), 'library_add_conditions:noSetpoint');
        end
    end

    methods (Access = private)
        function [csv, I] = write_dephased_profile(tc)
            cond = struct('label', 'geom_180', 'phase_shift_deg', [0 180 0 180], 'setpoint_ep_mm', 50.02);
            lib = library_add_conditions(tc.lib, tc.tran, cond);
            c = lib.conditions.geom_180;
            ann = tc.tran.transducer.annular;
            z = (2:0.5:90)';
            H = axial_transfer_matrix(tc.tran, z + ann.curv_radius_mm - ann.dist_geom_ep_mm, 1500, 994);
            I = abs(H * (c.elem_amp_weight(:) .* exp(1i * c.elem_phase_deg(:) * pi/180))).^2;
            I = 25 * I / max(I(z >= 10));
            csv = fullfile(tc.tmp, 'tran-CTX_250_001_ds-SC_203_035_acq-dephased_axial.csv');
            fid = fopen(csv, 'w');
            fprintf(fid, ',Focal depths (mm)\nDistance WRT Exit plane (mm),50.02\n');
            fprintf(fid, '%.2f,%.6f\n', [z, I]');
            fclose(fid);
            fid = fopen(strrep(csv, '.csv', '.yaml'), 'w');
            fprintf(fid, 'quantity: intensity\nunits: W/cm2\ncolumns:\n  - label: geom_180\n    phase_shift_deg: [0, 180, 0, 180]\n');
            fclose(fid);
        end
    end
end
