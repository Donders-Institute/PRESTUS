classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_resolve_geometric_entry < matlab.unittest.TestCase
% TEST_RESOLVE_GEOMETRIC_ENTRY  A transducer given by serial whose library entry
% is geometric (phases only, no drive) resolves: the phases come from the entry,
% the amplitude stays empty and FILL_ANNULAR_DRIVE sets the nominal drive. An
% entry with a calibrated drive keeps it. Uses the shipped equipment and
% calibration library; no k-Wave, no data.
%
%   results = runtests('tests/test_resolve_geometric_entry.m');

    properties
        Equip
        Lib
    end

    methods (TestClassSetup)
        function load_library(tc)
            tc.Equip = load_equipment_config();
            tc.Lib   = fullfile(get_prestus_path(), 'config', 'calibration');
        end
    end

    methods (Test)
        function test_geometric_entry_resolves_without_drive(tc)
            [~, found] = library_entry_file(tc.Lib, 'CTX_500');
            tc.assumeTrue(found, 'shipped geometric entry CTX_500 missing');
            n  = tc.Equip.trans.CTX_500.transducer.annular.elem_n;
            tr = struct('serial', 'CTX_500', 'focal_distance_ep', 60, 'target_isppa_wcm2', 30);
            out = resolve_transducer_from_serial(tr, tc.Equip, tc.Lib, 1);
            tc.verifyNumElements(out.annular.elem_phase_deg, n);
            tc.verifyTrue(~isfield(out.annular, 'elem_p_amp') || isempty(out.annular.elem_p_amp), ...
                'a geometric entry has no drive to copy');
            out = fill_annular_drive(out, 1500, 1);
            tc.verifyEqual(out.annular.elem_p_amp, repmat(100e3, 1, n));
        end

        function test_calibrated_entry_keeps_its_drive(tc)
            % Negative control: an entry with a drive still sets it
            tr = struct('serial', 'CTX_500_026', 'focal_distance_ep', 60, 'target_isppa_wcm2', 30);
            tr.combo.ds_serial = 'SC_203_035';
            out = resolve_transducer_from_serial(tr, tc.Equip, tc.Lib, 1);
            tc.verifyNotEmpty(out.annular.elem_p_amp);
            filled = fill_annular_drive(out, 1500, 1);
            tc.verifyEqual(filled.annular.elem_p_amp, out.annular.elem_p_amp);
        end
    end
end
