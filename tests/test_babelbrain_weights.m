classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_babelbrain_weights < matlab.unittest.TestCase
% TEST_BABELBRAIN_WEIGHTS  EXPORT_BABELBRAIN_WEIGHTS / IMPORT_BABELBRAIN_WEIGHTS
% against BabelBrain's OptimizedWeightsFile layout.
%
% BabelBrain writes and reads /CALIBRATION with BabelViscoFDTD's SaveToH5py /
% ReadFromH5py (TxCalibration.py, BabelIntegrationBASE.py): a 1-D h5py complex
% dataset (compound {r, i}) at the root, with a string attribute 'type' on every
% root item. The fixtures bb_calibration_ring.h5 (one weight per ring) and
% bb_calibration_subelem.h5 (one weight per sub-element, as TxCalibration fits)
% were written by fixtures/make_babelbrain_weights_fixture.py with the real
% SaveToH5py. No external data, Python or GPU needed.
%
%   Run with:   results = runtests('tests/test_babelbrain_weights.m');

    properties (Constant)
        RingPhaseDeg = [0 30 -45 120]   % phases encoded in the fixtures
        RingAmp      = [1 0.9 1.1 1]
    end

    properties
        Tran
        Params
        Dir
        Yaml
        RefDepth
    end

    methods (TestClassSetup)
        function setup(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(here, 'fixtures')));
            equip = load_equipment_config();
            serials = fieldnames(equip.trans);
            k = find(cellfun(@(s) isfield(equip.trans.(s).transducer, 'annular') ...
                && equip.trans.(s).n_elem == 4, serials), 1);
            tc.assumeNotEmpty(k, 'No shipped 4-ring annular transducer to match the fixtures.');
            tc.Tran = equip.trans.(serials{k});
            tc.RefDepth = tc.Tran.transducer.annular.dist_geom_ep_mm;
            tc.Params = struct();
            tc.Params.medium_properties.water.sound_speed = 1500;
            tc.Params.simulation.medium = 'water';
            tc.Params.transducer.annular.elem_n = 4;
            tc.Params.transducer.freq_hz = tc.Tran.transducer.freq_hz;
        end
    end

    methods (TestMethodSetup)
        function make_dir(tc)
            tc.Dir = tempname; mkdir(tc.Dir);
            tc.addTeardown(@() rmdir(tc.Dir, 's'));
            tc.Yaml = fullfile(tc.Dir, 'equipment.yaml');
            fid = fopen(tc.Yaml, 'w'); fprintf(fid, 'name: test\n'); fclose(fid);
        end
    end

    methods (Test, TestTags = {'unit', 'babelbrain'})

        function test_export_writes_babelbrain_layout(tc)
            phases = [0.1 1.2 -2.0 3.0];
            geo    = [0.0 0.5 0.5 -1.0];
            amp    = [2 1.8 2.2 2];
            out = fullfile(tc.Dir, 'w.h5');
            export_babelbrain_weights(phases, 0.1, tc.Params, out, geo, amp);

            info = h5info(out, '/CALIBRATION');
            tc.verifyEqual(info.Datatype.Class, 'H5T_COMPOUND');
            tc.verifyEqual({info.Datatype.Type.Member.Name}, {'r', 'i'});
            tc.verifyEqual(info.Dataspace.Size, 4);          % 1-D, as numpy writes it
            raw = h5read(out, '/CALIBRATION');
            w = double(raw.r) + 1i * double(raw.i);
            tc.verifyEqual(w(:).', (amp / amp(1)) .* exp(1i * (phases - geo)), 'AbsTol', 1e-6);

            % ReadFromH5py reads attrs['type'] of every item it visits.
            root = h5info(out, '/');
            names = [{root.Datasets.Name}, {root.Groups.Name}];
            for n = 1:numel(names)
                tc.verifyTrue(any(strcmp({h5info(out, ['/' strip(names{n}, 'left', '/')]).Attributes.Name}, 'type')), ...
                    sprintf('%s has no ''type'' attribute', names{n}));
            end
            tc.verifyEqual(h5readatt(out, '/CALIBRATION', 'type'), 'ndarray');
            tc.verifyEqual(h5readatt(out, '/metadata', 'type'), 'dict');
            tc.verifyEqual(h5readatt(out, '/metadata/freq_hz', 'type'), 'scalar');
            tc.verifyEqual(h5read(out, '/DRIVE_phases_rad').', phases);
        end

        function test_round_trip_export_import(tc)
            geo_deg   = set_real_phases(tc.Tran, tc.RefDepth, tc.Params);
            delta_deg = [0 25 -60 170];
            phases    = (geo_deg + delta_deg) * pi / 180;
            out = fullfile(tc.Dir, 'w.h5');
            export_babelbrain_weights(phases, 0.1, tc.Params, out, geo_deg * pi / 180);
            corr = import_babelbrain_weights(out, tc.Yaml, tc.Tran, tc.RefDepth, tc.Params);
            tc.verify_phase_equal(corr, delta_deg);
            d = yaml.loadFile(tc.Yaml, 'ConvertToArray', true);
            tc.verify_phase_equal(d.elem_phase_correction.deg, delta_deg);
            tc.verifyEqual(d.elem_phase_correction.ref_depth_ep_mm, tc.RefDepth);
        end

        function test_import_babelbrain_per_ring_file(tc)
            corr = import_babelbrain_weights(tc.fixture('bb_calibration_ring.h5'), ...
                tc.Yaml, tc.Tran, tc.RefDepth, tc.Params);
            tc.verify_phase_equal(corr, tc.RingPhaseDeg);
        end

        function test_import_babelbrain_per_subelement_file(tc)
            f = tc.fixture('bb_calibration_subelem.h5');
            tc.verifyError(@() import_babelbrain_weights(f, tc.Yaml, tc.Tran, tc.RefDepth, tc.Params), ?MException);
            corr = import_babelbrain_weights(f, tc.Yaml, tc.Tran, tc.RefDepth, tc.Params, [16 16 16 16]);
            tc.verify_phase_equal(corr, tc.RingPhaseDeg);
        end

        function test_import_legacy_split_layout(tc)
            w = exp(1i * tc.RingPhaseDeg * pi / 180).';
            out = fullfile(tc.Dir, 'legacy.h5');
            h5create(out, '/CALIBRATION/real', [4 1], 'Datatype', 'single');
            h5create(out, '/CALIBRATION/imag', [4 1], 'Datatype', 'single');
            h5write(out, '/CALIBRATION/real', single(real(w)));
            h5write(out, '/CALIBRATION/imag', single(imag(w)));
            corr = import_babelbrain_weights(out, tc.Yaml, tc.Tran, tc.RefDepth, tc.Params);
            tc.verify_phase_equal(corr, tc.RingPhaseDeg);
        end

    end

    methods (Access = private)
        function p = fixture(~, name)
            p = fullfile(fileparts(mfilename('fullpath')), 'fixtures', name);
        end

        function verify_phase_equal(tc, actual_deg, expected_deg)
            % Compare phases modulo 360 deg.
            err = mod(actual_deg(:).' - expected_deg(:).' + 180, 360) - 180;
            tc.verifyLessThan(max(abs(err)), 1e-3);
        end
    end
end
