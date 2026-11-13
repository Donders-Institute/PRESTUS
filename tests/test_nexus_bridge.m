classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_nexus_bridge < matlab.unittest.TestCase
% TEST_NEXUS_BRIDGE  The MATLAB side of the PRESTUS Nexus: NEXUS_BRIDGE (file protocol,
% allow-list, errors), NEXUS_NEURONAV on a synthetic Localite study, NEXUS_CALIBRATION on
% synthetic profiles of a shipped transducer, and NEURONAV_CHECK_TARGET_MAP.

    properties
        Dir
    end

    methods (TestMethodSetup)
        function make_dir(tc)
            tc.Dir = tempname; mkdir(tc.Dir);
            tc.addTeardown(@() rmdir(tc.Dir, 's'));
        end
    end

    methods (Test, TestTags = {'unit', 'nexus'})

        function test_bridge_serves_requests_and_refuses_others(tc)
            b = fullfile(tc.Dir, 'bridge');
            mkdir(fullfile(b, 'req'));
            req = @(id, fn, args) tc.write(fullfile(b, 'req', [id '.json']), ...
                jsonencode(struct('id', id, 'fn', fn, 'args', args)));
            req('001', 'nexus_ping', struct('x', [1 2]));
            req('002', 'system', struct('cmd', 'echo hi'));
            req('003', 'nexus_calibration', struct('action', 'view', 'session', 'none'));
            nexus_bridge(b, IdleTimeout=0.5, Poll=0.01);               % serves, then ends when idle
            tc.verifyTrue(isfile(fullfile(b, 'ready.json')));
            r1 = jsondecode(fileread(fullfile(b, 'resp', '001.json')));
            tc.verifyTrue(r1.ok);
            tc.verifyEqual(r1.result.echo.x, [1; 2]);
            r2 = jsondecode(fileread(fullfile(b, 'resp', '002.json')));
            tc.verifyFalse(r2.ok);
            tc.verifyEqual(r2.identifier, 'nexus_bridge:notAllowed');
            r3 = jsondecode(fileread(fullfile(b, 'resp', '003.json')));
            tc.verifyEqual(r3.identifier, 'nexus_calibration:noSession');
            tc.verifyEmpty(dir(fullfile(b, 'req', '*.json')));
        end

        function test_bridge_ends_without_heartbeat(tc)
            b = fullfile(tc.Dir, 'bridge');
            mkdir(b);
            fclose(fopen(fullfile(b, 'heartbeat'), 'w'));
            t = tic;
            nexus_bridge(b, HeartbeatTimeout=1, Poll=0.05);
            tc.verifyLessThan(toc(t), 30);
            tc.verifyFalse(isfile(fullfile(b, 'heartbeat')));
        end

        function test_check_target_map(tc)
            e = @(n, c) struct('series_index', 1, 'transducer_id', 0, 'target_name', n, 'color', c);
            tc.verifyEqual(neuronav_check_target_map([e('A', ''), e('B', '#00ff00')]), '');
            tc.verifySubstring(neuronav_check_target_map([]), 'Tick at least one');
            tc.verifySubstring(neuronav_check_target_map(e(' ', '')), 'needs a target name');
            tc.verifySubstring(neuronav_check_target_map([e('L-PUL', ''), e('L_PUL', '')]), 'unique');
            tc.verifySubstring(neuronav_check_target_map(e('A', 'red')), '#rrggbb');
        end

        function test_neuronav_detect_write_scan(tc)
            rz = @(d) [cosd(d) -sind(d) 0; sind(d) cosd(d) 0; 0 0 1];
            make_synthetic_localite_export(fullfile(tc.Dir, 'raw'), 'sub-001', 'ses-01', ...
                {[rz(20), [40; 90; 140]; 0 0 0 1], [rz(160), [140; 95; 135]; 0 0 0 1]}, [1; 1; 0], ...
                {'LEFT PUL', 'RIGHT PUL'});
            cfg = fullfile(tc.Dir, 'config_study.yaml');
            tc.write(cfg, sprintf("path:\n  localite_raw: '%s'\n  localite_post: '%s'\n", ...
                fullfile(tc.Dir, 'raw'), fullfile(tc.Dir, 'post')));
            a = struct('config', cfg, 'sub', 'sub-001', 'ses', 'ses-01');
            s = nexus_neuronav(setfield(a, 'action', 'sessions')); %#ok<SFLD>
            tc.verifyEqual(s.subjects{1}.sessions, {'ses-01'});
            d = nexus_neuronav(setfield(a, 'action', 'detect')); %#ok<SFLD>
            tc.verifyEqual(numel(d.series), 2);
            tc.verifyEqual(d.series{1}.name, 'LEFT PUL');
            tc.verifyEqual(d.series{1}.n_pulses, 80);
            tc.verifyEqual(numel(d.series{1}.positions), 80);
            tc.verifySubstring(d.series{1}.name_from, 'off axis');
            jsonencode(d, 'ConvertInfAndNaN', true);                   % must be sendable
            tm = jsondecode(['[{"series_index":1,"transducer_id":0,"target_name":"L_PUL","color":""},' ...
                '{"series_index":2,"transducer_id":0,"target_name":"L_PUL","color":""}]']);
            tc.verifyError(@() nexus_neuronav(setfield(setfield(a, 'action', 'write'), 'target_map', tm)), ...
                'nexus_neuronav:targetMap'); %#ok<SFLD>
            tm(2).target_name = 'R_PUL';
            w = nexus_neuronav(setfield(setfield(a, 'action', 'write'), 'target_map', tm)); %#ok<SFLD>
            tc.verifyTrue(isfile(w.xml) && isfile(w.json));
            tc.verifyEqual(w.n_targets, 2);
            sc = nexus_neuronav(struct('action', 'scan', 'config', cfg));
            tc.verifyTrue(sc.sessions{1}.ingested);
            tc.verifyEqual(numel(sc.labels), 2);
        end

        function test_calibration_open_view_fit_save(tc)
            equip = load_equipment_config();
            serials = fieldnames(equip.trans);
            serial = serials{find(cellfun(@(s) isfield(equip.trans.(s).transducer, 'annular'), serials), 1)};
            tran = equip.trans.(serial);
            d0 = tran.transducer.annular.dist_geom_ep_mm;
            z = (0:0.25:1.4 * d0)'; dm = round(d0 + [-10 0 10], 2); I = zeros(numel(z), 3);
            for j = 1:3
                m = fit_transducer_model(tran, struct('depths_ep_mm', dm(j), 'z_ep_mm', z, ...
                    'I', ones(numel(z), 1)), struct('EvaluateOnly', true, 'NearFieldSkipMm', 0));
                I(:, j) = 30 * interp1(m.z_ep_mm, m.I_model, z, 'linear', 0);
            end
            csv = fullfile(tc.Dir, 'prof.csv');
            write_axial_profiles(struct('depths_ep_mm', dm, 'z_ep_mm', z, 'I', I), csv);
            lib = fullfile(tc.Dir, 'lib');
            st = struct('n_starts', 2);
            o = nexus_calibration(struct('action', 'open', 'session', 't', 'serial', serial, ...
                'driving_system', '', 'library_path', lib, 'csv', csv));
            tc.verifyEqual(o.measured_depths, dm);
            tc.verifyNotEmpty(o.default_depths);
            v = nexus_calibration(struct('action', 'view', 'session', 't', 'predicted', false, 'settings', st));
            tc.verifyEqual(fieldnames(v.results), {'geometric'});
            tc.verifyEqual(numel(v.results.geometric.Imod), 3);           % one cell per column
            f = nexus_calibration(struct('action', 'fit', 'session', 't', 'settings', st, ...
                'mode', 'per_depth', 'free_params', {{'phase'; 'amp'}}, 'calibrate_at', 'measured', ...
                'measured_depths', dm(2)));
            tc.verifyEqual(f.params.kind, 'per_depth');
            tc.verifyEqual(f.fit_record.tool, 'nexus');
            jsonencode(f, 'ConvertInfAndNaN', true);
            s = nexus_calibration(struct('action', 'save', 'session', 't', 'out_dir', lib));
            tc.verifyTrue(isfile(s.file));
            w = nexus_calibration(struct('action', 'write_fit', 'session', 't', 'path', fullfile(tc.Dir, 'fit.mat')));
            tc.verifyEqual(w.depths_ep_mm, dm(2));
            x = load(w.path);
            tc.verifyEqual(x.fit.mode, 'per_depth');
            o = nexus_calibration(struct('action', 'open', 'session', 't', 'serial', serial, ...
                'driving_system', '', 'library_path', lib, 'csv', csv));
            tc.verifyEqual(o.library_file, s.file);
            tc.verifyEqual(char(o.fit_record.tool), 'nexus');
            nexus_calibration(struct('action', 'close', 'session', 't'));
            tc.verifyError(@() nexus_calibration(struct('action', 'view', 'session', 't')), ...
                'nexus_calibration:noSession');
        end
    end

    methods
        function write(~, f, txt)
            fid = fopen(f, 'w'); fwrite(fid, txt, 'char'); fclose(fid);
        end
    end
end
