classdef test_neuronav_io < matlab.unittest.TestCase
% TEST_NEURONAV_IO  Localite I/O: folder layout, path keys and pose selection.
%
% Self-contained: builds synthetic Localite TriggerMarkers exports with known
% 4×4 poses in a temporary folder, so no study data is needed.
%
%   Pose regression: the matrix handed to localite_matrix_to_positions must
%   have the recorded position in column 4 (a transpose used to put it at
%   the origin).

    properties
        Root    % temporary folder
        Poses   % {M1, M2}: the two recorded 4×4 poses
    end

    methods (TestMethodSetup)
        function make_root(tc)
            tc.Root = tempname;
            mkdir(tc.Root);
            tc.addTeardown(@() rmdir(tc.Root, 's'));
            R1 = rotz(20);  R2 = rotz(-35);
            tc.Poses = {[R1, [ 40; 90; 140]; 0 0 0 1], ...
                        [R2, [-45; 85; 135]; 0 0 0 1]};
        end
    end

    methods (Test, TestTags = {'unit', 'localite'})

        function test_file_mode_uses_column4_position(tc)
            f = tc.write_export(fullfile(tc.Root, 'raw'), 'sub-001', 'ses-01', 'loc');
            p = base_parameters();
            p.placement.localite.file = f;
            for k = 1:2
                p.placement.localite.position = k;
                M = neuronav_localite_pose(p);
                tc.verifyEqual(M(1:3, 4), tc.Poses{k}(1:3, 4), 'AbsTol', 1e-6);
                tc.verifyEqual(M(1:3, 1), tc.Poses{k}(1:3, 1), 'AbsTol', 1e-6);
            end
        end

        function test_position_transducer_localite_not_at_origin(tc)
            f = tc.write_export(fullfile(tc.Root, 'raw'), 'sub-001', 'ses-01', 'loc');
            p = base_parameters();
            hdr = struct('Transform', struct('T', eye(4)), 'ImageSize', [256 256 256], ...
                'PixelDimensions', [1 1 1]);
            [~, ~, trans_ras] = position_transducer_localite(f, hdr, p);
            M = tc.Poses{1};
            expected = M(1:3, 4)' + p.placement.localite.tracker_to_bowl_mm * M(1:3, 1)';
            tc.verifyEqual(trans_ras, expected, 'AbsTol', 1e-6);
        end

        function test_automatic_mode_ingests_and_selects(tc)
            for loc = {'loc', ''}   % '' = no Localite sub-folder level
                raw = fullfile(tc.Root, ['raw_' char(loc{1})]);
                tc.write_export(raw, 'sub-001', 'ses-02', loc{1});
                p = base_parameters();
                p.path.localite_raw  = raw;
                p.path.localite_post = fullfile(tc.Root, ['post_' char(loc{1})]);
                p.placement.localite.loc_folder = loc{1};
                p.placement.localite.session  = 2;
                p.placement.localite.position = 2;

                [M, info] = neuronav_localite_pose(p);
                tc.verifyEqual(info.mode, 'automatic');
                tc.verifyEqual(info.series_index, 2);
                tc.verifyEqual(M(1:3, 4), tc.Poses{2}(1:3, 4), 'AbsTol', 1e-6);
                tc.verifyTrue(isfile(fullfile(p.path.localite_post, 'sub-001', 'ses-02', ...
                    'InstrumentMarker_sub-001_ses-02.json')));

                % By name, from the (now existing) ingested file
                p.placement.localite.position = 1;
                p.placement.localite.target = 'series_2';
                M = neuronav_localite_pose(p);
                tc.verifyEqual(M(1:3, 4), tc.Poses{2}(1:3, 4), 'AbsTol', 1e-6);

                % A raw series that was not ingested must not be guessed
                p.placement.localite.target = '';
                tc.verifyError(@() neuronav_localite_pose(p), 'PRESTUS:localite:noSeries');
            end
        end

        function test_one_pose_per_transducer(tc)
            % placement.localite.position/target may list one entry per
            % transducer (bilateral setups): transducer i uses entry i.
            f = tc.write_export(fullfile(tc.Root, 'raw'), 'sub-001', 'ses-01', 'loc');
            anat = fullfile(tc.Root, 'anat');
            mkdir(anat);
            niftiwrite(zeros(20, 20, 20, 'single'), fullfile(anat, 'sub-001_T1w.nii'));
            p = base_parameters();
            p.path.anat = anat;
            p.path.t1_pattern = 'sub-%03d_T1w.nii';
            p.placement.mode = 'localite';
            p.placement.report = false;
            p.placement.localite.file = f;
            p.placement.localite.position = [2 1];
            p.transducer(2) = p.transducer(1);
            out = preproc_transducer_placement(p);
            hdr = niftiinfo(fullfile(anat, 'sub-001_T1w.nii'));
            for ti = 1:2
                M = tc.Poses{3 - ti};
                expected_ras = M(1:3, 4)' + p.placement.localite.tracker_to_bowl_mm * M(1:3, 1)';
                tc.verifyEqual(out.transducer(ti).trans_pos, ras_to_grid(expected_ras, hdr), ...
                    'AbsTol', 1, sprintf('transducer %d must use pose %d', ti, 3 - ti));
            end
            tc.verifyNotEqual(out.transducer(1).trans_pos, out.transducer(2).trans_pos);

            p.placement.localite.position = [1 2 1];
            tc.verifyError(@() preproc_transducer_placement(p), 'PRESTUS:localite:listLength');
        end

        function test_deprecated_path_alias(tc)
            p = struct('path', struct('localite', '/data/localite'));
            paths = tc.verifyWarning(@() neuronav_resolve_paths(p), ...
                'PRESTUS:localite:deprecatedPath');
            tc.verifyEqual(paths.raw, '/data/localite');
            tc.verifyEqual(paths.t1, '/data/localite');
            p.path.localite_raw = '/data/raw';
            p.path.localite_t1  = '/data/t1';
            paths = tc.verifyWarningFree(@() neuronav_resolve_paths(p));
            tc.verifyEqual(paths.raw, '/data/raw');
            tc.verifyEqual(paths.t1, '/data/t1');
        end

        function test_missing_source_errors(tc)
            p = base_parameters();
            tc.verifyError(@() neuronav_localite_pose(p), 'PRESTUS:localite:noSource');
        end
    end

    methods
        function f = write_export(tc, raw, sub, ses, loc)
            % Two tracked series of 80 pulses, 1 s apart, 60 s gap between series.
            d = fullfile(raw, sub, ses, loc, 'Session_20240101120000000', 'TMSTrigger');
            mkdir(d);
            f = fullfile(d, 'TriggerMarkers_Coil0_20240101123000000.xml');
            fid = fopen(f, 'w');
            fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n<TriggerMarkerList coordinateSpace="RAS">\n');
            t = 0;
            for s = 1:2
                M = tc.Poses{s};
                for k = 1:80
                    t = t + 1000;
                    fprintf(fid, '  <TriggerMarker description="" recordingTime="%d">\n    <Matrix4D', t);
                    for r = 0:3
                        for c = 0:3
                            fprintf(fid, ' data%d%d="%.10g"', r, c, M(r+1, c+1));
                        end
                    end
                    fprintf(fid, '/>\n  </TriggerMarker>\n');
                end
                t = t + 60000;
            end
            fprintf(fid, '</TriggerMarkerList>\n');
            fclose(fid);
        end
    end
end

function p = base_parameters()
    p = struct();
    p.subject_id = 1;
    p.placement.localite = struct('position', 1, 'tracker_to_bowl_mm', -10);
    p.transducer = struct('focal_distance_ep', 60, 'focal_distance_bowl', [], ...
        'focal_distance_offset', []);
    p.io.overwrite_files = 'always';
end

function R = rotz(deg)
    c = cosd(deg); s = sind(deg);
    R = [c -s 0; s c 0; 0 0 1];
end
