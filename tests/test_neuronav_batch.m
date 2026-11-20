classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_neuronav_batch < matlab.unittest.TestCase
% TEST_NEURONAV_BATCH  Multi-subject Localite ingestion: scan, label map, naming, batch.
%
% Synthetic study in a temp folder: 2 subjects × 2 sessions, each with a
% TriggerMarkers export (two tracked series) and an EntryTarget export whose
% targets lie on the recorded beam axes, labelled with two conventions
% ('LEFT PUL'/'RIGHT PUL' vs. 'L_PUL_1'/'R_PUL_1').

    properties
        Root
        Poses
    end

    methods (TestMethodSetup)
        function make_study(tc)
            tc.Root = tempname;
            mkdir(tc.Root);
            tc.addTeardown(@() rmdir(tc.Root, 's'));
            tc.Poses = {[rotz(20),  [ 40; 90; 140]; 0 0 0 1], ...
                        [rotz(160), [140; 95; 135]; 0 0 0 1]};
            conventions = {{'LEFT PUL', 'RIGHT PUL'}, {'L_PUL_1', 'R_PUL_1'}};
            subs = {'sub-001', 'sub-002'};
            for i = 1:2
                for j = 1:2
                    ses = sprintf('ses-%02d', j);
                    % shift each session a little so files differ (no de-duplication)
                    tc.write_export(subs{i}, ses, [i; j; 0], conventions{j});
                end
            end
        end
    end

    methods (Test, TestTags = {'unit', 'localite'})

        function test_map_label(tc)
            map = struct('target', {'L_PUL', 'R_PUL'}, 'labels', {{'LEFT PUL', 'l_pul*'}, {'RIGHT*'}});
            tc.verifyEqual(neuronav_map_label('LEFT PUL', map), 'L_PUL');
            tc.verifyEqual(neuronav_map_label('L_PUL_1', map), 'L_PUL');   % case-insensitive wildcard
            tc.verifyEqual(neuronav_map_label('right pul', map), 'R_PUL');
            [t, m] = neuronav_map_label('other', map);
            tc.verifyEqual(t, 'other'); tc.verifyFalse(m);
            tc.verifyEqual(neuronav_map_label('LEFT PUL', []), 'LEFT PUL');
        end

        function test_scan_lists_sessions_and_labels(tc)
            p = tc.params();
            [sessions, labels] = neuronav_scan_sessions(p);
            tc.verifyEqual(height(sessions), 4);
            tc.verifyEqual(sessions.n_trigger', [1 1 1 1]);
            tc.verifyEqual(sessions.n_entrytarget', [1 1 1 1]);
            tc.verifyEqual(sort(labels.label)', sort({'LEFT PUL', 'RIGHT PUL', 'L_PUL_1', 'R_PUL_1'}));
            tc.verifyEqual(labels.stem{strcmp(labels.label, 'L_PUL_1')}, 'L_PUL');
            tc.verifyEqual(labels.n_subjects', [2 2 2 2]);
            map = neuronav_label_map_from_scan(labels, 'stem');
            tc.verifyEqual(numel(map), 4);   % LEFT PUL, RIGHT PUL, L_PUL, R_PUL
        end

        function test_batch_names_via_label_map(tc)
            p = tc.params();
            p.neuronav.target_labels = struct('target', {'L_PUL', 'R_PUL'}, ...
                'labels', {{'LEFT PUL', 'L_PUL*'}, {'RIGHT PUL', 'R_PUL*'}});
            T = neuronav_ingest_batch(p);
            tc.verifyEqual(height(T), 8);                       % 4 sessions × 2 targets
            tc.verifyEqual(sort(unique(T.target))', {'L_PUL', 'R_PUL'});
            tc.verifyTrue(all(strcmp(T.rule, 'heuristic')));
            tc.verifyLessThan(max(T.axis_dist_mm), 2);
            % coordinates: pose column 4 of the series matched to L_PUL
            r = T(strcmp(T.sub, 'sub-001') & strcmp(T.ses, 'ses-01') & strcmp(T.target, 'L_PUL'), :);
            tc.verifyEqual([r.pose_x r.pose_y r.pose_z], (tc.Poses{1}(1:3, 4) + [1; 1; 0])', 'AbsTol', 1e-6);
            tc.verifyTrue(isfile(fullfile(p.path.localite_post, 'neuronav_ingest_summary.csv')));
            tc.verifyTrue(isfile(fullfile(p.path.localite_post, 'sub-002', 'ses-02', ...
                'InstrumentMarker_sub-002_ses-02.json')));
        end

        function test_batch_without_map_uses_labels(tc)
            p = tc.params();
            T = neuronav_ingest_batch(p, {'sub-001', 'ses-02'}, 'summary_file', '');
            tc.verifyEqual(sort(T.target)', {'L_PUL_1', 'R_PUL_1'});
        end

        function test_summary_merge_replaces_sessions(tc)
            p = tc.params();
            neuronav_ingest_batch(p);
            neuronav_ingest_batch(p, {'sub-001', 'ses-01'});
            f = fullfile(p.path.localite_post, 'neuronav_ingest_summary.csv');
            T = readtable(f, 'TextType', 'char');
            tc.verifyEqual(height(T), 8, sprintf( ...
                'Re-running one session must not duplicate its rows. Summary file:\n%s', fileread(f)));
        end

        function test_named_mode_and_defaults(tc)
            p = tc.params();
            p.neuronav.suggest_names = false;
            T = neuronav_ingest_batch(p, {'sub-001', 'ses-01'}, 'summary_file', '');
            tc.verifyEqual(sort(T.target)', {'series_1', 'series_2'});
            p.neuronav.batch_write = 'named';
            T = neuronav_ingest_batch(p, {'sub-001', 'ses-01'}, 'summary_file', '');
            tc.verifyEqual(height(T), 1);
            tc.verifyTrue(isempty(T.target{1}));
            tc.verifySubstring(T.status{1}, 'no series to write');
        end

        function test_name_series_coil_map_wins(tc)
            p = tc.params();
            p.neuronav.coil_map = struct('coil', 0, 'sequences', struct('target', {'MINE', ''}));
            p.neuronav.coil_map.sequences(2).target = 'coil0';   % default name -> heuristic
            D = neuronav_detect_series(p, 'sub-001', 'ses-01');
            N = neuronav_name_series(D, p);
            tc.verifyEqual({N.rule}, {'coil_map', 'heuristic'});
            tc.verifyEqual(N(1).name, 'MINE');
            tc.verifyEqual(N(2).name, 'RIGHT PUL');
        end

        function test_summary_to_mni(tc)
            % Linear synthetic field (MNI = segmentation RAS + [5 0 0]), so
            % its interpolation is exact; planning and segmentation headers
            % differ, as for a Localite planning T1 vs. final_tissues.
            p = tc.params();
            p.path.seg = fullfile(tc.Root, 'seg');
            p.path.localite_t1 = fullfile(tc.Root, 't1');
            m2m = fullfile(p.path.seg, 'm2m_sub-001');
            mkdir(fullfile(m2m, 'toMNI'));
            mkdir(fullfile(p.path.localite_t1, 'sub-001'));
            A_plan = [eye(3), zeros(3, 1); -10 -20 -30 1];   % row convention, 0-based voxels
            A_seg  = [eye(3), zeros(3, 1); 100 50 60 1];
            write_nii(fullfile(p.path.localite_t1, 'sub-001', 'sub-001_T1'), zeros(20, 20, 20, 'single'), A_plan);
            write_nii(fullfile(m2m, 'final_tissues'), zeros(20, 20, 20, 'single'), A_seg);
            [i, j, k] = ndgrid(0:19);
            ras = [i(:) j(:) k(:) ones(numel(i), 1)] * A_seg;
            write_nii(fullfile(m2m, 'toMNI', 'Conform2MNI_nonl'), ...
                single(reshape(ras(:, 1:3) + [5 0 0], [20 20 20 3])), A_seg);

            T = table({'sub-001'; 'sub-002'}, {'ses-01'; 'ses-01'}, {'L_PUL'; 'L_PUL'}, ...
                [-5; -5], [-12; -12], [-25; -25], [-2; -2], [-15; -15], [-22; -22], ...
                [NaN; NaN], [NaN; NaN], [NaN; NaN], 'VariableNames', {'sub', 'ses', 'target', ...
                'pose_x', 'pose_y', 'pose_z', 'trans_ras_x', 'trans_ras_y', 'trans_ras_z', ...
                'focus_ras_x', 'focus_ras_y', 'focus_ras_z'});
            M = neuronav_summary_to_mni(T, p);
            % planning [-5 -12 -25] -> voxel [5 8 5] -> segmentation [105 58 65] -> +5 x
            tc.verifyEqual([M.pose_mni_x(1) M.pose_mni_y(1) M.pose_mni_z(1)], [110 58 65], 'AbsTol', 1e-4);
            tc.verifyEqual([M.trans_mni_x(1) M.trans_mni_y(1) M.trans_mni_z(1)], [113 55 68], 'AbsTol', 1e-4);
            tc.verifyTrue(all(isnan([M.focus_mni_x; M.focus_mni_y; M.focus_mni_z])));
            tc.verifyEmpty(M.mni_status{1});
            tc.verifySubstring(M.mni_status{2}, 'm2m_sub-002');
            tc.verifyTrue(isnan(M.pose_mni_x(2)));
        end
    end

    methods
        function p = params(tc)
            p = struct();
            p.path.localite_raw  = fullfile(tc.Root, 'raw');
            p.path.localite_post = fullfile(tc.Root, 'post');
            p.placement.localite = struct('loc_folder', 'loc', 'tracker_to_bowl_mm', -10);
            p.io.overwrite_files = 'always';
            p.transducer = struct('focal_distance_ep', 60, 'focal_distance_bowl', [], ...
                'focal_distance_offset', []);
        end

        function write_export(tc, sub, ses, shift, names)
            make_synthetic_localite_export(fullfile(tc.Root, 'raw'), sub, ses, tc.Poses, shift, names);
        end
    end
end

function R = rotz(deg)
    c = cosd(deg); s = sind(deg);
    R = [c -s 0; s c 0; 0 0 1];
end

function write_nii(base, vol, A)
% Write base.nii.gz with the row-convention affine A.
    niftiwrite(vol, base);
    info = niftiinfo([base '.nii']);
    delete([base '.nii']);
    info.Transform = affine3d(A);
    info.TransformName = 'Sform';   % 'None' would drop the affine
    niftiwrite(vol, base, info, 'Compressed', true);
end
