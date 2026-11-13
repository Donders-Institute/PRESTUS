classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_placement_paths < matlab.unittest.TestCase
% TEST_PLACEMENT_PATHS  Placement output locations and run-record reuse.
%
% GET_PLACEMENT_DIR / PLACEMENT_RUN_DIR resolve <path.placement or
% path.sim/placement>/sub-XXX/<mode>/<target_name>; PREPROC_TRANSDUCER_PLACEMENT
% reuses a run record (placement_run.json, or the legacy placement_run.mat)
% per placement.reuse: auto reuses a matching record and runs again on a
% mismatch (an error under io.overwrite_files: never), always requires a
% record, never ignores it.
%
% Run with:   results = runtests('tests/test_placement_paths.m');

    properties
        TmpDir
    end

    methods (TestMethodSetup)
        function setup(tc)
            tc.TmpDir = tempname;
            mkdir(tc.TmpDir);
            tc.addTeardown(@() rmdir(tc.TmpDir, 's'));
        end
    end

    methods (Test)
        function test_default_root_is_under_path_sim(tc)
            p = base_params(tc.TmpDir);
            tc.verifyEqual(get_placement_dir(p), fullfile(tc.TmpDir, 'sims', 'placement', 'sub-007'));
            tc.verifyEqual(placement_run_dir(p, 'tus_entry'), ...
                fullfile(tc.TmpDir, 'sims', 'placement', 'sub-007', 'tus_entry', 'DLPFC'));
        end

        function test_path_placement_wins(tc)
            p = base_params(tc.TmpDir);
            p.path.placement = fullfile(tc.TmpDir, 'shared');
            tc.verifyEqual(placement_run_dir(p, 'tus_entry'), ...
                fullfile(tc.TmpDir, 'shared', 'sub-007', 'tus_entry', 'DLPFC'));
        end

        function test_work_dir_override_and_default_target(tc)
            p = base_params(tc.TmpDir);
            p.placement.tus_entry.work_dir = fullfile(tc.TmpDir, 'mine');
            tc.verifyEqual(placement_run_dir(p, 'tus_entry'), fullfile(tc.TmpDir, 'mine'));
            p.placement.plantus = struct();
            tc.verifyEqual(placement_run_dir(p, 'plantus'), ...
                fullfile(tc.TmpDir, 'sims', 'placement', 'sub-007', 'plantus', 'target'));
        end

        function test_errors(tc)
            p = base_params(tc.TmpDir);
            tc.verifyError(@() placement_run_dir(p, 'manual'), 'PRESTUS:placement:noRunDir');
            p.path = rmfield(p.path, 'sim');           % falls back to io.dir_output's parent
            tc.verifyEqual(get_placement_dir(p), fullfile(tc.TmpDir, 'sims', 'placement', 'sub-007'));
            p.io = rmfield(p.io, 'dir_output');
            tc.verifyError(@() get_placement_dir(p), 'get_placement_dir:missingPath');
        end

        function test_record_reused_under_never(tc)
            p = base_params(tc.TmpDir);
            write_record(p, placement_signature(p, 'tus_entry'));
            p.io.overwrite_files = 'never';
            % script_path / env_path are unset: running TUS_entry would fail,
            % so a result proves the record was reused
            p = preproc_transducer_placement(p);
            tc.verifyEqual(p.transducer(1).trans_pos, [10 20 30]);
            tc.verifyEqual(p.transducer(1).focus_pos, [40 50 60]);
            tc.verifyEqual(p.transducer(1).focus_pos_ras, [4 5 6]);
        end

        function test_record_mismatch_errors(tc)
            p = base_params(tc.TmpDir);
            write_record(p, placement_signature(p, 'tus_entry'));
            p.io.overwrite_files = 'never';
            p.placement.tus_entry.maximal_distance_cm = 6;
            tc.verifyError(@() preproc_transducer_placement(p), 'prestus:cache:signatureMismatch');
        end

        function test_json_record_round_trip(tc)
            p = base_params(tc.TmpDir);
            q = p;
            q.transducer = struct('curv_radius_mm', {63, 63}, ...
                'trans_pos', {[10 20 30], [11 21 31]}, 'focus_pos', {[40 50 60], [41 51 61]}, ...
                'trans_pos_ras', {[1 2 3], [1 2 4]}, 'focus_pos_ras', {[4 5 6], [4 5 7]});
            q.placement.tus_entry.target_ras = [7 8 9];
            file = write_placement_record(q, 'tus_entry', placement_run_dir(p, 'tus_entry'), ...
                placement_signature(p, 'tus_entry'));
            tc.verifyTrue(isfile(file));
            rec = read_placement_record(placement_run_dir(p, 'tus_entry'));
            tc.verifyEqual(rec.trans_pos, [10 20 30; 11 21 31]);
            tc.verifyEqual(rec.focus_pos_ras, [4 5 6; 4 5 7]);
            tc.verifyEqual(rec.target_ras, [7 8 9]);
        end

        function test_json_record_reused_by_default(tc)
            % placement.reuse defaults to auto: reused under overwrite_files: always too
            p = base_params(tc.TmpDir);
            write_json_record(p, placement_signature(p, 'tus_entry'));
            p = preproc_transducer_placement(p);
            tc.verifyEqual(p.transducer(1).trans_pos, [10 20 30]);
            tc.verifyEqual(p.transducer(1).focus_pos_ras, [4 5 6]);
        end

        function test_auto_runs_again_on_mismatch(tc)
            p = base_params(tc.TmpDir);
            write_json_record(p, placement_signature(p, 'tus_entry'));
            p.placement.tus_entry.maximal_distance_cm = 6;
            % TUS_entry is not installed: reaching it proves the record was not reused
            id = run_placement_id(p);
            tc.verifyNotEmpty(id);
            tc.verifyNotEqual(id, 'prestus:cache:signatureMismatch');
        end

        function test_never_ignores_record(tc)
            p = base_params(tc.TmpDir);
            write_json_record(p, placement_signature(p, 'tus_entry'));
            p.placement.reuse = 'never';
            tc.verifyNotEmpty(run_placement_id(p));
        end

        function test_always_requires_record(tc)
            p = base_params(tc.TmpDir);
            p.placement.reuse = 'always';
            tc.verifyError(@() preproc_transducer_placement(p), 'PRESTUS:placement:noRecord');
        end

        function test_always_refuses_mismatch(tc)
            p = base_params(tc.TmpDir);
            write_json_record(p, placement_signature(p, 'tus_entry'));
            p.placement.reuse = 'always';
            p.placement.tus_entry.maximal_distance_cm = 6;
            tc.verifyError(@() preproc_transducer_placement(p), 'prestus:cache:signatureMismatch');
        end

        function test_bad_reuse_value(tc)
            p = base_params(tc.TmpDir);
            p.placement.reuse = 'sometimes';
            tc.verifyError(@() preproc_transducer_placement(p), 'PRESTUS:placement:badReuse');
        end
    end
end

function id = run_placement_id(p)
% Identifier of the error preproc_transducer_placement stops with ('' if none).
    id = '';
    try
        preproc_transducer_placement(p);
    catch ME
        id = ME.identifier;
        if isempty(id); id = 'unidentified'; end
    end
end

function write_json_record(p, signature)
    q = p;
    q.transducer.trans_pos = [10 20 30];
    q.transducer.focus_pos = [40 50 60];
    q.transducer.trans_pos_ras = [1 2 3];
    q.transducer.focus_pos_ras = [4 5 6];
    write_placement_record(q, 'tus_entry', placement_run_dir(p, 'tus_entry'), signature);
end

function p = base_params(root)
    p = struct();
    p.subject_id = 7;
    p.path.sim  = fullfile(root, 'sims');
    p.path.seg  = fullfile(root, 'seg');
    p.path.anat = fullfile(root, 'anat');            % no T1: the overlay figure is skipped
    p.path.t1_pattern = 'sub-%03d_T1w.nii';
    p.io.dir_output = fullfile(root, 'sims', 'sub-007');
    p.io.output_affix = '';
    p.io.overwrite_files = 'always';
    p.placement.mode = 'tus_entry';
    p.placement.report = false;
    p.placement.tus_entry = struct('target_name', 'DLPFC', 'maximal_distance_cm', 7.46);
    p.transducer = struct('curv_radius_mm', 63, 'trans_pos', [], 'focus_pos', []);
end

function write_record(p, signature)
    run_dir = placement_run_dir(p, 'tus_entry');
    mkdir(run_dir);
    rec = struct('mode', 'tus_entry', 'subject_id', p.subject_id, 'created', 'test', ...
        'signature', signature, 'trans_pos', [10 20 30], 'focus_pos', [40 50 60], ...
        'trans_pos_ras', [1 2 3], 'focus_pos_ras', [4 5 6], 'target_ras', []);
    save(fullfile(run_dir, 'placement_run.mat'), '-struct', 'rec');
end
