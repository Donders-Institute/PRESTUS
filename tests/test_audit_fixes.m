classdef test_audit_fixes < matlab.unittest.TestCase
% TEST_AUDIT_FIXES  Regression tests for defects found in the 2026-09/10 code audit.
%
%   No SimNIBS, k-Wave or scheduler needed. HPC tests only write the job
%   script; submission fails without sbatch/qsub, which is expected.
%   Run with:   results = runtests('tests/test_audit_fixes.m');

    properties
        TmpDir
        RepoRoot
    end

    methods (TestClassSetup)
        function add_prestus_path(tc)
            tc.RepoRoot = fileparts(fileparts(mfilename('fullpath')));
            addpath(genpath(fullfile(tc.RepoRoot, 'functions')));
            addpath(genpath(fullfile(tc.RepoRoot, 'external')));
        end
    end

    methods (TestMethodSetup)
        function make_tmp(tc)
            tc.TmpDir = tempname;
            mkdir(tc.TmpDir);
        end
    end

    methods (TestMethodTeardown)
        function remove_tmp(tc)
            if isfolder(tc.TmpDir); rmdir(tc.TmpDir, 's'); end
        end
    end

    methods (Test, TestTags = {'audit'})

        function test_tutorial_config_positions_at_transducer_level(tc)
            cfg = yaml.loadFile(fullfile(tc.RepoRoot, 'config', 'config_tutorial.yaml'), ...
                'ConvertToArray', true);
            tc.verifyFalse(isfield(cfg.transducer, 'position'));
            tc.verifyEqual(double(cfg.transducer.trans_pos), [128 139 15]);
            tc.verifyEqual(double(cfg.transducer.focus_pos), [128 142 75]);
            f = cfg.transducer.matrix.matrix_shape.extract_from_file.file_path;
            tc.verifyTrue(isfile(fullfile(tc.RepoRoot, 'config', 'equipment', f)));
        end

        function test_affix_desc_prefix_not_doubled(tc)
            tc.verifyEqual(sanitize_output_affix('desc-foo'), '_desc-foo');
            tc.verifyEqual(sanitize_output_affix('_foo'), '_desc-foo');
            tc.verifyEqual(sanitize_output_affix('_desc-foo'), '_desc-foo');
            tc.verifyEqual(sanitize_output_affix('__'), '');
        end

        function test_matlab_escape_quotes_round_trip(tc)
            p = '/data/O''Brien/sub-001; delete(''x'')';
            tc.verifyEqual(eval(['''' matlab_escape_quotes(p) '''']), p);
        end

        function test_phase_expression_arithmetic_only(tc)
            tr.annular = struct('elem_id_mm', [0 10], 'elem_od_mm', [9 20], 'elem_n', 2, ...
                'curv_radius_mm', 60, 'dist_geom_ep_mm', 50, 'elem_p_amp', 1);
            tr.annular.elem_phase_deg = {0, '274*pi/180'};
            out = validate_annular_transducer(tr, 1);
            tc.verifyEqual(out.annular.elem_phase_deg, [0 274*pi/180], 'AbsTol', 1e-12);

            sentinel = fullfile(tc.TmpDir, 'PWNED');
            tr.annular.elem_phase_deg = {0, sprintf('fclose(fopen(''%s'',''w''))', sentinel)};
            tc.verifyError(@() validate_annular_transducer(tr, 1), ...
                'validate_annular_transducer:badPhaseExpression');
            tc.verifyFalse(isfile(sentinel));
        end

        function test_manifest_keeps_labels_and_serials_as_text(tc)
            f = fullfile(tc.TmpDir, 'manifest.tsv');
            fid = fopen(f, 'w');
            fprintf(fid, '# comment line\n');
            fprintf(fid, 'subject_id\tplacement.localite.target\ttransducer.serial\ttransducer.target_isppa_wcm2\n');
            fprintf(fid, '1\tno\t001\t30\n');
            fclose(fid);
            m = study_read_manifest(f);
            ov = m.overrides{1};
            tc.verifyEqual(ov{1, 2}, 'no');
            tc.verifyEqual(ov{2, 2}, '001');
            tc.verifyEqual(ov{3, 2}, 30);
        end

        function test_manifest_error_reports_file_line(tc)
            f = fullfile(tc.TmpDir, 'manifest.tsv');
            fid = fopen(f, 'w');
            fprintf(fid, '# comment\n\nsubject_id\tcondition\n1\tA\textra\n');
            fclose(fid);
            try
                study_read_manifest(f);
                tc.verifyFail('Expected tooManyCells error');
            catch ME
                tc.verifyEqual(ME.identifier, 'study_read_manifest:tooManyCells');
                tc.verifySubstring(ME.message, 'line 4');
            end
        end

        function test_script_json_cannot_close_script_block(tc)
            j = html_utils.script_json(struct('label', '</script><script>alert(1)</script>'));
            tc.verifyEmpty(strfind(j, '</'));
            tc.verifyEqual(jsondecode(j).label, '</script><script>alert(1)</script>');
            tc.verifyEqual(html_utils.escape('a''b'), 'a&#39;b');
        end

        function test_native_mni_warp_keeps_voxel_position(tc)
            % Identity warp on a 20^3 MNI grid (origin -10 mm); input grid
            % with origin -5 mm. A voxel at world (2,3,4) mm must land on
            % MNI voxel (13,14,15), 1-based.
            m2m = fullfile(tc.TmpDir, 'm2m_sub-001');
            mkdir(fullfile(m2m, 'toMNI'));
            n = 20;
            [i, j, k] = ndgrid(0:n-1, 0:n-1, 0:n-1);
            warp = single(cat(4, i - 10, j - 10, k - 10));
            write_nii(fullfile(m2m, 'toMNI', 'MNI2Conform_nonl'), warp, -10);
            write_nii(fullfile(m2m, 'toMNI', 'final_tissues_MNI'), zeros(n, n, n, 'uint8'), -10);
            vol = zeros(n, n, n, 'single');
            vol(8, 9, 10) = 1;
            in = write_nii(fullfile(tc.TmpDir, 'in'), vol, -5);
            out = fullfile(tc.TmpDir, 'out_MNI.nii.gz');
            convert_final_to_MNI_simnibs(in, m2m, out, struct(), ...
                'method', 'native', 'interpolation_order', 0);
            res = niftiread(out);
            [~, idx] = max(res(:));
            [a, b, c] = ind2sub(size(res), idx);
            tc.verifyEqual([a b c], [13 14 15]);
        end

        function test_qsub_script_walltime_and_cpu_request(tc)
            p = hpc_params();
            p.simulation.code_type = 'matlab_cpu';
            try %#ok<TRYNC>
                hpc_submit_job('qsub', 'job_m', p, tc.TmpDir);
            end
            txt = fileread(only_file(tc, 'temp_qsub_*.sh'));
            tc.verifySubstring(txt, 'walltime=04:00:00');
            tc.verifyEmpty(strfind(txt, 'gpus='));
        end

        function test_snellius_timelimit_check(tc)
            p = hpc_params();
            p.hpc.name = 'snellius';
            p.hpc.partition = 'gpu_a100';
            p.snellius.gpu_a100 = struct('memorylimit', 120, 'cores', 18, ...
                'timelimit', '120:00:00', 'n_gpu', 1);
            try %#ok<TRYNC>
                hpc_submit_job('slurm', 'job_m', p, tc.TmpDir);
            end
            tc.verifySubstring(fileread(only_file(tc, 'temp_slurm_*.sh')), '--time=04:00:00');

            p.hpc.timelimit = '121:00:00';
            tc.verifyError(@() hpc_submit_job('slurm', 'job_m', p, tc.TmpDir), ?MException);
        end

    end
end

function p = hpc_params()
    p.subject_id = 1;
    p.hpc = struct('name', 'default', 'partition', '', 'gpu', '', 'reservation', '', ...
        'timelimit', '04:00:00', 'memorylimit', 20, 'job_prefix', 'PS');
    p.simulation.code_type = 'matlab_gpu';
    p.io.output_affix = '';
end

function f = write_nii(base, vol, origin_mm)
% 1 mm isotropic NIfTI whose 0-based voxel [0 0 0] sits at origin_mm on each axis.
    niftiwrite(vol, base);
    info = niftiinfo([base '.nii']);
    T = eye(4);
    T(4, 1:3) = origin_mm;
    info.Transform = affine3d(T);
    info.TransformName = 'Sform';   % 'None' would drop the affine on write
    niftiwrite(vol, base, info, 'Compressed', true);
    delete([base '.nii']);
    f = [base '.nii.gz'];
end

function f = only_file(tc, pattern)
    d = dir(fullfile(tc.TmpDir, pattern));
    tc.assertNumElements(d, 1, sprintf('Expected one %s', pattern));
    f = fullfile(d.folder, d.name);
end
