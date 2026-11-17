classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_stage_sidecars < matlab.unittest.TestCase
% TEST_STAGE_SIDECARS  Starting a run from prior stage outputs via their sidecars.
%
% Covers the sidecar helpers (WRITE_SIDECAR, READ_SIDECAR, FIND_SIDECAR,
% CHECK_SIDECAR_SIGNATURE), the pseudo-CT check (CHECK_PCT_SIDECAR, incl. the
% legacy pseudoCT.json), and the acoustic boundary: WRITE_ACOUSTIC_SIDECAR
% describes a run's pressure NIfTI and LOAD_ACOUSTIC_INPUT starts a later run
% from it (settings check, rescaling to a new target, FOV crop). Pure logic on
% small synthetic volumes: no k-Wave, no data.
%
%   results = runtests('tests/test_stage_sidecars.m');

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

        % ── helpers ──────────────────────────────────────────────────────
        function test_sidecar_round_trip_and_header(tc)
            p = params(tc.TmpDir);
            f = write_sidecar(fullfile(tc.TmpDir, 'x.json'), 'demo', ...
                struct('value', 3, 'signature', struct('a', [1 2 3])), p);
            sc = read_sidecar(f);
            tc.verifyEqual(string(sc.kind), "demo");
            tc.verifyEqual(sc.subject_id, 1);
            tc.verifyEqual(string(sc.generated_by.name), "PRESTUS");
            tc.verifyTrue(all(isfield(sc.generated_by, {'version', 'commit', 'commit_date', 'dirty'})));
            tc.verifyTrue(check_sidecar_signature(f, sc.signature, struct('a', [1 2 3]), p));
            tc.verifyError(@() check_sidecar_signature(f, sc.signature, struct('a', [1 2 4]), p), ...
                'prestus:cache:signatureMismatch');
        end

        function test_sidecar_uses_run_version(tc)
            % The run's recorded version wins over the current checkout.
            p = params(tc.TmpDir);
            p.io.prestus_version = struct('version', 'v9.9.9-1-gabc1234', 'commit', 'abc1234', ...
                'commit_date', '2026-01-02', 'branch', 'x', 'dirty', 'yes');
            f = write_sidecar(fullfile(tc.TmpDir, 'v.json'), 'demo', struct('value', 1), p);
            sc = read_sidecar(f);
            tc.verifyEqual(string(sc.generated_by.commit), "abc1234");
            tc.verifyEqual(string(sc.generated_by.commit_date), "2026-01-02");
            tc.verifyEqual(string(sc.generated_by.dirty), "yes");
        end

        function test_prestus_version_info(tc)
            [hash, ver, info] = prestus_version();
            tc.verifyTrue(all(isfield(info, {'version', 'commit', 'commit_date', 'branch', 'dirty'})));
            tc.verifyEqual(info.commit, hash);
            tc.verifyEqual(info.version, ver);
            tc.verifyTrue(ismember(info.dirty, {'yes', 'no', 'unknown'}));
        end

        function test_find_sidecar(tc)
            nii = fullfile(tc.TmpDir, 'sub-001_pressure.nii.gz');
            [found, target] = find_sidecar(nii);
            tc.verifyEmpty(found);
            tc.verifyEqual(target, fullfile(tc.TmpDir, 'sub-001_pressure.json'));
            fclose(fopen(target, 'w'));
            tc.verifyEqual(find_sidecar(nii), target);
        end

        % ── pseudo-CT ────────────────────────────────────────────────────
        function test_pct_sidecar_match_and_mismatch(tc)
            p = params(tc.TmpDir);
            p.pct = struct('skull_mapping', 'kosciessa');
            write_sidecar(fullfile(tc.TmpDir, 'pseudoCT.json'), 'pseudoCT', ...
                struct('signature', pct_signature(p)), p);
            tc.verifyTrue(check_pct_sidecar(tc.TmpDir, p));
            p.pct.skull_mapping = 'carpino';
            tc.verifyError(@() check_pct_sidecar(tc.TmpDir, p), 'prestus:cache:signatureMismatch');
        end

        function test_pct_legacy_sidecar(tc)
            % pseudoCT.json with only SkullMapping (no signature)
            p = params(tc.TmpDir);
            p.pct = struct('skull_mapping', 'carpino');
            fid = fopen(fullfile(tc.TmpDir, 'pseudoCT.json'), 'w');
            fprintf(fid, '{"SkullMapping": "carpino", "SubjectID": 1}');
            fclose(fid);
            tc.verifyTrue(check_pct_sidecar(tc.TmpDir, p));
            p.pct.skull_mapping = 'kosciessa';
            tc.verifyError(@() check_pct_sidecar(tc.TmpDir, p), 'prestus:cache:signatureMismatch');
        end

        function test_pct_without_sidecar_warns(tc)
            p = params(tc.TmpDir);
            p.pct = struct('skull_mapping', 'kosciessa');
            tc.verifyWarning(@() check_pct_sidecar(tc.TmpDir, p), 'prestus:cache:unverified');
        end

        % ── acoustic pressure map ────────────────────────────────────────
        function test_acoustic_round_trip_rescales_to_new_target(tc)
            p = params(tc.TmpDir);
            [nii, vol] = write_pressure(p);
            % Source run: field scaled from a 15 W/cm² free-water drive to 30 W/cm²
            prov = struct('freefield_isppa_wcm2', 15, 'freefield_source', 'kwave_simulation', ...
                'input_signature', cache_signature(p, 'acoustic'));
            f = write_acoustic_sidecar(p, prov, struct('isppa_scale_p', sqrt(2)));
            tc.verifyTrue(isfile(f));
            sc = read_sidecar(f);
            tc.verifyEqual(sc.isppa.level_wcm2, 30);
            tc.verifyEqual(sc.isppa.freefield_baseline_wcm2, 15);

            % Later run at 60 W/cm²: pressure × sqrt(2), 4× the drive in intensity
            q = p;
            q.transducer.target_isppa_wcm2 = 60;
            q.path.acoustic_input = nii;
            [sd, prov2] = load_acoustic_input(nii, q, planimg_for(vol), ...
                struct('input_signature', cache_signature(q, 'acoustic')));
            tc.verifyEqual(double(sd.p_max_all), vol * sqrt(2), 'RelTol', 1e-5);
            tc.verifyEqual(sd.isppa_scale_p, 2, 'RelTol', 1e-9);
            tc.verifyEqual(prov2.freefield_isppa_wcm2, 15);
            tc.verifyEqual(prov2.acoustic_input, nii);
        end

        function test_acoustic_settings_mismatch(tc)
            p = params(tc.TmpDir);
            [nii, vol] = write_pressure(p);
            write_acoustic_sidecar(p, struct('freefield_isppa_wcm2', 30, ...
                'input_signature', cache_signature(p, 'acoustic')), struct());
            q = p;
            q.medium_properties.skull.sound_speed = 3000;
            tc.verifyError(@() load_acoustic_input(nii, q, planimg_for(vol), ...
                struct('input_signature', cache_signature(q, 'acoustic'))), ...
                'prestus:cache:signatureMismatch');
            % Thermal-only grid settings do not invalidate the map
            q = p;
            q.grid.thermal_resolution_mm = 1;
            sd = load_acoustic_input(nii, q, planimg_for(vol), ...
                struct('input_signature', cache_signature(q, 'acoustic')));
            tc.verifyEqual(double(sd.p_max_all), vol, 'RelTol', 1e-6);
        end

        function test_acoustic_without_sidecar_is_used_as_is(tc)
            p = params(tc.TmpDir);
            p.transducer.target_isppa_wcm2 = [];
            [nii, vol] = write_pressure(p);
            sd = tc.verifyWarning(@() load_acoustic_input(nii, p, planimg_for(vol), struct()), ...
                'prestus:cache:unverified');
            tc.verifyEqual(double(sd.p_max_all), vol, 'RelTol', 1e-6);
        end

        function test_acoustic_fov_crop(tc)
            p = params(tc.TmpDir);
            [nii, vol] = write_pressure(p);
            img = planimg_for(vol);
            img.fov_offset_ac = [2 3 1];
            p.grid.dims = [3 2 4];
            p.transducer.target_isppa_wcm2 = [];
            p.io.cache_signature_check = 'off';     % no sidecar needed for the crop itself
            sd = load_acoustic_input(nii, p, img, struct('full_ac_dims', size(vol)));
            tc.verifyEqual(double(sd.p_max_all), vol(2:4, 3:4, 1:4), 'RelTol', 1e-6);
        end

        function test_acoustic_input_must_be_t1_space(tc)
            p = params(tc.TmpDir);
            [nii, vol] = write_pressure(p);
            img = planimg_for(vol);
            img.t1_image_orig = zeros(5, 5, 5);
            tc.verifyError(@() load_acoustic_input(nii, p, img, struct()), ...
                'prestus:acousticInput:spaceMismatch');
        end

        function test_deprecated_key_is_read(tc)
            p = params(tc.TmpDir);
            p.io.external_acoustic_nifti = 'a.nii.gz';
            tc.verifyEqual(acoustic_input_file(p), 'a.nii.gz');
            p.path.acoustic_input = 'b.nii.gz';
            tc.verifyEqual(acoustic_input_file(p), 'b.nii.gz');
        end
    end
end

function p = params(root)
    p = struct();
    p.subject_id = 1;
    p.simulation = struct('medium', 'layered', 'interactive', 0);
    p.grid = struct('resolution_mm', 1, 'dims', [6 5 4]);
    p.layers = struct('skull', 4);
    p.headmodel = struct();
    p.pct = struct();
    p.medium_properties = struct('skull', struct('sound_speed', 2800));
    p.transducer = struct('freq_hz', 500e3, 'trans_pos', [3 3 1], 'focus_pos', [3 3 4], ...
        'target_isppa_wcm2', 30);
    p.path = struct('sim', root);
    p.io = struct('overwrite_files', 'always', 'cache_signature_check', 'error', ...
        'output_affix', '', 'dir_output', root, 'dir_nii_T1w', fullfile(root, 'nii', 'T1w'));
end

function [nii, vol] = write_pressure(p)
    vol = reshape(1:prod(p.grid.dims), p.grid.dims) * 1e3;
    if ~isfolder(p.io.dir_nii_T1w); mkdir(p.io.dir_nii_T1w); end
    stem = fullfile(p.io.dir_nii_T1w, sprintf('sub-%03d_%s_T1w%s_pressure', ...
        p.subject_id, p.simulation.medium, p.io.output_affix));
    niftiwrite(single(vol), stem, 'Compressed', true);
    nii = [stem '.nii.gz'];
end

function img = planimg_for(vol)
    img = struct('transf', eye(4), 't1_image_orig', zeros(size(vol)), 'fov_offset_ac', [1 1 1]);
end
