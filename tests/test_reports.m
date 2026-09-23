classdef test_reports < matlab.unittest.TestCase
% TEST_REPORTS  Unit/smoke tests for the HTML report generators.
%
%   Covers the report enhancements: the new safety metrics (TIC, ΔT-from-37),
%   the per-dimension image fix, the embedded config dump, and end-to-end
%   generation of both the per-subject and group reports from synthetic data
%   (no MRI / SimNIBS / GPU). The generators are pure HTML string builders
%   wrapped in per-section try/catch, so a successful run is asserted by the
%   absence of any "Section failed" marker plus the presence of the new
%   feature markup.
%
%   Run with:  results = runtests('tests/test_reports.m');

    properties
        tmp   % per-test temporary output directory
    end

    methods (TestMethodSetup)
        function setup(tc)
            tc.tmp = fullfile(tempdir, 'prestus_report_unittest');
            if exist(tc.tmp, 'dir'); rmdir(tc.tmp, 's'); end
            mkdir(tc.tmp);
        end
    end

    methods (TestMethodTeardown)
        function teardown(tc)
            if exist(tc.tmp, 'dir'); rmdir(tc.tmp, 's'); end
        end
    end

    methods (Test, TestTags = {'reports'})

        function test_risk_limits_has_new_metrics(tc)
            L = get_risk_limits(true);
            for f = {'riseT37_brain', 'riseT37_skull', 'riseT37_skin', 'TIC'}
                tc.verifyTrue(isfield(L, f{1}), sprintf('limit %s missing', f{1}));
            end
            tc.verifyEqual(L.riseT37_skull.limit, 2.0, 'AbsTol', 1e-9);
            tc.verifyTrue(isinf(L.TIC.limit), 'TIC should be informational (Inf limit)');
        end

        function test_compute_tic_valid(tc)
            tr = local_transducer();
            p  = struct('timing', struct('pd', 0.02, 'pri', 0.1));
            [v, info] = compute_tic(p, tr, 20);
            tc.verifyTrue(isfinite(v) && v > 0, 'TIC should be finite and positive');
            tc.verifyEqual(info.duty_cycle, 0.2, 'AbsTol', 1e-9);
            tc.verifyTrue(isfinite(info.Deq_cm) && info.Deq_cm > 0);
        end

        function test_compute_tic_missing_aperture_is_nan(tc)
            tr = struct('type', 'annular', 'annular', struct('curv_radius_mm', 63));
            p  = struct('timing', struct('pd', 0.02, 'pri', 0.1));
            v  = compute_tic(p, tr, 20);
            tc.verifyTrue(isnan(v), 'TIC should be NaN when aperture is unavailable');
        end

        function test_embed_image_dims_finds_slices(tc)
            base = 'sub-001_layered_maxT';
            suffix = '_test';
            for d = {'_x', '_y', '_z'}
                imwrite(uint8(zeros(4, 4, 3)), fullfile(tc.tmp, [base d{1} suffix '.png']));
            end
            html = html_utils.embed_image_dims(tc.tmp, base, suffix, 'alt', 'cap');
            tc.verifyEqual(numel(strfind(html, '<figure>')), 3, 'should embed all three slice images');
        end

        function test_embed_image_dims_no_dim_file(tc)
            base = 'sub-001_layered_maxT';
            suffix = '_test';
            imwrite(uint8(zeros(4, 4, 3)), fullfile(tc.tmp, [base suffix '.png']));
            html = html_utils.embed_image_dims(tc.tmp, base, suffix, 'alt', 'cap');
            tc.verifyEqual(numel(strfind(html, '<figure>')), 1, '2D no-dim file should embed once');
        end

        function test_config_dump_prunes_large_arrays(tc)
            p = struct('a', 1, 'big', rand(200, 1), 'nested', struct('b', 'hello'));
            html = build_config_dump(p);
            tc.verifyTrue(~isempty(html));
            tc.verifyTrue(~isempty(strfind(html, 'pruned')), 'large array should be pruned'); %#ok<STREMP>
        end

        function test_subject_report_generates(tc)
            p = local_make_subject_params(tc.tmp);
            rp = generate_simulation_report(p);
            tc.verifyTrue(~isempty(rp) && isfile(rp), 'no report produced');
            txt = fileread(rp);
            tc.verifyEqual(numel(strfind(txt, 'Section failed')), 0, 'a section builder errored');
            for m = {'<!DOCTYPE html>', 'id="safety"', 'id="other"', 'Mechanical', 'class="verdict', ...
                     'Cranial TI (TIC)', 'EP-to-focus', 'config-dump', 'theme-btn'}
                tc.verifyTrue(~isempty(strfind(txt, m{1})), sprintf('missing %s', m{1})); %#ok<STREMP>
            end
            % Structural split: an MI tile sits in the Safety section, TIC in Other.
            iSafety = strfind(txt, 'id="safety"'); iOther = strfind(txt, 'id="other"');
            iMI = strfind(txt, 'MI (transcranial)'); iTIC = strfind(txt, 'Cranial TI (TIC)');
            tc.verifyTrue(iMI(1) > iSafety(1) && iMI(1) < iOther(1), 'MI should be in the Safety section');
            tc.verifyTrue(iTIC(1) > iOther(1), 'TIC should be in the Other section');
        end

        function test_group_report_generates_and_flags(tc)
            ids = local_make_group_dir(tc.tmp);
            p = struct('path', struct('sim', tc.tmp), ...
                       'simulation', struct('medium', 'layered'), ...
                       'io', struct('output_affix', '_test'), ...
                       'modules', struct('run_heating_sims', 1));
            rp = generate_group_report(p, ids);
            tc.verifyTrue(~isempty(rp) && isfile(rp), 'no group report produced');
            txt = fileread(rp);
            tc.verifyEqual(numel(strfind(txt, 'Section failed')), 0, 'a section builder errored');
            tc.verifyTrue(~isempty(strfind(txt, 'Data Analysis')), 'analysis section missing'); %#ok<STREMP>
            tc.verifyTrue(~isempty(strfind(txt, 'PRESTUS_SUBJECTS')), 'analysis payload missing'); %#ok<STREMP>
            tc.verifyEqual(numel(strfind(txt, 'class="exceed"')), 1, 'exactly one subject should be flagged');
            % Structural split: Safety dashboard vs Other-metrics section.
            tc.verifyTrue(~isempty(strfind(txt, 'id="other-group"')), 'Other metrics section missing'); %#ok<STREMP>
            iDash = strfind(txt, 'id="dashboard"'); iOther = strfind(txt, 'id="other-group"');
            iMI = strfind(txt, 'MI (transcranial)'); iTIC = strfind(txt, 'Cranial TI (TIC)');
            tc.verifyTrue(iMI(1) > iDash(1) && iMI(1) < iOther(1), 'MI should be in the Safety dashboard');
            tc.verifyTrue(iTIC(1) > iOther(1), 'TIC should be in the Other section');
        end

    end
end

% =========================================================================
%  Local helpers
% =========================================================================
function tr = local_transducer()
    tr = struct('type', 'annular', 'name', 'CTX-500-026', 'freq_hz', 5e5);
    tr.annular = struct('elem_n', 4, 'curv_radius_mm', 63.2, 'elem_amp', 1.06e5, ...
        'elem_od_mm', [0 33.02 46.23 56.08], 'dist_geom_ep_mm', 51, 'elem_id_mm', [0 30 43 53]);
end

function row = local_csv_row(id, riseT_skull)
% One synthetic per-subject CSV row with the full safety-metric column set.
    row = table();
    row.subject_id = id; row.freq_Hz = 5e5;
    row.Isppa = 22.1; row.Ipa_target = 19.8;
    row.real_focal_distance_mm = 64.2; row.ep_focal_distance_mm = 12.8; row.ep_offset_mm = 51;
    row.MI_tc = 0.74; row.MI_brain = 0.71; row.MI_skull = 0.95; row.MI_skin = 0.42; row.TIC = 0.83;
    row.Isppa_brain = 21.4; row.Isppa_skull = 9.8; row.Isppa_skin = 9.0;
    row.maxT_skull = 38.8; row.maxT_brain = 38.2; row.maxT_skin = 37.4;
    row.riseT_skull = riseT_skull; row.riseT_brain = 1.2; row.riseT_skin = 0.4;
    row.riseT37_skull = riseT_skull; row.riseT37_brain = 1.2; row.riseT37_skin = 0.4;
    row.CEM43_skull = 11.4; row.CEM43_brain = 0.8; row.CEM43_skin = 5.9;
    row.CEM43iso_skull = 12; row.CEM43iso_brain = 0.9; row.CEM43iso_skin = 6;
end

function p = local_make_subject_params(out)
    csv = fullfile(out, 'sub-601_layered_test.csv');
    writetable(local_csv_row(601, 1.9), csv);

    p = struct();
    p.subject_id = 601;
    p.simulation = struct('medium', 'layered', 'debug', 0, 'code_type', 'matlab_gpu');
    p.io = struct('output_affix', '_test', 'dir_output', out, 'filename_table', csv);
    p.modules = struct('run_source_setup', 1, 'run_acoustic_sims', 1, 'run_heating_sims', 1, ...
        'run_posthoc_water_sims', 0, 'generate_report', 1);
    p.thermal = struct('cem43_iso', 1, 'pd', 0.02, 'pri', 0.1, 'ptd', 0.4, 'ptri', 0.8, 'ptrd', 0.4, 'post_ptri_dur', 1);
    p.thermal.temp_0 = struct('brain', 37, 'skull', 37, 'skin', 35.5);
    p.transducer = local_transducer();
    p.transducer.trans_pos = [100 100 30];
    p.transducer.focus_pos = [100 100 130];
    p.grid = struct('resolution_mm', 0.5, 'n_dims', 3, 'default_dims', [201 201 201], 'pml_size', 10, 'dims', [201 201 201]);
    p.timing = struct('pd', 0.02, 'pri', 0.1, 'dc', 0.2, 'ptd', 0.4, 'ptri', 0.8);
    p.analysis = struct('focus_area_radius', 5);
    p.layers = struct('skin', 1, 'skull', 2, 'brain', 3);
    p.medium_properties = struct();
    p.medium_properties.skin  = struct('sound_speed', 1610, 'density', 1090, 'thermal_conductivity', 0.4, 'specific_heat', 3500);
    p.medium_properties.skull = struct('sound_speed', 2800, 'density', 1900, 'thermal_conductivity', 0.32, 'specific_heat', 1300);
    p.medium_properties.brain = struct('sound_speed', 1560, 'density', 1040, 'thermal_conductivity', 0.51, 'specific_heat', 3630);
    p.pct = struct('enabled', 0);
end

function ids = local_make_group_dir(sim)
    ids = [601, 998, 11];           % sub-998 exceeds (riseT_skull = 2.3 > 2.0)
    rises = [1.9, 2.3, 1.4];
    for k = 1:numel(ids)
        d = fullfile(sim, sprintf('sub-%03d', ids(k)));
        if ~exist(d, 'dir'); mkdir(d); end
        row = local_csv_row(ids(k), rises(k));
        writetable(row, fullfile(d, sprintf('sub-%03d_layered_test.csv', ids(k))));
    end
end
