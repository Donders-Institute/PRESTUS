classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_reports < matlab.unittest.TestCase
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

        function test_compute_tic_source_power(tc)
            % W0 from the source: sum over rings of (scale_p*elem_p_amp)^2/(2 rho c) * ring area
            % * duty cycle; ISPPA scaling enters squared; pulse-train duty multiplies in
            tr = local_transducer();
            tr.annular = rmfield(tr.annular, 'elem_amp');
            tr.annular.elem_p_amp = 1e5;
            water = struct('density', 1000, 'sound_speed', 1500);
            p = struct('timing', struct('pd', 0.02, 'pri', 0.1), ...
                'medium_properties', struct('water', water));
            od = tr.annular.elem_od_mm(:); id = tr.annular.elem_id_mm(:);
            P_W = sum((1e5)^2 / (2 * 1000 * 1500) * pi/4 * (od.^2 - id.^2) * 1e-6);
            [~, info] = compute_tic(p, tr, 20);
            tc.verifyEqual(info.W0_mW, P_W * 0.2 * 1000, 'RelTol', 1e-12);
            tc.verifySubstring(info.note, 'source intensity');
            [~, info2] = compute_tic(p, tr, 20, 2);
            tc.verifyEqual(info2.W0_mW, 4 * info.W0_mW, 'RelTol', 1e-12);
            p.timing.ptd = 20; p.timing.ptri = 40;
            [~, info3] = compute_tic(p, tr, 20);
            tc.verifyEqual(info3.duty_cycle, 0.1, 'AbsTol', 1e-12);
        end

        function test_compute_tic_ramp_only_on_source_power(tc)
            % elem_p_amp is a peak amplitude: the ramp factor lowers W0; the Isppa
            % fallback is already a pulse average and is not scaled again
            tr = local_transducer();
            tr.annular = rmfield(tr.annular, 'elem_amp');
            tr.annular.elem_p_amp = 1e5;
            p = struct('timing', struct('pd', 0.02, 'pri', 0.1));
            [~, flat] = compute_tic(p, tr, 20);
            p.timing.pulse_ramp_dur = 0.004;
            [~, ramped] = compute_tic(p, tr, 20);
            tc.verifyEqual(ramped.W0_mW, flat.W0_mW * compute_ramp_factor(0.02, 0.004, 'linear'), 'RelTol', 1e-12);
            tr.annular = rmfield(tr.annular, 'elem_p_amp');
            [~, fb_flat] = compute_tic(rmfield(p, 'timing'), tr, 20);
            [~, fb_ramp] = compute_tic(p, tr, 20);
            tc.verifyEqual(fb_ramp.W0_mW, fb_flat.W0_mW * 0.2, 'RelTol', 1e-12);   % only the duty cycle
        end

        function test_compute_tic_fallback_without_amplitude(tc)
            tr = local_transducer();
            tr.annular = rmfield(tr.annular, 'elem_amp');
            p = struct('timing', struct('pd', 0.02, 'pri', 0.1));
            [v, info] = compute_tic(p, tr, 20);
            tc.verifyTrue(isfinite(v));
            tc.verifySubstring(info.note, 'fallback');
            Deq_cm = max(tr.annular.elem_od_mm) / 10;
            tc.verifyEqual(info.W0_mW, 20 * pi * (Deq_cm / 2)^2 * 0.2 * 1000, 'RelTol', 1e-12);
        end

        function test_skull_and_skin_mi_not_reported(tc)
            % kept in the CSV, dropped by the reader and absent from the verdict limits
            f = fullfile(tc.tmp, 'sub-001_layered.csv');
            writetable(table(1, 0.7, 0.9, 0.4, 0.6, 'VariableNames', ...
                {'subject_id', 'MI_tc', 'MI_skull', 'MI_skin', 'MI_brain'}), f);
            t = read_metrics_table(f);
            tc.verifyEmpty(intersect({'MI_skull', 'MI_skin'}, t.Properties.VariableNames));
            tc.verifyTrue(all(ismember({'MI_tc', 'MI_brain'}, t.Properties.VariableNames)));
            L = get_risk_limits(true);
            tc.verifyFalse(any(isfield(L, {'MI_skull', 'MI_skin'})));
            tc.verifyTrue(all(isfield(L, {'MI_tc', 'MI_brain'})));
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
                     'Cranial TI (TIC)', 'exit-plane-ref.', 'config-dump', 'theme-btn'}
                tc.verifyTrue(~isempty(strfind(txt, m{1})), sprintf('missing %s', m{1})); %#ok<STREMP>
            end
            % Structural split: an MI tile sits in the Safety section, TIC in Other.
            iSafety = strfind(txt, 'id="safety"'); iOther = strfind(txt, 'id="other"');
            iMI = strfind(txt, 'MI (transcranial)'); iTIC = strfind(txt, 'Cranial TI (TIC)');
            tc.verifyTrue(iMI(1) > iSafety(1) && iMI(1) < iOther(1), 'MI should be in the Safety section');
            tc.verifyTrue(iTIC(1) > iOther(1), 'TIC should be in the Other section');
        end

        function test_methods_text_is_manuscript_ready(tc)
            p = local_make_subject_params(tc.tmp);
            txt = fileread(generate_simulation_report(p));
            tc.verifyEqual(numel(strfind(txt, 'class="filter-btn methods-copy"')), 2, 'expected two copy buttons');
            iStart = strfind(txt, 'id="methods-text"');
            iRefs  = strfind(txt, 'id="methods-refs"');
            tc.verifyNotEmpty(iStart, 'methods text container missing');
            tc.verifyNotEmpty(iRefs, 'methods references block missing');
            methods = txt(iStart(1):iRefs(1));
            % Self-contained: no links back into the report.
            tc.verifyEmpty(strfind(methods, 'href="#'), 'methods text should not link to report sections');
            % Values stated inline (tissue properties, simulated exposure, software).
            for m = {'Simulation software.', 'c&#8201;=&#8201;1560&#8201;m/s', ...
                     'k&#8201;=&#8201;0.51', 'Simulated exposure.', '21.4&#8201;W/cm', ...
                     'MI<sub>tc</sub> was 0.74', 'CEM43, ISO form', '0.9&#8201;min (brain)'}
                tc.verifyNotEmpty(strfind(methods, m{1}), sprintf('methods text missing %s', m{1}));
            end
            tc.verifyNotEmpty(strfind(txt, '10.5281/zenodo.15095860'), 'PRESTUS reference missing');
            tc.verifyNotEmpty(strfind(txt, '10.1016/j.brs.2024.04.013'), 'Martin et al. 2024 DOI wrong');
            % ISO CEM43 form (cem43_iso = 1) cites its sources and lists ISO 14708-3.
            tc.verifyNotEmpty(strfind(methods, 'ISO&#8201;14708-3:2017, Formula&#8201;1'), 'ISO CEM43 form must cite ISO 14708-3');
            tc.verifyNotEmpty(strfind(methods, 'Aubry et al., 2025; Sapareto'), 'ISO CEM43 form must cite Aubry et al. 2025');
            tc.verifyNotEmpty(strfind(txt, 'iso.org/standard/60539'), 'ISO 14708-3 reference missing');
        end

        function test_methods_text_babelbrain_matches_bridge(tc)
            % BabelBrain runs fluid bone with PRESTUS properties injected
            % (run_babelbrain_case.py), so the text must not claim shear waves
            % and must state the PRESTUS tissue values.
            p = local_make_subject_params(tc.tmp);
            p.acoustic_backend = 'babelbrain';
            p.babelbrain = struct('ppw', 12, 'tx_system', 'CTX_500', 'computing_backend', 'metal', ...
                'ct_mapping', '', 'q_correction', [], 'water_amp_calibration', true);
            p.transducer.annular.elem_phase_rad = [0 0.4 0.9 1.5];
            txt = fileread(generate_simulation_report(p));
            iStart = strfind(txt, 'id="methods-text"');
            iRefs  = strfind(txt, 'id="methods-refs"');
            methods = txt(iStart(1):iRefs(1));
            tc.verifyEmpty(strfind(methods, 'mode conversion'), 'fluid-bone run must not claim mode conversion');
            tc.verifyEmpty(strfind(methods, 'viscoelastic'), 'fluid-bone run must not claim a viscoelastic model');
            for m = {'shear wave speed set to zero', 'c&#8201;=&#8201;1560&#8201;m/s', ...
                     'quality-factor correction 1)', 'Rayleigh&#8211;Sommerfeld', ...
                     'at least 12 points per wavelength', 'Metal (Apple Silicon)'}
                tc.verifyNotEmpty(strfind(methods, m{1}), sprintf('methods text missing %s', m{1}));
            end
            tc.verifyEmpty(strfind(txt, 'Goss, S. A.'), 'uncited BabelBrain table reference listed');
        end

        function test_babelbrain_thermal_reports_kwave_cem43(tc)
            % cem43_iso = 1 in the fixture; BabelBrain thermal has no ISO maps,
            % so the report must describe and cite BabelBrain's own CEM43.
            p = local_make_subject_params(tc.tmp);
            p.acoustic_backend = 'babelbrain';
            p.thermal_backend  = 'babelbrain';
            p.babelbrain = struct('ppw', 12, 'tx_system', 'CTX_500', 'computing_backend', 'metal', ...
                'ct_mapping', '', 'q_correction', [], 'water_amp_calibration', true);
            txt = fileread(generate_simulation_report(p));
            iStart = strfind(txt, 'id="methods-text"');
            iRefs  = strfind(txt, 'id="methods-refs"');
            methods = txt(iStart(1):iRefs(1));
            tc.verifyEmpty(strfind(methods, 'ISO'), 'BabelBrain thermal run must not claim the ISO CEM43 form');
            tc.verifyNotEmpty(strfind(methods, 'R&#8201;=&#8201;0.25 for T&#8201;&lt;&#8201;43'), 'BabelBrain CEM43 formula missing');
            tc.verifyEmpty(strfind(txt, 'iso.org/standard/60539'), 'ISO 14708-3 listed without the ISO form');
        end

        function test_methods_text_jwave_states_resolved_pml(tc)
            p = local_make_subject_params(tc.tmp);
            p.acoustic_backend = 'jwave';
            p.jwave = struct('solver', 'helmholtz', 'pml_size', 'auto');
            p.grid.pml_size = 'auto';
            p.grid.pml_size_effective = 24;
            txt = fileread(generate_simulation_report(p));
            tc.verifyNotEmpty(strfind(txt, 'perfectly matched layer thickness: 24 grid points'), ...
                'jwave methods text should state the resolved PML');
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

        function test_read_metrics_table_renames_apex_column(tc)
            % Earlier CSVs name the apex-referenced peak distance
            % real_focal_distance_mm; readers see real_focal_distance_bowl_mm
            f = fullfile(tc.tmp, 'sub-001_layered.csv');
            writetable(table(1, 61.5, 50.9, 'VariableNames', ...
                {'subject_id', 'real_focal_distance_mm', 'real_focal_distance_ep_mm'}), f);
            t = read_metrics_table(f);
            tc.verifyFalse(ismember('real_focal_distance_mm', t.Properties.VariableNames));
            tc.verifyEqual(t.real_focal_distance_bowl_mm, 61.5);
            tc.verifyEqual(t.real_focal_distance_ep_mm, 50.9);
        end

    end

    % html_utils escaping
    methods (Test, TestTags = {'reports'})

        function test_script_json_cannot_close_script_block(tc)
            j = html_utils.script_json(struct('label', '</script><script>alert(1)</script>'));
            tc.verifyEmpty(strfind(j, '</'));
            tc.verifyEqual(jsondecode(j).label, '</script><script>alert(1)</script>');
            tc.verifyEqual(html_utils.escape('a''b'), 'a&#39;b');
        end

    end

    % risk colours and style tokens
    methods (Test, TestTags = {'reports'})

        function test_risk_color_all_branches(tc)
            % [] errored before: isnan([]) gives an empty operand to ||
            tc.verifyEqual(risk_color([], 2), 'gray');
            tc.verifyEqual(risk_color([], Inf), 'gray');
            tc.verifyEqual(risk_color(NaN, 2), 'gray');
            tc.verifyEqual(risk_color(1, Inf), 'info');
            tc.verifyEqual(risk_color(3, 2), 'red');
            tc.verifyEqual(risk_color(1, 2), 'amber');
            tc.verifyEqual(risk_color(2, 2), 'amber');      % at the limit
            tc.verifyEqual(risk_color(0.5, 2), 'green');
        end

        function test_css_tokens_are_defined(tc)
            % every var(--name) in the report and group CSS resolves to a
            % token defined there (--border and --mono did not exist)
            root = fileparts(fileparts(which('risk_color')));
            files = [dir(fullfile(root, 'report', '*.m')); dir(fullfile(root, 'group', '*.m'))];
            src = cell(numel(files), 1);
            for k = 1:numel(files)
                src{k} = fileread(fullfile(files(k).folder, files(k).name));
            end
            t = regexp(strjoin(src, newline), '--([A-Za-z0-9-]+)\s*:', 'tokens');
            defined = unique([t{:}]);
            tc.assertNotEmpty(defined);
            for k = 1:numel(files)
                t = regexp(src{k}, 'var\(--([A-Za-z0-9-]+)', 'tokens');
                if isempty(t); continue; end
                used = unique([t{:}]);
                tc.verifyEmpty(setdiff(used, defined), ...
                    sprintf('%s uses undefined CSS tokens', files(k).name));
            end
            % negative control: an undefined token is caught
            tc.verifyNotEmpty(setdiff({'mono'}, defined));
        end

        function test_css_has_no_left_accent_bar(tc)
            % banners, callouts, flags, and cards are plain tinted boxes in
            % the base CSS and in every skin (the base had 4px coloured
            % left borders, a left inset shadow, and ::before bars)
            sheets = {'css_styles_base', css_styles_base()};
            skins = report_styles.names();
            for k = 1:numel(skins)
                sheets(end+1, :) = {skins{k}, ...
                    report_styles.css(struct('name', skins{k}, 'switcher', false))}; %#ok<AGROW>
            end
            sheets(end+1, :) = {'switcher', ...
                report_styles.css(struct('name', 'instrument', 'switcher', true))};
            for k = 1:size(sheets, 1)
                bars = local_left_bars(sheets{k, 2});
                tc.verifyEmpty(bars, sprintf('left accent bar in %s CSS: %s', ...
                    sheets{k, 1}, strjoin(bars, ' | ')));
            end
            % negative controls: each form of left bar is caught
            bad = {'.a{border-left:4px solid var(--red)}', ...
                   '.a{border-inline-start:2px solid var(--line)}', ...
                   '.a{border-left:1px solid var(--info)}', ...
                   '.a{border:0;border-left-color:var(--green)}', ...
                   '.a{border-left-width:3px}', ...
                   'tr td{box-shadow:inset 3px 0 0 var(--red)}', ...
                   '.a::before{content:"";position:absolute;left:0;top:0;bottom:0;width:6px}'};
            for k = 1:numel(bad)
                tc.verifyNotEmpty(local_left_bars(bad{k}), sprintf('not caught: %s', bad{k}));
            end
            % allowed: thin neutral rule, full ring, full tinted border
            ok = {'.a{border-left:1px solid var(--line)}', ...
                  '.a{box-shadow:0 0 0 2px var(--red) inset}', ...
                  '.a{border:1px solid color-mix(in srgb,var(--red) 30%,transparent)}'};
            for k = 1:numel(ok)
                tc.verifyEmpty(local_left_bars(ok{k}), sprintf('false positive: %s', ok{k}));
            end
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
    row.real_focal_distance_mm = 64.2; row.real_focal_distance_ep_mm = 12.8;
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

function bars = local_left_bars(css)
% Rules in css that draw a left accent bar: a border-left/border-inline-start
% wider than 1px or coloured other than var(--line), a left inset box-shadow,
% or an absolutely positioned narrow ::before/::after pinned to the left edge.
    bars = {};
    rules = regexp(css, '([^{}]+)\{([^{}]*)\}', 'tokens');
    for r = 1:numel(rules)
        sel = strtrim(rules{r}{1});
        body = rules{r}{2};
        decl = regexp(body, '(border-(?:left|inline-start)(?:-color|-width|-style)?)\s*:\s*([^;]*)', 'tokens');
        for d = 1:numel(decl)
            prop = decl{d}{1};
            val = strtrim(regexprep(decl{d}{2}, '\s*!important', ''));
            if endsWith(prop, '-style'); continue; end
            w = str2double(regexp(val, '(\d*\.?\d+)px', 'tokens', 'once'));
            colour = strtrim(regexprep(val, ...
                '(\d*\.?\d+px|\<(solid|dashed|dotted|double|none|hidden|thin|medium|thick|0)\>)', ''));
            if endsWith(prop, '-width')
                colour = '';
            elseif endsWith(prop, '-color')
                w = NaN;
                colour = val;
            end
            if any(strcmp(val, {'0', 'none'})); continue; end
            if (~isnan(w) && w > 1) || any(strcmp(val, {'medium', 'thick'})) ...
                    || (~isempty(colour) && ~strcmp(colour, 'var(--line)'))
                bars{end+1} = sprintf('%s{%s:%s}', sel, prop, val); %#ok<AGROW>
            end
        end
        sh = regexp(body, 'box-shadow\s*:\s*([^;]*)', 'tokens');
        for s = 1:numel(sh)
            parts = strsplit(sh{s}{1}, ',');
            for p = 1:numel(parts)
                if ~contains(parts{p}, 'inset'); continue; end
                xy = regexp(strtrim(strrep(parts{p}, 'inset', '')), ...
                    '^(-?\d*\.?\d+)(?:px)?\s+(-?\d*\.?\d+)(?:px)?\s+(-?\d*\.?\d+)(?:px)?', 'tokens', 'once');
                if ~isempty(xy) && str2double(xy{1}) > 0 && str2double(xy{2}) == 0
                    bars{end+1} = sprintf('%s{box-shadow:%s}', sel, sh{s}{1}); %#ok<AGROW>
                end
            end
        end
        if ~isempty(regexp(sel, '::?(before|after)', 'once')) ...
                && ~isempty(regexp(body, 'position\s*:\s*absolute', 'once')) ...
                && ~isempty(regexp(body, '(^|;)\s*left\s*:\s*0', 'once')) ...
                && ~isempty(regexp(body, '(^|;)\s*top\s*:\s*0', 'once')) ...
                && ~isempty(regexp(body, '(^|;)\s*bottom\s*:\s*0', 'once'))
            w = str2double(regexp(body, '(?:^|;)\s*width\s*:\s*(\d*\.?\d+)px', 'tokens', 'once'));
            if ~isnan(w) && w <= 10
                bars{end+1} = sprintf('%s{%s}', sel, body); %#ok<AGROW>
            end
        end
    end
end
