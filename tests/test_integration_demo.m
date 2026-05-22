classdef test_integration_demo < matlab.unittest.TestCase
% TEST_INTEGRATION_DEMO  Integration tests using real demo data (sub-009).
%
%   Uses a standardised test-data root outside the PRESTUS repository.
%   All tests skip gracefully when the root folder is absent.
%
%   Tags and expected runtime:
%     'demo_inputs'         — verify every expected input file exists    (~1 s)
%     'demo_config'         — load each placement-mode config            (~1 s)
%     'demo_localite'       — parse Localite XML → finite voxel coords   (~1 s)
%     'demo_pct'            — verify pre-computed pCT data is present    (~1 s)
%     'demo_head'           — head preprocessing end-to-end              (~1–5 min)
%     'demo_acoustic'       — full acoustic pipeline                     (~10–60 min)
%     'demo_thermal'        — thermal pipeline on cached acoustics
%
%   Set the data root before running:
%     setenv('PRESTUS_DEMO_DATA', '/path/to/prestus_testdata')
%
%   Run a single tag (from repo root):
%     runtests('tests/test_integration_demo.m', 'Tag', 'demo_inputs')
%
%   Run all tags up to and including head preprocessing:
%     run_all_tests('demo_head')

    properties (Constant)
        SUBJECT_ID = 9
    end

    properties
        data_root    % prestus_testdata root (from env var)
        repo_root    % PRESTUS repo root
    end

    % ------------------------------------------------------------------ %
    %% Class-level setup
    % ------------------------------------------------------------------ %
    methods (TestClassSetup)
        function resolve_paths(tc)
            test_dir     = fileparts(mfilename('fullpath'));
            tc.repo_root = fileparts(test_dir);
            addpath(genpath(fullfile(tc.repo_root, 'functions')));
            addpath(fullfile(tc.repo_root, 'config'));
            addpath(genpath(fullfile(tc.repo_root, 'external')));
        end
    end

    methods (TestMethodSetup)
        function require_demo_data(tc)
            root = getenv('PRESTUS_DEMO_DATA');
            tc.assumeFalse(isempty(root), ...
                'PRESTUS_DEMO_DATA not set. Skipping demo tests.');
            tc.assumeTrue(isfolder(root), ...
                sprintf('PRESTUS_DEMO_DATA folder not found: %s', root));
            tc.data_root = root;
        end
    end

    % ------------------------------------------------------------------ %
    %% Helpers
    % ------------------------------------------------------------------ %
    methods (Access = private)

        function p = load_demo_config(tc, config_name)
            % Load one of the three placement-mode demo configs and inject
            % the standard demo data paths plus subject ID.
            cfg_file = fullfile(tc.data_root, 'configs', config_name);
            tc.assumeTrue(isfile(cfg_file), ...
                sprintf('Demo config not found: %s', cfg_file));

            cfg_struct = yaml.loadFile(cfg_file, 'ConvertToArray', true);
            cfg_struct.simulation.debug = 0;
            p = load_parameters(cfg_struct);

            p.subject_id             = tc.SUBJECT_ID;
            p.simulation.interactive = 0;
            p.platform               = 'matlab';
            p.io.overwrite_files     = 'always';

            p.path.anat = fullfile(tc.data_root, 'bids', ...
                sprintf('sub-%03d', tc.SUBJECT_ID), 'anat');
            p.path.seg  = fullfile(tc.data_root, 'simnibs');
            p.path.sim  = fullfile(tc.data_root, 'sim_outputs');
        end

        function p = load_localite_config(tc)
            p = tc.load_demo_config('config_demo_localite.yaml');
            p.path.localite = fullfile(tc.data_root, 'localite');
        end

        function p = load_heuristic_config(tc)
            p = tc.load_demo_config('config_demo_heuristic.yaml');
            coord_file = fullfile(tc.data_root, 'coords', ...
                sprintf('sub-%03d_heuristic_target.json', tc.SUBJECT_ID));
            tc.assumeTrue(isfile(coord_file), ...
                sprintf('Heuristic coord file not found: %s', coord_file));
            coords = jsondecode(fileread(coord_file));
            p.placement.heuristic.mni_target_mm = coords.mni_target_mm(:)';
            p.placement.heuristic.target_name   = coords.target_name;
        end

    end

    % ------------------------------------------------------------------ %
    %% demo_inputs — file presence checks (no MATLAB processing)
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'demo_inputs'})

        function test_t1_exists(tc)
            f = fullfile(tc.data_root, 'bids', ...
                sprintf('sub-%03d', tc.SUBJECT_ID), 'anat', ...
                sprintf('sub-%03d_T1w.nii.gz', tc.SUBJECT_ID));
            tc.verifyTrue(isfile(f), ['T1 not found: ' f]);
        end

        function test_ute_exists(tc)
            f = fullfile(tc.data_root, 'bids', ...
                sprintf('sub-%03d', tc.SUBJECT_ID), 'anat', ...
                sprintf('sub-%03d_UTE.nii.gz', tc.SUBJECT_ID));
            tc.verifyTrue(isfile(f), ['UTE not found: ' f]);
        end

        function test_segmentation_exists(tc)
            f = fullfile(tc.data_root, 'simnibs', ...
                sprintf('m2m_sub-%03d', tc.SUBJECT_ID), 'final_tissues.nii.gz');
            tc.verifyTrue(isfile(f), ['Segmentation not found: ' f]);
        end

        function test_localite_xml_exists(tc)
            d = fullfile(tc.data_root, 'localite', ...
                sprintf('sub-%03d', tc.SUBJECT_ID), 'ses-02', 'localite');
            hits = dir(fullfile(d, '**', 'TriggerMarkers_Coil0*.xml'));
            tc.verifyGreaterThan(numel(hits), 0, ...
                ['No TriggerMarkers XML found under: ' d]);
        end

        function test_heuristic_coord_file_exists(tc)
            f = fullfile(tc.data_root, 'coords', ...
                sprintf('sub-%03d_heuristic_target.json', tc.SUBJECT_ID));
            tc.verifyTrue(isfile(f), ['Coord file not found: ' f]);
        end

        function test_heuristic_coord_file_valid(tc)
            f = fullfile(tc.data_root, 'coords', ...
                sprintf('sub-%03d_heuristic_target.json', tc.SUBJECT_ID));
            tc.assumeTrue(isfile(f));
            coords = jsondecode(fileread(f));
            tc.verifyTrue(isfield(coords, 'mni_target_mm'), ...
                'coord JSON missing field mni_target_mm');
            tc.verifyTrue(isfield(coords, 'target_name'), ...
                'coord JSON missing field target_name');
            tc.verifyEqual(numel(coords.mni_target_mm), 3, ...
                'mni_target_mm must have exactly 3 elements');
            tc.verifyTrue(ischar(coords.target_name) || isstring(coords.target_name), ...
                'target_name must be a string');
        end

        function test_pct_pseudoCT_exists(tc)
            f = fullfile(tc.data_root, 'simnibs', ...
                sprintf('m2m_sub-%03d', tc.SUBJECT_ID), 'pseudoCT.nii.gz');
            tc.verifyTrue(isfile(f), ['pseudoCT.nii.gz not found: ' f]);
        end

        function test_pct_ute_reg_exists(tc)
            f = fullfile(tc.data_root, 'simnibs', ...
                sprintf('m2m_sub-%03d', tc.SUBJECT_ID), 'UTE_reg.nii.gz');
            tc.verifyTrue(isfile(f), ['UTE_reg.nii.gz not found: ' f]);
        end

        function test_mni_transforms_exist(tc)
            d = fullfile(tc.data_root, 'simnibs', ...
                sprintf('m2m_sub-%03d', tc.SUBJECT_ID), 'toMNI');
            tc.verifyTrue(isfolder(d), ['toMNI folder not found: ' d]);
            tc.verifyTrue(isfile(fullfile(d, 'Conform2MNI_nonl.nii.gz')), ...
                'Conform2MNI_nonl.nii.gz not found in toMNI/');
        end

    end

    % ------------------------------------------------------------------ %
    %% demo_config — config loading for all three placement modes
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'demo_config'})

        function test_manual_config_loads(tc)
            p = tc.load_demo_config('config_demo_manual.yaml');
            tc.verifyClass(p, 'struct');
            tc.verifyEqual(p.placement.mode, 'manual');
            tc.verifyEqual(numel(p.transducer.annular.elem_amp), 10);
        end

        function test_localite_config_loads(tc)
            p = tc.load_localite_config();
            tc.verifyClass(p, 'struct');
            tc.verifyEqual(p.placement.mode, 'localite');
            tc.verifyTrue(p.placement.localite.enabled == 1);
        end

        function test_heuristic_config_loads(tc)
            p = tc.load_heuristic_config();
            tc.verifyClass(p, 'struct');
            tc.verifyEqual(p.placement.mode, 'heuristic');
            tc.verifyEqual(numel(p.placement.heuristic.mni_target_mm), 3);
            tc.verifyNotEmpty(p.placement.heuristic.target_name);
        end

        function test_t1_pattern_resolves(tc)
            p = tc.load_demo_config('config_demo_manual.yaml');
            t1_path = fullfile(p.path.anat, ...
                sprintf(p.path.t1_pattern, tc.SUBJECT_ID));
            % Strip glob wildcard for isfile check
            t1_path = strrep(t1_path, '*', '');
            hits = dir([t1_path '*']);
            tc.verifyGreaterThan(numel(hits), 0, ...
                ['T1 pattern does not resolve to an existing file: ' t1_path]);
        end

    end

    % ------------------------------------------------------------------ %
    %% demo_localite — Localite XML parsing → finite voxel positions
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'demo_localite'})

        function test_localite_xml_parses(tc)
            p = tc.load_localite_config();
            % Locate the TriggerMarkers XML
            xml_dir = fullfile(p.path.localite, ...
                sprintf('sub-%03d', tc.SUBJECT_ID), 'ses-02', 'localite');
            hits = dir(fullfile(xml_dir, '**', 'TriggerMarkers_Coil0*.xml'));
            tc.assumeGreaterThan(numel(hits), 0, ...
                'No TriggerMarkers XML to parse.');
            xml_file = fullfile(hits(1).folder, hits(1).name);

            t1_file = fullfile(p.path.anat, ...
                sprintf('sub-%03d_T1w.nii.gz', tc.SUBJECT_ID));
            tc.assumeTrue(isfile(t1_file));
            t1_hdr = niftiinfo(t1_file);

            [trans_pos, focus_pos, trans_ras, focus_ras] = ...
                position_transducer_localite(xml_file, t1_hdr, p);

            tc.verifyEqual(numel(trans_pos), 3, ...
                'trans_pos should have 3 elements');
            tc.verifyEqual(numel(focus_pos), 3, ...
                'focus_pos should have 3 elements');
            tc.verifyTrue(all(isfinite(trans_ras)), ...
                'trans_pos_ras contains non-finite values');
            tc.verifyTrue(all(isfinite(focus_ras)), ...
                'focus_pos_ras contains non-finite values');
            tc.verifyTrue(all(trans_pos > 0), ...
                'trans_pos voxel indices should be positive');
        end

        function test_neuronav_select_finds_session(tc)
            p = tc.load_localite_config();
            pn.data_postlocalite = p.path.localite;
            sub_id = sprintf('sub-%03d', tc.SUBJECT_ID);
            ses_id = p.placement.localite.session;

            localite = neuronav_select_localite(pn, sub_id, ses_id, 'TriggerMarkers');
            tc.verifyFalse(isempty(localite), ...
                'neuronav_select_localite returned empty — session not found.');
        end

    end

    % ------------------------------------------------------------------ %
    %% demo_pct — verify pCT inputs and pre-computed outputs
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'demo_pct'})

        function test_pct_enabled_config_loads(tc)
            p = tc.load_demo_config('config_demo_manual.yaml');
            p.pct.enabled = 1;
            tc.verifyTrue(p.pct.enabled == 1);
            tc.verifyTrue(isfield(p.pct, 'skull_mapping'));
        end

        function test_pseudoCT_nifti_readable(tc)
            f = fullfile(tc.data_root, 'simnibs', ...
                sprintf('m2m_sub-%03d', tc.SUBJECT_ID), 'pseudoCT.nii.gz');
            tc.assumeTrue(isfile(f));
            hdr = niftiinfo(f);
            tc.verifyEqual(numel(hdr.ImageSize), 3, ...
                'pseudoCT.nii.gz should be a 3-D volume');
            tc.verifyTrue(all(hdr.ImageSize > 0));
        end

        function test_ute_reg_nifti_readable(tc)
            f = fullfile(tc.data_root, 'simnibs', ...
                sprintf('m2m_sub-%03d', tc.SUBJECT_ID), 'UTE_reg.nii.gz');
            tc.assumeTrue(isfile(f));
            hdr = niftiinfo(f);
            tc.verifyEqual(numel(hdr.ImageSize), 3, ...
                'UTE_reg.nii.gz should be a 3-D volume');
        end

    end

    % ------------------------------------------------------------------ %
    %% demo_head — head preprocessing (segmentation + medium masks)
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'demo_head'})

        function test_preproc_head_manual(tc)
            p = tc.load_demo_config('config_demo_manual.yaml');
            p = tc.disable_sims(p);
            [p_out, medium_masks, segmentation, ~, ~] = preproc_head(p);
            tc.verifyNotEmpty(medium_masks);
            tc.verifyNotEmpty(segmentation);
            tc.verifyTrue(isfield(p_out, 'grid'));
            tc.verifyTrue(isfield(p_out.grid, 'dims'));
        end

        function test_medium_masks_label_range_manual(tc)
            p = tc.load_demo_config('config_demo_manual.yaml');
            p = tc.disable_sims(p);
            [~, medium_masks, ~, ~, ~] = preproc_head(p);
            n_layers = numel(fieldnames(p.medium_properties));
            tc.verifyGreaterThanOrEqual(min(medium_masks(:)), 0);
            tc.verifyLessThanOrEqual(max(medium_masks(:)), n_layers);
        end

        function test_preproc_head_pct(tc)
            p = tc.load_demo_config('config_demo_manual.yaml');
            p = tc.disable_sims(p);
            p.pct.enabled = 1;
            [p_out, medium_masks, ~, pseudoCT_crop, ~] = preproc_head(p);
            tc.verifyNotEmpty(medium_masks);
            tc.verifyNotEmpty(pseudoCT_crop, ...
                'pseudoCT crop should be non-empty when pct.enabled=1');
            tc.verifyTrue(isfield(p_out, 'grid'));
        end

    end

    % ------------------------------------------------------------------ %
    %% demo_acoustic — full acoustic pipeline
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'demo_acoustic'})

        function test_acoustic_output_created_manual(tc)
            p = tc.load_demo_config('config_demo_manual.yaml');
            p.modules.run_heating_sims       = 0;
            p.modules.run_thermal_analysis   = 0;
            p.modules.run_posthoc_water_sims = 0;
            p.modules.generate_report        = 0;

            prestus_pipeline(p);

            pat = fullfile(p.path.sim, ...
                sprintf('sub-%03d', tc.SUBJECT_ID), ...
                sprintf('sub-%03d_%s%s_intensity_orig_coord.nii.gz', ...
                    tc.SUBJECT_ID, p.simulation.medium, p.io.output_affix));
            tc.verifyTrue(isfile(pat), ...
                ['Acoustic output NIfTI not created: ' pat]);
        end

    end

    % ------------------------------------------------------------------ %
    %% demo_thermal — thermal pipeline (requires cached acoustic output)
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'demo_thermal'})

        function test_thermal_output_created_manual(tc)
            p = tc.load_demo_config('config_demo_manual.yaml');
            p.modules.run_acoustic_sims      = 0;  % reuse cached output
            p.modules.run_heating_sims       = 1;
            p.modules.run_posthoc_water_sims = 0;
            p.modules.generate_report        = 0;
            p.io.overwrite_files             = 'never';

            prestus_pipeline(p);

            pat = fullfile(p.path.sim, ...
                sprintf('sub-%03d', tc.SUBJECT_ID), ...
                sprintf('sub-%03d_%s%s_cem43_orig_coord.nii.gz', ...
                    tc.SUBJECT_ID, p.simulation.medium, p.io.output_affix));
            tc.verifyTrue(isfile(pat), ...
                ['Thermal CEM43 NIfTI not created: ' pat]);
        end

    end

    % ------------------------------------------------------------------ %
    %% Private helpers
    % ------------------------------------------------------------------ %
    methods (Access = private, Static)
        function p = disable_sims(p)
            p.modules.run_medium_setup       = 1;
            p.modules.run_source_setup       = 0;
            p.modules.run_acoustic_sims      = 0;
            p.modules.run_acoustic_analysis  = 0;
            p.modules.run_heating_sims       = 0;
            p.modules.run_thermal_analysis   = 0;
            p.modules.run_nifti_creation     = 0;
            p.modules.run_posthoc_water_sims = 0;
            p.modules.generate_report        = 0;
        end
    end

end
