function setup_demo_data(demo_root)
% SETUP_DEMO_DATA  Validate the prestus_testdata directory for demo/integration tests.
%
%   setup_demo_data(demo_root)
%
%   Checks that every file and folder required by test_integration_demo is
%   present under DEMO_ROOT, then prints a checklist.  No files are created
%   or modified — this is a read-only diagnostic.
%
%   After running, set the environment variable so the test suite can find
%   the data:
%
%     setenv('PRESTUS_DEMO_DATA', demo_root)
%     run_all_tests('demo_inputs')
%
%   Expected layout of demo_root (sub-009 as demo subject):
%
%     prestus_testdata/
%     ├── bids/sub-009/anat/
%     │   ├── sub-009_T1w.nii.gz        (symlink → simnibs/m2m_sub-009/T1.nii.gz)
%     │   └── sub-009_UTE.nii.gz        (symlink → simnibs/m2m_sub-009/UTE_reg.nii.gz)
%     ├── simnibs/
%     │   └── m2m_sub-009/              (symlink or full copy of SimNIBS output)
%     │       ├── final_tissues.nii.gz
%     │       ├── T1.nii.gz
%     │       ├── pseudoCT.nii.gz
%     │       ├── UTE_reg.nii.gz
%     │       └── toMNI/Conform2MNI_nonl.nii.gz
%     ├── localite/
%     │   └── sub-009/ses-02/localite/
%     │       └── Session_*/TMSTrigger/TriggerMarkers_Coil0*.xml
%     ├── coords/
%     │   └── sub-009_heuristic_target.json
%     ├── configs/
%     │   ├── config_demo_manual.yaml
%     │   ├── config_demo_localite.yaml
%     │   └── config_demo_heuristic.yaml
%     └── sim_outputs/                  (written by the pipeline; may be empty)

    if nargin < 1
        demo_root = getenv('PRESTUS_DEMO_DATA');
    end

    if isempty(demo_root)
        error('setup_demo_data:noRoot', ...
            'Provide demo_root or set the PRESTUS_DEMO_DATA environment variable.');
    end

    sub_id  = 9;
    sub_str = sprintf('sub-%03d', sub_id);
    ses_str = 'ses-02';

    checks = build_checklist(demo_root, sub_str, ses_str);

    fprintf('\n====================================================\n');
    fprintf('  PRESTUS demo-data validation\n');
    fprintf('  Root: %s\n', demo_root);
    fprintf('====================================================\n');

    n_pass = 0;  n_fail = 0;
    for i = 1:numel(checks)
        c = checks(i);
        if c.is_folder
            ok = isfolder(c.path);
        else
            ok = isfile(c.path);
        end

        if ok
            status = 'OK  ';
            n_pass = n_pass + 1;
        else
            status = 'MISS';
            n_fail = n_fail + 1;
        end
        fprintf('  [%s]  %s\n         %s\n', status, c.label, c.path);
    end

    % Extra check: coord JSON content
    coord_file = fullfile(demo_root, 'coords', ...
        sprintf('%s_heuristic_target.json', sub_str));
    if isfile(coord_file)
        try
            coords = jsondecode(fileread(coord_file));
            has_mm   = isfield(coords, 'mni_target_mm') && numel(coords.mni_target_mm) == 3;
            has_name = isfield(coords, 'target_name');
            if has_mm && has_name
                fprintf('  [OK  ]  coord JSON has mni_target_mm[3] and target_name\n');
                fprintf('         target="%s"  mni=[%.1f %.1f %.1f]\n', ...
                    coords.target_name, coords.mni_target_mm(1), ...
                    coords.mni_target_mm(2), coords.mni_target_mm(3));
                n_pass = n_pass + 1;
            else
                fprintf('  [FAIL]  coord JSON missing required fields\n');
                n_fail = n_fail + 1;
            end
        catch ME
            fprintf('  [FAIL]  coord JSON could not be parsed: %s\n', ME.message);
            n_fail = n_fail + 1;
        end
    end

    fprintf('----------------------------------------------------\n');
    fprintf('  Present: %d   Missing: %d\n', n_pass, n_fail);
    fprintf('====================================================\n\n');

    if n_fail == 0
        fprintf('All checks passed. Set the environment variable and run:\n');
        fprintf('  setenv(''PRESTUS_DEMO_DATA'', ''%s'')\n', demo_root);
        fprintf('  run_all_tests(''demo_inputs'')\n\n');
    else
        fprintf('Fix missing items above, then re-run setup_demo_data.\n\n');
    end

end

% -------------------------------------------------------------------------
function checks = build_checklist(root, sub_str, ses_str)
    anat_dir    = fullfile(root, 'bids',    sub_str, 'anat');
    m2m_dir     = fullfile(root, 'simnibs', sprintf('m2m_%s', sub_str));
    loc_dir     = fullfile(root, 'localite', sub_str, ses_str, 'localite');
    trig_hits   = dir(fullfile(loc_dir, '**', 'TriggerMarkers_Coil0*.xml'));

    entries = {
        % label                                     path                                            is_folder
        'BIDS anat folder'                          anat_dir                                        true
        'T1w NIfTI'                                 fullfile(anat_dir, [sub_str '_T1w.nii.gz'])     false
        'UTE NIfTI'                                 fullfile(anat_dir, [sub_str '_UTE.nii.gz'])     false
        'SimNIBS m2m folder'                        m2m_dir                                         true
        'final_tissues.nii.gz'                      fullfile(m2m_dir, 'final_tissues.nii.gz')       false
        'T1.nii.gz (in m2m)'                        fullfile(m2m_dir, 'T1.nii.gz')                  false
        'pseudoCT.nii.gz'                           fullfile(m2m_dir, 'pseudoCT.nii.gz')            false
        'UTE_reg.nii.gz'                            fullfile(m2m_dir, 'UTE_reg.nii.gz')             false
        'toMNI/ folder'                             fullfile(m2m_dir, 'toMNI')                      true
        'Conform2MNI_nonl.nii.gz'                   fullfile(m2m_dir, 'toMNI', 'Conform2MNI_nonl.nii.gz') false
        'Localite session folder'                   loc_dir                                         true
        'heuristic coord JSON'                      fullfile(root, 'coords', [sub_str '_heuristic_target.json']) false
        'config_demo_manual.yaml'                   fullfile(root, 'configs', 'config_demo_manual.yaml')   false
        'config_demo_localite.yaml'                 fullfile(root, 'configs', 'config_demo_localite.yaml') false
        'config_demo_heuristic.yaml'                fullfile(root, 'configs', 'config_demo_heuristic.yaml') false
        'sim_outputs/ folder'                       fullfile(root, 'sim_outputs')                   true
    };

    % TriggerMarkers XML (found via glob above)
    if ~isempty(trig_hits)
        trig_path = fullfile(trig_hits(1).folder, trig_hits(1).name);
    else
        trig_path = fullfile(loc_dir, 'Session_*/TMSTrigger/TriggerMarkers_Coil0*.xml (not found)');
    end

    n = size(entries, 1) + 1;
    checks(n) = struct('label', '', 'path', '', 'is_folder', false);
    for i = 1:size(entries, 1)
        checks(i).label     = entries{i, 1};
        checks(i).path      = entries{i, 2};
        checks(i).is_folder = entries{i, 3};
    end
    checks(n).label     = 'TriggerMarkers XML (Coil0)';
    checks(n).path      = trig_path;
    checks(n).is_folder = false;
end
