function make_bridge_fixtures()
% MAKE_BRIDGE_FIXTURES  Record nexus_calibration / nexus_neuronav answers as JSON fixtures
%
% The Nexus page tests (test_bridge_pages.py) replay these instead of running MATLAB, so
% they parse exactly what NEXUS_BRIDGE sends. Rerun after changing the façades:
%   matlab -batch "cd('<PRESTUS>/nexus/tests'); make_bridge_fixtures"
% Needs local/profiles/CTX_250_001 (hydrophone profiles) for the calibration answers.

    here = fileparts(mfilename('fullpath'));
    root = fileparts(fileparts(here));
    addpath(genpath(fullfile(root, 'functions')));
    addpath(genpath(fullfile(root, 'external', 'yaml')));
    addpath(fullfile(root, 'tests', 'fixtures'));
    out = fullfile(here, 'fixtures', 'bridge');
    if ~isfolder(out); mkdir(out); end

    % Calibration (CTX_250_001 / SC_105_010, shipped entry)
    st = struct('forward_model', 'rayleigh', 'near_field_skip_mm', 10, 'huber_delta', 0.8, ...
        'n_starts', 4, 'domain', 'intensity', 'isppa_wcm2', 30, 'align', 'label', ...
        'max_extrapolation_mm', 0, 'advanced', struct());
    a = struct('session', 'fx', 'settings', st);
    save_json(out, 'calibration_open', nexus_calibration(struct('action', 'open', 'session', 'fx', ...
        'serial', 'CTX_250_001', 'driving_system', 'SC_105_010', 'library_path', '')));
    save_json(out, 'calibration_view_measured', nexus_calibration(setf(a, 'action', 'view', 'predicted', false)));
    save_json(out, 'calibration_view_predicted', nexus_calibration(setf(a, 'action', 'view', ...
        'predicted', true, 'depths', [30; 45; 70])));
    save_json(out, 'calibration_fit_per_depth', nexus_calibration(setf(a, 'action', 'fit', ...
        'mode', 'per_depth', 'free_params', {'phase'; 'amp'}, 'calibrate_at', 'measured', ...
        'measured_depths', [31.11; 45.45])));
    save_json(out, 'calibration_fit_static', nexus_calibration(setf(a, 'action', 'fit', ...
        'mode', 'static', 'free_params', {'phase'; 'amp'}, 'calibrate_at', 'measured', ...
        'measured_depths', [31.11; 45.45; 54.09])));

    % Localite (synthetic study, as TEST_NEURONAV_BATCH)
    R = tempname; mkdir(R); c = onCleanup(@() rmdir(R, 's'));
    rz = @(deg) [cosd(deg) -sind(deg) 0; sind(deg) cosd(deg) 0; 0 0 1];
    poses = {[rz(20), [40; 90; 140]; 0 0 0 1], [rz(160), [140; 95; 135]; 0 0 0 1]};
    make_synthetic_localite_export(fullfile(R, 'raw'), 'sub-001', 'ses-01', poses, [1; 1; 0], {'LEFT PUL', 'RIGHT PUL'});
    cfg = fullfile(R, 'config_study.yaml');
    fid = fopen(cfg, 'w');
    fprintf(fid, "path:\n  localite_raw: '%s'\n  localite_post: '%s'\n", fullfile(R, 'raw'), fullfile(R, 'post'));
    fclose(fid);
    save_json(out, 'neuronav_detect', nexus_neuronav(struct('action', 'detect', 'config', cfg, ...
        'sub', 'sub-001', 'ses', 'ses-01')));
    tm = struct('series_index', {1, 2}, 'transducer_id', 0, 'target_name', {'L_PUL', 'R_PUL'}, 'color', '');
    save_json(out, 'neuronav_write', nexus_neuronav(struct('action', 'write', 'config', cfg, ...
        'sub', 'sub-001', 'ses', 'ses-01', 'target_map', tm)));
    save_json(out, 'neuronav_scan', nexus_neuronav(struct('action', 'scan', 'config', cfg)));
end

function s = setf(s, varargin)
    for k = 1:2:numel(varargin); s.(varargin{k}) = varargin{k + 1}; end
end

function save_json(folder, name, value)
% Machine paths are replaced by <PRESTUS> and <tmp>, so the fixtures can be committed
    txt = jsonencode(value, 'ConvertInfAndNaN', true);
    txt = strrep(txt, get_prestus_path(), '<PRESTUS>');
    txt = regexprep(txt, [regexptranslate('escape', strip(tempdir, 'right', filesep)) '/[^/"]+'], '<tmp>');
    fid = fopen(fullfile(folder, [name '.json']), 'w', 'n', 'UTF-8');
    fwrite(fid, txt, 'char');
    fclose(fid);
    fprintf('wrote %s.json\n', name);
end
