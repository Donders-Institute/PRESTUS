function out_files = export_babelbrain_inputs( ...
        parameters, medium_masks, focus_pos, planimg, equip_param)
% EXPORT_BABELBRAIN_INPUTS  Write all BabelBrain input files from a PRESTUS run
%
% Prepares every file needed to run BabelBrain on the current subject and
% target, using PRESTUS-derived geometry, segmentation, and calibrated phases.
% Output goes to parameters.io.dir_babelbrain (created if absent).
%
% Files written:
%   <target>_<freq>kHz_<PPW>PPW_BabelViscoInput.nii.gz
%       Tissue mask with BabelBrain integer labels (0=background, 1=skin,
%       2=cortical, 3=trabecular, 4=brain, 5=focal point).
%
%   <target>_CalibWeights.h5
%       OptimizedWeightsFile with per-element phasors encoding the PRESTUS
%       calibration (delta from geometric phases). Stored as split real/imag
%       float32 datasets (MATLAB cannot write native complex HDF5).
%       run_babelbrain_case.py auto-converts to BabelBrain's native complex64
%       format (_bb.h5) before passing to RunCases.
%
%   RunBabelBrain.py
%       Ready-to-run Python script that calls RUN_SIM().RunCases() with the
%       exact parameters derived from this PRESTUS run. Edit the PYTHON_ENV
%       variable at the top before executing.
%
%   babelbrain_config.yaml
%       Human-readable summary of all BabelBrain parameters for record-keeping.
%
% Use as:
%   out_files = export_babelbrain_inputs(parameters, medium_masks, focus_pos, planimg)
%   out_files = export_babelbrain_inputs(..., equip_param)
%
% Input:
%   parameters    - PRESTUS config with .babelbrain sub-struct and .io.dir_babelbrain
%   medium_masks  - [Nx x Ny x Nz] medium label volume from preproc_medium_mask
%   focus_pos     - [1x3] focus voxel in grid coordinates
%   planimg       - planimg struct from grid_tissue_setup
%   equip_param   - (optional) equipment config for transducer lookup
%
% Output:
%   out_files - struct with fields:
%     .mask_nii     - path to BabelViscoInput.nii.gz
%     .weights_h5   - path to CalibWeights.h5 (empty if no calibrated phases)
%     .run_script   - path to RunBabelBrain.py
%     .config_yaml  - path to babelbrain_config.yaml
%
% See also: BABELBRAIN_WRAPPER, WRITE_BABELBRAIN_INPUT, EXPORT_BABELBRAIN_WEIGHTS,
%           BABELBRAIN_TRANSDUCER_MAP

arguments
    parameters    (1,1) struct
    medium_masks  (:,:,:) {mustBeNumeric}
    focus_pos     (1,3)   {mustBeNumeric}
    planimg       (1,1)   struct
    equip_param         struct = struct()
end

%% Output directory
out_dir = '';
if isfield(parameters.io, 'dir_babelbrain') && ~isempty(parameters.io.dir_babelbrain)
    out_dir = parameters.io.dir_babelbrain;
elseif isfield(parameters.io, 'dir_output') && ~isempty(parameters.io.dir_output)
    out_dir = fullfile(parameters.io.dir_output, 'babelbrain');
else
    out_dir = fullfile(parameters.io.dir_cache, 'babelbrain');
end
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

out_files = struct('mask_nii', '', 'weights_h5', '', 'run_script', '', 'config_yaml', '');

%% Resolve target name
target = strrep(parameters.io.output_affix, '_', '');
if isempty(target)
    target = sprintf('sub%03d_target', parameters.subject_id);
end

%% Resolve BabelBrain preset
bb = struct();
if isfield(parameters, 'babelbrain')
    bb = parameters.babelbrain;
end
bb = apply_default(bb, 'tx_system',          '');
bb = apply_default(bb, 'computing_backend',   3);
bb = apply_default(bb, 'device_name',         'Apple');
bb = apply_default(bb, 'ppw',                 6);
bb = apply_default(bb, 'do_refocusing',       false);
bb = apply_default(bb, 'python_env',          '');

% Domain size: derived from parameters.grid.acoustic_fov_length_mm, not from babelbrain.*
fov_len_mm = [];
if isfield(parameters, 'grid') && isfield(parameters.grid, 'acoustic_fov_length_mm') && ...
        isnumeric(parameters.grid.acoustic_fov_length_mm)
    fov_len_mm = parameters.grid.acoustic_fov_length_mm;
end
tight_narrow_beam = ~isempty(fov_len_mm);
if tight_narrow_beam
    focal_depth_mm = 0;
    if isfield(parameters, 'transducer') && isfield(parameters.transducer(1), 'focal_distance_mm')
        focal_depth_mm = parameters.transducer(1).focal_distance_mm;
    elseif isfield(parameters, 'transducer') && isfield(parameters.transducer(1), 'focal_distance_ep')
        focal_depth_mm = parameters.transducer(1).focal_distance_ep;
    end
    z_beyond_focus_mm = max(40, fov_len_mm - focal_depth_mm);
else
    z_beyond_focus_mm = 40;
end
bb = apply_default(bb, 'tx_mech_adjust_z_m',  0.0);

if isempty(bb.tx_system)
    [bb.tx_system, ~] = babelbrain_transducer_map(parameters, equip_param, out_dir);
end

%% Compute PPW
freq_hz  = parameters.transducer(1).freq_hz;
freq_kHz = round(freq_hz / 1e3);
res_mm   = parameters.grid.resolution_mm;
c_min    = 1500;
if isfield(parameters, 'medium_properties') && isfield(parameters.medium_properties, 'brain')
    c_min = parameters.medium_properties.brain.sound_speed;
end
ppw_auto = round(c_min / (freq_hz * res_mm * 1e-3));
ppw_use  = max(bb.ppw, ppw_auto);
subj_id  = sprintf('sub-%03d', parameters.subject_id);

%% 1. BabelViscoInput NIfTI
fprintf('[export_babelbrain_inputs] Writing tissue mask...\n');
out_files.mask_nii = write_babelbrain_input( ...
    medium_masks, focus_pos, parameters, planimg, out_dir, target);

%% 2. OptimizedWeightsFile (calibrated phases)
tr = parameters.transducer(1);
weights_path = fullfile(out_dir, [target '_CalibWeights.h5']);

if isfield(tr, 'annular') && isfield(tr.annular, 'elem_phase_rad') && isfield(tr.annular, 'elem_p_amp')
    phases_rad     = tr.annular.elem_phase_rad;
    elem_p_amp       = tr.annular.elem_p_amp(1);
    geo_phases_rad = [];
    if isfield(tr.annular, 'geo_phase_rad')
        geo_phases_rad = tr.annular.geo_phase_rad;
    elseif isfield(tr.annular, 'elem_id_mm') && isfield(tr.annular, 'elem_od_mm') && ...
            isfield(tr.annular, 'curv_radius_mm') && isfield(tr.annular, 'dist_geom_ep_mm') && ...
            isfield(tr, 'focal_distance_ep') && isfield(tr, 'freq_hz') && ...
            isfield(parameters, 'medium_properties') && ...
            isfield(parameters.medium_properties, 'water') && ...
            isfield(parameters.medium_properties.water, 'sound_speed')
        try
            tran_geo.n_elem             = tr.annular.elem_n;
            tran_geo.manufact           = '';
            tran_geo.transducer.freq_hz = tr.freq_hz;
            tran_geo.transducer.annular = tr.annular;
            geo_deg = compute_phases( ...
                parameters.medium_properties.water.sound_speed, ...
                tran_geo, tr.focal_distance_ep);
            geo_phases_rad = geo_deg * pi / 180;
            fprintf('[export_babelbrain_inputs] Geometric phases auto-computed at %.1f mm.\n', ...
                tr.focal_distance_ep);
        catch ME
            warning('export_babelbrain_inputs:geoPhases', ...
                'Could not auto-compute geometric phases: %s', ME.message);
        end
    end
    fprintf('[export_babelbrain_inputs] Writing calibration weights...\n');
    export_babelbrain_weights(phases_rad, elem_p_amp, parameters, weights_path, geo_phases_rad);
    out_files.weights_h5 = weights_path;
else
    fprintf('[export_babelbrain_inputs] No elem_phase_rad — skipping weights file.\n');
end

%% 3. Ready-to-run Python script
run_script_path = fullfile(out_dir, 'RunBabelBrain.py');
out_files.run_script = run_script_path;

nii_abs = out_files.mask_nii;
weights_arg = '';
if ~isempty(out_files.weights_h5) && isfile(out_files.weights_h5)
    weights_arg = sprintf('    OptimizedWeightsFile  = r''%s'',\n', out_files.weights_h5);
end

% Locate BabelBrain root for the script
prestus_root = get_prestus_path();
bb_root_candidates = {
    fullfile(prestus_root, 'external', 'BabelBrain'), ...
};
bb_root_py = '';
for k = 1:numel(bb_root_candidates)
    cand = strrep(fullfile(bb_root_candidates{k}), '\', '/');
    if exist(fullfile(bb_root_candidates{k}, 'TranscranialModeling'), 'dir')
        bb_root_py = cand;
        break
    end
end

fid = fopen(run_script_path, 'w');
fprintf(fid, '#!/usr/bin/env python3\n');
fprintf(fid, '"""\n');
fprintf(fid, 'RunBabelBrain.py — auto-generated by PRESTUS export_babelbrain_inputs\n');
fprintf(fid, 'Subject: %s  Target: %s  Preset: %s  Frequency: %d kHz\n', ...
    subj_id, target, bb.tx_system, freq_kHz);
fprintf(fid, '\n');
fprintf(fid, 'Before running:\n');
fprintf(fid, '  1. Activate the BabelBrain conda environment:\n');
fprintf(fid, '       conda activate BabelBrain\n');
fprintf(fid, '  2. Edit PYTHON_ENV / BabelBrain_ROOT below if needed.\n');
fprintf(fid, '  3. Run: python RunBabelBrain.py\n');
fprintf(fid, '"""\n\n');
fprintf(fid, 'import sys, os\n\n');
fprintf(fid, '# ── Configuration ─────────────────────────────────────────────\n');
fprintf(fid, 'BabelBrain_ROOT = r''%s''\n', bb_root_py);
fprintf(fid, 'if BabelBrain_ROOT and BabelBrain_ROOT not in sys.path:\n');
fprintf(fid, '    sys.path.insert(0, BabelBrain_ROOT)\n\n');
fprintf(fid, 'from TranscranialModeling.BabelIntegrationANNULAR_ARRAY import RUN_SIM\n\n');
fprintf(fid, 'R = RUN_SIM()\n');
fprintf(fid, 'output_files = R.RunCases(\n');
fprintf(fid, '    targets        = [r''%s''],\n', target);
fprintf(fid, '    ID             = r''%s'',\n', subj_id);
fprintf(fid, '    basedir        = r''%s'',\n', strrep(out_dir, '\', '/'));
fprintf(fid, '    Frequencies    = [%g],\n', freq_hz);
fprintf(fid, '    basePPW        = [%d],\n', ppw_use);
fprintf(fid, '    COMPUTING_BACKEND           = %d,\n', bb.computing_backend);
fprintf(fid, '    deviceName                  = r''%s'',\n', bb.device_name);
fprintf(fid, '    TxMechanicalAdjustmentZ     = %.6f,\n', bb.tx_mech_adjust_z_m);
fprintf(fid, '    bUseRayleighForWater        = True,\n');
fprintf(fid, '    bTightNarrowBeamDomain      = %s,\n', mat2pyb(tight_narrow_beam));
fprintf(fid, '    zLengthBeyonFocalPointWhenNarrow = %.4f,  # m\n', z_beyond_focus_mm/1e3);
fprintf(fid, '    bDoRefocusing               = %s,\n', mat2pyb(bb.do_refocusing));
fprintf(fid, '    bDisplay                    = False,\n');
fprintf(fid, '%s', weights_arg);

% Geometry kwargs (from preset default.yaml if available)
prestus_root_str = get_prestus_path();
yaml_path = fullfile(bb_root_py, 'BabelBrain', ...
    ['Babel_' get_babel_folder(bb.tx_system)], 'default.yaml');
yaml_path = strrep(yaml_path, '/', filesep);
if isfile(yaml_path)
    cfg = yaml.loadFile(yaml_path, 'ConvertToArray', true);
    if isfield(cfg, 'TxDiam')
        fprintf(fid, '    Aperture       = %.6f,  # m\n', cfg.TxDiam);
    end
    if isfield(cfg, 'FocalLength')
        fprintf(fid, '    FocalLength    = %.6f,  # m\n', cfg.FocalLength);
    end
    if isfield(cfg, 'InDiameters')
        fprintf(fid, '    InDiameters    = [%s],\n', ...
            strjoin(arrayfun(@(x) sprintf('%.6f', x), cfg.InDiameters(:)', 'UniformOutput', false), ', '));
    end
    if isfield(cfg, 'OutDiameters')
        fprintf(fid, '    OutDiameters   = [%s],\n', ...
            strjoin(arrayfun(@(x) sprintf('%.6f', x), cfg.OutDiameters(:)', 'UniformOutput', false), ', '));
    end
end

fprintf(fid, ')\n\n');
fprintf(fid, 'for f in (output_files or []):\n');
fprintf(fid, '    print(''Output:'', f)\n');
fclose(fid);
fprintf('[export_babelbrain_inputs] Run script: %s\n', run_script_path);

%% 4. Config YAML summary
config_yaml_path = fullfile(out_dir, 'babelbrain_config.yaml');
out_files.config_yaml = config_yaml_path;

fid = fopen(config_yaml_path, 'w');
fprintf(fid, '# BabelBrain run parameters — auto-generated by PRESTUS export_babelbrain_inputs\n');
fprintf(fid, '# Subject: %s  Target: %s\n\n', subj_id, target);
fprintf(fid, 'subject_id:   %s\n', subj_id);
fprintf(fid, 'target:       %s\n', target);
fprintf(fid, 'tx_system:    %s\n', bb.tx_system);
fprintf(fid, 'freq_hz:      %g\n', freq_hz);
fprintf(fid, 'ppw:          %d\n', ppw_use);
fprintf(fid, 'res_mm:       %g\n', res_mm);
fprintf(fid, 'computing_backend: %d\n', bb.computing_backend);
fprintf(fid, 'device_name:  %s\n', bb.device_name);
fprintf(fid, 'tx_mech_adjust_z_m: %.6f\n', bb.tx_mech_adjust_z_m);
fprintf(fid, 'z_beyond_focus_mm: %.1f\n', z_beyond_focus_mm);
fprintf(fid, 'do_refocusing: %s\n', mat2pyb(bb.do_refocusing));
fprintf(fid, '\nfiles:\n');
fprintf(fid, '  mask_nii:    %s\n', out_files.mask_nii);
if ~isempty(out_files.weights_h5)
    fprintf(fid, '  weights_h5:  %s\n', out_files.weights_h5);
end
fprintf(fid, '  run_script:  %s\n', run_script_path);

% Timing for thermal (in case user wants to run BabelBrain thermal too)
if isfield(parameters, 'timing') && isfield(parameters.timing, 'pd')
    t = parameters.timing;
    fprintf(fid, '\nthermal_inputs:\n');
    fprintf(fid, '  # Map to BabelBrain: DutyCycle = pd/pri, PRF = 1/pri, DurationUS = ptrd\n');
    if isfield(t, 'pd') && isfield(t, 'pri') && t.pri > 0
        fprintf(fid, '  DutyCycle:   %.4f  # pd=%.3f s / pri=%.3f s\n', t.pd/t.pri, t.pd, t.pri);
        fprintf(fid, '  PRF:         %.1f  # Hz = 1/pri\n', 1/t.pri);
    end
    if isfield(t, 'ptrd')
        fprintf(fid, '  DurationUS:  %.1f  # ptrd [s]\n', t.ptrd);
    end
end
fclose(fid);
fprintf('[export_babelbrain_inputs] Config YAML: %s\n', config_yaml_path);

fprintf('[export_babelbrain_inputs] Done. Files in: %s\n', out_dir);

end

%% Helpers

function s = apply_default(s, field, val)
    if ~isfield(s, field) || isempty(s.(field))
        s.(field) = val;
    end
end

function py_bool = mat2pyb(val)
    if val; py_bool = 'True'; else; py_bool = 'False'; end
end

function folder = get_babel_folder(tx_system)
    MAP = struct( ...
        'CTX_500',   'CTX500',   'CTX_250',   'CTX250',   'CTX_250_2ch', 'CTX250_2ch', ...
        'R15287',    'R15287',   'R15473',    'R15473',   'R15646',       'R15646', ...
        'DPX_500',   'DPX500',   'DPXPC_300', 'DPXPC300', 'R15148',       'R15148', ...
        'DomeTx',    'DomeTx',   'H317',      'H317',     'H246',         'H246', ...
        'H301',      'H301',     'Single',    'SingleTx', 'BSonix',       'SingleTx', ...
        'REMOPD',    'REMOPD',   'I12378',    'I12378',   'IGT64_500',    'IGT64_500', ...
        'ATAC',      'ATAC' ...
    );
    if isfield(MAP, tx_system)
        folder = MAP.(tx_system);
    else
        folder = tx_system;
    end
end
