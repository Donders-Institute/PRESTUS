function out_files = export_babelbrain_inputs( ...
        parameters, medium_masks, focus_pos, planimg)
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
%
% Input:
%   parameters    - PRESTUS config with .babelbrain sub-struct and .io.dir_babelbrain
%   medium_masks  - [Nx x Ny x Nz] medium label volume from preproc_medium_mask
%   focus_pos     - [1x3] focus voxel in grid coordinates
%   planimg       - planimg struct from grid_tissue_setup
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
bb = apply_default(bb, 'z_beyond_focus_mm',  40);
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
    if isfield(parameters, 'transducer') && isfield(parameters.transducer(1), 'focal_distance_ep') && ...
            ~isempty(parameters.transducer(1).focal_distance_ep)
        focal_depth_mm = parameters.transducer(1).focal_distance_ep;
    end
    z_beyond_focus_mm = max(bb.z_beyond_focus_mm, fov_len_mm - focal_depth_mm);
else
    z_beyond_focus_mm = bb.z_beyond_focus_mm;
end
bb = apply_default(bb, 'tx_mech_adjust_z_m',  0.0);

% 'auto' / named backends -> BabelBrain integer (RunCases and the YAML need a number)
[bb.computing_backend, bb.device_name] = resolve_bb_backend(bb.computing_backend, bb.device_name);

%% Compute PPW
freq_hz  = parameters.transducer(1).freq_hz;
freq_kHz = round(freq_hz / 1e3);
res_mm   = parameters.grid.resolution_mm;
% PPW target: BabelBrain's own spatial step is SmallestSOS/(freq*PPW), where
% SmallestSOS is floored by its built-in Pichardo shear-speed reference at
% density=1000 kg/m3 (1102.515 m/s) under the fluid-bone model -- see
% write_babelbrain_input.m for the full derivation. Using plain brain/water
% sound speed here (as this used to) targets the wrong reference and
% inflates the resampling write_babelbrain_input.m has to do.
pichardo_shear_floor_ms = 1000.0 * 0.422 + 680.515;
c_min = pichardo_shear_floor_ms;
if isfield(parameters, 'medium_properties')
    mp_names_ = fieldnames(parameters.medium_properties);
    for mi_ = 1:numel(mp_names_)
        mp_ = parameters.medium_properties.(mp_names_{mi_});
        if isfield(mp_, 'sound_speed') && mp_.sound_speed > 0
            c_min = min(c_min, mp_.sound_speed);
        end
    end
end
ppw_auto = round(c_min / (freq_hz * res_mm * 1e-3));
ppw_use  = max(bb.ppw, ppw_auto);
subj_id  = sprintf('sub-%03d', parameters.subject_id);

%% 1. BabelViscoInput NIfTI
fprintf('[export_babelbrain_inputs] Writing tissue mask...\n');
% Pass trans_pos so write_babelbrain_input can orient the depth flip to the
% actual beam direction (focus_pos(3) vs trans_pos(3)); empty falls back to the
% beam-+z assumption.
trans_pos_export = [];
if isfield(parameters.transducer(1), 'trans_pos') && ~isempty(parameters.transducer(1).trans_pos)
    trans_pos_export = parameters.transducer(1).trans_pos;
end
out_files.mask_nii = write_babelbrain_input( ...
    medium_masks, focus_pos, parameters, planimg, out_dir, target, ppw_use, [], trans_pos_export);

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

% Simulated geometry (PRESTUS transducer, or a named BabelBrain preset);
% its focal length also drives the placement correction.
if isempty(bb_root_py)
    error('export_babelbrain_inputs: BabelBrain source root not found (external/BabelBrain).');
end
bb_geom = babelbrain_tx_geometry(parameters.transducer(1), bb.tx_system, strrep(bb_root_py, '/', filesep));
bb.tx_system = bb_geom.tx_system;
preset_focal_length_m = bb_geom.focal_length_m;

% Geometry keyword arguments for RunCases: preset default.yaml, overridden by
% the PRESTUS values (PRESTUS geometry, or Single/BSonix aperture and focal length).
cfg = bb_geom.cfg;
geom_py = struct();
if isfield(cfg, 'TxDiam'),       geom_py.Aperture     = cfg.TxDiam;       end
if isfield(cfg, 'FocalLength'),  geom_py.FocalLength  = cfg.FocalLength;  end
if isfield(cfg, 'InDiameters'),  geom_py.InDiameters  = cfg.InDiameters;  end
if isfield(cfg, 'OutDiameters'), geom_py.OutDiameters = cfg.OutDiameters; end
flag_keys = struct('aperture', 'Aperture', 'focal_length', 'FocalLength', ...
                   'in_diameters', 'InDiameters', 'out_diameters', 'OutDiameters');
for k = 1:2:numel(bb_geom.flags)
    geom_py.(flag_keys.(erase(bb_geom.flags{k}, '--'))) = str2double(strsplit(bb_geom.flags{k+1}, ','));
end
if strcmp(bb_geom.tx_system, 'PRESTUS')
    run_sim_module = 'TranscranialModeling.babel_integration.integration_templates.babel_integration_focused_annular_array';
else
    run_sim_module = sprintf('TranscranialModeling.babel_integration.%s.babel_integration_%s', ...
        bb_geom.type, bb_geom.tx_system);
end

% Mechanical placement correction so BabelBrain's apex lands on trans_pos
% (same rule as babelbrain_wrapper; see babelbrain_placement_correction).
[geom_corr_m, apex_to_target_mm, placement_src] = babelbrain_placement_correction( ...
    parameters.transducer(1), focus_pos, trans_pos_export, res_mm, preset_focal_length_m);
tx_mech_adjust_total_m = bb.tx_mech_adjust_z_m + geom_corr_m;
if ~strcmp(placement_src, 'none')
    fprintf(['[export_babelbrain_inputs] Tx placement correction: %+.1f mm ' ...
        '(FL %.1f mm - apex-to-target %.2f mm from %s)\n'], ...
        geom_corr_m*1e3, preset_focal_length_m*1e3, apex_to_target_mm, placement_src);
else
    warning('export_babelbrain_inputs:noPlacementCorr', ...
        ['Preset FocalLength or apex-to-target distance unavailable; TxMechanicalAdjustmentZ ' ...
         'omits the placement correction and the BabelBrain field may be shifted axially.']);
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
fprintf(fid, 'BabelBrain_ROOT = ''%s''\n', py_escape(bb_root_py));
fprintf(fid, 'if BabelBrain_ROOT and BabelBrain_ROOT not in sys.path:\n');
fprintf(fid, '    sys.path.insert(0, BabelBrain_ROOT)\n');
fprintf(fid, '    sys.path.append(os.path.join(BabelBrain_ROOT, ''BabelBrain''))\n\n');
fprintf(fid, 'import numpy as np\n');
fprintf(fid, 'from %s import RUN_SIM\n\n', run_sim_module);
phased = any(strcmp(bb_geom.type, {'focused_array', 'flat_array_2D'}));
if phased
    % Element tables are long; the script reads them from the preset yaml.
    fprintf(fid, 'import yaml\n');
    fprintf(fid, 'with open(''%s'') as _f:\n', py_escape(strrep(bb_geom.yaml, '\', '/')));
    fprintf(fid, '    PRESET = yaml.safe_load(_f)\n\n');
end
fprintf(fid, 'R = RUN_SIM()\n');
fprintf(fid, 'output_files = R.RunCases(\n');
fprintf(fid, '    TxSystem       = ''%s'',\n', py_escape(bb_geom.tx_system));
fprintf(fid, '    targets        = [''%s''],\n', py_escape(target));
fprintf(fid, '    ID             = ''%s'',\n', py_escape(subj_id));
fprintf(fid, '    basedir        = ''%s'',\n', py_escape(strrep(out_dir, '\', '/')));
fprintf(fid, '    Frequencies    = [%g],\n', freq_hz);
fprintf(fid, '    basePPW        = [%d],\n', ppw_use);
fprintf(fid, '    COMPUTING_BACKEND           = %d,\n', bb.computing_backend);
fprintf(fid, '    deviceName                  = ''%s'',\n', py_escape(bb.device_name));
fprintf(fid, '    TxMechanicalAdjustmentZ     = %.6f,  # m: %.1f mm placement corr + %.1f mm user\n', ...
    tx_mech_adjust_total_m, geom_corr_m*1e3, bb.tx_mech_adjust_z_m*1e3);
if isfield(bb, 'use_rayleigh_for_water') && ~isempty(bb.use_rayleigh_for_water)
    use_rayleigh_water = logical(bb.use_rayleigh_for_water);
else
    use_rayleigh_water = false;   % PRESTUS default: FDTD water for backend comparability
end
fprintf(fid, '    bUseRayleighForWater        = %s,\n', mat2pyb(use_rayleigh_water));
fprintf(fid, '    bTightNarrowBeamDomain      = %s,\n', mat2pyb(tight_narrow_beam));
fprintf(fid, '    zLengthBeyonFocalPointWhenNarrow = %.4f,  # m\n', z_beyond_focus_mm/1e3);
fprintf(fid, '    bDoRefocusing               = %s,\n', mat2pyb(bb.do_refocusing));
fprintf(fid, '    bDisplay                    = False,\n');
fprintf(fid, '%s', weights_arg);

% Geometry kwargs [m]
if isfield(geom_py, 'Aperture')
    fprintf(fid, '    Aperture       = %.7f,\n', geom_py.Aperture);
end
if isfield(geom_py, 'FocalLength')
    fprintf(fid, '    FocalLength    = %.7f,\n', geom_py.FocalLength);
end
for f = {'InDiameters', 'OutDiameters'}
    if isfield(geom_py, f{1})
        fprintf(fid, '    %-14s = np.array([%s]),\n', f{1}, ...
            strjoin(arrayfun(@(x) sprintf('%.7f', x), geom_py.(f{1})(:)', 'UniformOutput', false), ', '));
    end
end
if phased
    fprintf(fid, '    elements          = PRESET[''Elements''],\n');
    fprintf(fid, '    num_elements      = PRESET[''NumElements''],\n');
    fprintf(fid, '    element_size      = PRESET[''ElementSize''],\n');
    fprintf(fid, '    coordinate_system = PRESET.get(''CoordinateSystem'', ''cartesian''),\n');
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
fprintf(fid, 'tx_mech_adjust_z_m: %.6f  # total = %.1f mm placement corr (%s) + %.1f mm user\n', ...
    tx_mech_adjust_total_m, geom_corr_m*1e3, placement_src, bb.tx_mech_adjust_z_m*1e3);
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

function s = py_escape(s)
% Text for a single-quoted Python string literal: no quote or newline in a
% value can end the literal and run the rest as Python.
    s = char(string(s));
    s = strrep(s, '\', '\\');
    s = strrep(s, '''', '\''');
    s = strrep(s, newline, '\n');
    s = strrep(s, char(13), '\r');
end
