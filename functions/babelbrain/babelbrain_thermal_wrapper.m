function results_heating = babelbrain_thermal_wrapper( ...
        parameters, sensor_data, planimg, kwave_medium, medium_masks)
% BABELBRAIN_THERMAL_WRAPPER  Run BabelBrain thermal simulation and read results
%
% Calls BabelBrain's CalculateThermalProcess via the Python bridge, reads the
% output HDF5, and returns a results_heating struct compatible with PRESTUS
% thermal_analysis and nifti_thermal so that all downstream post-processing
% (focal metrics, CEM43, NIfTI export, HTML report) runs unchanged.
%
% The thermal simulation uses the *DataForSim.h5 acoustic output from
% babelbrain_wrapper (sensor_data.babelbrain_source). BabelBrain applies its
% own Pennes BHTE solver on the same domain as the acoustic simulation.
%
% Timing parameters are read from parameters.timing and converted to the
% BabelBrain convention:
%   DutyCycle   = pd / pri
%   PRF         = 1  / pri        [Hz]
%   DurationUS  = ptd             [s]  (ON duration of one pulse train)
%   DurationOff = ptri - ptd      [s]  (gap after each train)
%   Repetitions = n_ptri_reps
%   PauseBetweenGroupedSonications = post_ptri_dur [s] (cooldown after the
%                 last train; NumberGroupedSonications = 1)
% See babelbrain_thermal_timing.
%
% After the simulation, the 3D temperature and CEM43 maps are resampled from
% BabelBrain's simulation grid back into the PRESTUS grid space using the
% affine stored in the DataForSim.h5 file, so that PRESTUS thermal analysis
% can operate on them.
%
% Use as:
%   results_heating = babelbrain_thermal_wrapper(parameters, sensor_data, planimg)
%   results_heating = babelbrain_thermal_wrapper(..., kwave_medium)
%   results_heating = babelbrain_thermal_wrapper(..., kwave_medium, medium_masks)
%
% Input:
%   parameters    - PRESTUS config; uses .timing, .thermal, .babelbrain
%   sensor_data   - struct from babelbrain_wrapper; must contain
%                   .babelbrain_source (path to *DataForSim.h5)
%   planimg       - planimg struct (affine; used for NIfTI resampling)
%   kwave_medium  - (optional) PRESTUS medium struct on the ACOUSTIC grid
%                   (same grid as sensor_data.p_max_all); its sound_speed and
%                   density give the rho*c used for the in-situ ISPPA
%   medium_masks  - (optional) PRESTUS label volume on the acoustic grid; used
%                   for the brain mask when sensor_data carries no
%                   babelbrain_material_map
%
% ISPPA semantics: BabelBrain's BHTE rescales the field so that the peak
% ISPPA INSIDE THE BRAIN equals the Isppa it receives
% (PressureRatio = sqrt(Isppa*2*rho*c) / max(p_amp(brain)) in
% CalculateTemperatureEffects.py). Unless babelbrain.isppa_w_per_cm2 is set
% explicitly, this wrapper therefore passes the in-brain peak ISPPA of
% sensor_data.p_max_all (the same calibrated/scaled field k-Wave thermal
% consumes) -- NOT the free-water calibration ISPPA, which would drop the
% transcranial loss and overestimate heating.
%
% Output:
%   results_heating - struct with fields matching thermal_simulation output
%                     (see READ_BABELBRAIN_THERMAL), on the PRESTUS grid:
%     .p_amp            peak pressure [Pa]
%     .maxT             temperature at end of sonication [°C]
%     .heating_endT     final temperature (after cooling) [°C]
%     .CEM43, .CEM43_end  thermal dose at end of sonication / after cooling [min]
%     .CEM43_iso, .CEM43_iso_end  NaN (BabelBrain computes no ISO CEM43)
%     .focal_planeT, .focal_planeCEM43, .focal_planeCEM43_iso  empty
%     .timeseries       temperature time courses (.TempProfileTarget,
%                       .TimeProfileTarget, per-tissue .T and .Tdiff)
%     .MaxIsppa, .MaxIspta, .Isppa, .Ispta  BabelBrain intensity metrics [W/cm²]
%     .dx_step          grid voxel size [m]
%     .babelbrain_thermal_source  path to thermal HDF5 output
%
% See also: BABELBRAIN_WRAPPER, THERMAL_SIMULATION, THERMAL_ANALYSIS

arguments
    parameters   (1,1) struct
    sensor_data  (1,1) struct
    planimg      (1,1) struct
    kwave_medium (1,1) struct = struct()
    medium_masks = []
end

%% Validate acoustic source
if ~isfield(sensor_data, 'babelbrain_source') || ...
        isempty(sensor_data.babelbrain_source) || ...
        strcmp(sensor_data.babelbrain_source, 'dry_run')
    error(['babelbrain_thermal_wrapper: sensor_data.babelbrain_source is missing or ' ...
        'is a dry-run placeholder. Run babelbrain_wrapper first with dry_run=false.']);
end

acoustic_h5 = sensor_data.babelbrain_source;
if ~isfile(acoustic_h5)
    error('babelbrain_thermal_wrapper: acoustic HDF5 not found: %s', acoustic_h5);
end

%% Resolve BabelBrain config
bb = struct();
if isfield(parameters, 'babelbrain')
    bb = parameters.babelbrain;
end
if isfield(bb, 'medium') && strcmpi(char(bb.medium), 'prestus')
    % The BHTE reads tissue regions and thermal properties from BabelBrain's
    % tissue labels, which babelbrain.medium 'prestus' replaces by material classes
    error('babelbrain_thermal_wrapper:prestusMedium', ...
        ['babelbrain.medium ''prestus'' field: BabelBrain''s thermal solver needs its tissue labels. ' ...
         'Use thermal_backend ''kwave'' for heating.']);
end
bb = apply_default(bb, 'computing_backend',  'auto');
bb = apply_default(bb, 'device_name',        'auto');
bb = apply_default(bb, 'python_env',         '');
bb = apply_default(bb, 'root',               '');
if isfield(bb, 'isppa_w_per_cm2') && ~isempty(bb.isppa_w_per_cm2)
    fprintf(['[babelbrain_thermal_wrapper] Using explicit babelbrain.isppa_w_per_cm2 = %.2f W/cm² ' ...
        '(interpreted by BabelBrain as the in-situ peak ISPPA in brain).\n'], bb.isppa_w_per_cm2);
else
    bb.isppa_w_per_cm2 = babelbrain_insitu_isppa(parameters, sensor_data, kwave_medium, medium_masks);
    ff_isppa = get_tx_isppa(parameters.transducer(1));
    if isfinite(ff_isppa) && ff_isppa > 0
        fprintf(['[babelbrain_thermal_wrapper] In-situ brain peak ISPPA = %.2f W/cm² ' ...
            '(free-water %.2f W/cm²; transmission %.1f%%) -> BabelBrain Isppa.\n'], ...
            bb.isppa_w_per_cm2, ff_isppa, 100 * bb.isppa_w_per_cm2 / ff_isppa);
    else
        fprintf('[babelbrain_thermal_wrapper] In-situ brain peak ISPPA = %.2f W/cm² -> BabelBrain Isppa.\n', ...
            bb.isppa_w_per_cm2);
    end
end

[bb.computing_backend, bb.device_name] = resolve_backend(bb.computing_backend, bb.device_name);

%% Locate Python binary
python_bin = resolve_python_env(bb.python_env, 'python');
ensure_hdf5_plugin_path(python_bin);

%% Locate runner script
this_dir  = fileparts(mfilename('fullpath'));
bb_runner = fullfile(this_dir, 'run_babelbrain_thermal.py');
if ~isfile(bb_runner)
    error(['babelbrain_thermal_wrapper: run_babelbrain_thermal.py not found at:\n  %s\n' ...
        'This file should be in functions/babelbrain/ alongside babelbrain_thermal_wrapper.m.'], ...
        bb_runner);
end

%% Convert PRESTUS timing → BabelBrain convention
% ptd -> DurationUS, ptri-ptd -> DurationOff, n_ptri_reps -> Repetitions,
% post_ptri_dur -> PauseBetweenGroupedSonications (single group); see
% babelbrain_thermal_timing.
bt = babelbrain_thermal_timing(parameters.timing);

if ~isfinite(bb.isppa_w_per_cm2) || bb.isppa_w_per_cm2 <= 0
    error(['babelbrain_thermal_wrapper: ISPPA=%.4f W/cm² is invalid (<= 0 or non-finite).\n' ...
        'Check sensor_data.p_max_all or set babelbrain.isppa_w_per_cm2.'], bb.isppa_w_per_cm2);
end

% Select pressure field: use refocused map when refocusing was requested.
if isfield(bb, 'do_refocusing') && bb.do_refocusing
    sel_p = 'p_amp_refocus';
else
    sel_p = 'p_amp';
end

baseline_temp = 37.0;
if isfield(parameters, 'thermal') && isfield(parameters.thermal, 'temp_0') && ...
        isfield(parameters.thermal.temp_0, 'brain')
    baseline_temp = parameters.thermal.temp_0.brain;
end

%% Per-run log file
log_dir = fullfile(fileparts(acoustic_h5), 'log');
if ~isfolder(log_dir), mkdir(log_dir); end
[~, acoustic_stem, ~] = fileparts(strrep(acoustic_h5, '.h5', ''));
log_file = fullfile(log_dir, [acoustic_stem '_thermal.log']);
fprintf('[babelbrain_thermal_wrapper] Log file: %s\n', log_file);

%% Build Python command
% run_command invokes python directly via argv, no shell at all -- every
% path/string reaches the process as a literal argument regardless of its
% contents.
thermal_argv = {python_bin, bb_runner, ...
    '--acoustic_h5', acoustic_h5, ...
    '--isppa', sprintf('%.4f', bb.isppa_w_per_cm2), ...
    '--duty_cycle', sprintf('%.6f', bt.duty_cycle), ...
    '--prf', sprintf('%.2f', bt.prf_hz), ...
    '--duration_s', sprintf('%.4f', bt.duration_s), ...
    '--duration_off_s', sprintf('%.4f', bt.duration_off_s), ...
    '--repetitions', num2str(bt.repetitions), ...
    '--pause_after_s', sprintf('%.4f', bt.pause_after_s), ...
    '--sel_p', sel_p, ...
    '--computing_backend', num2str(bb.computing_backend), ...
    '--device_name', bb.device_name, ...
    '--baseline_temp', sprintf('%.2f', baseline_temp), ...
    '--log_file', log_file};

if ~isempty(bb.root)
    thermal_argv = [thermal_argv, {'--babelbrain_root', bb.root}];
end

if isfield(parameters, 'transducer') && ~isempty(parameters.transducer) && ...
        isfield(parameters.transducer(1), 'freq_hz')
    thermal_argv = [thermal_argv, {'--frequency', sprintf('%.0f', parameters.transducer(1).freq_hz)}];
end

fprintf(['[babelbrain_thermal_wrapper] Timing: DC=%.3f  PRF=%.0f Hz  duration=%.1f s  ' ...
    'duration_off=%.1f s  repetitions=%d  pause_after=%.1f s\n'], ...
    bt.duty_cycle, bt.prf_hz, bt.duration_s, bt.duration_off_s, bt.repetitions, bt.pause_after_s);
fprintf('[babelbrain_thermal_wrapper] ISPPA: %.2f W/cm²\n', bb.isppa_w_per_cm2);
fprintf('[babelbrain_thermal_wrapper] Command:\n  %s\n', strjoin(thermal_argv, ' '));

t0     = tic;
[status, cmd_out] = run_command(thermal_argv, 'Echo', true);
elapsed = toc(t0);

if status ~= 0
    error('babelbrain_thermal_wrapper: Python process exited with status %d.\nOutput:\n%s', status, cmd_out);
end
fprintf('[babelbrain_thermal_wrapper] BabelBrain thermal finished in %.1f s.\n', elapsed);

%% Locate thermal output HDF5
% BabelBrain names the thermal file <acoustic stem>-ThermalField_AllCombinations.h5
% (the _DataForSim part is kept, e.g. x_300kHz_6PPW_DataForSim-ThermalField_...).
% Try that exact name first; fall back to a glob + newest-by-mtime if the
% convention has changed.
[acous_dir_, acous_stem_] = fileparts(acoustic_h5);
stem_prefix_ = regexprep(acous_stem_, '_DataForSim$', '');
stem_match_  = fullfile(acous_dir_, [acous_stem_ '-ThermalField_AllCombinations.h5']);
thermal_h5   = '';
if isfile(stem_match_)
    thermal_h5 = stem_match_;
else
    thermal_candidates = dir(fullfile(acous_dir_, [stem_prefix_ '*ThermalField*.h5']));
    if isempty(thermal_candidates)
        % Wider fallback: any ThermalField file in the same directory
        thermal_candidates = dir(fullfile(acous_dir_, '*ThermalField*.h5'));
    end
    if ~isempty(thermal_candidates)
        [~, newest_idx] = max([thermal_candidates.datenum]);
        thermal_h5 = fullfile(thermal_candidates(newest_idx).folder, thermal_candidates(newest_idx).name);
        if numel(thermal_candidates) > 1
            warning('babelbrain_thermal_wrapper:ambiguousThermalFile', ...
                ['%d *ThermalField*.h5 files found; using newest:\n  %s\n' ...
                'Delete stale files to avoid ambiguity.'], numel(thermal_candidates), thermal_h5);
        end
    end
end

if isempty(thermal_h5)
    error(['babelbrain_thermal_wrapper: thermal output HDF5 not found in %s.\n' ...
        'Expected a file matching *ThermalField*.h5'], acous_dir_);
end
fprintf('[babelbrain_thermal_wrapper] Reading thermal output: %s\n', thermal_h5);

% BabelBrain writes thermal HDF5 with blosc compression; rewrite to gzip so
% MATLAB's h5read can load it without a blosc filter plugin.
rewriter = fullfile(this_dir, 'rewrite_hdf5_no_blosc.py');
if isfile(rewriter)
    % run_command invokes python directly via argv, no shell at all.
    [st_rw, out_rw] = run_command({python_bin, rewriter, thermal_h5});
    if st_rw ~= 0
        warning('babelbrain_thermal_wrapper:bloscRewrite', ...
            'Thermal HDF5 blosc-rewrite failed (exit %d):\n%s\n  h5read may fail.', ...
            st_rw, out_rw);
    end
end

%% Read thermal output and resample to PRESTUS grid
% The thermal volumes are resampled onto the THERMAL grid (parameters =
% parameters_th from thermal_stage), so the anchor positions and step must be
% in thermal-grid coordinates. sensor_data.focus_pos / trans_pos are on the
% ACOUSTIC grid (acoustic-FOV crop coordinates when an acoustic FOV was used,
% which thermal_stage restores to the full grid for the thermal run; a separate
% thermal grid may also be offset / at a different step), so they are used only
% as a fallback. focal_length_mm (fallback z anchor) comes from sensor_data.
if isfield(sensor_data, 'focal_length_mm')
    thermal_focal_length_mm = sensor_data.focal_length_mm;
else
    error(['babelbrain_thermal_wrapper: sensor_data is missing focal_length_mm.\n' ...
        'It is set by babelbrain_wrapper — re-run babelbrain_wrapper to refresh sensor_data.']);
end
[thermal_focus_pos, thermal_trans_pos, params_rd] = thermal_grid_anchor(parameters, sensor_data);
results_heating = read_babelbrain_thermal(thermal_h5, params_rd, planimg, thermal_focus_pos, thermal_focal_length_mm, thermal_trans_pos);

end

%% Helpers

function [focus_pos, trans_pos, params_rd] = thermal_grid_anchor(parameters, sensor_data)
% Focus / transducer positions and grid step of the thermal grid.
% thermal_grid_setup stores the unrounded positions for a separate thermal grid
% (grid.thermal_focus_pos_exact / thermal_trans_pos_exact); otherwise the
% thermal grid is the acoustic grid in thermal_stage's coordinates, whose
% transducer(1).focus_pos / trans_pos are the right ones.
params_rd = parameters;
g = parameters.grid;
if isfield(g, 'thermal_focus_pos_exact') && ~isempty(g.thermal_focus_pos_exact)
    focus_pos = g.thermal_focus_pos_exact;
    trans_pos = g.thermal_trans_pos_exact;
    if isfield(g, 'thermal_resolution_mm') && ~isempty(g.thermal_resolution_mm)
        params_rd.grid.resolution_mm = g.thermal_resolution_mm;
    end
elseif isfield(parameters.transducer(1), 'focus_pos') && ~isempty(parameters.transducer(1).focus_pos)
    focus_pos = parameters.transducer(1).focus_pos;
    trans_pos = parameters.transducer(1).trans_pos;
elseif isfield(sensor_data, 'focus_pos')
    warning('babelbrain_thermal_wrapper:acousticAnchor', ...
        ['No thermal-grid focus_pos; anchoring on the acoustic-grid sensor_data.focus_pos. ' ...
         'Correct only when the thermal grid equals the acoustic grid.']);
    focus_pos = sensor_data.focus_pos;
    trans_pos = [];
    if isfield(sensor_data, 'trans_pos'); trans_pos = sensor_data.trans_pos; end
else
    error('babelbrain_thermal_wrapper:noAnchor', 'No focus position available for the thermal resample.');
end
focus_pos = double(focus_pos(:)');
trans_pos = double(trans_pos(:)');
end

function [backend, device] = resolve_backend(backend, device)
% Resolve 'auto' computing_backend / device_name to concrete integers/strings.
% Priority: Metal (macOS) > CUDA (nvidia-smi) > error.
%
% BabelBrain's thermal engine (CalculateThermalProcess.py) maps backends 1..4
% to 'CUDA','OpenCL','Metal','MLX' and has no CPU option (unlike the acoustic
% engine), so backend 4 means MLX (Apple Silicon only). Without a GPU this
% errors rather than falling back to a CPU path that does not exist.
    if isstring(backend); backend = char(backend); end
    if ischar(backend) && ~strcmpi(backend, 'auto')
        % Named backend ('cuda' | 'opencl' | 'metal' | 'cpu') -> 1-4, as for
        % the acoustic run; the device is left to the logic below.
        backend = resolve_bb_backend(backend, 'unused');
    end
    if ~(ischar(backend) && strcmpi(backend, 'auto'))
        if backend == 4
            error('babelbrain_thermal_wrapper:noCpuBackend', ...
                ['computing_backend=4 requested for the BabelBrain thermal engine, but BabelBrain''s ' ...
                 'thermal Backend list (CalculateThermalProcess.py) has no CPU option -- index 4 is MLX ' ...
                 '(Apple Silicon only). Use 1=CUDA, 2=OpenCL, or 3=Metal, matching available hardware.']);
        end
        if ischar(device) && strcmp(device, 'auto')
            switch backend
                case 3;  device = 'Apple';
                otherwise; device = 'CPU';
            end
        end
        return
    end
    if ismac
        backend = 3;
        device  = 'Apple';
        fprintf('[babelbrain_thermal_wrapper] Auto-detected backend: Metal (macOS)\n');
        return
    end
    [status, out] = run_command({'nvidia-smi', '--query-gpu=name', '--format=csv,noheader'});
    if status == 0
        lines = strtrim(strsplit(strtrim(out), newline));
        lines = lines(~cellfun(@isempty, lines));
        if ~isempty(lines)
            backend = 1;
            device  = lines{1};
            fprintf('[babelbrain_thermal_wrapper] Auto-detected backend: CUDA, device="%s"\n', device);
            return
        end
    end
    error('babelbrain_thermal_wrapper:noBackendAvailable', ...
        ['No GPU detected (not macOS, no CUDA device via nvidia-smi) and BabelBrain''s thermal engine ' ...
         'has no CPU backend -- it cannot run here. Set babelbrain.computing_backend explicitly if a ' ...
         'supported device (CUDA/OpenCL/Metal) is actually available but was not auto-detected.']);
end

function s = apply_default(s, field, val)
    if ~isfield(s, field) || isempty(s.(field))
        s.(field) = val;
    end
end

function flag = bb_root_flag(bb_root)
    if ~isempty(bb_root)
        flag = sprintf('--babelbrain_root "%s"', bb_root);
    else
        flag = '';
    end
end

function isppa = get_tx_isppa(tr)
% Read ISPPA from tx_calibration.isppa_wcm2 with fallback to legacy field.
    if isfield(tr, 'tx_calibration') && isfield(tr.tx_calibration, 'isppa_wcm2') && ...
            isscalar(tr.tx_calibration.isppa_wcm2) && isfinite(tr.tx_calibration.isppa_wcm2)
        isppa = tr.tx_calibration.isppa_wcm2;
    elseif isfield(tr, 'set_intensity_w_per_cm2') && isscalar(tr.set_intensity_w_per_cm2)
        isppa = tr.set_intensity_w_per_cm2;
    else
        isppa = NaN;
    end
end
