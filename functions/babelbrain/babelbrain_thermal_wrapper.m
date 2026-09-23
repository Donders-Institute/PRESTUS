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
%   DutyCycle = pd / pri
%   PRF       = 1  / pri   [Hz]
%   DurationUS = ptrd       [s]  (pulse train repetition duration)
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
% explicitly, this wrapper therefore passes the in-brain peak Isptp of
% sensor_data.p_max_all (the same calibrated/scaled field k-Wave thermal
% consumes) -- NOT the free-water calibration ISPPA, which would drop the
% transcranial loss and overestimate heating.
%
% Output:
%   results_heating - struct with fields matching thermal_simulation output:
%     .dx_step          grid voxel size [m]
%     .p_amp            [Nx x Ny x Nz] peak pressure [Pa], resampled to PRESTUS grid
%     .T                [Nx x Ny x Nz] temperature rise at end of sonication [°C]
%     .T_end            [Nx x Ny x Nz] final temperature (after cooling) [°C]
%     .CEM43            [Nx x Ny x Nz] thermal dose [CEM43]
%     .CEM43_end        [Nx x Ny x Nz] final CEM43
%     .TempProfileTarget [1 x Nt]     temperature at focal target over time [°C]
%     .TimeProfileTarget [1 x Nt]     time axis [s]
%     .MaxIsppa         peak ISPPA in brain [W/cm²]
%     .MaxIspta         peak ISPTA in brain [W/cm²]
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
% PRESTUS timing hierarchy (see thermal_parameters.m): pd/pri (single pulse
% within a train) -> ptd (ON duration of ONE pulse train) -> ptri (period of
% one pulse-train repetition, ON+OFF) -> ptrd (TOTAL span across ALL
% repetitions, = n_ptri_reps * ptri) -> post_ptri_dur (extra cooldown after
% the whole repeated sequence, not currently threaded through to BabelBrain
% -- see warning below).
t = parameters.timing;
if ~isfield(t,'pd') || ~isfield(t,'pri') || t.pri <= 0
    error('babelbrain_thermal_wrapper: parameters.timing.pd and .pri must be set (pri > 0).');
end
if ~isfield(t, 'ptd')
    error('babelbrain_thermal_wrapper: parameters.timing.ptd must be set (single pulse-train ON duration).');
end
duty_cycle = t.pd / t.pri;
prf_hz     = 1 / t.pri;
duration_s = t.ptd;   % ON duration of a single pulse-train repetition (NOT ptrd, the total across all reps)

% Validate converted timing parameters before passing to BabelBrain
if duty_cycle <= 0 || duty_cycle > 1
    error(['babelbrain_thermal_wrapper: duty_cycle=%.4f is outside (0, 1].\n' ...
        'Check parameters.timing.pd=%.4f s and .pri=%.4f s.'], duty_cycle, t.pd, t.pri);
end
if prf_hz <= 0 || ~isfinite(prf_hz)
    error(['babelbrain_thermal_wrapper: PRF=%.2f Hz is invalid (<= 0 or non-finite).\n' ...
        'Check parameters.timing.pri=%.4f s.'], prf_hz, t.pri);
end
if duration_s <= 0
    error(['babelbrain_thermal_wrapper: duration=%.4f s is <= 0.\n' ...
        'Check parameters.timing.ptd.'], duration_s);
end
if ~isfinite(bb.isppa_w_per_cm2) || bb.isppa_w_per_cm2 <= 0
    error(['babelbrain_thermal_wrapper: ISPPA=%.4f W/cm² is invalid (<= 0 or non-finite).\n' ...
        'Check sensor_data.p_max_all or set babelbrain.isppa_w_per_cm2.'], bb.isppa_w_per_cm2);
end

% Cooling interval between repetitions (BabelBrain DurationOff).
% In PRESTUS, ptri is the period of one pulse-train repetition (ON+OFF); the
% off-time is the remainder after the ON duration (ptd). If ptri is not set
% or is not longer than ptd, there is no inter-train gap.
if isfield(t, 'ptri') && t.ptri > t.ptd
    duration_off_s = t.ptri - t.ptd;
else
    duration_off_s = 0.0;
end

% Number of pulse-train repetitions: PRESTUS computes this as
% n_ptri_reps = round(ptrd / ptri) (thermal_parameters.m). BabelBrain's
% Repetitions parameter runs this many (DurationUS on / DurationOff off) cycles,
% so cooling between trains is modelled.
repetitions = 1;
if isfield(t, 'n_ptri_reps') && t.n_ptri_reps > 1
    repetitions = round(t.n_ptri_reps);
end

% NOTE: parameters.timing.post_ptri_dur (cooldown after the whole sequence in
% the native thermal path) is not passed to BabelBrain. The candidate options,
% NumberGroupedSonications/PauseBetweenGroupedSonications, are left at their
% defaults (1, 0.0) until their semantics are confirmed, so FinalTemp/FinalDose
% reflect the end of the last repetition only.
if isfield(t, 'post_ptri_dur') && t.post_ptri_dur > 0
    warning('PRESTUS:BabelBrain:postPtriDurNotApplied', ...
        ['parameters.timing.post_ptri_dur (%.1f s) is set but not currently passed to BabelBrain -- ' ...
         'its thermal output reflects cooldown only up to the end of the last repetition, not the ' ...
         'full post_ptri_dur cooldown window the native thermal path would track.'], t.post_ptri_dur);
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
    '--duty_cycle', sprintf('%.6f', duty_cycle), ...
    '--prf', sprintf('%.2f', prf_hz), ...
    '--duration_s', sprintf('%.4f', duration_s), ...
    '--duration_off_s', sprintf('%.4f', duration_off_s), ...
    '--repetitions', num2str(repetitions), ...
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

fprintf('[babelbrain_thermal_wrapper] Timing: DC=%.3f  PRF=%.0f Hz  duration=%.1f s  duration_off=%.1f s  repetitions=%d\n', ...
    duty_cycle, prf_hz, duration_s, duration_off_s, repetitions);
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
% BabelBrain names the thermal file by replacing the DataForSim suffix with
% -ThermalField_AllCombinations.h5.  Try an exact stem match first; fall back
% to a glob + newest-by-mtime if the convention has changed.
[acous_dir_, acous_stem_] = fileparts(acoustic_h5);
stem_prefix_ = regexprep(acous_stem_, '_DataForSim$', '');
stem_match_  = fullfile(acous_dir_, [stem_prefix_ '-ThermalField_AllCombinations.h5']);
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
% focus_pos and focal_length_mm are stored in sensor_data by babelbrain_wrapper
% so that the thermal reader can anchor BB's internal coordinate system.
if isfield(sensor_data, 'focus_pos') && isfield(sensor_data, 'focal_length_mm')
    thermal_focus_pos       = sensor_data.focus_pos;
    thermal_focal_length_mm = sensor_data.focal_length_mm;
else
    error(['babelbrain_thermal_wrapper: sensor_data is missing focus_pos or focal_length_mm.\n' ...
        'These are set by babelbrain_wrapper — re-run babelbrain_wrapper to refresh sensor_data.']);
end
thermal_trans_pos = [];
if isfield(sensor_data, 'trans_pos')
    thermal_trans_pos = sensor_data.trans_pos;   % beam-axis direction for the depth resample
end
results_heating = read_babelbrain_thermal(thermal_h5, parameters, planimg, thermal_focus_pos, thermal_focal_length_mm, thermal_trans_pos);

end

%% Helpers

function [backend, device] = resolve_backend(backend, device)
% Resolve 'auto' computing_backend / device_name to concrete integers/strings.
% Priority: Metal (macOS) > CUDA (nvidia-smi) > error.
%
% BabelBrain's thermal engine (CalculateThermalProcess.py) maps backends 1..4
% to 'CUDA','OpenCL','Metal','MLX' and has no CPU option (unlike the acoustic
% engine), so backend 4 means MLX (Apple Silicon only). Without a GPU this
% errors rather than falling back to a CPU path that does not exist.
    if ~(ischar(backend) && strcmp(backend, 'auto'))
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
