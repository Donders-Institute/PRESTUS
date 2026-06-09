function results_heating = babelbrain_thermal_wrapper( ...
        parameters, sensor_data, planimg, kwave_medium)
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
%
% Input:
%   parameters    - PRESTUS config; uses .timing, .thermal, .babelbrain
%   sensor_data   - struct from babelbrain_wrapper; must contain
%                   .babelbrain_source (path to *DataForSim.h5)
%   planimg       - planimg struct (affine; used for NIfTI resampling)
%   kwave_medium  - (optional) PRESTUS medium struct; unused here but kept
%                   for interface compatibility with thermal_simulation
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
    kwave_medium (1,1) struct = struct()  %#ok<INUSA>
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
bb = apply_default(bb, 'computing_backend',  3);
bb = apply_default(bb, 'device_name',        'Apple');
bb = apply_default(bb, 'python_env',         '');
bb = apply_default(bb, 'root',               '');
bb = apply_default(bb, 'isppa_w_per_cm2',    parameters.transducer(1).set_intensity_w_per_cm2);

%% Locate Python binary
if ~isempty(bb.python_env)
    python_bin = fullfile(bb.python_env, 'bin', 'python');
    if ~isfile(python_bin)
        python_bin = fullfile(bb.python_env, 'python');
    end
else
    python_bin = 'python3';
end

%% Locate runner script
this_dir  = fileparts(mfilename('fullpath'));
bb_runner = fullfile(this_dir, 'run_babelbrain_thermal.py');
if ~isfile(bb_runner)
    error(['babelbrain_thermal_wrapper: run_babelbrain_thermal.py not found at:\n  %s\n' ...
        'This file should be in functions/babelbrain/ alongside babelbrain_thermal_wrapper.m.'], ...
        bb_runner);
end

%% Convert PRESTUS timing → BabelBrain convention
t = parameters.timing;
if ~isfield(t,'pd') || ~isfield(t,'pri') || t.pri <= 0
    error('babelbrain_thermal_wrapper: parameters.timing.pd and .pri must be set (pri > 0).');
end
duty_cycle = t.pd / t.pri;
prf_hz     = 1 / t.pri;
duration_s = t.ptrd;   % pulse train repetition duration

% Cooling interval between repetitions (BabelBrain DurationOff).
% In PRESTUS, ptri is the period of one pulse-train block; the off-time is
% the remainder after the on-duration (ptrd). If ptri is not set or is not
% longer than ptrd, there is no inter-train gap.
if isfield(t, 'ptri') && t.ptri > t.ptrd
    duration_off_s = t.ptri - t.ptrd;
else
    duration_off_s = 0.0;
end

% Warn for multi-PTRI protocols: DurationUS maps to one block (ptrd), not
% the total duration across all repetitions.
if isfield(t, 'n_ptri_reps') && t.n_ptri_reps > 1
    warning('PRESTUS:BabelBrain:multiPTRI', ...
        ['BabelBrain DurationUS maps to one pulse-train block (ptrd = %.1f s). ' ...
         'Multi-PTRI protocol (%d reps) uses BabelBrain Repetitions parameter ' ...
         'for the full sonication — verify the timing in the BabelBrain log.'], ...
        t.ptrd, t.n_ptri_reps);
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

%% Build Python command
cmd = sprintf([ ...
    '"%s" "%s" '  ...
    '--acoustic_h5 "%s" '  ...
    '--isppa %.4f '  ...
    '--duty_cycle %.6f '  ...
    '--prf %.2f '  ...
    '--duration_s %.4f '  ...
    '--duration_off_s %.4f '  ...
    '--sel_p "%s" '  ...
    '--computing_backend %d '  ...
    '--device_name "%s" '  ...
    '%s '], ...
    python_bin, bb_runner, ...
    acoustic_h5, ...
    bb.isppa_w_per_cm2, duty_cycle, prf_hz, duration_s, duration_off_s, sel_p, ...
    bb.computing_backend, bb.device_name, ...
    bb_root_flag(bb.root));

if ~isempty(bb.root)
    cmd = [cmd sprintf('--babelbrain_root "%s" ', bb.root)];
end

fprintf('[babelbrain_thermal_wrapper] Timing: DC=%.3f  PRF=%.0f Hz  duration=%.1f s  duration_off=%.1f s\n', ...
    duty_cycle, prf_hz, duration_s, duration_off_s);
fprintf('[babelbrain_thermal_wrapper] ISPPA: %.2f W/cm²\n', bb.isppa_w_per_cm2);
fprintf('[babelbrain_thermal_wrapper] Command:\n  %s\n', cmd);

t0     = tic;
status = system(cmd);
elapsed = toc(t0);

if status ~= 0
    error('babelbrain_thermal_wrapper: Python process exited with status %d.', status);
end
fprintf('[babelbrain_thermal_wrapper] BabelBrain thermal finished in %.1f s.\n', elapsed);

%% Locate thermal output HDF5
% BabelBrain names output: same prefix as DataForSim but with thermal suffix
thermal_candidates = dir(fileparts(acoustic_h5));
thermal_h5 = '';
for k = 1:numel(thermal_candidates)
    if contains(thermal_candidates(k).name, 'ThermalField') && ...
            endsWith(thermal_candidates(k).name, '.h5')
        thermal_h5 = fullfile(thermal_candidates(k).folder, thermal_candidates(k).name);
        break
    end
end

if isempty(thermal_h5)
    error(['babelbrain_thermal_wrapper: thermal output HDF5 not found in %s.\n' ...
        'Expected a file matching *ThermalField*.h5'], fileparts(acoustic_h5));
end
fprintf('[babelbrain_thermal_wrapper] Reading thermal output: %s\n', thermal_h5);

%% Read thermal output and resample to PRESTUS grid
results_heating = read_babelbrain_thermal(thermal_h5, parameters, planimg);

end

%% Helpers

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
