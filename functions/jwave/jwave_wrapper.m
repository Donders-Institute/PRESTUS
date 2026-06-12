function sensor_data = jwave_wrapper(parameters, kgrid, medium, source, source_labels)
% JWAVE_WRAPPER  Run acoustic simulation via jwave (JAX-based Python backend)
%
% Translates PRESTUS kgrid/medium/source structs to HDF5, invokes
% run_jwave_sim.py, and reads back the pressure field in the same
% sensor_data format returned by acoustic_simulation.m.
%
% Only the direct-geometry source path (use_kWaveArray == 0) is supported.
% Axisymmetric grids are not supported. Thermal simulation is not supported.
%
% Use as:
%   sensor_data = jwave_wrapper(parameters, kgrid, medium, source, source_labels)
%
% Input:
%   parameters    - PRESTUS config; must contain parameters.jwave (see below)
%   kgrid         - kWaveGrid with time axis set (kgrid.t_array ready)
%   medium        - k-Wave medium struct (sound_speed, density, alpha_coeff,
%                   alpha_power); thermal fields are silently ignored
%   source        - struct with p_mask [grid_dims logical] and p [N_src x Nt]
%   source_labels - [grid_dims] int array (element IDs, 0 = inactive); used
%                   for provenance only, not required for the simulation
%
% parameters.jwave fields (all optional with defaults):
%   .python_env     - path to CPU Python venv/conda env root or bin dir, or '' for system python3
%   .python_env_gpu - path to GPU (CUDA) Python env bin dir; leave empty to always use CPU
%   .use_gpu        - [] auto-detect via nvidia-smi (default) | true force GPU | false force CPU
%   .solver         - 'time_domain' (default) | 'helmholtz'
%   .pml_size       - PML thickness in grid points (default: 20)
%   .gmres_restart  - GMRES restart parameter for Helmholtz solver (default: 40)
%   .gmres_maxiter  - Max GMRES outer iterations for Helmholtz solver (default: 300)
%   .precision      - 'single' (default) | 'double'
%   .dry_run        - true to write inputs and exit without simulating (default false)
%
% Output:
%   sensor_data.p_max_all   [grid_dims] double [Pa]  peak positive pressure
%   sensor_data.p_min_all   [grid_dims] double [Pa]  peak negative pressure
%   sensor_data.p_final     [grid_dims] double [Pa]  pressure at final time step
%   sensor_data.jwave_input_h5   char  path to HDF5 input written
%   sensor_data.jwave_output_h5  char  path to HDF5 output read
%
% See also: BABELBRAIN_WRAPPER, ACOUSTIC_SIMULATION, PRESTUS_PIPELINE

arguments
    parameters    (1,1) struct
    kgrid         (1,1)
    medium        (1,1) struct
    source        (1,1) struct
    source_labels = []
end

% -------------------------------------------------------------------------
% 0. Ensure HDF5 plugin path is valid before any h5read/h5write call
% -------------------------------------------------------------------------
ensure_hdf5_plugin_path();

% -------------------------------------------------------------------------
% 1. Resolve parameters
% -------------------------------------------------------------------------

jw = struct();
if isfield(parameters, 'jwave')
    jw = parameters.jwave;
end

python_env     = get_field(jw, 'python_env',     '');
python_env_gpu = get_field(jw, 'python_env_gpu', '');
use_gpu        = get_field(jw, 'use_gpu',        []);   % [] = auto-detect
solver         = get_field(jw, 'solver',          'time_domain');
pml_size       = get_field(jw, 'pml_size',        20);
gmres_restart  = get_field(jw, 'gmres_restart',   40);
gmres_maxiter  = get_field(jw, 'gmres_maxiter',   300);
precision      = get_field(jw, 'precision',       'single');
dry_run        = get_field(jw, 'dry_run',         false);

% Resolve python executable, selecting GPU env when appropriate
[python_bin, using_gpu_env] = resolve_python_env(python_env, python_env_gpu, use_gpu);

% Locate the Python bridge script (same folder as this .m file)
this_dir   = fileparts(mfilename('fullpath'));
bridge_py  = fullfile(this_dir, 'run_jwave_sim.py');
if ~isfile(bridge_py)
    error('jwave_wrapper: run_jwave_sim.py not found at:\n  %s\n%s', ...
        bridge_py, 'This file must be in functions/jwave/ alongside jwave_wrapper.m.');
end

% Validate solver choice
valid_solvers = {'time_domain', 'helmholtz'};
if ~ismember(solver, valid_solvers)
    error('jwave_wrapper: Unknown solver ''%s''. Choose from: %s', ...
        solver, strjoin(valid_solvers, ', '));
end

% -------------------------------------------------------------------------
% 2. Prepare output paths
% -------------------------------------------------------------------------

out_dir    = char(parameters.io.dir_cache);
affix      = char(parameters.io.output_affix);
subj_id    = parameters.subject_id;
medium_tag = parameters.simulation.medium;

input_h5  = fullfile(out_dir, sprintf('sub-%03d_%s_jwave_input%s.h5',  subj_id, medium_tag, affix));
output_h5 = fullfile(out_dir, sprintf('sub-%03d_%s_jwave_output%s.h5', subj_id, medium_tag, affix));

% -------------------------------------------------------------------------
% 3. Write HDF5 input
% -------------------------------------------------------------------------

fprintf('[jwave_wrapper] Writing simulation input: %s\n', input_h5);
write_jwave_input(input_h5, kgrid, medium, source, parameters, pml_size, precision);

if dry_run
    fprintf('[jwave_wrapper] Dry run — skipping simulation.\n');
    sensor_data.p_max_all      = [];
    sensor_data.p_min_all      = [];
    sensor_data.p_final        = [];
    sensor_data.jwave_input_h5 = input_h5;
    sensor_data.jwave_output_h5 = '';
    return;
end

% -------------------------------------------------------------------------
% 4. Invoke Python bridge
% -------------------------------------------------------------------------

if using_gpu_env
    jax_prefix = '';
else
    jax_prefix = 'JAX_PLATFORMS=cpu ';
end

% Append Python stdout+stderr directly to the main PRESTUS log (no separate log file).
if isfield(parameters, 'io') && isfield(parameters.io, 'log_file') && ~isempty(parameters.io.log_file)
    log_redirect = sprintf(' >> "%s" 2>&1', parameters.io.log_file);
else
    log_redirect = '';
end

gmres_flags = '';
if strcmp(solver, 'helmholtz')
    gmres_flags = sprintf(' --gmres_restart %d --gmres_maxiter %d', gmres_restart, gmres_maxiter);
end
cmd = sprintf('%s"%s" "%s" --input "%s" --output "%s" --solver %s%s%s', ...
    jax_prefix, python_bin, bridge_py, input_h5, output_h5, solver, gmres_flags, log_redirect);

fprintf('[jwave_wrapper] Launching jwave simulation...\n');
fprintf('[jwave_wrapper] Command:\n  %s\n', cmd);

t_start = tic;
[status, cmd_out] = system(cmd);
elapsed = toc(t_start);

if ~isempty(cmd_out)
    disp(cmd_out);
end

if status ~= 0
    error('jwave_wrapper: Python process exited with status %d.\n%s', status, cmd_out);
end

fprintf('[jwave_wrapper] jwave finished in %.1f s.\n', elapsed);

% -------------------------------------------------------------------------
% 5. Read HDF5 output
% -------------------------------------------------------------------------

if ~isfile(output_h5)
    error('jwave_wrapper: Output HDF5 not found.\n  Expected: %s', output_h5);
end

fprintf('[jwave_wrapper] Reading output: %s\n', output_h5);
sensor_data = read_jwave_output(output_h5, parameters);

sensor_data.jwave_input_h5  = input_h5;
sensor_data.jwave_output_h5 = output_h5;

% Verify output dimensions match PRESTUS grid
expected_dims = parameters.grid.dims(:)';
for chk_ = {'p_max_all', 'p_min_all', 'p_final'}
    arr = sensor_data.(chk_{1});
    if ~isempty(arr) && ~isequal(size(arr), expected_dims)
        error(['jwave_wrapper: output %s has size [%s] but expected PRESTUS grid [%s].\n' ...
            'Check dimension transpositions in run_jwave_sim.py.'], ...
            chk_{1}, num2str(size(arr)), num2str(expected_dims));
    end
end

end


% =========================================================================
% Local helper: write_jwave_input
% =========================================================================

function write_jwave_input(h5_path, kgrid, medium, source, parameters, pml_size, precision)
% Write kgrid, medium, and source to HDF5 for run_jwave_sim.py.

if isfile(h5_path)
    delete(h5_path);
end

cast_fn = @(x) cast(gather(x), precision);

% --- Grid ---
dims = numel(parameters.grid.dims);
h5create(h5_path, '/grid/N',  [1 dims]);
h5write( h5_path, '/grid/N',  int32(parameters.grid.dims(:)'));
h5create(h5_path, '/grid/dx', [1 dims]);
h5write( h5_path, '/grid/dx', cast_fn(parameters.grid.resolution_mm / 1e3 * ones(1, dims)));

% --- Time axis ---
Nt = kgrid.Nt;
dt = kgrid.dt;
h5create(h5_path, '/time/Nt', [1 1]);
h5write( h5_path, '/time/Nt', int32(Nt));
h5create(h5_path, '/time/dt', [1 1]);
h5write( h5_path, '/time/dt', cast_fn(dt));

% --- Medium ---
thermal_fields = {'thermal_conductivity', 'specific_heat', 'perfusion_coeff'};
for f = thermal_fields
    if isfield(medium, f{1}), medium = rmfield(medium, f{1}); end
end

grid_sz = parameters.grid.dims(:)';
write_medium_field(h5_path, 'sound_speed', medium.sound_speed, grid_sz, cast_fn);
write_medium_field(h5_path, 'density',     medium.density,     grid_sz, cast_fn);
write_medium_field(h5_path, 'alpha_coeff', medium.alpha_coeff, grid_sz, cast_fn);

alpha_power = get_alpha_power(medium, grid_sz);
write_medium_field(h5_path, 'alpha_power', alpha_power, grid_sz, cast_fn);

% --- Source ---
% p_mask: store linear indices of active source voxels
src_indices = int32(find(source.p_mask(:))' - 1);  % 0-based for Python
h5create(h5_path, '/source/mask_indices', size(src_indices));
h5write( h5_path, '/source/mask_indices', src_indices);

% Source signal: [N_src x Nt]
h5create(h5_path, '/source/p', size(source.p));
h5write( h5_path, '/source/p', cast_fn(source.p));

% Record start index (1-based MATLAB → 0-based Python)
record_start = 1;
if isfield(parameters, 'grid') && isfield(parameters.grid, 'record_start_index')
    record_start = parameters.grid.record_start_index;
end
h5create(h5_path, '/source/record_start_index', [1 1]);
h5write( h5_path, '/source/record_start_index', int32(record_start - 1));

% --- PML ---
h5create(h5_path, '/pml/size', [1 1]);
h5write( h5_path, '/pml/size', int32(pml_size));

% --- Metadata ---
h5writeatt(h5_path, '/', 'prestus_version', 'PRESTUS');
h5writeatt(h5_path, '/', 'dims',            int32(dims));
h5writeatt(h5_path, '/', 'freq_hz',         parameters.transducer(1).freq_hz);

end


function write_medium_field(h5_path, name, val, grid_sz, cast_fn)
% Write a medium field that may be scalar or a full spatial array.
dset = ['/medium/' name];
if isscalar(val)
    h5create(h5_path, dset, [1 1]);
    h5write( h5_path, dset, cast_fn(val));
else
    h5create(h5_path, dset, grid_sz);
    h5write( h5_path, dset, cast_fn(val));
end
end


function alpha_power = get_alpha_power(medium, grid_sz)
if isfield(medium, 'alpha_power')
    alpha_power = medium.alpha_power;
else
    alpha_power = 2.0;  % k-Wave fractional Laplacian default
end
% Broadcast scalar to match expectations on Python side
if isscalar(alpha_power)
    alpha_power = alpha_power * ones(grid_sz);
end
end


% =========================================================================
% Local helper: read_jwave_output
% =========================================================================

function sensor_data = read_jwave_output(h5_path, parameters)
% Read p_max_all, p_min_all, p_final from jwave HDF5 output.

sensor_data = struct();

fields = {'p_max_all', 'p_min_all', 'p_final'};
for f = fields
    key = ['/pressure/' f{1}];
    try
        sensor_data.(f{1}) = double(h5read(h5_path, key));
    catch ME
        warning('jwave_wrapper: Could not read %s from %s:\n  %s', key, h5_path, ME.message);
        sensor_data.(f{1}) = [];
    end
end

end


% =========================================================================
% Utility
% =========================================================================

function [python_bin, is_gpu] = resolve_python_env(cpu_env, gpu_env, use_gpu)
% Select CPU or GPU Python environment and return the full path to the
% python3 executable, and whether the GPU env was chosen.
%
% use_gpu == true  : force GPU env (error if gpu_env is empty)
% use_gpu == false : force CPU env
% use_gpu == []    : auto-detect via nvidia-smi; use GPU env if a device is
%                    found and gpu_env is non-empty, otherwise fall back to
%                    CPU env with a warning

    if islogical(use_gpu) && use_gpu
        % Explicit GPU request
        if isempty(gpu_env)
            error('jwave_wrapper: parameters.jwave.use_gpu = true but python_env_gpu is not set.');
        end
        fprintf('[jwave_wrapper] GPU env selected (use_gpu = true).\n');
        chosen_env = gpu_env;
        is_gpu = true;

    elseif islogical(use_gpu) && ~use_gpu
        % Explicit CPU request
        chosen_env = cpu_env;
        is_gpu = false;

    else
        % Auto-detect: probe nvidia-smi, but also honour CUDA_VISIBLE_DEVICES.
        % nvidia-smi lists all node GPUs regardless of scheduler allocation;
        % CUDA_VISIBLE_DEVICES="" or "NoDevFiles" means no GPU was granted to
        % this job even if hardware is present.
        [smi_status, smi_out] = system('nvidia-smi -L 2>/dev/null');
        cvd = getenv('CUDA_VISIBLE_DEVICES');
        cvd_blocked = ~isempty(cvd) && (isempty(strtrim(cvd)) || strcmpi(strtrim(cvd), 'nodevfiles'));
        smi_out_trim = strtrim(smi_out);
        % Strip trailing punctuation so "No devices found." matches the check
        smi_out_norm = regexprep(smi_out_trim, '[.!]+$', '');
        smi_lists_device = (smi_status == 0) && ~isempty(smi_out_norm) && ...
            ~strcmpi(smi_out_norm, 'no devices found') && ...
            ~strcmpi(smi_out_norm, 'no devices were found') && ...
            ~contains(lower(smi_out_norm), 'no devices');
        gpu_found = smi_lists_device && ~cvd_blocked;

        if gpu_found && ~isempty(gpu_env)
            fprintf('[jwave_wrapper] GPU detected (%s) — using GPU env.\n', ...
                smi_out_trim(1:min(end,60)));
            chosen_env = gpu_env;
            is_gpu = true;
        else
            if gpu_found && isempty(gpu_env)
                warning('jwave_wrapper:noGpuEnv', ...
                    'GPU detected but python_env_gpu is not set — falling back to CPU env.');
            end
            chosen_env = cpu_env;
            is_gpu = false;
        end
    end

    % Resolve chosen_env path to a python3 executable
    if isempty(chosen_env)
        python_bin = 'python3';
    elseif isfile(chosen_env)
        % caller passed the binary directly
        python_bin = chosen_env;
    elseif isfile(fullfile(chosen_env, 'bin', 'python3'))
        % venv / conda root
        python_bin = fullfile(chosen_env, 'bin', 'python3');
    elseif isfile(fullfile(chosen_env, 'python3'))
        % bin directory passed directly
        python_bin = fullfile(chosen_env, 'python3');
    else
        warning('jwave_wrapper:pythonEnvNotFound', ...
            'Could not find python3 under ''%s'' — falling back to system python3.', chosen_env);
        python_bin = 'python3';
    end
end


function val = get_field(s, name, default)
if isfield(s, name)
    val = s.(name);
else
    val = default;
end
end
