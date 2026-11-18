function [sensor_data] = acoustic_simulation(kgrid, medium, source, sensor, input_args, parameters, p_complex_freqs)
% ACOUSTIC_SIMULATION  Dispatch a k-Wave acoustic simulation across supported compute backends
%
% Selects and calls the appropriate k-Wave solver function based on
% parameters.simulation.code_type ('matlab_cpu', 'matlab_gpu', 'cpp_cpu',
% 'cpp_gpu', 'cpp_mlx') and grid dimensionality (2D, 2D axisymmetric, or 3D). With
% acoustic_backend 'kwave_elastic', hands off to ACOUSTIC_SIMULATION_ELASTIC
% (pstdElastic2D/3D, shear waves in the skull). Strips
% thermal-only medium fields before passing to k-Wave, sets DataCast/precision
% for GPU/C++ paths, and configures HDF5 file I/O for C++ backends.
%
% With p_complex_freqs (io.save_p_complex), also returns the complex
% steady-state amplitude at those frequencies (STEADY_STATE_COMPLEX_AMPLITUDE)
% as sensor_data.p_complex_fit [n_points x n_freq], and the recording window
% as sensor_data.p_complex_window. The same method on every k-Wave build and
% the elastic backend:
%   one frequency  - the pressure time series is recorded over the last
%                    quarter period only (memory: ~ppp/4 samples per point
%                    instead of steady_state_num_periods periods). That window
%                    also bounds p_max_all / p_min_all, so they are replaced
%                    by +/- the fitted amplitude.
%   several        - the steady-state window as configured (tones are only
%                    separated over several periods); p_max_all / p_min_all
%                    stay k-Wave's; stops above 8 GB of recording
%                    (PRESTUS:pComplexTooLarge).
%
% Use as:
%   sensor_data = acoustic_simulation(kgrid, medium, source, sensor, input_args, parameters)
%   sensor_data = acoustic_simulation(kgrid, medium, source, sensor, input_args, parameters, p_complex_freqs)
%
% Input:
%   kgrid      - simulation grid object
%   medium     - k-Wave medium with sound_speed, density, alpha_coeff,
%                thermal_conductivity, specific_heat, perfusion_coeff (thermal fields
%                are removed internally before k-Wave is called)
%   source     - k-Wave source definition (p_mask, p)
%   sensor     - k-Wave sensor definition (mask, record)
%   input_args - additional k-Wave key-value options (PMLInside, PMLSize, etc.)
%   parameters - PRESTUS config with simulation.code_type, simulation.interactive,
%                grid.dims, io.dir_output, transducer, simulation.precision
%   p_complex_freqs - optional frequencies [Hz] for sensor_data.p_complex_fit ([] = none)
%
% Output:
%   sensor_data - k-Wave sensor output (p_max_all, p_final, etc.) [Pa]
%
% See also: ACOUSTIC_WRAPPER, ACOUSTIC_SIMULATION_ELASTIC, KSPACEFIRSTORDER3D,
%           KSPACEFIRSTORDER2D, KSPACEFIRSTORDERAS

arguments
    kgrid      (1,1)
    medium     (1,1) struct
    source     (1,1) struct
    sensor     (1,1) struct
    input_args (1,1) struct
    parameters (1,1) struct
    p_complex_freqs double = []
end

% Select which k-Wave code to run
%   1: MATLAB CPU code 'matlab_cpu'
%   2: MATLAB GPU code 'matlab_gpu'
%   3: C++ CPU code 'cpp_cpu'
%   4: C++ GPU code 'cpp_gpu'
%   5: C++ Apple Silicon GPU code 'cpp_mlx' (kspaceFirstOrder-MLX)

% Disable k-Wave's PlotSim when not interactive or when using C++ backends
% (which run out-of-process). Interactive MATLAB runs keep the real-time
% animation.
cpp_backend = ismember(parameters.simulation.code_type, {'cpp_cpu', 'cpp_gpu', 'cpp_mlx'});
if ~parameters.simulation.interactive || cpp_backend
    input_args.PlotSim = false;
    if parameters.simulation.interactive
        disp('Interactive acoustic visualization is not available for C++ backends (cpp_cpu/cpp_gpu/cpp_mlx).');
    end
end

% Remove fields that are not recognized for acoustic simulations
thermal_fields = {'thermal_conductivity', 'specific_heat', 'perfusion_coeff'};
for f = thermal_fields
    if isfield(medium, f{1})
        medium = rmfield(medium, f{1});
    end
end

% Define scale for plotting
input_args.PlotScale = [-1, 1] * parameters.transducer(1).(parameters.transducer(1).type).elem_p_amp(1);

% Complex amplitude (io.save_p_complex): record the pressure time series
p_window = [];
if ~isempty(p_complex_freqs)
    if numel(p_complex_freqs) == 1
        % last quarter period (at least 4 samples), ending at the last step
        ppp = 1 / (p_complex_freqs * kgrid.dt);
        n_rec = min(kgrid.Nt, max(4, ceil(ppp / 4) + 1));
        sensor.record_start_index = kgrid.Nt - n_rec + 1;
        p_window = struct('method', 'quarter_period', 'record_start_index', sensor.record_start_index, ...
            'n_samples', n_rec, 'points_per_period', ppp, 'p_max_from_fit', true);
    else
        record_start = 1;
        if isfield(sensor, 'record_start_index') && ~isempty(sensor.record_start_index)
            record_start = double(sensor.record_start_index);
        end
        n_rec = kgrid.Nt - record_start + 1;
        record_gb = nnz(sensor.mask) * n_rec * 4 / 1e9;   % single precision
        max_record_gb = 8;
        if record_gb > max_record_gb
            error('PRESTUS:pComplexTooLarge', ...
                ['The complex pressure fields of a multi-frequency setup need the steady-state window of the pressure ' ...
                 'time series in memory: %.1f GB (%d points x %d time steps); the limit is %d GB. Use a partial acoustic FOV ' ...
                 '(grid.acoustic_fov_diameter_mm, grid.acoustic_fov_length_mm), a coarser grid or fewer periods ' ...
                 '(grid.steady_state_num_periods).'], record_gb, nnz(sensor.mask), n_rec, max_record_gb);
        end
        p_window = struct('method', 'steady_state_window', 'record_start_index', record_start, ...
            'n_samples', n_rec, 'points_per_period', 1 ./ (p_complex_freqs * kgrid.dt), 'p_max_from_fit', false);
    end
    if ~ismember('p', sensor.record)
        sensor.record{end+1} = 'p';
    end
end

% Elastic (shear-wave) skull model: pstdElastic2D/3D instead of kspaceFirstOrder*.
% Also taken by the free-water baseline, so both runs share one forward model.
if isfield(parameters, 'acoustic_backend') && strcmpi(parameters.acoustic_backend, 'kwave_elastic')
    sensor_data = acoustic_simulation_elastic(kgrid, medium, source, sensor, input_args, parameters);
    sensor_data = fit_complex(sensor_data, kgrid, p_complex_freqs, p_window);
    return
end

% k-Wave C++ backends write large scratch .h5 files (input + output) that are
% removed on success (DeleteData=true). Their names are deterministic
% (subject + medium + output_affix), so two jobs sharing an affix — e.g. a
% re-launched run while a prior one is still in flight — would otherwise fight
% over the same file and fail with an HDF5 "cannot create file" iostream error.
% Prefer per-job node-local scratch: io.dir_scratch (HPC_SCRATCH_SETUP, with
% hpc.tmp_gb > 0), else $TMPDIR (set by SLURM/qsub). Both are unique per job
% and keep the 100s-of-MB I/O off the shared (network) project mount. Falls
% back to the cache dir for local runs where neither is set.
kwave_scratch_dir = char(get_nonempty_field(parameters, {'io', 'dir_scratch'}, ''));
if isempty(kwave_scratch_dir) || ~isfolder(kwave_scratch_dir)
    kwave_scratch_dir = getenv('TMPDIR');
end
if isempty(kwave_scratch_dir) || ~isfolder(kwave_scratch_dir)
    kwave_scratch_dir = char(parameters.io.dir_cache);
end

% Select submission based on code type
switch parameters.simulation.code_type
    case 'cpp_cpu'

      % Force precision for C++ HDF5 compatibility
      medium  = cast_struct(medium, parameters.simulation.precision);
      source  = cast_struct(source, parameters.simulation.precision);
      sensor  = cast_struct(sensor, parameters.simulation.precision);

      % Pathname for the input and output files (used only for non-interactive computations)
      h5_dir = kwave_scratch_dir;
      input_args.SaveToDisk = char(fullfile(h5_dir, ...
         sprintf('sub-%03d_%s_input%s.h5', parameters.subject_id, ...
         parameters.simulation.medium, parameters.io.output_affix)));

      input_args.DataName = sprintf('kwave_sub-%03d%s', parameters.subject_id, parameters.io.output_affix);
      input_args.DataPath = h5_dir;
      input_args.DeleteData = true;

      if numel(parameters.grid.dims) == 3
         input_args.FunctionName = 'kspaceFirstOrder3D';
      elseif numel(parameters.grid.dims) == 2 && isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1
         input_args.FunctionName = 'kspaceFirstOrderAS';
         % kspaceFirstOrder3DC must know the grid is axisymmetric so it correctly
         % crops only the outer-edge PML from the radial dimension (y1=1, not y1=pml+1).
         input_args.Axisymmetric = true;
      else
         input_args.FunctionName = 'kspaceFirstOrder2D';
      end

      n_threads = cpu_thread_count();
      if ~isempty(n_threads), input_args.NumThreads = n_threads; end

      input_args_cell = zip_fields(input_args);
      sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args_cell{:});

   case 'cpp_gpu'

      if numel(parameters.grid.dims) == 2 && isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1
         % No CUDA binary supports axisymmetric geometry; fall back to matlab_gpu path
         % (kspaceFirstOrderAS with GPU DataCast) which preserves GPU acceleration.
         % Do NOT cast medium/source/sensor to CPU precision first — DataCast handles
         % both precision conversion and GPU transfer in one step, avoiding an
         % intermediate in-memory copy of every medium array at full precision.
         warn('kspaceFirstOrder-CUDA does not support axisymmetric simulations. Falling back to matlab_gpu (kspaceFirstOrderAS with gpuArray).');
         as_args = zip_fields(rmfield(input_args, intersect(fieldnames(input_args), ...
             {'SaveToDisk','DataName','DataPath','DeleteData','FunctionName','DeviceNum','NumThreads'})));
         as_args_cast = [as_args, {'DataCast', ['gpuArray-', char(parameters.simulation.precision)]}];
         sensor_data = kspaceFirstOrderAS(kgrid, medium, source, sensor, as_args_cast{:}, 'RadialSymmetry', 'WSWA-FFT');
      else
         % Force precision for C++ HDF5 compatibility (only needed for the C++ path;
         % the axisymmetric fallback above uses DataCast instead).
         medium  = cast_struct(medium, parameters.simulation.precision);
         source  = cast_struct(source, parameters.simulation.precision);
         sensor  = cast_struct(sensor, parameters.simulation.precision);

         % Pathname for the input and output files
         h5_dir = kwave_scratch_dir;
         input_args.SaveToDisk = char(fullfile(h5_dir, ...
           sprintf('sub-%03d_%s_input%s.h5', parameters.subject_id, ...
           parameters.simulation.medium, parameters.io.output_affix)));
         input_args.DataName = sprintf('kwave_sub-%03d%s', parameters.subject_id, parameters.io.output_affix);
         input_args.DataPath = h5_dir;
         input_args.DeleteData = true;
         input_args.FunctionName = 'kspaceFirstOrder3D';
         gpu_id = str2double(getenv('SLURM_LOCALID'));
         input_args.DeviceNum = gpu_id;
         input_args.NumThreads = 1;
         input_args_cell = zip_fields(input_args);
         sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args_cell{:});
      end
      
   case 'cpp_mlx'

      if numel(parameters.grid.dims) == 2 && isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1
         % kspaceFirstOrder-MLX is a port of the CUDA solver, which has no axisymmetric
         % geometry; use the MATLAB solver on the CPU (no gpuArray on macOS).
         warn('kspaceFirstOrder-MLX does not support axisymmetric simulations. Falling back to matlab_cpu (kspaceFirstOrderAS).');
         as_args = zip_fields(rmfield(input_args, intersect(fieldnames(input_args), ...
             {'SaveToDisk','DataName','DataPath','DeleteData','FunctionName','DeviceNum','NumThreads','BinaryName'})));
         as_args_cast = [as_args, {'DataCast', char(parameters.simulation.precision)}];
         sensor_data = kspaceFirstOrderAS(kgrid, medium, source, sensor, as_args_cast{:}, 'RadialSymmetry', 'WSWA-FFT');
      else
         % The user supplies the binary, as for the other C++ backends: both files
         % must sit in k-Wave's binaries folder (MLX loads mlx.metallib from the
         % folder of the executable).
         mlx_dir = getkWavePath('binaries');
         for f = {'kspaceFirstOrder-MLX', 'mlx.metallib'}
            if ~isfile(fullfile(mlx_dir, f{1}))
               error('prestus:mlxBinaryMissing', ...
                  ['code_type cpp_mlx needs %s in %s. Build kspaceFirstOrder-MLX for Apple ' ...
                   'Silicon and copy it and mlx.metallib there (see doc_kwave.md).'], f{1}, mlx_dir);
            end
         end

         % Force precision for C++ HDF5 compatibility
         medium  = cast_struct(medium, parameters.simulation.precision);
         source  = cast_struct(source, parameters.simulation.precision);
         sensor  = cast_struct(sensor, parameters.simulation.precision);

         h5_dir = kwave_scratch_dir;
         input_args.SaveToDisk = char(fullfile(h5_dir, ...
           sprintf('sub-%03d_%s_input%s.h5', parameters.subject_id, ...
           parameters.simulation.medium, parameters.io.output_affix)));
         input_args.DataName = sprintf('kwave_sub-%03d%s', parameters.subject_id, parameters.io.output_affix);
         input_args.DataPath = h5_dir;
         input_args.DeleteData = true;
         if numel(parameters.grid.dims) == 3
            input_args.FunctionName = 'kspaceFirstOrder3D';
         else
            input_args.FunctionName = 'kspaceFirstOrder2D';
         end
         % One Apple GPU: no device number or thread count (those are CUDA / OpenMP options).
         input_args.BinaryName = 'kspaceFirstOrder-MLX';
         input_args_cell = zip_fields(input_args);
         sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args_cell{:});
      end

   case 'matlab_gpu'

      input_args.DataCast = ['gpuArray-', char(parameters.simulation.precision)];
      input_args_cell = zip_fields(input_args);
      
      if numel(parameters.grid.dims) == 3
         sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args_cell{:});
      elseif numel(parameters.grid.dims) == 2 && isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1
         sensor_data = kspaceFirstOrderAS(kgrid, medium, source, sensor, input_args_cell{:}, 'RadialSymmetry', 'WSWA-FFT');
      else % 2D simulation
         sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args_cell{:});
      end
      
   case 'matlab_cpu'

      n_threads = cpu_thread_count();
      if ~isempty(n_threads), maxNumCompThreads(n_threads); end
      input_args.DataCast = char(parameters.simulation.precision);
      input_args_cell = zip_fields(input_args);
      
      if numel(parameters.grid.dims) == 3
         sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args_cell{:});
      elseif numel(parameters.grid.dims) == 2 && isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1
         sensor_data = kspaceFirstOrderAS(kgrid, medium, source, sensor, input_args_cell{:}, 'RadialSymmetry', 'WSWA-FFT');
      else % 2D simulation
         sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args_cell{:});
      end
      
   otherwise

      error('Unsupported code_type: %s. Supported options: matlab_cpu, matlab_gpu, cpp_cpu, cpp_gpu, cpp_mlx', parameters.simulation.code_type);

end

sensor_data = fit_complex(sensor_data, kgrid, p_complex_freqs, p_window);

end

function n = cpu_thread_count()
% Threads for the k-Wave CPU builds: the cores SLURM allocated to the job
% (hpc.cores -> --cpus-per-task). Outside SLURM, or without --cpus-per-task,
% [] keeps k-Wave's and MATLAB's own default.
    n = str2double(getenv('SLURM_CPUS_PER_TASK'));
    if ~(isfinite(n) && n >= 1)
        n = [];
        fprintf('[acoustic_simulation] CPU threads: default\n');
    else
        fprintf('[acoustic_simulation] CPU threads: %d\n', n);
    end
end

function sensor_data = fit_complex(sensor_data, kgrid, p_complex_freqs, p_window)
% Complex amplitude from the recorded window; the time series is dropped.
% With the quarter-period window, p_max_all / p_min_all (bounded by the same
% window) become +/- the fitted amplitude.
    if isempty(p_complex_freqs) || ~isstruct(sensor_data) || ~isfield(sensor_data, 'p')
        return
    end
    sensor_data.p_complex_fit = steady_state_complex_amplitude(sensor_data.p, kgrid.dt, ...
        p_complex_freqs, p_window.record_start_index);
    sensor_data = rmfield(sensor_data, 'p');
    if p_window.p_max_from_fit
        amp = abs(sensor_data.p_complex_fit(:, 1));
        for f = {'p_max_all', 'p_min_all'}
            if ~isfield(sensor_data, f{1}); continue; end
            if numel(sensor_data.(f{1})) ~= numel(amp)
                error('PRESTUS:pComplexMask', ...
                    '%s has %d values but the sensor recorded %d points: the complex field needs a whole-grid sensor mask.', ...
                    f{1}, numel(sensor_data.(f{1})), numel(amp));
            end
            sgn = 1 - 2 * strcmp(f{1}, 'p_min_all');
            sensor_data.(f{1}) = cast(sgn * reshape(amp, size(sensor_data.(f{1}))), class(sensor_data.(f{1})));
        end
        fprintf(['[acoustic_simulation] p_complex: fitted over the last quarter period (%d samples, %.1f per period); ' ...
                 'p_max_all / p_min_all are +/- the fitted amplitude.\n'], p_window.n_samples, p_window.points_per_period);
    end
    sensor_data.p_complex_window = p_window;
end
