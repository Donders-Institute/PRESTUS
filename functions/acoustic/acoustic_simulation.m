function [sensor_data] = acoustic_simulation(kgrid, medium, source, sensor, input_args, parameters)
% ACOUSTIC_SIMULATION  Dispatch a k-Wave acoustic simulation across supported compute backends
%
% Selects and calls the appropriate k-Wave solver function based on
% parameters.simulation.code_type ('matlab_cpu', 'matlab_gpu', 'cpp_cpu',
% 'cpp_gpu') and grid dimensionality (2D, 2D axisymmetric, or 3D). Strips
% thermal-only medium fields before passing to k-Wave, sets DataCast/precision
% for GPU/C++ paths, and configures HDF5 file I/O for C++ backends.
%
% Use as:
%   sensor_data = acoustic_simulation(kgrid, medium, source, sensor, input_args, parameters)
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
%                grid.dims, io.output_dir, transducer, simulation.precision
%
% Output:
%   sensor_data - k-Wave sensor output (p_max_all, p_final, etc.) [Pa]
%
% See also: ACOUSTIC_WRAPPER, KSPACEFIRSTORDER3D, KSPACEFIRSTORDER2D, KSPACEFIRSTORDERAS

arguments
    kgrid      (1,1)
    medium     (1,1) struct
    source     (1,1) struct
    sensor     (1,1) struct
    input_args (1,1) struct
    parameters (1,1) struct
end

% Select which k-Wave code to run
%   1: MATLAB CPU code 'matlab_cpu'
%   2: MATLAB GPU code 'matlab_gpu'
%   3: C++ CPU code 'cpp_cpu'
%   4: C++ GPU code 'cpp_gpu'

% Only actively plot each timepoint of the simulations if it is interactive
if ~parameters.simulation.interactive
   input_args.PlotSim = false;
end

% Remove fields that are not recognized for acoustic simulations
thermal_fields = {'thermal_conductivity', 'specific_heat', 'perfusion_coeff'};
for f = thermal_fields
    if isfield(medium, f{1})
        medium = rmfield(medium, f{1});
    end
end

% Define scale for plotting
input_args.PlotScale = [-1, 1] * parameters.transducer(1).(parameters.transducer(1).type).elem_amp(1);

% Select submission based on code type
switch parameters.simulation.code_type
    case 'cpp_cpu'

      % Force precision for C++ HDF5 compatibility
      medium  = cast_struct(medium, parameters.simulation.precision);
      source  = cast_struct(source, parameters.simulation.precision);
      sensor  = cast_struct(sensor, parameters.simulation.precision);

      % Pathname for the input and output files (used only for non-interactive computations)
      h5_dir = char(parameters.io.cache_dir);
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

      input_args_cell = zip_fields(input_args);
      sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args_cell{:});

   case 'cpp_gpu'

      % Force precision for C++ HDF5 compatibility
      medium  = cast_struct(medium, parameters.simulation.precision);
      source  = cast_struct(source, parameters.simulation.precision);
      sensor  = cast_struct(sensor, parameters.simulation.precision);

      input_args_cell = zip_fields(input_args);

      if numel(parameters.grid.dims) == 2 && isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1
         % No CUDA binary supports axisymmetric geometry; fall back to matlab_gpu path
         % (kspaceFirstOrderAS with GPU data cast) which preserves GPU acceleration.
         warning('kspaceFirstOrder-CUDA does not support axisymmetric simulations. Falling back to matlab_gpu (kspaceFirstOrderAS with gpuArray).');
         as_args = zip_fields(rmfield(input_args, intersect(fieldnames(input_args), ...
             {'SaveToDisk','DataName','DataPath','DeleteData','FunctionName','DeviceNum','NumThreads'})));
         as_args_cast = [as_args, {'DataCast', ['gpuArray-', char(parameters.simulation.precision)]}];
         sensor_data = kspaceFirstOrderAS(kgrid, medium, source, sensor, as_args_cast{:}, 'RadialSymmetry', 'WSWA-FFT');
      else
         % Pathname for the input and output files
         h5_dir = char(parameters.io.cache_dir);
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

      error('Unsupported code_type: %s. Supported options: matlab_cpu, matlab_gpu, cpp_cpu, cpp_gpu', parameters.simulation.code_type);

end

end