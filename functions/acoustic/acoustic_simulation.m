function [sensor_data] = acoustic_simulation(kgrid, medium, source, sensor, input_args, parameters)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                  Run acoustic simulation                          %
%                                                                   %
% Loads the parameters into a k-wave function.                      %
% Changes some input parameters based on what cpu, gpu and which    %
% which mode it is run in (interactive vs. non-interactive)         %
%                                                                   %
% Some notes:                                                       %
% - This script only loads the parameters into k-wave, for details  %
% on the simulations itself, look at k-wave documentation for the   %
% functions 'kspaceFirstOrder3D' and 'kspaceFirstOrder2D'.          %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
% Select which k-Wave code to run
%   1: MATLAB CPU code 'matlab_cpu'
%   2: MATLAB GPU code 'matlab_gpu'
%   3: C++ CPU code 'cpp_cpu'
%   4: C++ GPU code 'cpp_gpu'

% Only actively plot each timepoint of the simulations if it is interactive
if ~parameters.interactive 
   input_args.PlotSim = false;
end

% Remove fields that are not recognized for acoustic simulations
medium = rmfield(medium,'thermal_conductivity');
medium = rmfield(medium,'specific_heat');
medium = rmfield(medium,'perfusion_coeff');

% Define scale for plotting
input_args.PlotScale = [-1, 1] * parameters.transducer(1).source_amp(1);

% Select submission based on code type
switch parameters.code_type
    case 'cpp_cpu'

      % Force precision for C++ HDF5 compatibility
      medium  = cast_struct(medium, parameters.precision); 
      source  = cast_struct(source, parameters.precision);
      sensor  = cast_struct(sensor, parameters.precision);
   
      % Pathname for the input and output files (used only for non-interactive computations)
      input_args.SaveToDisk = fullfile(parameters.output_dir, ...
         sprintf('sub-%03d_%s_input%s.h5', parameters.subject_id, ...
         parameters.simulation_medium, parameters.results_filename_affix));
      
      input_args.DataName = sprintf('kwave_sub-%03d%s', parameters.subject_id, parameters.results_filename_affix);
      input_args.DataPath = parameters.output_dir;
      input_args.DeleteData = true;

      if parameters.n_sim_dims == 3
         input_args.FunctionName = 'kspaceFirstOrder3D';
      elseif parameters.n_sim_dims == 2 && isfield(parameters, 'axisymmetric') && parameters.axisymmetric == 1
         input_args.FunctionName = 'kspaceFirstOrderAS';
      else
         input_args.FunctionName = 'kspaceFirstOrder2D';
      end

      input_args_cell = zip_fields(input_args);
      sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args_cell{:});
      
   case 'cpp_gpu'

      % Force precision for C++ HDF5 compatibility
      medium  = cast_struct(medium, parameters.precision); 
      source  = cast_struct(source, parameters.precision);
      sensor  = cast_struct(sensor, parameters.precision);

      % Pathname for the input and output files (used only for non-interactive computations)
      input_args.SaveToDisk = fullfile(parameters.output_dir, ...
        sprintf('sub-%03d_%s_input%s.h5', parameters.subject_id, ...
        parameters.simulation_medium, parameters.results_filename_affix));

      input_args.DataName = sprintf('kwave_sub-%03d%s', parameters.subject_id, parameters.results_filename_affix);
      input_args.DataPath = parameters.output_dir;
      input_args.DeleteData = true;

      if parameters.n_sim_dims == 3
         input_args.FunctionName = 'kspaceFirstOrder3D';
      elseif parameters.n_sim_dims == 2 && isfield(parameters, 'axisymmetric') && parameters.axisymmetric == 1
         input_args.FunctionName = 'kspaceFirstOrderAS';
      else
         input_args.FunctionName = 'kspaceFirstOrder2D';
      end

      gpu_id = str2double(getenv('SLURM_LOCALID'));
      input_args.DeviceNum = gpu_id;
      input_args.NumThreads = 1;

      input_args_cell = zip_fields(input_args);
      sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args_cell{:});
      
   case 'matlab_gpu'

      input_args.DataCast = ['gpuArray-', char(parameters.precision)];
      input_args_cell = zip_fields(input_args);
      
      if parameters.n_sim_dims == 3
         sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args_cell{:});
      elseif parameters.n_sim_dims == 2 && isfield(parameters, 'axisymmetric') && parameters.axisymmetric == 1
         sensor_data = kspaceFirstOrderAS(kgrid, medium, source, sensor, input_args_cell{:}, 'RadialSymmetry', 'WSWA-FFT');
      else % 2D simulation
         sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args_cell{:});
      end
      
   case 'matlab_cpu'

      input_args.DataCast = parameters.precision;
      input_args_cell = zip_fields(input_args);
      
      if parameters.n_sim_dims == 3
         sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args_cell{:});
      elseif parameters.n_sim_dims == 2 && isfield(parameters, 'axisymmetric') && parameters.axisymmetric == 1
         sensor_data = kspaceFirstOrderAS(kgrid, medium, source, sensor, input_args_cell{:}, 'RadialSymmetry', 'WSWA-FFT');
      else % 2D simulation
         sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args_cell{:});
      end
      
   otherwise

      error('Unsupported code_type: %s. Supported options: matlab_cpu, matlab_gpu, cpp_cpu, cpp_gpu', parameters.code_type);

end

end