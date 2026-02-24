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

      if parameters.n_sim_dims ~= 3
         error("C++ option only supported for 3D acoustic simulations. Please choose a different code_type.");
      end

      input_args_cell = zip_fields(input_args);
      sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args_cell{:});
      
   case 'cpp_gpu'

      if parameters.n_sim_dims ~= 3
         error("CUDA option only supported for 3D acoustic simulations. Please choose a different code_type.");
      end
      input_args_cell = zip_fields(input_args);
      sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args_cell{:});
      
   case 'matlab_gpu'

      input_args.DataCast = 'gpuArray-single';
      input_args_cell = zip_fields(input_args);
      
      if parameters.n_sim_dims == 3
         sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args_cell{:});
      elseif parameters.n_sim_dims == 2 && isfield(parameters, 'axisymmetric') && parameters.axisymmetric == 1
         sensor_data = kspaceFirstOrderAS(kgrid, medium, source, sensor, input_args_cell{:}, 'RadialSymmetry', 'WSWA-FFT');
      else % 2D simulation
         sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args_cell{:});
      end
      
   case 'matlab_cpu'

      input_args.DataCast = 'single';
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