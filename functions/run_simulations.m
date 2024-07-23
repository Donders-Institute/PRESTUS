    function  [sensor_data] = run_simulations(kgrid, medium, source, sensor, input_args, parameters)

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %                           Run simulations                         %
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
    %   3: C++ code (Interactive) 'cpp_interactive'
    %   4: C++ code (Non-Interactive) 'cpp_noninteractive'
    %   5: CUDA code 'cuda'
   
    % Only actively plot each timepoint of the simulations if it is interactive
    input_args_cell = zip_fields(input_args);
    if ~parameters.interactive 
        input_args.PlotSim = false;
    end

    % Runs simulations in C++ in interactive mode or on
    % non-interactive code for Nvidia gpu's.
    if strcmp(parameters.code_type, 'cpp_interactive')
       sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args_cell{:});
    elseif strcmp(parameters.code_type, 'cuda')
       sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args_cell{:});
    else
        
    input_args.PlotScale = [-1, 1] * parameters.transducer.source_amp(1);
   
   if strcmp(parameters.code_type,'matlab_gpu')
      input_args.DataCast = 'gpuArray-single';
   else
      input_args.DataCast = 'single';       
   end
   
   % If using C++ in non-interactive mode without an Nvidia gpu
   if strcmp(parameters.code_type,'cpp_noninteractive')
       input_args.SaveToDisk = parameters.kwave_input_filename;
   end
   
   % Save new input arguments for simulations
   input_args_cell = zip_fields(input_args);
   
   % Debugging info for non-interactive simulations on the Donders cluster
   medium = rmfield(medium,'thermal_conductivity');
   medium = rmfield(medium,'specific_heat');

   % Runs simulations on the CPU only in 3 or 2 dimensions
   if parameters.n_sim_dims == 3
       sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args_cell{:});
   else
       sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args_cell{:});
   end

end