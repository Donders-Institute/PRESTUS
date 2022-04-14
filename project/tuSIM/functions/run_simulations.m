    function  [sensor_data] = run_simulations(kgrid, medium, source, sensor, input_args, parameters)

    % select which k-Wave code to run
    %   1: MATLAB CPU code 'matlab_cpu'
    %   2: MATLAB GPU code 'matlab_gpu'
    %   3: C++ code (Interactive) 'cpp_interactive'
    %   4: C++ code (Non-Interactive) 'cpp_noninteractive'
    %   5: CUDA code 'cuda'

    % run code
%     if isfield(parameters,'run_simulations_with_qsub') && parameters.run_simulations_with_qsub == 1
%        disp('Source')
%        disp(source)
%        disp('Medium')
%        disp(medium)
%        disp('Sensor')
%        disp(sensor)
%        disp('Grid')
%        disp(kgrid)
%        disp('Parameters')
%        disp(parameters)
%     end
%     if isfield(parameters,'paths_to_add')
%        path(path, parameters.paths_to_add)
%     end
   
    input_args_cell = zip_fields(input_args);
    if ~parameters.interactive 
        input_args.PlotSim = false;
    end
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
   
   if strcmp(parameters.code_type,'cpp_noninteractive')
       input_args.SaveToDisk = parameters.kwave_input_filename;
   end
   
   input_args_cell = zip_fields(input_args);
   
   % debugging info for non-interactive simulations on the cluster

   
   sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args_cell{:});

   % the code below is for preparing commands for non-interactive
   % processing in terminal
   % currently, not used
%    if strcmp(parameters.code_type,'cpp_noninteractive') || strcmp(parameters.code_type, 'cuda')
%         % display the required syntax to run the C++ simulation
% %         disp('~..~..~..~..~..~.~..~..~..~..~..~.~..~..~..~..~..~.~..~..~..~..~..~.~..~..~..~.'); % line spacing
% %         disp(['Using a terminal window, navigate to the ' filesep 'binaries folder of the k-Wave Toolbox']);
% %         disp('Then, use the syntax shown below to run the simulation:');
% %         disp('||');
% %         disp('||');
% %         disp('\/');
%         if ~isfield(parameters, 'kwave_binaries_path')
%             error("For non-interactive computations, parameters should have a 'kwave_binaries_path' field set")
%         end
%         
%         cmd = '';
%         
%         if strcmp(parameters.code_type,'cpp_noninteractive')
%             kwave_binary = 'kspaceFirstOrder-CUDA';
%         else
%             kwave_binary = 'kspaceFirstOrder-OMP';
%         end
%         if isunix
%             cmd = 'module load gcc; ';
%         else            
%             kwave_binary = [kwave_binary '.exe'];
%         end
%         
%         path_to_kwave = fullfile(parameters.kwave_binaries_path, kwave_binary);
%         
%         cmd = [cmd sprintf('%s -i %s  -o  %s --p_max_all', path_to_kwave, parameters.kwave_input_filename, parameters.kwave_output_filename)];
%         
%         if isunix
%             % 1) make logs directory if does not exist
%     
%             log_dir = fullfile(parameters.data_path, 'batch_job_logs');
%             if ~exist(log_dir, 'file' )
%                 mkdir(log_dir)
%             end
% 
%             subj_id_string = sprintf('sub-%03d', parameters.subject_id);
% 
%             qsub_call = sprintf('qsub -l "nodes=1:gpus=1,feature=cuda,walltime=06:00:00,mem=20gb,reqattr=cudacap>=5.0" -o %s -e %s -w %s', ...
%                 fullfile(log_dir, sprintf('%s_qsub_output.log', subj_id_string)),...
%                 fullfile(log_dir, sprintf('%s_qsub_error.log', subj_id_string)),...
%                 data_path);
% 
%             full_cmd = sprintf('echo "cd %s; %s" | %s', data_path, cmd, qsub_call);
% 
%             fprintf('Running kWave with command \n%s', full_cmd)
%             [res, out] = system(full_cmd)
%             
%             fprintf('Now wait for the job with the id listed above to finish')
%             exit();
%         else
%             fprintf('Running the kwave command %s now...\n', cmd)
%             [res, out] = system(cmd)
%             fprintf('Trying to get the outputs from %s\n', parameters.kwave_output_filename)
%             sensor_data.p_max_all = h5read(parameters.kwave_output_filename, '/p_max_all');
%         end
%         
%     end 
end
