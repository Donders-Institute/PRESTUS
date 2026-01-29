function [sensor_data] = acoustic_wrapper(parameters, kgrid, kwave_medium, source, sensor, medium_masks, filename_sensor_data)
    
    disp('Specifying and Starting acoustic simulations...')

    % Pathname for the input and output files (used only for non-interactive computations)
    parameters.kwave_input_filename  = fullfile(parameters.output_dir, ...
        sprintf('sub-%03d_%s_input%s.h5', parameters.subject_id, ...
        parameters.simulation_medium, parameters.results_filename_affix));
    parameters.kwave_output_filename = fullfile(parameters.output_dir, ...
        sprintf('sub-%03d_%s_output%s.h5', parameters.subject_id, ...
        parameters.simulation_medium, parameters.results_filename_affix));

    % Note: The edge of the simulation grid = the edge of the PML layer
    plminside = true;

    kwave_input_args = struct('PMLInside', plminside, ...
        'PMLSize', parameters.pml_size, ...
        'PlotPML', true);

    if contains(parameters.simulation_medium, {'layered'}) && ...
            any(ismember(fieldnames(parameters.layers), {'skull'}))
        % Extract the skull edge ...
        mask = tissuemask_binary(parameters, medium_masks);
        skull_edge = edge3(ismember(medium_masks, mask.skull), 'approxcanny', 0.1);
        % ... to set as display mask
        kwave_input_args.DisplayMask = skull_edge;
    end

    if parameters.run_source_setup==0
        error('Source setup not requested. Not able to proceed with acoustic simulation.')
    end

    sensor_data = acoustic_simulation(kgrid, kwave_medium, source, sensor, ...
        kwave_input_args, parameters);

    % keep 'parameters' in info so not to confuse future runs when saving
    acoustic_info.parameters = parameters;
    acoustic_info.kwave_input_args = kwave_input_args;

    if isfield(parameters, 'savemat') && parameters.savemat==0
        disp("Not saving acoustic output matrices ...")
    else
        % keep 'parameters' as a copy so not to confuse future runs
        acoustic_info.parameters = parameters;
        acoustic_info.kwave_input_args = kwave_input_args;
        save(filename_sensor_data, ...
            'sensor_data', ...
            'kgrid', ...
            'kwave_medium', ...
            'source', ...
            'sensor', ...
            'acoustic_info' ,'-v7.3')
    end